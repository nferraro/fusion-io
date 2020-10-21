#include "fusion_io.h"
#include <iostream>

int m3dc1_source::open(const char* filename)
{
  if(!file.open(filename))
    return FIO_FILE_ERROR;

  // read some relevant data from the file
  file.read_parameter("extsubtract", &extsubtract);
  file.read_parameter("icomplex", &icomplex);
  file.read_parameter("3d", &i3d);
  file.read_parameter("eqsubtract", &eqsubtract);
  file.read_parameter("linear", &linear);
  file.read_parameter("ion_mass", &ion_mass);
  file.read_parameter("bzero", &bzero);
  file.read_parameter("rzero", &rzero);
  file.read_parameter("n0_norm", &n0);
  file.read_parameter("l0_norm", &L0);
  file.read_parameter("b0_norm", &B0);
  file.read_parameter("version", &version);
  file.read_parameter("itor", &itor);
  file.read_parameter("ntime", &ntime);
  file.read_parameter("ntor", &ntor);

  if(version < 23) 
    file.read_parameter("zeff", &z_ion);
  else
    file.read_parameter("z_ion", &z_ion);

  if(version >= 23)
    file.read_parameter("kprad_z", &kprad_z);
  else
    kprad_z = 0;

  const double c = 3.e10;
  const double m_p = 1.6726e-24;

  // define some normalization quantities in cgs
  p0 = B0*B0/(4.*M_PI);
  J0 = B0*c/(4.*M_PI*L0);
  v0 = B0/sqrt(4.*M_PI*ion_mass*m_p*n0);
  t0 = L0/v0;
  Phi0 = L0*v0*B0/c;
  temp0 = p0/n0;

  // convert normalization quantities to mks
  n0 /= 1.e-6;
  B0 /= 1e4;
  L0 /= 100.;
  p0 /= 10.;
  J0 /= (c*1e-5);
  v0 /= 100.;
  Phi0 /= (1.e8/c);
  temp0 /= 1.6022e-12;    // Temperature is output in eV

  if(itor==1) {
    period = 2.*M_PI;
  } else {
    period = 2.*M_PI*rzero*L0;
  }

  // determine ion species (assume one proton and no electrons)
  ion_species = fio_species(ion_mass, 1, 0);

  return FIO_SUCCESS;
}

int m3dc1_source::allocate_search_hint(void** s)
{
  *s = new int;
  *((int*)(*s)) = -1;
  //  std::cerr << "Allocated hint at " << *s << std::endl;
  return FIO_SUCCESS;
}
int m3dc1_source::deallocate_search_hint(void** s)
{
  delete (int*)(*s);
  return FIO_SUCCESS;
}


int m3dc1_source::close()
{
  if(!file.close())
    return FIO_FILE_ERROR;
  return FIO_SUCCESS;
}

int m3dc1_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  opt->add_option(FIO_TIMESLICE, 0);
  opt->add_option(FIO_LINEAR_SCALE, 1.);
  opt->add_option(FIO_PHASE, 0.);
  opt->add_option(FIO_PART, FIO_TOTAL);
  opt->add_option(FIO_SPECIES, 0);

  return FIO_SUCCESS;
}

int m3dc1_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_DENSITY);
  fields->push_back(FIO_CURRENT_DENSITY);
  fields->push_back(FIO_FLUID_VELOCITY);
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_POLOIDAL_FLUX);
  fields->push_back(FIO_POLOIDAL_FLUX_NORM);
  fields->push_back(FIO_PRESSURE);
  fields->push_back(FIO_TEMPERATURE);
  fields->push_back(FIO_TOTAL_PRESSURE);
  fields->push_back(FIO_VECTOR_POTENTIAL);
  fields->push_back(FIO_ELECTRIC_FIELD);

  return FIO_SUCCESS;
}

int m3dc1_source::get_int_parameter(const parameter_type t, int* p) const
{
  switch(t) {

  case FIO_NUM_TIMESLICES:
    *p = ntime;
    return FIO_SUCCESS;

  case FIO_GEOMETRY:
    *p = (itor==1 ? FIO_CYLINDRICAL : FIO_CARTESIAN);
    return FIO_SUCCESS;

  case FIO_TOROIDAL_MODE:
    *p = ntor;
    return FIO_SUCCESS;

  default:    return FIO_UNSUPPORTED;
  }
}

int m3dc1_source::get_real_parameter(const parameter_type t, double* p) const
{
  switch(t) {

  case FIO_PERIOD:
    *p = period;
    return FIO_SUCCESS;

  default:    return FIO_UNSUPPORTED;
  }
}

int m3dc1_source::get_field(const field_type t, fio_field** f,
			    const fio_option_list* opt)
{
  *f = 0;
  m3dc1_fio_field* mf(0);
  fio_series *psi0, *psi1;
  double psi0_val, psi1_val;
  bool unneeded_species = false;
  int s, result;

  result = FIO_SUCCESS;

  opt->get_option(FIO_SPECIES, &s);
  if(s==FIO_MAIN_ION) s = ion_species;

  switch(t) {
  case(FIO_DENSITY):
    if(s==fio_electron) {
      if(version < 23) 
	mf = new m3dc1_scalar_field(this, "den", n0*z_ion);
      else
	mf = new m3dc1_scalar_field(this, "ne", n0);
    } else if(s==ion_species) {
      mf = new m3dc1_scalar_field(this, "den", n0);
    } else {
      if((fio_species(s).atomic_number() != kprad_z) || (kprad_z==0)) {
	std::cerr << "Error: requesting density of species with atomic number "
		  << fio_species(s).atomic_number()
		  << ".  Impurity species in simulation has "
		  << "atomic number " << kprad_z << std::endl;
	result = FIO_BAD_SPECIES;
      } else {
	char name[11];
	snprintf(name, 11, "kprad_n_%02d", fio_species(s).charge());
	mf = new m3dc1_scalar_field(this, name, n0);
      }
    }
    break;

  case(FIO_TEMPERATURE):
    if(s==fio_electron) {
      mf = new m3dc1_scalar_field(this, "te", temp0);
    } else if(s==ion_species) {
      mf = new m3dc1_scalar_field(this, "ti", temp0);
    } else {
      result = FIO_BAD_SPECIES;
    }
    break;

  case(FIO_CURRENT_DENSITY):
    mf = new m3dc1_current_density(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_ELECTRIC_FIELD):
    mf = new m3dc1_electric_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_FLUID_VELOCITY):
    mf = new m3dc1_fluid_velocity(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_MAGNETIC_FIELD):
    mf = new m3dc1_magnetic_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_POLOIDAL_FLUX):
    mf = new m3dc1_scalar_field(this, "psi", 2.*M_PI*B0*L0*L0);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_POLOIDAL_FLUX_NORM):
    int t;
    double slice_time;
    opt->get_option(FIO_TIMESLICE, &t);
    if((result = get_slice_time(t, &slice_time)) != FIO_SUCCESS)
      break;

    if((result = get_series(FIO_MAGAXIS_PSI, &psi0)) != FIO_SUCCESS)
      break;
    if((result = get_series(FIO_LCFS_PSI,    &psi1)) != FIO_SUCCESS)
      break;
    if((result = psi0->eval(slice_time, &psi0_val)) != FIO_SUCCESS)
       break;
    if((result = psi1->eval(slice_time, &psi1_val)) != FIO_SUCCESS)
       break;
    delete(psi0);
    delete(psi1);

    double offset;
    int part;
    opt->get_option(FIO_PART, &part);
    if(part==FIO_PERTURBED_ONLY) {
      offset = 0.;
    } else {
      offset = psi0_val;
    }
    mf = new m3dc1_scalar_field(this, "psi",
				1./(psi1_val - psi0_val), 
				offset);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_PRESSURE):
    if(s==fio_electron) {
      mf = new m3dc1_scalar_field(this, "Pe", p0);
    } else if(s==ion_species) {
      mf = new m3dc1_pi_field(this);
    } else {
      result = FIO_BAD_SPECIES;
    }
    break;

  case(FIO_SCALAR_POTENTIAL):
    mf = new m3dc1_phi_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_TOTAL_PRESSURE):
    mf = new m3dc1_scalar_field(this, "P", p0);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_VECTOR_POTENTIAL):
    mf = new m3dc1_vector_potential(this);
    if(s!=0) unneeded_species = true;
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  if(result==FIO_BAD_SPECIES) {
    std::cerr << "Unsupported species: " << fio_species(s).name() << std::endl;
    std::cerr << "Main ions: " << ion_species.name() << std::endl;
    return result;
  }

  result = mf->load(opt);
  if(result == FIO_SUCCESS) {
    *f = mf;
  } else {
    delete(mf);
  }
  return result;
}

int m3dc1_source::get_series(const series_type t,fio_series** s)
{
  m3dc1_fio_series* ms;
  int result;

  *s = 0;

  result = FIO_SUCCESS;

  switch(t) {
  case(FIO_MAGAXIS_PSI):
    ms = new m3dc1_fio_series(this, "psimin", B0*L0*L0);
    break;

  case(FIO_LCFS_PSI):
    ms = new m3dc1_fio_series(this, "psi_lcfs", B0*L0*L0);
    break;

  case(FIO_MAGAXIS_R):
    ms = new m3dc1_fio_series(this, "xmag", L0);
    break;

  case(FIO_MAGAXIS_Z):
    ms = new m3dc1_fio_series(this, "zmag", L0);
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  result = ms->load();
  if(result == FIO_SUCCESS) {
    *s = ms;
  } else {
    delete(ms);
  }
  return result;
}

int m3dc1_source::get_slice_time(const int slice, double* time)
{
  if(!file.get_slice_time(slice, time)) {
    return FIO_OUT_OF_BOUNDS;
  }
  *time *= t0;
  return FIO_SUCCESS;
}

