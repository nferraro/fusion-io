#ifdef FUSIONIO_ENABLE_NIMROD

#include "nimrod_source.h"
#include "nimrod_field.h"
#include "fusion_io.h"
#include <iostream>

extern "C" {
  void c_nimfield_init(const char*, const char*);
}

nimrod_source::nimrod_source()
{
}

nimrod_source::~nimrod_source()
{
}

int nimrod_source::open(const char* filename)
{
  dump_file = filename;
  // Initialize nimfield with all possible fields we might want
  // b=magnetic field, v=velocity, n=density, t=temperature, p=pressure
  c_nimfield_init(filename, "bvntp");
  return FIO_SUCCESS;
}

int nimrod_source::close()
{
  return FIO_SUCCESS;
}

int nimrod_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_FLUID_VELOCITY);
  fields->push_back(FIO_DENSITY);
  fields->push_back(FIO_PRESSURE);
  fields->push_back(FIO_TEMPERATURE);
  return FIO_SUCCESS;
}

int nimrod_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();
  opt->add_option(FIO_PART, FIO_TOTAL);
  opt->add_option(FIO_SPECIES, FIO_MAIN_ION);
  return FIO_SUCCESS;
}

int nimrod_source::get_field(const field_type t, fio_field** f, const fio_option_list* opt)
{
  int part = FIO_TOTAL;
  if(opt) opt->get_option(FIO_PART, &part);

  int species = FIO_MAIN_ION;
  if(opt) opt->get_option(FIO_SPECIES, &species);

  int eq;
  switch(part) {
    case FIO_TOTAL:             eq = 2; break;
    case FIO_PERTURBED_ONLY:    eq = 0; break;
    case FIO_EQUILIBRIUM_ONLY:  eq = 1; break;
    case FIO_AXISYMMETRIC_ONLY: eq = 3; break;
    default:
      std::cerr << "Error: Unsupported NIMROD FIO_PART." << std::endl;
      return FIO_UNSUPPORTED;
      break;
  }

  switch(t) {
    case FIO_MAGNETIC_FIELD:
      *f = new nimrod_field("b", 3, eq);
      break;
    case FIO_FLUID_VELOCITY:
      *f = new nimrod_field("v", 3, eq);
      break;
    case FIO_DENSITY:
      *f = new nimrod_field("n", 1, eq);
      break;
    case FIO_PRESSURE:
      if(species == FIO_ELECTRON)
        *f = new nimrod_field("pe", 1, eq);
      else if(species == FIO_MAIN_ION)
        *f = new nimrod_field("pi", 1, eq);
      else
        *f = new nimrod_field("p", 1, eq);
      break;
    case FIO_TEMPERATURE:
      if(species == FIO_ELECTRON)
        *f = new nimrod_field("te", 1, eq);
      else
        *f = new nimrod_field("ti", 1, eq);
      break;
    default:
      std::cerr << "Error: Unsupported NIMROD field." << std::endl;
      return FIO_UNSUPPORTED;
  }

  return FIO_SUCCESS;
}

int nimrod_source::get_int_parameter(const parameter_type t, int* p) const
{
  switch(t) {
    case FIO_GEOMETRY:
      *p = FIO_CYLINDRICAL;
      return FIO_SUCCESS;
    default:
      return FIO_UNSUPPORTED;
  }
}

int nimrod_source::get_real_parameter(const parameter_type t, double* p) const
{
  switch(t) {
    case FIO_PERIOD:
      *p = 2.*M_PI;
      return FIO_SUCCESS;
    default:
      return FIO_UNSUPPORTED;
  }
}

int nimrod_source::allocate_search_hint(fio_hint* s)
{
  double* h = new double[2];
  h[0] = -1.0;
  h[1] = -1.0;
  *s = (fio_hint)h;
  return FIO_SUCCESS;
}

int nimrod_source::deallocate_search_hint(fio_hint* s)
{
  double* h = (double*)(*s);
  delete[] h;
  *s = 0;
  return FIO_SUCCESS;
}

#endif // FUSIONIO_ENABLE_NIMROD
