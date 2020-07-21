#include "fusion_io.h"
#include "interpolate.h"

#include <iostream>
#include <fstream>
#include <iomanip>

static double pow(double x, int p)
{
  double r = 1.;
  if(p < 0) return 0.;
  for(int i=0; i<p; i++) r = r*x;
  return r;
}


gpec_source::gpec_source()
  : b0(false, 3), b1(true, 3), bx(true, 3)
{
}

gpec_source::~gpec_source()
{
}

int gpec_source::gpec_field_data::allocate(const int nr0, const int nz0)
{
  nr = nr0;
  nz = nz0;

  r = new double[nr];
  z = new double[nz];
  v_real = new double**[n_comp];
  if(is_complex) v_imag = new double**[n_comp];
  for(int k=0; k<n_comp; k++) {
    v_real[k] = new double*[nr];
    if(is_complex) v_imag[k] = new double*[nr];
    for(int i=0; i<nr; i++) {
      v_real[k][i] = new double[nz];
      if(is_complex) v_imag[k][i] = new double[nz];
    }
  }
  return FIO_SUCCESS;
}

int gpec_source::gpec_field_data::free()
{
  if(!is_populated)
    return FIO_SUCCESS;

  delete[] r;
  delete[] z;

  for(int k=0; k<n_comp; k++) {
    for(int i=0; i<nr; i++) {
      delete[] v_real[k][i];
      if(is_complex) delete[] v_imag[k][i];
    }
    delete[] v_real[k];
    if(is_complex) delete[] v_imag[k];
  }
  delete[] v_real;
  if(is_complex) delete[] v_imag;

  r = 0;
  z = 0;
  v_real = 0;
  v_imag = 0;
}

int gpec_source::read_field_data(const char* filename, gpec_field_data* d)
{
  if(d->is_populated)
    return FIO_SUCCESS;

  std::ifstream file;
  std::string dum;

  dum = path + "/" + filename;
  file.open(dum.c_str(), std::fstream::in);

  if(!file) {
    std::cerr << "Error: cannot open file " << filename << std::endl;
    return FIO_FILE_ERROR;
  }

  int nr, nz;

  std::getline(file, dum);  std::cerr << dum << std::endl;
  std::getline(file, dum);  std::cerr << dum << std::endl;
  std::getline(file, dum);  std::cerr << dum << std::endl;
  file >> dum >> dum >> d->ntor;
  std::cerr << "Toroidal Mode Number = " << d->ntor << std::endl;
  file >> dum >> dum >> nr >> dum >> dum >> nz;
  std::cerr << "Dimensions: " << nr << " x " << nz << std::endl;

  std::getline(file, dum); std::cerr << dum << std::endl;
  std::getline(file, dum); std::cerr << dum << std::endl;
  std::getline(file, dum); std::cerr << dum << std::endl;

  int result;
  std::cerr << "Allocating gpec field data" << std::endl;
  result = d->allocate(nr, nz);
  if(result != FIO_SUCCESS) return result;

  std::cerr << "Reading gpec field data" << std::endl;
  for(int i=0; i<nr; i++) {
    for(int j=0; j<nz; j++) {
      file >> dum >> d->r[i] >> d->z[j];
      for(int k=0; k<d->n_comp; k++) {
	file >> d->v_real[k][i][j];
	if(d->is_complex)
	  file >> d->v_imag[k][i][j];
      }
    }
  }

  std::cerr << "nr = " << nr << std::endl;
  std::cerr << "nz = " << nz << std::endl;
  std::cerr << d->r[0] << " < r < " << d->r[d->nr-1] << std::endl;
  std::cerr << d->z[0] << " < z < " << d->z[d->nz-1] << std::endl;

  file.close();

  return 0;
}

int gpec_source::open(const char* p)
{
  path = p;
  return FIO_SUCCESS;
}

int gpec_source::close()
{
  return FIO_SUCCESS;
}

int gpec_source::center(double *r0, double *z0) const
{
  return FIO_SUCCESS;
}

int gpec_source::extent(double* r0, double* r1, double* z0, double* z1)
  const
{
  return FIO_SUCCESS;
}

int gpec_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();
  opt->add_option(FIO_LINEAR_SCALE, 1.);
  opt->add_option(FIO_PHASE, 0.);
  opt->add_option(FIO_PART, FIO_TOTAL);
  opt->add_option(FIO_TIMESLICE, 1);

  return FIO_SUCCESS;
}

int gpec_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_MAGNETIC_FIELD);

  return FIO_SUCCESS;
}

int gpec_source::get_field(const field_type t,fio_field** f,
			   const fio_option_list* opt)
{
  gpec_field* gf;
  *f = 0;

  switch(t) {

  case(FIO_MAGNETIC_FIELD):
    gf = new gpec_magnetic_field(this);
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  int result = gf->load(opt);
  if(result != FIO_SUCCESS) {
    delete(gf);
    return result;
  }
 
  *f = gf;

  return FIO_SUCCESS;
}

int gpec_source::get_series(const series_type t,fio_series** s)
{
  switch(t) {
  case(FIO_MAGAXIS_R):
    *s = new fio_scalar_series(1.875);
    break;

  case(FIO_MAGAXIS_Z):
    *s = new fio_scalar_series(0.0);
    break;

  default:
    *s = 0;
    return FIO_UNSUPPORTED;
  };
  
  return FIO_SUCCESS;
}


int gpec_source::gpec_field_data::interpolate(const double r0, const double z0,
					      double* vr, double* vi) const
{
  double a[4][4];

  double dx = (r[nr-1]-r[0]) / (nr - 1.);
  double dz = (z[nr-1]-z[0]) / (nz - 1.);

  double p = (r0 - r[0])/dx;
  double q = (z0 - z[0])/dz;
  int i = (int)p;
  int j = (int)q;

  if(i < 1 || i > nr) return FIO_OUT_OF_BOUNDS;
  if(j < 1 || j > nz) return FIO_OUT_OF_BOUNDS;

  /*
  for(int k=0; k<n_comp; k++) {
    vr[k] = v_real[k][i][j];
    if(is_complex) vi[k] = v_imag[k][i][j];
  }
  */

  bool result;
  double si[6];
  bool imag = false;
  int loops = is_complex ? 1 : 0;

  for(int l=0; l<=loops; l++) {
    for(int k=0; k<n_comp; k++) {
      if(l==0) {
	result = 
	  bicubic_interpolation_coeffs((const double**)(v_real[k]),nr,nz,p,q,a);
      } else { 
	result = 
	  bicubic_interpolation_coeffs((const double**)(v_imag[k]),nr,nz,p,q,a);
      }
      if(!result)
	return FIO_OUT_OF_BOUNDS;
      
      for(int n=0; n<6; n++)
	si[n] = 0.;
      
      //double temp;
      for(int n=0; n<4; n++) {
	for(int m=0; m<4; m++) {
	  si[0] += a[m][n]      *pow(p-i,n  )*pow(q-j,m  );
	  /*	  
	  temp = a[m][n]      *n*pow(p-i,n-1)*pow(q-j,m  );
	  si[1] += temp/dx;
	  
	  temp = a[m][n]      *m*pow(p-i,n  )*pow(q-j,m-1);
	  si[2] += temp/dz;
	  
	  temp = a[m][n]*n*(n-1)*pow(p-i,n-2)*pow(q-j,m  );
	  si[3] += temp/(dx*dx);
	  
	  temp = a[m][n]    *m*n*pow(p-i,n-1)*pow(q-j,m-1);
	  si[4] += temp/(dx*dz);
	  
	  temp = a[m][n]*m*(m-1)*pow(p-i,n  )*pow(q-j,m-2);
	  si[5] += temp/(dz*dz);
	  */
	}
      }
      if(l==0) {
	vr[k] = si[0];
      } else {
	vi[k] = si[0];
      }    
    }
  }

  return FIO_SUCCESS;
}
