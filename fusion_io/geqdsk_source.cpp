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


geqdsk_source::geqdsk_source()
  : psi(0), psirz(0), fpol(0), ffprime(0), press(0)
{
}

geqdsk_source::~geqdsk_source()
{
  close();
}

int geqdsk_source::open(const char* filename)
{
  std::fstream gfile;
  std::string dum;

  // Read eqdsk file
  gfile.open(filename, std::fstream::in);

  if(!gfile) {
    std::cerr << "Error: cannot open file " << filename << std::endl;
    return FIO_FILE_ERROR;
  }

  double rdim, zdim, rcentr;
  double bcentr, current;
    
  for(int i=0; i<5; i++) gfile >> dum;
  gfile >> std::setw(16) >> nw >> nh;
  gfile >> std::setw(16) >> rdim >> zdim >> rcentr >> rleft >> zmid;
  gfile >> std::setw(16) >> rmaxis >> zmaxis >> simag >> sibry >> bcentr;
  gfile >> std::setw(16) >> current >> simag >> dum >> rmaxis >> dum;
  gfile >> std::setw(16) >> zmaxis >> dum >> sibry >> dum >> dum;
  zbottom = zmid - zdim/2.;

  std::cout << "nw = " << nw << ", nh = " << nh << std::endl;
  std::cout << "rmaxis = " << rmaxis << ", zmaxis = " << zmaxis << std::endl;
  std::cout << "rleft = " << rleft << ", zmid = " << zmid << std::endl;
  std::cout << "rdim = " << rdim << ", zdim = " << zdim << std::endl;
  std::cout << "simag = " << simag << std::endl;
  std::cout << "dum = " << dum << std::endl;

  psirz = new double*[nh];
  for(int i=0; i<nh; i++)
    psirz[i] = new double[nw];
  ffprime = new double[nw];
  fpol = new double[nw];
  press = new double[nw];

  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> fpol[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> press[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> ffprime[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // pprime
  for(int i=0; i<nh; i++)
    for(int j=0; j<nw; j++)
      gfile >> std::setw(16) >> psirz[i][j];

  gfile.close();

  // calculate spacing of grid points
  dx = rdim/(double)(nw - 1);
  dz = zdim/(double)(nh - 1);

  // set up array of flux values
  psi = new double[nw];
  for(int i=0; i<nw; i++) 
    psi[i] = (sibry - simag)*(double)i/(double)(nw - 1) + simag;

  return FIO_SUCCESS;
}

int geqdsk_source::close()
{
  if(psirz) {
    for(int i=0; i<nh; i++)
      delete[] psirz[i];
    delete[] psirz;
  }
  if(fpol) delete[] fpol;
  if(ffprime) delete[] ffprime;
  if(press) delete[] press;
  if(psi) delete[] psi;

  psirz = 0;
  fpol = 0;
  ffprime = 0;
  press = 0;
  psi = 0;

  return FIO_SUCCESS;
}

int geqdsk_source::center(double *r0, double *z0) const
{
  *r0 = rmaxis;
  *z0 = zmaxis;
  return FIO_SUCCESS;
}

int geqdsk_source::extent(double* r0, double* r1, double* z0, double* z1)
  const
{
  *r0 = rleft;
  *r1 = rleft + dx*(double)(nw-1);
  *z0 = zbottom;
  *z1 = zbottom + dz*(double)(nh-1);

  return FIO_SUCCESS;
}

int geqdsk_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  return FIO_SUCCESS;
}

int geqdsk_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_CURRENT_DENSITY);
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_POLOIDAL_FLUX);
  fields->push_back(FIO_POLOIDAL_FLUX_NORM);
  fields->push_back(FIO_TOTAL_PRESSURE);

  return FIO_SUCCESS;
}

int geqdsk_source::get_field(const field_type t,fio_field** f,
			     const fio_option_list* opt)
{
  switch(t) {
  case(FIO_CURRENT_DENSITY):
    *f = new geqdsk_current_density(this);
    break;

  case(FIO_MAGNETIC_FIELD):
    *f = new geqdsk_magnetic_field(this);
    break;

  case(FIO_POLOIDAL_FLUX):
    *f = new geqdsk_psi_field(this, 2.*M_PI);
    break;

  case(FIO_POLOIDAL_FLUX_NORM):
    *f = new geqdsk_psi_field(this, 1./(sibry-simag), simag);
    break;

  case(FIO_TOTAL_PRESSURE):
    *f = new geqdsk_pressure_field(this);
    break;

  default:
    *f = 0;
    return FIO_UNSUPPORTED;
  };
 
  return FIO_SUCCESS;
}

int geqdsk_source::get_series(const series_type t,fio_series** s)
{
  switch(t) {
  case(FIO_MAGAXIS_PSI):
    *s = new fio_scalar_series(simag);
    break;

  case(FIO_LCFS_PSI):
    *s = new fio_scalar_series(sibry);
    break;

  case(FIO_MAGAXIS_R):
    *s = new fio_scalar_series(rmaxis);
    break;

  case(FIO_MAGAXIS_Z):
    *s = new fio_scalar_series(zmaxis);
    break;

  default:
    *s = 0;
    return FIO_UNSUPPORTED;
  };

  return FIO_SUCCESS;
}


int geqdsk_source::interpolate_psi(const double r0, const double z0,
				   double* si) const
{
  double a[4][4];

  double p = (r0-rleft)/dx;
  double q = (z0-zbottom)/dz;
  int i = (int)p;
  int j = (int)q;

  bool result = 
    bicubic_interpolation_coeffs((const double**)psirz,nh,nw,q,p,a);

  if(!result)
    return FIO_OUT_OF_BOUNDS;

  for(int n=0; n<6; n++)
    si[n] = 0.;

  double temp;
  for(int n=0; n<4; n++) {
    for(int m=0; m<4; m++) {
      si[0] += a[m][n]      *pow(p-i,n  )*pow(q-j,m  );
      
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
    }
  }

  return FIO_SUCCESS;
}
