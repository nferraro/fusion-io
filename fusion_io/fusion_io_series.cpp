#include "fusion_io.h"
#include "fusion_io_series.h"

int fio_scalar_series::bounds(double* tmin, double* tmax) const 
{
  *tmin = 0.;
  *tmax = 0.;

  return FIO_SUCCESS;
}

int fio_scalar_series::eval(const double t, double* x)
{
  *x = data;
  return FIO_SUCCESS;
}

fio_array_series::fio_array_series(const int n_in, const double* d_in)
{
  n = n_in;
  data = new double[n];
  for(int i=0; i<n; i++)
    data[i] = d_in[i];
}

fio_array_series::~fio_array_series()
{
  if(data) delete[] data;
}

int fio_array_series::eval(const double t, double* x)
{
  int i = (int)t;
  if(i<0 || i>=n) return FIO_OUT_OF_BOUNDS;
  *x = data[i];
  return FIO_SUCCESS;
}

int fio_array_series::bounds(double* tmin, double* tmax) const
{
  *tmin = 0.;
  *tmax = (double)(n-1);
  return FIO_SUCCESS;
}
