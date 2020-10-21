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
