#include "fusion_io.h"

#include <iostream>

fio_field& fio_field::operator+(const fio_field& f)
{
  return *(new fio_field_sum(this, &f));
}

fio_field& fio_field::operator*(const fio_field& f)
{
  return *(new fio_field_product(this, &f));
}


int fio_field::find_val_on_line(const double f0,
			 const double* x0_in, const double *x1,
			 double* x_val, fio_hint hint, const double tol)
{
  const int maxits = 100;
  bool checklen;
  double x0[3];
  double last_l;

  double vec[3];
  if(x0_in[0]==x1[0] && x0_in[1]==x1[1] && x0_in[2]==x1[2]) {
    checklen = false;
    x0[0] = 0.;
    x0[1] = 0.;
    x0[2] = 0.;
  } else {
    checklen = true;
    x0[0] = x0_in[0];
    x0[1] = x0_in[1];
    x0[2] = x0_in[2];
  }

  vec[0] = x1[0] - x0[0];
  vec[1] = x1[1] - x0[1];
  vec[2] = x1[2] - x0[2];

  double len = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

  vec[0] /= len;
  vec[1] /= len;
  vec[2] /= len;

  // Initial guess defined by x_val
  double l =
    (x_val[0] - x0[0])*vec[0] +
    (x_val[1] - x0[1])*vec[1] +
    (x_val[2] - x0[2])*vec[2];

  for(int i=0; i<maxits; i++) {
    if(l < 0  ) l=0.;
    if(l > len && checklen) l=len;

    x_val[0] = l*vec[0] + x0[0];
    x_val[1] = l*vec[1] + x0[1];
    x_val[2] = l*vec[2] + x0[2];

    int result;
    double fval;
    double dfval[3];

    if((result = eval(x_val, &fval, hint)) != FIO_SUCCESS) {
      if(result==FIO_OUT_OF_BOUNDS && i>0) {
	l = (l + last_l)/2.;
	continue;
      }
      std::cerr << "Error in find_val_on_line: "
		<< "Can't evaluate field" << std::endl;
      return result;
    }

    if(fabs(fval - f0) < tol) {
      std::cerr << "Found field = " << fval << " at ( "
		<<  x_val[0] << ", " << x_val[1] << ", " << x_val[2]
		<< ") after " << i << " iterations." << std::endl;
      return FIO_SUCCESS;
    }

    if((result = eval_deriv(x_val, dfval, hint)) != FIO_SUCCESS) {
      std::cerr << "Error in find_val_on_line: "
		<< "Can't evaluate field deriv" << std::endl;
      return result;
    }

    double dfdl = (dfval[0]*vec[0] + dfval[1]*vec[1] + dfval[2]*vec[2])/len;
    last_l = l;
    l -= (fval-f0) / dfdl;
  }

  std::cerr << "find_val_on_line failed to converge after " << maxits
	    << " iterations" << std::endl;
  return FIO_DIVERGED;
}
