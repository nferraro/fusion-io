#include <iostream>
#include <deque>
#include <math.h>

#include "fusion_io_field.h"
#include "interpolate.h"

int fio_find_val_2d(fio_field* f, const double val, double* x,
		    const double tol, const double max_step, fio_hint h=0)
{
  if(f->dimension() != 1) {
    std::cerr << "Error in fio_find_val_2d: field is not a scalar field"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  const int max_its = 20;
  int result;
  double last_x[3], dx[3], dv[3], v, dl, mod_dv;
  
  for(int i=0; i<max_its; i++) {
    // Evaluate the field
    result = f->eval(x, &v, h);
    /*    
    std::cerr << "f(" << x[0] << ", " << x[1] << ", " << x[2] << ") = "
	      << v << std::endl;
    std::cerr << val - v << std::endl;
    */
    if(result == FIO_OUT_OF_BOUNDS) {
      if(i==0) { 
	std::cerr << "Error in fio_find_val_2d: initial guess is out of bounds"
		  << std::endl;
	return FIO_OUT_OF_BOUNDS;
      } else {
	// If we've moved out of bounds, cut step size in half and try again
	dx[0] /= 2;
	dx[2] /= 2;
	x[0] = last_x[0] + dx[0];
	x[2] = last_x[2] + dx[2];
	continue;
      }
    } else if(result != FIO_SUCCESS) {
      return result;
    }

    // Check if we are within tolerance
    if(fabs(v - val) < tol) {
      return FIO_SUCCESS;
    }

    // Evaluate the derivative of the field and take a step
    result = f->eval_deriv(x, dv, h);
    if(result != FIO_SUCCESS) return result;

    mod_dv = sqrt(dv[0]*dv[0] + dv[2]*dv[2]);
    dl = (val - v) / mod_dv;

    if(dl==0.) {
      std::cerr << "Error in fio_find_val_2d: zero gradient" << std::endl;
      return FIO_DIVERGED;
    }
    if(dl > max_step) {
      dl = max_step;
    } else if(dl < -max_step) {
      dl = -max_step;
    }
    dx[0] = dl*dv[0]/mod_dv;
    dx[2] = dl*dv[2]/mod_dv;

    last_x[0] = x[0];
    last_x[2] = x[2];
    x[0] += dx[0];
    x[2] += dx[2];
  }

  std::cerr << "Error in fio_find_val_2d: did not converge"
	    << std::endl;
  return FIO_DIVERGED;
}

int fio_isosurface_2d(fio_field* f, const double val, const double* guess,
		      const double* enclose, const double dl, const double tol,
		      const double max_step, int* n, double*** path,
		      fio_hint h=0)
{
  const int max_pts = 1000;

  double x[3];
  x[0] = guess[0];
  x[1] = guess[1];
  x[2] = guess[2];
  
  std::deque<double> p[3];
  p[0].clear();
  p[1].clear();
  p[2].clear();

  double dtheta = 0;
  *n = 0;

  for(int i=0; i<max_pts; i++) {
    int result;
    double t[3], e[3], dv[3];
    double dt, de2;

    result = fio_find_val_2d(f, val, x, tol, max_step, h);
    if(result != FIO_SUCCESS) return result;
    
    p[0].push_back(x[0]);
    p[1].push_back(x[1]);
    p[2].push_back(x[2]);
    (*n)++;

    //    std::cout << x[0] << ", " << x[2] << ", " << dtheta << std::endl;
    
    result = f->eval_deriv(x, dv, h);
    if(result != FIO_SUCCESS) return result;
    
    t[0] = dv[2];
    t[2] = -dv[0];
    dt = sqrt(t[0]*t[0] + t[2]*t[2]);

    if(dt==0.) {
      std::cerr << "Error in fio_isosurface_2d: no gradient" << std::endl;
      return FIO_DIVERGED;
    }

    t[0] *= dl/dt;
    t[2] *= dl/dt;

    // calculate the approximate additional angle of rotation around "enclose"
    e[0] = x[0] - enclose[0];
    e[2] = x[2] - enclose[2];
    de2 = e[0]*e[0] + e[2]*e[2];
    dtheta += (t[2]*e[0] - t[0]*e[2]) / de2;
    //    std::cerr << "dtheta = " << dtheta << std::endl;
    
    if(fabs(dtheta) > 2.*M_PI) {
      (*path) = new double*[3];
      (*path)[0] = new double[*n];
      (*path)[1] = new double[*n];
      (*path)[2] = new double[*n];
      for(int j=0; j<*n; j++) {
	(*path)[0][j] = p[0][j];
	(*path)[1][j] = p[1][j];
	(*path)[2][j] = p[2][j];
      }
      return FIO_SUCCESS;
    }

    // take a step
    x[0] += t[0];
    x[2] += t[2];
  }

  std::cerr << "Error in fio_isosurface_2d: surface not closed after max_pts"
	    << std::endl;

  return FIO_DIVERGED;
}

int fio_isosurface(fio_field* f, const double val, const double* guess, 
		   const double* enclose)
{
  return FIO_SUCCESS;
}

int fio_gridify_surface_2d(const int m0, double** path0, const double* axis, 
			   const int m, double** path, double* theta)
{
  // Calculate theta of original points
  double* theta0 = new double[m0];
  double min_theta;
  int min_index = 0;
  for(int i=0; i<m0; i++) {
    theta0[i] = atan2(path0[2][i]-axis[2],path0[0][i]-axis[0]);
    if(theta0[i] < 0.) theta0[i] += 2.*M_PI;
    if(i==0 || theta0[i] < min_theta) {
      min_theta = theta0[i];
      min_index = i;
    }
  }

  // shift arrays so that theta0 is monotonic
  shift_array(m0, theta0, min_index);
  shift_array(m0, path0[0], min_index);
  shift_array(m0, path0[1], min_index);
  shift_array(m0, path0[2], min_index);
  /*
  for(int i=0; i<m0; i++) {
    std::cout << path[0][i] << ", " << path[2][i] << ", " 
	      << theta[i] << std::endl;
  }
  */
  // interpolate path onto regular theta grid
  for(int i=0; i<m; i++) {
    theta[i] = 2.*M_PI*i/m;
    cubic_interpolation(m0, theta0, theta[i], path0[0], &(path[0][i]));
    cubic_interpolation(m0, theta0, theta[i], path0[1], &(path[1][i]));
    cubic_interpolation(m0, theta0, theta[i], path0[2], &(path[2][i]));
  }

  delete[] theta0;

  return FIO_SUCCESS;
}

