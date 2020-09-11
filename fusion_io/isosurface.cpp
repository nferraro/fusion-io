#include <iostream>
#include <deque>
#include <math.h>
#include <string.h>

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

  const int max_its = 50;
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

// fio_isosurface
// ~~~~~~~~~~~~~~
// finds an isosurface in field f
//
// inputs: 
//    val: value of isosurface (i.e. find surface where f = val)
//    guess: initial guess for first point on isosurface
//    enclose[3]:  R, phi, Z coordinates of point that surface encloases
//    dl: approximate distance between points on isosurface
//    tol: max acceptable different between f and val on isosurface
//    max_step: max step size of Newton iteration
//    h: hint for first point on path
//
// outputs:
//    n: number of points on isosurface
//    path[3][n]: R, phi, Z coordinates of points on path
//    h: hint for last point on path
//
// Note: the path will be traced until either max_pts is reached, or
//       the path has accumulated 2*pi rotation around "enclose"
int fio_isosurface(fio_field* f, const double val, const double* guess,
		   const double* enclose, const double dl, const double tol,
		   const double max_step, const double dphi, 
		   int* n, double*** path,
		   fio_hint h=0)
{
  const int max_pts = 100000;
  const double toroidal_extent = 2.*M_PI; // TODO: generalize this

  double x[3];
  x[0] = guess[0];
  x[1] = guess[1];
  x[2] = guess[2];
  
  std::deque<double> p[3];
  p[0].clear();
  p[1].clear();
  p[2].clear();

  double offset = 0.;
  double dtheta = 0.;
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
      std::cerr << "Error in fio_isosurface: no gradient" << std::endl;
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
      
      // if we've gone around once poloidally, take a step toroidally
      x[0] = guess[0];
      x[1] += dphi;
      x[2] = guess[2];
      if(x[1] >= toroidal_extent) {
	x[1] -= toroidal_extent;
	offset += toroidal_extent;
      }
      if(x[1] < 0.) {
	x[1] += toroidal_extent;
	offset -= toroidal_extent;
      }
      dtheta = 0.;

      // if we've gone around toroidally, convert data to arrays and exit
      if(x[1]+offset - guess[1] >= toroidal_extent) {
      
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
    }

    // take a step
    x[0] += t[0];
    x[2] += t[2];
  }

  std::cerr << "Error in fio_isosurface: surface not closed after max_pts"
	    << std::endl;

  return FIO_DIVERGED;
}


// fio_gridify_surface
// ~~~~~~~~~~~~~~~~~~~
// interpolates a path in (R, Z) onto a points regularly spaced in the
// geometric angle.
//
// input:
//   m0: number of points in original path 
//   path0[3][m0]:  arrays containing R, phi, Z coordinates of original path
//   axis[3]:       R, phi, Z coordinates of axis for calculating angle
//   m: number of points for regularly-spaced (gridded) path
// output: 
//   path[3][m]:    arrays containing R, phi, Z coordinates of gridded path
//   theta[m]:      array containing theta coordinates of gridded path
//
// Note: "path" and "theta" must be allocated before calling this function
// Note: this function may fail if path0 is not single-valued in theta

int fio_gridify_surface_2d(const int m0, double** path0, const double* axis, 
			   const int m, double** path, double* theta)
{

  const double toroidal_period = 2.*M_PI;    // TODO: generalize this
  double* theta0 = new double[m0+3];         // extra space for interpolation

  // Calculate theta of original points 
  double min_theta;
  int min_index = 0;
  int n_ascending = 0;
  for(int i=0; i<m0; i++) {
    theta0[i+1] = atan2(path0[2][i]-axis[2],path0[0][i]-axis[0]);
    if(theta0[i+1] < 0.) theta0[i+1] += 2.*M_PI;
    if(i==0 || theta0[i+1] < min_theta) {
      min_theta = theta0[i+1];
      min_index = i;
    }
    if((i>0) && (theta0[i+1] > theta0[i]))
      n_ascending++;
  }

  // if poloidal angle is decreasing, reverse everything
  if(n_ascending <= 2) {
    std::cerr << "Reversing poloidal angle..." << std::endl;
    reverse_array(m0, &(theta0[1]));
    reverse_array(m0, path[0]);
    reverse_array(m0, path[1]);
    reverse_array(m0, path[2]);
    min_index = m0-1-min_index;
  }

  // shift arrays so that theta0 is monotonic
  shift_array(m0, &(theta0[1]), min_index);
  shift_array(m0, path0[0], min_index);
  shift_array(m0, path0[1], min_index);
  shift_array(m0, path0[2], min_index);
  
  // add one point before beginning and two points after end to improve interpolation
  double* tmp_path[3];
  tmp_path[0] = new double[m0+3];
  tmp_path[1] = new double[m0+3];
  tmp_path[2] = new double[m0+3];
  memcpy(&(tmp_path[0][1]), path0[0], sizeof(double)*m0);
  memcpy(&(tmp_path[1][1]), path0[1], sizeof(double)*m0);
  memcpy(&(tmp_path[2][1]), path0[2], sizeof(double)*m0);

  theta0[0   ] = theta0[m0]-toroidal_period;
  theta0[m0+1] = theta0[1 ]+toroidal_period;
  theta0[m0+2] = theta0[2 ]+toroidal_period;
  for(int i=0; i<3; i++) {
    tmp_path[i][0   ] = tmp_path[i][m0];
    tmp_path[i][m0+1] = tmp_path[i][1 ];
    tmp_path[i][m0+2] = tmp_path[i][2 ];
  }

  /*
  std::cerr << "theta = ";
  for(int i=0; i<m0+3; i++) {
    std::cerr << theta0[i] << " ";
  }
  std::cerr << std::endl;
  std::cerr << "tmp_path = ";
  for(int i=0; i<m0+3; i++) {
    std::cerr << tmp_path[0][i] << " ";
  }
  std::cerr << std::endl;
  */

  // interpolate path onto regular theta grid
  for(int i=0; i<m; i++) {
    theta[i] = 2.*M_PI*i/m;
    cubic_interpolation(m0+3, theta0, theta[i], tmp_path[0], &(path[0][i]));
    cubic_interpolation(m0+3, theta0, theta[i], tmp_path[1], &(path[1][i]));
    cubic_interpolation(m0+3, theta0, theta[i], tmp_path[2], &(path[2][i]));
  }

  delete[] theta0;
  delete[] tmp_path[0];
  delete[] tmp_path[1];
  delete[] tmp_path[2];

  return FIO_SUCCESS;
}

int fio_gridify_surface(const int m0, double** path0, const double* axis, 
			const int nphi, const int ntheta, 
			double** path, double* phi, double* theta)
{
  int istart, result, iphi;
  double* tmp_path0[3];
  double* tmp_path[3];

  istart=0;
  iphi = 0;
  for(int i=1; i<m0; i++) {
    int m;

    m = 0;

    if(i==m0-1) {
      m = i-istart+1;
    } else if(path0[1][i] != path0[1][i-1]) {
      m = i-istart;
    }

    if(m > 0) {
      phi[iphi] = path0[1][i-1];
      /*
      std::cerr << "iphi = " << iphi << std::endl;
      std::cerr << "istart = " << istart << std::endl;
      std::cerr << "Found " << m << " points at phi = " << phi[iphi]
		<< std::endl;
      */
      tmp_path0[0] = &(path0[0][istart]);
      tmp_path0[1] = &(path0[1][istart]);
      tmp_path0[2] = &(path0[2][istart]);
      tmp_path[0] = &(path[0][ntheta*iphi]);
      tmp_path[1] = &(path[1][ntheta*iphi]);
      tmp_path[2] = &(path[2][ntheta*iphi]);
      result = fio_gridify_surface_2d(m, tmp_path0, axis, 
				      ntheta, tmp_path, theta);
      if(result != FIO_SUCCESS)
	return result;
      
      istart = i;
      iphi++;
    }

  }
  if(iphi != nphi) {
    std::cerr 
      << "Error in fio_gridify_surface: wrong number of toroidal points" 
      << std::endl;
  }

  return FIO_SUCCESS;
}
