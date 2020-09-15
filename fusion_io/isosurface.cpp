#include <iostream>
#include <deque>
#include <math.h>
#include <string.h>

#include "fusion_io_field.h"
#include "interpolate.h"

int fio_find_val_1d_brute(fio_field* f, const double val, double* x,
		    const double tol, const double max_step, 
		    const double* axis, fio_hint h=0)
{
  const int max_its = 100;
  double d0 = x[0] - axis[0];
  double d2 = x[2] - axis[2];
  double d = sqrt(d0*d0 + d2*d2);
  double co = d0/d;
  double sn = d2/d;
  double v;
  double error, last_error;
  double last_d;
  double step = max_step;

  double x0[3];
  x0[0] = x[0];
  x0[1] = x[1];
  x0[2] = x[2];
  
  d = step;

  for(int i=0; i<max_its; i++) {
    x[0] = d*co + axis[0];
    x[2] = d*sn + axis[2];
    
    int result = f->eval(x, &v, h);

    if(result == FIO_OUT_OF_BOUNDS)
      std::cerr << "No sign flip found." << std::endl;

    if(result != FIO_SUCCESS)
      break;

    error = v - val;
    if(fabs(error) < tol)
      return FIO_SUCCESS;

    if(i > 0) {
      if(error*last_error < 0) {
	d = last_d;
	error = last_error;
	step /= 2;
	continue;
      }
    }

    last_d = d;
    last_error = error;
    d += step;
  }

  x[0] = x0[0];
  x[1] = x0[1];
  x[2] = x0[2];

  std::cerr << "Brute force failed" << std::endl;
  return FIO_DIVERGED;
}

int fio_find_val_2d(fio_field* f, const double val, double* x,
		    const double tol, const double max_step, 
		    const double* axis, fio_hint h=0)
{
  if(f->dimension() != 1) {
    std::cerr << "Error in fio_find_val_2d: field is not a scalar field"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  const int max_its = 50;
  int result;
  double last_x[3], dx[3], dv[3], v, dl, mod_dv, co, sn;
  double x0[3];
  x0[0] = x[0];
  x0[1] = x[1];
  x0[2] = x[2];

  if(axis) {
    double d0 = x[0] - axis[0];
    double d2 = x[2] - axis[2];
    double d = sqrt(d0*d0 + d2*d2);
    co = d0/d;
    sn = d2/d;
  }

  double x2[3], v2;
  const double del = 1e-5;
  double error;
  
  for(int i=0; i<max_its; i++) {
    // Evaluate the field
    result = f->eval(x, &v, h);

    if(result == FIO_OUT_OF_BOUNDS) {
      if(i==0) { 
	std::cerr << "Error in fio_find_val_2d: initial guess is out of bounds"
		  << std::endl;
	return FIO_OUT_OF_BOUNDS;
      } else {
	// If we've moved out of bounds, cut step size in half and try again
	/*
	std::cerr << "Out of bounds.  re-stepping." << std::endl;
	std::cerr << dx[0] << ", " << dx[2] << ", " << dl << std::endl;
	std::cerr << x[0] << ", " << x[2] << std::endl;
	*/
	dx[0] /= 2.;
	dx[2] /= 2.;
	x[0] = last_x[0] + dx[0];
	x[2] = last_x[2] + dx[2];
	continue;
      }
    } else if(result != FIO_SUCCESS) {
      return result;
    }

    // Check if we are within tolerance
    error = v - val;
    if(fabs(error) < tol) {
      return FIO_SUCCESS;
    }

    // Evaluate the derivative of the field and take a step
    if(axis) {
      x2[0] = x[0] + del*co;
      x2[1] = x[1];
      x2[2] = x[2] + del*sn;
      result = f->eval(x2, &v2, h);
    } else {
      result = f->eval_deriv(x, dv, h);
    }

    if(result != FIO_SUCCESS) return result;

    if(axis) {
      // if angle is provided, then maintain angle
      //      dl = -error / (dv[0]*co + dv[2]*sn);
      dl = -error / (v2 - v) * del;
    } else {
      // if no angle is provided, then follow path of fastest descent
      mod_dv = sqrt(dv[0]*dv[0] + dv[2]*dv[2]);
      dl = -error / mod_dv;
    }

    if(dl==0.) {
      std::cerr << "Error in fio_find_val_2d: zero gradient" << std::endl;
      return FIO_DIVERGED;
    }
    if(dl > max_step) {
      dl = max_step;
    } else if(dl < -max_step) {
      dl = -max_step;
    }

    if(axis) {
      dx[0] = dl*co;
      dx[2] = dl*sn;
    } else {
      dx[0] = dl*dv[0]/mod_dv;
      dx[2] = dl*dv[2]/mod_dv;
    }

    last_x[0] = x[0];
    last_x[2] = x[2];
    x[0] += dx[0];
    x[2] += dx[2];
  }

  x[0] = x0[0];
  x[1] = x0[1];
  x[2] = x0[2];

  if(axis) {
    // If Newton iterations don't converge, try brute force search
    result = fio_find_val_1d_brute(f, val, x, tol, max_step, axis, h);
    if(result==FIO_SUCCESS)
      return FIO_SUCCESS;
  }

  std::cerr << "Error: fio_find_val_2d did not converge after " << max_its
	    << " iterations" << std::endl;

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

    result = fio_find_val_2d(f, val, x, tol, max_step, 0, h);
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

// fio_gridded_isosurface
// ~~~~~~~~~~~~~~~~~~~~~~
// finds an isosurface in field f and returns path on uniform theta, phi grid
// where theta and phi are the geometric poloidal and toroidal angles.
//
// inputs: 
//    val: value of isosurface (i.e. find surface where f = val)
//    guess: initial guess for first point on isosurface
//    axis[3]:  R, phi, Z coordinates of poloidal axis
//    tol: max acceptable different between f and val on isosurface
//    nphi: number of toroidal points
//    ntheta: number ot poloidal points
//
// outputs:
//    path[3][ntheta*nphi]: R, phi, Z coordinates of points on path
//    h: hint for last point on path
//
int fio_gridded_isosurface(fio_field* f, const double val, const double* guess,
			   const double* axis, const double tol,
			   const double max_step,
			   const double nphi, const double ntheta, 
			   double*** path, fio_hint h=0)
{
  const double toroidal_period = 2.*M_PI; // TODO: generalize this

  double x[3];
  double dx = guess[0] - axis[0];
  double dz = guess[2] - axis[2];
  double l = sqrt(dx*dx + dz*dz);
  double phi, theta;
  int result, k;

  for(int i=0; i<nphi; i++) {
    phi = toroidal_period*i/nphi;
    for(int j=0; j<ntheta; j++) {
      k = i*ntheta + j;
      theta = 2.*M_PI*j/ntheta;
      
      x[0] = l*cos(theta) + axis[0];
      x[1] = phi;
      x[2] = l*sin(theta) + axis[2];

      /*
      std::cerr << "theta, l = " << theta << ", " << l << std::endl;
      std::cerr << "Guess: (" << x[0] << ", " << x[1] << ", " << x[2] << ")"
		<< std::endl;
      */
      result = fio_find_val_2d(f, val, x, tol, max_step, axis, h);
      /*
      if(result != FIO_SUCCESS)
	return result;
      */

      (*path)[0][k] = x[0];
      (*path)[1][k] = x[1];
      (*path)[2][k] = x[2];

      // update guess for distance from axis
      dx = x[0] - axis[0];
      dz = x[2] - axis[2];
      l = sqrt(dx*dx + dz*dz);
    }
  }

  return FIO_SUCCESS;
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


int fio_q_at_surface(fio_field* f, const int n, double** x, double* q,
		     double* bpol, fio_hint h=0)
{
  int result;

  *q = 0.;

  if(f->dimension() != 3) {
    std::cerr << "Error in fio_q_at_surface: field is not a vector field"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  for(int i=0; i<n; i++) {
    double b[3], p[3];
  
    p[0] = x[0][i];
    p[1] = x[1][i];
    p[2] = x[2][i];
    result = f->eval(p, b, h);
    if(result !=FIO_SUCCESS) 
      return result;
    
    int im, ip;
    if(i==0  ) im = n-1; else im = i-1;
    if(i==n-1) ip = 0;   else ip = i+1;

    double dx[3];
    dx[0] = (x[0][ip] - x[0][im]) / 2.;
    dx[1] = (x[1][ip] - x[1][im]) / 2.;
    dx[2] = (x[2][ip] - x[2][im]) / 2.;

    double dl, bp;
    dl = sqrt(dx[0]*dx[0] + 0.*dx[1]*dx[1] + dx[2]*dx[2]);
    if(dl==0.) continue;
    bp = (b[0]*dx[0] + 0.*b[1]*dx[1] + b[2]*dx[2]) / dl;

    //    std::cerr << p[0] << ", " << bp << ", " << b[1] << ", " << dl << std::endl;

    *q += b[1]/(p[0]*bp) * dl;

    if(bpol)
      bpol[i] = bp;
  }

  *q /= (2.*M_PI);

  return FIO_SUCCESS;
}



int fio_surface_average(fio_field* f, const int n, double** x, double* a,
			double* bpol, fio_hint h=0)
{
  int result;

  if(f->dimension() != 1) {
    std::cerr << "Error in fio_q_at_surface: field is not a vector field"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  *a = 0.;
  double denom = 0.;

  for(int i=0; i<n; i++) {
    double v, p[3];
  
    p[0] = x[0][i];
    p[1] = x[1][i];
    p[2] = x[2][i];
    result = f->eval(p, &v, h);
    if(result !=FIO_SUCCESS) 
      return result;

    int im, ip;
    if(i==0  ) im = n-1; else im = i-1;
    if(i==n-1) ip = 0;   else ip = i+1;

    double dx[3];
    dx[0] = (x[0][ip] - x[0][im]) / 2.;
    dx[1] = (x[1][ip] - x[1][im]) / 2.;
    dx[2] = (x[2][ip] - x[2][im]) / 2.;

    double dl;
    dl = sqrt(dx[0]*dx[0] + 0.*dx[1]*dx[1] + dx[2]*dx[2]);
    if(dl==0.) continue;

    *a += v/bpol[i] * dl;
    denom += 1./bpol[i] * dl;
  }
  *a /= denom;

  return FIO_SUCCESS;
}
