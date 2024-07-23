#include <iostream>
#include <fstream>
#include <iomanip>
#include <deque>
#include <math.h>
#include <string.h>

#include "fusion_io_field.h"
#include "interpolate.h"

// find x such that f(x) = val
// dim = 1:
// such that x on the line passing through the initial guess for x
// in direction norm
// dim = 2:
// such that x is in the plane passing through the initial guess for x
// with normal norm
int fio_find_val(fio_field* f, const double val, double* x,
		 const double tol, const double max_step, 
		 const int dim, const double* norm, fio_hint h=0)
{
  const double toroidal_extent = 2.*M_PI; // TODO: generalize this

  if(f->dimension() != 1) {
    std::cerr << "Error in fio_find_val: field is not a scalar field"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  const int max_its = 100;
  int result;
  double last_x[3], dx[3], dv[3], v, dl, mod_dv, n[3];
  double x0[3];
  x0[0] = x[0];
  x0[1] = x[1];
  x0[2] = x[2];

  // calculate unit vector from norm
  double d = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
  n[0] = norm[0] / d;
  n[1] = norm[1] / d;
  n[2] = norm[2] / d;

  double s[3];
  double error, last_error, df;

  for(int i=0; i<max_its; i++) {
    // Evaluate the field
    result = f->eval(x, &v, h);

    if(result == FIO_OUT_OF_BOUNDS) {
      if(i==0) {
	std::cerr << "Error in fio_find_val: initial guess is out of bounds"
		  << std::endl;
	std::cerr << "( " << x[0] << ", " << x[1] << ", " << x[2]
		  <<  " )" << std::endl;
	return FIO_OUT_OF_BOUNDS;
      } else {
	// If we've moved out of bounds, cut step size in half and try again

	std::cerr << "Out of bounds.  re-stepping." << std::endl;

	dx[0] /= 2.;
	dx[1] /= 2.;
	dx[2] /= 2.;
	x[0] = last_x[0] + dx[0];
	x[1] = last_x[1] + dx[1];
	x[2] = last_x[2] + dx[2];
	continue;
      }
    } else if(result != FIO_SUCCESS) {
      return result;
    }

    // Check if we are within tolerance
    error = v - val;
    if(fabs(error) < tol) {
      //      std::cerr << "Found." << std::endl;
      return FIO_SUCCESS;
    }

    if(i>10 && last_error < error) {
      if(error*last_error < 0) {
	// spiralling; split the difference
	/*
	std::cerr << "Spiralling.  re-stepping. " << i << std::endl;
	std::cerr << dx[0] << ", " << dx[2] << ", " << dl << std::endl;
	std::cerr << x[0] << ", " << x[2] << std::endl;
	*/
	dx[0] /= 2.;
	dx[1] /= 2.;
	dx[2] /= 2.;
	x[0] = last_x[0] + dx[0];
	x[1] = last_x[1] + dx[1];
	x[2] = last_x[2] + dx[2];
	continue;
      } else {
	// moving in wrong direction; back up
	//	std::cerr << "Moving in wrong direction.  re-stepping. " << i << std::endl;
	dx[0] /= 2.;
	dx[1] /= 2.;
	dx[2] /= 2.;
	x[0] = last_x[0] + dx[0];
	x[1] = last_x[1] + dx[1];
	x[2] = last_x[2] + dx[2];
	continue;
      }
    } else {
      last_error = error;
    }

    // Evaluate the derivative of the field and take a step
    result = f->eval_deriv(x, dv, h);
    if(result != FIO_SUCCESS) return result;

    // define the direction to move
    if(dim==1) {
      // move in direction supplied by user
      s[0] = n[0];
      s[1] = n[1];
      s[2] = n[2];
    } else {
      // move in direction of steepest gradient but normal to n
      s[0] = dv[0]*(1.-n[0]);
      s[1] = dv[1]*(1.-n[1]);
      s[2] = dv[2]*(1.-n[2]);

      double sn = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
      s[0] /= sn;
      s[1] /= sn;
      s[2] /= sn;
    }

    // derivative in the direction to move
    df = dv[0]*s[0] + dv[1]*s[1] + dv[2]*s[2];

    if(fabs(df) < tol) {
      std::cerr << "Error in fio_find_val: zero gradient" << std::endl;
      return FIO_DIVERGED;
    }

    double dl = -error / df;
    if(dl > max_step) {
      dl = max_step;
    } else if(dl < -max_step) {
      dl = -max_step;
    }

    dx[0] = dl*s[0];
    dx[1] = dl*s[1];
    dx[2] = dl*s[2];
    /*
    std::cerr << "dv: " << dv[0] << ", " << dv[1] << ", " << dv[2] << std::endl;
    std::cerr << "dx: " << dx[0] << ", " << dx[1] << ", " << dx[2] << std::endl;
    */
    last_x[0] = x[0];
    last_x[1] = x[1];
    last_x[2] = x[2];
    x[0] += dx[0];
    x[1] += dx[1];
    x[2] += dx[2];

    if(x[1] < 0.) x[1] += toroidal_extent;
    if(x[1] >= toroidal_extent) x[1] -= toroidal_extent;
  }

  /*
  if(axis) {
    // If Newton iterations don't converge, try brute force search
    result = fio_find_val_1d_brute(f, val, x, tol, max_step, axis, h);
    if(result==FIO_SUCCESS)
      return FIO_SUCCESS;
  }
  */
  /*
  std::cerr << "Error: fio_find_val did not converge after " << max_its
	    << " iterations" << std::endl;
  */
  /*
  std::cerr << " df = " << df << std::endl;
  std::cerr << " dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
  std::cerr << " s = " << s[0] << " " << s[1] << " " << s[2] << std::endl;
  */
  /*
  for(int i=0; i<max_its; i++) {
    std::cerr << "f(" << x_hist[0][i] << "," << x_hist[1][i] << "," << x_hist[2][i]
	      << ")_ " << h_hist[i] << " = " << val_hist[i] << std::endl;
  }
  std::cerr << "looking for " << val << std::endl;

  x[0] = x0[0];
  x[1] = x0[1];
  x[2] = x0[2];
  */
  return FIO_DIVERGED;


  //  return FIO_SUCCESS;
}


// find local maximum of f
// dim = 1: only search on line passing through x in direction norm
// dim = 2: only search in plane passing through x with normal norm
// dim = 3: fully 3D search
int fio_find_max(fio_field* f, double* val, double* x,
		 const double tol, const double max_step, 
		 const int dim, const double* norm, fio_hint h=0)
{
  const double toroidal_extent = 2.*M_PI; // TODO: generalize this

  if(f->dimension() != 1) {
    std::cerr << "Error in fio_find_val: field is not a scalar field"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  const int max_its = 100;
  int result;
  double last_x[3], dx[3], dv[3], v, dl, mod_dv, n[3];
  double x0[3];
  x0[0] = x[0];
  x0[1] = x[1];
  x0[2] = x[2];

  // calculate unit vector from norm
  double d = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
  n[0] = norm[0] / d;
  n[1] = norm[1] / d;
  n[2] = norm[2] / d;


  double s[3];
  double error, last_error;

  for(int i=0; i<max_its; i++) {
    // Evaluate the field
    result = f->eval_deriv(x, dv, h);

    if(result == FIO_OUT_OF_BOUNDS) {
      if(i==0) {
	std::cerr << "Error in fio_find_max: initial guess is out of bounds"
		  << std::endl;
	std::cerr << "( " << x[0] << ", " << x[1] << ", " << x[2]
		  <<  " )" << std::endl;
	return FIO_OUT_OF_BOUNDS;
      } else {
	// If we've moved out of bounds, cut step size in half and try again

	std::cerr << "fio_find_max: Out of bounds; re-stepping." << std::endl;
	/*
	std::cerr << dx[0] << ", " << dx[2] << ", " << dl << std::endl;
	std::cerr << x[0] << ", " << x[2] << std::endl;
	*/
	dx[0] /= 2.;
	dx[1] /= 2.;
	dx[2] /= 2.;
	x[0] = last_x[0] + dx[0];
	x[1] = last_x[1] + dx[1];
	x[2] = last_x[2] + dx[2];
	continue;
      }
    } else if(result != FIO_SUCCESS) {
      return result;
    }

    // restrict gradient to search direction / plane
    if(dim==1) {
      dv[0] = dv[0]*n[0];
      dv[1] = dv[1]*n[1];
      dv[2] = dv[2]*n[2];
    } else if(dim==2) {
      dv[0] = dv[0]*(1. - n[0]);
      dv[1] = dv[1]*(1. - n[1]);
      dv[2] = dv[2]*(1. - n[2]);
    }
    double df = sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);

    result = f->eval(x, val, h);
    /*
    std::cerr << "Te( " << x[0] << ", " << x[1] << ", " << x[2] << " ) = "
	      << *val << std::endl;
    std::cerr << " grad(Te) = ( "
	      << dv[0] << ", " << dv[1] << ", " << dv[2] << " )"
	      << std::endl;
    */
    if(fabs(df) < tol) {
      /*
      std::cerr << "Found maximum at ("
		<< x[0] << ", " << x[1] << ", " << x[2]
		<< " with value " << *val << std::endl;
      */
      return FIO_SUCCESS;
    }

    // find the derivative of |grad(f)| in the direction of steepest descent
    double x1[3], dv1[3], df1;
    const double del = 1e-5;
    x1[0] = x[0] + dv[0]*del/df;
    x1[1] = x[1] + dv[1]*del/df;
    x1[2] = x[2] + dv[2]*del/df;

    result = f->eval_deriv(x1, dv1, h);
    if(dim==1) {
      dv1[0] = dv1[0]*n[0];
      dv1[1] = dv1[1]*n[1];
      dv1[2] = dv1[2]*n[2];
    } else if(dim==2) {
      dv1[0] = dv1[0]*(1. - n[0]);
      dv1[1] = dv1[1]*(1. - n[1]);
      dv1[2] = dv1[2]*(1. - n[2]);
    }
    df1 = sqrt(dv1[0]*dv1[0] + dv1[1]*dv1[1] + dv1[2]*dv1[2]);

    double ddf = (df1 - df) / del;

    double dl = -df / ddf;

    // take a newton step in the direction of steepest descent
    dx[0] = dv[0]*dl/df;
    dx[1] = dv[1]*dl/df;
    dx[2] = dv[2]*dl/df;

    if(dl > max_step) {
      dx[0] *= max_step/dl;
      dx[1] *= max_step/dl;
      dx[2] *= max_step/dl;
    }
    /*
    std::cerr << " x1 = ( " << x1[0] << ", " << x1[1] << ", " << x1[2] << " )"
	      << std::endl;
    std::cerr << " df = " << df << "; ddf = " << ddf << "; dl = " << dl << std::endl;
    std::cerr << " dx = ( " << dx[0] << ", " << dx[1] << ", " << dx[2] << " )"
	      << std::endl;
    */
    last_x[0] = x[0];
    last_x[1] = x[1];
    last_x[2] = x[2];
    x[0] += dx[0];
    x[1] += dx[1];
    x[2] += dx[2];

    if(x[1] < 0.) x[1] += toroidal_extent;
    if(x[1] >= toroidal_extent) x[1] -= toroidal_extent;
  }

  std::cerr << "Error: fio_find_max did not converge after " << max_its
	    << " iterations" << std::endl;

  /*
  x[0] = x0[0];
  x[1] = x0[1];
  x[2] = x0[2];
  */
  return FIO_DIVERGED;
}



int fio_isosurface_2d(fio_field* f, const double val, const double* guess,
		      const double* axis, const double* norm, const double* dd, 
		      const double dl, const double tol, const double max_step,
		      int* n, double*** path, fio_hint h=0)
{
  const int max_pts = 1000;

  int nn;
  double x[3];
  x[0] = guess[0];
  x[1] = guess[1];
  x[2] = guess[2];

  std::deque<double> p[3];
  p[0].clear();
  p[1].clear();
  p[2].clear();

  *n = 0;
  nn = 0;

  while(nn < max_pts) {
    int result;
    double dv[3], dx[3];
    double t[3], dt, de2;

    result = fio_find_val(f, val, x, tol, max_step, 2, norm, h);

    if(result == FIO_OUT_OF_BOUNDS) {
      std::cerr << "fio_find_val failed at nn = " << nn << "\n"
		<< " x = " << x[0] << " " << x[1] << " " << x[2] << "\n"
		<< " norm = " << norm[0] << " " << norm[1] << " " << norm[2] << "\n"
		<< std::endl;
      std::cout << "\n";
      for(int j=0; j<nn; j++) {
	std::cout << p[0][j] << " " << p[1][j] << " " << p[2][j] << std::endl;
      }
      return result;
    }
    
    p[0].push_back(x[0]);
    p[1].push_back(x[1]);
    p[2].push_back(x[2]);
    nn++;

    //    std::cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << std::endl;

    // t = norm x grad(f)
    if(dd) {
      t[0] = dd[0];
      t[1] = dd[1];
      t[2] = dd[2];
    } else {
      result = f->eval_deriv(x, dv, h);
      //      dv[1] /= x[0];  // Convert from partial derivs to gradient
      if(result != FIO_SUCCESS) return result;

      t[0] = (norm[1]*dv[2]-norm[2]*dv[1]);
      t[1] = (norm[2]*dv[0]-norm[0]*dv[2]);
      t[2] = (norm[0]*dv[1]-norm[1]*dv[0]);
    }
    dt = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);

    if(dt==0.) {
      std::cerr << "Error in fio_isosurface: no gradient" << std::endl;
      return FIO_DIVERGED;
    }

    dx[0] = x[0]*cos(x[1]) - p[0][0]*cos(p[1][0]);
    dx[1] = x[0]*sin(x[1]) - p[0][0]*sin(p[1][0]);
    dx[2] = x[2] - p[2][0];
    de2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    
    bool looped;
    looped = ((de2 < 4.*dl*dl) && nn > 3);

    if(looped) {
      /*
      std::cerr << "Start: "
		<< p[0][0] << ", " << p[1][0] << ", " << p[2][0]
		<< std::endl;
      std::cerr << "t: "
		<< t[0] << ", " << t[1] << ", " << t[2]
		<< std::endl;
      std::cerr << "End: "
		<< x[0] << ", " << x[1] << ", " << x[2]
		<< std::endl;
      std::cerr << "Dist: "
		<< dx[0] << ", " << dx[1] << ", " << dx[2]
		<< std::endl;
      std::cerr << "dt = " << dt << std::endl;
      */
      (*path) = new double*[3];
      (*path)[0] = new double[nn];
      (*path)[1] = new double[nn];
      (*path)[2] = new double[nn];
      for(int j=0; j<nn; j++) {
	(*path)[0][j] = p[0][j];
	(*path)[1][j] = p[1][j];
	(*path)[2][j] = p[2][j];
      }
      *n = nn;
	
      return FIO_SUCCESS;
    }
    /*
    std::cerr << "dl = " << dl << std::endl;
    std::cerr << "x[0] = " << x[0] << std::endl;
    std::cerr << "t[1] / dt = " << t[1] / dt << std::endl;
    std::cerr << "dphi = " << dl * t[1] / dt / x[0] << std::endl;
    */
    x[0] += dl * t[0] / dt;
    x[1] += dl * t[1] / dt;
    x[2] += dl * t[2] / dt;
  }

  std::cerr << "Error in fio_isosurface: surface not closed after max_pts"
	    << std::endl;

  std::ofstream dump_file;
  std::string filename = "dump_" + std::to_string(guess[1])
    + "_" + std::to_string(guess[0]) + "_" + std::to_string(guess[2]);
  std::cerr << "Dumping to " << filename << std::endl;

  dump_file.open(filename, std::ofstream::out|std::ofstream::trunc);
  for(int i=0; i<nn; i++) {
    dump_file << std::setprecision(10) << std::setw(16) << p[1][i]
	      << std::setprecision(10) << std::setw(16) << p[0][i]
	      << std::setprecision(10) << std::setw(16) << p[2][i]
	      << std::endl;
  }

  dump_file.close();

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
		   const double max_step, const int nphi, 
		   int* n, double*** path,
		   fio_hint h=0)
{
  const int max_pts = 1000000;
  const double toroidal_extent = 2.*M_PI; // TODO: generalize this

  double pplane[3];
  pplane[0] = 0.;
  pplane[1] = 1.;
  pplane[2] = 0.;
  double x[3];
  x[0] = guess[0];
  x[1] = guess[1];
  x[2] = guess[2];

  double dx = x[0] - enclose[0];
  double dz = x[2] - enclose[2];
  double d = sqrt(dx*dx+dz*dz);
  double co = dx/d;
  double sn = dz/d;

  
  std::deque<double> p[3];
  p[0].clear();
  p[1].clear();
  p[2].clear();

  double offset = 0.;
  double dtheta = 0.;
  *n = 0;

  int i_start = 0;
  int iphi = 0;

  while(*n < max_pts) {
    int result;
    double t[3], e[3], dv[3];
    double dt, de2;

    if(*n==i_start) {
      // For first point, keep angle fixed
      result = fio_find_val(f, val, x, tol, max_step, 1, enclose, h);
    } else {
      result = fio_find_val(f, val, x, tol, max_step, 2, pplane, h);
    }
    if(result != FIO_SUCCESS) return result;
    /*    
    if(*n >= 1) {
      dx = x[0] - p[0][*n-1];
      dz = x[2] - p[2][*n-1];
      if(dx*dx + dz*dz > 16.*dl*dl) {
	std::cerr << "Warning: big jump in position ("
		  << p[0][*n-1] << ", " << p[2][*n-1] << ") -> ("
		  << x[0] << ", " << x[2] << std::endl;
      }
    }
    */
    
    p[0].push_back(x[0]);
    p[1].push_back(x[1]);
    p[2].push_back(x[2]);
    (*n)++;

    //    std::cout << x[0] << ", " << x[2] << ", " << dtheta << std::endl;
    
    result = f->eval_deriv(x, dv, h);
    if(result != FIO_SUCCESS) return result;
    
    // t = \hat{phi} x grad(f)
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

    dx = x[0] - p[0][i_start];
    dz = x[2] - p[2][i_start];
    de2 = dx*dx + dz*dz;
    
    bool looped;
    /*
    looped = (fabs(dtheta) >= 2.*M_PI)
      || ((de2 < dl*dl) && (*n - i_start) > 3);
    */
    looped =  ((de2 < dl*dl) && (*n - i_start) > 3);

    if(looped) {
      /*
      e[0] = x[0]-p[0][i_start];
      e[2] = x[2]-p[2][i_start];
      de2 = e[0]*e[0] + e[2]*e[2];
      if(de2 > dl*dl) {
	std::cerr << "Warning: path does not close" << std::endl;
	std::cerr << p[0][i_start] << ", " << p[2][i_start] << std::endl;
	std::cerr << p[0][*n] << ", " << p[2][*n] << std::endl;
	std::cerr << x[0] << ", " << x[2] << std::endl;
	std::cerr << dtheta / (2.*M_PI) << std::endl;
      }
      */

      // close path
      /*
      x[0] = p[0][i_start];
      x[1] = p[1][i_start];
      x[2] = p[2][i_start];
      p[0].push_back(x[0]);
      p[1].push_back(x[1]);
      p[2].push_back(x[2]);
      (*n)++;
      */

      // if we've gone around toroidally, convert data to arrays and exit
      if(iphi==nphi-1) {      
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

      // if we've gone around once poloidally, take a step toroidally
      // new guess should maintain poloidal angle but change distance to 
      // stay on isosurface
      //      0 = dl*(dv/dx * co + dv/dz * sn) + dv/dphi * dphi
      dx = p[0][i_start] - enclose[0];
      dz = p[2][i_start] - enclose[2];
      d = sqrt(dx*dx + dz*dz);
      double dphi = toroidal_extent/nphi;
      double dl = -dv[1]*dphi / (dv[0]*co + dv[2]*sn);
      x[0] = enclose[0] + (d + dl)*co;
      x[1] = guess[1] + toroidal_extent * (double)(iphi+1) / (double)nphi;
      x[2] = enclose[2] + (d + dl)*sn;
      if(x[1] >= toroidal_extent) {
	x[1] -= toroidal_extent;
	offset += toroidal_extent;
      }
      if(x[1] < 0.) {
	x[1] += toroidal_extent;
	offset -= toroidal_extent;
      }
      dtheta = 0.;
      i_start = *n+1;
      iphi++;

    }

    // take a step
    x[0] += t[0];
    x[2] += t[2];
  }

  std::cerr << "Error in fio_isosurface: surface not closed after max_pts"
	    << std::endl;

  return FIO_DIVERGED;
}

// fio_geom_isosurface
// ~~~~~~~~~~~~~~~~~~~
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

int fio_geom_isosurface(fio_field* f, const double val, const double* guess,
			const double* axis, const double tol,
			const double max_step,
			const int nphi, const int ntheta, 
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


//      std::cerr << "theta, l = " << theta << ", " << l << std::endl;
//      std::cerr << "Guess: (" << x[0] << ", " << x[1] << ", " << x[2] << ")"
//		<< std::endl;
      result = fio_find_val(f, val, x, tol, max_step, 1, axis, h);
      if(result != FIO_SUCCESS)
	return result;

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


// fio_gridify_loop
// ~~~~~~~~~~~~~~~~
// interpolates a path in (R, Z) onto a points regularly spaced in the
// geometric angle.
//
// input:
//   m0: number of points in original path 
//   path0[3][m0]:  arrays containing R, phi, Z coordinates of original path
//   axis[3]:       R, phi, Z coordinates of axis for calculating angle
//   m: number of points for regularly-spaced (gridded) path
//   parameter:   1 = poloidal angle, 2 = toroidal angle, 
//                other = normlized arclength
// output: 
//   path[3][m]:    arrays containing R, phi, Z coordinates of gridded path
//   theta[m]:      array containing theta coordinates of gridded path
//
// Note: "path" and "theta" must be allocated before calling this function
// Note: this function may fail if path0 is not single-valued in theta

int fio_gridify_loop(const int m0, double** path0, const double* axis, 
		     const int m, double** path, double* theta, 
		     const int param=0)
{ 
  int dir = 0;
  double* l0 = new double[m0+1];

  for(int i=0; i<=m0; i++) {
    int ii = (i < m0) ? i : 0;
    if(param==1) {
      // parameterize by poloidal angle
      double dx[3];
      dx[0] = path0[0][ii] - axis[0];
      dx[2] = path0[2][ii] - axis[2];
      l0[i] = atan2(dx[2],dx[0]);
      if(l0[i] < 0.) l0[i] += 2.*M_PI;
    } else if(param==2) {
      // parameterize by toroidal angle
      l0[i] = path0[1][ii];
      /*
      if(l0[i] < 0.) {
	std::cerr << "error: toroidal angle < 0" << std::endl;
	std::cerr << "i = " << i << std::endl;
	std::cerr << "l0[i] = " << l0[i] << std::endl;
      }
      */
    } else {

      // parameterize by arclength
      if(i==0)
	l0[i] = 0.;
      else {
	double dx[3];
	dx[0] = path0[0][ii]*cos(path0[1][ii]) 
	  - path0[0][i-1]*cos(path0[1][i-1]);
	dx[1] = path0[0][ii]*sin(path0[1][ii])
	  - path0[0][i-1]*sin(path0[1][i-1]);
	dx[2] = path0[2][ii]-path0[2][i-1];
	
	double dl = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
	l0[i] = l0[i-1] + dl;
      }
    }

    if(i > 0) {
      double d = l0[i] - l0[i-1];
      
      if(i==1)
	dir = (d < 0) ? -1 : 1;
      else if (d*dir < 0) {
	std::cerr << "Error in fio_gridify_surface_2d: d < 0\n"
		  << i << " " << l0[i] << " " << l0[i-1] << " " 
		  << path0[1][ii] << std::endl;
	return FIO_DIVERGED;
      }
    }
  }

  //  std::cerr << "l0[m0-1] = " << l0[m0-1] << std::endl;

  // interpolate path onto regular arc length grid with m points
  path[0][0] = path0[0][0];
  path[1][0] = path0[1][0];
  path[2][0] = path0[2][0];
  theta[0] = 0.;
  int j=0;

  for(int i=1; i<m; i++) {
    double l;
    if(param==1 || param==2) {
      theta[i] = (double)i/(double)(m) * 2.*M_PI;
      l = theta[i];
    } else {
      theta[i] = (double)i/(double)(m);
      l = theta[i]*l0[m0];
    }
    /*
    path[0][i] = path0[0][i];
    path[1][i] = path0[1][i];
    path[2][i] = path0[2][i];
    */

    while(l0[j+1] < l) {
      j++;
      if(j==m0) {
	std::cerr << " j not found! " << std::endl;
	std::cerr << "l = " << l << std::endl;
	std::cerr << "l0[m0] = " << l0[m0] << std::endl;
	/*
	std::cerr << "gridify_loop called with\n";
	std::cerr << "m0 = " << m0 << '\n';
	std::cerr << "m = " << m << '\n';
	std::cerr << "param = " << param << std::endl;
	*/
	j = m0-1;
	break;
      } 
    }
    double delta_l = l0[j+1] - l0[j];
    double dl = l - l0[j];
    int jp = (j+1<m0) ? j+1 : 0;
    path[0][i] = path0[0][j]*(1.-dl/delta_l) + path0[0][jp]*(dl/delta_l);
    path[1][i] = path0[1][j]*(1.-dl/delta_l) + path0[1][jp]*(dl/delta_l);
    path[2][i] = path0[2][j]*(1.-dl/delta_l) + path0[2][jp]*(dl/delta_l);

    /*
    bool r = true;
    r = cubic_interpolation(m0, l0, l, path0[0], &(path[0][i])) && r;
    r = cubic_interpolation(m0, l0, l, path0[1], &(path[1][i])) && r;
    r = cubic_interpolation(m0, l0, l, path0[2], &(path[2][i])) && r;

    if(!r)
      std::cerr << "Error interpolating " << i << std::endl;
    */
  }

  delete[] l0;

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
      
      result = fio_gridify_loop(m, tmp_path0, axis, 
				ntheta, tmp_path, theta, 0);
      if(result != FIO_SUCCESS)
	return result;
      
      istart = i;
      iphi++;
    }

  }
  if(iphi != nphi) {
    std::cerr 
      << "Error in fio_gridify_surface: wrong number of toroidal points " 
      << iphi << ", " << nphi << std::endl;
  }

  return FIO_SUCCESS;
}


// fio_gridded_isosurface
// ~~~~~~~~~~~~~~~~~~~~~~
// finds an isosurface in field f and returns path on uniform theta, phi grid
// where theta is the fractional distance around the surface
//
// inputs: 
//    val: value of isosurface (i.e. find surface where f = val)
//    guess: initial guess for first point on isosurface
//    axis[3]:  R, phi, Z coordinates of poloidal axis
//    dl: initial approximate spacing between points when finding surface
//    tol: max acceptable different between f and val on isosurface
//    nphi: number of toroidal points
//    ntheta: number ot poloidal points
//
// outputs:
//    path[3][ntheta*nphi]: R, phi, Z coordinates of points on path
//    phi[nphi]: values of phi
//    theta
//    h: hint for last point on path
//
int fio_gridded_isosurface(fio_field* f, const double val, const double* guess,
			   double** axis, const double dl_tor, const double dl_pol,
			   const double tol, const double max_step,
			   const int nphi, const int ntheta, 
			   double* phi, double* theta,
			   double** path, const char* label, fio_hint h=0)
{
  int result;
  int n;

  std::ofstream file;
  std::string filename;

  double norm[3];
  norm[0] = 0.;
  norm[1] = 1.;
  norm[2] = 0.;

  if(label) {
    filename = "pol_loop_raw_" + std::string(label) + ".dat";
    file.open(filename, std::ofstream::out|std::ofstream::trunc);
  }

  double t[3];
  double x[3];
  x[0] = guess[0];

  double a[3];
  a[0] = axis[0][0];
  a[1] = axis[1][0];
  a[2] = axis[2][0];

  // Loop over each toroidal point and calculate poloidal path
  for(int i=0; i<nphi; i++) {
    double** path_pol0;

    phi[i] = 2.*M_PI*i/nphi;

    // find axis in phi = phi[i] plane
    double te_max;
    t[0] = 0.;
    t[1] = 1.;
    t[2] = 0.;
    a[1] = phi[i];
    
    result = fio_find_max(f, &te_max, a, tol, max_step, 2, t, h);
    if(result != FIO_SUCCESS) {
      std::cerr << "Error finding magnetic axis at phi = "
		<< phi[i] << std::endl;
      return result;
    }
    axis[0][i] = a[0];
    axis[1][i] = a[1];
    axis[2][i] = a[2];

    // find initial point in horizontal line from axis
    t[0] = 1.;
    t[1] = 0.;
    t[2] = 0.;
    x[1] = phi[i];
    x[2] = a[2];
    result = fio_find_val(f, val, x, tol, 0.1, 1, t, h);

    if(result != FIO_SUCCESS) {
      std::cerr << "Error finding initial point" << std::endl;
      break;
    }
    /*
    std::cerr << "Poloidal loop " << i << " at Te = " << val
	      << " phi = " << x[1] << std::endl;
    */

    double ax[3];
    ax[0] = axis[0][i];
    ax[1] = axis[1][i];
    ax[2] = axis[2][i];
    result = fio_isosurface_2d(f, val, x, ax, norm, 0, dl_pol, tol, max_step,
			       &n, &path_pol0, h);

    if(result != FIO_SUCCESS) {
      std::cerr << "Error finding poloidal loop at phi = " << phi[i]
		<< std::endl;
      break;
    }
    x[0] = path_pol0[0][n-1];
    x[1] = path_pol0[1][n-1];
    x[2] = path_pol0[2][n-1];

    if(label) {
      for(int j=0; j<n; j++) {
	file << std::setprecision(10) << std::setw(16) << path_pol0[1][j] 
	     << std::setprecision(10) << std::setw(16) << path_pol0[0][j] 
	     << std::setprecision(10) << std::setw(16) << path_pol0[2][j] 
	     << std::endl;
      }
    }

    double* path_pol[3];
    path_pol[0] = &(path[0][i*ntheta]);
    path_pol[1] = &(path[1][i*ntheta]);
    path_pol[2] = &(path[2][i*ntheta]);

    result = fio_gridify_loop(n, path_pol0, 0, ntheta, path_pol, theta, 0);
    
    if(result != FIO_SUCCESS) {
      std::cerr << "Error gridifying poloidal path " << std::endl;
      break;
    }
  }

  //  std::cerr << "Writing final path data" << std::endl;
  if(label) {
    file.close();

    filename = "surface_" + std::string(label) + ".dat";
    file.open(filename, std::ofstream::out|std::ofstream::trunc);

    int k=0;
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {
	file << std::setprecision(10) << std::setw(16) << path[1][k] 
	     << std::setprecision(10) << std::setw(16) << path[0][k] 
	     << std::setprecision(10) << std::setw(16) << path[2][k] 
	     << std::endl;
	k++;
      }
      file << std::endl;
    }
    file.close();
  }

  return result;
}


int fio_q_at_plane(fio_field* f, const int n, double** x, double* q,
		   double* bpol, double** bout, fio_hint h=0)
{
  int result;

  *q = 0.;

  if(f->dimension() != 3) {
    std::cerr << "Error in fio_q_at_plane: field is not a vector field"
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
    if(bout) {
      bout[0][i] = b[0];
      bout[1][i] = b[1];
      bout[2][i] = b[2];
    }
  }

  *q /= (2.*M_PI);

  return FIO_SUCCESS;
}

int fio_eval_on_path(fio_field* f, const int n, double** x,
		     double** fout, fio_hint h=0)
{
  int result;
  int dim = f->dimension();

  double* b = new double[dim];
  
  for(int i=0; i<n; i++) {
    double p[3];
    
    p[0] = x[0][i];
    p[1] = x[1][i];
    p[2] = x[2][i];

    result = f->eval(p, b, h);
    if(result != FIO_SUCCESS) 
      return result;

    for(int j=0; j<dim; j++)
      fout[j][i] = b[j];
  }

  delete[] b;

  return FIO_SUCCESS;
}


int fio_surface_average(fio_field* f, int ns, int nphi, int ntheta,
			double** x, double* jac,
			double* fout, fio_hint h=0)
{
  if(f->dimension() > 1) {
    std::cerr << "Error: can only surface-average scalar fields"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  int di = ntheta*nphi;
  int dj = ntheta;
  int dk = 1;

  for(int i=0; i<ns; i++) {
    fout[i] = 0.;
    double dV = 0.;
    for(int j=0; j<nphi; j++) {
      for(int k=0; k<ntheta; k++) {
	int ijk = k*dk + j*dj + i*di;
	
	double p[3];
	p[0] = x[0][ijk];
	p[1] = x[1][ijk];
	p[2] = x[2][ijk];

	double b;
	int result = f->eval(p, &b, h);
	if(result != FIO_SUCCESS) 
	  return result;

	fout[i] += b*jac[ijk];
	dV += jac[ijk];
      }
    }
    fout[i] /= dV;
  }

  return FIO_SUCCESS;
}

int fio_isosurface_jacobian(int ns, int nphi, int ntheta,
			    double** x, double* jac)
{
  int result;

  if(ns==1) {
    std::cerr << "Error: cannot calculate jacobian with one surface"
	      << std::endl;
    return FIO_UNSUPPORTED;
  }

  int di = ntheta*nphi;
  int dj = ntheta;
  int dk = 1;

  double* R = x[0];
  double* Z = x[2];
  double dphidj = 2.*M_PI / nphi;

  int jac_sign = 0;
  int ierr = FIO_SUCCESS;

  for(int i=0; i<ns; i++ ) {
    for(int j=0; j<nphi; j++) {
      for(int k=0; k<ntheta; k++) {
	double dRdi, dRdj, dRdk, dZdi, dZdj, dZdk;
	int ijk = k*dk + j*dj + i*di;

	if(i==0) {
	  dRdi = R[ijk + di] - R[ijk];
	  dZdi = Z[ijk + di] - Z[ijk];
	} else if(i==ns-1) {
	  dRdi = R[ijk] - R[ijk - di];
	  dZdi = Z[ijk] - Z[ijk - di];
	} else {
	  dRdi = (R[ijk+di] - R[ijk-di])/2.;
	  dZdi = (Z[ijk+di] - Z[ijk-di])/2.;
	}
	
	int jm = (j==0 ? nphi-1 : j-1);
	int jp = (j==nphi-1 ? 0 : j+1);
	dRdj = (R[k*dk + jp*dj + i*di] - R[k*dk + jm*dj + i*di])/2.;
	dZdj = (Z[k*dk + jp*dj + i*di] - Z[k*dk + jm*dj + i*di])/2.;

	int km = (k==0 ? ntheta-1 : k-1);
	int kp = (k==ntheta-1 ? 0 : k+1);
	dRdk = (R[kp*dk + j*dj + i*di] - R[km*dk + j*dj + i*di])/2.;
	dZdk = (Z[kp*dk + j*dj + i*di] - Z[km*dk + j*dj + i*di])/2.;
	
	// di . (dk x dj)
	// where i is the radial index
	//       k is the poloidal index
	//       j is the toroidal index
	jac[ijk] = (dRdi*dZdk - dRdk*dZdi)*R[ijk]*dphidj;

	if(jac_sign==0) {
	  jac_sign = (jac[ijk] > 0) ? 1 : -1;
	} else {
	  if(jac[ijk]*jac_sign <= 0.) {
	    //	    std::cerr << "Error: Jacobian flips sign. "
	    //		      << i << ", " << j << ", " << k << std::endl;
	    ierr = 1.;
	  }
	}

      }
    }
  }

  if(ierr==1) {
    std::cerr << "Error: Jacobian flips sign" << std::endl;
  }

  return ierr;
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

int fio_eval_on_surface(fio_field* f, const int n, double** x, double** a,
			fio_hint h=0)
{
  int result, dim;

  dim = f->dimension();

  for(int i=0; i<n; i++) {
    double v[dim], p[3];

    p[0] = x[0][i];
    p[1] = x[1][i];
    p[2] = x[2][i];
    result = f->eval(p, v, h);
    if(result !=FIO_SUCCESS)
      return result;

    for(int j=0; j<dim; j++)
      a[i][j] = v[j];
  }

  return FIO_SUCCESS;
}
