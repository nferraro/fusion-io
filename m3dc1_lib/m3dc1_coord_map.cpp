#include "m3dc1_coord_map.h"
#include "output_stream.h"

#include <iostream>

bool m3dc1_mesh_element::is_in_element(const double r, const double phi, const double z,
		     double *xi_frac, double *zi_frac, double *eta_frac) const
{
  const double tol = 1e-6;

  // determine whether phi falls inside element
  if(phi < Phi[0]-tol) return false;
  if(phi > Phi[1]+tol) return false;
  double d = (phi - Phi[0]) / (Phi[1] - Phi[0]);

  // define the triangle in the plane of phi
  double R0 = R[3]*d + R[0]*(1.-d);
  double R1 = R[4]*d + R[1]*(1.-d);
  double R2 = R[5]*d + R[2]*(1.-d);
  double Z0 = Z[3]*d + Z[0]*(1.-d);
  double Z1 = Z[4]*d + Z[1]*(1.-d);
  double Z2 = Z[5]*d + Z[2]*(1.-d);

  if((r - R0)*(Z1 - Z0) - (z - Z0)*(R1 - R0) > tol) return false;
  if((r - R1)*(Z2 - Z1) - (z - Z1)*(R2 - R1) > tol) return false;
  if((r - R2)*(Z0 - Z2) - (z - Z2)*(R0 - R2) > tol) return false;

  // return "fractional" local coordinates

  double h2 = (R1-R0)*(R1-R0) + (Z1-Z0)*(Z1-Z0);
  double co_h = (R1 - R0)/h2;
  double sn_h = (Z1 - Z0)/h2;
  double b_h = ( (R2 - R0)*(R1 - R0) + (Z2 - Z0)*(Z1 - Z0))/h2;
  double c_h = (-(R2 - R0)*(Z1 - Z0) + (Z2 - Z0)*(R1 - R0))/h2;
  *xi_frac =    (r - R0)*co_h + (z - Z0)*sn_h - b_h;
  *eta_frac = (-(r - R0)*sn_h + (z - Z0)*co_h) / c_h;
  *zi_frac = d;

  return true;
}


m3dc1_coord_map::m3dc1_coord_map(m3dc1_3d_mesh* m)
{
  mesh = m;
  elm = new m3dc1_mesh_element[mesh->nelms];
}

m3dc1_coord_map::~m3dc1_coord_map()
{
  delete[] elm;
}

bool m3dc1_coord_map::load(m3dc1_field* Rf, m3dc1_field* Zf)
{  
  R_field = Rf;
  Z_field = Zf;

  // Populate elm array
  for(int i=0; i<mesh->nelms; i++) {
    double xi[6], zi[6], eta[6];

    xi[0] = -mesh->b[i];
    xi[1] =  mesh->a[i];
    xi[2] = 0.;
    xi[3] = -mesh->b[i];
    xi[4] =  mesh->a[i];
    xi[5] = 0.;
    zi[0] = 0.;
    zi[1] = 0.;
    zi[2] = 0.;
    zi[3] = mesh->d[i];
    zi[4] = mesh->d[i];
    zi[5] = mesh->d[i];
    eta[0] = 0.;
    eta[1] = 0.;
    eta[2] = mesh->c[i];
    eta[3] = 0.;
    eta[4] = 0.;
    eta[5] = mesh->c[i];

    for(int j=0; j<6; j++) {
      int e;
      double x, phi, y, r[m3dc1_field::OP_NUM], z[m3dc1_field::OP_NUM];

      mesh->local_to_global(i,xi[j],zi[j],eta[j],&x,&phi,&y);

      e = i;

      if(!R_field->eval(x,phi,y,m3dc1_field::GET_VAL,r,&e)) {
	std::cerr << "Error reading R" << std::endl;
	return false;
      }
      if(!Z_field->eval(x,phi,y,m3dc1_field::GET_VAL,z,&e)) {
	std::cerr << "Error reading Z" << std::endl;
	return false;
      }
      /*
      if(e != i) {
	std::cerr << "Warning: node found in different element" << std::endl;
	std::cerr << e << " " << i << std::endl;
      }
      */
      elm[i].R[j] = r[m3dc1_field::OP_1];
      elm[i].Z[j] = z[m3dc1_field::OP_1];
      if(j==0) {
	elm[i].Phi[0] = phi;
      } else if(j==3) {
	elm[i].Phi[1] = phi;
      }
      if(elm[i].Phi[0] > elm[i].Phi[1])
	elm[i].Phi[1] += mesh->period / mesh->nperiods;
    }
  }
  std::cerr << "Successfully loaded map!" << std::endl;

  return true;
}

bool m3dc1_coord_map::find_element(const double R, const double Phi, const double Z,
			      double *xi_frac, double *zi_frac, double *eta_frac, int* e) const
{
  bool found;

  double p = mesh->period / mesh->nperiods;
  double phi = Phi;
  while(phi <  0.) phi += p;
  while(phi >= p ) phi -= p;

  if(*e >= 0) {
    found = elm[*e].is_in_element(R, phi, Z, xi_frac, zi_frac, eta_frac);
    if(found) return true;
  }

  // Test neighbors
  for(int n=0; n<mesh->nneighbors[*e]; n++) {
    int m = mesh->neighbor[*e][n];
    if(elm[m].is_in_element(R,phi,Z,xi_frac,zi_frac,eta_frac)) {
      *e = m;
      return true;
    }
  }

  // Test neighbors' neighbors
  for(int n=0; n<mesh->nneighbors[*e]; n++) {
    int m = mesh->neighbor[*e][n];
    for(int nn=0; nn<mesh->nneighbors[m]; nn++) {
      int l = mesh->neighbor[m][nn];
      if(elm[l].is_in_element(R,phi,Z,xi_frac,zi_frac,eta_frac)) {
	*e = l;
	return true;
      }
    }
  }

  // Test all
  for(int i=0; i<mesh->nelms; i++) {
    found = elm[i].is_in_element(R, phi, Z, xi_frac, zi_frac, eta_frac);
    if(found) {
      *e = i;
      return true;
    }
  }

  /* std::cerr << "Failed; guess was " << *e << std::endl; */
  output::qerr << "Failed; guess was " << *e << std::endl;

  return false;
}

bool m3dc1_coord_map::find_coordinates(const double R, const double Phi, const double Z,
				       double *x, double* phi, double* y,
				       int* e, int refine) const
{
  double xi_frac, zi_frac, eta_frac;
  double xi, zi, eta;

  // find element containing the R, Phi, Z point
  if(!find_element(R, Phi, Z, &xi_frac, &zi_frac, &eta_frac, e)) {
    /* std::cerr << "failed find_element containing (" */
	      /* << R << ", " << Phi << ", " << Z << ")" << std::endl; */
    output::qerr << "failed find_element containing ("
	      << R << ", " << Phi << ", " << Z << ")" << std::endl;
    return false;
  }

  // ensure that local coordinates are inside triangle
  double h = mesh->a[*e] + mesh->b[*e];
  if(xi_frac >  mesh->a[*e]/h) xi_frac =  mesh->a[*e]/h;
  if(xi_frac < -mesh->b[*e]/h) xi_frac = -mesh->b[*e]/h;
  if(eta_frac > 1.+(h/mesh->b[*e])*xi_frac) eta_frac = 1.+(h/mesh->b[*e])*xi_frac;
  if(eta_frac > 1.-(h/mesh->a[*e])*xi_frac) eta_frac = 1.-(h/mesh->a[*e])*xi_frac;
  if(eta_frac < 0.) eta_frac = 0.;
  if(zi_frac < 0.) zi_frac = 0.;
  if(zi_frac > 1.) zi_frac = 1.;

  // calculate local coordinates in element
  // these coordinates have been calculated using a linear interpolation of
  // the (R, Phi, Z) values at the element nodes and should be refined
  xi = xi_frac*(mesh->a[*e] + mesh->b[*e]);
  zi = zi_frac*mesh->d[*e];
  eta = eta_frac*mesh->c[*e];

  // calculate global logical coordinates
  mesh->local_to_global(*e, xi, zi, eta, x, phi, y);


  // refine the coordinates through Newton iterations
  //  int e0 =  *e;
  double R0[m3dc1_field::OP_NUM], Z0[m3dc1_field::OP_NUM];
  const double min_det = 1.e-10;
  const double tol = 1e-5;
  int op = m3dc1_field::GET_VAL | m3dc1_field::GET_DVAL;
  double x_last = *x;
  double y_last = *y;
  int e_last = *e;
  for(int i=0; i<refine; i++) {

    if(!R_field->eval(*x,*phi,*y,(m3dc1_field::m3dc1_get_op)op,R0,e)) {
      if(i==0) return false;

      // If we were inside the domain on the last step but not in this step,
      // then backtrack halfway and try again.
      //      std::cerr << "backtracking..." << std::endl;
      *x = (*x + x_last) / 2.;
      *y = (*y + y_last) / 2.;
      *e = e_last;
      continue;
    }

    if(!Z_field->eval(*x,*phi,*y,(m3dc1_field::m3dc1_get_op)op,Z0,e))
      return false;

    double dR = R - R0[m3dc1_field::OP_1];
    double dZ = Z - Z0[m3dc1_field::OP_1];
    if(abs(dR) < tol && abs(dZ) < tol)
      return true;

    // store last successful evaluation
    e_last = *e;
    x_last = *x;
    y_last = *y;

    /*
    R(x,y) = R(x0+dx,y0+dy) = R0 + (dR/dx)dx + (dR/dy)dy
    Z(x,y) = ...            = Z0 + (dZ/dx)dx + (dZ/dy)dy
    ( R - R0 ) = ( dR/dx  dR/dy ) ( dx )
    ( Z - Z0 ) = ( dZ/dx  dZ/dy ) ( dy )

    ( dx ) = _______________1________________ (  dZ/dy  -dR/dy ) ( R - R0 )
    ( dy ) = (dR/dx)(dZ/dy) - (dR/dy)(dZ/dx)) ( -dZ/dx   dR/dx ) ( Z - Z0 )
    */
    double det = R0[m3dc1_field::OP_DR]*Z0[m3dc1_field::OP_DZ]
      -          R0[m3dc1_field::OP_DZ]*Z0[m3dc1_field::OP_DR];
    if(abs(det) < min_det) {
      std::cerr << "determinant is too small.  aborting Newton iterations"
		<< std::endl;
      break;
    }

    double dx = ( Z0[m3dc1_field::OP_DZ]*dR - R0[m3dc1_field::OP_DZ]*dZ)/det;
    double dy = (-Z0[m3dc1_field::OP_DR]*dR + R0[m3dc1_field::OP_DR]*dZ)/det;

    double max_step = (mesh->a[*e] + mesh->b[*e] + mesh->c[*e])/12.;

    if(dx > max_step) {
      dy *= max_step/dx;
      dx  = max_step;
    } else if(dx < -max_step) {
      dy *= -max_step/dx;
      dx  = -max_step;
    }
    if(dy > max_step) {
      dx *= max_step/dy;
      dy  = max_step;
    } else if(dy < -max_step) {
      dx *= -max_step/dy;
      dy  = -max_step;
    }

    *x += dx;
    *y += dy;
  }

  // Make sure the last step was a successful one
  *x = x_last;
  *y = y_last;
  *e = e_last;

  if(!R_field->eval(*x,*phi,*y,m3dc1_field::GET_VAL,R0,e))
    return false;
  if(!Z_field->eval(*x,*phi,*y,m3dc1_field::GET_VAL,Z0,e))
    return false;

  /*
  std::cerr << "( " << R  << ", " << Z  << ") /"
	    << "( " << R0[m3dc1_field::OP_1] << ", "
	    << Z0[m3dc1_field::OP_1] << ")"
	    << std::endl;
  */
  return true;
}

bool m3dc1_coord_map::eval_map_deriv(const double x, const double phi, const double y, double* R, double* Z, int* e) const
{
  if(!R_field->eval(x,phi,y,m3dc1_field::GET_ALL,R,e))
    return false;
  if(!Z_field->eval(x,phi,y,m3dc1_field::GET_ALL,Z,e))
    return false;

  return true;
}
