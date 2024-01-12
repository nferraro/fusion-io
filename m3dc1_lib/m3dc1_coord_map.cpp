#include "m3dc1_coord_map.h"

#include <iostream>

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
    //    std::cerr << "i = " << i << std::endl;

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
      double x, phi, y, r, z;

      mesh->local_to_global(i,xi[j],zi[j],eta[j],&x,&phi,&y);
      
      e = i;

      std::cerr << "Before eval " << std::endl;
      if(!R_field->eval(x,phi,y,m3dc1_field::qGET_VAL,&r,&e)) {
	std::cerr << "failed;  returning" << std::endl;
	return false;
      }
      if(!Z_field->eval(x,phi,y,m3dc1_field::GET_VAL,&z,&e))
	return false;

      if(e != i) {
	std::cerr << "Warning: node found in different element" << std::endl;
	std::cerr << e << " " << i << std::endl;
      }
      /*
      r = x;
      z = y;
      */
      elm[i].R[j] = r;
      elm[i].Z[j] = z;
      if(j==0) {
	elm[i].Phi[0] = phi;
      } else if(j==3) {
	elm[i].Phi[1] = phi;
      }
    }

  }

  return true;
}

bool m3dc1_coord_map::find_element(const double R, const double Phi, const double Z,
			      double *xi_frac, double *zi_frac, double *eta_frac, int* e) const
{
  bool found;
  
  if(*e >= 0) {
    found = elm[*e].is_in_element(R, Phi, Z, xi_frac, zi_frac, eta_frac);
    if(found) return true;
  }

  for(int i=0; i<mesh->nelms; i++) {
    found = elm[i].is_in_element(R, Phi, Z, xi_frac, zi_frac, eta_frac);
    if(found) {
      *e = i;
      return true;
    }
  }

  return false;
}

bool m3dc1_coord_map::find_coordinates(const double R, const double Phi, const double Z,
				       double *x, double* phi, double* y, int* e, int refine) const
{
  double xi_frac, zi_frac, eta_frac;
  double xi, zi, eta;

  // find element containing the R, Phi, Z point
  if(!find_element(R, Phi, Z, &xi_frac, &zi_frac, &eta_frac, e)) {
    return false;
  }

  // calculate local coordinates in element
  // these coordinates have been calculated using a linear interpolation of
  // the (R, Phi, Z) values at the element nodes and should be refined
  xi = xi_frac*(mesh->a[*e] + mesh->b[*e]);
  zi = zi_frac*mesh->d[*e];
  eta = eta_frac*mesh->c[*e];

  // calculate global logical coordinates
  mesh->local_to_global(*e, xi, zi, eta, x, phi, y);
  
  // refine the coordinates through Newton iterations
  int e0 = *e;
  double R0, Z0;
  for(int i=0; i<refine; i++) {
    
    if(!R_field->eval(*x,*phi,*y,m3dc1_field::GET_VAL,&R0,&e0))
      return false;
    if(!Z_field->eval(*x,*phi,*y,m3dc1_field::GET_VAL,&Z0,&e0))
      return false;
  }

  std::cerr << "( " << R  << ", " << Z  << ") /"
	    << "( " << R0 << ", " << Z0 << ")"
	    << std::endl;
    
  return true;
}
