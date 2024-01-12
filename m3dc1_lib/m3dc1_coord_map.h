#ifndef M3DC1_COORD_MAP_H
#define M3DC1_COORD_MAP_H

#include "m3dc1_mesh.h"
#include "m3dc1_field.h"

struct m3dc1_mesh_element {
  double R[6];
  double Phi[2];
  double Z[6];
  bool is_in_element(const double r, const double phi, const double z,
		     double *xi_frac, double *zi_frac, double *eta_frac) {
    // determine whether phi falls inside element
    if(phi < Phi[0]) return false;
    if(phi > Phi[1]) return false;
    double d = (phi - Phi[0]) / (Phi[1] - Phi[0]);

    // define the triangle in the plane of phi
    double R0 = R[3]*d - R[0]*(1.-d);
    double R1 = R[4]*d - R[1]*(1.-d);
    double R2 = R[5]*d - R[2]*(1.-d);
    double Z0 = Z[3]*d - Z[0]*(1.-d);
    double Z1 = Z[4]*d - Z[1]*(1.-d);
    double Z2 = Z[5]*d - Z[2]*(1.-d);

    const double tol = 1e-6;
    if((r - R0)*(Z1 - Z0) - (z - Z0)*(R1 - R0) > tol) return false;
    if((r - R1)*(Z2 - Z1) - (z - Z1)*(R2 - R1) > tol) return false;
    if((r - R2)*(Z0 - Z1) - (z - Z2)*(R0 - R2) > tol) return false;

    // return "fractional" local coordinates
    *xi_frac = ((r - R2)*(R1 - R0) + (z - Z2)*(Z1 - Z0)) /
      ((R1 - R0)*(R1 - R0) + (Z1 - Z0)*(Z1 - Z0));
    *zi_frac = d;
    *eta_frac = ((r - R0)*(Z1 - Z0) - (z - Z0)*(R1 - R0)) /
      ((R2 - R0)*(Z1 - Z0) - (Z2 - Z0)*(R1 - R0));

    return true;
  }
  
};

class m3dc1_coord_map {
  m3dc1_field* R_field;
  m3dc1_field* Z_field;

  m3dc1_3d_mesh* mesh;
  m3dc1_mesh_element* elm;     // data 

 public:
  m3dc1_coord_map(m3dc1_3d_mesh*);
  virtual ~m3dc1_coord_map();

  bool load(m3dc1_field* R, m3dc1_field* Z);

  bool find_coordinates(const double R, const double Phi, const double Z,
			double *x, double* phi, double* y, int* e, int refine) const;
  bool find_element(const double R, const double Phi, const double Z,
		    double *xi_frac, double *zi_frac, double *eta_frac, int* e) const;
    
};

#endif
