#ifndef M3DC1_COORD_MAP_H
#define M3DC1_COORD_MAP_H

#include "m3dc1_mesh.h"
#include "m3dc1_field.h"

struct m3dc1_mesh_element {
  double R[6];
  double Phi[2];
  double Z[6];
  bool is_in_element(const double r, const double phi, const double z,
		     double *xi_frac, double *zi_frac, double *eta_frac) const;
  
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
