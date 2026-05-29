#ifndef M3DC1_H
#define M3DC1_H

#include <hdf5.h>
#include <map>
#include <vector>
#include <string>

#include "m3dc1_mesh.h"
#include "m3dc1_field.h"
#include "m3dc1_coord_map.h"

class m3dc1_stell_field : public m3dc1_3d_field {
  m3dc1_coord_map* map;
  
 public:
  m3dc1_stell_field(m3dc1_mesh*, m3dc1_coord_map*);

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element=0);
  virtual bool extent(double* X0, double* X1,
		      double* Phi0, double* Phi1,
		      double* Z0, double* Z1) const;
};

#endif
