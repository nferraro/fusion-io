#include "m3dc1_stell_field.h"

#include <iostream>

m3dc1_stell_field::m3dc1_stell_field(m3dc1_mesh* mesh_in, m3dc1_coord_map* map_in)
  : m3dc1_3d_field(mesh_in)
{
  map = map_in;
}

bool m3dc1_stell_field::eval(const double R, const double Phi, const double Z, 
			     const m3dc1_field::m3dc1_get_op op, double* val, 
			     int* element)
{
  double x, phi, y;
  int refine=1;
  
  // find logical coordinates associated with R, Phi, Z
  if(!map->find_coordinates(R,Phi,Z,&x,&phi,&y,element,refine))
    return false;

  // evaluate field on logical mesh
  if(!m3dc1_3d_field::eval(x,phi,y,op,val,element))
    return false;
  
  return true;
}
