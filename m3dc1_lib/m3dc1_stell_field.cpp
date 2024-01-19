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
  double R0[OP_NUM], Z0[OP_NUM];
  double x, phi, y;
  int refine=20;
  
  // find logical coordinates associated with R, Phi, Z
  if(!map->find_coordinates(R,Phi,Z,&x,&phi,&y,element,refine))
    return false;

  // evaluate field on logical mesh
  if(!m3dc1_3d_field::eval(x,phi,y,op,val,element))
    return false;

  if(op==GET_VAL)
    return true;

  // translate derivatives (if requested)
  // This follows Zhou et al 2021 Nucl. Fusion 61 086015

  if(!map->eval_map_deriv(x,phi,y,R0,Z0,element)) {
    std::cerr << "Error evaluating map derivatives" << std::endl;
    return false;
  }

  // derivatives of v, R0, and Z0 here are wrt to logical coordinates
  double v[OP_NUM];
  for(int i=0; i<OP_NUM; i++) {
    v[i] = val[i];
    if(i != OP_1) val[i] = 0;
  }

  double A = R0[OP_DZ]*Z0[OP_DP] - R0[OP_DP]*Z0[OP_DZ];
  double B = R0[OP_DP]*Z0[OP_DR] - R0[OP_DR]*Z0[OP_DP];
  double D = R0[OP_DR]*Z0[OP_DZ] - R0[OP_DZ]*Z0[OP_DR];

  // First derivatives
  val[OP_DR] = ( Z0[OP_DZ]*v[OP_DR] - Z0[OP_DR]*v[OP_DZ])/D;
  val[OP_DZ] = (-R0[OP_DZ]*v[OP_DR] + R0[OP_DR]*v[OP_DZ])/D;
  val[OP_DP] = v[OP_DP] + (A*v[OP_DR] + B*v[OP_DZ])/D;

  // Second derivatives
  double F = R0[OP_DR]*(R0[OP_DRZ]*Z0[OP_DZ ] - R0[OP_DZZ]*Z0[OP_DR ] +
			R0[OP_DR ]*Z0[OP_DZZ] - R0[OP_DZ ]*Z0[OP_DRZ])
    -        R0[OP_DZ]*(R0[OP_DRR]*Z0[OP_DZ ] - R0[OP_DRZ]*Z0[OP_DR ] +
			R0[OP_DR ]*Z0[OP_DRZ] - R0[OP_DZ ]*Z0[OP_DRR]);
  double G = Z0[OP_DZ]*(R0[OP_DRR]*Z0[OP_DZ ] - R0[OP_DRZ]*Z0[OP_DR ] +
			R0[OP_DR ]*Z0[OP_DRZ] - R0[OP_DZ ]*Z0[OP_DRR])
    -        Z0[OP_DR]*(R0[OP_DRZ]*Z0[OP_DZ ] - R0[OP_DZZ]*Z0[OP_DR ] +
			R0[OP_DR ]*Z0[OP_DZZ] - R0[OP_DZ ]*Z0[OP_DRZ]);

  // note that the following reference some real-space derivatives (val).
  val[OP_DRR] = (Z0[OP_DZ]*Z0[OP_DZ]*v[OP_DRR] + Z0[OP_DR]*Z0[OP_DR]*v[OP_DZZ]
		 - 2.*Z0[OP_DR]*Z0[OP_DZ]*v[OP_DRZ] - G*val[OP_DR]
		 + (Z0[OP_DZ]*Z0[OP_DRZ] - Z0[OP_DR]*Z0[OP_DZZ])*v[OP_DR])/(D*D);
  val[OP_DZZ] = (R0[OP_DZ]*R0[OP_DZ]*v[OP_DRR] + R0[OP_DR]*R0[OP_DR]*v[OP_DZZ]
		 - 2.*R0[OP_DR]*R0[OP_DZ]*v[OP_DRZ] - F*val[OP_DZ]
		 + (R0[OP_DZ]*R0[OP_DRZ] - R0[OP_DR]*R0[OP_DZZ])*v[OP_DR])/(D*D);
  val[OP_DRZ] = ((R0[OP_DR]*Z0[OP_DZ] + R0[OP_DZ]*Z0[OP_DR])*v[OP_DRZ]
		 - R0[OP_DZ]*Z0[OP_DZ]*v[OP_DRR] - R0[OP_DR]*Z0[OP_DR]*v[OP_DZZ]
		 - G*val[OP_DZ]
		 -(Z0[OP_DZ]*R0[OP_DRZ] - Z0[OP_DR]*R0[OP_DZZ])*v[OP_DR]
		 -(Z0[OP_DR]*R0[OP_DRZ] - Z0[OP_DZ]*R0[OP_DRR])*v[OP_DZ])/(D*D);

  double AP = R0[OP_DZP]*Z0[OP_DP ] - R0[OP_DPP]*Z0[OP_DZ ]
    +         R0[OP_DZ ]*Z0[OP_DPP] - R0[OP_DP ]*Z0[OP_DZP];
  double BP = R0[OP_DPP]*Z0[OP_DR ] - R0[OP_DRP]*Z0[OP_DP ]
    +         R0[OP_DP ]*Z0[OP_DRP] - R0[OP_DR ]*Z0[OP_DPP];
  double DP = R0[OP_DRP]*Z0[OP_DZ ] - R0[OP_DZP]*Z0[OP_DR ]
    +         R0[OP_DR ]*Z0[OP_DZP] - R0[OP_DZ ]*Z0[OP_DRP];

  // note that the following reference some real-space derivatives (val).
  val[OP_DRP] = -R0[OP_DP]*val[OP_DRR] - Z0[OP_DP]*val[OP_DRZ]
    + ( Z0[OP_DZP]*v[OP_DR ] - Z0[OP_DRP]*v[OP_DZ ]
      + Z0[OP_DZ ]*v[OP_DRP] - Z0[OP_DR ]*v[OP_DZP])/D
    - ( Z0[OP_DZ ]*v[OP_DR ] - Z0[OP_DR ]*v[OP_DZ ])*DP/(D*D);
  val[OP_DZP] = -R0[OP_DP]*val[OP_DRZ] - Z0[OP_DP]*val[OP_DZZ]
    + ( Z0[OP_DZP]*v[OP_DR ] - Z0[OP_DRP]*v[OP_DZ ]
      + Z0[OP_DZ ]*v[OP_DRP] - Z0[OP_DR ]*v[OP_DZP])/D
    - ( Z0[OP_DZ ]*v[OP_DR ] - Z0[OP_DR ]*v[OP_DZ])*DP/(D*D);
  val[OP_DPP] = -R0[OP_DP]*val[OP_DRP] - Z0[OP_DP]*val[OP_DZP]
    + v[OP_DPP] + (A*v[OP_DRP] + B*v[OP_DZP] + AP*v[OP_DR] + BP*v[OP_DZ])/D;

  // Second derivatives
  if(((op & GET_DDVAL) == GET_DDVAL) &&
     ((op & GET_PPVAL) == GET_PPVAL)) {
    std::cerr << "Warning: cannot evaluate both poloidal and toroidal second derivatives in stellarator geometry" << std::endl;
  }

  return true;
}
