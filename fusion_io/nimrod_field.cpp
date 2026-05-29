#ifdef FUSIONIO_ENABLE_NIMROD

#include "nimrod_field.h"
#include "fusion_io.h"
#include <iostream>
#include <cmath>

extern "C" {
  void c_eval_field_dm(const double* RZP, const char* fname, double* fval, 
                       const int* nqty, double* xyguess, int* ifail, const int* dmode, const int* eq);
}

int nimrod_field::eval(const double* coords, double* v, fio_hint s)
{
  double RZP[3];
  RZP[0] = coords[0];    // R
  RZP[1] = coords[2];    // Z
  RZP[2] = -coords[1];   // Phi (NIMROD phi = -fusion phi)

  double val[3];
  int ifail = 0;
  int dmode = 0;
  double* hint = (double*)s;
  double dummy_hint[2] = {-1.0, -1.0};
  if(!hint) hint = dummy_hint;

  c_eval_field_dm(RZP, name.c_str(), val, &nqty, hint, &ifail, &dmode, &eq);

  if(ifail != 0) return FIO_OUT_OF_BOUNDS;

  if(nqty == 3) {
    v[0] = val[0];         // R
    v[1] = -val[2];        // Phi (sign swap because phi_nim = -phi_fio)
    v[2] = val[1];         // Z
  } else {
    v[0] = val[0];
  }

  return FIO_SUCCESS;
}

int nimrod_field::eval_deriv(const double* coords, double* v, fio_hint s)
{
  double RZP[3];
  RZP[0] = coords[0];
  RZP[1] = coords[2];
  RZP[2] = -coords[1];

  double val[12]; // nqty * 4 for dmode=1
  int ifail = 0;
  int dmode = 1;
  double* hint = (double*)s;
  double dummy_hint[2] = {-1.0, -1.0};
  if(!hint) hint = dummy_hint;

  c_eval_field_dm(RZP, name.c_str(), val, &nqty, hint, &ifail, &dmode, &eq);

  if(ifail != 0) return FIO_OUT_OF_BOUNDS;

  // val[0...nqty-1] is value
  // val[nqty...2nqty-1] is d/dR
  // val[2nqty...3nqty-1] is d/dZ
  // val[3nqty...4nqty-1] is d/dPhi_nimrod (which is (1/R) d/dphi_nimrod)

  const double R = coords[0];

  if(nqty == 1) {
    // v[FIO_DR], v[FIO_DPHI], v[FIO_DZ]
    v[0] = val[1];               // d/dR
    v[1] = -R * val[3];          // d/dphi_fio = d/dphi_nim * dphi_nim/dphi_fio = (R*val[3])*(-1)
    v[2] = val[2];               // d/dZ
  } else if(nqty == 3) {
    // Partial derivatives: first index = partial deriv, second index = coord
    // v[0] = d/dR (R)
    // v[1] = d/dR (Phi)
    // v[2] = d/dR (Z)
    // v[3] = d/dPhi (R)
    // v[4] = d/dPhi (Phi)
    // v[5] = d/dPhi (Z)
    // v[6] = d/dZ (R)
    // v[7] = d/dZ (Phi)
    // v[8] = d/dZ (Z)

    // d/dR
    v[0] = val[3];               // d/dR (R)
    v[1] = -val[5];              // d/dR (Phi)
    v[2] = val[4];               // d/dR (Z)

    // d/dPhi_fio = -d/dPhi_nimrod
    // Note: val[9...11] is (1/R) d/dphi_nimrod
    v[3] = -R * val[9];          // d/dPhi (R)
    v[4] =  R * val[11];         // d/dPhi (Phi) -> -(-R*val[11]) = R*val[11]
    v[5] = -R * val[10];         // d/dPhi (Z)

    // d/dZ
    v[6] = val[6];               // d/dZ (R)
    v[7] = -val[8];              // d/dZ (Phi)
    v[8] = val[7];               // d/dZ (Z)
  }

  return FIO_SUCCESS;
}

#endif // FUSIONIO_ENABLE_NIMROD
