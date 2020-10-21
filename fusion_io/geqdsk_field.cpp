#include "fusion_io.h"
#include "interpolate.h"


int geqdsk_current_density::eval(const double* x, double* j, void*)
{
  double psi[6];
  int ierr;

  ierr = source->interpolate_psi(x[0], x[2], psi);
  if(ierr != FIO_SUCCESS) return ierr;

  j[1] =  (psi[3] + psi[5] - psi[1]/x[0])/x[0];

  double f, ffp;
  cubic_interpolation(source->nw,source->psi,psi[0],source->fpol,
		       &f  );
  cubic_interpolation(source->nw,source->psi,psi[0],source->ffprime, 
		       &ffp);

  double fp = ffp/f;
  j[0] = -fp*psi[2]/x[0];
  j[2] =  fp*psi[1]/x[0]; 

  j[0] /= M_PI*4e-7;
  j[1] /= M_PI*4e-7;
  j[2] /= M_PI*4e-7;

  return FIO_SUCCESS;
}

int geqdsk_magnetic_field::eval(const double* x, double* b, void*)
{
  double psi[6];
  int ierr;

  ierr = source->interpolate_psi(x[0], x[2], psi);
  if(ierr != FIO_SUCCESS) return ierr;

  b[0] =  psi[2]/x[0];
  b[2] = -psi[1]/x[0];

  double f;
  cubic_interpolation(source->nw, source->psi, psi[0], source->fpol, &f);
  b[1] =  f/x[0];
  
  return FIO_SUCCESS;
}

int geqdsk_pressure_field::eval(const double* x, double* p, void*)
{
  double psi[6];
  int ierr;

  ierr = source->interpolate_psi(x[0], x[2], psi);
  if(ierr != FIO_SUCCESS) return ierr;

  cubic_interpolation(source->nw, source->psi, psi[0], source->press, p);
  return FIO_SUCCESS;
}

int geqdsk_psi_field::eval(const double* x, double* b, void*)
{
  double psi[6];
  int ierr;

  ierr = source->interpolate_psi(x[0], x[2], psi);
  if(ierr != FIO_SUCCESS) return ierr;

  psi[0] = (psi[0] - offset) * factor;
  
  *b = psi[0];
  
  return FIO_SUCCESS;
}

int geqdsk_psi_field::eval_deriv(const double* x, double *db, void*)
{
  double psi[6];
  int ierr;

  ierr = source->interpolate_psi(x[0], x[2], psi);
  if(ierr != FIO_SUCCESS) return ierr;

  db[FIO_DR  ] = psi[1]*factor;
  db[FIO_DPHI] = 0.;
  db[FIO_DZ  ] = psi[2]*factor;
  
  return FIO_SUCCESS;
}


