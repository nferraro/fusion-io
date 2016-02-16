#include "fusion_io.h"
#include <iostream>

int mars_field::load(const fio_option_list* opt)
{
  opt->get_option(FIO_TIMESLICE, &time);
  opt->get_option(FIO_LINEAR_SCALE, &linfac);
  opt->get_option(FIO_PART, &part);
  opt->get_option(FIO_PHASE, &phase);

  if(linfac != 1. && part==FIO_EQUILIBRIUM_ONLY) {
    std::cerr << "Linear scale linfac is ignored for equilibrium data." 
	      << std::endl;
    linfac = 1.;
  }
  

  return FIO_SUCCESS;
}

int mars_field::get_indices(const double* x, int* n, int* m, 
			    double* dn, double *dm)
{
  //  std::cerr << "Point ( " << x[0] << ", " << x[1] << ", " << x[2] << ")";

  double det, drdi, drdj, dzdi, dzdj;

  for(int i=0; i<source->nrp1-1; i++) {
    for(int j=0; j<source->nchi; j++) {
      int jj = (j+1) % source->nchi;
      if(x[0] < source->eqdata->req[i][j] &&
	 x[0] < source->eqdata->req[i][jj] &&
	 x[0] < source->eqdata->req[i+1][j] &&
	 x[0] < source->eqdata->req[i+1][jj])
	continue;
      if(x[0] > source->eqdata->req[i][j] &&
	 x[0] > source->eqdata->req[i][jj] &&
	 x[0] > source->eqdata->req[i+1][j] &&
	 x[0] > source->eqdata->req[i+1][jj])
	continue;
      if(x[2] < source->eqdata->zeq[i][j] &&
	 x[2] < source->eqdata->zeq[i][jj] &&
	 x[2] < source->eqdata->zeq[i+1][j] &&
	 x[2] < source->eqdata->zeq[i+1][jj])
	continue;
      if(x[2] > source->eqdata->zeq[i][j] &&
	 x[2] > source->eqdata->zeq[i][jj] &&
	 x[2] > source->eqdata->zeq[i+1][j] &&
	 x[2] > source->eqdata->zeq[i+1][jj])
	continue;
      //      std::cerr << " found in cell " << i << ", " << j << std::endl;

      *n = i;
      *m = j;

      drdi = ((source->eqdata->req[i+1][j ]-source->eqdata->req[i  ][j ])
             +(source->eqdata->req[i+1][jj]-source->eqdata->req[i  ][jj]))/2.;
      drdj = ((source->eqdata->req[i  ][jj]-source->eqdata->req[i  ][j ])
	     +(source->eqdata->req[i+1][jj]-source->eqdata->req[i+1][j ]))/2.;
      dzdi = ((source->eqdata->zeq[i+1][j ]-source->eqdata->zeq[i  ][j ])
	     +(source->eqdata->zeq[i+1][jj]-source->eqdata->zeq[i  ][jj]))/2.;
      dzdj = ((source->eqdata->zeq[i  ][jj]-source->eqdata->zeq[i  ][j ])
	     +(source->eqdata->zeq[i+1][jj]-source->eqdata->zeq[i+1][j ]))/2.;
      
      det = drdi*dzdj - drdj*dzdi;

      *dn = (dzdj*(x[0] - source->eqdata->req[i][j]) -
	     drdj*(x[2] - source->eqdata->zeq[i][j])) / det;
      *dm = (-dzdi*(x[0] - source->eqdata->req[i][j]) -
	      drdi*(x[2] - source->eqdata->zeq[i][j])) / det;
      //      std::cerr << " (n, m) =  " << *n << ", " << *m << std::endl;
      //      std::cerr << " (dn, dm) =  " << *dn << ", " << *dm << std::endl;

      return FIO_SUCCESS;
    }
  }

  //  std::cerr << " not found" << std::endl;
  
  return FIO_OUT_OF_BOUNDS;
}


int mars_magnetic_field::load(const fio_option_list* opt)
{
  mars_field::load(opt);

  int result;

  if(part != FIO_EQUILIBRIUM_ONLY) {
    result = source->load_eigdata("BPLASMA.OUT",&(source->bplasma),false);
    if(result != FIO_SUCCESS) return result;
  }

  return FIO_SUCCESS;
}


int mars_magnetic_field::eval(const double* x, double* v)
{  
  v[0] = 0.;
  v[1] = 0.;
  v[2] = 0.;

  int n, m;
  double dn, dm;
  double tmp[3];
  int result = get_indices(x, &n, &m, &dn, &dm);
  if(result != FIO_SUCCESS)
    return result;

  double chi = ((double)m + dm) * (2.*M_PI) / source->nchi;

  /*
  // add the equilibrium part
  if(part != FIO_PERTURBED_ONLY) {
    double t;
    if(n < source->bplasma->nrp-1) {
      t = source->eqdata->t[n]*(1.-dn)
	+ source->eqdata->t[n+1]*dn;
    } else {
      t = source->eqdata->t[n];
    }
    v[1] += t/x[0];
  }
  */

  // add the perturbed part
  if(part != FIO_EQUILIBRIUM_ONLY) {
    tmp[0] = 0.;
    tmp[1] = 0.;
    tmp[2] = 0.;
    for(int i=0; i<source->bplasma->maxm-1; i++) {
      for(int k=0; k<3; k++) { 
	double val_r, val_i;
	val_r = source->bplasma->val_r[k][i][n  ]*(1.-dn) 
	  +     source->bplasma->val_r[k][i][n+1]*dn; 
	val_i = source->bplasma->val_i[k][i][n  ]*(1.-dn)
	  +     source->bplasma->val_i[k][i][n+1]*dn;

	tmp[k] += ( val_r*cos(chi*source->bplasma->mnum_r[k][i])
                   -val_i*sin(chi*source->bplasma->mnum_r[k][i]));
      }
      
      v[0] += tmp[0]*linfac;
      v[1] += tmp[1]*linfac;
      v[2] += tmp[2]*linfac;
    }
  }

  return FIO_SUCCESS;
}

int mars_fluid_velocity::load(const fio_option_list* opt)
{
  mars_field::load(opt);

  return FIO_SUCCESS;
}

int mars_fluid_velocity::eval(const double* x, double* v)
{  
  return FIO_SUCCESS;
}
