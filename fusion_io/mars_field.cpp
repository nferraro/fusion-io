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
			    double* dn, double *dm, 
			    double* drdio, double* drdjo,
			    double* dzdio, double* dzdjo)
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

      if(drdio != NULL) *drdio = drdi;
      if(drdjo != NULL) *drdjo = drdj;
      if(dzdio != NULL) *dzdio = dzdi;
      if(dzdjo != NULL) *dzdjo = dzdj;

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


int mars_magnetic_field::eval(const double* x, double* v, void*)
{  
  v[0] = 0.;
  v[1] = 0.;
  v[2] = 0.;

  int n, m;
  double dn, dm, drdi, drdj, dzdi, dzdj;
  double tmp[3];
  int result = get_indices(x, &n, &m, &dn, &dm, &drdi, &drdj, &dzdi, &dzdj);
  if(result != FIO_SUCCESS)
    return result;

  int m1 = (m+1) % source->nchi;

  double dthetadj = 2.*M_PI/source->nchi;

  double dRds = drdi/
    (source->eqdata->cse[n+1] - source->eqdata->cse[n]);
  double dRdchi = drdj/dthetadj;
  double dZds = dzdi/
    (source->eqdata->cse[n+1] - source->eqdata->cse[n]);
  double dZdchi = dzdj/dthetadj;

  double g11    = (1.-dm)*(source->eqdata->g11l[n][m]*(1.-dn) + 
			   source->eqdata->g11l[n+1][m]*dn) + 
    +             dm*(source->eqdata->g11l[n][m1]*(1.-dn) +
		      source->eqdata->g11l[n+1][m1]*dn);
  double g22    = (1.-dm)*(source->eqdata->g22l[n][m]*(1.-dn) + 
			   source->eqdata->g22l[n+1][m]*dn) + 
    +             dm*(source->eqdata->g22l[n][m1]*(1.-dn) +
		      source->eqdata->g22l[n+1][m1]*dn);
  double g33    = (1.-dm)*(source->eqdata->g33l[n][m]*(1.-dn) + 
			   source->eqdata->g33l[n+1][m]*dn) + 
    +             dm*(source->eqdata->g33l[n][m1]*(1.-dn) +
		      source->eqdata->g33l[n+1][m1]*dn);
  double jac    = (1.-dm)*(source->eqdata->rja[n][m]*(1.-dn) + 
			   source->eqdata->rja[n+1][m]*dn) + 
    +             dm*(source->eqdata->rja[n][m1]*(1.-dn) +
		      source->eqdata->rja[n+1][m1]*dn);

  // add the equilibrium part
  if(part != FIO_PERTURBED_ONLY) {
    double t = source->eqdata->t[n]*(1.-dn)
      +        source->eqdata->t[n+1]*dn;
    v[1] += t/x[0];

    double dpsids = source->eqdata->dpsids[n]*(1.-dn)
      +             source->eqdata->dpsids[n+1]*dn;
    v[0] += (dpsids*dRdchi)/jac;
    v[2] += (dpsids*dZdchi)/jac;
  }

  // add the perturbed part
  if(part != FIO_EQUILIBRIUM_ONLY) {
    double chi = ((double)m + dm) * dthetadj;

    tmp[0] = 0.;
    tmp[1] = 0.;
    tmp[2] = 0.;

    // Sum over poloidal fourier modes
    for(int i=0; i<source->bplasma->maxm-1; i++) {

      // Sum over components
      for(int k=0; k<3; k++) { 
	double val_r, val_i;
	val_r = source->bplasma->val_r[k][i][n  ]*(1.-dn) 
	  +     source->bplasma->val_r[k][i][n+1]*dn; 
	val_i = source->bplasma->val_i[k][i][n  ]*(1.-dn)
	  +     source->bplasma->val_i[k][i][n+1]*dn;

	tmp[k] += ( val_r*cos(chi*source->bplasma->mnum_r[k][i])
		   -val_i*sin(chi*source->bplasma->mnum_r[k][i]));
      }
    }

    tmp[0] *= sqrt(g11)/jac;
    tmp[1] *= sqrt(g22)/jac;
    tmp[2] *= sqrt(g33)/jac;

    v[0] += linfac*(tmp[0]*dRds + tmp[1]*dRdchi)/jac;
    v[1] += linfac*tmp[2]*source->eqdata->cse[0]/jac;
    v[2] += linfac*(tmp[0]*dZds + tmp[1]*dZdchi)/jac;
  }

  // convert to SI
  v[0] *= source->b0exp;
  v[1] *= source->b0exp;
  v[2] *= source->b0exp;

  /*
  std::cerr << "Jacobian: " 
	    << (-dRdchi*dZds + dRds*dZdchi)*x[0] << "\t"
	    << jac << std::endl;
  std::cerr << "sqrt(g)               =\t"
	    << sqrt(source->eqdata->g11l[n][m]) << ",\t"
	    << sqrt(source->eqdata->g22l[n][m]) << ",\t"
	    << sqrt(source->eqdata->g33l[n][m]) << std::endl;
  std::cerr << "sqrt(g)               =\t"
	    << sqrt(g11) << ",\t"
	    << sqrt(g22) << ",\t"
	    << sqrt(g33) << std::endl;
  std::cerr << "sqrt(dRds^2 + dZds^2) =\t" 
  	    << sqrt(dRds*dRds + dZds*dZds) << ",\t" 
  	    << sqrt(dRdchi*dRdchi + dZdchi*dZdchi) << ",\t"
  	    << source->eqdata->req[n][m] << std::endl;
  */

  return FIO_SUCCESS;
}

int mars_fluid_velocity::load(const fio_option_list* opt)
{
  mars_field::load(opt);

  return FIO_SUCCESS;
}

int mars_fluid_velocity::eval(const double* x, double* v, void*)
{  
  return FIO_SUCCESS;
}
