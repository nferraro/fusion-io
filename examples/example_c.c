//==========================================================================
// fio_example
// ~~~~~~~~~~~
// In this example, several files from linear calculations are loaded.
// The fields are evaluated at several points, and the results are summed.
//==========================================================================

#include <fusion_io_defs.h>
#include <fusion_io_c.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int isrc;
int ipres, ine, ini;
int imag, ij, ia;
void* hint;

void free_sources();

int main()
{
  // factor by which to multiply perturbed (linear) part of solution
  const double factor = 1.;

  // time slice to read
  const int timeslice = 0;

  char filename[256];

  double p, ne, ni, b[3], x[3], x2[3], curr[3], db[9], db_n[9], a[3];
  double t1, t3[3], t32[3], t9[9];
  const double epsilon = 1e-6;

  const int npts = 10;
  double R0, R1, Z0, Z1, phi0, phi1;
  int i, j, ierr;

  double psi_axis, psi_lcfs, r0, z0, tmin, tmax;
  int ipsi_axis, ipsi_lcfs, ir0, iz0;

  int isz;

  strcpy(filename, "C1.h5");

  // read files and fields
  printf("Reading %s\n", filename);
  ierr = fio_open_source(FIO_M3DC1_SOURCE, filename, &(isrc));
  if(ierr != 0) {
    printf("Error reading %s\n", filename);
    free_sources();
    return 1;
  }
    
  isz = fio_sizeof_search_hint(isrc);
  hint = calloc(1, isz);

  // Set options appropriate to this source
  ierr = fio_get_options(isrc);
  ierr = fio_set_int_option(FIO_TIMESLICE, timeslice);
  ierr = fio_set_real_option(FIO_LINEAR_SCALE, factor);

  // For first file, read full fields (equilibrium + perturbed)
  // For subsequent file, read only perturbed parts
  // if(i > 0) ierr = fio_set_int_option(FIO_PART, FIO_PERTURBED_ONLY);
  ierr = fio_set_int_option(FIO_PART, FIO_PERTURBED_ONLY);

  // read fields
  // magnetic field and total pressure are species-independent
  ierr = fio_get_field(isrc, FIO_TOTAL_PRESSURE, &(ipres));
  ierr = fio_get_field(isrc, FIO_MAGNETIC_FIELD, &(imag));
  ierr = fio_get_field(isrc, FIO_CURRENT_DENSITY, &(ij));
  ierr = fio_get_field(isrc, FIO_VECTOR_POTENTIAL, &(ia));

  // density is species dependent; specify electrons
  ierr = fio_set_int_option(FIO_SPECIES, FIO_ELECTRON);
  ierr = fio_get_field(isrc, FIO_DENSITY, &(ine));
    
  ierr = fio_set_int_option(FIO_SPECIES, FIO_MAIN_ION);
  ierr = fio_get_field(isrc, FIO_DENSITY, &(ini));
 
  ierr = fio_get_series(isrc, FIO_MAGAXIS_PSI, &ipsi_axis);
  ierr = fio_get_series(isrc, FIO_LCFS_PSI, &ipsi_lcfs);
  ierr = fio_get_series(isrc, FIO_MAGAXIS_R, &ir0);
  ierr = fio_get_series(isrc, FIO_MAGAXIS_Z, &iz0);
  printf("ipsi_axis: %d\n", ipsi_axis);
  printf("ipsi_lcfs: %d\n", ipsi_lcfs);

  ierr = fio_eval_series(ipsi_axis, 0., &psi_axis);
  ierr = fio_eval_series(ipsi_lcfs, 0., &psi_lcfs);
  ierr = fio_eval_series(ir0, 0., &r0);
  ierr = fio_eval_series(iz0, 0., &z0);
  printf("Psi at magnetic axis: %g\n", psi_axis);
  printf("Psi at lcfs: %g\n", psi_lcfs);
  printf("R at magnetic axis: %g\n", r0);
  printf("Z at magnetic axis: %g\n", z0);

  ierr = fio_get_series_bounds(ipsi_axis, &tmin, &tmax);
  printf("tmin, tmax: %g %g\n", tmin, tmax);

  ierr = fio_close_series(ipsi_axis);
  ierr = fio_close_series(ipsi_lcfs);
  ierr = fio_close_series(ir0);
  ierr = fio_close_series(iz0);


  R0 = 1.6;
  R1 = 1.6;
  Z0 = 0.4;
  Z1 = 0.4;
  phi0 = 0.;
  phi1 = 350.;

  for(i=0; i<npts; i++) {
    x[0] = R0 + (R1-R0)*(double)i/(double)(npts-1);
    x[1] = phi0 + (phi1-phi0)*(double)i/(double)(npts-1);
    x[2] = Z0 + (Z1-Z0)*(double)i/(double)(npts-1);

    printf("(%g, %g, %g)\n", x[0], x[1], x[2]);
    
    p = 0.;
    ne = 0.;
    ni = 0.;
    for(j=0; j<3; j++) {
      a[j] = 0.;
      b[j] = 0.;
      curr[j] = 0.;
    }
    for(j=0; j<9; j++) { 
      db[j] = 0.;
      db_n[j] = 0.;
    }

    ierr = fio_eval_field(ipres, x, &t1, hint);
    p += t1;
    ierr = fio_eval_field(ine, x, &t1, hint);
    ne += t1;
    ierr = fio_eval_field(ini, x, &t1, hint);
    ni += t1;
    ierr = fio_eval_field(imag, x, t3, hint);
    b[0] += t3[0];
    b[1] += t3[1];
    b[2] += t3[2];
    ierr = fio_eval_field(ia, x, t3, hint);
    a[0] += t3[0];
    a[1] += t3[1];
    a[2] += t3[2];
    
    // calculate numerical derivatives of B
    // dR
    x2[0] = x[0] + epsilon;  x2[1] = x[1];  x2[2] = x[2];
    ierr = fio_eval_field(imag, x2, t32, hint);
    db_n[FIO_DR_R] += (t32[0] - t3[0]) / epsilon;
    db_n[FIO_DR_PHI] += (t32[1] - t3[1]) / epsilon;
    db_n[FIO_DR_Z] += (t32[2] - t3[2]) / epsilon;
    // dPHI
    x2[0] = x[0];  x2[1] = x[1] + epsilon;  x2[2] = x[2];
    ierr = fio_eval_field(imag, x2, t32, hint);
    db_n[FIO_DPHI_R] += (t32[0] - t3[0]) / epsilon;
    db_n[FIO_DPHI_PHI] += (t32[1] - t3[1]) / epsilon;
    db_n[FIO_DPHI_Z] += (t32[2] - t3[2]) / epsilon;
    // dZ
    x2[0] = x[0];  x2[1] = x[1];  x2[2] = x[2] + epsilon;
    ierr = fio_eval_field(imag, x2, t32, hint);
    db_n[FIO_DZ_R] += (t32[0] - t3[0]) / epsilon;
    db_n[FIO_DZ_PHI] += (t32[1] - t3[1]) / epsilon;
    db_n[FIO_DZ_Z] += (t32[2] - t3[2]) / epsilon;     
    
    ierr = fio_eval_field(ij, x, t3, hint);
    curr[0] += t3[0];
    curr[1] += t3[1];
    curr[2] += t3[2];
    
    ierr = fio_eval_field_deriv(imag, x, t9, hint);
    for(j=0; j<9; j++)
      db[j] += t9[j];
    
  }

  printf("\tpressure = %g\n", p);
  printf("\telectron density = %g\n", ne);
  printf("\tion density = %g\n", ni);
  printf("\ttotal B = (%g, %g, %g)\n", b[0], b[1], b[2]);
  printf("\ttotal J = (%g, %g, %g)\n", curr[0], curr[1], curr[2]);
  printf("\ttotal A = (%g, %g, %g)\n", a[0], a[1], a[2]);
    /*
    printf("Interpolated derivatives: \n");
    printf("\tdB/dR = (%g, %g, %g)\n", 
	   db[FIO_DR_R], db[FIO_DR_PHI], db[FIO_DR_Z]);
    printf("\tdB/dphi = (%g, %g, %g)\n", 
	   db[FIO_DPHI_R], db[FIO_DPHI_PHI], db[FIO_DPHI_Z]);
    printf("\tdB/dZ = (%g, %g, %g)\n", 
	   db[FIO_DZ_R], db[FIO_DZ_PHI], db[FIO_DZ_Z]);

    printf("Numerical derivatives: \n");
    printf("\tdB/dR = (%g, %g, %g)\n", 
	   db_n[FIO_DR_R], db_n[FIO_DR_PHI], db_n[FIO_DR_Z]);
    printf("\tdB/dphi = (%g, %g, %g)\n", 
	   db_n[FIO_DPHI_R], db_n[FIO_DPHI_PHI], db_n[FIO_DPHI_Z]);
    printf("\tdB/dZ = (%g, %g, %g)\n", 
	   db_n[FIO_DZ_R], db_n[FIO_DZ_PHI], db_n[FIO_DZ_Z]);
    */

  free_sources();

  return 0;
}

void free_sources()
{
  int i, ierr;
  free(hint);
  ierr = fio_close_field(ine);
  ierr = fio_close_field(ini);
  ierr = fio_close_field(ipres);
  ierr = fio_close_field(imag);
  ierr = fio_close_field(ij);
  ierr = fio_close_field(ia);
  ierr = fio_close_source(isrc);
}
