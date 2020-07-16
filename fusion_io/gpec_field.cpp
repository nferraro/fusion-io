#include "fusion_io.h"
#include "interpolate.h"


int gpec_series::bounds(double* tmin, double* tmax) const 
{
  *tmin = 0.;
  *tmax = 0.;

  return FIO_SUCCESS;
}

int gpec_series::eval(const double t, double* x)
{
  *x = data;
  return FIO_SUCCESS;
}

int gpec_magnetic_field::load(const fio_option_list* opt)
{
  int result;

  opt->get_option(FIO_LINEAR_SCALE, &linfac);
  opt->get_option(FIO_PART, &ilin);

  if(ilin != FIO_PERTURBED_ONLY) {
    result = source->read_field_data("gpec_eqbrzphi_n1.out", &source->b0);
    if(result != FIO_SUCCESS) return result;
  }
  if(ilin != FIO_EQUILIBRIUM_ONLY) {
    result = source->read_field_data("gpec_brzphi_n1.out", &source->b1);
    if(result != FIO_SUCCESS) return result;
  }
}

int gpec_magnetic_field::eval(const double* x, double* b, void*)
{
  double psi[6];
  int ierr;
  double vr[3], vi[3];
  double br[3], bi[3];

  for(int i=0; i<3; i++) {
    br[i] = 0.;
    bi[i] = 0.;
  }

  // Read equilibrium part
  if(ilin != FIO_PERTURBED_ONLY) {
    ierr = source->b0.interpolate(x[0], x[2], vr, vi);
    if(ierr != FIO_SUCCESS) return ierr;
    
    // GPEC data is stored in (BR, BZ, BPhi) order, so we re-order it here
    b[0] = vr[0];
    b[1] = vr[2];
    b[2] = vr[1];
  }

  // Read perturbed part
  if(ilin != FIO_EQUILIBRIUM_ONLY) {
    ierr = source->b1.interpolate(x[0], x[2], vr, vi);
    if(ierr != FIO_SUCCESS) return ierr;

    double phase = source->b1.ntor*x[1]*M_PI/180.;

    // GPEC data is stored in (BR, BZ, BPhi) order, so we re-order it here
    b[0] += linfac*(vr[0]*cos(phase) - vi[0]*sin(phase));
    b[1] += linfac*(vr[2]*cos(phase) - vi[2]*sin(phase));
    b[2] += linfac*(vr[1]*cos(phase) - vi[1]*sin(phase));
  }

  return FIO_SUCCESS;
}



