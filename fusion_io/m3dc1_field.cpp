#include "fusion_io.h"
#include <iostream>

int m3dc1_fio_series::load()
{
  data = source->file.read_scalar(name.c_str());
  if(!data) return -1;
  time = source->file.read_scalar("time");
  if(!time) return -1;

  if(time->size() != data->size()) {
    std::cerr << "Error: series " << name
	      << "has different length than time series." << std::endl;
    return -1;
  }

  return FIO_SUCCESS;
}

int m3dc1_fio_series::bounds(double* tmin, double* tmax) const
{
  if(time->size() > 0) {
    *tmin = time->at(0) * source->t0;
    *tmax = time->at(time->size()-1) * source->t0;
    return FIO_SUCCESS;
  } else {
    return FIO_NO_DATA;
  }
}

int m3dc1_fio_series::eval(const double t, double* x)
{
  if(time->size()==0)
    return FIO_NO_DATA;
  
  if(time->size()==1) {
    *x = factor * data->at(0);
    return FIO_SUCCESS;
  }

  if(t < time->at(0) || t > time->at(time->size()-1))
    return FIO_OUT_OF_BOUNDS;

  // linearly interpolate data
  int i;
  for(i=0; i<time->size()-1; i++) {
    if(time->at(i) * source->t0 <= t && time->at(i+1) * source->t0 >= t)
      break;
  }
  *x = factor * (
		 data->at(i) + 
    (data->at(i+1)-data->at(i))*
		 (t / source->t0 - time->at(i))/(time->at(i+1)-time->at(i)));
  
  return FIO_SUCCESS;
}

int m3dc1_fio_field::load(const fio_option_list* opt)
{
  int ilin;

  opt->get_option(FIO_TIMESLICE, &time);
  opt->get_option(FIO_LINEAR_SCALE, &linfac);
  opt->get_option(FIO_PART, &ilin);
  opt->get_option(FIO_PHASE, &phase);

  if(ilin==FIO_EQUILIBRIUM_ONLY) {
    if(source->eqsubtract==1) 
      time = -1;
    else
      std::cerr << "Equilibrium Only is ignored when eqsubtract==0" 
		<< std::endl;
  }

  extsub = (source->extsubtract==1 && source->version < 8);
  eqsub = (source->eqsubtract==1) && (ilin != FIO_PERTURBED_ONLY);
  use_f = (source->i3d==1 || source->icomplex==1);

  if(time==-1) {
    eqsub = false;   // equilibrium fields need not be added in
    use_f = false;   // equilibrium is assumed axisymmetric
  }

  if(linfac != 1.) {
    if(source->linear == 0) {
      std::cerr << "Linear scale linfac is ignored for nonlinear data."
		<< std::endl;
      linfac = 1.;
    }
    if(time==-1) {
      std::cerr << "Linear scale linfac is ignored for equilibrium data." 
		<< std::endl;
      linfac = 1.;
    }
  }

  return FIO_SUCCESS;
}

int m3dc1_fio_field::get_real_parameter(const field_parameter t, double* p)
{
  switch(t) {
  case FIO_TIME:  
    if(!source->file.get_slice_time(time, p))
      return FIO_NO_DATA;
    *p *= source->t0;
    return FIO_SUCCESS;

  default:
    return FIO_UNSUPPORTED;
  }
}

int m3dc1_scalar_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  f1 = source->file.load_field(name.c_str(), time);
  if(!f1) return 1;

  if(eqsub) {
    f0 = source->file.load_field(name.c_str(), -1);
    if(!f0) return 1;
  }

  return FIO_SUCCESS;
}


int m3dc1_scalar_field::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op get = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  double val[m3dc1_field::OP_NUM];

  if(!f1->eval(x[0], x[1]-phase, x[2], get, val, (int*)s)) {
    return FIO_OUT_OF_BOUNDS;
  }
  *v = linfac*val[m3dc1_field::OP_1];

  if(eqsub) {
    if(!f0->eval(x[0], x[1]-phase, x[2], get, val, (int*)s)) {
      return FIO_OUT_OF_BOUNDS;
    }
    *v += val[m3dc1_field::OP_1];
  }

  *v *= factor;

  return FIO_SUCCESS;
}

int m3dc1_pi_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);
  
  pe = new m3dc1_scalar_field(source, "Pe", source->p0);
  int result = pe->load(opt);
  if(result != FIO_SUCCESS) return result;

  p = new m3dc1_scalar_field(source, "P", source->p0);  
  result = p->load(opt);
  if(result != FIO_SUCCESS) return result;

  return FIO_SUCCESS;
}


int m3dc1_pi_field::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op get = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  double val;
  int result;

  result = p->eval(x, v, s);
  if(result != FIO_SUCCESS) return result;

  result = pe->eval(x, &val, (int*)s);
  if(result != FIO_SUCCESS) return result;
  
  *v = *v - val;

  return FIO_SUCCESS;
}

m3dc1_pi_field::~m3dc1_pi_field()
{
  if(p) delete(p);
  if(pe) delete(pe);
}


int m3dc1_phi_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  return FIO_SUCCESS;
}

int m3dc1_phi_field::eval(const double* x, double* v, void* s)
{
  *v = 0.;

  return FIO_SUCCESS;
}

int m3dc1_electric_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  if(!(E1[0] = source->file.load_field("E_R", time)))
     return 1;
  if(!(E1[1] = source->file.load_field("E_PHI", time)))
     return 1;
  if(!(E1[2] = source->file.load_field("E_Z", time)))
     return 1;

  if(eqsub) {
    if(!(E0[0] = source->file.load_field("E_R", -1)))
      return 1;
    if(!(E0[1] = source->file.load_field("E_PHI", -1)))
      return 1;
    if(!(E0[2] = source->file.load_field("E_Z", -1)))
      return 1;
  }
     
  return FIO_SUCCESS;
}

int m3dc1_electric_field::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op get = 
    (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_VAL);

  double val[m3dc1_field::OP_NUM];

  for(int i=0; i<3; i++) {
    if(!E1[i]->eval(x[0], x[1]-phase, x[2], get, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[i] = linfac*val[m3dc1_field::OP_1];

    if(eqsub) {
      if(!E0[i]->eval(x[0], x[1]-phase, x[2], get, val, (int*)s))
	return FIO_OUT_OF_BOUNDS;
      v[i] += val[m3dc1_field::OP_1];
    }

    v[i] *= source->Phi0 / source->L0;
  }

  return FIO_SUCCESS;
}

int m3dc1_fluid_velocity::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  phi1 = source->file.load_field("phi", time);
  if(!phi1) return 1;
  w1 = source->file.load_field("V", time);
  if(!w1) return 1;
  chi1 = source->file.load_field("chi", time);
  if(!chi1) return 1;

  if(eqsub) {
    phi0 = source->file.load_field("phi", -1);
    if(!phi0) return 1;
    w0 = source->file.load_field("V", -1);
    if(!w0) return 1;
    chi0 = source->file.load_field("chi", -1);
    if(!chi0) return 1;
  }

  return FIO_SUCCESS;
}


int m3dc1_fluid_velocity::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op phiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op vget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op chiget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  double val[m3dc1_field::OP_NUM];
  double r = (source->itor==1) ? x[0] : 1.;

  // V_R   = -(dphi/dZ)*R + dchi/dR 
  // V_Z   =  (dphi/dR)*R + dchi/dZ
  // V_Phi =  V*R

  if(!phi1->eval(x[0], x[1], x[2], phiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;
  v[0] = -linfac*val[m3dc1_field::OP_DZ]*r;
  v[2] =  linfac*val[m3dc1_field::OP_DR]*r;

  if(!w1->eval(x[0], x[1], x[2], vget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;
  v[1] =  linfac*val[m3dc1_field::OP_1]*r;

  if(!chi1->eval(x[0], x[1], x[2], chiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;
  v[0] += linfac*val[m3dc1_field::OP_DR]/(r*r);
  v[2] += linfac*val[m3dc1_field::OP_DZ]/(r*r);

  if(eqsub) {
    if(!phi0->eval(x[0], x[1], x[2], phiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[0] = -val[m3dc1_field::OP_DZ]*r;
    v[2] =  val[m3dc1_field::OP_DR]*r;
    
    if(!w0->eval(x[0], x[1], x[2], vget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[1] =  val[m3dc1_field::OP_1]*r;

    if(!chi0->eval(x[0], x[1], x[2], chiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[0] += val[m3dc1_field::OP_DR]/(r*r);
    v[2] += val[m3dc1_field::OP_DZ]/(r*r);
  }

  v[0] *= source->v0;
  v[1] *= source->v0;
  v[2] *= source->v0;
  
  return FIO_SUCCESS;
}

int m3dc1_vector_potential::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  if(source->igeometry==1) {
    rst = source->file.load_field("rst", time);
    if(!rst) return 1;
    zst = source->file.load_field("zst", time);
    if(!zst) return 1;
  }

  return FIO_SUCCESS;
}


int m3dc1_vector_potential::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op fget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op rstget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op zstget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  double val[m3dc1_field::OP_NUM];
  double r = (source->itor==1) ? x[0] : 1.;

  if((source->igeometry==1)&&(source->itor==1)) {
    if(!rst->eval(x[0], x[1]-phase, x[2], rstget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    r = val[m3dc1_field::OP_1];
  }

  // A_R   =  R (df/dZ)
  // A_Z   = -R (df/dR) - F0 ln(R)
  // A_Phi = psi/R

  if(!psi1->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;

  v[1] = linfac*val[m3dc1_field::OP_1]/r;

  if(use_f) {
    if(!f1->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] =  linfac*r*val[m3dc1_field::OP_DZ];
    v[2] = -linfac*r*val[m3dc1_field::OP_DR];
  } else {
    v[0] = v[2] = 0.;
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[1] += val[m3dc1_field::OP_1]/r;

    if(source->itor==1) 
      v[2] -= source->bzero*source->rzero*log(r);
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[1] += linfac*val[m3dc1_field::OP_1]/r;

    if(use_f) {
      if(!fx->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
        return FIO_OUT_OF_BOUNDS;

      v[0] += linfac*r*val[m3dc1_field::OP_DZ];
      v[2] -= linfac*r*val[m3dc1_field::OP_DR];
    }
  }

  // convert to mks
  v[0] *= source->B0*source->L0;
  v[1] *= source->B0*source->L0;
  v[2] *= source->B0*source->L0;
  
  return FIO_SUCCESS;
}

int m3dc1_magnetic_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  i1 = source->file.load_field("I", time);
  if(!i1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
    i0 = source->file.load_field("I", -1);
    if(!i0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    ix = source->file.load_field("I_ext", time);
    if(!ix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_magnetic_field::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op gget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op fget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);

  double val[m3dc1_field::OP_NUM];
  double r = (source->itor==1) ? x[0] : 1.;

  // B_R   = -(dpsi/dZ)/R - (d2f/dRdphi)
  // B_Z   =  (dpsi/dR)/R - (d2f/dZdphi)
  // B_Phi =  F/R

  if(!psi1->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;

  v[0] = -linfac*val[m3dc1_field::OP_DZ]/r;
  v[2] =  linfac*val[m3dc1_field::OP_DR]/r;

  if(!i1->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;
  v[1] =  linfac*val[m3dc1_field::OP_1]/r;

  if(use_f) {
    if(!f1->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= linfac*val[m3dc1_field::OP_DRP];
    v[2] -= linfac*val[m3dc1_field::OP_DZP];
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= val[m3dc1_field::OP_DZ]/r;
    v[2] += val[m3dc1_field::OP_DR]/r;

    if(!i0->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[1] += val[m3dc1_field::OP_1]/r;
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= linfac*val[m3dc1_field::OP_DZ]/r;
    v[2] += linfac*val[m3dc1_field::OP_DR]/r;

    if(!ix->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[1] += linfac*val[m3dc1_field::OP_1]/r;

    if(use_f) {
      if(!fx->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
        return FIO_OUT_OF_BOUNDS;

      v[0] -= linfac*val[m3dc1_field::OP_DRP];
      v[1] -= linfac*val[m3dc1_field::OP_DZP];
    }
  }

  v[0] *= source->B0;
  v[1] *= source->B0;
  v[2] *= source->B0;
  
  return FIO_SUCCESS;
}

int m3dc1_magnetic_field::eval_deriv(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_DDVAL);

  const m3dc1_field::m3dc1_get_op gget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL | m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);

  const m3dc1_field::m3dc1_get_op fget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | 
     m3dc1_field::GET_DDVAL | m3dc1_field::GET_PPVAL);

  double val[m3dc1_field::OP_NUM];
  double r = (source->itor==1) ? x[0] : 1.;

  // B_R   = -(dpsi/dZ)/R - (d2f/dRdphi)
  // B_Z   =  (dpsi/dR)/R - (d2f/dZdphi)
  // B_Phi =  F/R

  if(!psi1->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;

  v[FIO_DR_R  ] = -linfac*val[m3dc1_field::OP_DRZ]/r;
  if(source->itor==1)
    v[FIO_DR_R] += linfac*val[m3dc1_field::OP_DZ]/(r*r);
  v[FIO_DPHI_R] = -linfac*val[m3dc1_field::OP_DZP]/r;
  v[FIO_DZ_R  ] = -linfac*val[m3dc1_field::OP_DZZ]/r;

  v[FIO_DR_Z  ] =  linfac*(val[m3dc1_field::OP_DRR]/r);
  if(source->itor==1) 
    v[FIO_DR_Z] -= linfac*val[m3dc1_field::OP_DR]/(r*r);
  v[FIO_DPHI_Z] =  linfac*val[m3dc1_field::OP_DRP]/r;
  v[FIO_DZ_Z  ] =  linfac*val[m3dc1_field::OP_DRZ]/r;

  if(!i1->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;
  v[FIO_DR_PHI  ] =  linfac*val[m3dc1_field::OP_DR]/r;
  if(source->itor==1)
    v[FIO_DR_PHI] -= linfac*val[m3dc1_field::OP_1]/(r*r);
  v[FIO_DPHI_PHI] =  linfac*val[m3dc1_field::OP_DP]/r;
  v[FIO_DZ_PHI  ] =  linfac*val[m3dc1_field::OP_DZ]/r;

  if(use_f) {
    if(!f1->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[FIO_DR_R  ] -= linfac*val[m3dc1_field::OP_DRRP];
    v[FIO_DPHI_R] -= linfac*val[m3dc1_field::OP_DRPP];
    v[FIO_DZ_R  ] -= linfac*val[m3dc1_field::OP_DRZP];

    v[FIO_DR_Z  ] -= linfac*val[m3dc1_field::OP_DRZP];
    v[FIO_DPHI_Z] -= linfac*val[m3dc1_field::OP_DZPP];
    v[FIO_DZ_Z  ] -= linfac*val[m3dc1_field::OP_DZZP];
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[FIO_DR_R  ] -= val[m3dc1_field::OP_DRZ]/r;
    if(source->itor==1) 
      v[FIO_DR_R] += val[m3dc1_field::OP_DZ]/(r*r);
    v[FIO_DPHI_R] -= val[m3dc1_field::OP_DZP]/r;
    v[FIO_DZ_R  ] -= val[m3dc1_field::OP_DZZ]/r;

    v[FIO_DR_Z  ] += val[m3dc1_field::OP_DRR]/r;
    if(source->itor==1) 
      v[FIO_DR_Z] -= val[m3dc1_field::OP_DR]/(r*r);
    v[FIO_DPHI_Z] += val[m3dc1_field::OP_DRP]/r;
    v[FIO_DZ_Z  ] += val[m3dc1_field::OP_DRZ]/r;

    if(!i0->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[FIO_DR_PHI  ] += val[m3dc1_field::OP_DR]/r;
    if(source->itor==1)
      v[FIO_DR_PHI] -= val[m3dc1_field::OP_1]/(r*r);
    v[FIO_DPHI_PHI] += val[m3dc1_field::OP_DP]/r;
    v[FIO_DZ_PHI  ] += val[m3dc1_field::OP_DZ]/r;
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[FIO_DR_R  ] -= linfac*val[m3dc1_field::OP_DRZ]/r;
    if(source->itor==1) 
      v[FIO_DR_R] += linfac*val[m3dc1_field::OP_DZ]/(r*r);
    v[FIO_DPHI_R] -= linfac*val[m3dc1_field::OP_DZP]/r;
    v[FIO_DZ_R  ] -= linfac*val[m3dc1_field::OP_DZZ]/r;

    v[FIO_DR_Z  ] += linfac*val[m3dc1_field::OP_DRR]/r;
    if(source->itor==1)
      v[FIO_DR_Z] -= linfac*val[m3dc1_field::OP_DR]/(r*r);
    v[FIO_DPHI_Z] += linfac*val[m3dc1_field::OP_DRP]/r;
    v[FIO_DZ_Z  ] += linfac*val[m3dc1_field::OP_DRZ]/r;

    if(!ix->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[FIO_DR_PHI  ] += linfac*val[m3dc1_field::OP_DR]/r; 
    if(source->itor==1)
      v[FIO_DR_PHI] -= linfac*val[m3dc1_field::OP_1]/(r*r);
    v[FIO_DPHI_PHI] += linfac*val[m3dc1_field::OP_DP]/r;
    v[FIO_DZ_PHI  ] += linfac*val[m3dc1_field::OP_DZ]/r;

    if(use_f) {
      if(!fx->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
        return FIO_OUT_OF_BOUNDS;

      v[FIO_DR_R  ] -= linfac*val[m3dc1_field::OP_DRRP];
      v[FIO_DPHI_R] -= linfac*val[m3dc1_field::OP_DRPP];
      v[FIO_DZ_R  ] -= linfac*val[m3dc1_field::OP_DRZP];

      v[FIO_DR_Z  ] -= linfac*val[m3dc1_field::OP_DRZP];
      v[FIO_DPHI_Z] -= linfac*val[m3dc1_field::OP_DZPP];
      v[FIO_DZ_Z  ] -= linfac*val[m3dc1_field::OP_DZZP];
    }
  }

  v[FIO_DR_R    ] *= source->B0;
  v[FIO_DR_PHI  ] *= source->B0;
  v[FIO_DR_Z    ] *= source->B0;
  v[FIO_DPHI_R  ] *= source->B0;
  v[FIO_DPHI_PHI] *= source->B0;
  v[FIO_DPHI_Z  ] *= source->B0;
  v[FIO_DZ_R    ] *= source->B0;
  v[FIO_DZ_PHI  ] *= source->B0;
  v[FIO_DZ_Z    ] *= source->B0;
  
  return FIO_SUCCESS;
}


int m3dc1_current_density::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  i1 = source->file.load_field("I", time);
  if(!i1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
    i0 = source->file.load_field("I", -1);
    if(!i0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    ix = source->file.load_field("I_ext", time);
    if(!ix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_current_density::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_DDVAL);

  const m3dc1_field::m3dc1_get_op gget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op fget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_PPVAL);

  double val[m3dc1_field::OP_NUM];
  double r = (source->itor==1) ? x[0] : 1.;

  // J_R   = -(d(F+f'')/dZ)/R + (d(psi')/dR)/R^2
  // J_Z   =  (d(F+f'')/dR)/R + (d(psi')/dZ)/R^2
  // J_Phi = -Del*[psi]/R

  if(!psi1->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;

  v[0] =  linfac*val[m3dc1_field::OP_DRP]/(r*r);
  v[2] =  linfac*val[m3dc1_field::OP_DZP]/(r*r);
  v[1] = -linfac*(val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ])/r;
  if(source->itor==1) v[1] += linfac*val[m3dc1_field::OP_DR]/(r*r); 

  if(!i1->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;
  v[0] -= linfac*val[m3dc1_field::OP_DZ]/r;
  v[2] += linfac*val[m3dc1_field::OP_DR]/r;

  if(use_f) {
    if(!f1->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= linfac*val[m3dc1_field::OP_DZPP]/r;
    v[2] += linfac*val[m3dc1_field::OP_DRPP]/r;
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] +=  val[m3dc1_field::OP_DRP]/(r*r);
    v[2] +=  val[m3dc1_field::OP_DZP]/(r*r);
    v[1] -= (val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ])/r;
    if(source->itor==1) v[1] += val[m3dc1_field::OP_DR]/(r*r); 

    if(!i0->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[0] -= val[m3dc1_field::OP_DZ]/r;
    v[2] += val[m3dc1_field::OP_DR]/r;
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] += linfac*val[m3dc1_field::OP_DRP]/(r*r);
    v[2] += linfac*val[m3dc1_field::OP_DZP]/(r*r);
    v[1] -= linfac*(val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ])/r;
    if(source->itor==1) v[1] += linfac*val[m3dc1_field::OP_DR]/(r*r); 

    if(!ix->eval(x[0], x[1]-phase, x[2], gget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;
    v[0] -= linfac*val[m3dc1_field::OP_DZ]/r;
    v[2] += linfac*val[m3dc1_field::OP_DR]/r;

    if(use_f) {
      if(!fx->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
        return FIO_OUT_OF_BOUNDS;

      v[0] -= linfac*val[m3dc1_field::OP_DZPP]/r;
      v[2] += linfac*val[m3dc1_field::OP_DRPP]/r;
    }
  }

  v[0] *= source->J0;
  v[1] *= source->J0;
  v[2] *= source->J0;
  
  return FIO_SUCCESS;
}



int m3dc1_grad_vector_potential::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  i1 = source->file.load_field("I", time);
  if(!i1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
    i0 = source->file.load_field("I", -1);
    if(!i0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    ix = source->file.load_field("I_ext", time);
    if(!ix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_grad_vector_potential::eval(const double* x, double* v, void* s)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL | m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);

  const m3dc1_field::m3dc1_get_op fget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_DDVAL);

  double val[m3dc1_field::OP_NUM];

  if(!psi1->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
    return FIO_OUT_OF_BOUNDS;

  v[1] = linfac*(val[m3dc1_field::OP_DR] - val[m3dc1_field::OP_1]/x[0])/x[0];
  v[4] = linfac*val[m3dc1_field::OP_DP]/x[0];
  v[7] = linfac*val[m3dc1_field::OP_DZ]/x[0];

  if(use_f) {
    if(!f1->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[0] =  linfac*(x[0]*val[m3dc1_field::OP_DRZ] + val[m3dc1_field::OP_DZ]);
    v[2] = -linfac*(x[0]*val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DR]);
    v[3] =  linfac*(x[0]*val[m3dc1_field::OP_DZP]);
    v[5] = -linfac*(x[0]*val[m3dc1_field::OP_DRP]);
    v[6] =  linfac*(x[0]*val[m3dc1_field::OP_DZZ]);
    v[8] = -linfac*(x[0]*val[m3dc1_field::OP_DRZ]);
  } else {
    v[0] = v[2] = v[3] = v[5] = v[6] = v[8] = 0.;
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[1] += (val[m3dc1_field::OP_DR]-val[m3dc1_field::OP_1]/x[0])/x[0];
    v[4] += val[m3dc1_field::OP_DP]/x[0];
    v[7] += val[m3dc1_field::OP_DZ]/x[0];

    v[2] -= source->bzero*source->rzero/x[0];
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1]-phase, x[2], psiget, val, (int*)s))
      return FIO_OUT_OF_BOUNDS;

    v[1] += linfac*(val[m3dc1_field::OP_DR]-val[m3dc1_field::OP_1]/x[0])/x[0];
    v[4] += linfac*val[m3dc1_field::OP_DP]/x[0];
    v[7] += linfac*val[m3dc1_field::OP_DZ]/x[0];

    if(use_f) {
      if(!fx->eval(x[0], x[1]-phase, x[2], fget, val, (int*)s))
        return FIO_OUT_OF_BOUNDS;

      v[0] += linfac*(x[0]*val[m3dc1_field::OP_DRZ] + val[m3dc1_field::OP_DZ]);
      v[2] -= linfac*(x[0]*val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DR]);
      v[3] += linfac*(x[0]*val[m3dc1_field::OP_DZP]);
      v[5] -= linfac*(x[0]*val[m3dc1_field::OP_DRP]);
      v[6] += linfac*(x[0]*val[m3dc1_field::OP_DZZ]);
      v[8] -= linfac*(x[0]*val[m3dc1_field::OP_DRZ]);
    }
  }

  // convert to mks
  v[0] *= source->B0;
  v[1] *= source->B0;
  v[2] *= source->B0;
  v[3] *= source->B0*source->L0;
  v[4] *= source->B0*source->L0;
  v[5] *= source->B0*source->L0;
  v[6] *= source->B0;
  v[7] *= source->B0;
  v[8] *= source->B0;
  
  return FIO_SUCCESS;
}
