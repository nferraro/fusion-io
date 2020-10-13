#include "trace_integrator.h"

#include <iostream>
#include <iomanip>
#include <math.h>

trace_integrator::trace_integrator()
{
  reverse = false;
  toroidal = false;
  period = 2.*M_PI;
  plane = 0.;
  nplanes = 1;
  tpts = 1;
}

trace_integrator::~trace_integrator()
{
  close_file();
}


bool trace_integrator::eval(const double r, const double phi, const double z, 
			    double* br, double* bphi, double* bz)
{
  double m[3], x[3];

  x[0] = r;
  x[1] = phi;
  x[2] = z;

  *br = 0.;
  *bphi = 0.;
  *bz = 0.;
      
  trace_source_list::iterator i = sources.begin();

  while(i != sources.end()) {
    int result = i->field->eval(x, m, i->hint);
    if(result != FIO_SUCCESS)
      return false;
    *br   += m[0];
    *bphi += m[1];
    *bz   += m[2];
    i++;
  }

  return true;
}

bool trace_integrator::eval_psi(const double r, const double phi, const double z, 
			    double* psi)
{
  if(!sources[0].psi_norm)
    return false;
  
  double x[3];
  x[0] = r;
  x[1] = phi;
  x[2] = z;

  int result = sources[0].psi_norm->eval(x, psi, sources[0].hint);
  if(result != FIO_SUCCESS)
    return false;

  return true;
}

bool trace_integrator::open_file(const char* filename)
{
  if(file.is_open()) close_file();

  file.open(filename, std::fstream::out | std::fstream::trunc);

  if(!file) return false;
  return true;
}

bool trace_integrator::close_file()
{
  file.close();

  return true;
}

void trace_integrator::set_reverse(const bool r)
{
  reverse = r;
}

bool trace_integrator::set_pos(const double r,const double phi,const double z)
{
  R = r;
  Phi = phi;
  Z = z;

  return true;
}

bool trace_integrator::center(double* R0, double* Z0) const
{
  if(sources.size()==0)
    return false;

  *R0 = sources[0].magaxis[0];
  *Z0 = sources[0].magaxis[1];

  return true;
}


bool trace_integrator::get_surface(const double r0, const double phi0, 
				   const double z0, const double ds, 
				   double** r, double** z, int* n)
{
  if(sources.size()==0)
    return false;

  return true;
  //  return sources[0]->get_surface(r0, phi0, z0, ds, r, z, n);
}


bool trace_integrator::integrate(int transits, int steps_per_transit, 
				 integrator_data* data)
{
  double dphi;
  bool plot, result;
  int i, k=0;
  int steps = transits*steps_per_transit + 1;
  bool ptrans, add = false;
  double avg_steps_per_pol_transit;
  double steps_since_pol_transit;
  double pd = (reverse ? -period : period);

  if(data) {
    data->toroidal_transits = 0;
    data->poloidal_transits = 0;
    data->q = 0;
    data->distance = 0.;
    avg_steps_per_pol_transit = 0;
    steps_since_pol_transit = 0;
  }

  double R0, Z0;
  center(&R0,&Z0);

  dphi = pd/(double)steps_per_transit;

  for(i=0; i<steps; i++) {

    double next_Phi = Phi + dphi;
    if(reverse) {
      if(next_Phi <= pd) {
	next_Phi -= pd;
	Phi -= pd;
      }
    } else {
      if(next_Phi >= pd) {
	next_Phi -= pd;
	Phi -= pd;
      }
    }

    double last_R = R;
    double last_Phi = Phi;
    double last_Z = Z;
    double pl;

    // if Phi will pass through the plotting plane on this step, plot intercept
    if(nplanes <=1) {
      if(reverse) {
	plot = (Phi > plane && next_Phi <= plane);
      } else  {
	plot = (Phi < plane && next_Phi >= plane);
      }
      pl = plane;
    } else {
      for(int j=0; j<nplanes; j++) {
	pl = pd*(double)j/(double)nplanes + plane;
	if(reverse) {
	  plot = (Phi > pl && next_Phi <= pl);
	} else {
	  plot = (Phi < pl && next_Phi >= pl);
	}
	if(plot) break;
      }
    }

      
    // shift derivatives at previous timesteps for predictor-correctors
    for(int j=3; j>0; j--) {
      dr[j] = dr[j-1];
      dz[j] = dz[j-1];
    }

    // take timestep
    //result = step_euler(dphi);
    result = step_rk4(dphi);
    //    result = step_rk3(dphi);
    //result = (i < 4) ? step_rk4(dphi) : step_predcorr(dphi);

    if(!result) return false;

    if(data) {
      if(toroidal) {
	data->distance += 
	  sqrt((R-last_R)*(R-last_R) + (Z-last_Z)*(Z-last_Z)
	       + R*last_R*dphi*dphi);
      } else {
	data->distance += 
	  sqrt((R-last_R)*(R-last_R) + (Z-last_Z)*(Z-last_Z)
	       + dphi*dphi);
      }

      // count toroidal transits
      k++;
      if(k==steps_per_transit) {
	k = 0;
	data->toroidal_transits++;
      }

      // count poloidal transits
      ptrans = false;
      if     (last_Z < Z0 && Z >= Z0) ptrans = true;
      else if(last_Z > Z0 && Z <= Z0) ptrans = true;
      if(ptrans) {
	if(add) {
	  data->poloidal_transits++;
	  data->q = (double)data->toroidal_transits/
	    (double)data->poloidal_transits;
	  avg_steps_per_pol_transit = 
	    (double)(avg_steps_per_pol_transit * (data->poloidal_transits - 1.)
	     + steps_since_pol_transit) / (double)(data->poloidal_transits);
	  steps_since_pol_transit = 0;
	}
	add = !add;
      } else {
	steps_since_pol_transit++;
      }
    }
    if(plot) {
      double f = (pl - last_Phi)/dphi;

      double R_plot = R*f + last_R*(1.-f);
      double Z_plot = Z*f + last_Z*(1.-f);
      double theta_plot = atan2(Z_plot - Z0, R_plot - R0)*180./M_PI;
      double psi_plot;
      if(!eval_psi(R_plot, plane, Z_plot, &psi_plot))
	psi_plot = 0;

      file << std::setiosflags(std::ios::scientific)
	   << std::setw(20) << std::setprecision(12) 
	   << pl*180./M_PI
	   << std::setiosflags(std::ios::scientific)
	   << std::setw(20) << std::setprecision(12) 
	   << R_plot 
	   << std::setiosflags(std::ios::scientific) 
	   << std::setw(20) << std::setprecision(12)
	   << Z_plot
	   << std::setiosflags(std::ios::scientific)
	   << std::setw(20) << std::setprecision(12)
	   << theta_plot
	   << std::setiosflags(std::ios::scientific) 
	   << std::setw(20) << std::setprecision(12)
	   << psi_plot << std::endl;
    }
  }

  if(data) {
    double denom = data->poloidal_transits 
      + (double)steps_since_pol_transit/avg_steps_per_pol_transit;
    data->q = (double)data->toroidal_transits/denom;
  }
    
  return true;
}

bool trace_integrator::step_euler(double dphi)
{
  double b_r, b_phi, b_z;

  if(!eval(R,Phi,Z,&b_r,&b_phi,&b_z))
    return false;

  Phi += dphi;
  dr[0] = R*b_r/b_phi;
  dz[0] = R*b_z/b_phi;

  R += dphi*dr[0];
  Z += dphi*dz[0];

  return true;
}

bool trace_integrator::step_rk3(double dphi)
{
  double b_r, b_phi, b_z, dR, dZ;

  if(!eval(R,Phi,Z,&b_r,&b_phi,&b_z))
    return false;

  Phi += dphi/2.;
  double k1_R = R*dphi*b_r/b_phi;
  double k1_Z = R*dphi*b_z/b_phi;
   
  if(!eval(R + k1_R/2.,Phi,Z + k1_Z/2.,&b_r,&b_phi,&b_z))
    return false;

  Phi += dphi/2.;
  double k2_R = (R+k1_R/2.)*dphi*b_r/b_phi;
  double k2_Z = (R+k1_R/2.)*dphi*b_z/b_phi;
  
  if(!eval(R - k1_R + 2.*k2_R,Phi,Z - k1_Z + 2.*k2_Z,&b_r,&b_phi,&b_z))
    return false;

  double k3_R = (R-k1_R+2.*k2_R)*dphi*b_r/b_phi;
  double k3_Z = (R-k1_R+2.*k2_R)*dphi*b_z/b_phi;

  dR = k1_R/6. + 2.*k2_R/3. + k3_R/6.;
  dZ = k1_Z/6. + 2.*k2_Z/3. + k3_Z/6.;

  R += dR;
  Z += dZ;

  // store the derivative at the current step
  // for use with predictor-corrector methods
  dr[0] = dR/dphi;
  dz[0] = dZ/dphi;

  return true;
}


bool trace_integrator::step_rk4(double dphi)
{
  double b_r, b_phi, b_z, dR, dZ;
  double RR = toroidal ? R : 1.;

  if(!eval(R,Phi,Z,&b_r,&b_phi,&b_z))
    return false;

  Phi += dphi/2.;

  double k1_R = RR*dphi*b_r/b_phi;
  double k1_Z = RR*dphi*b_z/b_phi;
   
  if(!eval(R + k1_R/2.,Phi,Z + k1_Z/2.,&b_r,&b_phi,&b_z))
    return false;

  RR = toroidal ? (R+k1_R/2.) : 1.;
  double k2_R = RR*dphi*b_r/b_phi;
  double k2_Z = RR*dphi*b_z/b_phi;

  if(!eval(R + k2_R/2.,Phi,Z + k2_Z/2.,&b_r,&b_phi,&b_z))
    return false;

  Phi += dphi/2.;
  RR = toroidal ? (R+k2_R/2.) : 1.;
  double k3_R = RR*dphi*b_r/b_phi;
  double k3_Z = RR*dphi*b_z/b_phi;
  
  if(!eval(R + k3_R,Phi,Z + k3_Z,&b_r,&b_phi,&b_z))
    return false;

  RR = toroidal ? (R+k3_R) : 1.;
  double k4_R = RR*dphi*b_r/b_phi;
  double k4_Z = RR*dphi*b_z/b_phi;

  dR = k1_R/6. + k2_R/3. + k3_R/3. + k4_R/6.;
  dZ = k1_Z/6. + k2_Z/3. + k3_Z/3. + k4_Z/6.;

  R += dR;
  Z += dZ;

  // store the derivative at the current step
  // for use with predictor-corrector methods
  dr[0] = dR/dphi;
  dz[0] = dZ/dphi;

  return true;
}

bool trace_integrator::step_predcorr(double dphi)
{
  double b_r, b_phi, b_z;
  double h12 = dphi/12.;

  // Adams-Bashforth predictor step
  double dR = h12*(23.*dr[1] - 16.*dr[2] + 5.*dr[3]);
  double dZ = h12*(23.*dz[1] - 16.*dz[2] + 5.*dz[3]);
  Phi += dphi;

  // Derivative evaulation
  if(!eval(R+dR,Phi,Z+dZ,&b_r,&b_phi,&b_z))
    return false;

  dr[0] = (R+dR)*b_r/b_phi;
  dz[0] = (R+dR)*b_z/b_phi;

  // Adams-Moulton corrector step
  R += h12*(5.*dr[0] + 8.*dr[1] - dr[2]);
  Z += h12*(5.*dz[0] + 8.*dz[1] - dz[2]);

  return true;
}
