#include "trace_integrator_s.h"

#include <iostream>
#include <iomanip>
#include <math.h>

trace_integrator_s::trace_integrator_s()
{
  reverse = false;
  toroidal = true;
  period = 2.*M_PI;
  plane = 0.;
  nplanes = 1;
  tpts = 1;
}

trace_integrator_s::~trace_integrator_s()
{
  close_file();
}

bool trace_integrator_s::rhs(const double r, const double phi, const double z, 
			    double* dx, double* dy)
{
  double m[3], x[3], br, bphi, bz;
  double dr[3], dz[3], d, rr, r0, rx, ry, rz, zx, zy, zz;

  x[0] = r;
  x[1] = phi;
  x[2] = z;

  br = 0.;
  bphi = 0.;
  bz = 0.;
  rx = 0.;
  ry = 0.;
  rz = 0.;
  zx = 0.;
  zy = 0.;
  zz = 0.;
  rr = 0.;
      
  trace_source_list::iterator i = sources.begin();
  int result = i->field->eval(x, m, i->hint);
  if(toroidal) {
    result = i->rst->eval(x, &r0, i->hint);
    if(result != FIO_SUCCESS)
      return false;
    rr += r0;
  } else {
    rr = 1.;
  }
  result = i->rst->eval_deriv(x, dr, i->hint);
  if(result != FIO_SUCCESS)
    return false;
  rx += dr[0];
  rz += dr[1];
  ry += dr[2];
  result = i->zst->eval_deriv(x, dz, i->hint);
  if(result != FIO_SUCCESS)
    return false;
  zx += dz[0];
  zz += dz[1];
  zy += dz[2];

  while(i != sources.end()) {
    int result = i->field->eval(x, m, i->hint);
    if(result != FIO_SUCCESS)
      return false;
    br   += m[0];
    bphi += m[1];
    bz   += m[2];
    i++;
  }
  d = rx*zy - ry*zx; 
  /*std::cerr << "Br " << br << std::endl;
  std::cerr << "Bphi " << bphi << std::endl;
  std::cerr << "Bz " << bz << std::endl;
  std::cerr << "R " << rr << std::endl;
  std::cerr << "x " << r << std::endl;
  std::cerr << "y " << z << std::endl;
  std::cerr << "phi " << phi << std::endl;
  //d *= 5; 
  if(fabs(phi)<1e-6) {
    std::cerr << "Br " << br << std::endl;
    std::cerr << "Bphi " << bphi << std::endl;
    std::cerr << "Bz " << bz << std::endl;
    std::cerr << "R " << rr << std::endl;
    std::cerr << "x " << r << std::endl;
    std::cerr << "y " << z << std::endl;
    std::cerr << "phi " << phi << std::endl;
    std::cerr << "Rx " << rx << std::endl;
    std::cerr << "Ry " << ry << std::endl;
    std::cerr << "Rz " << rz << std::endl;
    std::cerr << "Zx " << zx << std::endl;
    std::cerr << "Zy " << zy << std::endl;
    std::cerr << "Zz " << zz << std::endl;
    std::cerr << "D " << d << std::endl;
  }*/
  *dx = rr*(zy*br - ry*bz)/(d*bphi)+(ry*zz-rz*zy)/d;
  *dy = rr*(rx*bz - zx*br)/(d*bphi)+(rz*zx-rx*zz)/d;
  //*dx = rr*br/(bphi);//+(ry*zz-rz*zy)/d;
  //*dy = rr*bz/(bphi);//+(rz*zx-rx*zz)/d;

  return true;
}

bool trace_integrator_s::eval(const double r, const double phi, const double z, 
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

bool trace_integrator_s::eval_psi(const double r, const double phi, const double z, 
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

bool trace_integrator_s::open_file(const char* filename)
{
  if(file.is_open()) close_file();

  file.open(filename, std::fstream::out | std::fstream::trunc);

  if(!file) return false;
  return true;
}

bool trace_integrator_s::close_file()
{
  file.close();

  return true;
}

void trace_integrator_s::set_reverse(const bool r)
{
  reverse = r;
}

bool trace_integrator_s::set_pos(const double r,const double phi,const double z)
{
  R = r;
  Phi = phi;
  Z = z;

  double z0, x[3] = {R, Phi, Z};
  Zp = 0;
  trace_source_list::iterator i = sources.begin();

  int result = i->zst->eval(x, &z0, i->hint);
  if(result != FIO_SUCCESS)
    return false;
  Zp += z0;
  return true;
}

bool trace_integrator_s::center(double* R0, double* Z0) const
{
  if(sources.size()==0)
    return false;

  //*R0 = sources[0].magaxis[0];
  //*Z0 = sources[0].magaxis[1];
  *R0 = 0.;
  *Z0 = 0.;

  return true;
}


bool trace_integrator_s::get_surface(const double r0, const double phi0, 
				   const double z0, const double ds, 
				   double** r, double** z, int* n)
{
  if(sources.size()==0)
    return false;

  return true;
  //  return sources[0]->get_surface(r0, phi0, z0, ds, r, z, n);
}


bool trace_integrator_s::integrate(int transits, int steps_per_transit, 
				 integrator_data* data)
{
  double dphi;
  bool plot, result;
  int i, k=0, iresult;
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
    double last_Zp = Zp;
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

      double zz0, xx[3] = {R, Phi, Z};
      Zp = 0;
      trace_source_list::iterator m = sources.begin();
      iresult = m->zst->eval(xx, &zz0, m->hint);
      if(iresult != FIO_SUCCESS)
        return false;
      Zp += zz0;

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

      trace_source_list::iterator i = sources.begin();
      double Rout = 0.;
      double Zout = 0.;
      double r0;
      double z0;
      double x[3] = {R_plot, pl, Z_plot};

      iresult = i->rst->eval(x, &r0, i->hint);
      if(iresult != FIO_SUCCESS)
        return false;
      Rout += r0;
      iresult = i->zst->eval(x, &z0, i->hint);
      if(iresult != FIO_SUCCESS)
        return false;
      Zout += z0;

      file << std::setiosflags(std::ios::scientific)
	   << std::setw(20) << std::setprecision(12) 
	   << pl*180./M_PI
	   << std::setiosflags(std::ios::scientific)
	   << std::setw(20) << std::setprecision(12) 
	   << Rout 
	   << std::setiosflags(std::ios::scientific) 
	   << std::setw(20) << std::setprecision(12)
	   << Zout
	   << std::setiosflags(std::ios::scientific)
	   << std::setw(20) << std::setprecision(12)
	   << theta_plot
	   << std::setiosflags(std::ios::scientific) 
	   << std::setw(20) << std::setprecision(12)
	   << psi_plot << std::endl;
    }
  }

  if(data) {
    std::cerr << data->poloidal_transits << " pol. transits" << std::endl;
    double denom = data->poloidal_transits 
      + (double)steps_since_pol_transit/avg_steps_per_pol_transit;
    data->q = (double)data->toroidal_transits/denom;
  }
    
  return true;
}

bool trace_integrator_s::step_euler(double dphi)
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

bool trace_integrator_s::step_rk3(double dphi)
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


bool trace_integrator_s::step_rk4(double dphi)
{
  double dx, dy, dR, dZ;
 
  if(!rhs(R,Phi,Z,&dx,&dy))
    return false;

  Phi += dphi/2.;

  double k1_R = dx*dphi;
  double k1_Z = dy*dphi;
   
  if(!rhs(R + k1_R/2.,Phi,Z + k1_Z/2.,&dx,&dy))
    return false;

  double k2_R = dx*dphi;
  double k2_Z = dy*dphi;

  if(!rhs(R + k2_R/2.,Phi,Z + k2_Z/2.,&dx,&dy))
    return false;

  Phi += dphi/2.;
  double k3_R = dx*dphi;
  double k3_Z = dy*dphi;
  
  if(!rhs(R + k3_R,Phi,Z + k3_Z,&dx,&dy))
    return false;

  double k4_R = dx*dphi;
  double k4_Z = dy*dphi;

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

bool trace_integrator_s::step_predcorr(double dphi)
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
