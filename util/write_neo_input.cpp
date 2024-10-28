#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <compound_field.h>

#include <iostream>
#include <cstring>

#include <netcdf.h>

int nr = 10;      // number of radial points for surfaces
int ntheta = 400; // number of poloidal points
int nphi = 32;    // number of toroidal points
double dl_tor = 0.01;     // step size when finding surfaces
double dl_pol = 0.01;     // step size when finding surfaces
double max_step = 0.02;   // Maximum step size for Newton iterations
double tol = 0.1;         // Tolerance for Te when finding isosurface
double dR0 = 0.0;         // Guess for offset from magnetic axis
double R0 = 0.0;          // Guess for R coordinate of magnetic axis
double psi_start = -1;
double psi_end = -1;
double te_start = -1;
double te_end = -1;
std::deque<fio_source*> sources;
fio_compound_field electron_density, electron_temperature;
fio_compound_field ion_density, ion_temperature, psin, mag;
double axis[3];
double psi0, psi1;

void print_usage();
int process_command_line(int argc, char* argv[]);
int process_line(const std::string& opt, const int argc, 
		  const std::string argv[]);
void delete_sources();


int main(int argc, char* argv[])
{
  int result;

  if(argc < 2) {
    print_usage();
    return 1;
  }
  
  if(process_command_line(argc, argv) != FIO_SUCCESS) {
    print_usage();
    return 1;
  }

  if(sources.size()==0) {
    std::cerr << "Error, no source is loaded" << std::endl;
    print_usage();
    return 1;
  }

  if(psi_start <= 0. || psi_start > 1.) psi_start = 0.;
  if(psi_end <= 0. || psi_end > 1.) psi_end = 1.;
  if(psi_start <= 0.) psi_start = (psi_end-psi_start)/nr;
  if(psi_end==1.) psi_end = 1. - (psi_end-psi_start)/nr;

  std::cerr << "Input parameters\n=======================\n";
  std::cerr << "Maximum Newton iteration step (max_step) = " << max_step
	    << '\n';
  std::cerr << "Radial grid points (nr) = " << nr << '\n';
  std::cerr << "Toroidal grid points (nphi) = " << nphi << '\n';
  std::cerr << "Poloidal grid points (ntheta) = " << ntheta << '\n';
  std::cerr << "First value of psi_norm (psi_start) = "
	    << psi_start << '\n';
  std::cerr << "Last value of psi_norm (psi_end) = "
	    << psi_end << '\n';
  std::cerr << "Tolerance for finding isosurface (tol) = " << tol << " eV\n";
  std::cerr << "=======================" << std::endl;

  // Find magnetic axis
  axis[1] = 0.;

  fio_source* src0 = sources[0];

  void* h;
  src0->allocate_search_hint(&h);

  double x[3];


  // Calculate 3D surfaces
  double** path = new double*[3];
  path[0] = new double[nr*nphi*ntheta];
  path[1] = new double[nr*nphi*ntheta];
  path[2] = new double[nr*nphi*ntheta];

  double* theta = new double[ntheta];
  double* phi = new double[nphi];
  double** path_plane = new double*[3];  // contains path of one toroidal plane
  double** path_surf = new double*[3];   // contains path of one surface
  //  double* bpol = new double[nphi*ntheta];
  //  double* btor = new double[nr*nphi*ntheta];
  double** b = new double*[3];
  double** b_surf = new double*[3];
  b[0] = new double[nr*nphi*ntheta];
  b[1] = new double[nr*nphi*ntheta];
  b[2] = new double[nr*nphi*ntheta];

  double* q = new double[nr];
  double* psi_surf = new double[nr];
  double* ne = new double[nr];
  double* te = new double[nr];
  double* ni = new double[nr];
  double* ti = new double[nr];
  double* psi_t = new double[nr];
  double* dpsi_t = new double[nr];
  double* psi_p = new double[nr];
  double* dpsi_p = new double[nr];
  double* axis_3d[3];
  axis_3d[0] = new double[nphi];
  axis_3d[1] = new double[nphi];
  axis_3d[2] = new double[nphi];
  
  // Surfaces
  std::ofstream gplot, splot;
  gplot.open("plot_surfaces", std::ofstream::out|std::ofstream::trunc);
  splot.open("splot_surfaces", std::ofstream::out|std::ofstream::trunc);
  
  gplot << "plot ";
  splot << "set hidden3d\nsplot ";

  // Find magnetic axis
  double n[3];
  double te_max;
  n[0] = 0.;
  n[1] = 1.;
  n[2] = 0.;

  if(R0 > 0.0) {
    x[0] = R0;
  } else {		 
    x[0] = axis[0] + dR0 - 0.01;
  }
  x[1] = 0.;
  x[2] = 0.;

  result = fio_find_max(&electron_temperature, &te_max, x,
			tol, max_step, 2, n, h);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error finding magnetic axis " << std::endl;
    return result;
  } else {
    std::cerr << "Found magnetic axis" << std::endl;
    std::cerr << "TE MAX = " << te_max << " AT ("
	      << x[0] << ", " << x[1] << ", " << x[2] << ") " 
	      << std::endl;

    axis[0] = x[0];
    axis[1] = 0;
    axis[2] = x[2];
  }

  double psi_norm;

  std::cerr << "Calculating surfaces..." << std::endl;

  for(int s=0; s<nr; s++) {
    std::cerr << "Surface " << s+1 << " of " << nr;

    x[0] = axis[0] - 0.01;
    x[1] = 0.;
    x[2] = axis[2];

    if((te_end>0.) && (te_start>0.)) {
      // use Te as radial coordinate
      double te0 = (te_end - te_start)*s/(nr - 1.) + te_start;
      std::cerr << " ( Te = " << te0 << " )" << std::endl;
      result = fio_find_val(&electron_temperature, te0, x, tol, max_step, 1, axis, h);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error finding Te = " << te0  << std::endl;
	return result;
      }
      psi_norm = (te_max-te0)/te_max;

    } else {
      // use psi_norm as radial coordinate
      if(nr==1) {
	psi_norm = psi_start;
      } else {
	psi_norm = (psi_end-psi_start)*s/(nr-1.) + psi_start;
      }
      std::cerr << " ( psi_norm = " << psi_norm << " )"
		<< std::endl;

      result = fio_find_val(&psin, psi_norm, x, 1e-4, 0.1, 1, axis, h);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error finding psi_norm = " << psi_norm
		  << std::endl;
	return result;
      }
    }

    // Find electron temperature on surface
    double temp;
    result = electron_temperature.eval(x, &temp, h);
    if(result != FIO_SUCCESS) {
      std::cerr << "Error evaluating electron temperature" << std::endl;
      break;
    }
    std::cerr << "Found Te = " << temp
	      << ", psi_norm = " << psi_norm << std::endl;

    char* label;
    char surf_label[256];

    sprintf(surf_label, "%d", s);
    label = surf_label;

    path_surf[0] = &(path[0][s*nphi*ntheta]);
    path_surf[1] = &(path[1][s*nphi*ntheta]);
    path_surf[2] = &(path[2][s*nphi*ntheta]);

    axis_3d[0][0] = axis[0];
    axis_3d[1][0] = 0;
    axis_3d[2][0] = axis[2];
    result = fio_gridded_isosurface(&electron_temperature, temp, x,
				    axis_3d, dl_tor, dl_pol, tol, max_step,
				    nphi, ntheta, 
				    phi, theta,
				    path_surf, label, h);

    if(result!=FIO_SUCCESS) {
      std::cerr << "Error finding surface at Te = " << temp << std::endl;
      return result;
    }

    te[s] = temp;
    psi_surf[s] = psi_norm*(psi1 - psi0) + psi0;
    
    gplot << "'surface_" << std::to_string(s) << ".dat' u 2:3 w l";
    if(s < nr-1)
      gplot << ", \\\n";
  }

  gplot << std::endl;  
  gplot.close();
  splot.close();

  // Calculate Jacobian
  std::cerr << "Calculating the Jacobian..." << std::endl;
  double* jac = new double[nr*nphi*ntheta];
  result = fio_isosurface_jacobian(nr, nphi, ntheta, path, jac);
  if(result != FIO_SUCCESS)
    std::cerr << "Error evaluating jacobian" << std::endl;

  // calculate B at each point
  std::cerr << "Calculating B on the surfaces..." << std::endl;
  result = fio_eval_on_path(&mag, nr*nphi*ntheta, path, b, h);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error finding B" << std::endl;
  }

  // Calculate flux surface averages for profiles
  std::cerr << "Calculating flux surface averages..." << std::endl;  
  result = fio_surface_average(&electron_density, nr, nphi, ntheta, path,
			       jac, ne, h);
  result = fio_surface_average(&ion_temperature, nr, nphi, ntheta, path,
			       jac, ti, h);
  result = fio_surface_average(&ion_density, nr, nphi, ntheta, path,
			       jac, ni, h);

  // swap indices
  double* R = new double[nr*nphi*ntheta];
  double* Z = new double[nr*nphi*ntheta];
  double* BR   = new double[nr*nphi*ntheta];
  double* BPhi = new double[nr*nphi*ntheta];
  double* BZ   = new double[nr*nphi*ntheta];
  double* Jac = new double[nr*nphi*ntheta];

  for(int i=0; i<nr; i++) {
    for(int j=0; j<nphi; j++) {
      for(int k=0; k<ntheta; k++) {
	Jac[i + j*nr + k*nr*nphi] = jac[k + j*ntheta + i*ntheta*nphi];
	R[i + j*nr + k*nr*nphi] = path[0][k + j*ntheta + i*ntheta*nphi];
	Z[i + j*nr + k*nr*nphi] = path[2][k + j*ntheta + i*ntheta*nphi];
	BR  [i + j*nr + k*nr*nphi] = b[0][k + j*ntheta + i*ntheta*nphi];
	BPhi[i + j*nr + k*nr*nphi] = b[1][k + j*ntheta + i*ntheta*nphi];
	BZ  [i + j*nr + k*nr*nphi] = b[2][k + j*ntheta + i*ntheta*nphi];
      }
    }
  }


  // Toroidal Flux
  std::cerr << "Calculating toroidal flux" << std::endl;
  double* dPsids = new double[nr*nphi];
  
  for(int i=0; i<nr; i++) {
    for(int j=0; j<nphi; j++) {
      dPsids[i+j*nr] = 0.;
      for(int k=0; k<ntheta; k++) {
	double dRdi, dZdi, dRdk, dZdk;
	int ijk = i + j*nr + k*nr*nphi;

	if(i==0) {
	  dRdi = R[ijk + 1] - R[ijk];
	  dZdi = Z[ijk + 1] - Z[ijk];
	} else if(i==nr-1) {
	  dRdi = R[ijk] - R[ijk-1];
	  dZdi = Z[ijk] - Z[ijk-1];
	} else {
	  dRdi = (R[ijk+1] - R[ijk-1])/2.;
	  dZdi = (Z[ijk+1] - Z[ijk-1])/2.;
	};
	int km = (k==0 ? ntheta-1 : k-1);
	int kp = (k==ntheta-1 ? 0 : k+1);
	dRdk = (R[i + j*nr + kp*nr*nphi] - R[i + j*nr + km*nr*nphi])/2.;
	dZdk = (Z[i + j*nr + kp*nr*nphi] - Z[i + j*nr + km*nr*nphi])/2.;

	// flux through dtheta x dr
	dPsids[i+j*nr] += BPhi[ijk]*(dRdi*dZdk - dRdk*dZdi);
      }
    }
  }

  // calculate toroidal flux
  double fluxt = 0.;
  for(int i=0; i<nr; i++) {
    dpsi_t[i] = 0.;

    // calculate mean toroidal flux
    for(int j=0; j<nphi; j++)
      dpsi_t[i] += dPsids[i+j*nr];
    dpsi_t[i] /= nphi;
    fluxt += dpsi_t[i];
    psi_t[i] = fluxt;

    // calculate variance in toroidal flux
    double fluxt_var = 0.;
    for(int j=0; j<nphi; j++) {
      double fluxt_j = 0.;
      for(int ii=0; ii<=i; ii++)
	fluxt_j += dPsids[ii+j*nr];

      fluxt_var += (fluxt_j - psi_t[i])*(fluxt_j - psi_t[i]);
    }
    fluxt_var /= nphi;

    std::cerr << "Toroidal flux at surface " << i << " = "
	      << psi_t[i] << " +/- " << sqrt(fluxt_var) << std::endl;
  }


  // Poloidal flux
  std::cerr << "Calculating poloidal flux" << std::endl;
  double* dpsids = new double[nr*ntheta];
  double* dpsix  = new double[nr*ntheta];
  double* dAds = new double[nr*ntheta];
  double* dVds = new double[nr*ntheta];

  for(int i=0; i<nr; i++) {
    for(int k=0; k<ntheta; k++) {
      dpsids[i+k*nr] = 0.;
      dpsix[i+k*nr] = 0.;
      dAds[i+k*nr] = 0.;

      // Do integral over phi
      for(int j=0; j<nphi; j++) {
	double dRdi, dZdi, dRdj, dZdj, dRdk, dZdk, dphidj;
	int ijk = i + j*nr + k*nr*nphi;

	if(i==0) {
	  dRdi = R[ijk + 1] - R[ijk];
	  dZdi = Z[ijk + 1] - Z[ijk];
	} else if(i==nr-1) {
	  dRdi = R[ijk] - R[ijk - 1];
	  dZdi = Z[ijk] - Z[ijk - 1];
	} else {
	  dRdi = (R[ijk+1] - R[ijk-1])/2.;
	  dZdi = (Z[ijk+1] - Z[ijk-1])/2.;
	};
	int jm = (j==0 ? nphi-1 : j-1);
	int jp = (j==nphi-1 ? 0 : j+1);
	dRdj = (R[i + jp*nr + k*nr*nphi] - R[i + jm*nr + k*nr*nphi])/2.;
	dZdj = (Z[i + jp*nr + k*nr*nphi] - Z[i + jm*nr + k*nr*nphi])/2.;
	dphidj = 2.*M_PI / nphi;

	int km = (k==0 ? ntheta-1 : k-1);
	int kp = (k==ntheta-1 ? 0 : k+1);
	dRdk = (R[i + j*nr + kp*nr*nphi] - R[i + j*nr + km*nr*nphi])/2.;
	dZdk = (Z[i + j*nr + kp*nr*nphi] - Z[i + j*nr + km*nr*nphi])/2.;

	// B . dphi x dr (poloidal flux)
	dpsids[i + k*nr] -= R[ijk]*dphidj*
	  ( BZ[ijk] * dRdi - BR[ijk] * dZdi) - 
	    BPhi[ijk] * (dRdi*dZdj - dRdj*dZdi);
	
	// B . dtheta x dphi (normal flux)
	dpsix[i + k*nr] -= R[ijk]*dphidj*
	  ( BR[ijk] * dZdk - BZ[ijk] * dRdk) -
	  BPhi[ijk] * (dRdj*dZdk - dRdk*dZdj);

	dAds[i + k*nr] += Jac[ijk] / sqrt(dRdk*dRdk + dZdk*dZdk);
	dVds[i + k*nr] += Jac[ijk];
      }
    }
  }

  // calculate poloidal flux
  double fluxp = 0.;
  for(int i=0; i<nr; i++) {
    dpsi_p[i] = 0.;

    // calculate mean poloidal flux
    for(int k=0; k<ntheta; k++)
      dpsi_p[i] += dpsids[i+k*nr];
    dpsi_p[i] /= ntheta;
    fluxp += dpsi_p[i];   
    psi_p[i] = fluxp;

    // calculate variance in poloidal flux
    double fluxp_var = 0.;
    for(int k=0; k<ntheta; k++) {
      double fluxp_k = 0.;
      for(int ii=0; ii<=i; ii++)
	fluxp_k += dpsids[ii+k*nr];
      
      fluxp_var += (fluxp_k - psi_p[i])*(fluxp_k - psi_p[i]);
    }
    fluxp_var /= ntheta;

    std::cerr << "Poloidal flux at surface " << i << " = "
	      << psi_p[i] << " +/- " << sqrt(fluxp_var) << std::endl;
  }

  // Sanity check: calculate normal flux through surfaces
  // (should be zero)
  for(int i=0; i<nr; i++) {
    double flux_norm = 0;
    for(int k=0; k<ntheta; k++) {
      flux_norm += dpsix[i + k*nr];
    }
    std::cerr << "Net flux through surface " << i << ": "
	      << flux_norm << std::endl;
  }

  // Calculate q
  for(int i=0; i<nr; i++) {
    q[i] = dpsi_t[i] / dpsi_p[i];
    std::cerr << "q at surface " << i << ": "
	      << q[i] << std::endl;
  }


  // Calculate volume inside each surface
  double volume = 0;
  for(int i=0; i<nr; i++) {
    for(int k=0; k<ntheta; k++) {
      volume += dVds[i + k*nr];
    }
 
    std::cerr << "Net volume inside surface " << i << ": "
	      << volume << std::endl;
  }


  double ion_mass = 2.;  // TODO: generalize this
  int version = 4;

  // version 3
  // * Includes calculation of B_tor and Jacobian

  // version 4
  // * Includes 3D axis info

  // convert phi to degrees
  for(int i=0; i<nphi; i++) { 
    phi[i] = 180.*phi[i]/M_PI;
    //    std::cerr << "phi[" << i << "] = " << phi[i] << std::endl;
  }

  // Create netcdf file
  std::cerr << "Writing netcdf file" << std::endl;

  int ncid;
  result = nc_create("neo_input.nc", NC_CLOBBER, &ncid);
  if(result != 0) {
    std::cerr << "Error opening netcdf file\n" << nc_strerror(result)
	      << std::endl;
  }

  //  std::cerr << " attributes..." << std::endl;  
  nc_put_att_int(ncid, NC_GLOBAL, "version", NC_SHORT, 1, &version);
  nc_put_att_double(ncid, NC_GLOBAL, "ion_mass", NC_FLOAT, 1, &ion_mass);
  nc_put_att_double(ncid, NC_GLOBAL, "psi_0", NC_FLOAT, 1, &psi0);
  nc_put_att_double(ncid, NC_GLOBAL, "psi_1", NC_FLOAT, 1, &psi1);

  //  std::cerr << " defining variables..." << std::endl;  
  int nr_dimid, np_dimid, nt_dimid, npsi_dimid;
  nc_def_dim(ncid, "npsi", nr,   &npsi_dimid);
  nc_def_dim(ncid, "nr", nr,     &nr_dimid);
  nc_def_dim(ncid, "np", ntheta, &np_dimid);
  nc_def_dim(ncid, "nt", nphi,   &nt_dimid);

  int phi_id;
  nc_def_var(ncid, "Phi", NC_FLOAT, 1, &nt_dimid, &phi_id);

  int axis_r_id;
  int axis_z_id;
  nc_def_var(ncid, "axis_R", NC_FLOAT, 1, &nt_dimid, &axis_r_id); // ver 4
  nc_def_var(ncid, "axis_Z", NC_FLOAT, 1, &nt_dimid, &axis_z_id); // ver 4
  
  int q_id, psi_id, psi_t_id, dpsi_t_id, psi_p_id, dpsi_p_id;
  nc_def_var(ncid, "q", NC_FLOAT, 1, &nr_dimid, &q_id);
  nc_def_var(ncid, "psi", NC_FLOAT, 1, &nr_dimid, &psi_id);
  nc_def_var(ncid, "Psi_t", NC_FLOAT, 1, &nr_dimid, &psi_t_id);    // ver 3
  nc_def_var(ncid, "dPsi_t", NC_FLOAT, 1, &nr_dimid, &dpsi_t_id);  // ver 3
  nc_def_var(ncid, "Psi_p", NC_FLOAT, 1, &nr_dimid, &psi_p_id);    // ver 4
  nc_def_var(ncid, "dPsi_p", NC_FLOAT, 1, &nr_dimid, &dpsi_p_id);  // ver 4  
  int ne_id, te_id, ni_id, ti_id, psi0_id;
  nc_def_var(ncid, "psi0", NC_FLOAT, 1, &npsi_dimid, &psi0_id);
  nc_def_var(ncid, "ne0",  NC_FLOAT, 1, &npsi_dimid, &ne_id);
  nc_def_var(ncid, "Te0",  NC_FLOAT, 1, &npsi_dimid, &te_id);
  nc_def_var(ncid, "ni0",  NC_FLOAT, 1, &npsi_dimid, &ni_id);
  nc_def_var(ncid, "Ti0",  NC_FLOAT, 1, &npsi_dimid, &ti_id);

  int r_id, z_id, jac_id, br_id, bphi_id, bz_id;
  int dims[3];
  dims[2] = nr_dimid;
  dims[1] = nt_dimid;
  dims[0] = np_dimid;
  nc_def_var(ncid, "R", NC_FLOAT, 3, dims, &r_id);
  nc_def_var(ncid, "Z", NC_FLOAT, 3, dims, &z_id);
  nc_def_var(ncid, "Jac", NC_FLOAT, 3, dims, &jac_id);    // added in version 3
  nc_def_var(ncid, "B_R",   NC_FLOAT, 3, dims, &br_id);   // added in version 3
  nc_def_var(ncid, "B_Phi", NC_FLOAT, 3, dims, &bphi_id); // added in version 3
  nc_def_var(ncid, "B_Z",   NC_FLOAT, 3, dims, &bz_id);   // added in version 3

  nc_enddef(ncid);

  //  std::cerr << " writing variables..." << std::endl;  
  nc_put_var_double(ncid, phi_id, phi);
  nc_put_var_double(ncid, axis_r_id, axis_3d[0]);
  nc_put_var_double(ncid, axis_z_id, axis_3d[2]);
  nc_put_var_double(ncid, q_id, q);
  nc_put_var_double(ncid, psi_id, psi_surf);
  nc_put_var_double(ncid, psi0_id, psi_surf);
  nc_put_var_double(ncid, ne_id, ne);
  nc_put_var_double(ncid, te_id, te);
  nc_put_var_double(ncid, ni_id, ni);
  nc_put_var_double(ncid, ti_id, ti);
  nc_put_var_double(ncid, r_id, R);
  nc_put_var_double(ncid, z_id, Z);
  nc_put_var_double(ncid, jac_id, Jac);
  nc_put_var_double(ncid, br_id, BR);
  nc_put_var_double(ncid, bphi_id, BPhi);
  nc_put_var_double(ncid, bz_id, BZ);
  nc_put_var_double(ncid, psi_t_id, psi_t);
  nc_put_var_double(ncid, dpsi_t_id, dpsi_t);
  nc_put_var_double(ncid, psi_p_id, psi_p);
  nc_put_var_double(ncid, dpsi_p_id, dpsi_p);

  nc_close(ncid);

  std::cerr << "Done writing netcdf file." << std::endl;
  std::cerr << "Freeing memory." << std::endl;

  delete[] R;
  delete[] Z;
  delete[] BR;
  delete[] BPhi;
  delete[] BZ;
  delete[] Jac;
  delete[] dPsids;

  delete[] jac;
  delete[] theta;
  delete[] phi;
  delete[] q;
  delete[] psi_surf;
  delete[] psi_t;
  delete[] dpsi_t;
  delete[] psi_p;
  delete[] dpsi_p;
  delete[] ne;
  delete[] te;
  delete[] ni;
  delete[] ti;
  delete[] b[0];
  delete[] b[1];
  delete[] b[2];
  delete[] b;
  delete[] b_surf;
  delete[] path_plane;
  delete[] path_surf;
  delete[] path[0];
  delete[] path[1];
  delete[] path[2];
  delete[] path;
  delete[] axis_3d[0];
  delete[] axis_3d[1];
  delete[] axis_3d[2];

  // Deallocate fields
  src0->deallocate_search_hint(&h);

  /*
  electron_density.close();
  electron_temperature.close();
  ion_density.close();
  ion_temperature.close();
  psin.close();
  mag.close();
  */
  delete_sources();

  std::cerr << "Done." << std::endl;
  return 0;
}


bool create_source(const int type, const int argc, const std::string argv[]) 
{
  fio_source* src;
  fio_option_list fopt;
  int result;
  double slice_time = 0.;
  int timeslice = 0;

  std::cerr << "Opening source\n";

  switch(type) {
  case(FIO_M3DC1_SOURCE):
    std::cerr << " Type: M3DC1\n";
    src = new m3dc1_source();
    if(argc>=1) {
      result = src->open(argv[0].c_str());
    } else {
      result = src->open("C1.h5");
    }
    if(result != FIO_SUCCESS) {
      std::cerr << "Error opening file" << std::endl;
      delete(src);
      return false;
    };
    // set options for fields obtained from this source
    src->get_field_options(&fopt);
    fopt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
    if(argc<2) fopt.set_option(FIO_TIMESLICE,-1);  // default to timeslice=-1
    break;

  case(FIO_GEQDSK_SOURCE):
    std::cerr << " Type: GEQDSK\n";
    src = new geqdsk_source();
    if(argc>=1) {
      if((result = src->open(argv[0].c_str())) != FIO_SUCCESS) {
	std::cerr << "Error opening file" << std::endl;
	delete(src);
	return false;
      };
    } else {
      std::cerr << "Filename must be provided for geqdsk files." << std::endl;
      delete(src);
      return false;
    }
    src->get_field_options(&fopt);
    break;

  case(FIO_GPEC_SOURCE):
    std::cerr << " Type: GPEC\n";
    src = new gpec_source();
    if(argc>=1) {
      if((result = src->open(argv[0].c_str())) != FIO_SUCCESS) {
	std::cerr << "Error opening file" << std::endl;
	delete(src);
	return false;
      };
    } else {
      std::cerr << "Directory must be provided for gpec files." << std::endl;
      delete(src);
      return false;
    }
    src->get_field_options(&fopt);
    break;

  default:
    return false;
  }

  if(argc>=2) {
    timeslice = atoi(argv[1].c_str());
    std::cerr << " Time slice: " << timeslice << std::endl;
    fopt.set_option(FIO_TIMESLICE, timeslice);
    if(src->get_slice_time(timeslice, &slice_time)==FIO_SUCCESS) {
      std::cerr << "  Slice time = " << slice_time << std::endl;
    }
  }
  if(argc>=3) {
    double scale = atof(argv[2].c_str());
    std::cerr << " Linear scale factor: " << scale << std::endl;
    fopt.set_option(FIO_LINEAR_SCALE, scale);
  }
  if(argc>=4) {
    double phase = atof(argv[3].c_str())*M_PI/180.;
    std::cerr << " Phase: " << phase << std::endl;
    fopt.set_option(FIO_PHASE, phase);
  }

  // Allocate search hint
  fio_hint hint;
  result = src->allocate_search_hint(&hint);
  if(result != FIO_SUCCESS) {
    std::cerr << "Warning: couldn't allocate search hint" << std::endl;
  }



  if(sources.size()==0) {
    // If this is the first source being read, use this source to 
    // determine equilibrium fields and magnetic axis.

    fio_series *magaxis[2];
    fio_series *magaxis_psi, *lcfs_psi;

    // Find magnetic axis
    if((result = src->get_series(FIO_MAGAXIS_R, &magaxis[0])) != FIO_SUCCESS) {
      std::cerr << "Couldn't load MAGAXIS_R" << std::endl;
      delete(src);
      return result;
    }
    if((result = src->get_series(FIO_MAGAXIS_Z, &magaxis[1])) != FIO_SUCCESS) {
      std::cerr << "Couldn't load MAGAXIS_Z" << std::endl;
      delete(src);
      return result;
    }
    if((result = magaxis[0]->eval(slice_time, &(axis[0]))) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate MAGAXIS_R" << std::endl;
	delete(src);
	return result;
    }
    if((result = magaxis[1]->eval(slice_time, &(axis[2]))) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate MAGAXIS_Z" << std::endl;
	delete(src);
	return result;
      }
    delete(magaxis[0]);
    delete(magaxis[1]);


    std::cerr << "Magnetic axis found at (" 
	      << axis[0] << ", " << axis[2] << ")." << std::endl;


    if((result = src->get_series(FIO_MAGAXIS_PSI, &magaxis_psi)) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't load MAGAXIS_PSI" << std::endl;
	delete(src);
	return result;
      }
    if((result = src->get_series(FIO_LCFS_PSI, &lcfs_psi)) != FIO_SUCCESS) 
      {
	std::cerr << "Couldn't load LCFS_PSI" << std::endl;
	delete(src);
	return result;
      }
    if((result = magaxis_psi->eval(slice_time, &psi0)) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate MAGAXIS_PSI" << std::endl;
	delete(src);
	return result;
      }
    if((result = lcfs_psi->eval(slice_time, &psi1)) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate LCFS_PSI" << std::endl;
	delete(src);
	return result;
      }
    delete(magaxis_psi);
    delete(lcfs_psi);
    
    std::cerr << "psi0 = " << psi0 << std::endl;
    std::cerr << "psi1 = " << psi1 << std::endl;


    fopt.set_option(FIO_PART, FIO_TOTAL);
    
  } else {
    fopt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
  }

  // Read fields
  fio_field* field;

  result = src->get_field(FIO_MAGNETIC_FIELD, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening magnetic field field" << std::endl;
    delete(src);
    return result;
  }
  mag.add_field(field, FIO_ADD, 1., hint);
  
  result = src->get_field(FIO_POLOIDAL_FLUX_NORM, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening psi norm field" << std::endl;
    delete(src);
    return result;
  }
  psin.add_field(field, FIO_ADD, 1., hint);


  // open fields
  fopt.set_option(FIO_SPECIES, FIO_ELECTRON);
  result = src->get_field(FIO_TEMPERATURE, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron temperature field" << std::endl;
    delete(src);
    return result;
  }
  electron_temperature.add_field(field, FIO_ADD, 1., hint);

  result = src->get_field(FIO_DENSITY, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron density field" << std::endl;
    delete(src);
    return result;
  }
  electron_density.add_field(field, FIO_ADD, 1., hint);

  fopt.set_option(FIO_SPECIES, FIO_MAIN_ION);
  result = src->get_field(FIO_DENSITY, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening ion density field" << std::endl;
    delete(src);
    return result;
  }
  ion_density.add_field(field, FIO_ADD, 1., hint);

  result = src->get_field(FIO_TEMPERATURE, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening ion temperature field" << std::endl;
    delete(src);
    return result;
  }
  ion_temperature.add_field(field, FIO_ADD, 1., hint);
  

  // Add source to list
  sources.push_back(src);
    

  return FIO_SUCCESS;
}




int process_command_line(int argc, char* argv[])
{
  const int max_args = 4;
  const int num_opts = 15;
  std::string arg_list[num_opts] = 
    { "-dl", "-dR0", "-m3dc1", "-max_step", "-nphi", "-nphi", "-npsi", "-nr",
      "-ntheta", "-psi_end", "-psi_start", "-R0", "-te_start", "-te_end", "-tol" };
  std::string opt = "";
  std::string arg[max_args];
  int args = 0;
  bool is_opt;
  bool processed = true;
  int result;

  for(int i=1; i<argc; i++) {
    // determine if current cl arg is an option
    is_opt = false;
    for(int j=0; j<num_opts; j++) {
      if(argv[i]==arg_list[j]) {
	is_opt = true;
	break;
      }
    }
    
    if(is_opt) {     // if so, process current option
      if(!processed)
	if((result = process_line(opt, args, arg)) != FIO_SUCCESS) 
	  return result;

      opt = argv[i];
      args = 0;
      processed = false;
      
    } else {         // otherwise, add argument
      if(args >= max_args)
	std::cerr << "Too many arguments for " << opt << std::endl;
      else 
	arg[args++] = argv[i];
    }
  }

  if(!processed)
    if((result = process_line(opt, args, arg)) != FIO_SUCCESS)
      return result;

  return FIO_SUCCESS;

}

int process_line(const std::string& opt, const int argc, const std::string argv[])
{
  bool argc_err = false;

  if(opt=="-dl") {
    if(argc==1) dl_pol = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dR0") {
    if(argc==1) dR0 = atof(argv[0].c_str());
    else argc_err = true;    
  } else if(opt=="-m3dc1") {
    return create_source(FIO_M3DC1_SOURCE, argc, argv);
  } else if(opt=="-max_step") {
    if(argc==1) max_step = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-nphi") {
    if(argc==1) nphi = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-nr") {
    if(argc==1) nr = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-ntheta") {
    if(argc==1) ntheta = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-psi_end") {
    if(argc==1) psi_end = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-psi_start") {
    if(argc==1) psi_start = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-R0") {
    if(argc==1) R0 = atof(argv[0].c_str());
    else argc_err = true;    
  } else if(opt=="-te_end") {
    if(argc==1) te_end = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-te_start") {
    if(argc==1) te_start = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-tol") {
    if(argc==1) tol = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-npsi") {
    std::cerr << "Error: option -npsi is deprecated\n"
	      << " profiles are now only calculated at surface locations"
	      << std::endl;
    return FIO_UNSUPPORTED;
  } else if(opt=="-pert_prof") {
    std::cerr << "Error: option -pert_prof is deprecated\n"
	      << " profiles are now always calculated on perturbed surfaces"
	      << std::endl;
    return FIO_UNSUPPORTED;
  } else {
    std::cerr << "Unrecognized option " << opt << std::endl;
    return FIO_UNSUPPORTED;
  }

  if(argc_err) {
    std::cerr << "Incorrect number of arguments for option " 
	      << opt << std::endl;
    return FIO_UNSUPPORTED;
  }

  return FIO_SUCCESS;
}


void print_usage()
{
  std::cerr << "write_neo_input"
	    << " -dl <dl>"
	    << " -dR0 <dR0>"
	    << " -m3dc1 <m3dc1_source> <time> <scale> <phase>"
	    << " -nphi <nphi>"
	    << " -nr <nr>"
	    << " -ntheta <ntheta>"
	    << " -psi_end <psi_end>"
	    << " -psi_start <psi_start>"
	    << std::endl;

  std::cerr
    << "<dl>:           poloidal step size when tracing isosurface\n"
    << "<dR0>:          offset to major radius\n"
    << "<m3dc1_source>: filename of M3D-C1 source file\n"
    << "<nphi>:         number of toroidal points per surface\n"
    << "<nr>:           number of surfaces\n"
    << "<ntheta>:       number of poloidal points per surface\n"
    << "<phase>:        phase factor to apply to fields\n"
    << "<psi_end>:      psi_norm of outermost surface\n"
    << "<psi_start>:    psi_norm of innermost surface\n"
    << "<scale>:        scale factor to apply to linear perturbation\n"
    << "<time>:         timeslice of fields to read\n"
    << "<tol>:          tolerance for finding Te isourface (in eV)\n";
}

void delete_sources()
{
  for(auto it = sources.cbegin(); it != sources.end(); it++) {
    (*it)->close();
    delete(*it);
  }

  sources.clear();
}
