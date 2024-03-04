#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <compound_field.h>

#include <iostream>
#include <cstring>

#include <netcdf.h>

int npsi = 50;   // number of psi points for profiles
int nr = 10;      // number of radial points for surfaces
int ntheta = 400; // number of poloidal points
int nphi = 32;    // number of toroidal points
double dl_tor = 0.01;     // step size when finding surfaces
double dl_pol = 0.01;     // step size when finding surfaces
double max_step = 0.02;   // Maximum step size for Newton iterations
double tol = 0.1;         // Tolerance for Te when finding isosurface
double dR0 = 0.0;         // Guess for offset from magnetic axis
double psi_start = -1;
double psi_end = -1;
double te_start = -1;
double te_end = -1;
std::deque<fio_source*> sources;
fio_compound_field electron_density, electron_temperature;
fio_compound_field ion_density, ion_temperature, psin, mag;
fio_field *psin0, *te0, *te_prof;
double axis[3];
double psi0, psi1;

bool pert_prof = true; // Use perturbed surfaces for calculating profiles

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
  std::cerr << "Radial profile points (npsi) = " << npsi << '\n';
  std::cerr << "First value of psi_norm (psi_start) = "
	    << psi_start << '\n';
  std::cerr << "Last value of psi_norm (psi_end) = "
	    << psi_end << '\n';
  std::cerr << "Use perturbed surfaces for profiles (pert_prof)? "
	    << pert_prof << '\n';
  std::cerr << "Tolerance for finding isosurface (tol) = " << tol << " eV\n";
  std::cerr << "=======================" << std::endl;

  // Find magnetic axis
  axis[1] = 0.;

  fio_source* src0 = sources[0];

  if(!pert_prof) {
    fio_option_list opt;
    src0->get_field_options(&opt);

    opt.set_option(FIO_PART, FIO_EQUILIBRIUM_ONLY);
    result = src0->get_field(FIO_POLOIDAL_FLUX_NORM, &psin0, &opt);
    if(result != FIO_SUCCESS) {
      std::cerr << "Error opening equilirbium psi norm field" << std::endl;
      psin0 = 0;
    }
    opt.set_option(FIO_SPECIES, FIO_ELECTRON);
    result = src0->get_field(FIO_TEMPERATURE, &te0, &opt);
    if(result != FIO_SUCCESS) {
      std::cerr << "Error opening equilibrium Te field" << std::endl;
      te0 = 0;
    }
  }

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
  double* bpol = new double[nphi*ntheta];
  //  double* btor = new double[nr*nphi*ntheta];
  double** b = new double*[3];
  double** b_plane = new double*[3];
  b[0] = new double[nr*nphi*ntheta];
  b[1] = new double[nr*nphi*ntheta];
  b[2] = new double[nr*nphi*ntheta];

  double* q = new double[nr];
  double* psi_surf = new double[nr];
  double* psi_prof = new double[npsi];
  double* ne = new double[npsi];
  double* te = new double[npsi];
  double* ni = new double[npsi];
  double* ti = new double[npsi];
  double** axis_3d = new double*[nphi];
  for(int i=0; i<nphi; i++)
    axis_3d[i] = new double[3];

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

  x[0] = axis[0] + dR0 - 0.01;
  x[1] = 0.;
  x[2] = 0.;

  result = fio_find_max(&electron_temperature, &te_max, x,
			1., 0.1, 2, n, h);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error finding magnetic axis " << std::endl;
  } else {
    std::cerr << "Found magnetic axis" << std::endl;
    std::cerr << "TE MAX = " << te_max << " AT ("
	    << x[0] << ", " << x[1] << ", " << x[2] << std::endl;

    axis[0] = x[0];
    axis[1] = 0;
    axis[2] = x[2];
  }


  for(int surfaces=0; surfaces<2; surfaces++) {
  
    int s_max;
    double s_start, s_end;
    double psi_norm;
    
    if(surfaces==1) {
      std::cerr << "Calculating surfaces..." << std::endl;
      s_max = nr;
      s_start = psi_start;
      s_end = psi_end;
    } else {
      std::cerr << "Calculating profiles..." << std::endl; 
      s_max = npsi;
      s_start = 1./npsi;
      s_end = 1. - 1./npsi;
    }

    if(surfaces==0 && !pert_prof) {
      te_prof = te0;
    } else {
      te_prof = &electron_temperature;
    }

    for(int s=0; s<s_max; s++) {
      
      std::cerr << "Surface " << s+1 << " of " << s_max << std::endl;

      x[0] = axis[0] - 0.01;
      x[1] = 0.;
      x[2] = axis[2];

      if((te_end>0.) && (te_start>0.)) {
	// use Te as radial coordinate
	double te0 = (te_end - te_start)*s/(s_max - 1.) + te_start;
	std::cerr << "Searching for Te = " << te0 << " at "
		  << "( " << x[0] << ", " << x[1] << ", " << x[2] << " )"
		  << std::endl;
	result = fio_find_val(&electron_temperature, te0, x, 1., 0.1, 1, axis, h);
	if(result != FIO_SUCCESS) {
	  std::cerr << "Error finding Te = " << te0  << std::endl;
	  break;
	} else {
	  std::cerr << "Found Te = " << te0 << std::endl;
	}
	psi_norm = (te_max-te0)/te_max;

      } else {
	// use psi_norm as radial coordinate
	if(s_max==1) {
	  psi_norm = s_start;
	} else {
	  psi_norm = (s_end-s_start)*s/(s_max-1.) + s_start;
	}

	if(surfaces==0 && !pert_prof) {
	  result = fio_find_val(psin0, psi_norm, x, 1e-4, 0.1, 1, axis, h);
	} else {
	  result = fio_find_val(&psin, psi_norm, x, 1e-4, 0.1, 1, axis, h);
	}
	if(result != FIO_SUCCESS) {
	  std::cerr << "Error finding psi_norm = " << psi_norm
		    << std::endl;
	  break;
	}

      }

      // Find electron temperature on surface
      double temp;
      result = te_prof->eval(x, &temp, h);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error evaluating electron temperature" << std::endl;
	break;
      }
      
      char* label;
      char surf_label[256];
      if(surfaces==1) {
	sprintf(surf_label, "%d", s);
	label = surf_label;
      } else {
	label = 0;
      }
      if(surfaces==1) {
	path_surf[0] = &(path[0][s*nphi*ntheta]);
	path_surf[1] = &(path[1][s*nphi*ntheta]);
	path_surf[2] = &(path[2][s*nphi*ntheta]);
      } else {
	path_surf[0] = path[0];
	path_surf[1] = path[1];
	path_surf[2] = path[2];
      }
      result = fio_gridded_isosurface(te_prof, temp, x,
				      axis_3d, dl_tor, dl_pol, tol, max_step,
				      nphi, ntheta, 
				      phi, theta,
				      path_surf, label, h);

      if(result!=FIO_SUCCESS) {
	std::cerr << "Error finding surface at Te = " << temp << std::endl;
	continue;
      }


      // Estimate q by averaging q over each toroidal plane
      double q_plane;
      if(surfaces==1) {
        q[s] = 0.;
      } else {
        ne[s] = 0.;
      }
      for(int i=0; i<nphi; i++) {
        path_plane[0] = &(path_surf[0][i*ntheta]);
        path_plane[1] = &(path_surf[1][i*ntheta]);
        path_plane[2] = &(path_surf[2][i*ntheta]);
	if(surfaces==1) {
	  b_plane[0] = &(b[0][s*nphi*ntheta+i*ntheta]);
	  b_plane[1] = &(b[1][s*nphi*ntheta+i*ntheta]);
	  b_plane[2] = &(b[2][s*nphi*ntheta+i*ntheta]);
	  result = fio_q_at_surface(&mag, ntheta, path_plane, &q_plane,
				    &(bpol[i*ntheta]), b_plane, h);
	} else {
	  result = fio_q_at_surface(&mag, ntheta, path_plane, &q_plane,
				    &(bpol[i*ntheta]), NULL, h);
	}
        //      std::cerr << "q_plane = " << q_plane << std::endl;              
        if(surfaces==1) {
          q[s] += q_plane;
        }
      }

      if(surfaces==1) {
	// Calculate q
        psi_surf[s] = psi_norm*(psi1 - psi0) + psi0;
        q[s] /= nphi;
	std::cerr <<  " q = " <<  q[s] << std::endl;                    
      } else {

	// Calculate flux surface averages for profiles
        psi_prof[s] = psi_norm*(psi1 - psi0) + psi0;
        te[s] = temp;
        result = fio_surface_average(&electron_density, nphi*ntheta, path_surf,
                                     &(ne[s]), bpol, h);
        result = fio_surface_average(&ion_temperature, nphi*ntheta, path_surf,
                                     &(ti[s]), bpol, h);
        result = fio_surface_average(&ion_density, nphi*ntheta, path_surf,
                                     &(ni[s]), bpol, h);
      }

      if(surfaces==1) {
	gplot << "'surface_" << std::to_string(s) << ".dat' u 2:3 w l";
	if(s < s_max-1)
	  gplot << ", \\\n";
      }
    }

  }

  gplot << std::endl;  
  gplot.close();
  splot.close();


  // Create netcdf file
  std::cerr << "Writing netcdf file" << std::endl;

  int ncid;
  result = nc_create("neo_input.nc", NC_CLOBBER, &ncid);
  if(result != 0) {
    std::cerr << "Error opening netcdf file\n" << nc_strerror(result)
	      << std::endl;
  }

  // swap indices
  double* R = new double[nr*nphi*ntheta];
  double* Z = new double[nr*nphi*ntheta];
  double* BR   = new double[nr*nphi*ntheta];
  double* BPhi = new double[nr*nphi*ntheta];
  double* BZ   = new double[nr*nphi*ntheta];

  for(int i=0; i<nr; i++) {
    for(int j=0; j<nphi; j++) {
      for(int k=0; k<ntheta; k++) {
	R[i + j*nr + k*nr*nphi] = path[0][k + j*ntheta + i*ntheta*nphi];
	Z[i + j*nr + k*nr*nphi] = path[2][k + j*ntheta + i*ntheta*nphi];
	BR  [i + j*nr + k*nr*nphi] = b[0][k + j*ntheta + i*ntheta*nphi];
	BPhi[i + j*nr + k*nr*nphi] = b[1][k + j*ntheta + i*ntheta*nphi];
	BZ  [i + j*nr + k*nr*nphi] = b[2][k + j*ntheta + i*ntheta*nphi];
      }
    }
  }

  // Calculate Jacobian
  double* Jac = new double[nr*nphi*ntheta];
  double* dPsids = new double[nr*nphi];
  
  for(int i=0; i<nr; i++) {
    for(int j=0; j<nphi; j++) {
      dPsids[i+j*nr] = 0.;
      for(int k=0; k<ntheta; k++) {
	double dRdi, dZdj, dRdj, dZdi;
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

	dRdj = (R[i + j*nr + kp*nr*nphi] - R[i + j*nr + km*nr*nphi])/2.;
	dZdj = (Z[i + j*nr + kp*nr*nphi] - Z[i + j*nr + km*nr*nphi])/2.;

	Jac[ijk] = dRdi*dZdj - dRdj*dZdi;

	// note that btor is stored in [i][j][k] whereas Jac is in [k][j][i]
	dPsids[i+j*nr] = dPsids[i+j*nr] + Jac[ijk]*b[1][k + j*ntheta + i*ntheta*nphi];
      }
    }
  }
 
  double ion_mass = 2.;  // TODO: generalize this
  int version = 3;

  // version 3
  // * Includes calculation of B_tor and Jacobian

  // convert phi to degrees
  for(int i=0; i<nphi; i++) { 
    phi[i] = 180.*phi[i]/M_PI;
    //    std::cerr << "phi[" << i << "] = " << phi[i] << std::endl;
  }

  // sanity check of toroidal flux
  double *psi_t = new double[nr];
  double *dpsi_t = new double[nr];
  for(int j=0; j<nphi; j++) {
    double flux = 0.;
    for(int i=0; i<nr; i++) {
      flux = flux + dPsids[i + j*nr];
      if(j==0) {
	psi_t[i] = flux;
	dpsi_t[i] = dPsids[i + j*nr];
	std::cerr << "Flux inside surface " << i << " at phi=0: " << psi_t[i] << std::endl;
      } else {
	if(abs(flux - psi_t[i]) > 0.01*abs(psi_t[i])) {
	  std::cerr << "Warning: toroidal flux at " << phi[j] << ": " << flux
		    << "         toroidal flux at " << phi[0] << ": " << psi_t[i]
		    << std::endl;
	}
      }
    }
    std::cerr << "Total toroidal flux at phi = " << phi[j]
	      << ": " << flux << std::endl;
  }

  
  nc_put_att_int(ncid, NC_GLOBAL, "version", NC_SHORT, 1, &version);
  nc_put_att_double(ncid, NC_GLOBAL, "ion_mass", NC_FLOAT, 1, &ion_mass);
  nc_put_att_double(ncid, NC_GLOBAL, "psi_0", NC_FLOAT, 1, &psi0);
  nc_put_att_double(ncid, NC_GLOBAL, "psi_1", NC_FLOAT, 1, &psi1);

  int nr_dimid, np_dimid, nt_dimid, npsi_dimid;
  nc_def_dim(ncid, "npsi", npsi, &npsi_dimid);
  nc_def_dim(ncid, "nr", nr,     &nr_dimid);
  nc_def_dim(ncid, "np", ntheta, &np_dimid);
  nc_def_dim(ncid, "nt", nphi,   &nt_dimid);

  int phi_id;
  nc_def_var(ncid, "Phi", NC_FLOAT, 1, &nt_dimid, &phi_id);

  int q_id, psi_id, psi_t_id, dpsi_t_id;
  nc_def_var(ncid, "q", NC_FLOAT, 1, &nr_dimid, &q_id);
  nc_def_var(ncid, "psi", NC_FLOAT, 1, &nr_dimid, &psi_id);
  nc_def_var(ncid, "Psi_t", NC_FLOAT, 3, &nr_dimid, &psi_t_id);   // added in version 3
  nc_def_var(ncid, "dPsi_t", NC_FLOAT, 3, &nr_dimid, &dpsi_t_id);   // added in version 3

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

  nc_put_var_double(ncid, phi_id, phi);
  nc_put_var_double(ncid, q_id, q);
  nc_put_var_double(ncid, psi_id, psi_surf);
  nc_put_var_double(ncid, psi0_id, psi_prof);
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

  nc_close(ncid);

  delete[] R;
  delete[] Z;
  delete[] BR;
  delete[] BPhi;
  delete[] BZ;
  delete[] Jac;
  delete[] dPsids;

  delete[] theta;
  delete[] phi;
  delete[] q;
  delete[] psi_surf;
  delete[] psi_t;
  delete[] dpsi_t;
  delete[] ne;
  delete[] te;
  delete[] ni;
  delete[] ti;
  delete[] psi_prof;
  delete[] bpol;
  //  delete[] btor;
  delete[] b[0];
  delete[] b[1];
  delete[] b[2];
  delete[] b;
  delete[] b_plane;
  delete[] path_plane;
  delete[] path_surf;
  delete[] path[0];
  delete[] path[1];
  delete[] path[2];
  delete[] path;
  for(int i=0; i<nphi; i++)
    delete[] axis_3d[i];
  delete[] axis_3d;

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
  if(!pert_prof) {
    delete(te0);
    delete(psin0);
  }
  delete_sources();

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
  const int num_opts = 14;
  std::string arg_list[num_opts] = 
    { "-dR0", "-m3dc1", "-max_step", "-nphi", "-nphi", "-npsi", "-nr",
      "-ntheta", "-pert_prof", "-psi_end", "-psi_start", "-te_start", "-te_end", "-tol" };
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

  if(opt=="-dR0") {
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
  } else if(opt=="-npsi") {
    if(argc==1) npsi = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-nr") {
    if(argc==1) nr = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-ntheta") {
    if(argc==1) ntheta = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-pert_prof") {
    if(argc==1) pert_prof = (atoi(argv[0].c_str()) != 0);
    else argc_err = true;
  } else if(opt=="-psi_end") {
    if(argc==1) psi_end = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-psi_start") {
    if(argc==1) psi_start = atof(argv[0].c_str());
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
	    << " -dR0 <dR0>"
	    << " -m3dc1 <m3dc1_source> <time> <scale> <phase>"
	    << " -nphi <nphi>"
	    << " -npsi <npsi>"
	    << " -nr <nr>"
	    << " -ntheta <ntheta>"
	    << " -psi_end <psi_end>"
	    << " -psi_start <psi_start>"
	    << std::endl;

  std::cerr
    << "<dR0>:          offset to major radius\n"
    << "<m3dc1_source>: filename of M3D-C1 source file\n"
    << "<nphi>:         number of toroidal points per surface\n"
    << "<npsi>:         number of radial points for T, n profiles\n"
    << "<nr>:           number of surfaces\n"
    << "<ntheta>:       number of poloidal points per surface\n"
    << "<pert_prof>:    if 1, use perturbed surfaces for T, n profiles\n"
    << "                if 0, use unperturbed surfaces\n"
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
