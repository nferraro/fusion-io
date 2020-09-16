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
double psi_start = -1;
double psi_end = -1;
double scalefac = 1.;
fio_source* src;

bool parse_args(int argc, char* argv[]);

int main(int argc, char* argv[])
{
  int result;
  const int timeslice = 1;
  fio_field *electron_density, *electron_temperature;
  fio_field *ion_density, *ion_temperature;
  fio_field *psin0, *mag0;
  fio_option_list opt;

  /*
  if(argc < 2) {
    std::cerr << "Usage: example <source_type>\n"
	      << " where <source_type> is one of \n"
	      << " m3dc1, geqdsk, gpec" << std::endl;
    return 1;
  }


  std::string source_type(argv[1]);
  std::cerr << source_type << std::endl;

  if( source_type == "gpec" ) {
    std::cerr << "source_type = gpec" << std::endl;
  }
  std::cerr << source_type << std::endl;

  if(source_type == "m3dc1") {
    // Open an m3dc1 source
    */
    /*
  } else if(source_type == "geqdsk") {
    result = fio_open_source(&src, FIO_GEQDSK_SOURCE, "../examples/data/geqdsk/g158115.04701");

  } else if(source_type == "gpec") {
    result = fio_open_source(&src, FIO_GPEC_SOURCE, "../examples/data/gpec");

  } else {
    std::cerr << "Error: source type " << argv[1]
	      << " not recognized" << std::endl;
    return 1;
  };
    */
  
  if(!parse_args(argc, argv))
    return 1;

  if(psi_start <= 0. || psi_start > 1.) psi_start = 0.;
  if(psi_end <= 0. || psi_end > 1.) psi_end = 1.;
  if(psi_start <= 0.) psi_start = (psi_end-psi_start)/nphi;
  if(psi_end==1.) psi_end = 1. - (psi_end-psi_start)/nphi;

  std::cerr << "Input parameters\n=======================\n";
  std::cerr << "Scale factor (scale) = " << scalefac << '\n';
  std::cerr << "Radial grid points (nr) = " << nr << '\n';
  std::cerr << "Toroidal grid points (nphi) = " << nphi << '\n';
  std::cerr << "Poloidal grid points (ntheta) = " << ntheta << '\n';
  std::cerr << "Radial profile points (npsi) = " << npsi << '\n';
  std::cerr << "First value of psi_norm (psi_start) = "
	    << psi_start << '\n';
  std::cerr << "Last value of psi_norm (psi_end) = "
	    << psi_end << '\n';
  std::cerr << "=======================" << std::endl;

  double slice_time = 0;
  src->get_slice_time(timeslice, &slice_time);
  
  std::cerr << "Slice time = " << slice_time << std::endl;

  // Find magnetic axis
  fio_series* magaxis[2];
  double axis[3];
  axis[1] = 0.;

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
  if((result = magaxis[0]->eval(slice_time, &axis[0])) != FIO_SUCCESS)
    {
      std::cerr << "Couldn't evaluate MAGAXIS_R" << std::endl;
      delete(src);
      return result;
    }
  if((result = magaxis[1]->eval(slice_time, &axis[2])) != FIO_SUCCESS)
    {
      std::cerr << "Couldn't evaluate MAGAXIS_Z" << std::endl;
      delete(src);
      return result;
    }
  delete(magaxis[0]);
  delete(magaxis[1]);

  std::cerr << "Magnetic axis found at (" 
	    << axis[0] << ", " << axis[2] << ")." << std::endl;


  // Find magnetic axis
  fio_series *magaxis_psi, *lcfs_psi;
  double psi0, psi1;

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


  // set options for fields obtained from this source
  src->get_field_options(&opt);
  opt.set_option(FIO_TIMESLICE, timeslice);
  opt.set_option(FIO_PART, FIO_TOTAL);
  opt.set_option(FIO_LINEAR_SCALE, scalefac);

  // open fields
  result = src->get_field(FIO_MAGNETIC_FIELD, &mag0, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening magnetic field field" << std::endl;
    mag0 = 0;
  };
  result = src->get_field(FIO_POLOIDAL_FLUX_NORM, &psin0, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening psi norm field" << std::endl;
    psin0 = 0;
  };


  opt.set_option(FIO_SPECIES, FIO_ELECTRON);
  result = src->get_field(FIO_DENSITY, &electron_density, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening density field" << std::endl;
    electron_density = 0;
  };
  result = src->get_field(FIO_TEMPERATURE, &electron_temperature, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron temperature field" << std::endl;
    electron_temperature = 0;
  };
  opt.set_option(FIO_SPECIES, FIO_MAIN_ION);
  result = src->get_field(FIO_DENSITY, &ion_density, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening density field" << std::endl;
    ion_density = 0;
  };
  result = src->get_field(FIO_TEMPERATURE, &ion_temperature, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron temperature field" << std::endl;
    ion_temperature = 0;
  };
  

  void* h;
  src->allocate_search_hint(&h);

  double x[3];

  // Calculate 3D surfaces
  double** path = new double*[3];
  double* theta = new double[ntheta];
  double* phi = new double[nphi];
  path[0] = new double[nr*nphi*ntheta];
  path[1] = new double[nr*nphi*ntheta];
  path[2] = new double[nr*nphi*ntheta];
  double** path_plane = new double*[3];  // contains path of one toroidal plane
  double** path_surf = new double*[3];   // contains path of one surface
  double* bpol = new double[nphi*ntheta];

  double* q = new double[nr];
  double* psi_surf = new double[nr];
  double* psi_prof = new double[npsi];
  double* ne = new double[npsi];
  double* te = new double[npsi];
  double* ni = new double[npsi];
  double* ti = new double[npsi];
  /*
  for(int i=0; i<ntheta; i++) 
    theta[i] = 2.*M_PI*i/ntheta;
  for(int i=0; i<nphi; i++)
    phi[i] = 360.*i/nphi;
  */

  // Surfaces
  std::ofstream gplot, splot;
  gplot.open("plot_surfaces", std::ofstream::out|std::ofstream::trunc);
  splot.open("splot_surfaces", std::ofstream::out|std::ofstream::trunc);
  
  gplot << "plot ";
  splot << "set hidden3d\nsplot ";

  for(int surfaces=0; surfaces<2; surfaces++) {
    int s_max;
    double psi_norm;
    
    if(surfaces==1) {
      std::cerr << "Calculating surfaces..." << std::endl;
      s_max = nr;
    } else {
      std::cerr << "Calculating profiles..." << std::endl; 
      s_max = npsi;
    }
  
    x[0] = axis[0] + 0.01;

    for(int s=0; s<s_max; s++) {
      
      if(s_max==1) {
	psi_norm = psi_start;
      } else {
	psi_norm = (psi_end-psi_start)*s/(s_max-1.) + psi_start;
      }
      
      x[1] = axis[1];
      x[2] = axis[2];

      result = fio_find_val_2d(psin0, psi_norm, x, 1e-4, 0.1, axis, h);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error finding psi_norm = " << psi_norm
		  << std::endl;
	std::cerr << " psi_norm = " << psi_norm << std::endl;
	break;
      }
      
      // Find electron temperature on surface
      double temp;
      result = electron_temperature->eval(x, &temp, h);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error evaluating electron temperature" << std::endl;
	break;
      }
      
      // Calculate 3D surface
      int ind;
      if(surfaces==1) {
	ind = s*nphi*ntheta;
      } else {
	ind = 0;
      }
      path_surf[0] = &(path[0][ind]);
      path_surf[1] = &(path[1][ind]);
      path_surf[2] = &(path[2][ind]);

      result = fio_gridded_isosurface(electron_temperature, temp, x,
				      axis, 0.01, 1., 0.01,
				      nphi, ntheta, phi, theta,
				      path_surf, h);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error constructing surface at psi_norm = "
		  << psi_norm << std::endl;
      }

      /*
      result = fio_gridded_isosurface(electron_temperature, temp, x,
				      axis, 1., 0.1,
				      nphi, ntheta, &path_surf, h);
      */
	/*
      double** temp_path;
      int temp_n;
      //      std::cerr << "finding isosurface..." << std::endl;
      result = fio_isosurface(electron_temperature, temp, x,
			      axis, 0.01, 1., 0.1, 
			      2.*M_PI/nphi, &temp_n, &temp_path, h);

      if(result != FIO_SUCCESS) {
	std::cerr << "Error finding surface" << std::endl;
	std::cerr << " psi_norm = " << psi_norm << std::endl;
	break;
      }
      //      std::cerr << "Found isosurface with " << temp_n  << " points" << std::endl;

      //      std::cerr << "Gridifying isosurface..." << std::endl;
      result = fio_gridify_surface(temp_n, temp_path, axis, 
				   nphi, ntheta, 
				   path_surf, phi, theta);
            
      if(result != FIO_SUCCESS) {
	std::cerr << "Error gridifying surface" << std::endl;
	std::cerr << " psi_norm = " << psi_norm << std::endl;
	break;
      }
      //      std::cerr << "Successfully gridified surface" << std::endl;

      delete[] temp_path[0];
      delete[] temp_path[1];
      delete[] temp_path[2];
      delete[] temp_path;
	*/

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
	result = fio_q_at_surface(mag0, ntheta, path_plane, &q_plane, 
				  &(bpol[i*ntheta]), h);
	//	std::cerr << "q_plane = " << q_plane << std::endl;


	if(surfaces==1) {
	  q[s] += q_plane;
	}
      }

      if(surfaces==1) {
	psi_surf[s] = psi_norm*(psi1 - psi0) + psi0;
	q[s] /= nphi;
	//	std::cerr <<  " q = " <<  q[s] << std::endl;
      } else {
	psi_prof[s] = psi_norm*(psi1 - psi0) + psi0;
	te[s] = temp;
	result = fio_surface_average(electron_density, nphi*ntheta, path_surf, 
				     &(ne[s]), bpol, h);
	result = fio_surface_average(ion_temperature, nphi*ntheta, path_surf, 
				     &(ti[s]), bpol, h);
	result = fio_surface_average(ion_density, nphi*ntheta, path_surf, 
				     &(ni[s]), bpol, h);
      }

      
      if(surfaces==1) {
	// Write surface data
	
	if(result==FIO_SUCCESS) {	  
	  std::ofstream file;
	  std::string surface_filename;
	  
	  surface_filename = "surface_" + std::to_string(s) + ".dat";
	  file.open(surface_filename, std::ofstream::out|std::ofstream::trunc);
	  int k=0;
	  for(int j=0; j<nphi; j++) {
	    if(j > 0) file << '\n';
	    for(int i=0; i<ntheta; i++) {
	      file << path_surf[1][k] << ", " 
		   << path_surf[0][k] << ", " 
		   << path_surf[2][k] << '\n';;
	      k++;
	    }
	  }
	  file.close();
	  
	  gplot << " '" << surface_filename << "' u 2:3 w l";
	  splot << " '" << surface_filename << "' w l";
	  if(s<nr-1) {
	    gplot << ", \\\n";
	    splot << ", \\\n";
	  }
	}
      }
    }
  }

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

  double ion_mass = 2.;  // TODO: generalize this
  int version = 2;

  // convert phi to degrees
  for(int i=0; i<nphi; i++)
    phi[i] = 180.*phi[i]/M_PI;


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

  int q_id, psi_id;
  nc_def_var(ncid, "q", NC_FLOAT, 1, &nr_dimid, &q_id);
  nc_def_var(ncid, "psi", NC_FLOAT, 1, &nr_dimid, &psi_id);

  int ne_id, te_id, ni_id, ti_id, psi0_id;
  nc_def_var(ncid, "psi0", NC_FLOAT, 1, &npsi_dimid, &psi0_id);
  nc_def_var(ncid, "ne0",  NC_FLOAT, 1, &npsi_dimid, &ne_id);
  nc_def_var(ncid, "Te0",  NC_FLOAT, 1, &npsi_dimid, &te_id);
  nc_def_var(ncid, "ni0",  NC_FLOAT, 1, &npsi_dimid, &ni_id);
  nc_def_var(ncid, "Ti0",  NC_FLOAT, 1, &npsi_dimid, &ti_id);

    int r_id, z_id;
  int dims[3];
  dims[0] = nr_dimid;
  dims[1] = nt_dimid;
  dims[2] = np_dimid;
  nc_def_var(ncid, "R", NC_FLOAT, 3, dims, &r_id);
  nc_def_var(ncid, "Z", NC_FLOAT, 3, dims, &z_id);

  nc_enddef(ncid);

  nc_put_var_double(ncid, phi_id, phi);
  nc_put_var_double(ncid, q_id, q);
  nc_put_var_double(ncid, psi_id, psi_surf);
  nc_put_var_double(ncid, psi0_id, psi_prof);
  nc_put_var_double(ncid, ne_id, ne);
  nc_put_var_double(ncid, te_id, te);
  nc_put_var_double(ncid, ni_id, ni);
  nc_put_var_double(ncid, ti_id, ti);
  nc_put_var_double(ncid, r_id, path[0]);
  nc_put_var_double(ncid, z_id, path[2]);

  nc_close(ncid);

  delete[] theta;
  delete[] phi;
  delete[] q;
  delete[] psi_surf;
  delete[] ne;
  delete[] te;
  delete[] ni;
  delete[] ti;
  delete[] psi_prof;
  delete[] bpol;
  delete[] path_plane;
  delete[] path_surf;
  delete[] path[0];
  delete[] path[1];
  delete[] path[2];
  delete[] path;


  // Deallocate fields
  src->deallocate_search_hint(&h);

  fio_close_field(&electron_density);
  fio_close_field(&electron_temperature);
  fio_close_field(&ion_density);
  fio_close_field(&ion_temperature);
  fio_close_field(&psin0);
  fio_close_field(&mag0);
  fio_close_source(&src);

  return 0;
}



bool parse_args(int argc, char* argv[])
{
  for(int i=1; i<argc; i++) {
    if(strcmp(argv[i],"-psi_start")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -psi_start requires an argument" << std::endl;
	return false;
      }
      psi_start = atof(argv[i+1]);
    }
    if(strcmp(argv[i],"-psi_end")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -psi_end requires an argument" << std::endl;
	return false;
      }
      psi_end = atof(argv[i+1]);
    }
    if(strcmp(argv[i],"-scale")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -scale requires an argument" << std::endl;
	return false;
      }
      scalefac = atof(argv[i+1]);
    }
    if(strcmp(argv[i],"-npsi")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -npsi requires an argument" << std::endl;
	return false;
      }
      npsi = atoi(argv[i+1]);
    }
    if(strcmp(argv[i],"-nr")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -nr requires an argument" << std::endl;
	return false;
      }
      nr = atoi(argv[i+1]);
    }
    if(strcmp(argv[i],"-ntheta")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -ntheta requires an argument" << std::endl;
	return false;
      }
      ntheta = atoi(argv[i+1]);
    }
    if(strcmp(argv[i],"-nphi")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -nphi requires an argument" << std::endl;
	return false;
      }
      nphi = atoi(argv[i+1]);
    }
    if(strcmp(argv[i],"-m3dc1")==0) {
      if(i+1 >= argc) {
	std::cerr << "Error: -m3dc1 requires an argument" << std::endl;
	return false;
      }
      int result = fio_open_source(&src, FIO_M3DC1_SOURCE, argv[i+1]);
      if(result != FIO_SUCCESS) {
	std::cerr << "Error opening m3dc1 file " << argv[1];
	return false;
      }
    }
  }

  return true;
}
