#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <compound_field.h>

#include <iostream>


int main(int argc, char* argv[])
{
  int result;
  const int timeslice = 1;
  fio_source* src;
  fio_field *pressure, *density, *electron_temperature;
  fio_field *psin0, *mag0;
  fio_option_list opt;

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
    result = fio_open_source(&src, FIO_M3DC1_SOURCE, "../examples/data/m3dc1/C1.h5");
    
  } else if(source_type == "geqdsk") {
    result = fio_open_source(&src, FIO_GEQDSK_SOURCE, "../examples/data/geqdsk/g158115.04701");

  } else if(source_type == "gpec") {
    result = fio_open_source(&src, FIO_GPEC_SOURCE, "../examples/data/gpec");

  } else {
    std::cerr << "Error: source type " << argv[1]
	      << " not recognized" << std::endl;
    return 1;
  };


  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening file" << std::endl;
    delete(src);
    return result;
  };


  // set options for fields obtained from this source
  src->get_field_options(&opt);
  opt.set_option(FIO_TIMESLICE, timeslice);
  opt.set_option(FIO_PART, FIO_TOTAL);
  opt.set_option(FIO_LINEAR_SCALE, 5.);

  // open fields
  result = src->get_field(FIO_TOTAL_PRESSURE, &pressure, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening pressure field" << std::endl;
    pressure = 0;
  };

  opt.set_option(FIO_SPECIES, FIO_ELECTRON);
  result = src->get_field(FIO_DENSITY, &density, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening density field" << std::endl;
    density = 0;
  };
  result = src->get_field(FIO_TEMPERATURE, &electron_temperature, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron temperature field" << std::endl;
    density = 0;
  };
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
  

  void* h;
  src->allocate_search_hint(&h);

  double x[3];

  /*
  int npts = 10;
  double R0 = 1.6;
  double R1 = 2.1;
  double Z0 = 0.0;
  double Z1 = 0.0;
  double phi0 = 0.;
  double phi1 = 0.;

  double p, n, te;
    
  for(int i=0; i<npts; i++) {
    x[0] = R0 + (R1-R0)*i/(npts-1);
    x[1] = phi0 + (phi1-phi0)*i/(npts-1);
    x[2] = Z0 + (Z1-Z0)*i/(npts-1);
    
    std::cout << "(" << x[0] << ", " << x[1] << ", " << x[2] << "):\n";

    if(pressure) {
      result = pressure->eval(x, &p);
      std::cout << "\tpressure = " << p << "\n";
    }

    if(density) {
      result = density->eval(x, &n);
      std::cout << "\tdensity = " << n << "\n";
    }

    if(electron_temperature) {
      result = electron_temperature->eval(x, &te);
      std::cout << "\telectron temperature = " << te << "\n";
    }
  }
  */

  double axis[3];
  axis[0] = 1.7;
  axis[1] = 0.;
  axis[2] = 0.;

  /*
  x[0] = 2.17526;
  x[1] = 0.;
  x[2] = -0.312184;
  result = fio_find_val_2d(electron_temperature, 1578.1, x, 0.1, 0.1, axis);
  if(result==FIO_SUCCESS) {
    std::cerr << "Found at " << x[0] << ", " << x[1] << ", " << x[2] 
	      << std::endl;
  }

  return 0;
  */

  x[0] = 1.8;
  x[1] = 0.;
  x[2] = 0.4;

  // calcualte 2D profiles
  //  int npsi = 10;    // number of psi points for profiles
  //result = fio_find_val_2d(electron_temperature, 374.928, x, 1., 0.01, axis);


  int nphi = 1;
  int ntheta = 400;
  double** path;
  path = new double*[3];
  path[0] = new double[nphi*ntheta];
  path[1] = new double[nphi*ntheta];
  path[2] = new double[nphi*ntheta];

  result = fio_gridded_isosurface(psin0, 0.5, x,
				  axis, 1e-4, 0.01,
				  nphi, ntheta, &path, h);

  if(result==FIO_SUCCESS) {
    std::ofstream file;
    file.open("surface_2d.dat", std::ofstream::out|std::ofstream::trunc);
    for(int i=0; i<ntheta; i++) {
      file << path[1][i] << ", " << path[0][i] << ", " << path[2][i] 
	   << '\n';
    }
    file.close();

    double q;
    result = fio_q_at_surface(mag0, ntheta, path, &q, 0, h);

    std::cerr << "q = " << q << std::endl;
  }
  delete[] path[0];
  delete[] path[1];
  delete[] path[2];
  delete[] path;


  // Calculate 3D surfaces
  nphi = 16;
  path = new double*[3];
  path[0] = new double[nphi*ntheta];
  path[1] = new double[nphi*ntheta];
  path[2] = new double[nphi*ntheta];
  double** path_plane = new double*[3];
  double* bpol = new double[nphi*ntheta];

  int nsurf = 10;
  double* q = new double[nsurf];
  double* psi_norm_surf = new double[nsurf];


  for(int s=0; s<nsurf; s++) {
    std::cerr << "Surface " << s << std::endl;
    psi_norm_surf[s] = (s+1.)/(nsurf+1.);
    std::cerr << "  psi_norm = " << psi_norm_surf[s] << std::endl;

    x[0] = 1.8;
    x[1] = 0.;
    x[2] = 0.;

    result = fio_find_val_2d(psin0, psi_norm_surf[s], x, 
			     1e-4, 0.1, axis, h);
    if(result != FIO_SUCCESS) {
      std::cerr << "Error finding psi_norm = " << psi_norm_surf[s] << std::endl;
      break;
    }
    std::cerr << "  found at " << x[0] << ", " <<  x[1] << ", " << x[2]
	      << std::endl;

    // Find electron temperature on surface
    double te;
    result = electron_temperature->eval(x, &te, h);
    if(result != FIO_SUCCESS) {
      std::cerr << "Error evaluating electron temperature" << std::endl;
      break;
    }
    std::cerr << "  Te = " << te << std::endl;

    // Calculate 3D surface
    result = fio_gridded_isosurface(electron_temperature, te, x,
				    axis, 1., 0.1,
				    nphi, ntheta, &path, h);

    if(result != FIO_SUCCESS) {
      std::cerr << "Error finding surface" << std::endl;
      //break;
    }

    if(result==FIO_SUCCESS) {
      // Estimate q by averaging q over each toroidal plane
      double q_plane;
      q[s] = 0.;
      for(int i=0; i<nphi; i++) {
	path_plane[0] = &(path[0][i*ntheta]);
	path_plane[1] = &(path[1][i*ntheta]);
	path_plane[2] = &(path[2][i*ntheta]);
	result = fio_q_at_surface(mag0, ntheta, path_plane, &q_plane, 
				  &(bpol[i*ntheta]), h);
	//	std::cerr << "q_plane = " << q_plane << std::endl;
	q[s] += q_plane;
      }
      q[s] /= nphi;
      std::cerr <<  " q = " <<  q[s] << std::endl;

      std::ofstream file;
      std::string surface_filename;

      surface_filename = "surface_" + std::to_string(s) + ".dat";
      std::cerr << " writing surface data to " << surface_filename 
		<< std::endl;
      file.open(surface_filename, std::ofstream::out|std::ofstream::trunc);
      int k=0;
      for(int j=0; j<nphi; j++) {
	if(j > 0) file << '\n';
	for(int i=0; i<ntheta; i++) {
	  file << path[1][k] << ", " << path[0][k] << ", " << path[2][k] 
	       << '\n';;
	  k++;
	}
      }
      file.close();
      

    }

  }

  delete[] psi_norm_surf;
  delete[] bpol;
  delete[] path_plane;
  delete[] path[0];
  delete[] path[1];
  delete[] path[2];
  delete[] path;


  // Deallocate fields
  src->deallocate_search_hint(&h);

  fio_close_field(&pressure);
  fio_close_field(&density);
  fio_close_field(&electron_temperature);
  fio_close_field(&psin0);
  fio_close_field(&mag0);
  fio_close_source(&src);

  return 0;
}



