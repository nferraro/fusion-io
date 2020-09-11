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

  x[0] = 1.8;
  x[1] = 0.;
  x[2] = 0.;
  /*
  result = fio_find_val_2d(electron_temperature, 3900., x, 1e-2, 0.1);
  if(result==FIO_SUCCESS) {
    std::cerr << "Found at " << x[0] << ", " << x[1] << ", " << x[2] 
	      << std::endl;
	      }
  */
  double enclose[3];
  enclose[0] = 1.7;
  enclose[1] = 0.;
  enclose[2] = 0.;
  double** path0;
  int m0;
  const int nphi = 16;
  result = fio_isosurface(electron_temperature, 3900., x,
			  enclose, 1e-2, 1e-2,
			  0.1, 2.*M_PI/nphi, &m0, &path0);

  if(result == FIO_SUCCESS) {
    std::cerr << "Found path with " << m0 << " points." << std::endl;


    std::ofstream file;

    file.open("surface_raw.dat", std::ofstream::out | std::ofstream::trunc);
    for(int i=0; i<m0; i++) {
      if(i > 0) {
	if(path0[1][i] != path0[1][i-1])
	  file << '\n';
      }
      file << path0[1][i] << ", " << path0[0][i] << ", " << path0[2][i] 
	   << '\n';;
    }
    file.close();
    
    
    const int ntheta = 100;
    double* theta = new double[ntheta];
    double* phi = new double[nphi];
    double** path = new double*[3];
    path[0] = new double[ntheta*nphi];
    path[1] = new double[ntheta*nphi];
    path[2] = new double[ntheta*nphi];
    result = 
      fio_gridify_surface(m0, path0, enclose, nphi, ntheta, path, phi, theta);

    file.open("surface_gridded.dat", std::ofstream::out|std::ofstream::trunc);
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
    
    std::cerr << "Deallocating.." << std::endl;
    delete[] path0[0];
    delete[] path0[1];
    delete[] path0[2];
    delete[] path0;

    delete[] theta;
    delete[] path[0];
    delete[] path[1];
    delete[] path[2];
    delete[] path;
  }

  fio_close_field(&pressure);
  fio_close_field(&density);
  fio_close_field(&electron_temperature);
  fio_close_source(&src);

  return 0;
}



