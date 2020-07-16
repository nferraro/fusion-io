#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <options.h>
#include <compound_field.h>

#include <iostream>

int main(int argc, char* argv[])
{
  int result;
  fio_source* src;
  fio_field *pressure, *density, *magnetic_field;
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
    result = fio_open_source(&src, FIO_M3DC1_SOURCE, "data/m3dc1/C1.h5");
    
  } else if(source_type == "geqdsk") {
    result = fio_open_source(&src, FIO_GEQDSK_SOURCE, "data/geqdsk/g158115.04701");

  } else if(source_type == "gpec") {
    result = fio_open_source(&src, FIO_GPEC_SOURCE, "data/gpec");

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
  opt.set_option(FIO_TIMESLICE, 1);
  opt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
  opt.set_option(FIO_LINEAR_SCALE, 10.);

  // open fields
  result = src->get_field(FIO_MAGNETIC_FIELD, &magnetic_field, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening magnetic field" << std::endl;
    magnetic_field = 0;
  };

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
  
  int npts = 10;
  double R0 = 1.6;
  double R1 = 2.1;
  double Z0 = 0.0;
  double Z1 = 0.0;
  double phi0 = 0.;
  double phi1 = 0.;
  double x[3];
  double p, n, b[3];
  
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

    if(magnetic_field) {
      result = magnetic_field->eval(x, b);
      std::cout << "\tB = (" << b[0] << ", " << b[1] << ", " << b[2] << "):\n";
    }
  }

  fio_close_field(&pressure);
  fio_close_field(&density);
  fio_close_field(&magnetic_field);
  fio_close_source(&src);

  return 0;
}


