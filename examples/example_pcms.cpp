#include <fusion_io.h>
#include "pcms_lib.h"

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: example <source_type> <isRdv>\n"
      << "where <source_type> is one of m3dc1, geqdsk, gpec\n" 
      << "where <isRdv> is 0 or 1\n" << std::endl;
    return 1;
  }

  auto mpi = Omega_h::Library(&argc, &argv);
  PCMS_Library lib(argc, argv, mpi.world()->get_impl());

  int result;
  fio_source *src;
  fio_field *pressure, *density, *magnetic_field;
  fio_option_list opt;

  std::string source_type(argv[1]);
  std::cout << "source_type = " << source_type << std::endl;

  if(source_type == "m3dc1") {
    // Open an m3dc1 source
    result = fio_open_source(&src, FIO_M3DC1_SOURCE, "data/m3dc1/C1.h5", lib);
      
  } else if(source_type == "geqdsk") {
    result = fio_open_source(&src, FIO_GEQDSK_SOURCE, "data/geqdsk/g158115.04701", lib);

  } else if(source_type == "gpec") {
    result = fio_open_source(&src, FIO_GPEC_SOURCE, "data/gpec", lib);

  } else {
    std::cerr << "Error: source type " << argv[1]
      << " not recognized" << std::endl;
    return 1;
  };

  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening file" << std::endl;
    exit(1);
  };

  // set options for fields obtained from this source
  src->get_field_options(&opt);
  opt.set_option(FIO_TIMESLICE, 1);
  opt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
  opt.set_option(FIO_LINEAR_SCALE, 100.);
  opt.set_option(FIO_PHASE, 180.);

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
      auto value = pressure->evaluate(x, &p);
      printf("\tpressure = %.5e\n", value[0]);
    }

    if(density) {
      auto value = density->evaluate(x, &n);
      printf("\tdensity = %.5e\n", value[0]);
    }

    if(magnetic_field) {
      auto value = magnetic_field->evaluate(x, b);
      printf("\tB = (%.5e, %.5e, %.5e):\n", value[0], value[1], value[2]);
    }
  }

  fio_close_field(&pressure);
  fio_close_field(&density);
  fio_close_field(&magnetic_field);
  fio_close_source(&src);

  return 0;
}