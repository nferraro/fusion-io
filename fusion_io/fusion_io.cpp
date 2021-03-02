#include "fusion_io.h"
#include <iostream>

int fio_open_source(fio_source** src, const int type, const char* filename)
{
  int ierr;
  *src = 0;

  switch(type) {
  case(FIO_GATO_SOURCE):
    *src = new gato_source();
    ierr = (*src)->open(filename);
    break;

  case(FIO_GEQDSK_SOURCE):
    *src = new geqdsk_source();
    ierr = (*src)->open(filename);
    break;

  case(FIO_GPEC_SOURCE):
    *src = new gpec_source();
    ierr = (*src)->open(filename);
    break;

  case(FIO_M3DC1_SOURCE):
    *src = new m3dc1_source();
    ierr = (*src)->open(filename);
    break;

  case(FIO_MARS_SOURCE):
    *src = new mars_source();
    ierr = (*src)->open(filename);
    break;

  default:
    std::cerr << "Source type " << type << " unsupported." << std::endl;
    return FIO_UNSUPPORTED;
  }

  if(ierr != FIO_SUCCESS) {
    delete(*src);
    return ierr;
  }
  
  return FIO_SUCCESS;
}

int fio_close_source(fio_source** source)
{
  if(*source) {
    int result = (*source)->close();
    delete(*source);
    *source = 0;
    return result;
  } else return 1;
}

int fio_close_field(fio_field** field)
{
  if(*field) {
    delete(*field);
    *field = 0;
    return FIO_SUCCESS;
  } else return 1;
}

int fio_get_source_name(const int i, std::string* s)
{
  switch(i) {
  case(FIO_GATO_SOURCE):       *s = "GATO";             break;
  case(FIO_GEQDSK_SOURCE):     *s = "GEQDSK";           break;
  case(FIO_M3DC1_SOURCE):      *s = "M3D-C1";           break;
  case(FIO_NIMROD_SOURCE):     *s = "NIMROD";           break;
  default:
    *s = "Unknown source";
    return FIO_UNSUPPORTED;
  }

  return FIO_SUCCESS;
}

int fio_get_field_name(const field_type f, std::string* s)
{
  switch(f) {
  case(FIO_ALPHA):             *s = "alpha";            break;
  case(FIO_CURRENT_DENSITY):   *s = "current density";  break;
  case(FIO_DENSITY):           *s = "density";          break;
  case(FIO_ELECTRIC_FIELD):    *s = "electric field";   break;
  case(FIO_FLUID_VELOCITY):    *s = "fluid velocity";   break; 
  case(FIO_MAGNETIC_FIELD):    *s = "magnetic field";   break;
  case(FIO_POLOIDAL_FLUX):     *s = "poloidal flux";    break;
  case(FIO_POLOIDAL_FLUX_NORM):*s = "normalized poloidal flux"; break;
  case(FIO_PRESSURE):          *s = "pressure";         break;
  case(FIO_SCALAR_POTENTIAL):  *s = "scalar potential"; break;
  case(FIO_TEMPERATURE):       *s = "temperature";      break;
  case(FIO_TOTAL_PRESSURE):    *s = "total pressure";   break;
  case(FIO_VECTOR_POTENTIAL):  *s = "vector potential"; break;
  case(FIO_VELOCITY):          *s = "velocity";         break; 
  case(FIO_RSTELLA):           *s = "rst";              break;
  case(FIO_ZSTELLA):           *s = "zst";              break;
  case(FIO_JPHI):              *s = "jphi";             break;
  case(FIO_FPRIME):            *s = "fp";               break;
  case(FIO_WALL_DIST):         *s = "wall distance";    break;
  default:
    *s = "Unnamed field";
    return FIO_UNSUPPORTED;
  }

  return FIO_SUCCESS; 
}

int fio_get_option_name(const int opt, std::string* s)
{
  switch(opt) {
  case(FIO_TIMESLICE):    *s = "Time slice";   break;
  case(FIO_PART):         *s = "Part";         break;
  case(FIO_SPECIES):      *s = "Species";      break;
  case(FIO_LINEAR_SCALE): *s = "Linear Scale"; break;
  case(FIO_PHASE):        *s = "Phase";        break;
  default:
    *s = "Unknown Option";
    return FIO_UNSUPPORTED;
  }

  return FIO_SUCCESS; 
}
