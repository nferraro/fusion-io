#ifndef FUSION_IO_SOURCE_H
#define FUSION_IO_SOURCE_H

#include "fusion_io_species.h"
#include "fusion_io_field.h"
#include "options.h"

#include <vector>
#include <math.h>

typedef int field_type;
typedef int series_type;

class fio_source {
 public:
  typedef std::vector<field_type> fio_field_list;

  virtual ~fio_source()
    { }
  virtual int open(const char*) = 0;
  virtual int close() = 0;

  virtual int get_field_options(fio_option_list*) const = 0;
  virtual int get_field(const field_type, fio_field**, const fio_option_list*)
    = 0;
  virtual int get_series(const series_type, fio_series**)
  { return FIO_UNSUPPORTED; }
  virtual int get_available_fields(fio_field_list*) const = 0;
  virtual int get_coordinate_system(int* cs) const
  { 
    *cs = FIO_CYLINDRICAL;
    return FIO_SUCCESS; 
  }
  virtual int get_period(double* p) const
  {
    *p = 2.*M_PI;
    return FIO_SUCCESS;
  }
};

#endif 
