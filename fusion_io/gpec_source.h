#ifndef GPEC_SOURCE_H
#define GPEC_SOURCE_H

#include "fusion_io_source.h"

#include <string>

class gpec_source : public fio_source {
 public:

  struct gpec_field_data {
    bool is_populated;
    bool is_complex;
    int ntor, nr, nz, n_comp;
    double *r, *z;
    double ***v_real, ***v_imag;
    gpec_field_data(const bool ic, const int nc) {
      is_populated = false;
      r = 0;
      z = 0;
      v_real = 0;
      v_imag = 0;
      is_complex = ic;
      n_comp = nc;
    }
    virtual ~gpec_field_data() {
      free();
    }
    int allocate(const int nr0, const int nz0);
    int free();
    int interpolate(const double r0, const double z0,
		    double* vr, double* vi) const;

  };

  gpec_field_data b0;  // axisymmetric field
  gpec_field_data b1;  // total perturbed field
  gpec_field_data bx;  // vacuum perturbed field

  int read_field_data(const char* filename, gpec_field_data* d);

 public:
  std::string path;
  
  gpec_source();
  virtual ~gpec_source();

  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);
  int get_available_fields(fio_field_list*) const;
  int get_series(const series_type, fio_series**);

  int open(const char*);
  int close();

  int center(double* r0, double* z0) const;
  int extent(double* r0, double* r1, double* z0, double* z1) const;
};

#endif
