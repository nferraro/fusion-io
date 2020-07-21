#ifndef GEQDSK_SOURCE_H
#define GEQDSK_SOURCE_H

#include "fusion_io_source.h"

#include <string>

class geqdsk_source : public fio_source {
 public:
  double rmaxis, zmaxis;

  int nw, nh;
  double dx, dz;
  double rleft, zmid, zbottom;
  double simag, sibry;

  double* psi;
  double** psirz;
  double* fpol;
  double* ffprime;
  double* press;

 public:
  std::string filename;
  
  geqdsk_source();
  virtual ~geqdsk_source();

  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);
  int get_available_fields(fio_field_list*) const;
  int get_series(const series_type, fio_series**);

  int open(const char*);
  int close();

  int center(double* r0, double* z0) const;
  int extent(double* r0, double* r1, double* z0, double* z1) const;

  int interpolate_psi(const double r0, const double z0, 
		      double* psi) const;
};

#endif
