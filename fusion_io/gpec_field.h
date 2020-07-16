#ifndef GPEC_FIELD_H
#define GPEC_FIELD_H

#include "gpec_source.h"

class gpec_series : public fio_series {
  double data;

 public:
  gpec_series(const double d) { data = d; }
  virtual ~gpec_series()
    { }

  virtual int eval(const double, double*);
  virtual int bounds(double*, double*) const;
};


class gpec_field : public fio_field {
 protected:
  gpec_source* source;
  double linfac;
  int ilin;

 public:
  gpec_field(gpec_source* s) 
    { source = s; }
  virtual ~gpec_field()
    { }

  virtual int load(const fio_option_list*) = 0;
  virtual int dimension() const = 0;
  virtual int eval(const double*, double*, void* =0) = 0;
};

class gpec_magnetic_field : public gpec_field {
 public:
  gpec_magnetic_field(gpec_source* s)
    : gpec_field(s) { }
  gpec_magnetic_field* clone() const { return new gpec_magnetic_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*, void* =0);
};


#endif
