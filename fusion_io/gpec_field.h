#ifndef GPEC_FIELD_H
#define GPEC_FIELD_H

#include "gpec_source.h"


class gpec_field : public fio_field {
 protected:
  gpec_source* source;
  double linfac;
  double phase;
  int ilin;
  int time;

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
  gpec_source::gpec_field_data *b0, *b1;

 public:
  gpec_magnetic_field(gpec_source* s)
    : gpec_field(s) { }
  gpec_magnetic_field* clone() const { return new gpec_magnetic_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*, void* =0);
};


#endif
