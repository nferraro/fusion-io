#ifndef MARS_FIELD_H
#define MARS_FIELD_H

#include "fusion_io_field.h"
#include "mars_source.h"
#include "options.h"


class mars_field : public fio_field {
 protected:
  int part, time;
  double linfac, phase;
  mars_source* source;
 public:
  mars_field(mars_source* s) 
    { source = s; }
  virtual int load(const fio_option_list*) = 0;

  int get_indices(const double*, int*, int*, double*, double*,
		  double*, double*, double*, double*);
};



// V
class mars_fluid_velocity : public mars_field {
 public:
  mars_fluid_velocity(mars_source* s) 
    : mars_field(s) { }
  mars_fluid_velocity* clone() const 
  { return new mars_fluid_velocity(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*, void* =0);
};


// B
class mars_magnetic_field : public mars_field {
  double **b[3];
 public:
  mars_magnetic_field(mars_source* s) 
    : mars_field(s) { }
  mars_magnetic_field* clone() const 
  { return new mars_magnetic_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*, void* =0);
};


#endif
