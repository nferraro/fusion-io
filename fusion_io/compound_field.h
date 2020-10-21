#ifndef COMPOUND_FIELD_H
#define COMPOUND_FIELD_H

#include <deque>
#include "fusion_io_field.h"
#include "fusion_io_defs.h"

class fio_compound_field : public fio_field {
  struct component_field {
    fio_field* field;
    int        op;
    double     factor;
    fio_hint   hint;
    component_field(fio_field* f, const int o, const int d,
		    const fio_hint h=0);
    component_field(const component_field& f);
  };
  typedef std::deque<component_field> field_list;
  field_list fields;
  int dim;

 public:
  fio_compound_field();
  fio_compound_field* clone() const { return new fio_compound_field(*this); }
  int dimension() const;
  int eval(const double*, double*, fio_hint=0);
  int eval_deriv(const double*, double*, fio_hint=0);
  int add_field(fio_field*, const int, const double, const fio_hint h=0);
};

#endif
