#include "compound_field.h"
#include "fusion_io_defs.h"

#include <iostream>

fio_compound_field::
component_field::component_field(fio_field* f, const int o, const int d,
				 const fio_hint h)
{
  field = f;
  op = o;
  factor = d;
  hint = h;
}

fio_compound_field::
component_field::component_field(const component_field& f)
{
  field = f.field;
  op = f.op;
  factor = f.factor;
  hint = f.hint;
}

fio_compound_field::fio_compound_field()
{
  dim = 0;
}

int fio_compound_field::dimension() const
{
  return dim;
}

int fio_compound_field::add_field(fio_field* f, const int op, const double d,
				  const fio_hint h)
{
  if(fields.empty()) {
    dim = f->dimension();

    if(op != FIO_ADD) {
      std::cerr << "Error: operation for first field must be FIO_ADD"
		<< std::endl;
      return 1;
    }
  } else {
    if(f->dimension() != dimension()) {
      std::cerr << "Error: adding field of incompatible dimension" 
		<< std::endl;
      return 1;
    }
  }

  fields.push_back(component_field(f, op, d, h));
  return FIO_SUCCESS;
}

int fio_compound_field::eval(const double* x, double* v, fio_hint h)
{
  double* z = new double[dimension()];

  for(int j=0; j<dimension(); j++)
    v[j] = 0.;

  field_list::iterator i = fields.begin();
  while(i != fields.end()) {
    int result = i->field->eval(x, z, i->hint);
    if(result != FIO_SUCCESS) return result;

    for(int j=0; j<dimension(); j++)
      switch(i->op) {
      case(FIO_ADD):       v[j] += i->factor*z[j]; break;
      case(FIO_MULTIPLY):  v[j] *= i->factor*z[j]; break;
      case(FIO_DIVIDE):    v[j] /= i->factor*z[j]; break;
      }     
    i++;
  }

  delete[] z;

  return FIO_SUCCESS;
}

int fio_compound_field::eval_deriv(const double* x, double* v, fio_hint h)
{
  double* z = new double[3*dimension()];

  for(int j=0; j<3*dimension(); j++)
    v[j] = 0.;

  field_list::iterator i = fields.begin();
  while(i != fields.end()) {
    int result = i->field->eval_deriv(x, z, i->hint);
    if(result != FIO_SUCCESS) return result;

    for(int j=0; j<3*dimension(); j++)
      switch(i->op) {
      case(FIO_ADD):       v[j] += i->factor*z[j]; break;
      default:
	std::cerr << "Eval deriv for op != FIO_ADD unsupported" << std::endl;
	return FIO_UNSUPPORTED;
      }     
    i++;
  }

  delete[] z;

  return FIO_SUCCESS;
}

