#ifndef NIMROD_SOURCE_H
#define NIMROD_SOURCE_H

#ifdef FUSIONIO_ENABLE_NIMROD

#include "fusion_io_source.h"
#include <string>

class nimrod_source : public fio_source {
  std::string dump_file;

 public:
  nimrod_source();
  virtual ~nimrod_source();

  int open(const char*);
  int close();

  int get_available_fields(fio_field_list*) const;
  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);
  int get_int_parameter(const parameter_type, int*) const;
  int get_real_parameter(const parameter_type, double*) const;

  int sizeof_search_hint() const
  { return 2*sizeof(double); }
  int allocate_search_hint(fio_hint* s);
  int deallocate_search_hint(fio_hint* s);
};

#endif // FUSIONIO_ENABLE_NIMROD
#endif // NIMROD_SOURCE_H
