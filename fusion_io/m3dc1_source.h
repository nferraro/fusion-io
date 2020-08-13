#ifndef M3DC1_SOURCE_H
#define M3DC1_SOURCE_H

#include <m3dc1_file.h>

#include "fusion_io_source.h"

class m3dc1_source : public fio_source {
 public:
  m3dc1_file file;
  double bzero, rzero, z_ion, ion_mass, period;
  double n0, L0, B0, p0, t0, v0, J0, Phi0, temp0;
  int linear, eqsubtract, extsubtract, icomplex, i3d, version, itor, ntor, ntime;
  int kprad_z;
  fio_species ion_species;

 public:
  int open(const char*);
  int close();

  int get_available_fields(fio_field_list*) const;
  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);
  int get_series(const series_type, fio_series**);
  int get_int_parameter(const parameter_type, int*) const;
  int get_real_parameter(const parameter_type, double*) const;
  int get_slice_time(const int, double*);

  int sizeof_search_hint() const
  { return sizeof(int); }
  int allocate_search_hint(void** s);
  int deallocate_search_hint(void** s);
};


#endif
