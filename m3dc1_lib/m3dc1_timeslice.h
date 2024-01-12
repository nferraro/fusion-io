#ifndef M3DC1_TIMESLICE_H
#define M3DC1_TIMESLICE_H

#include "m3dc1_mesh.h"
#include "m3dc1_coord_map.h"

class m3dc1_timeslice {
 public:
  typedef std::map<std::string, m3dc1_field*> m3dc1_field_map;
  m3dc1_field_map field_map;

  int ntor;
  int is_3d;
  int is_complex;
  int is_stell;

  double time;
  m3dc1_mesh* mesh;
  m3dc1_coord_map* map;

  m3dc1_timeslice();
  ~m3dc1_timeslice();

  bool get_field(const char*, m3dc1_field**) const;
};

#endif
