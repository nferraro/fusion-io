#include "m3dc1_timeslice.h"

m3dc1_timeslice::m3dc1_timeslice()
  : mesh(0), map(0)
{
  is_3d = 0;
  is_stell = 0;
}

m3dc1_timeslice::~m3dc1_timeslice()
{
  std::map<std::string, m3dc1_field*>::iterator i = field_map.begin();
  while(i != field_map.end()) {
    delete(i->second);
    i++;
  }

  if(map) delete(map);
  if(mesh) delete(mesh);
}
