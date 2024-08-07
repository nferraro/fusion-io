#ifndef M3DC1_FILE_H
#define M3DC1_FILE_H

#include "m3dc1_field.h"
#include "m3dc1_stell_field.h"
#include "m3dc1_timeslice.h"
#include <map>
#include <string>
#include <vector>

typedef std::vector<double> m3dc1_scalar_list;

class m3dc1_file {
  int version;
  hid_t file;

  typedef std::map<int, m3dc1_timeslice> m3dc1_timeslice_map;
  m3dc1_timeslice_map timeslice_map;

  typedef std::map<std::string, m3dc1_scalar_list> m3dc1_scalar_map;
  m3dc1_scalar_map scalar_map;

  bool time_name(int time, char* name);
  hid_t open_timeslice(int time);
  m3dc1_timeslice* load_timeslice(const int time);
  m3dc1_mesh* read_mesh(int time);
  void load_adjacency(m3dc1_mesh* mesh, const hid_t mesh_group);

 public:
  const static int M3DC1_ADD_EQUILIBRIUM = 1;
  const static int M3DC1_ADD_EXTERNAL    = 2;
  const static int M3DC1_TOTAL           = 3;
  const static int M3DC1_LOGICAL_COORDS  = 8;

 public:
  m3dc1_file();
  ~m3dc1_file();

  bool open(const char* filename);
  bool close();

  bool get_slice_time(const int, double*);
  m3dc1_scalar_list* get_slice_times();
  m3dc1_scalar_list* read_scalar(const char*);
  m3dc1_field* load_field(const char* name, const int time, 
			  const int options=0);
  bool unload_field(const char* name, const int time);
  bool read_parameter(const char*, int*); 
  bool read_parameter(const char*, double*);
  bool read_adjacency(m3dc1_mesh* mesh, const hid_t mesh_group);

  bool extent(int, double*, double*, double*, double*, double*, double*);
};

#endif
