#include "m3dc1_file.h"
#include <math.h>
#include <libgen.h>

#include <iostream>

m3dc1_file::m3dc1_file()
{
  file = -1;
}

m3dc1_file::~m3dc1_file()
{
  close();
}

bool m3dc1_file::open(const char* filename)
{
  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  path = dirname(strdup(filename));

  if(file < 0) {
    std::cerr << "Error: couldn't open file " << filename << std::endl;
    return false;
  }

  return true; 
}

bool m3dc1_file::close()
{
  if(file < 0) return 0;
  if(H5Fclose(file) < 0) return false;
  file = -1;
  return true;
}

bool m3dc1_file::time_name(int time, char* name)
{
  if(time < 0)
    sprintf(name, "equilibrium");
  else
    sprintf(name, "time_%03d", time);
  return true;
}

hid_t m3dc1_file::open_timeslice(int t)
{
  char name[256];

  time_name(t, name);
  hid_t time_group = H5Gopen(file, name, H5P_DEFAULT);

  return time_group;
}

m3dc1_timeslice* m3dc1_file::load_timeslice(const int t)
{
  char name[256];
  m3dc1_timeslice_map::iterator i;

  i = timeslice_map.find(t);
  if(i != timeslice_map.end()) {
    //std::cerr << "Timeslice " << t << " has already been read." << std::endl;
    return &i->second;
  }

  hid_t time_group = open_timeslice(t);
  if(time_group < 0) {
    std::cerr << "Error opening time slice " << name << std::endl;
    return 0;
  }

  std::pair<m3dc1_timeslice_map::iterator,bool> tsi;
  tsi = timeslice_map.insert(
    m3dc1_timeslice_map::value_type(t, m3dc1_timeslice()));

  if(!tsi.second) {
    std::cerr << "Error inserting timeslice" << std::endl;
    return 0;
  }
  m3dc1_timeslice* ts = &tsi.first->second;

  if(t >= 0) read_parameter("ntor", &ts->ntor);
  else ts->ntor = 0;
  //  std::cerr << "ntor = " << ts->ntor << std::endl;

  read_parameter("icomplex", &ts->is_complex);
  read_parameter("3d", &ts->is_3d);
  //  std::cerr << "3d = " << ts->is_3d << std::endl;

  hid_t attr_id = H5Aopen(time_group, "time", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_DOUBLE, (void*)(&ts->time));
  H5Aclose(attr_id);
  
  ts->mesh = read_mesh(t);
  if(!ts->mesh) { 
    std::cerr << "Error reading mesh" << std::endl;
    return 0;
  }

  H5Gclose(time_group);

  return ts;
}

m3dc1_mesh* m3dc1_file::read_mesh(const int t)
{
  int nelms;

  hid_t time_group = open_timeslice(t);
  if(time_group < 0) {
    std::cerr << "Error opening time group" << std::endl;
    return 0;
  }

  hid_t mesh_group = H5Gopen(time_group, "mesh", H5P_DEFAULT);
  if(mesh_group < 0) {
    std::cerr << "Error opening mesh" << std::endl;
    return 0;
  }

  hid_t nelms_attr = H5Aopen(mesh_group, "nelms", H5P_DEFAULT);
  H5Aread(nelms_attr, H5T_NATIVE_INT, &nelms);
  H5Aclose(nelms_attr);
  if(nelms < 1) {
    std::cerr << "Error: too few elements" << std::endl;
    return 0;
  }

  int is_3d = 0;
  read_parameter("3d", &is_3d);  
  //  std::cerr << "is_3d = " << is_3d << std::endl;

  int version = 0;
  read_parameter("version", &version);

  int nfields;
  if(is_3d) nfields = 9;
  else nfields = 7;
  if(version>=15) nfields++;

  double* data = new double[nfields*nelms];

  hid_t mesh_dataset = H5Dopen(mesh_group, "elements", H5P_DEFAULT);
  H5Dread(mesh_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	  (void*)data);
  H5Dclose(mesh_dataset);

  m3dc1_mesh* mesh;
  if(is_3d) mesh = new m3dc1_3d_mesh(nelms);
  else mesh = new m3dc1_mesh(nelms);

  // set default mesh memory depth for searching
  // (this is being phased out in favor of connectivity tree)
  //  mesh->set_memory_depth(3);
  mesh->set_memory_depth(0);    

  int offset;
  if(version >= 15)
    offset = 8;
  else 
    offset = 7;

  for(int i=0; i<nelms; i++) {
    mesh->a[i]     =      data[i*nfields  ];
    mesh->b[i]     =      data[i*nfields+1];
    mesh->c[i]     =      data[i*nfields+2];
    mesh->co[i]    =  cos(data[i*nfields+3]);
    mesh->sn[i]    =  sin(data[i*nfields+3]);
    mesh->x[i]     =      data[i*nfields+4];
    mesh->z[i]     =      data[i*nfields+5];
    mesh->bound[i] = (int)data[i*nfields+6];
    if(version >= 15)
      mesh->region[i] = (int)data[i*nfields+7];
    else
      mesh->region[i] = 0;

    if(is_3d) {
      ((m3dc1_3d_mesh*)mesh)->d[i]   =      data[i*nfields+offset];
      ((m3dc1_3d_mesh*)mesh)->phi[i] =      data[i*nfields+offset+1];
    }
  }
  delete[] data;
  
  // Determine geometry
  double rzero = 0.;
  int itor = 0;
  read_parameter("itor", &itor);
  mesh->toroidal = (itor==1);
  if(itor==0) {
    read_parameter("rzero", &rzero);
    mesh->period = 2.*M_PI*rzero;
  } else {
    mesh->period = 2.*M_PI;
  }

  read_parameter("nplanes", &(mesh->nplanes));

  // Calculate connectivity tree
  std::string stash = path + "/.stash" + std::to_string(t);
  mesh->find_neighbors(stash.c_str());

  H5Gclose(mesh_group);
  H5Gclose(time_group);

  return mesh;
}


bool m3dc1_file::read_parameter(const char* n, double* d)
{
  hid_t attr_id = H5Aopen(file, n, H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_DOUBLE, d);
  H5Aclose(attr_id);

  return true;  
}

bool m3dc1_file::read_parameter(const char* n, int* d)
{
  hid_t attr_id = H5Aopen(file, n, H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, d);
  H5Aclose(attr_id);

  return true;  
}

bool m3dc1_file::get_slice_time(const int i, double* t)
{
  m3dc1_timeslice* ts = load_timeslice(i);
  if(ts==NULL) return false;
  
  *t = ts->time;
  return true;
}

m3dc1_scalar_list* m3dc1_file::get_slice_times()
{
  const char name[] = "slice_times";

  m3dc1_scalar_map::iterator i = scalar_map.find(name);
  if(i != scalar_map.end())
    return &i->second;

  int ntimes;
  if(!read_parameter("ntime", &ntimes)) return NULL;

  m3dc1_scalar_list* s = &(scalar_map[name] = m3dc1_scalar_list(ntimes));

  for(int t=0; t<ntimes; t++) {
    m3dc1_timeslice* ts = load_timeslice(t);
    if(ts==NULL) {
      std::cerr << "Error reading timeslice " << t << std::endl;
      return NULL;
    }
    (*s)[t] = ts->time;
  }

  return s;
}

m3dc1_scalar_list* m3dc1_file::read_scalar(const char* name)
{
  m3dc1_scalar_map::iterator i;
  i = scalar_map.find(name);
  if(i != scalar_map.end()) {
    // std::cerr << "Scalar " << name << " has already been read" << std::endl;
    return &i->second;
  }

  // Open fields group
  hid_t scalar_group = H5Gopen(file, "scalars", H5P_DEFAULT);
  if(scalar_group < 0) {
    std::cerr << "Error opening scalars group" << std::endl;
    return 0;
  }

  // Open field dataset
  hid_t scalar_dataset = H5Dopen(scalar_group, name, H5P_DEFAULT);
  if(scalar_dataset < 0) {
    std::cerr << "Error opening scalar " << name << std::endl;
    H5Gclose(scalar_group);
    return 0;
  }

  hid_t scalar_space = H5Dget_space(scalar_dataset);
  hssize_t n = H5Sget_simple_extent_npoints(scalar_space);
  H5Sclose(scalar_space);

  if(n==0) {
    std::cerr << "Error: no elements in dataset " << name << std::endl;
    return 0;
  }

  m3dc1_scalar_list* s = &(scalar_map[name] = m3dc1_scalar_list(n));

  // Read Field Data
  H5Dread(scalar_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  	  (void*)(&((*s)[0])));

  H5Dclose(scalar_dataset);
  H5Gclose(scalar_group);

  return s;
}

m3dc1_field* m3dc1_file::load_field(const char* n, const int t, 
				    int options)
{
  hid_t fields_group, field_dataset, fieldi_dataset;
  std::map<std::string, m3dc1_field*>::iterator i;
  std::string name = n;
  std::string name_i = name + "_i";

  //  std::cerr << "Reading field " << n << " at timelice " << t << std::endl;

  m3dc1_timeslice* ts = load_timeslice(t);
  if(!ts) {
    std::cerr << "Error reading timeslice" << std::endl;
    return 0;
  }

  // Check if field has already been read
  i = ts->field_map.find(name);
  if(i != ts->field_map.end()) {
    // std::cerr << "Field " << name << " has already been read" << std::endl;
    return i->second;
  }

  m3dc1_compound_field* cfield = 0;

  //  std::cerr << "options = " << options << std::endl;

  // Check if this will be a compund field
  // Are we adding in the equilibrium?
  if((options & M3DC1_ADD_EQUILIBRIUM)==M3DC1_ADD_EQUILIBRIUM) {
    int eqsubtract;
    read_parameter("eqsubtract", &eqsubtract);
    std::cerr << "eqsubtract = " << eqsubtract << std::endl;
    if(eqsubtract==1 && t !=-1) {
      std::cerr << "Adding equilibrium to compound field" << std::endl;
      cfield = new m3dc1_compound_field;
      m3dc1_field* temp_field = load_field(n, -1);
      if(!temp_field) {
	std::cerr << "Error reading equilibrium field" << std::endl;
	delete(cfield);
	return 0;
      }
      cfield->subfield.push_back(temp_field);
    }
  }
  // Are we adding in external fields?
  if((options & M3DC1_ADD_EXTERNAL)==M3DC1_ADD_EXTERNAL) {
    int extsubtract;
    read_parameter("extsubtract", &extsubtract);
    if(extsubtract==1) {
      if(name=="psi" || name=="f" || name=="I") {
	std::cerr << "Adding external field to compound field" << std::endl;
	std::string name_ext = name + "_ext";
	if(!cfield) cfield = new m3dc1_compound_field;
	m3dc1_field* temp_field = load_field(name_ext.c_str(), t);
	if(!temp_field) {
	  std::cerr << "Error reading external field" << std::endl;
	  delete(cfield);
	  return 0;
	}
	cfield->subfield.push_back(temp_field);
      }
    }
  }
 
  bool is_3d = (ts->is_3d==1);
  bool is_complex = (ts->is_complex==1);

  // open time group
  hid_t time_group = open_timeslice(t);
  if(time_group < 0) {
    std::cerr << "Error opening timeslice" << std::endl;
    return 0;
  }

  // Open fields group
  fields_group = H5Gopen(time_group, "fields", H5P_DEFAULT);
  if(fields_group < 0) {
    std::cerr << "Error opening fields group" << std::endl;
    return 0;
  }

  // Open field dataset
  field_dataset = H5Dopen(fields_group, name.data(), H5P_DEFAULT);
  if(field_dataset < 0) {
    std::cerr << "Error opening field " << name << std::endl;
    H5Gclose(fields_group);
    return 0;
  }

  // For complex fields, open imaginary field dataset
  if(is_complex) {
    fieldi_dataset = H5Dopen(fields_group, name_i.data(), H5P_DEFAULT);
    if(fieldi_dataset < 0) {
      std::cerr << "Error opening field " << name_i << std::endl;
      is_complex = false;
    }
  }

  // Create new field
  m3dc1_field* field;
  if(is_3d)
    field = new m3dc1_3d_field(ts->mesh);
  else if(is_complex)
    field = new m3dc1_complex_field(ts->mesh,ts->ntor);
  else
    field = new m3dc1_field(ts->mesh);

  // Read Field Data
  H5Dread(field_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  	  (void*)field->data);
  H5Dclose(field_dataset);

  if(is_complex) {
    H5Dread(fieldi_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	    (void*)((m3dc1_complex_field*)field)->data_i);
    H5Dclose(fieldi_dataset);
  }
  
  H5Gclose(fields_group);
  H5Gclose(time_group);

  field->time = ts->time;

  if(cfield) {
    cfield->subfield.push_back(field);
    // Add field to list of fields that have been read
    std::string name_new = name + "_base";
    ts->field_map.insert(
      m3dc1_timeslice::m3dc1_field_map::value_type(name_new, field));
    ts->field_map.insert(
      m3dc1_timeslice::m3dc1_field_map::value_type(name, cfield));
    return cfield;
  } else {
    ts->field_map.insert(
      m3dc1_timeslice::m3dc1_field_map::value_type(name, field));
    return field;
  }
}

bool m3dc1_file::unload_field(const char* n, int t)
{
  m3dc1_timeslice_map::iterator its;

  // Check if timeslice is loaded
  its = timeslice_map.find(t);
  if(its == timeslice_map.end()) {
    std::cerr << "Timeslice " << t << " is not loaded."  << std::endl;
    return false;
  }
  m3dc1_timeslice* ts = &its->second;

  // Check if field is loaded
  std::string name = n;
  std::map<std::string, m3dc1_field*>::iterator i;
  i = ts->field_map.find(name);
  if(i == ts->field_map.end()) {
    std::cerr << "Field " << name << " is not loaded." << std::endl;
    return false;
  }
  
  // Unload field
  ts->field_map.erase(i);
  return true;
}

bool m3dc1_file::extent(const int t, double* r0, double* r1, 
			double* phi0, double* phi1, 
			double* z0, double* z1)
{
  m3dc1_timeslice* ts = load_timeslice(t);

  if(!ts) return false;

  ts->mesh->extent(r0, r1, phi0, phi1, z0, z1); 
  return true;
}


