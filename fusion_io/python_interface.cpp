#include "fusion_io.h"
#include "fusion_io_c.h"
#include "m3dc1_numpy_source.h"
#include <Python.h>

#include <iostream>
#include <cstring>
#include <cmath>

extern "C" {
  static PyObject* fio_add_field_py(PyObject*, PyObject*);
  static PyObject* fio_allocate_hint_py(PyObject*, PyObject*);
  static PyObject* fio_close_field_py(PyObject*, PyObject*);
  static PyObject* fio_close_series_py(PyObject*, PyObject*);
  static PyObject* fio_close_source_py(PyObject*, PyObject*);
  static PyObject* fio_create_compound_field_py(PyObject*, PyObject*);
  static PyObject* fio_eval_scalar_field_py(PyObject*, PyObject*);
  static PyObject* fio_eval_scalar_field_deriv_py(PyObject*, PyObject*);
  static PyObject* fio_eval_series_py(PyObject*, PyObject*);
  static PyObject* fio_eval_vector_field_py(PyObject*, PyObject*);
  static PyObject* fio_eval_vector_field_deriv_py(PyObject*, PyObject*);
  static PyObject* fio_eval_tensor_field_py(PyObject*, PyObject*);
  static PyObject* fio_free_hint_py(PyObject*, PyObject*);
  static PyObject* fio_get_available_fields_py(PyObject*, PyObject*);
  static PyObject* fio_get_field_py(PyObject*, PyObject*);
  static PyObject* fio_get_field_name_py(PyObject*, PyObject*);
  static PyObject* fio_get_option_name_py(PyObject*, PyObject*);
  static PyObject* fio_get_options_py(PyObject*, PyObject*);
  static PyObject* fio_get_series_py(PyObject*, PyObject*);
  static PyObject* fio_get_series_bounds_py(PyObject*, PyObject*);
  static PyObject* fio_open_source_py(PyObject*, PyObject*);
  static PyObject* fio_set_int_option_py(PyObject*, PyObject*);
  static PyObject* fio_set_str_option_py(PyObject*, PyObject*);
  static PyObject* fio_set_real_option_py(PyObject*, PyObject*);
  static PyObject* fio_get_int_parameter_py(PyObject*, PyObject*);
  static PyObject* fio_get_real_parameter_py(PyObject*, PyObject*);
  static PyObject* fio_get_real_field_parameter_py(PyObject*, PyObject*);
  static PyObject* fio_set_quiet_option_py(PyObject*, PyObject*);
  // numpy injection
  static PyObject* fio_open_numpy_source_py(PyObject*, PyObject*);
  static PyObject* fio_set_numpy_mesh_py(PyObject*, PyObject*);
  static PyObject* fio_set_numpy_field_py(PyObject*, PyObject*);
  static PyObject* fio_set_numpy_scalar_py(PyObject*, PyObject*);
  static PyObject* fio_set_numpy_parameters_py(PyObject*, PyObject*);

  static PyMethodDef fio_methods[] = {
    {"add_field",                fio_add_field_py,                METH_VARARGS, ""},
    {"allocate_hint",            fio_allocate_hint_py,            METH_VARARGS, ""},
    {"close_field",              fio_close_field_py,              METH_VARARGS, ""},
    {"close_series",             fio_close_series_py,             METH_VARARGS, ""},
    {"close_source",             fio_close_source_py,             METH_VARARGS, ""},
    {"create_compound_field",    fio_create_compound_field_py,    METH_VARARGS, ""},
    {"eval_scalar_field",        fio_eval_scalar_field_py,        METH_VARARGS, ""},
    {"eval_scalar_field_deriv",  fio_eval_scalar_field_deriv_py,  METH_VARARGS, ""},
    {"eval_series",              fio_eval_series_py,              METH_VARARGS, ""},
    {"eval_vector_field",        fio_eval_vector_field_py,        METH_VARARGS, ""},
    {"eval_vector_field_deriv",  fio_eval_vector_field_deriv_py,  METH_VARARGS, ""},
    {"eval_tensor_field",        fio_eval_tensor_field_py,        METH_VARARGS, ""},
    {"free_hint",                fio_free_hint_py,                METH_VARARGS, ""},
    {"get_available_fields",     fio_get_available_fields_py,     METH_VARARGS, ""},
    {"get_options",              fio_get_options_py,              METH_VARARGS, ""},
    {"get_field",                fio_get_field_py,                METH_VARARGS, ""},
    {"get_field_name",           fio_get_field_name_py,           METH_VARARGS, ""},
    {"get_option_name",          fio_get_option_name_py,          METH_VARARGS, ""},
    {"get_series",               fio_get_series_py,               METH_VARARGS, ""},
    {"get_series_bounds",        fio_get_series_bounds_py,        METH_VARARGS, ""},
    {"open_source",              fio_open_source_py,              METH_VARARGS, ""},
    {"set_int_option",           fio_set_int_option_py,           METH_VARARGS, ""},
    {"set_str_option",           fio_set_str_option_py,           METH_VARARGS, ""},
    {"set_real_option",          fio_set_real_option_py,          METH_VARARGS, ""},
    {"get_int_parameter",        fio_get_int_parameter_py,        METH_VARARGS, ""},
    {"get_real_parameter",       fio_get_real_parameter_py,       METH_VARARGS, ""},
    {"get_real_field_parameter", fio_get_real_field_parameter_py, METH_VARARGS, ""},
    {"set_quiet_option",         fio_set_quiet_option_py,         METH_VARARGS, ""},
    // numpy injection API
    {"open_numpy_source",        fio_open_numpy_source_py,        METH_VARARGS,
     "open_numpy_source() -> int\n"
     "Create a numpy-backed M3D-C1 source. Returns an integer handle\n"
     "compatible with all other fio_py functions."},
    {"set_numpy_mesh",           fio_set_numpy_mesh_py,           METH_VARARGS,
     "set_numpy_mesh(handle, nelms, a, b, c, co, sn, x, z,\n"
     "               bound, region, nplanes, period, toroidal) -> None\n"
     "Inject mesh arrays into the numpy source.\n"
     "a,b,c,co,sn,x,z: bytes of float64 arrays of length nelms.\n"
     "bound,region:     bytes of int32  arrays of length nelms.\n"
     "co and sn must be cos/sin of the element angle (not the angle itself)."},
    {"set_numpy_field",          fio_set_numpy_field_py,          METH_VARARGS,
     "set_numpy_field(handle, timeslice, field_name,\n"
     "                data_real, data_imag_or_None, ntor) -> None\n"
     "Inject field coefficient arrays for one field at one timeslice.\n"
     "data_real: bytes of float64, length nelms*20.\n"
     "data_imag: bytes of float64, same length, or None for real fields.\n"
     "ntor: toroidal mode number (0 if not complex)."},
    {"set_numpy_scalar",         fio_set_numpy_scalar_py,         METH_VARARGS,
     "set_numpy_scalar(handle, scalar_name, values_bytes) -> None\n"
     "Inject a scalar time series (bytes of float64) into the source.\n"
     "Use scalar_name='slice_times' for the time axis."},
    {"set_numpy_parameters",     fio_set_numpy_parameters_py,     METH_VARARGS,
     "set_numpy_parameters(handle,\n"
     "                     bzero, rzero, ntime, ntor,\n"
     "                     ion_mass, z_ion, itor, numvar, version,\n"
     "                     n0, L0, B0) -> None\n"
     "Set physical normalization parameters (n0, L0, B0 in CGS as stored\n"
     "in the M3D-C1 HDF5 file). Must be called before get_field()."},
    {NULL, NULL, 0, NULL}
  };

  static PyObject* fio_error;

#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef fio_py_module =
  {
      PyModuleDef_HEAD_INIT,
      "fio_py",
      "",
      -1,
      fio_methods
  };

  PyMODINIT_FUNC PyInit_fio_py()
#else
  PyMODINIT_FUNC initfio_py()
#endif
  {
    PyObject *m, *d;
#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&fio_py_module);
#else
    m = Py_InitModule("fio_py", fio_methods);
#endif

    if(m == NULL)
#if PY_MAJOR_VERSION >= 3
      return NULL;
#else
      return;
#endif

    fio_error = PyErr_NewException("fio.error", NULL, NULL);
    Py_INCREF(fio_error);
    PyModule_AddObject(m, "error", fio_error);

    d = PyModule_GetDict(m);
#if PY_MAJOR_VERSION >= 3
    PyDict_SetItemString(d, "test", PyLong_FromLong(12));
#else
    PyDict_SetItemString(d, "test", PyInt_FromLong(12));
#endif

#undef FUSION_IO_DEFS_H
#define PYTHON
#define PYTHON_DICT d
#include "fusion_io_defs.h"
#if PY_MAJOR_VERSION >= 3
    return m;
#endif
  }
} // extern "C"

// ---------------------------------------------------------------------------
// Helper: retrieve m3dc1_numpy_source* from a handle, or set Python exception
// ---------------------------------------------------------------------------

static m3dc1_numpy_source* get_numpy_source(int handle)
{
    fio_source* src = nullptr;
    if(fio_get_source(handle, &src) != FIO_SUCCESS || src == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid source handle");
        return nullptr;
    }
    m3dc1_numpy_source* ns = dynamic_cast<m3dc1_numpy_source*>(src);
    if(!ns) {
        PyErr_SetString(PyExc_TypeError,
                        "Source handle does not refer to a numpy source");
        return nullptr;
    }
    return ns;
}

// ---------------------------------------------------------------------------
// Existing bindings (unchanged from original python_interface.cpp)
// ---------------------------------------------------------------------------

PyObject* fio_add_field_py(PyObject* self, PyObject *args)
{
  int icfield, ifield, op;
  double fac;
  if(!PyArg_ParseTuple(args, "iiid", &icfield, &ifield, &op, &fac))
    return NULL;
  int ierr = fio_add_field(icfield, ifield, op, fac);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_allocate_hint_py(PyObject* self, PyObject *args)
{
  int isrc;
  if(!PyArg_ParseTuple(args, "i", &isrc)) return NULL;
  int sz = fio_sizeof_search_hint(isrc);
  if(sz == 0) return NULL;
  void* hint = calloc(1, sz);
  return PyLong_FromVoidPtr(hint);
}

PyObject* fio_free_hint_py(PyObject* self, PyObject *args)
{
  PyObject* hint_py;
  if(!PyArg_ParseTuple(args, "O", &hint_py)) return NULL;
  free(PyLong_AsVoidPtr(hint_py));
  Py_RETURN_NONE;
}

PyObject* fio_close_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  if(!PyArg_ParseTuple(args, "i", &ifield)) return NULL;
  int ierr = fio_close_field(ifield);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_close_source_py(PyObject* self, PyObject *args)
{
  int isrc;
  if(!PyArg_ParseTuple(args, "i", &isrc)) return NULL;
  int ierr = fio_close_source(isrc);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_close_series_py(PyObject* self, PyObject *args)
{
  int isrc;
  if(!PyArg_ParseTuple(args, "i", &isrc)) return NULL;
  int ierr = fio_close_series(isrc);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_create_compound_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  int ierr = fio_create_compound_field(&ifield);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("i", ifield);
}

PyObject* fio_eval_scalar_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;
  if(!PyArg_ParseTuple(args, "i(ddd)O",
                       &ifield, &x[0], &x[1], &x[2], &hint_py))
    return NULL;
  void* hint = PyLong_AsVoidPtr(hint_py);
  double v;
  int ierr = fio_eval_field(ifield, x, &v, hint);
  if(ierr != FIO_SUCCESS) return Py_BuildValue("");
  hint_py = PyLong_FromVoidPtr(hint);
  return Py_BuildValue("d", v);
}

PyObject* fio_eval_scalar_field_deriv_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;
  if(!PyArg_ParseTuple(args, "i(ddd)O",
                       &ifield, &x[0], &x[1], &x[2], &hint_py))
    return NULL;
  void* hint = PyLong_AsVoidPtr(hint_py);
  double v[3];
  int ierr = fio_eval_field_deriv(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS) return Py_BuildValue("(2)", Py_None);
  hint_py = PyLong_FromVoidPtr(hint);
  return Py_BuildValue("(ddd)", v[0], v[1], v[2]);
}

PyObject* fio_eval_series_py(PyObject* self, PyObject *args)
{
  int iseries;
  double x;
  if(!PyArg_ParseTuple(args, "id", &iseries, &x)) return NULL;
  double v;
  int ierr = fio_eval_series(iseries, x, &v);
  if(ierr != FIO_SUCCESS) return Py_BuildValue("");
  return Py_BuildValue("d", v);
}

PyObject* fio_eval_vector_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;
  if(!PyArg_ParseTuple(args, "i(ddd)O",
                       &ifield, &x[0], &x[1], &x[2], &hint_py))
    return NULL;
  void* hint = PyLong_AsVoidPtr(hint_py);
  double v[3];
  int ierr = fio_eval_field(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS) return Py_BuildValue("(2)", Py_None);
  hint_py = PyLong_FromVoidPtr(hint);
  return Py_BuildValue("(ddd)", v[0], v[1], v[2]);
}

PyObject* fio_eval_vector_field_deriv_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;
  if(!PyArg_ParseTuple(args, "i(ddd)O",
                       &ifield, &x[0], &x[1], &x[2], &hint_py))
    return NULL;
  void* hint = PyLong_AsVoidPtr(hint_py);
  double v[9];
  int ierr = fio_eval_field_deriv(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS) return Py_BuildValue("(8)", Py_None);
  hint_py = PyLong_FromVoidPtr(hint);
  return Py_BuildValue("(ddddddddd)",
                       v[0], v[1], v[2], v[3], v[4],
                       v[5], v[6], v[7], v[8]);
}

PyObject* fio_eval_tensor_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;
  if(!PyArg_ParseTuple(args, "i(ddd)O",
                       &ifield, &x[0], &x[1], &x[2], &hint_py))
    return NULL;
  void* hint = PyLong_AsVoidPtr(hint_py);
  double v[9];
  int ierr = fio_eval_field(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS) return Py_BuildValue("(8)", Py_None);
  hint_py = PyLong_FromVoidPtr(hint);
  return Py_BuildValue("(ddddddddd)",
                       v[0], v[1], v[2], v[3], v[4],
                       v[5], v[6], v[7], v[8]);
}

PyObject* fio_get_available_fields_py(PyObject* self, PyObject *args)
{
  int isrc;
  if(!PyArg_ParseTuple(args, "i", &isrc)) return NULL;
  field_type* f;
  int n;
  int ierr = fio_get_available_fields(isrc, &n, &f);
  if(ierr != FIO_SUCCESS) return NULL;
  PyObject* list = PyList_New(n);
  for(int i = 0; i < n; i++)
#if PY_MAJOR_VERSION >= 3
    PyList_SET_ITEM(list, i, PyLong_FromLong(f[i]));
#else
    PyList_SET_ITEM(list, i, PyInt_FromLong(f[i]));
#endif
  return list;
}

PyObject* fio_get_field_py(PyObject* self, PyObject *args)
{
  int isrc, itype, handle;
  if(!PyArg_ParseTuple(args, "ii", &isrc, &itype)) return NULL;
  int ierr = fio_get_field(isrc, itype, &handle);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("i", handle);
}

PyObject* fio_get_field_name_py(PyObject* self, PyObject *args)
{
  int itype;
  if(!PyArg_ParseTuple(args, "i", &itype)) return NULL;
  std::string name;
  int ierr = fio_get_field_name(itype, &name);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("s", name.c_str());
}

PyObject* fio_get_option_name_py(PyObject* self, PyObject *args)
{
  int iopt;
  if(!PyArg_ParseTuple(args, "i", &iopt)) return NULL;
  std::string name;
  int ierr = fio_get_option_name(iopt, &name);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("s", name.c_str());
}

PyObject* fio_get_options_py(PyObject* self, PyObject *args)
{
  int isrc;
  if(!PyArg_ParseTuple(args, "i", &isrc)) return NULL;
  int ierr = fio_get_options(isrc);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_get_int_parameter_py(PyObject* self, PyObject *args)
{
  int isrc, t, p;
  if(!PyArg_ParseTuple(args, "ii", &isrc, &t)) return NULL;
  int ierr = fio_get_int_parameter(isrc, t, &p);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("i", p);
}

PyObject* fio_get_real_parameter_py(PyObject* self, PyObject *args)
{
  int isrc, t;
  double p;
  if(!PyArg_ParseTuple(args, "ii", &isrc, &t)) return NULL;
  int ierr = fio_get_real_parameter(isrc, t, &p);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("d", p);
}

PyObject* fio_get_real_field_parameter_py(PyObject* self, PyObject *args)
{
  int ifield, t;
  double p;
  if(!PyArg_ParseTuple(args, "ii", &ifield, &t)) return NULL;
  int ierr = fio_get_real_field_parameter(ifield, t, &p);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("d", p);
}

PyObject* fio_get_series_py(PyObject* self, PyObject *args)
{
  int isrc, itype, handle;
  if(!PyArg_ParseTuple(args, "ii", &isrc, &itype)) return NULL;
  int ierr = fio_get_series(isrc, itype, &handle);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("i", handle);
}

PyObject* fio_get_series_bounds_py(PyObject* self, PyObject *args)
{
  int iseries;
  double tmin, tmax;
  if(!PyArg_ParseTuple(args, "i", &iseries)) return NULL;
  int ierr = fio_get_series_bounds(iseries, &tmin, &tmax);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("dd", tmin, tmax);
}

PyObject* fio_open_source_py(PyObject* self, PyObject *args)
{
  int type;
  char* filename;
  if(!PyArg_ParseTuple(args, "is", &type, &filename)) return NULL;
  int handle;
  int ierr = fio_open_source(type, filename, &handle);
  if(ierr != FIO_SUCCESS) return NULL;
  return Py_BuildValue("i", handle);
}

PyObject* fio_set_int_option_py(PyObject* self, PyObject *args)
{
  int iopt, v;
  if(!PyArg_ParseTuple(args, "ii", &iopt, &v)) return NULL;
  int ierr = fio_set_int_option(iopt, v);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_set_str_option_py(PyObject* self, PyObject *args)
{
  int iopt;
  char* v;
  if(!PyArg_ParseTuple(args, "is", &iopt, &v)) return NULL;
  int ierr = fio_set_str_option(iopt, v);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_set_real_option_py(PyObject* self, PyObject *args)
{
  int iopt;
  double v;
  if(!PyArg_ParseTuple(args, "id", &iopt, &v)) return NULL;
  int ierr = fio_set_real_option(iopt, v);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_set_quiet_option_py(PyObject* self, PyObject *args)
{
  bool v;
  if(!PyArg_ParseTuple(args, "p", &v)) return NULL;
  int ierr = fio_set_quiet_option(v);
  if(ierr != FIO_SUCCESS) return NULL;
  Py_RETURN_NONE;
}

// ---------------------------------------------------------------------------
// New numpy-injection bindings
// ---------------------------------------------------------------------------

/*
 * open_numpy_source() -> int handle
 */
PyObject* fio_open_numpy_source_py(PyObject* self, PyObject* args)
{
    m3dc1_numpy_source* ns = new m3dc1_numpy_source();
    int handle;
    int ierr = fio_register_source(ns, &handle);
    if(ierr != FIO_SUCCESS) {
        delete ns;
        PyErr_SetString(PyExc_RuntimeError,
                        "fio_open_numpy_source: failed to register source");
        return nullptr;
    }
    return Py_BuildValue("i", handle);
}

/*
 * set_numpy_parameters(handle,
 *                      bzero, rzero, ntime, ntor,
 *                      ion_mass, z_ion, itor, numvar, version,
 *                      n0, L0, B0) -> None
 *
 * Total: 13 arguments (handle + 12 parameters).
 * Format string: "iddiiddiiiddd"
 *   i  = handle
 *   d  = bzero
 *   d  = rzero
 *   i  = ntime
 *   i  = ntor
 *   d  = ion_mass
 *   d  = z_ion
 *   i  = itor
 *   i  = numvar
 *   i  = version
 *   d  = n0
 *   d  = L0
 *   d  = B0
 */
PyObject* fio_set_numpy_parameters_py(PyObject* self, PyObject* args)
{
    int handle, ntime, ntor, itor, numvar, version;
    double bzero, rzero, ion_mass, z_ion, n0, L0, B0;

    if(!PyArg_ParseTuple(args, "iddiiddiiiddd",
                         &handle,
                         &bzero, &rzero,
                         &ntime, &ntor,
                         &ion_mass, &z_ion,
                         &itor, &numvar, &version,
                         &n0, &L0, &B0))
        return nullptr;

    m3dc1_numpy_source* ns = get_numpy_source(handle);
    if(!ns) return nullptr;

    ns->bzero    = bzero;
    ns->rzero    = rzero;
    ns->ntime    = ntime;
    ns->ntor_global = ntor;
    ns->ion_mass = ion_mass;
    ns->z_ion    = z_ion;
    ns->itor     = itor;
    ns->numvar   = numvar;
    ns->version  = version;

    // Reproduce normalization from m3dc1_source::open() exactly.
    const double c   = 3.e10;
    const double m_p = 1.6726e-24;

    double B0_cgs = B0;
    double L0_cgs = L0;
    double n0_cgs = n0;

    ns->p0   = B0_cgs * B0_cgs / (4. * M_PI);
    ns->J0   = B0_cgs * c / (4. * M_PI * L0_cgs);
    ns->v0   = B0_cgs / std::sqrt(4. * M_PI * ion_mass * m_p * n0_cgs);
    ns->t0   = L0_cgs / ns->v0;
    ns->Phi0 = L0_cgs * ns->v0 * B0_cgs / c;
    ns->temp0 = ns->p0 / n0_cgs;
    ns->eta0  = 4. * M_PI * ns->t0 * ns->v0 * ns->v0 / (c * c);

    // Convert to MKS.
    ns->n0   = n0_cgs  / 1.e-6;
    ns->B0   = B0_cgs  / 1e4;
    ns->L0   = L0_cgs  / 100.;
    ns->p0  /= 10.;
    ns->J0  /= (c * 1e-5);
    ns->v0  /= 100.;
    ns->Phi0 /= (1.e8 / c);
    ns->temp0 /= 1.6022e-12;
    ns->eta0  /= 1e11 / (c * c);

    if(itor == 1)
        ns->period = 2. * M_PI;
    else
        ns->period = 2. * M_PI * rzero * ns->L0;

    ns->ion_species = fio_species((int)ion_mass, 1, 0);

    Py_RETURN_NONE;
}

/*
 * set_numpy_mesh(handle, nelms,
 *                a, b, c, co, sn, x, z, bound, region,
 *                nplanes, period, toroidal) -> None
 */
PyObject* fio_set_numpy_mesh_py(PyObject* self, PyObject* args)
{
    int handle, nelms, nplanes, toroidal;
    double period;
    Py_buffer a_buf, b_buf, c_buf, co_buf, sn_buf, x_buf, z_buf;
    Py_buffer bound_buf, region_buf;

    if(!PyArg_ParseTuple(args, "iiy*y*y*y*y*y*y*y*y*idi",
                         &handle, &nelms,
                         &a_buf, &b_buf, &c_buf,
                         &co_buf, &sn_buf,
                         &x_buf, &z_buf,
                         &bound_buf, &region_buf,
                         &nplanes, &period, &toroidal))
        return nullptr;

    m3dc1_numpy_source* ns = get_numpy_source(handle);
    if(!ns) {
        PyBuffer_Release(&a_buf);   PyBuffer_Release(&b_buf);
        PyBuffer_Release(&c_buf);   PyBuffer_Release(&co_buf);
        PyBuffer_Release(&sn_buf);  PyBuffer_Release(&x_buf);
        PyBuffer_Release(&z_buf);   PyBuffer_Release(&bound_buf);
        PyBuffer_Release(&region_buf);
        return nullptr;
    }

    // Replace existing mesh.
    delete ns->mesh;
    ns->mesh = new m3dc1_mesh(nelms);

    const double* a      = (const double*)a_buf.buf;
    const double* b      = (const double*)b_buf.buf;
    const double* c      = (const double*)c_buf.buf;
    const double* co     = (const double*)co_buf.buf;
    const double* sn     = (const double*)sn_buf.buf;
    const double* xv     = (const double*)x_buf.buf;
    const double* zv     = (const double*)z_buf.buf;
    const int*    bound  = (const int*)bound_buf.buf;
    const int*    region = (const int*)region_buf.buf;

    for(int i = 0; i < nelms; i++) {
        ns->mesh->a[i]      = a[i];
        ns->mesh->b[i]      = b[i];
        ns->mesh->c[i]      = c[i];
        ns->mesh->co[i]     = co[i];
        ns->mesh->sn[i]     = sn[i];
        ns->mesh->x[i]      = xv[i];
        ns->mesh->z[i]      = zv[i];
        ns->mesh->bound[i]  = bound[i];
        ns->mesh->region[i] = region[i];
    }
    ns->mesh->nplanes  = nplanes;
    ns->mesh->period   = period;
    ns->mesh->toroidal = (toroidal != 0);

    PyBuffer_Release(&a_buf);   PyBuffer_Release(&b_buf);
    PyBuffer_Release(&c_buf);   PyBuffer_Release(&co_buf);
    PyBuffer_Release(&sn_buf);  PyBuffer_Release(&x_buf);
    PyBuffer_Release(&z_buf);   PyBuffer_Release(&bound_buf);
    PyBuffer_Release(&region_buf);

    // Allocate neighbor arrays and compute connectivity.
    ns->mesh->nneighbors = new int[nelms];
    ns->mesh->neighbor   = new int*[nelms];
    for(int i = 0; i < nelms; i++) {
        ns->mesh->nneighbors[i] = 0;
        ns->mesh->neighbor[i]   = new int[ns->mesh->max_neighbors()];
    }
    ns->mesh->find_neighbors();

    Py_RETURN_NONE;
}

/*
 * set_numpy_field(handle, timeslice, field_name,
 *                 data_real, data_imag_or_None, ntor) -> None
 */
PyObject* fio_set_numpy_field_py(PyObject* self, PyObject* args)
{
    int handle, timeslice, ntor;
    const char* field_name;
    Py_buffer real_buf;
    PyObject* imag_obj;

    if(!PyArg_ParseTuple(args, "iisy*Oi",
                         &handle, &timeslice, &field_name,
                         &real_buf, &imag_obj, &ntor))
        return nullptr;

    m3dc1_numpy_source* ns = get_numpy_source(handle);
    if(!ns) {
        PyBuffer_Release(&real_buf);
        return nullptr;
    }

    m3dc1_numpy_field_data fd;
    fd.ntor       = ntor;
    fd.is_complex = (imag_obj != Py_None);

    size_t n = real_buf.len / sizeof(double);
    fd.real.resize(n);
    std::memcpy(fd.real.data(), real_buf.buf, real_buf.len);
    PyBuffer_Release(&real_buf);

    if(fd.is_complex) {
        Py_buffer imag_buf;
        if(PyObject_GetBuffer(imag_obj, &imag_buf, PyBUF_SIMPLE) < 0)
            return nullptr;
        size_t ni = imag_buf.len / sizeof(double);
        fd.imag.resize(ni);
        std::memcpy(fd.imag.data(), imag_buf.buf, imag_buf.len);
        PyBuffer_Release(&imag_buf);
    }

    auto key = m3dc1_numpy_field_key(timeslice, std::string(field_name));
    ns->field_store[key] = std::move(fd);

    Py_RETURN_NONE;
}

/*
 * set_numpy_scalar(handle, scalar_name, values_bytes) -> None
 */
PyObject* fio_set_numpy_scalar_py(PyObject* self, PyObject* args)
{
    int handle;
    const char* scalar_name;
    Py_buffer val_buf;

    if(!PyArg_ParseTuple(args, "isy*", &handle, &scalar_name, &val_buf))
        return nullptr;

    m3dc1_numpy_source* ns = get_numpy_source(handle);
    if(!ns) {
        PyBuffer_Release(&val_buf);
        return nullptr;
    }

    m3dc1_numpy_scalar sc;
    size_t n = val_buf.len / sizeof(double);
    sc.values.resize(n);
    std::memcpy(sc.values.data(), val_buf.buf, val_buf.len);
    PyBuffer_Release(&val_buf);

    ns->scalar_store[std::string(scalar_name)] = std::move(sc);

    Py_RETURN_NONE;
}
