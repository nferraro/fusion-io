
#include "fusion_io.h"
#include <Python.h>

#include <iostream>

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

  static PyMethodDef fio_methods[] = {
    {"add_field", fio_add_field_py, METH_VARARGS, ""},
    {"allocate_hint", fio_allocate_hint_py, METH_VARARGS, ""},
    {"close_field", fio_close_field_py, METH_VARARGS, ""},
    {"close_series", fio_close_series_py, METH_VARARGS, ""},
    {"close_source", fio_close_source_py, METH_VARARGS, ""},
    {"create_compound_field", fio_create_compound_field_py, METH_VARARGS, ""},
    {"eval_scalar_field", fio_eval_scalar_field_py, METH_VARARGS, ""},
    {"eval_scalar_field_deriv", fio_eval_scalar_field_deriv_py, METH_VARARGS, ""},
    {"eval_series", fio_eval_series_py, METH_VARARGS, ""},
    {"eval_vector_field", fio_eval_vector_field_py, METH_VARARGS, ""},
    {"eval_vector_field_deriv", fio_eval_vector_field_deriv_py, METH_VARARGS, ""},
    {"eval_tensor_field", fio_eval_tensor_field_py, METH_VARARGS, ""},
    {"free_hint", fio_free_hint_py, METH_VARARGS, ""},
    {"get_available_fields", fio_get_available_fields_py, METH_VARARGS, ""},
    {"get_options", fio_get_options_py, METH_VARARGS, ""},
    {"get_field", fio_get_field_py, METH_VARARGS, ""},
    {"get_field_name", fio_get_field_name_py, METH_VARARGS, ""},
    {"get_option_name", fio_get_option_name_py, METH_VARARGS, ""},
    {"get_series", fio_get_series_py, METH_VARARGS, ""},
    {"get_series_bounds", fio_get_series_py, METH_VARARGS, ""},
    {"open_source", fio_open_source_py, METH_VARARGS, ""},
    {"set_int_option", fio_set_int_option_py, METH_VARARGS, ""},
    {"set_str_option", fio_set_str_option_py, METH_VARARGS, ""},
    {"set_real_option", fio_set_real_option_py, METH_VARARGS, ""},
    {"get_int_parameter", fio_get_int_parameter_py, METH_VARARGS, ""},
    {"get_real_parameter", fio_get_real_parameter_py, METH_VARARGS, ""},
    {"get_real_field_parameter", fio_get_real_field_parameter_py, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
  };

  static PyObject* fio_error;

#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef fio_py =
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
    m = PyModule_Create(&fio_py);
#else
    m = Py_InitModule("fio_py", fio_methods);
#endif

    if(m==NULL)
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
}

PyObject* fio_add_field_py(PyObject* self, PyObject *args)
{
  int icfield;
  int ifield;
  int op;
  double fac;

  if(!PyArg_ParseTuple(args, "iiid", &icfield, &ifield, &op, &fac))
    return NULL;

  int ierr = fio_add_field(icfield, ifield, op, fac);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_allocate_hint_py(PyObject* self, PyObject *args)
{
  int isrc;

  if(!PyArg_ParseTuple(args, "i", &isrc))
    return NULL;

  int sz = fio_sizeof_search_hint(isrc);

  if(sz==0)
    return NULL;

  void* hint = calloc(1, sz);
  return PyLong_FromVoidPtr(hint);
}

PyObject* fio_free_hint_py(PyObject* self, PyObject *args)
{
  PyObject* hint_py;
  void* hint;

  if(!PyArg_ParseTuple(args, "O", &hint_py))
    return NULL;

  hint = PyLong_AsVoidPtr(hint_py);

  free(hint);

  Py_RETURN_NONE;
}

PyObject* fio_close_field_py(PyObject* self, PyObject *args)
{
  int ifield;

  if(!PyArg_ParseTuple(args, "i", &ifield))
    return NULL;

  int ierr = fio_close_field(ifield);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_close_source_py(PyObject* self, PyObject *args)
{
  int isrc;

  if(!PyArg_ParseTuple(args, "i", &isrc))
    return NULL;

  int ierr = fio_close_source(isrc);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_close_series_py(PyObject* self, PyObject *args)
{
  int isrc;

  if(!PyArg_ParseTuple(args, "i", &isrc))
    return NULL;

  int ierr = fio_close_series(isrc);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_create_compound_field_py(PyObject* self, PyObject *args)
{
  int ifield;

  int ierr = fio_create_compound_field(&ifield);
  if(ierr != FIO_SUCCESS)
    return NULL;
  return Py_BuildValue("i", ifield);
}

PyObject* fio_eval_scalar_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;

  if(!PyArg_ParseTuple(args, "i(ddd)O", &ifield, &(x[0]), &(x[1]), &(x[2]), &hint_py))
    return NULL;

  void* hint = PyLong_AsVoidPtr(hint_py);

  double v;
  int ierr = fio_eval_field(ifield, x, &v, hint);
  if(ierr != FIO_SUCCESS)
    /* return the buildvalue of an empty string translates to python's None!*/
    return Py_BuildValue("");
    /*return NULL; */

  hint_py = PyLong_FromVoidPtr(hint);

  return Py_BuildValue("d", v);
}

PyObject* fio_eval_scalar_field_deriv_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;

  if(!PyArg_ParseTuple(args, "i(ddd)O", &ifield, &(x[0]), &(x[1]), &(x[2]), &hint_py))
    return NULL;

  void* hint = PyLong_AsVoidPtr(hint_py);

  double v[3];
  /*int ierr = fio_eval_field(ifield, x, &v, hint);*/
  int ierr = fio_eval_field_deriv(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS)
    /* return the buildvalue of an empty string translates to python's None!*/
    return Py_BuildValue("(2)", Py_None);
    /*return NULL; */

  hint_py = PyLong_FromVoidPtr(hint);

  return Py_BuildValue("(ddd)", v[0], v[1], v[2]);
}

PyObject* fio_eval_series_py(PyObject* self, PyObject *args)
{
  int iseries;
  double x;

  if(!PyArg_ParseTuple(args, "id", &iseries, &x))
    return NULL;

  double v;
  int ierr = fio_eval_series(iseries, x, &v);
  if(ierr != FIO_SUCCESS)
    /* return the buildvalue of an empty string translates to python's None!*/
    return Py_BuildValue("");
    /*return NULL; */

  return Py_BuildValue("d", v);
}

PyObject* fio_eval_vector_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;

  if(!PyArg_ParseTuple(args, "i(ddd)O", &ifield, &(x[0]), &(x[1]), &(x[2]), &hint_py))
    return NULL;

  void* hint = PyLong_AsVoidPtr(hint_py);

  double v[3];
  /* Do the evaluation!*/
  int ierr = fio_eval_field(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS)
    /* return a 3-tuple of nulls if it couldn't read*/
    return Py_BuildValue("(2)", Py_None);

  hint_py = PyLong_FromVoidPtr(hint);

  return Py_BuildValue("(ddd)", v[0], v[1], v[2]);
}

PyObject* fio_eval_vector_field_deriv_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;

  if(!PyArg_ParseTuple(args, "i(ddd)O", &ifield, &(x[0]), &(x[1]), &(x[2]), &hint_py))
    return NULL;

  void* hint = PyLong_AsVoidPtr(hint_py);

  double v[9];
  /* Do the evaluation!*/
  int ierr = fio_eval_field_deriv(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS)
    /* return a 9-tuple of nulls if it couldn't read*/
    return Py_BuildValue("(8)", Py_None);

  hint_py = PyLong_FromVoidPtr(hint);

  return Py_BuildValue("(ddddddddd)", v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
}

PyObject* fio_eval_tensor_field_py(PyObject* self, PyObject *args)
{
  int ifield;
  double x[3];
  PyObject* hint_py;

  if(!PyArg_ParseTuple(args, "i(ddd)O", &ifield, &(x[0]), &(x[1]), &(x[2]), &hint_py))
    return NULL;

  void* hint = PyLong_AsVoidPtr(hint_py);

  double v[9];
  /* Do the evaluation!*/
  int ierr = fio_eval_field(ifield, x, v, hint);
  if(ierr != FIO_SUCCESS)
    /* return a 9-tuple of nulls if it couldn't read*/
    return Py_BuildValue("(8)", Py_None);

  hint_py = PyLong_FromVoidPtr(hint);

  return Py_BuildValue("(ddddddddd)", v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
}

PyObject* fio_get_available_fields_py(PyObject* self, PyObject *args)
{
  int isrc;

  if(!PyArg_ParseTuple(args, "i", &isrc))
    return NULL;

  field_type* f;
  int n;
  int ierr = fio_get_available_fields(isrc, &n, &f);
  if(ierr != FIO_SUCCESS)
    return NULL;

  PyObject* list = PyList_New(n);
  for(int i=0; i<n; i++)
#if PY_MAJOR_VERSION >= 3
    PyList_SET_ITEM(list, i, PyLong_FromLong(f[i]));
#else
    PyList_SET_ITEM(list, i, PyInt_FromLong(f[i]));
#endif
  return list;
}

PyObject* fio_get_field_py(PyObject* self, PyObject *args)
{
  int isrc;
  int itype;
  int handle;

  if(!PyArg_ParseTuple(args, "ii", &isrc, &itype))
    return NULL;

  int ierr = fio_get_field(isrc, itype, &handle);
  if(ierr != FIO_SUCCESS)
    return NULL;
  return Py_BuildValue("i", handle);
}

PyObject* fio_get_field_name_py(PyObject* self, PyObject *args)
{
  int itype;

  if(!PyArg_ParseTuple(args, "i", &itype))
    return NULL;

  std::string name;
  int ierr = fio_get_field_name(itype, &name);
  if(ierr != FIO_SUCCESS)
    return NULL;

  return Py_BuildValue("s", name.c_str());
}

PyObject* fio_get_option_name_py(PyObject* self, PyObject *args)
{
  int iopt;

  if(!PyArg_ParseTuple(args, "i", &iopt))
    return NULL;

  std::string name;
  int ierr = fio_get_option_name(iopt, &name);
  if(ierr != FIO_SUCCESS)
    return NULL;

  return Py_BuildValue("s", name.c_str());
}

PyObject* fio_get_options_py(PyObject* self, PyObject *args)
{
  int isrc;
  if(!PyArg_ParseTuple(args, "i", &isrc))
    return NULL;

  int ierr = fio_get_options(isrc);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_get_int_parameter_py(PyObject* self, PyObject *args)
{
  int isrc, t;
  int p;

  if(!PyArg_ParseTuple(args, "ii", &isrc, &t))
    return NULL;

  int ierr = fio_get_int_parameter(isrc, t, &p);
  if(ierr != FIO_SUCCESS)
    return NULL;

  return Py_BuildValue("i", p);
}

PyObject* fio_get_real_parameter_py(PyObject* self, PyObject *args)
{
  int isrc, t;
  double p;

  if(!PyArg_ParseTuple(args, "ii", &isrc, &t))
    return NULL;

  int ierr = fio_get_real_parameter(isrc, t, &p);
  if(ierr != FIO_SUCCESS)
    return NULL;

  return Py_BuildValue("d", p);
}

PyObject* fio_get_real_field_parameter_py(PyObject* self, PyObject *args)
{
  int ifield, t;
  double p;

  if(!PyArg_ParseTuple(args, "ii", &ifield, &t))
    return NULL;

  int ierr = fio_get_real_field_parameter(ifield, t, &p);
  if(ierr != FIO_SUCCESS)
    return NULL;

  return Py_BuildValue("d", p);
}

PyObject* fio_get_series_py(PyObject* self, PyObject *args)
{
  int isrc;
  int itype;
  int handle;

  if(!PyArg_ParseTuple(args, "ii", &isrc, &itype))
    return NULL;

  int ierr = fio_get_series(isrc, itype, &handle);
  if(ierr != FIO_SUCCESS)
    return NULL;
  return Py_BuildValue("i", handle);
}

PyObject* fio_get_series_bounds_py(PyObject* self, PyObject *args)
{
  int iseries;
  double tmin, tmax;

  if(!PyArg_ParseTuple(args, "i", &iseries, &iseries))
    return NULL;

  int ierr = fio_get_series_bounds(iseries, &tmin, &tmax);
  if(ierr != FIO_SUCCESS)
    return NULL;
  return Py_BuildValue("dd", tmin, tmax);
}

PyObject* fio_open_source_py(PyObject* self, PyObject *args)
{
  int type;
  char* filename;
  if(!PyArg_ParseTuple(args, "is", &type, &filename))
    return NULL;

  int handle;
  int ierr = fio_open_source(type, filename, &handle);
  if(ierr != FIO_SUCCESS)
    return NULL;
  return Py_BuildValue("i", handle);
}

PyObject* fio_set_int_option_py(PyObject* self, PyObject *args)
{
  int iopt;
  int v;
  if(!PyArg_ParseTuple(args, "ii", &iopt, &v))
    return NULL;

  int ierr = fio_set_int_option(iopt, v);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_set_str_option_py(PyObject* self, PyObject *args)
{
  int iopt;
  char v[256];
  if(!PyArg_ParseTuple(args, "is", &iopt, &v))
    return NULL;

  int ierr = fio_set_str_option(iopt, v);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}

PyObject* fio_set_real_option_py(PyObject* self, PyObject *args)
{
  int iopt;
  double v;
  if(!PyArg_ParseTuple(args, "id", &iopt, &v))
    return NULL;

  int ierr = fio_set_real_option(iopt, v);
  if(ierr != FIO_SUCCESS)
    return NULL;
  Py_RETURN_NONE;
}
