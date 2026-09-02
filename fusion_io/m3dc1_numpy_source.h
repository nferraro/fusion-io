#ifndef M3DC1_NUMPY_SOURCE_H
#define M3DC1_NUMPY_SOURCE_H

/*
 * m3dc1_numpy_source.h
 *
 * An fio_source implementation that is populated entirely from numpy arrays
 * injected from Python, rather than by reading an HDF5 file.
 *
 * This allows the Python layer to stream only the field coefficient arrays
 * it needs from a remote host (via SSH + h5py), inject them here, and then
 * use the full fusion-io finite element evaluation machinery without
 * transferring complete HDF5 files.
 *
 * Workflow from Python:
 *   handle = fio_py.open_numpy_source()
 *   fio_py.set_numpy_parameters(handle, bzero, rzero, ntime, ntor,
 *                               ion_mass, z_ion, itor, numvar, version,
 *                               n0, L0, B0)
 *   fio_py.set_numpy_mesh(handle, nelms,
 *                         a, b, c, co, sn, x, z, bound, region,
 *                         nplanes, period, toroidal)
 *   fio_py.set_numpy_scalar(handle, 'slice_times', times_bytes)
 *   fio_py.set_numpy_field(handle, timeslice, 'psi',
 *                          data_real_bytes, data_imag_bytes_or_None, ntor)
 *   # Then use exactly as a normal source:
 *   fio_py.get_options(handle)
 *   fio_py.set_int_option(fio_py.FIO_TIMESLICE, t)
 *   ifield = fio_py.get_field(handle, fio_py.FIO_POLOIDAL_FLUX)
 *   val = fio_py.eval_scalar_field(ifield, (R, phi, Z), hint)
 */

#include "fusion_io_source.h"
#include "fusion_io_field.h"
#include "fusion_io_series.h"
#include "fusion_io_species.h"
#include "m3dc1_field.h"
#include "m3dc1_mesh.h"

#include <map>
#include <string>
#include <vector>
#include <cmath>

// Key for the field map: (timeslice, field_name)
typedef std::pair<int, std::string> m3dc1_numpy_field_key;

// One stored field: real coefficients, optional imaginary, ntor for complex
struct m3dc1_numpy_field_data {
    std::vector<double> real;   // size nelms * nbasis
    std::vector<double> imag;   // empty if not complex
    int ntor;
    bool is_complex;
};

// One stored scalar time series
struct m3dc1_numpy_scalar {
    std::vector<double> values;
};

class m3dc1_numpy_source;

/*
 * fio_series backed by a scalar array injected into an m3dc1_numpy_source.
 * Evaluates by finding the nearest time in 'slice_times'.
 */
class m3dc1_numpy_series : public fio_series {
public:
    const m3dc1_numpy_source* src;
    std::string name;
    double scale;

    m3dc1_numpy_series(const m3dc1_numpy_source* s,
                       const std::string& n, double sc)
        : src(s), name(n), scale(sc) {}
    virtual ~m3dc1_numpy_series() {}

    virtual int eval(double t, double* v);
    virtual int bounds(double* tmin, double* tmax) const;
};

/*
 * fio_field wrapper that owns an m3dc1_field* created from injected numpy
 * data and evaluates it using the standard m3dc1_field::eval().
 */
class m3dc1_numpy_fio_field : public fio_field {
public:
    m3dc1_field* owned_field;   // owned by this object
    double scale;
    int dim;

    m3dc1_numpy_fio_field(m3dc1_field* f, double sc, int d)
        : owned_field(f), scale(sc), dim(d) {}
    virtual ~m3dc1_numpy_fio_field() { delete owned_field; }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return dim; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0);
};

class m3dc1_numpy_source : public fio_source {
public:
    // Mesh (shared across all time slices, owned by this object)
    m3dc1_mesh* mesh;

    // Field coefficient storage: (timeslice, name) -> data
    std::map<m3dc1_numpy_field_key, m3dc1_numpy_field_data> field_store;

    // Scalar time series storage: name -> values
    std::map<std::string, m3dc1_numpy_scalar> scalar_store;

    // Physical parameters (mirrors m3dc1_source members)
    double bzero, rzero;
    double n0, L0, B0, p0, J0, v0, t0, Phi0, temp0, eta0;
    int ntime, ntor_global, itor, numvar, version;
    double ion_mass, z_ion;
    int kprad_z;
    double period;
    fio_species ion_species;
    int ifprime;

    m3dc1_numpy_source();
    virtual ~m3dc1_numpy_source();

    // fio_source interface
    virtual int open(const char*) { return FIO_SUCCESS; }
    virtual int close()           { return FIO_SUCCESS; }
    virtual int get_field_options(fio_option_list*) const;
    virtual int get_available_fields(fio_field_list*) const;
    virtual int get_int_parameter(const parameter_type, int*) const;
    virtual int get_real_parameter(const parameter_type, double*) const;
    virtual int get_field(const field_type, fio_field**,
                          const fio_option_list*);
    virtual int get_series(const series_type, fio_series**);
    virtual int allocate_search_hint(fio_hint*);
    virtual int deallocate_search_hint(fio_hint*);
    virtual int sizeof_search_hint() const { return sizeof(int); }

    // Internal helpers
    m3dc1_field* make_field(const std::string& name, int timeslice) const;
    int get_slice_time(int slice, double* time) const;
};

#endif // M3DC1_NUMPY_SOURCE_H
