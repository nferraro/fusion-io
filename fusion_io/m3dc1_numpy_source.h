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
 * Used for simple scalar fields (psi, ne, te, ti, P, Pe, den, ...).
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

/*
 * Composite fio_field for the vector potential A.
 * Assembles A from injected 'psi' and 'f' component fields, using the
 * same formulas as m3dc1_vector_potential::eval() in m3dc1_field.cpp:
 *
 *   A_R   =  R * (df/dZ)
 *   A_Z   = -R * (df/dR)
 *   A_Phi =  psi / R
 *
 * Owns the component m3dc1_field* objects.
 */
class m3dc1_numpy_vector_potential : public fio_field {
public:
    m3dc1_field* psi;   // owned
    m3dc1_field* f;     // owned, may be nullptr
    double scale;       // B0 * L0
    int itor;           // 1 for cylindrical

    m3dc1_numpy_vector_potential(m3dc1_field* psi_, m3dc1_field* f_,
                                  double sc, int itor_)
        : psi(psi_), f(f_), scale(sc), itor(itor_) {}
    virtual ~m3dc1_numpy_vector_potential()
    {
        delete psi;
        delete f;
    }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return 3; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0)
    {
        return FIO_UNSUPPORTED;
    }
};

/*
 * Composite fio_field for the magnetic field B.
 * Assembles B from injected 'psi', 'I', and 'f'/'fp' component fields,
 * using the same formulas as m3dc1_magnetic_field::eval():
 *
 *   B_R   = -(dpsi/dZ) / R
 *   B_Z   =  (dpsi/dR) / R
 *   B_Phi =  I / R
 *
 * (The f/fp correction terms apply only for 3D/complex runs.)
 * Owns the component m3dc1_field* objects.
 */
class m3dc1_numpy_magnetic_field : public fio_field {
public:
    m3dc1_field* psi;   // owned
    m3dc1_field* I;     // owned
    m3dc1_field* f;     // owned, may be nullptr
    double scale;       // B0
    int itor;
    int ifprime;        // 0 -> use f, 1 -> use fp
    bool use_f;         // true if f is available

    m3dc1_numpy_magnetic_field(m3dc1_field* psi_, m3dc1_field* I_,
                                m3dc1_field* f_, double sc, int itor_,
                                int ifprime_, bool use_f_)
        : psi(psi_), I(I_), f(f_), scale(sc), itor(itor_),
          ifprime(ifprime_), use_f(use_f_) {}
    virtual ~m3dc1_numpy_magnetic_field()
    {
        delete psi;
        delete I;
        delete f;
    }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return 3; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0)
    {
        return FIO_UNSUPPORTED;
    }
};

/*
 * Composite fio_field for the current density J.
 * Assembles J from injected 'psi', 'I', and 'f'/'fp' component fields,
 * using the same formulas as m3dc1_current_density::eval().
 * Owns the component m3dc1_field* objects.
 */
class m3dc1_numpy_current_density : public fio_field {
public:
    m3dc1_field* psi;
    m3dc1_field* I;
    m3dc1_field* f;
    double scale;       // J0
    int itor;
    int ifprime;
    bool use_f;

    m3dc1_numpy_current_density(m3dc1_field* psi_, m3dc1_field* I_,
                                 m3dc1_field* f_, double sc, int itor_,
                                 int ifprime_, bool use_f_)
        : psi(psi_), I(I_), f(f_), scale(sc), itor(itor_),
          ifprime(ifprime_), use_f(use_f_) {}
    virtual ~m3dc1_numpy_current_density()
    {
        delete psi;
        delete I;
        delete f;
    }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return 3; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0)
    {
        return FIO_UNSUPPORTED;
    }
};

/*
 * Composite fio_field for the fluid velocity v.
 * Assembles v from injected 'phi', 'V', 'chi' component fields,
 * using the same formulas as m3dc1_fluid_velocity::eval():
 *
 *   V_R   = -(dphi/dZ)*R + dchi/dR / R^2
 *   V_Z   =  (dphi/dR)*R + dchi/dZ / R^2
 *   V_Phi =  V * R
 *
 * Owns the component m3dc1_field* objects.
 */
class m3dc1_numpy_fluid_velocity : public fio_field {
public:
    m3dc1_field* phi;
    m3dc1_field* V;
    m3dc1_field* chi;
    double scale;       // v0
    int itor;

    m3dc1_numpy_fluid_velocity(m3dc1_field* phi_, m3dc1_field* V_,
                                m3dc1_field* chi_, double sc, int itor_)
        : phi(phi_), V(V_), chi(chi_), scale(sc), itor(itor_) {}
    virtual ~m3dc1_numpy_fluid_velocity()
    {
        delete phi;
        delete V;
        delete chi;
    }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return 3; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0)
    {
        return FIO_UNSUPPORTED;
    }
};

/*
 * Composite fio_field for the electric field E.
 * Assembles E from injected 'E_R', 'E_PHI', 'E_Z' component fields,
 * using the same formulas as m3dc1_electric_field::eval().
 * Owns the component m3dc1_field* objects.
 */
class m3dc1_numpy_electric_field : public fio_field {
public:
    m3dc1_field* E_R;
    m3dc1_field* E_PHI;
    m3dc1_field* E_Z;
    double scale;       // Phi0 / L0

    m3dc1_numpy_electric_field(m3dc1_field* er, m3dc1_field* ep,
                                m3dc1_field* ez, double sc)
        : E_R(er), E_PHI(ep), E_Z(ez), scale(sc) {}
    virtual ~m3dc1_numpy_electric_field()
    {
        delete E_R;
        delete E_PHI;
        delete E_Z;
    }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return 3; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0)
    {
        return FIO_UNSUPPORTED;
    }
};

/*
 * Composite fio_field for ion pressure pi = P - Pe.
 * Owns the component m3dc1_field* objects.
 */
class m3dc1_numpy_ion_pressure : public fio_field {
public:
    m3dc1_field* P;
    m3dc1_field* Pe;
    double scale;       // p0

    m3dc1_numpy_ion_pressure(m3dc1_field* P_, m3dc1_field* Pe_, double sc)
        : P(P_), Pe(Pe_), scale(sc) {}
    virtual ~m3dc1_numpy_ion_pressure()
    {
        delete P;
        delete Pe;
    }

    virtual fio_field* clone() const { return nullptr; }
    virtual int dimension() const { return 1; }
    virtual int eval(const double* x, double* v, fio_hint hint = 0);
    virtual int eval_deriv(const double* x, double* v, fio_hint hint = 0)
    {
        return FIO_UNSUPPORTED;
    }
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
