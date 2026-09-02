#include "m3dc1_numpy_source.h"
#include "fusion_io_defs.h"
#include "options.h"

#include <cmath>
#include <cstring>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------------------------
// m3dc1_numpy_series
// ---------------------------------------------------------------------------

int m3dc1_numpy_series::eval(double t, double* v)
{
    auto it = src->scalar_store.find(name);
    if(it == src->scalar_store.end()) return FIO_NOT_FOUND;

    auto tit = src->scalar_store.find("slice_times");
    if(tit == src->scalar_store.end()) return FIO_NOT_FOUND;

    const auto& times = tit->second.values;
    const auto& vals  = it->second.values;

    if(times.empty() || vals.empty()) return FIO_NOT_FOUND;

    // Find index of nearest time value.
    int idx = 0;
    double best = std::abs(times[0] - t);
    for(int i = 1; i < (int)times.size(); i++) {
        double diff = std::abs(times[i] - t);
        if(diff < best) { best = diff; idx = i; }
    }
    if(idx >= (int)vals.size()) return FIO_OUT_OF_BOUNDS;
    *v = vals[idx] * scale;
    return FIO_SUCCESS;
}

int m3dc1_numpy_series::bounds(double* tmin, double* tmax) const
{
    auto tit = src->scalar_store.find("slice_times");
    if(tit == src->scalar_store.end()) {
        *tmin = *tmax = 0.;
        return FIO_NOT_FOUND;
    }
    const auto& times = tit->second.values;
    if(times.empty()) {
        *tmin = *tmax = 0.;
        return FIO_NOT_FOUND;
    }
    *tmin = times.front();
    *tmax = times.back();
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_fio_field
// ---------------------------------------------------------------------------

int m3dc1_numpy_fio_field::eval(const double* x, double* v, fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    bool ok = owned_field->eval(x[0], x[1], x[2],
                                m3dc1_field::GET_VAL, val, element);
    if(!ok) return FIO_OUT_OF_BOUNDS;
    for(int i = 0; i < dim; i++)
        v[i] = val[i] * scale;
    return FIO_SUCCESS;
}

int m3dc1_numpy_fio_field::eval_deriv(const double* x, double* v, fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    bool ok = owned_field->eval(x[0], x[1], x[2],
                                m3dc1_field::GET_DVAL, val, element);
    if(!ok) return FIO_OUT_OF_BOUNDS;
    if(dim == 1) {
        v[0] = val[m3dc1_field::OP_DR] * scale;
        v[1] = 0.;
        v[2] = val[m3dc1_field::OP_DZ] * scale;
    } else {
        return FIO_UNSUPPORTED;
    }
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_source
// ---------------------------------------------------------------------------

m3dc1_numpy_source::m3dc1_numpy_source()
    : mesh(nullptr),
      bzero(0.), rzero(0.),
      n0(0.), L0(0.), B0(0.), p0(0.), J0(0.), v0(0.), t0(0.),
      Phi0(0.), temp0(0.), eta0(0.),
      ntime(0), ntor_global(0), itor(1), numvar(0), version(0),
      ion_mass(1.0), z_ion(1.0), kprad_z(0),
      period(2. * M_PI), ifprime(0)
{}

m3dc1_numpy_source::~m3dc1_numpy_source()
{
    delete mesh;
}

int m3dc1_numpy_source::allocate_search_hint(fio_hint* s)
{
    *s = new int;
    *((int*)(*s)) = -1;
    return FIO_SUCCESS;
}

int m3dc1_numpy_source::deallocate_search_hint(fio_hint* s)
{
    delete (int*)(*s);
    return FIO_SUCCESS;
}

int m3dc1_numpy_source::get_field_options(fio_option_list* opt) const
{
    opt->clear();
    opt->add_option(FIO_TIMESLICE,    0);
    opt->add_option(FIO_LINEAR_SCALE, 1.);
    opt->add_option(FIO_PHASE,        0.);
    opt->add_option(FIO_PART,         FIO_TOTAL);
    opt->add_option(FIO_SPECIES,      0);
    opt->add_option(FIO_FIELD_NAME,   std::string("P"));
    return FIO_SUCCESS;
}

int m3dc1_numpy_source::get_available_fields(fio_field_list* fields) const
{
    fields->clear();
    fields->push_back(FIO_DENSITY);
    fields->push_back(FIO_CURRENT_DENSITY);
    fields->push_back(FIO_FLUID_VELOCITY);
    fields->push_back(FIO_MAGNETIC_FIELD);
    fields->push_back(FIO_POLOIDAL_FLUX);
    fields->push_back(FIO_PRESSURE);
    fields->push_back(FIO_TEMPERATURE);
    fields->push_back(FIO_TOTAL_PRESSURE);
    fields->push_back(FIO_TOTAL_RADIATION);
    fields->push_back(FIO_VECTOR_POTENTIAL);
    fields->push_back(FIO_ELECTRIC_FIELD);
    fields->push_back(FIO_SCALAR_FIELD);
    return FIO_SUCCESS;
}

int m3dc1_numpy_source::get_int_parameter(const parameter_type t, int* p) const
{
    switch(t) {
    case FIO_NUM_TIMESLICES:
        *p = ntime;
        return FIO_SUCCESS;
    case FIO_GEOMETRY:
        *p = (itor == 1 ? FIO_CYLINDRICAL : FIO_CARTESIAN);
        return FIO_SUCCESS;
    case FIO_TOROIDAL_MODE:
        *p = ntor_global;
        return FIO_SUCCESS;
    default:
        return FIO_UNSUPPORTED;
    }
}

int m3dc1_numpy_source::get_real_parameter(const parameter_type t,
                                            double* p) const
{
    switch(t) {
    case FIO_PERIOD:
        *p = period;
        return FIO_SUCCESS;
    default:
        return FIO_UNSUPPORTED;
    }
}

int m3dc1_numpy_source::get_slice_time(int slice, double* time) const
{
    auto it = scalar_store.find("slice_times");
    if(it == scalar_store.end()) return FIO_NOT_FOUND;
    const auto& times = it->second.values;
    if(slice < 0 || slice >= (int)times.size()) return FIO_OUT_OF_BOUNDS;
    *time = times[slice] * t0;
    return FIO_SUCCESS;
}

/*
 * make_field
 * ~~~~~~~~~~
 * Create an m3dc1_field (or m3dc1_complex_field) from the injected data
 * for the given field name and timeslice.
 * Returns nullptr if the data has not been injected yet.
 * The caller takes ownership of the returned pointer.
 */
m3dc1_field* m3dc1_numpy_source::make_field(const std::string& name,
                                             int timeslice) const
{
    if(!mesh) return nullptr;

    auto key = m3dc1_numpy_field_key(timeslice, name);
    auto it = field_store.find(key);
    if(it == field_store.end()) return nullptr;

    const m3dc1_numpy_field_data& fd = it->second;

    if(fd.is_complex && !fd.imag.empty()) {
        m3dc1_complex_field* cf = new m3dc1_complex_field(mesh, fd.ntor);
        std::memcpy(cf->data,   fd.real.data(),
                    fd.real.size() * sizeof(double));
        std::memcpy(cf->data_i, fd.imag.data(),
                    fd.imag.size() * sizeof(double));
        return cf;
    } else {
        m3dc1_field* f = new m3dc1_field(mesh);
        std::memcpy(f->data, fd.real.data(),
                    fd.real.size() * sizeof(double));
        return f;
    }
}

/*
 * get_field
 * ~~~~~~~~~
 * Mirrors m3dc1_source::get_field() but uses make_field() instead of
 * m3dc1_file::load_field(). Normalization constants must have been set
 * via set_numpy_parameters() before this is called.
 */
int m3dc1_numpy_source::get_field(const field_type t, fio_field** f,
                                   const fio_option_list* opt)
{
    *f = nullptr;
    if(!mesh) {
        std::cerr << "m3dc1_numpy_source::get_field: mesh not injected"
                  << std::endl;
        return FIO_FILE_ERROR;
    }

    int timeslice = 0;
    opt->get_option(FIO_TIMESLICE, &timeslice);

    int s = 0;
    opt->get_option(FIO_SPECIES, &s);
    if(s == FIO_MAIN_ION) s = (int)ion_species;

    std::string hdf5_name;
    double scale = 1.0;
    int dim = 1;

    switch(t) {
    case FIO_DENSITY:
        hdf5_name = (s == (int)fio_electron) ? "ne" : "den";
        scale = n0; dim = 1;
        break;

    case FIO_TEMPERATURE:
        hdf5_name = (s == (int)fio_electron) ? "te" : "ti";
        scale = temp0; dim = 1;
        break;

    case FIO_POLOIDAL_FLUX:
        hdf5_name = "psi";
        scale = 2. * M_PI * B0 * L0 * L0;
        dim = 1;
        break;

    case FIO_PRESSURE:
        hdf5_name = (s == (int)fio_electron) ? "Pe" : "pi";
        scale = p0; dim = 1;
        break;

    case FIO_TOTAL_PRESSURE:
        hdf5_name = "P"; scale = p0; dim = 1;
        break;

    case FIO_TOTAL_RADIATION:
        hdf5_name = "kprad_rad";
        scale = (t0 > 0.) ? p0 / t0 : 1.;
        dim = 1;
        break;

    case FIO_RESISTIVITY:
        hdf5_name = "eta"; scale = eta0; dim = 1;
        break;

    case FIO_MAGNETIC_FIELD:
        hdf5_name = "B"; scale = B0; dim = 3;
        break;

    case FIO_VECTOR_POTENTIAL:
        hdf5_name = "A"; scale = B0 * L0; dim = 3;
        break;

    case FIO_ELECTRIC_FIELD:
        hdf5_name = "E"; scale = B0 * v0; dim = 3;
        break;

    case FIO_CURRENT_DENSITY:
        hdf5_name = "J"; scale = J0; dim = 3;
        break;

    case FIO_FLUID_VELOCITY:
        hdf5_name = "v"; scale = v0; dim = 3;
        break;

    case FIO_SCALAR_FIELD: {
        std::string field_name;
        opt->get_option(FIO_FIELD_NAME, &field_name);
        hdf5_name = field_name;
        scale = 1.0; dim = 1;
        break;
    }

    default:
        return FIO_UNSUPPORTED;
    }

    m3dc1_field* mf = make_field(hdf5_name, timeslice);
    if(!mf) {
        std::cerr << "m3dc1_numpy_source::get_field: field '"
                  << hdf5_name << "' not injected for timeslice "
                  << timeslice << std::endl;
        return FIO_NOT_FOUND;
    }

    *f = new m3dc1_numpy_fio_field(mf, scale, dim);
    return FIO_SUCCESS;
}

int m3dc1_numpy_source::get_series(const series_type t, fio_series** s)
{
    *s = nullptr;
    switch(t) {
    case FIO_MAGAXIS_PSI:
        *s = new m3dc1_numpy_series(this, "psimin",    B0 * L0 * L0);
        return FIO_SUCCESS;
    case FIO_LCFS_PSI:
        *s = new m3dc1_numpy_series(this, "psi_lcfs",  B0 * L0 * L0);
        return FIO_SUCCESS;
    case FIO_MAGAXIS_R:
        *s = new m3dc1_numpy_series(this, "xmag",      L0);
        return FIO_SUCCESS;
    case FIO_MAGAXIS_Z:
        *s = new m3dc1_numpy_series(this, "zmag",      L0);
        return FIO_SUCCESS;
    default:
        return FIO_UNSUPPORTED;
    }
}
