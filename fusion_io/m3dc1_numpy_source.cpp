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
// m3dc1_numpy_fio_field  (simple scalar wrapper)
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
// m3dc1_numpy_vector_potential  eval
// ---------------------------------------------------------------------------
// Mirrors m3dc1_vector_potential::eval() from m3dc1_field.cpp:
//   A_R   =  R * (df/dZ)
//   A_Z   = -R * (df/dR)
//   A_Phi =  psi / R
// All converted to MKS by multiplying by scale = B0*L0.

int m3dc1_numpy_vector_potential::eval(const double* x, double* v,
                                        fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    double r = (itor == 1) ? x[0] : 1.;

    // A_Phi = psi / R
    if(!psi->eval(x[0], x[1], x[2], m3dc1_field::GET_VAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[1] = val[m3dc1_field::OP_1] / r;

    if(f) {
        // A_R =  R * (df/dZ)
        // A_Z = -R * (df/dR)
        if(!f->eval(x[0], x[1], x[2], m3dc1_field::GET_DVAL, val, element))
            return FIO_OUT_OF_BOUNDS;
        v[0] =  r * val[m3dc1_field::OP_DZ];
        v[2] = -r * val[m3dc1_field::OP_DR];
    } else {
        v[0] = v[2] = 0.;
    }

    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_vector_potential  eval_deriv
// ---------------------------------------------------------------------------
// Mirrors m3dc1_vector_potential::eval_deriv() from m3dc1_field.cpp exactly.
//
// Unit conversion (from m3dc1_field.cpp):
//   v[DR_*]   *= B0       (spatial R derivative, dimension 1/L0 cancels)
//   v[DPHI_*] *= B0*L0    (toroidal derivative)
//   v[DZ_*]   *= B0       (spatial Z derivative)
// With scale = B0*L0 and scale_L = L0:
//   B0 = scale / scale_L
//   B0*L0 = scale

int m3dc1_numpy_vector_potential::eval_deriv(const double* x, double* v,
                                              fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];

    const m3dc1_field::m3dc1_get_op psiget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_VAL |
                                    m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL);

    const m3dc1_field::m3dc1_get_op fget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL |
                                    m3dc1_field::GET_DDVAL);

    // Initialise all 9 output components to zero.
    for(int i = 0; i < 9; i++) v[i] = 0.;

    if(!psi->eval(x[0], x[1], x[2], psiget, val, element))
        return FIO_OUT_OF_BOUNDS;

    // dA_Phi/dR = (dpsi/dR - psi/R) / R
    v[FIO_DR_PHI  ] = (val[m3dc1_field::OP_DR] -
                       val[m3dc1_field::OP_1] / x[0]) / x[0];
    // dA_Phi/dphi = dpsi/dphi / R
    v[FIO_DPHI_PHI] = val[m3dc1_field::OP_DP] / x[0];
    // dA_Phi/dZ = dpsi/dZ / R
    v[FIO_DZ_PHI  ] = val[m3dc1_field::OP_DZ] / x[0];

    if(f) {
        if(!f->eval(x[0], x[1], x[2], fget, val, element))
            return FIO_OUT_OF_BOUNDS;

        // dA_R/dR =  R * d2f/dRdZ + df/dZ
        v[FIO_DR_R  ] =  x[0] * val[m3dc1_field::OP_DRZ] +
                         val[m3dc1_field::OP_DZ];
        // dA_Z/dR = -(R * d2f/dR2 + df/dR)
        v[FIO_DR_Z  ] = -(x[0] * val[m3dc1_field::OP_DRR] +
                          val[m3dc1_field::OP_DR]);
        // dA_R/dphi =  R * d2f/dZdphi
        v[FIO_DPHI_R] =  x[0] * val[m3dc1_field::OP_DZP];
        // dA_Z/dphi = -R * d2f/dRdphi
        v[FIO_DPHI_Z] = -x[0] * val[m3dc1_field::OP_DRP];
        // dA_R/dZ =  R * d2f/dZ2
        v[FIO_DZ_R  ] =  x[0] * val[m3dc1_field::OP_DZZ];
        // dA_Z/dZ = -R * d2f/dRdZ
        v[FIO_DZ_Z  ] = -x[0] * val[m3dc1_field::OP_DRZ];
    }

    // Unit conversion — mirrors m3dc1_vector_potential::eval_deriv():
    //   spatial (R,Z) derivatives: *= B0  = scale / scale_L
    //   toroidal (phi) derivatives: *= B0*L0 = scale
    double B0 = scale / scale_L;
    v[FIO_DR_R    ] *= B0;
    v[FIO_DR_PHI  ] *= B0;
    v[FIO_DR_Z    ] *= B0;
    v[FIO_DPHI_R  ] *= scale;
    v[FIO_DPHI_PHI] *= scale;
    v[FIO_DPHI_Z  ] *= scale;
    v[FIO_DZ_R    ] *= B0;
    v[FIO_DZ_PHI  ] *= B0;
    v[FIO_DZ_Z    ] *= B0;

    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_magnetic_field  eval
// ---------------------------------------------------------------------------
// Mirrors m3dc1_magnetic_field::eval() from m3dc1_field.cpp:
//   B_R   = -(dpsi/dZ) / R
//   B_Z   =  (dpsi/dR) / R
//   B_Phi =  I / R
// Plus f/fp correction for 3D/complex runs.
// All converted to MKS by multiplying by scale = B0.

int m3dc1_numpy_magnetic_field::eval(const double* x, double* v,
                                      fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    double r = (itor == 1) ? x[0] : 1.;

    // B from psi
    if(!psi->eval(x[0], x[1], x[2], m3dc1_field::GET_DVAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[0] = -val[m3dc1_field::OP_DZ] / r;
    v[2] =  val[m3dc1_field::OP_DR] / r;

    // B_Phi from I
    if(!I->eval(x[0], x[1], x[2], m3dc1_field::GET_VAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[1] = val[m3dc1_field::OP_1] / r;

    // f/fp correction (only for 3D/complex runs where use_f is true)
    if(use_f && f) {
        const m3dc1_field::m3dc1_get_op fget =
            (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_DVAL |
                                        m3dc1_field::GET_PVAL);
        if(!f->eval(x[0], x[1], x[2], fget, val, element))
            return FIO_OUT_OF_BOUNDS;
        if(ifprime == 0) {
            v[0] -= val[m3dc1_field::OP_DRP];
            v[2] -= val[m3dc1_field::OP_DZP];
        } else {
            v[0] -= val[m3dc1_field::OP_DR];
            v[2] -= val[m3dc1_field::OP_DZ];
        }
    }

    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_magnetic_field  eval_deriv
// ---------------------------------------------------------------------------
// Mirrors m3dc1_magnetic_field::eval_deriv() from m3dc1_field.cpp exactly.
//
// Unit conversion (from m3dc1_field.cpp):
//   v[DR_*]   *= B0 / L0   (spatial R derivative)
//   v[DPHI_*] *= B0        (toroidal derivative, cylindrical itor==1)
//   v[DPHI_*] *= B0 / L0   (toroidal derivative, slab itor==0)
//   v[DZ_*]   *= B0 / L0   (spatial Z derivative)
// With scale = B0 and scale_L = L0:
//   B0/L0 = scale / scale_L

int m3dc1_numpy_magnetic_field::eval_deriv(const double* x, double* v,
                                            fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    double r = (itor == 1) ? x[0] : 1.;

    const m3dc1_field::m3dc1_get_op psiget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL |
                                    m3dc1_field::GET_DDVAL);

    const m3dc1_field::m3dc1_get_op gget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_VAL |
                                    m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL);

    const m3dc1_field::m3dc1_get_op fget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL |
                                    m3dc1_field::GET_DDVAL |
                                    m3dc1_field::GET_PPVAL);

    // Initialise all 9 output components to zero.
    for(int i = 0; i < 9; i++) v[i] = 0.;

    // --- psi contribution ---
    if(!psi->eval(x[0], x[1], x[2], psiget, val, element))
        return FIO_OUT_OF_BOUNDS;

    v[FIO_DR_R  ] = -val[m3dc1_field::OP_DRZ] / r;
    if(itor == 1)
        v[FIO_DR_R] += val[m3dc1_field::OP_DZ] / (r * r);
    v[FIO_DPHI_R] = -val[m3dc1_field::OP_DZP] / r;
    v[FIO_DZ_R  ] = -val[m3dc1_field::OP_DZZ] / r;

    v[FIO_DR_Z  ] =  val[m3dc1_field::OP_DRR] / r;
    if(itor == 1)
        v[FIO_DR_Z] -= val[m3dc1_field::OP_DR] / (r * r);
    v[FIO_DPHI_Z] =  val[m3dc1_field::OP_DRP] / r;
    v[FIO_DZ_Z  ] =  val[m3dc1_field::OP_DRZ] / r;

    // --- I contribution ---
    if(!I->eval(x[0], x[1], x[2], gget, val, element))
        return FIO_OUT_OF_BOUNDS;

    v[FIO_DR_PHI  ] =  val[m3dc1_field::OP_DR] / r;
    if(itor == 1)
        v[FIO_DR_PHI] -= val[m3dc1_field::OP_1] / (r * r);
    v[FIO_DPHI_PHI] =  val[m3dc1_field::OP_DP] / r;
    v[FIO_DZ_PHI  ] =  val[m3dc1_field::OP_DZ] / r;

    // --- f/fp contribution ---
    if(use_f && f) {
        if(!f->eval(x[0], x[1], x[2], fget, val, element))
            return FIO_OUT_OF_BOUNDS;

        if(ifprime == 0) {
            v[FIO_DR_R  ] -= val[m3dc1_field::OP_DRRP];
            v[FIO_DPHI_R] -= val[m3dc1_field::OP_DRPP];
            v[FIO_DZ_R  ] -= val[m3dc1_field::OP_DRZP];

            v[FIO_DR_Z  ] -= val[m3dc1_field::OP_DRZP];
            v[FIO_DPHI_Z] -= val[m3dc1_field::OP_DZPP];
            v[FIO_DZ_Z  ] -= val[m3dc1_field::OP_DZZP];
        } else {
            v[FIO_DR_R  ] -= val[m3dc1_field::OP_DRR];
            v[FIO_DPHI_R] -= val[m3dc1_field::OP_DRP];
            v[FIO_DZ_R  ] -= val[m3dc1_field::OP_DRZ];

            v[FIO_DR_Z  ] -= val[m3dc1_field::OP_DRZ];
            v[FIO_DPHI_Z] -= val[m3dc1_field::OP_DZP];
            v[FIO_DZ_Z  ] -= val[m3dc1_field::OP_DZZ];
        }
    }

    // --- Unit conversion ---
    // Mirrors m3dc1_magnetic_field::eval_deriv() in m3dc1_field.cpp:
    //   spatial (R,Z) derivatives: *= B0 / L0  = scale / scale_L
    //   toroidal (phi) derivatives: *= B0       (itor==1)
    //                            or *= B0 / L0  (itor==0)
    double B0overL0 = scale / scale_L;
    v[FIO_DR_R    ] *= B0overL0;
    v[FIO_DR_PHI  ] *= B0overL0;
    v[FIO_DR_Z    ] *= B0overL0;
    if(itor == 1) {
        v[FIO_DPHI_R  ] *= scale;
        v[FIO_DPHI_PHI] *= scale;
        v[FIO_DPHI_Z  ] *= scale;
    } else {
        v[FIO_DPHI_R  ] *= B0overL0;
        v[FIO_DPHI_PHI] *= B0overL0;
        v[FIO_DPHI_Z  ] *= B0overL0;
    }
    v[FIO_DZ_R    ] *= B0overL0;
    v[FIO_DZ_PHI  ] *= B0overL0;
    v[FIO_DZ_Z    ] *= B0overL0;

    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_current_density
// ---------------------------------------------------------------------------

int m3dc1_numpy_current_density::eval(const double* x, double* v,
                                       fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    double r = (itor == 1) ? x[0] : 1.;

    const m3dc1_field::m3dc1_get_op psiget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL |
                                    m3dc1_field::GET_DDVAL);
    const m3dc1_field::m3dc1_get_op gget = m3dc1_field::GET_DVAL;
    const m3dc1_field::m3dc1_get_op fget =
        (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_DVAL |
                                    m3dc1_field::GET_PVAL |
                                    m3dc1_field::GET_PPVAL);

    if(!psi->eval(x[0], x[1], x[2], psiget, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[0] =  val[m3dc1_field::OP_DRP] / (r * r);
    v[2] =  val[m3dc1_field::OP_DZP] / (r * r);
    v[1] = -(val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ]) / r;
    if(itor == 1) v[1] += val[m3dc1_field::OP_DR] / (r * r);

    if(!I->eval(x[0], x[1], x[2], gget, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[0] -= val[m3dc1_field::OP_DZ] / r;
    v[2] += val[m3dc1_field::OP_DR] / r;

    if(use_f && f) {
        if(!f->eval(x[0], x[1], x[2], fget, val, element))
            return FIO_OUT_OF_BOUNDS;
        if(ifprime == 0) {
            v[0] -= val[m3dc1_field::OP_DZPP] / r;
            v[2] += val[m3dc1_field::OP_DRPP] / r;
        } else {
            v[0] -= val[m3dc1_field::OP_DZP] / r;
            v[2] += val[m3dc1_field::OP_DRP] / r;
        }
    }

    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_fluid_velocity
// ---------------------------------------------------------------------------

int m3dc1_numpy_fluid_velocity::eval(const double* x, double* v,
                                      fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];
    double r = (itor == 1) ? x[0] : 1.;

    if(!phi->eval(x[0], x[1], x[2], m3dc1_field::GET_DVAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[0] = -val[m3dc1_field::OP_DZ] * r;
    v[2] =  val[m3dc1_field::OP_DR] * r;

    if(!V->eval(x[0], x[1], x[2], m3dc1_field::GET_VAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[1] = val[m3dc1_field::OP_1] * r;

    if(!chi->eval(x[0], x[1], x[2], m3dc1_field::GET_DVAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    v[0] += val[m3dc1_field::OP_DR] / (r * r);
    v[2] += val[m3dc1_field::OP_DZ] / (r * r);

    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_electric_field
// ---------------------------------------------------------------------------

int m3dc1_numpy_electric_field::eval(const double* x, double* v,
                                      fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];

    m3dc1_field* comps[3] = { E_R, E_PHI, E_Z };
    for(int i = 0; i < 3; i++) {
        if(!comps[i]->eval(x[0], x[1], x[2],
                           m3dc1_field::GET_VAL, val, element))
            return FIO_OUT_OF_BOUNDS;
        v[i] = val[m3dc1_field::OP_1] * scale;
    }
    return FIO_SUCCESS;
}

// ---------------------------------------------------------------------------
// m3dc1_numpy_ion_pressure  (pi = P - Pe)
// ---------------------------------------------------------------------------

int m3dc1_numpy_ion_pressure::eval(const double* x, double* v,
                                    fio_hint hint)
{
    int* element = (int*)hint;
    double val[m3dc1_field::OP_NUM];

    if(!P->eval(x[0], x[1], x[2], m3dc1_field::GET_VAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    double pval = val[m3dc1_field::OP_1];

    if(!Pe->eval(x[0], x[1], x[2], m3dc1_field::GET_VAL, val, element))
        return FIO_OUT_OF_BOUNDS;
    double peval = val[m3dc1_field::OP_1];

    *v = (pval - peval) * scale;
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

    bool use_f = (timeslice != -1);

    switch(t) {

    case FIO_DENSITY: {
        std::string nm = (s == (int)fio_electron) ? "ne" : "den";
        m3dc1_field* mf = make_field(nm, timeslice);
        if(!mf) {
            std::cerr << "m3dc1_numpy_source::get_field: field '"
                      << nm << "' not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        *f = new m3dc1_numpy_fio_field(mf, n0, 1);
        return FIO_SUCCESS;
    }

    case FIO_TEMPERATURE: {
        std::string nm = (s == (int)fio_electron) ? "te" : "ti";
        m3dc1_field* mf = make_field(nm, timeslice);
        if(!mf) {
            std::cerr << "m3dc1_numpy_source::get_field: field '"
                      << nm << "' not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        *f = new m3dc1_numpy_fio_field(mf, temp0, 1);
        return FIO_SUCCESS;
    }

    case FIO_POLOIDAL_FLUX: {
        m3dc1_field* mf = make_field("psi", timeslice);
        if(!mf) {
            std::cerr << "m3dc1_numpy_source::get_field: field 'psi'"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        *f = new m3dc1_numpy_fio_field(mf, 2. * M_PI * B0 * L0 * L0, 1);
        return FIO_SUCCESS;
    }

    case FIO_TOTAL_PRESSURE: {
        m3dc1_field* mf = make_field("P", timeslice);
        if(!mf) {
            std::cerr << "m3dc1_numpy_source::get_field: field 'P'"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        *f = new m3dc1_numpy_fio_field(mf, p0, 1);
        return FIO_SUCCESS;
    }

    case FIO_TOTAL_RADIATION: {
        m3dc1_field* mf = make_field("kprad_rad", timeslice);
        if(!mf) {
            std::cerr << "m3dc1_numpy_source::get_field: field 'kprad_rad'"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        double sc = (t0 > 0.) ? p0 / t0 : 1.;
        *f = new m3dc1_numpy_fio_field(mf, sc, 1);
        return FIO_SUCCESS;
    }

    case FIO_PRESSURE: {
        if(s == (int)fio_electron) {
            m3dc1_field* mf = make_field("Pe", timeslice);
            if(!mf) {
                std::cerr << "m3dc1_numpy_source::get_field: field 'Pe'"
                          << " not injected for timeslice "
                          << timeslice << std::endl;
                return FIO_NOT_FOUND;
            }
            *f = new m3dc1_numpy_fio_field(mf, p0, 1);
        } else {
            m3dc1_field* mP  = make_field("P",  timeslice);
            m3dc1_field* mPe = make_field("Pe", timeslice);
            if(!mP || !mPe) {
                delete mP;
                delete mPe;
                std::cerr << "m3dc1_numpy_source::get_field: fields 'P' or 'Pe'"
                          << " not injected for timeslice "
                          << timeslice << std::endl;
                return FIO_NOT_FOUND;
            }
            *f = new m3dc1_numpy_ion_pressure(mP, mPe, p0);
        }
        return FIO_SUCCESS;
    }

    case FIO_VECTOR_POTENTIAL: {
        m3dc1_field* mPsi = make_field("psi", timeslice);
        if(!mPsi) {
            std::cerr << "m3dc1_numpy_source::get_field: field 'psi'"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        m3dc1_field* mF = nullptr;
        if(use_f) {
            mF = make_field("f", timeslice);
        }
        // Pass B0*L0 as scale and L0 as scale_L.
        *f = new m3dc1_numpy_vector_potential(mPsi, mF, B0 * L0, L0, itor);
        return FIO_SUCCESS;
    }

    case FIO_MAGNETIC_FIELD: {
        m3dc1_field* mPsi = make_field("psi", timeslice);
        m3dc1_field* mI   = make_field("I",   timeslice);
        if(!mPsi || !mI) {
            delete mPsi;
            delete mI;
            std::cerr << "m3dc1_numpy_source::get_field: fields 'psi' or 'I'"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        m3dc1_field* mF = nullptr;
        if(use_f) {
            mF = make_field(ifprime ? "fp" : "f", timeslice);
        }
        // Pass B0 as scale and L0 as scale_L.
        *f = new m3dc1_numpy_magnetic_field(mPsi, mI, mF, B0, L0, itor,
                                             ifprime, use_f && mF != nullptr);
        return FIO_SUCCESS;
    }

    case FIO_CURRENT_DENSITY: {
        m3dc1_field* mPsi = make_field("psi", timeslice);
        m3dc1_field* mI   = make_field("I",   timeslice);
        if(!mPsi || !mI) {
            delete mPsi;
            delete mI;
            std::cerr << "m3dc1_numpy_source::get_field: fields 'psi' or 'I'"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        m3dc1_field* mF = nullptr;
        if(use_f) {
            mF = make_field(ifprime ? "fp" : "f", timeslice);
        }
        *f = new m3dc1_numpy_current_density(mPsi, mI, mF, J0, itor,
                                              ifprime, use_f && mF != nullptr);
        return FIO_SUCCESS;
    }

    case FIO_FLUID_VELOCITY: {
        m3dc1_field* mPhi = make_field("phi", timeslice);
        m3dc1_field* mV   = make_field("V",   timeslice);
        m3dc1_field* mChi = make_field("chi", timeslice);
        if(!mPhi || !mV || !mChi) {
            delete mPhi;
            delete mV;
            delete mChi;
            std::cerr << "m3dc1_numpy_source::get_field: fields 'phi', 'V',"
                      << " or 'chi' not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        *f = new m3dc1_numpy_fluid_velocity(mPhi, mV, mChi, v0, itor);
        return FIO_SUCCESS;
    }

    case FIO_ELECTRIC_FIELD: {
        m3dc1_field* mER   = make_field("E_R",   timeslice);
        m3dc1_field* mEPHI = make_field("E_PHI", timeslice);
        m3dc1_field* mEZ   = make_field("E_Z",   timeslice);
        if(!mER || !mEPHI || !mEZ) {
            delete mER;
            delete mEPHI;
            delete mEZ;
            std::cerr << "m3dc1_numpy_source::get_field: E-field components"
                      << " not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        double sc = (L0 > 0.) ? Phi0 / L0 : 1.;
        *f = new m3dc1_numpy_electric_field(mER, mEPHI, mEZ, sc);
        return FIO_SUCCESS;
    }

    case FIO_SCALAR_FIELD: {
        std::string field_name;
        opt->get_option(FIO_FIELD_NAME, &field_name);
        m3dc1_field* mf = make_field(field_name, timeslice);
        if(!mf) {
            std::cerr << "m3dc1_numpy_source::get_field: field '"
                      << field_name << "' not injected for timeslice "
                      << timeslice << std::endl;
            return FIO_NOT_FOUND;
        }
        *f = new m3dc1_numpy_fio_field(mf, 1.0, 1);
        return FIO_SUCCESS;
    }

    default:
        return FIO_UNSUPPORTED;
    }
}

int m3dc1_numpy_source::get_series(const series_type t, fio_series** s)
{
    *s = nullptr;
    switch(t) {
    case FIO_MAGAXIS_PSI:
        *s = new m3dc1_numpy_series(this, "psimin",   B0 * L0 * L0);
        return FIO_SUCCESS;
    case FIO_LCFS_PSI:
        *s = new m3dc1_numpy_series(this, "psi_lcfs", B0 * L0 * L0);
        return FIO_SUCCESS;
    case FIO_MAGAXIS_R:
        *s = new m3dc1_numpy_series(this, "xmag",     L0);
        return FIO_SUCCESS;
    case FIO_MAGAXIS_Z:
        *s = new m3dc1_numpy_series(this, "zmag",     L0);
        return FIO_SUCCESS;
    default:
        return FIO_UNSUPPORTED;
    }
}
