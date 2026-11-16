#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <compound_field.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <cmath>
#include <deque>
#include <algorithm>
#include <mpi.h>
#include <netcdf.h> 

// Helper to check NetCDF errors
#define NC_CHECK(e) { if((e) != NC_NOERR) { \
  std::cerr << "NetCDF Error: " << nc_strerror(e) << std::endl; \
  MPI_Abort(MPI_COMM_WORLD, 2); } }

// ========================================================
// Global Variables & Configuration
// ========================================================
int nR = 20;      // number of radial points
int nZ = 20;      // number of Z points
int nphi = 18;    // number of toroidal points
double Rmin = 5.0, Rmax = 10.0;
double Zmin = -4.0, Zmax = 4.0;
int ibootstrap = 0;
int ivecpot = 0;

std::deque<fio_source*> sources;
fio_compound_field electron_density, electron_temperature, pressure;
fio_compound_field current_density, vector_potential, mag;
fio_compound_field ion_density, ion_temperature, JpdotB;
fio_compound_field JpdotB_dndpsi, JpdotB_dtedpsi, JpdotB_dtidpsi;
fio_compound_field JBS_r, JBS_phi, JBS_z;
fio_compound_field eta_fld;

// MPI Globals
int world_rank = 0;
int world_size = 1;

void print_usage();
int process_command_line(int argc, char* argv[]);
int process_line(const std::string& opt, const int argc, const std::string argv[]);
void delete_sources();
bool create_source(const int type, const int argc, const std::string argv[]);

// ========================================================
// MAIN
// ========================================================
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if(argc < 2) {
    if(world_rank == 0) print_usage();
    MPI_Finalize();
    return 1;
  }
  
  if(process_command_line(argc, argv) != FIO_SUCCESS) {
    if(world_rank == 0) print_usage();
    MPI_Finalize();
    return 1;
  }

  if(sources.size() == 0) {
    if(world_rank == 0) {
      std::cerr << "Error, no source is loaded" << std::endl;
      print_usage();
    }
    MPI_Finalize();
    return 1;
  }

  if(world_rank == 0) {
    std::cerr << "Input parameters\n=======================\n";
    std::cerr << "Radial grid points (nR) = " << nR << " [" << Rmin << ", " << Rmax << "]\n";
    std::cerr << "Z grid points (nZ)      = " << nZ << " [" << Zmin << ", " << Zmax << "]\n";
    std::cerr << "Toroidal grid points (nphi) = " << nphi << '\n';
    std::cerr << "Bootstrap output = " << ibootstrap << '\n';
    std::cerr << "Vector Potential output = " << ivecpot << '\n';
    std::cerr << "=======================" << std::endl;
  }

  void* h;
  sources[0]->allocate_search_hint(&h);

  // -------------------------------------------------------------------------
  // ULTRA-OPTIMIZED PARALLEL FIELD EVALUATION 
  // -------------------------------------------------------------------------
  if(world_rank == 0) std::cerr << "Evaluating fields on Eulerian grid..." << std::endl;

  long long total_work = (long long)nR * nphi * nZ;
  long long pts_per_rank = total_work / world_size;
  long long rem_pts = total_work % world_size;

  long long my_pts_start = world_rank * pts_per_rank + std::min((long long)world_rank, rem_pts);
  long long my_pts_count = pts_per_rank + (world_rank < rem_pts ? 1 : 0);
  long long my_pts_end = my_pts_start + my_pts_count;

  // Local Buffers: Only allocate what this rank needs
  double* loc_Jac = new double[my_pts_count]();
  double* loc_Mask = new double[my_pts_count](); // Computational Boundary Mask
  double* loc_BR = new double[my_pts_count]();
  double* loc_BPhi = new double[my_pts_count]();
  double* loc_BZ = new double[my_pts_count]();
  double* loc_JR = new double[my_pts_count]();
  double* loc_JPhi = new double[my_pts_count]();
  double* loc_JZ = new double[my_pts_count]();
  double* loc_VecPot_R = new double[my_pts_count]();
  double* loc_VecPot_Phi = new double[my_pts_count]();
  double* loc_VecPot_Z = new double[my_pts_count]();
  double* loc_Te = new double[my_pts_count]();
  double* loc_Ti = new double[my_pts_count]();
  double* loc_Ne = new double[my_pts_count]();
  double* loc_Ni = new double[my_pts_count]();
  double* loc_P = new double[my_pts_count]();
  double* loc_JpdotB = new double[my_pts_count]();
  double* loc_JpdotB_dn = new double[my_pts_count]();
  double* loc_JpdotB_dte = new double[my_pts_count]();
  double* loc_JpdotB_dti = new double[my_pts_count]();
  double* loc_JBSR   = new double[my_pts_count]();
  double* loc_JBSPhi = new double[my_pts_count]();
  double* loc_JBSZ   = new double[my_pts_count]();
  double* loc_Eta = new double[my_pts_count]();

  double dR = (nR > 1) ? (Rmax - Rmin) / (nR - 1) : 0;
  double dZ = (nZ > 1) ? (Zmax - Zmin) / (nZ - 1) : 0;
  double dPhi = 2.0 * M_PI / nphi;

  // =========================================================================
  // SNAKE-PATH CACHE OPTIMIZATION
  // Sort the local evaluation order to zig-zag through space.
  // This ensures the spatial jump is exactly 1 grid cell every single iteration,
  // making the fio_hint cache nearly 100% successful and eliminating global searches!
  // =========================================================================
  std::vector<long long> eval_order(my_pts_count);
  for(long long step = 0; step < my_pts_count; step++) {
      eval_order[step] = step;
  }

  std::sort(eval_order.begin(), eval_order.end(), [&](long long a, long long b) {
      long long idx_a = my_pts_start + a;
      long long idx_b = my_pts_start + b;
      
      int ia = idx_a % nR; int ja = (idx_a / nR) % nphi; int ka = idx_a / (nR * nphi);
      int ib = idx_b % nR; int jb = (idx_b / nR) % nphi; int kb = idx_b / (nR * nphi);
      
      if (ka != kb) return ka < kb;
      // Zig-zag Phi based on Z
      if (ka % 2 == 1) { if (ja != jb) return ja > jb; } else { if (ja != jb) return ja < jb; }
      // Zig-zag R based on Phi
      if (ja % 2 == 1) return ia > ib;
      return ia < ib;
  });

  long long count_processed = 0;
  long long print_interval = my_pts_count / 10; 
  if(print_interval == 0) print_interval = 1;

  for(long long step = 0; step < my_pts_count; step++) {
      if(count_processed % print_interval == 0) {
           std::cerr << "[Rank " << world_rank << "] Progress: " 
                     << (count_processed * 100 / my_pts_count) << "%" << std::flush << std::endl;
      }
      
      long long loc_idx = eval_order[step]; 
      long long idx = my_pts_start + loc_idx;

      int i = idx % nR;
      int j = (idx / nR) % nphi;
      int k = idx / (nR * nphi);

      double p[3];
      p[0] = Rmin + i * dR;      // R
      p[1] = j * dPhi;           // Phi
      p[2] = Zmin + k * dZ;      // Z

      double v3[3] = {0.0, 0.0, 0.0};
      loc_Jac[loc_idx] = p[0];

      // =========================================================================
      // OUT-OF-BOUNDS SKIP
      // Evaluate Mag Field first. If it fails, the point is outside the plasma mesh.
      // We instantly skip evaluating the other variables, saving massive CPU time.
      // =========================================================================
      if (mag.eval(p, v3, h) == FIO_SUCCESS) {
          
          loc_Mask[loc_idx] = 1.0; // Mark as inside the computational boundary
          loc_BR[loc_idx] = v3[0]; loc_BPhi[loc_idx] = v3[1]; loc_BZ[loc_idx] = v3[2];

          // Evaluate Current Density
          current_density.eval(p, v3, h);
          loc_JR[loc_idx] = v3[0]; loc_JPhi[loc_idx] = v3[1]; loc_JZ[loc_idx] = v3[2];
          
          if(ivecpot) {
              vector_potential.eval(p, v3, h);
              loc_VecPot_R[loc_idx] = v3[0]; loc_VecPot_Phi[loc_idx] = v3[1]; loc_VecPot_Z[loc_idx] = v3[2];
          }
          
          double v_ne = 0.0, v_ni = 0.0, v_ti = 0.0, v_te = 0.0, v_p = 0.0;
          electron_density.eval(p, &v_ne, h);      loc_Ne[loc_idx] = v_ne;
          ion_density.eval(p, &v_ni, h);           loc_Ni[loc_idx] = v_ni;
          ion_temperature.eval(p, &v_ti, h);       loc_Ti[loc_idx] = v_ti;
          electron_temperature.eval(p, &v_te, h);  loc_Te[loc_idx] = v_te;
          pressure.eval(p, &v_p, h);               loc_P[loc_idx] = v_p;
          double v_eta = 0.0;
          eta_fld.eval(p, &v_eta, h);              loc_Eta[loc_idx] = v_eta;
          
          if(ibootstrap >= 1) {
              double v_jb = 0.0, v_dn = 0.0, v_dte = 0.0, v_dti = 0.0;
              JpdotB.eval(p, &v_jb, h);          loc_JpdotB[loc_idx] = v_jb;
              JpdotB_dndpsi.eval(p, &v_dn, h);   loc_JpdotB_dn[loc_idx] = v_dn;
              JpdotB_dtedpsi.eval(p, &v_dte, h); loc_JpdotB_dte[loc_idx] = v_dte;
              JpdotB_dtidpsi.eval(p, &v_dti, h); loc_JpdotB_dti[loc_idx] = v_dti;
              double v_r=0., v_ph=0., v_zz=0.;
              JBS_r.eval(p,&v_r,h);    loc_JBSR[loc_idx]=v_r;
              JBS_phi.eval(p,&v_ph,h); loc_JBSPhi[loc_idx]=v_ph;
              JBS_z.eval(p,&v_zz,h);   loc_JBSZ[loc_idx]=v_zz;
          }
      }
      count_processed++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank == 0) std::cerr << "Eval done. Gathering data..." << std::endl;
  // -------------------------------------------------------------------------
  // PER-FIELD GATHER + WRITE (NetCDF and binary VTK), rank 0 holds ONE array
  // -------------------------------------------------------------------------
  struct Fld { const char* name; double* loc; int ncvar; };
  std::vector<Fld> flds = {
    {"Comp_Boundary", loc_Mask}, {"Jac", loc_Jac},
    {"B_R", loc_BR}, {"B_Phi", loc_BPhi}, {"B_Z", loc_BZ},
    {"J_R", loc_JR}, {"J_Phi", loc_JPhi}, {"J_Z", loc_JZ}
  };
  if(ivecpot) {
    flds.push_back({"VecPot_R", loc_VecPot_R});
    flds.push_back({"VecPot_Phi", loc_VecPot_Phi});
    flds.push_back({"VecPot_Z", loc_VecPot_Z});
  }
  flds.push_back({"te_3d", loc_Te}); flds.push_back({"ti_3d", loc_Ti});
  flds.push_back({"ne_3d", loc_Ne}); flds.push_back({"ni_3d", loc_Ni});
  flds.push_back({"p_3d", loc_P});   flds.push_back({"eta_3d", loc_Eta});
  if(ibootstrap >= 1) {
    flds.push_back({"JpdotB_3d", loc_JpdotB});
    flds.push_back({"JpdotB_dndpsi_3d", loc_JpdotB_dn});
    flds.push_back({"JpdotB_dtedpsi_3d", loc_JpdotB_dte});
    flds.push_back({"JpdotB_dtidpsi_3d", loc_JpdotB_dti});
    flds.push_back({"JBS_R", loc_JBSR});
    flds.push_back({"JBS_Phi", loc_JBSPhi});
    flds.push_back({"JBS_Z", loc_JBSZ});
  }

  std::vector<int> flat_counts(world_size), flat_displs(world_size);
  { long long d = 0;
    for(int r=0; r<world_size; r++) {
      long long c = pts_per_rank + (r < rem_pts ? 1 : 0);
      flat_counts[r] = (int)c; flat_displs[r] = (int)d; d += c;
    } }
  int my_pts = (int)my_pts_count;

  // big-endian float writer for legacy binary VTK
  auto write_be_float = [](std::ofstream& f, float v) {
    uint32_t u; std::memcpy(&u, &v, 4);
    u = __builtin_bswap32(u);
    f.write(reinterpret_cast<const char*>(&u), 4);
  };

  std::ofstream vtk; int ncid = -1;
  double* buf = nullptr;
  long long vtk_points = (long long)nR * (nphi + 1) * nZ;

  if(world_rank == 0) {
    buf = new double[total_work];

    // ---- NetCDF definitions ----
    NC_CHECK(nc_create("neo_eulerian.nc", NC_CLOBBER | NC_NETCDF4, &ncid));
    int nr_d, nphi_d, nz_d;
    NC_CHECK(nc_def_dim(ncid, "Z", nZ, &nz_d));
    NC_CHECK(nc_def_dim(ncid, "Phi", nphi, &nphi_d));
    NC_CHECK(nc_def_dim(ncid, "R", nR, &nr_d));
    int rid, phid, zid;
    NC_CHECK(nc_def_var(ncid, "R", NC_FLOAT, 1, &nr_d, &rid));
    NC_CHECK(nc_def_var(ncid, "Phi", NC_FLOAT, 1, &nphi_d, &phid));
    NC_CHECK(nc_def_var(ncid, "Z", NC_FLOAT, 1, &nz_d, &zid));
    int dims[3] = {nz_d, nphi_d, nr_d};
    for(auto& f : flds) {
      NC_CHECK(nc_def_var(ncid, f.name, NC_FLOAT, 3, dims, &f.ncvar));
      NC_CHECK(nc_def_var_deflate(ncid, f.ncvar, 1, 1, 1)); // shuffle+zlib level 1
    }
    NC_CHECK(nc_enddef(ncid));

    std::vector<double> rax(nR), phax(nphi), zax(nZ);
    for(int i=0;i<nR;i++)   rax[i]  = Rmin + i*dR;
    for(int j=0;j<nphi;j++) phax[j] = j*dPhi;
    for(int k=0;k<nZ;k++)   zax[k]  = Zmin + k*dZ;
    NC_CHECK(nc_put_var_double(ncid, rid, rax.data()));
    NC_CHECK(nc_put_var_double(ncid, phid, phax.data()));
    NC_CHECK(nc_put_var_double(ncid, zid, zax.data()));

    // ---- VTK header + points (BINARY) ----
    std::cerr << "Opening VTK (binary)..." << std::endl;
    vtk.open("neo_eulerian.vtk", std::ios::binary);
    vtk << "# vtk DataFile Version 3.0\nEulerian 3D fields\nBINARY\nDATASET STRUCTURED_GRID\n";
    vtk << "DIMENSIONS " << nR << " " << nphi + 1 << " " << nZ << "\n";
    vtk << "POINTS " << vtk_points << " float\n";
    for(int k=0;k<nZ;k++) for(int j=0;j<=nphi;j++) for(int i=0;i<nR;i++) {
      double r = Rmin + i*dR, p = j*dPhi, z = Zmin + k*dZ;
      write_be_float(vtk, (float)(r*cos(p)));
      write_be_float(vtk, (float)(r*sin(p)));
      write_be_float(vtk, (float)z);
    }
    vtk << "\nPOINT_DATA " << vtk_points << "\n";
  }

  for(auto& f : flds) {
    // all ranks participate; only rank 0 has a valid recv buffer
    MPI_Gatherv(f.loc, my_pts, MPI_DOUBLE, buf, flat_counts.data(),
                flat_displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(world_rank == 0) {
      NC_CHECK(nc_put_var_double(ncid, f.ncvar, buf));
      vtk << "SCALARS " << f.name << " float 1\nLOOKUP_TABLE default\n";
      for(int k=0;k<nZ;k++) for(int j=0;j<=nphi;j++) {
        int rj = (j == nphi) ? 0 : j;
        for(int i=0;i<nR;i++)
          write_be_float(vtk, (float)buf[i + rj*nR + (long long)k*nR*nphi]);
      }
      vtk << "\n";
      std::cerr << "Wrote " << f.name << std::endl;
    }
  }

  if(world_rank == 0) {
    vtk.close();
    NC_CHECK(nc_close(ncid));
    delete[] buf;
    std::cerr << "Successfully wrote output files!" << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);   // all ranks leave together
  if(world_rank == 0) std::cerr << "Done." << std::endl;

  delete[] loc_Jac; delete[] loc_Mask; delete[] loc_BR; delete[] loc_BPhi; delete[] loc_BZ;
  delete[] loc_JR; delete[] loc_JPhi; delete[] loc_JZ;
  delete[] loc_VecPot_R; delete[] loc_VecPot_Phi; delete[] loc_VecPot_Z;
  delete[] loc_Te; delete[] loc_Ti; delete[] loc_Ne; delete[] loc_Ni; delete[] loc_P; delete[] loc_Eta;
  delete[] loc_JpdotB; delete[] loc_JpdotB_dn; delete[] loc_JpdotB_dte; delete[] loc_JpdotB_dti;
  delete[] loc_JBSR; delete[] loc_JBSPhi; delete[] loc_JBSZ;
  sources[0]->deallocate_search_hint(&h);
  delete_sources();
  MPI_Finalize();
  return 0;
}

// ========================================================
// Source Handlers
// ========================================================
// Load a field into a compound field; warn on rank 0 if it can't be loaded.
static bool load_field(fio_source* src, int type, fio_option_list& fopt,
                       fio_compound_field& dst, fio_hint hint, const char* name)
{
  fio_field* f = nullptr;
  if(src->get_field(type, &f, &fopt) == FIO_SUCCESS) {
    dst.add_field(f, FIO_ADD, 1., hint);
    return true;
  }
  if(world_rank == 0)
    std::cerr << "WARNING: could not load field '" << name
              << "' (output will be zeros)" << std::endl;
  return false;
}


bool create_source(const int type, const int argc, const std::string argv[]) 
{
  fio_source* src;
  fio_option_list fopt;
  int result;

  if (type == FIO_M3DC1_SOURCE) {
    src = new m3dc1_source();
    if(argc >= 1) result = src->open(argv[0].c_str());
    else result = src->open("C1.h5");

    if(result != FIO_SUCCESS) {
      std::cerr << "Error opening file" << std::endl;
      delete(src);
      return false;
    }
    src->get_field_options(&fopt);
    fopt.set_option(FIO_PART, FIO_TOTAL);
  } else {
    return false;
  }

  if(argc >= 2) fopt.set_option(FIO_TIMESLICE, atoi(argv[1].c_str()));
  if(argc >= 3) fopt.set_option(FIO_LINEAR_SCALE, atof(argv[2].c_str()));
  if(argc >= 4) fopt.set_option(FIO_PHASE, atof(argv[3].c_str()) * M_PI / 180.);

  fio_hint hint;
  src->allocate_search_hint(&hint);

    // species-independent fields
  load_field(src, FIO_MAGNETIC_FIELD,   fopt, mag,              hint, "B");
  load_field(src, FIO_CURRENT_DENSITY,  fopt, current_density,  hint, "J");
  load_field(src, FIO_VECTOR_POTENTIAL, fopt, vector_potential, hint, "A");
  load_field(src, FIO_TOTAL_PRESSURE,   fopt, pressure,         hint, "P");
  load_field(src, FIO_RESISTIVITY,      fopt, eta_fld,          hint, "eta");

  // bootstrap fields
  load_field(src, FIO_JBS,         fopt, JpdotB,         hint, "JpdotB");
  load_field(src, FIO_JBS_dndpsi,  fopt, JpdotB_dndpsi,  hint, "JpdotB_dndpsi");
  load_field(src, FIO_JBS_dtedpsi, fopt, JpdotB_dtedpsi, hint, "JpdotB_dtedpsi");
  load_field(src, FIO_JBS_dtidpsi, fopt, JpdotB_dtidpsi, hint, "JpdotB_dtidpsi");
  load_field(src, FIO_JBS_R,       fopt, JBS_r,          hint, "Jp_BS_r");
  load_field(src, FIO_JBS_PHI,     fopt, JBS_phi,        hint, "Jp_BS_phi");
  load_field(src, FIO_JBS_Z,       fopt, JBS_z,          hint, "Jp_BS_z");

  // electron species
  fopt.set_option(FIO_SPECIES, FIO_ELECTRON);
  load_field(src, FIO_TEMPERATURE, fopt, electron_temperature, hint, "te");
  load_field(src, FIO_DENSITY,     fopt, electron_density,     hint, "ne");

  // main ion species
  fopt.set_option(FIO_SPECIES, FIO_MAIN_ION);
  load_field(src, FIO_DENSITY,     fopt, ion_density,     hint, "den");
  load_field(src, FIO_TEMPERATURE, fopt, ion_temperature, hint, "ti");

  sources.push_back(src);
  return FIO_SUCCESS;
}

int process_command_line(int argc, char* argv[])
{
  const int max_args = 4;
  const int num_opts = 10;
  std::string arg_list[num_opts] = { "-bootstrap", "-vecpot", "-m3dc1", "-nR", "-nZ", "-nphi", "-Rmin", "-Rmax", "-Zmin", "-Zmax"};
  std::string opt = "";
  std::string arg[max_args];
  int args = 0;
  bool is_opt;
  bool processed = true;
  int result;

  for(int i=1; i<argc; i++) {
    is_opt = false;
    for(int j=0; j<num_opts; j++) {
      if(argv[i] == arg_list[j]) { is_opt = true; break; }
    }
    
    if(is_opt) {     
      if(!processed && (result = process_line(opt, args, arg)) != FIO_SUCCESS) return result;
      opt = argv[i]; args = 0; processed = false;
    } else {         
      if(args < max_args) arg[args++] = argv[i];
    }
  }

  if(!processed && (result = process_line(opt, args, arg)) != FIO_SUCCESS) return result;
  return FIO_SUCCESS;
}

int process_line(const std::string& opt, const int argc, const std::string argv[])
{
  if(opt=="-bootstrap") {
    if(argc==1) ibootstrap = atoi(argv[0].c_str());
  } else if(opt=="-vecpot") {
    if(argc==1) ivecpot = atoi(argv[0].c_str());
  } else if(opt=="-m3dc1") {
    return create_source(FIO_M3DC1_SOURCE, argc, argv);
  } else if(opt=="-nR") {
    if(argc==1) nR = atoi(argv[0].c_str());
  } else if(opt=="-nZ") {
    if(argc==1) nZ = atoi(argv[0].c_str());
  } else if(opt=="-nphi") {
    if(argc==1) nphi = atoi(argv[0].c_str());
  } else if(opt=="-Rmin") {
    if(argc==1) Rmin = atof(argv[0].c_str());
  } else if(opt=="-Rmax") {
    if(argc==1) Rmax = atof(argv[0].c_str());
  } else if(opt=="-Zmin") {
    if(argc==1) Zmin = atof(argv[0].c_str());
  } else if(opt=="-Zmax") {
    if(argc==1) Zmax = atof(argv[0].c_str());
  } else {
    return FIO_UNSUPPORTED;
  }
  return FIO_SUCCESS;
}

void print_usage()
{
  std::cerr << "write_neo_Eulerian\n"
      << " -nR <nR>\n"
      << " -nZ <nZ>\n"
      << " -nphi <nphi>\n"
      << " -Rmin <Rmin> -Rmax <Rmax>\n"
      << " -Zmin <Zmin> -Zmax <Zmax>\n"
      << " -bootstrap <0/1/2/3>\n"
      << " -vecpot <0/1>\n"
      << " -m3dc1 <m3dc1_source> <time> <scale> <phase>\n";
}

void delete_sources()
{
  for(auto it = sources.cbegin(); it != sources.end(); it++) {
    (*it)->close();
    delete(*it);
  }
  sources.clear();
}