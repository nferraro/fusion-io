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
          
          if(ibootstrap >= 1) {
              double v_jb = 0.0, v_dn = 0.0, v_dte = 0.0, v_dti = 0.0;
              JpdotB.eval(p, &v_jb, h);          loc_JpdotB[loc_idx] = v_jb;
              JpdotB_dndpsi.eval(p, &v_dn, h);   loc_JpdotB_dn[loc_idx] = v_dn;
              JpdotB_dtedpsi.eval(p, &v_dte, h); loc_JpdotB_dte[loc_idx] = v_dte;
              JpdotB_dtidpsi.eval(p, &v_dti, h); loc_JpdotB_dti[loc_idx] = v_dti;
          }
      }
      count_processed++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank == 0) std::cerr << "Eval done. Gathering data..." << std::endl;

  // -------------------------------------------------------------------------
  // GATHER DATA TO RANK 0
  // -------------------------------------------------------------------------
  std::vector<int> flat_counts(world_size);
  std::vector<int> flat_displs(world_size);
  long long disp = 0;
  for(int r=0; r<world_size; r++) {
      long long r_cnt = pts_per_rank + (r < rem_pts ? 1 : 0);
      flat_counts[r] = (int)r_cnt; 
      flat_displs[r] = (int)disp;
      disp += r_cnt;
  }

  double *out_Jac=nullptr, *out_Mask=nullptr;
  double *out_BR=nullptr, *out_BPhi=nullptr, *out_BZ=nullptr;
  double *out_JR=nullptr, *out_JPhi=nullptr, *out_JZ=nullptr;
  double *out_VecPot_R=nullptr, *out_VecPot_Phi=nullptr, *out_VecPot_Z=nullptr;
  double *out_Te=nullptr, *out_Ti=nullptr, *out_Ne=nullptr, *out_Ni=nullptr, *out_P=nullptr;
  double *out_JpdotB=nullptr, *out_JpdotB_dn=nullptr, *out_JpdotB_dte=nullptr, *out_JpdotB_dti=nullptr;

  if(world_rank == 0) {
      out_Jac = new double[total_work]; out_Mask = new double[total_work];
      out_BR = new double[total_work]; out_BPhi = new double[total_work]; out_BZ = new double[total_work];
      out_JR = new double[total_work]; out_JPhi = new double[total_work]; out_JZ = new double[total_work];
      if(ivecpot) { out_VecPot_R = new double[total_work]; out_VecPot_Phi = new double[total_work]; out_VecPot_Z = new double[total_work]; }
      out_Te = new double[total_work]; out_Ti = new double[total_work]; out_Ne = new double[total_work]; out_Ni = new double[total_work]; out_P = new double[total_work];
      if(ibootstrap >= 1) { out_JpdotB = new double[total_work]; out_JpdotB_dn = new double[total_work]; out_JpdotB_dte = new double[total_work]; out_JpdotB_dti = new double[total_work]; }
  }

  int my_pts = (int)my_pts_count;

  #define OPT_GATHER(loc_arr, global_arr) \
      MPI_Gatherv(loc_arr, my_pts, MPI_DOUBLE, global_arr, flat_counts.data(), flat_displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  OPT_GATHER(loc_Jac, out_Jac);
  OPT_GATHER(loc_Mask, out_Mask);
  OPT_GATHER(loc_BR, out_BR); OPT_GATHER(loc_BPhi, out_BPhi); OPT_GATHER(loc_BZ, out_BZ);
  OPT_GATHER(loc_JR, out_JR); OPT_GATHER(loc_JPhi, out_JPhi); OPT_GATHER(loc_JZ, out_JZ);
  
  if(ivecpot) { OPT_GATHER(loc_VecPot_R, out_VecPot_R); OPT_GATHER(loc_VecPot_Phi, out_VecPot_Phi); OPT_GATHER(loc_VecPot_Z, out_VecPot_Z); }

  OPT_GATHER(loc_Te, out_Te); OPT_GATHER(loc_Ti, out_Ti); OPT_GATHER(loc_Ne, out_Ne); OPT_GATHER(loc_Ni, out_Ni); OPT_GATHER(loc_P, out_P);
  
  if(ibootstrap >= 1) {
      OPT_GATHER(loc_JpdotB, out_JpdotB); OPT_GATHER(loc_JpdotB_dn, out_JpdotB_dn);
      OPT_GATHER(loc_JpdotB_dte, out_JpdotB_dte); OPT_GATHER(loc_JpdotB_dti, out_JpdotB_dti);
  }

  // -------------------------------------------------------------------------
  // VTK AND NETCDF FILE OUTPUT (Rank 0 Only) 
  // -------------------------------------------------------------------------
  if(world_rank == 0) {
      // 1. VTK OUTPUT (For ParaView, explicitly closes torus)
      std::cerr << "Writing VTK file for ParaView..." << std::endl;
      std::ofstream vtk("neo_eulerian.vtk");
      
      long long vtk_points = (long long)nR * (nphi + 1) * nZ;
      vtk << "# vtk DataFile Version 3.0\nEulerian 3D fields\nASCII\nDATASET STRUCTURED_GRID\n";
      vtk << "DIMENSIONS " << nR << " " << nphi + 1 << " " << nZ << "\n";
      vtk << "POINTS " << vtk_points << " float\n";
      
      for(int k = 0; k < nZ; k++) {
          for(int j = 0; j <= nphi; j++) {
              for(int i = 0; i < nR; i++) {
                  double r = Rmin + i * dR;
                  double p = j * dPhi;
                  double z = Zmin + k * dZ;
                  vtk << r * cos(p) << " " << r * sin(p) << " " << z << "\n";
              }
          }
      }

      vtk << "POINT_DATA " << vtk_points << "\n";
      auto write_vtk_scalar = [&](const char* name, double* data) {
          vtk << "SCALARS " << name << " float 1\nLOOKUP_TABLE default\n";
          for(int k = 0; k < nZ; k++) {
              for(int j = 0; j <= nphi; j++) {
                  int read_j = (j == nphi) ? 0 : j; 
                  for(int i = 0; i < nR; i++) {
                      long long idx = i + read_j * nR + k * nR * nphi;
                      vtk << data[idx] << "\n";
                  }
              }
          }
      };

      write_vtk_scalar("Comp_Boundary", out_Mask);
      write_vtk_scalar("Jac", out_Jac);
      write_vtk_scalar("B_R", out_BR); write_vtk_scalar("B_Phi", out_BPhi); write_vtk_scalar("B_Z", out_BZ);
      write_vtk_scalar("J_R", out_JR); write_vtk_scalar("J_Phi", out_JPhi); write_vtk_scalar("J_Z", out_JZ);
      if(ivecpot == 1) {
          write_vtk_scalar("VecPot_R", out_VecPot_R); write_vtk_scalar("VecPot_Phi", out_VecPot_Phi); write_vtk_scalar("VecPot_Z", out_VecPot_Z);
      }
      write_vtk_scalar("te_3d", out_Te); write_vtk_scalar("ti_3d", out_Ti); 
      write_vtk_scalar("ne_3d", out_Ne); write_vtk_scalar("ni_3d", out_Ni); write_vtk_scalar("p_3d", out_P);
      if(ibootstrap >= 1) {
          write_vtk_scalar("JpdotB_3d", out_JpdotB); write_vtk_scalar("JpdotB_dndpsi_3d", out_JpdotB_dn);
          write_vtk_scalar("JpdotB_dtedpsi_3d", out_JpdotB_dte); write_vtk_scalar("JpdotB_dtidpsi_3d", out_JpdotB_dti);
      }
      vtk.close();

      // 2. NETCDF OUTPUT (For Python/Xarray, strict array mapping)
      std::cerr << "Writing NetCDF file for Python..." << std::endl;
      int ncid;
      NC_CHECK(nc_create("neo_eulerian.nc", NC_CLOBBER | NC_NETCDF4, &ncid));

      int nr_dimid, nphi_dimid, nz_dimid;
      NC_CHECK(nc_def_dim(ncid, "Z", nZ, &nz_dimid));
      NC_CHECK(nc_def_dim(ncid, "Phi", nphi, &nphi_dimid)); 
      NC_CHECK(nc_def_dim(ncid, "R", nR, &nr_dimid));
      
      int r_coord_id, phi_coord_id, z_coord_id;
      NC_CHECK(nc_def_var(ncid, "R", NC_FLOAT, 1, &nr_dimid, &r_coord_id));
      NC_CHECK(nc_def_var(ncid, "Phi", NC_FLOAT, 1, &nphi_dimid, &phi_coord_id));
      NC_CHECK(nc_def_var(ncid, "Z", NC_FLOAT, 1, &nz_dimid, &z_coord_id));

      int dims[3] = {nz_dimid, nphi_dimid, nr_dimid};

      int jac_id, mask_id, br_id, bphi_id, bz_id, jr_id, jphi_id, jz_id, vecpot_r_id, vecpot_phi_id, vecpot_z_id;
      NC_CHECK(nc_def_var(ncid, "Jac", NC_FLOAT, 3, dims, &jac_id));
      NC_CHECK(nc_def_var(ncid, "Comp_Boundary", NC_FLOAT, 3, dims, &mask_id));
      NC_CHECK(nc_def_var(ncid, "B_R", NC_FLOAT, 3, dims, &br_id));
      NC_CHECK(nc_def_var(ncid, "B_Phi", NC_FLOAT, 3, dims, &bphi_id));
      NC_CHECK(nc_def_var(ncid, "B_Z", NC_FLOAT, 3, dims, &bz_id));
      NC_CHECK(nc_def_var(ncid, "J_R", NC_FLOAT, 3, dims, &jr_id));
      NC_CHECK(nc_def_var(ncid, "J_Phi", NC_FLOAT, 3, dims, &jphi_id));
      NC_CHECK(nc_def_var(ncid, "J_Z", NC_FLOAT, 3, dims, &jz_id));
      if(ivecpot == 1) {  
          NC_CHECK(nc_def_var(ncid, "VecPot_R", NC_FLOAT, 3, dims, &vecpot_r_id));
          NC_CHECK(nc_def_var(ncid, "VecPot_Phi", NC_FLOAT, 3, dims, &vecpot_phi_id));
          NC_CHECK(nc_def_var(ncid, "VecPot_Z", NC_FLOAT, 3, dims, &vecpot_z_id));
      }
      
      int te3d_id, ti3d_id, ne3d_id, ni3d_id, p3d_id;
      int JpdotB3d_id, JpdotB_dn_id, JpdotB_dte_id, JpdotB_dti_id;
      NC_CHECK(nc_def_var(ncid, "te_3d", NC_FLOAT, 3, dims, &te3d_id));
      NC_CHECK(nc_def_var(ncid, "ti_3d", NC_FLOAT, 3, dims, &ti3d_id));
      NC_CHECK(nc_def_var(ncid, "ne_3d", NC_FLOAT, 3, dims, &ne3d_id));
      NC_CHECK(nc_def_var(ncid, "ni_3d", NC_FLOAT, 3, dims, &ni3d_id));
      NC_CHECK(nc_def_var(ncid, "p_3d", NC_FLOAT, 3, dims, &p3d_id));
      
      if (ibootstrap >= 1){
          NC_CHECK(nc_def_var(ncid, "JpdotB_3d", NC_FLOAT, 3, dims, &JpdotB3d_id));
          NC_CHECK(nc_def_var(ncid, "JpdotB_dndpsi_3d", NC_FLOAT, 3, dims, &JpdotB_dn_id));
          NC_CHECK(nc_def_var(ncid, "JpdotB_dtedpsi_3d", NC_FLOAT, 3, dims, &JpdotB_dte_id));
          NC_CHECK(nc_def_var(ncid, "JpdotB_dtidpsi_3d", NC_FLOAT, 3, dims, &JpdotB_dti_id));
      }
      
      NC_CHECK(nc_enddef(ncid));

      double* r_ax = new double[nR];   for(int i=0; i<nR; i++) r_ax[i] = Rmin + i * dR;
      double* phi_ax = new double[nphi]; for(int j=0; j<nphi; j++) phi_ax[j] = j * dPhi; 
      double* z_ax = new double[nZ];   for(int k=0; k<nZ; k++) z_ax[k] = Zmin + k * dZ;
      NC_CHECK(nc_put_var_double(ncid, r_coord_id, r_ax));
      NC_CHECK(nc_put_var_double(ncid, phi_coord_id, phi_ax));
      NC_CHECK(nc_put_var_double(ncid, z_coord_id, z_ax));

      NC_CHECK(nc_put_var_double(ncid, jac_id, out_Jac));
      NC_CHECK(nc_put_var_double(ncid, mask_id, out_Mask));
      NC_CHECK(nc_put_var_double(ncid, br_id, out_BR));
      NC_CHECK(nc_put_var_double(ncid, bphi_id, out_BPhi));
      NC_CHECK(nc_put_var_double(ncid, bz_id, out_BZ));
      NC_CHECK(nc_put_var_double(ncid, jr_id, out_JR));
      NC_CHECK(nc_put_var_double(ncid, jphi_id, out_JPhi));
      NC_CHECK(nc_put_var_double(ncid, jz_id, out_JZ));
      if(ivecpot == 1) {
          NC_CHECK(nc_put_var_double(ncid, vecpot_r_id, out_VecPot_R));
          NC_CHECK(nc_put_var_double(ncid, vecpot_phi_id, out_VecPot_Phi));
          NC_CHECK(nc_put_var_double(ncid, vecpot_z_id, out_VecPot_Z));
      }
      NC_CHECK(nc_put_var_double(ncid, te3d_id, out_Te));
      NC_CHECK(nc_put_var_double(ncid, ti3d_id, out_Ti));
      NC_CHECK(nc_put_var_double(ncid, ne3d_id, out_Ne));
      NC_CHECK(nc_put_var_double(ncid, ni3d_id, out_Ni));
      NC_CHECK(nc_put_var_double(ncid, p3d_id, out_P));
      if(ibootstrap >= 1) {
          NC_CHECK(nc_put_var_double(ncid, JpdotB3d_id, out_JpdotB));
          NC_CHECK(nc_put_var_double(ncid, JpdotB_dn_id, out_JpdotB_dn));
          NC_CHECK(nc_put_var_double(ncid, JpdotB_dte_id, out_JpdotB_dte));
          NC_CHECK(nc_put_var_double(ncid, JpdotB_dti_id, out_JpdotB_dti));
      }
      NC_CHECK(nc_close(ncid));
      
      delete[] r_ax; delete[] phi_ax; delete[] z_ax;
      std::cerr << "Successfully wrote output files!" << std::endl;
  }

  // -------------------------------------------------------------------------
  // CLEANUP MEMORY ON ALL RANKS
  // -------------------------------------------------------------------------
  delete[] loc_Jac; delete[] loc_Mask; delete[] loc_BR; delete[] loc_BPhi; delete[] loc_BZ;
  delete[] loc_JR; delete[] loc_JPhi; delete[] loc_JZ;
  delete[] loc_VecPot_R; delete[] loc_VecPot_Phi; delete[] loc_VecPot_Z;
  delete[] loc_Te; delete[] loc_Ti; delete[] loc_Ne; delete[] loc_Ni; delete[] loc_P; 
  delete[] loc_JpdotB; delete[] loc_JpdotB_dn; delete[] loc_JpdotB_dte; delete[] loc_JpdotB_dti;

  if(world_rank == 0) {
      delete[] out_Jac; delete[] out_Mask; delete[] out_BR; delete[] out_BPhi; delete[] out_BZ;
      delete[] out_JR; delete[] out_JPhi; delete[] out_JZ;
      if(ivecpot) { delete[] out_VecPot_R; delete[] out_VecPot_Phi; delete[] out_VecPot_Z; }
      delete[] out_Te; delete[] out_Ti; delete[] out_Ne; delete[] out_Ni; delete[] out_P; 
      if(ibootstrap >= 1) { delete[] out_JpdotB; delete[] out_JpdotB_dn; delete[] out_JpdotB_dte; delete[] out_JpdotB_dti; }
  }

  sources[0]->deallocate_search_hint(&h);
  delete_sources();

  if(world_rank == 0) std::cerr << "Done." << std::endl;
  MPI_Finalize();
  return 0;
}

// ========================================================
// Source Handlers
// ========================================================

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

  fio_field* field;

  // Unconditionally load bootstrap fields so arg-parsing order doesn't break initialization.
  if(src->get_field(FIO_JBS, &field, &fopt) == FIO_SUCCESS)
      JpdotB.add_field(field, FIO_ADD, 1., hint);
      
  if(src->get_field(FIO_JBS_dndpsi, &field, &fopt) == FIO_SUCCESS)
      JpdotB_dndpsi.add_field(field, FIO_ADD, 1., hint);
      
  if(src->get_field(FIO_JBS_dtedpsi, &field, &fopt) == FIO_SUCCESS)
      JpdotB_dtedpsi.add_field(field, FIO_ADD, 1., hint);
      
  if(src->get_field(FIO_JBS_dtidpsi, &field, &fopt) == FIO_SUCCESS)
      JpdotB_dtidpsi.add_field(field, FIO_ADD, 1., hint);

  // Load everything else
  if(src->get_field(FIO_MAGNETIC_FIELD, &field, &fopt) == FIO_SUCCESS) mag.add_field(field, FIO_ADD, 1., hint);
  if(src->get_field(FIO_CURRENT_DENSITY, &field, &fopt) == FIO_SUCCESS) current_density.add_field(field, FIO_ADD, 1., hint);
  if(src->get_field(FIO_VECTOR_POTENTIAL, &field, &fopt) == FIO_SUCCESS) vector_potential.add_field(field, FIO_ADD, 1., hint);

  fopt.set_option(FIO_SPECIES, FIO_ELECTRON);
  if(src->get_field(FIO_TEMPERATURE, &field, &fopt) == FIO_SUCCESS) electron_temperature.add_field(field, FIO_ADD, 1., hint);
  if(src->get_field(FIO_DENSITY, &field, &fopt) == FIO_SUCCESS) electron_density.add_field(field, FIO_ADD, 1., hint);

  fopt.set_option(FIO_SPECIES, FIO_MAIN_ION);
  if(src->get_field(FIO_DENSITY, &field, &fopt) == FIO_SUCCESS) ion_density.add_field(field, FIO_ADD, 1., hint);
  if(src->get_field(FIO_TEMPERATURE, &field, &fopt) == FIO_SUCCESS) ion_temperature.add_field(field, FIO_ADD, 1., hint);  
  if(src->get_field(FIO_TOTAL_PRESSURE, &field, &fopt) == FIO_SUCCESS) pressure.add_field(field, FIO_ADD, 1., hint);

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