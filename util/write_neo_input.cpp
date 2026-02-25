#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <compound_field.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <deque>
#include <algorithm>
#include <netcdf.h>

// ========================================================
// Global Variables & Configuration
// ========================================================
int nr = 10;      // number of radial points for surfaces
int ntheta = 400; // number of poloidal points
int nphi = 32;    // number of toroidal points
int ibootstrap;
int ivecpot;
int hybrid = 0;       // 0: use regular gridded isosurface, 1: use hybrid method (for surfaces closer to LCFS))
double dl_tor = 0.01;     // step size when finding surfaces
double dl_pol = 0.01;     // step size when finding surfaces
double max_step = 0.02;   // Maximum step size for Newton iterations
double tol = 0.1;         // Tolerance for Te when finding isosurface
double dR0 = 0.0;         // Guess for offset from magnetic axis
double R0 = 0.0;          // Guess for R coordinate of magnetic axis
double psi_start = -1;
double psi_end = -1;
double te_start = -1;
double te_end = -1;

std::deque<fio_source*> sources;
fio_compound_field electron_density, electron_temperature, pressure;
fio_compound_field current_density, vector_potential;
fio_compound_field ion_density, ion_temperature, psin, mag, JpdotB, JpdotB_dndpsi, JpdotB_dtedpsi, JpdotB_dtidpsi;
fio_compound_field JpdotB_L31, JpdotB_L32, JpdotB_L34, JpdotB_alpha;
double axis[3];
double psi0, psi1;

// ========================================================
// Flags & Parameters for Calculating Bootstrap Coeffs
// ========================================================
int UseSauter0Redl1 = 1;   // 0: Sauter, 1: Redl
int nfp = 1;               // Stellarator field periods
int Nzp = 0;               // Toroidal mode number
int UseBfore = 1;
int lambda_count = 100;
double tol_trap = 1e-3;
double Zcharge_e_input = 1.0;
double Zcharge_eff = 1.0;  // For Redl terms

// MPI Globals
int world_rank = 0;
int world_size = 1;

void print_usage();
int process_command_line(int argc, char* argv[]);
int process_line(const std::string& opt, const int argc, 
		  const std::string argv[]);
void delete_sources();

// Data Structures

struct CurrentIntegrals {
    std::vector<double> Izp_sumj;
    std::vector<double> Izt_sumk;
    std::vector<double> GplusiI_fa;
    std::vector<double> Gbar_by_iminusN;
};
struct TrapResult {
    std::vector<double> ftrap;
    double Zcharge_e;
    double Zcharge_i;
    double Zeff;
};



struct CollisionResults {
    std::vector<double> nu_e_star, nu_i_star, epsilon, qR;
};

struct RedlResults {
    std::vector<double> L31, L32, L34, alpha_nu_i;
    // std::vector<double> ft31, ft32ee, ft32ei, ft34; // Optional: for debugging
};

struct SauterResults {
    std::vector<double> L31, L32, L34, alpha_nu_i;
   // std::vector<double> ft31, ft32ee, ft32ei, ft34; // Optional: for debugging
};

CurrentIntegrals calculate_current_integrals(
    int nr, int nphi, int ntheta, int nfp, int Nzp,
    const double* R, const double* Z, 
    const double* BR, const double* BPhi, const double* BZ,
    const double* q);

std::vector<double> calc_ftrap(
    int nr, int nphi, int ntheta, int lambda_count, 
    const double* Bmax_fpsi, const double* B2_fa, const double* Bmag, 
    const double* Jac, const double* psi_eval, double tol);

void calculate_Zcharge(const double* ne_fa, const double* ni_fa, int nr, double Zinput, 
                       double& Z_e, double& Z_i, double& Zeff_out);

CollisionResults calculate_collision_frequencies(
    int nr, int nphi, int ntheta,
    const double* ne_fa, const double* ni_fa, 
    const double* Te_fa, const double* Ti_fa,
    double Ze, double Zi, const double* Bmax, const double* Bmin,
    const double* q, const double* GplusiI, const double* Binv_fa,
    const double* R,
    int UseBfore, int Nzp);


RedlResults calculate_redl_terms(int nr, double Zeff, 
                                 const std::vector<double>& ftrap, 
                                 const std::vector<double>& nue, 
                                 const std::vector<double>& nui);

SauterResults calculate_sauter_terms(int nr, double Zeff, const std::vector<double>& ftrap, 
                                 const std::vector<double>& nue, const std::vector<double>& nui);


struct NeoCoeffs {
    std::vector<double> L31, L32, L34, alpha_nu_i;
};


// ========================================================
// MAIN
// ========================================================
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int result;

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

  if(sources.size()==0) {
    if(world_rank == 0) {
      std::cerr << "Error, no source is loaded" << std::endl;
      print_usage();
    }
    MPI_Finalize();
    return 1;
  }

  //Defaults
  if(psi_start <= 0. || psi_start > 1.) psi_start = 0.;
  if(psi_end <= 0. || psi_end > 1.) psi_end = 1.;
  if(psi_start <= 0.) psi_start = (psi_end-psi_start)/nr;
  if(psi_end==1.) psi_end = 1. - (psi_end-psi_start)/nr;

   // Rank 0 Setup and Prints
  if(world_rank == 0) {

    std::cerr << "Input parameters\n=======================\n";
    std::cerr << "Maximum Newton iteration step (max_step) = " << max_step
        << '\n';
    std::cerr << "Radial grid points (nr) = " << nr << '\n';
    std::cerr << "Toroidal grid points (nphi) = " << nphi << '\n';
    std::cerr << "Poloidal grid points (ntheta) = " << ntheta << '\n';
    std::cerr << "First value of psi_norm (psi_start) = "
        << psi_start << '\n';
    std::cerr << "Last value of psi_norm (psi_end) = "
        << psi_end << '\n';
    std::cerr << "Tolerance for finding isosurface (tol) = " << tol << " eV\n";
    std::cerr << "=======================" << std::endl;
  }
  


  // Calculate 3D surfaces
  double** path = new double*[3];
  path[0] = new double[nr*nphi*ntheta];
  path[1] = new double[nr*nphi*ntheta];
  path[2] = new double[nr*nphi*ntheta];

  double* theta = new double[ntheta];
  double* phi = new double[nphi];
  double** path_plane = new double*[3];  // contains path of one toroidal plane
  double** path_surf = new double*[3];   // contains path of one surface
  //  double* bpol = new double[nphi*ntheta];
  //  double* btor = new double[nr*nphi*ntheta];
  double** b = new double*[3];
  double** b_surf = new double*[3];
  b[0] = new double[nr*nphi*ntheta];
  b[1] = new double[nr*nphi*ntheta];
  b[2] = new double[nr*nphi*ntheta];

  double** currden = new double*[3];
  currden[0] = new double[nr*nphi*ntheta];
  currden[1] = new double[nr*nphi*ntheta];
  currden[2] = new double[nr*nphi*ntheta];

  double** vecpot = new double*[3];
  vecpot[0] = new double[nr*nphi*ntheta];
  vecpot[1] = new double[nr*nphi*ntheta];
  vecpot[2] = new double[nr*nphi*ntheta];
  
  

  double* telec = new double[nr*nphi*ntheta];
  double* tion = new double[nr*nphi*ntheta];
  double* nelec= new double[nr*nphi*ntheta];
  double* nion = new double[nr*nphi*ntheta];
  double* press = new double[nr*nphi*ntheta];
  double* jdotb = new double[nr*nphi*ntheta];
  double* jdotB_dndpsi = new double[nr*nphi*ntheta];
  double* jdotB_dtedpsi = new double[nr*nphi*ntheta];
  double* jdotB_dtidpsi = new double[nr*nphi*ntheta];
  double* jdotB_L31 = new double[nr*nphi*ntheta];
  double* jdotB_L32 = new double[nr*nphi*ntheta];
  double* jdotB_L34 = new double[nr*nphi*ntheta];
  double* jdotB_alpha = new double[nr*nphi*ntheta];


  
  // 1D Profile Arrays
  double* q = new double[nr];
  double* psi_surf = new double[nr];
  double* ne = new double[nr];
  double* te = new double[nr];
  double* ni = new double[nr];
  double* ti = new double[nr];
  double* JpdotB_fluxavg = new double[nr];
  double* JpdotB_dndpsi_fluxavg = new double[nr];
  double* JpdotB_dtedpsi_fluxavg = new double[nr];
  double* JpdotB_dtidpsi_fluxavg = new double[nr];
  double* JpdotB_L31_fluxavg = new double[nr];
  double* JpdotB_L32_fluxavg = new double[nr];
  double* JpdotB_L34_fluxavg = new double[nr];
  double* JpdotB_alpha_fluxavg = new double[nr];
  
  double* psi_t = new double[nr];
  double* dpsi_t = new double[nr];
  double* psi_p = new double[nr];
  double* dpsi_p = new double[nr];

  double* bx_fa = new double[nr];
  double* by_fa = new double[nr];
  double* bz_fa = new double[nr];
  double* bmag_fa = new double[nr];
  double* b2_fa = new double[nr];
  double* jx_fa = new double[nr];
  double* jy_fa = new double[nr];
  double* jz_fa = new double[nr];

  
  double* axis_3d[3];
  axis_3d[0] = new double[nphi];
  axis_3d[1] = new double[nphi];
  axis_3d[2] = new double[nphi];

  double psi_norm[nr];
  double te_max;

  // Find magnetic axis
  axis[1] = 0.;
  fio_source* src0 = sources[0];
  void* h;
  src0->allocate_search_hint(&h);
  double x[3];

  // -------------------------------------------------------------------------
  // FIND MAGNETIC AXIS (Rank 0 only)
  // ---------------------
    if(world_rank == 0) {
      axis[1] = 0.;
      double n[3];
      
      n[0] = 0.;
      n[1] = 1.;
      n[2] = 0.;

      if(R0 > 0.0) {
        x[0] = R0;
      } else {		 
        x[0] = axis[0] + dR0 - 0.01;
      }
      x[1] = 0.;
      x[2] = 0.;
      result = fio_find_max(&electron_temperature, &te_max, x, tol, max_step, 2, n, h);
      if(result != FIO_SUCCESS) {
        std::cerr << "Error finding magnetic axis " << std::endl;
        return result;
      } else {
        std::cerr << "Found magnetic axis" << std::endl;
        std::cerr << "TE MAX = " << te_max << " AT ("
            << x[0] << ", " << x[1] << ", " << x[2] << ") " 
            << std::endl;

        axis[0] = x[0];
        axis[1] = 0;
        axis[2] = x[2];
      }
    }
  

  // Broadcast Axis and TeMax
  MPI_Bcast(axis, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&te_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Rank 0 handles the gplot files creation (the index file)
  std::ofstream gplot, splot;
  if(world_rank == 0) {
      gplot.open("plot_surfaces", std::ofstream::out|std::ofstream::trunc);
      splot.open("splot_surfaces", std::ofstream::out|std::ofstream::trunc);
      splot << "set hidden3d\nsplot ";
      std::cerr << "Calculating surfaces..." << std::endl;
  }

    // Decompose work
  int surfaces_per_rank = nr / world_size;
  int remainder = nr % world_size;
  int my_start = world_rank * surfaces_per_rank + std::min(world_rank, remainder);
  int my_count = surfaces_per_rank + (world_rank < remainder ? 1 : 0);
  int my_end = my_start + my_count;

  // Temp buffers for fio_gridded_isosurface pointers
  double psi_norm_val;

  for(int s=my_start; s<my_end; s++) {
    double te0 = 0.0;
    double target_val_for_solver = 0.0;
    bool seed_found = false;

    x[0] = axis[0] - 0.01;
    x[1] = 0.;
    x[2] = axis[2];

    if((te_end >= 0.) && (te_start >= 0.)) {
      // use Te as radial coordinate
      te0 = (te_end - te_start)*s/(nr - 1.) + te_start;
      target_val_for_solver = te0;
      psi_norm[s] = (te_max - te0) / te_max;

      // Optimization: Only try "Fast Search" if we are NOT at the edge.
      // T_e=0 is numerically unstable for Newton, so skip straight to Robust.
      if (te0 > 50.0) { 
             result = fio_find_val(&electron_temperature, te0, x, tol, max_step, 1, axis, h);
             if (result == FIO_SUCCESS) seed_found = true;
      } else {
        //use previous surface as x for quicker convergenvce
        double te0_prev = (te_end - te_start)*(s-1)/(nr - 1.) + te_start;        
        result = fio_find_val(&electron_temperature, te0_prev, x, tol, max_step, 1, axis, h);
        seed_found = false;
        target_val_for_solver = te0;
        psi_norm[s] = (te_max - te0) / te_max;
      }
    
    }
      else {
      // use psi_norm as radial coordinate
      if(nr==1) {
	      psi_norm[s] = psi_start;
      } else {
	      psi_norm[s] = (psi_end-psi_start)*s/(nr-1.) + psi_start;
      }
      result = fio_find_val(&psin, psi_norm[s], x, 1e-4, 0.1, 1, axis, h);

      if(result == FIO_SUCCESS) {
          double temp_on_surf;
          electron_temperature.eval(x, &temp_on_surf, h);
          target_val_for_solver = temp_on_surf;
          seed_found = true;
          
      } else {
          target_val_for_solver = te_max * (1.0 - psi_norm[s]);
          seed_found = false;
          x[0] = axis[0] + 0.05;
          x[1] = 0.0;
          x[2] = axis[2];
        }
    }


    std::cerr << "[Rank " << world_rank << "] Solving Surface " << s+1 << " of " << nr 
              << " | Target Te=" << target_val_for_solver << "|Psi_norm="<<psi_norm[s]
              << " | Method: " << (seed_found ? "Fast Seed" : "Radial with x from prev surface") << std::endl;

    char* label;
    char surf_label[256];
    sprintf(surf_label, "%d", s);
    label = surf_label;

    path_surf[0] = &(path[0][s*nphi*ntheta]);
    path_surf[1] = &(path[1][s*nphi*ntheta]);
    path_surf[2] = &(path[2][s*nphi*ntheta]);

    //Setup Axis for solver
    axis_3d[0][0] = axis[0];
    axis_3d[1][0] = 0;
    axis_3d[2][0] = axis[2];
    if(hybrid == 0) {
      result = fio_gridded_isosurface(&electron_temperature, target_val_for_solver, x,
                    axis_3d, dl_tor, dl_pol, tol, max_step,
                    nphi, ntheta, 
                    phi, theta,
                    path_surf, label, h);
      if(world_rank == 0) {
            std::cerr << "Using regular gridded isosurface" << std::endl;
      }
    } else if (hybrid == 1) {
      result = fio_gridded_isosurface_hybrid(&electron_temperature, target_val_for_solver, x,
                    axis_3d, dl_tor, dl_pol, tol, max_step,
                    nphi, ntheta, 
                    phi, theta,
                    path_surf, label, h);
      if(world_rank == 0) {
            std::cerr << "Using hybrid version of gridded isosurface" << std::endl;
      }
    } else {
      std::cerr << "Error: Invalid value for -hybrid option. Use 0 for regular gridded isosurface, 1 for hybrid method." << std::endl;
      MPI_Finalize();
      return 1;
    }
    te[s] = target_val_for_solver;
    psi_surf[s] = psi_norm[s]*(psi1 - psi0) + psi0;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Reconstruct gplot commands on Rank 0
  if(world_rank == 0) {
      for(int s = 0; s < nr; s++) {
          // Add 2D plot command (R vs Z)
          gplot << "' surface_" << s << ".dat' u 2:3 w l";
          
          // Add 3D plot command (Phi, R, Z)
          splot << "'surface_" << s << ".dat' u 1:2:3 w l";
          
          if(s < nr - 1) {
              gplot << ", \\\n ";
              splot << ", \\\n ";
          }
      }
      gplot << std::endl;
      splot << std::endl;
      
      gplot.close();
      splot.close();
      std::cerr << "Gnuplot scripts 'plot_surfaces' and 'splot_surfaces' generated." << std::endl;
  }

 
  // -------------------------------------------------------------------------
  // gather GEOMETRY ACROSS ALL RANKS
  // -------------------------------------------------------------------------
  std::vector<int> count_bytes(world_size);
  std::vector<int> disp_bytes(world_size);
  std::vector<int> count_dbls(world_size); 
  std::vector<int> disp_dbls(world_size);

  int current_disp = 0;
  int current_disp_dbl = 0;
  for(int r=0; r<world_size; r++) {
      int cnt = (nr / world_size) + (r < remainder ? 1 : 0);
      count_dbls[r] = cnt;
      disp_dbls[r] = current_disp_dbl;
      
      count_bytes[r] = cnt * nphi * ntheta;
      disp_bytes[r] = current_disp;
      
      current_disp_dbl += cnt;
      current_disp += count_bytes[r];
  }

  long long my_geo_offset = (long long)my_start * nphi * ntheta;
  int my_geo_count = count_bytes[world_rank];
  std::vector<double> send_buf_geo(my_geo_count > 0 ? my_geo_count : 1);

 
  if(my_geo_count > 0) std::memcpy(send_buf_geo.data(), path[0] + my_geo_offset, my_geo_count * sizeof(double));
  MPI_Allgatherv(send_buf_geo.data(), my_geo_count, MPI_DOUBLE, 
                 path[0], count_bytes.data(), disp_bytes.data(), MPI_DOUBLE, MPI_COMM_WORLD);


  if(my_geo_count > 0) std::memcpy(send_buf_geo.data(), path[1] + my_geo_offset, my_geo_count * sizeof(double));
  MPI_Allgatherv(send_buf_geo.data(), my_geo_count, MPI_DOUBLE, 
                 path[1], count_bytes.data(), disp_bytes.data(), MPI_DOUBLE, MPI_COMM_WORLD);


  if(my_geo_count > 0) std::memcpy(send_buf_geo.data(), path[2] + my_geo_offset, my_geo_count * sizeof(double));
  MPI_Allgatherv(send_buf_geo.data(), my_geo_count, MPI_DOUBLE, 
                 path[2], count_bytes.data(), disp_bytes.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  
  
  int my_prof_count = count_dbls[world_rank];
  std::vector<double> send_buf_prof(my_prof_count > 0 ? my_prof_count : 1);

  if(my_prof_count > 0) std::memcpy(send_buf_prof.data(), te + my_start, my_prof_count * sizeof(double));
  MPI_Allgatherv(send_buf_prof.data(), my_prof_count, MPI_DOUBLE, 
                 te, count_dbls.data(), disp_dbls.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  if(my_prof_count > 0) std::memcpy(send_buf_prof.data(), psi_surf + my_start, my_prof_count * sizeof(double));
  MPI_Allgatherv(send_buf_prof.data(), my_prof_count, MPI_DOUBLE, 
                 psi_surf, count_dbls.data(), disp_dbls.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  

  // -------------------------------------------------------------------------
  //  CALCULATE JACOBIAN on rank 0
  // -------------------------------------------------------------------------
  
  double* jac = new double[nr*nphi*ntheta];
  if(world_rank == 0) {
      std::cerr << "Calculating the Jacobian..." << std::endl;
      result = fio_isosurface_jacobian(nr, nphi, ntheta, path, jac);
      if(result != FIO_SUCCESS) std::cerr << "Error evaluating jacobian" << std::endl;
  }
  // Broadcast Jacobian to all ranks
  MPI_Bcast(jac, nr*nphi*ntheta, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // -------------------------------------------------------------------------
  // PARALLEL FIELD EVALUATION 
  // -------------------------------------------------------------------------
  if(world_rank == 0) std::cerr << "Evaluating fields on surfaces..." << std::endl;

  // New Decomposition by TOTAL POINTS for load balancing
  long long total_work = (long long)nr * nphi * ntheta;
  long long pts_per_rank = total_work / world_size;
  long long rem_pts = total_work % world_size;

  long long my_pts_start = world_rank * pts_per_rank + std::min((long long)world_rank, rem_pts);
  long long my_pts_count = pts_per_rank + (world_rank < rem_pts ? 1 : 0);
  long long my_pts_end = my_pts_start + my_pts_count;

  // Local Partial Accumulators 
  std::vector<double> part_te(nr, 0.0), part_ne(nr, 0.0), part_ni(nr, 0.0), part_ti(nr, 0.0);
  std::vector<double> part_bx(nr, 0.0), part_by(nr, 0.0), part_bz(nr, 0.0);
  std::vector<double> part_bmag(nr, 0.0), part_b2(nr, 0.0);
  std::vector<double> part_jx(nr, 0.0), part_jy(nr, 0.0), part_jz(nr, 0.0);
  std::vector<double> part_jb(nr, 0.0), part_jb_dn(nr, 0.0), part_jb_dte(nr, 0.0), part_jb_dti(nr, 0.0);
  std::vector<double> part_L31(nr, 0.0), part_L32(nr, 0.0), part_L34(nr, 0.0), part_alpha(nr, 0.0);
  std::vector<double> part_dV(nr, 0.0);

  long long count_processed = 0;
  long long print_interval = my_pts_count / 10; 
  if(print_interval == 0) print_interval = 1;

  for(long long idx = my_pts_start; idx < my_pts_end; idx++) {
      
      // Progress Print
      if(count_processed % print_interval == 0) {
           std::cerr << "[Rank " << world_rank << "] Progress: " 
                     << (count_processed * 100 / my_pts_count) << "%" << std::flush << std::endl;
      }
      count_processed++;

      // Map global point index 'idx' back to surface 's'
      int s = idx / (nphi * ntheta);

      double p[3] = { path[0][idx], path[1][idx], path[2][idx] };
      double v3[3], v1;
      double J = jac[idx];

      // Evaluate B Field
      mag.eval(p, v3, h);
      b[0][idx] = v3[0]; b[1][idx] = v3[1]; b[2][idx] = v3[2];
      
      // Evaluate Current
      current_density.eval(p, v3, h);
      currden[0][idx] = v3[0]; currden[1][idx] = v3[1]; currden[2][idx] = v3[2];

      // Evaluate VecPot
      if(ivecpot) {
          vector_potential.eval(p, v3, h);
          vecpot[0][idx] = v3[0]; vecpot[1][idx] = v3[1]; vecpot[2][idx] = v3[2];
      }
      
      // Evaluate Scalars
      electron_density.eval(p, &v1, h); nelec[idx] = v1;
      ion_density.eval(p, &v1, h);      nion[idx] = v1;
      ion_temperature.eval(p, &v1, h);  tion[idx] = v1;
      electron_temperature.eval(p, &v1, h);  telec[idx] = v1;
      pressure.eval(p, &v1, h);         press[idx] = v1;
      

      //Calculate Partial Sums for Flux Averages here
      part_dV[s] += J;
      part_ne[s] += nelec[idx] * J;
      part_ni[s] += nion[idx] * J;
      part_ti[s] += tion[idx] * J;
      part_te[s] += telec[idx] * J;

      part_bx[s] += b[0][idx] * J; 
      part_by[s] += b[1][idx] * J; 
      part_bz[s] += b[2][idx] * J;
      double bmod = sqrt(b[0][idx]*b[0][idx] + b[1][idx]*b[1][idx] + b[2][idx]*b[2][idx]);
      part_bmag[s] += bmod * J;
      part_b2[s]   += (bmod*bmod) * J;
      
      part_jx[s] += currden[0][idx] * J;
      part_jy[s] += currden[1][idx] * J;
      part_jz[s] += currden[2][idx] * J;

      if(ibootstrap >= 1) {
          JpdotB.eval(p, &v1, h);         jdotb[idx] = v1; part_jb[s] += v1 * J;
          JpdotB_dndpsi.eval(p, &v1, h);  part_jb_dn[s] += v1 * J;
          JpdotB_dtedpsi.eval(p, &v1, h); part_jb_dte[s] += v1 * J;
          JpdotB_dtidpsi.eval(p, &v1, h); part_jb_dti[s] += v1 * J;
      }
      if(ibootstrap == 3) {
          JpdotB_L31.eval(p, &v1, h);   jdotB_L31[idx] = v1; part_L31[s] += v1 * J;
          JpdotB_L32.eval(p, &v1, h);   jdotB_L32[idx] = v1; part_L32[s] += v1 * J;
          JpdotB_L34.eval(p, &v1, h);   jdotB_L34[idx] = v1; part_L34[s] += v1 * J;
          JpdotB_alpha.eval(p, &v1, h); jdotB_alpha[idx] = v1; part_alpha[s] += v1 * J;
      }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank == 0) std::cerr << "Eval done. Gathering data..." << std::endl;

  
  std::vector<int> flat_counts(world_size);
  std::vector<int> flat_displs(world_size);
  long long disp = 0;
  for(int r=0; r<world_size; r++) {
      long long r_cnt = pts_per_rank + (r < rem_pts ? 1 : 0);
      flat_counts[r] = (int)r_cnt; 
      flat_displs[r] = (int)disp;
      disp += r_cnt;
  }

  // --- SAFE GATHER: Copy to temp buffer first ---
  int my_pts = (int)my_pts_count;
  std::vector<double> send_buf_3d(my_pts > 0 ? my_pts : 1);

 #define SAFE_GATHER_3D(arr) \
      if(my_pts > 0) std::memcpy(send_buf_3d.data(), arr + my_pts_start, my_pts * sizeof(double)); \
      MPI_Gatherv(send_buf_3d.data(), my_pts, MPI_DOUBLE, \
                  arr, flat_counts.data(), flat_displs.data(), MPI_DOUBLE, \
                  0, MPI_COMM_WORLD);

  SAFE_GATHER_3D(b[0]); 
  SAFE_GATHER_3D(b[1]); 
  SAFE_GATHER_3D(b[2]);
  
  SAFE_GATHER_3D(currden[0]); 
  SAFE_GATHER_3D(currden[1]); 
  SAFE_GATHER_3D(currden[2]);
  
  if(ivecpot) { 
      SAFE_GATHER_3D(vecpot[0]); 
      SAFE_GATHER_3D(vecpot[1]); 
      SAFE_GATHER_3D(vecpot[2]); 
  }

  SAFE_GATHER_3D(telec); 
  SAFE_GATHER_3D(nelec); 
  SAFE_GATHER_3D(tion); 
  SAFE_GATHER_3D(nion); 
  SAFE_GATHER_3D(press);
  
  if(ibootstrap >= 1) SAFE_GATHER_3D(jdotb);


  // Reduce 1D Flux Averages (Sum partials from all ranks)
  std::vector<double> total_dV(nr);
  MPI_Reduce(part_dV.data(), total_dV.data(), nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(part_ne.data(), ne, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_ni.data(), ni, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_ti.data(), ti, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(part_bx.data(), bx_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_by.data(), by_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_bz.data(), bz_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_bmag.data(), bmag_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_b2.data(), b2_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(part_jx.data(), jx_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_jy.data(), jy_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_jz.data(), jz_fa, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ibootstrap >= 1) {
      MPI_Reduce(part_jb.data(), JpdotB_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(part_jb_dn.data(), JpdotB_dndpsi_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(part_jb_dte.data(), JpdotB_dtedpsi_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(part_jb_dti.data(), JpdotB_dtidpsi_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  if(ibootstrap == 3) {
      MPI_Reduce(part_L31.data(), JpdotB_L31_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(part_L32.data(), JpdotB_L32_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(part_L34.data(), JpdotB_L34_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(part_alpha.data(), JpdotB_alpha_fluxavg, nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  // Normalize Averages
  if(world_rank == 0) {
      for(int i=0; i<nr; i++) {
          if(total_dV[i] > 0) {
              ne[i] /= total_dV[i]; ni[i] /= total_dV[i]; ti[i] /= total_dV[i];
              bx_fa[i] /= total_dV[i]; by_fa[i] /= total_dV[i]; bz_fa[i] /= total_dV[i];
              bmag_fa[i] /= total_dV[i]; b2_fa[i] /= total_dV[i];
              jx_fa[i] /= total_dV[i]; jy_fa[i] /= total_dV[i]; jz_fa[i] /= total_dV[i];
              
              if(ibootstrap >= 1) {
                  JpdotB_fluxavg[i] /= total_dV[i];
                  JpdotB_dndpsi_fluxavg[i] /= total_dV[i];
                  JpdotB_dtedpsi_fluxavg[i] /= total_dV[i];
                  JpdotB_dtidpsi_fluxavg[i] /= total_dV[i];
              }
              if(ibootstrap == 3) {
                  JpdotB_L31_fluxavg[i] /= total_dV[i];
                  JpdotB_L32_fluxavg[i] /= total_dV[i];
                  JpdotB_L34_fluxavg[i] /= total_dV[i];
                  JpdotB_alpha_fluxavg[i] /= total_dV[i];
              }
          }
      }
  }

  // -------------------------------------------------------------------------
  // POST-PROCESSING for neo write (Rank 0 Only) 
  // -------------------------------------------------------------------------
  if(world_rank == 0) {
    // swap indices
    double* R = new double[nr*nphi*ntheta];
    double* Z = new double[nr*nphi*ntheta];
    double* BR   = new double[nr*nphi*ntheta];
    double* BPhi = new double[nr*nphi*ntheta];
    double* BZ   = new double[nr*nphi*ntheta];
    double* Jac = new double[nr*nphi*ntheta];

    double* Te_3d = new double[nr*nphi*ntheta];
    double* Ti_3d = new double[nr*nphi*ntheta];
    double* Ne_3d = new double[nr*nphi*ntheta];
    double* Ni_3d = new double[nr*nphi*ntheta];
    double* P_3d = new double[nr*nphi*ntheta];
    double* JpdotB_3d = new double[nr*nphi*ntheta];
    double* VecPot_R = new double[nr*nphi*ntheta];
    double* VecPot_Phi = new double[nr*nphi*ntheta];
    double* VecPot_Z = new double[nr*nphi*ntheta];
  
    for(int i=0; i<nr; i++) {
      for(int j=0; j<nphi; j++) {
        for(int k=0; k<ntheta; k++) {
          Jac[i + j*nr + k*nr*nphi] = jac[k + j*ntheta + i*ntheta*nphi];
          R[i + j*nr + k*nr*nphi] = path[0][k + j*ntheta + i*ntheta*nphi];
          Z[i + j*nr + k*nr*nphi] = path[2][k + j*ntheta + i*ntheta*nphi];
          BR  [i + j*nr + k*nr*nphi] = b[0][k + j*ntheta + i*ntheta*nphi];
          BPhi[i + j*nr + k*nr*nphi] = b[1][k + j*ntheta + i*ntheta*nphi];
          BZ  [i + j*nr + k*nr*nphi] = b[2][k + j*ntheta + i*ntheta*nphi];
          Te_3d  [i + j*nr + k*nr*nphi] = telec[k + j*ntheta + i*ntheta*nphi];
          Ti_3d  [i + j*nr + k*nr*nphi] = tion[k + j*ntheta + i*ntheta*nphi];
          Ne_3d  [i + j*nr + k*nr*nphi] = nelec[k + j*ntheta + i*ntheta*nphi];
          Ni_3d  [i + j*nr + k*nr*nphi] = nion[k + j*ntheta + i*ntheta*nphi];
          P_3d  [i + j*nr + k*nr*nphi] = press[k + j*ntheta + i*ntheta*nphi];
          if(ivecpot==1){
            VecPot_R[i + j*nr + k*nr*nphi] = vecpot[0][k + j*ntheta + i*ntheta*nphi];
            VecPot_Phi[i + j*nr + k*nr*nphi] = vecpot[1][k + j*ntheta + i*ntheta*nphi];
            VecPot_Z[i + j*nr + k*nr*nphi] = vecpot[2][k + j*ntheta + i*ntheta*nphi];
          }
          if (ibootstrap==1 || ibootstrap==2 || ibootstrap==3){
            JpdotB_3d  [i + j*nr + k*nr*nphi] = jdotb[k + j*ntheta + i*ntheta*nphi];
          }
        }
      }
    }

    // Toroidal Flux
    std::cerr << "Calculating toroidal flux" << std::endl;
    double* dPsids = new double[nr*nphi];
    
    for(int i=0; i<nr; i++) {
      for(int j=0; j<nphi; j++) {
        dPsids[i+j*nr] = 0.;
        for(int k=0; k<ntheta; k++) {
          double dRdi, dZdi, dRdk, dZdk;
          int ijk = i + j*nr + k*nr*nphi;

          if(i==0) {
            dRdi = R[ijk + 1] - R[ijk];
            dZdi = Z[ijk + 1] - Z[ijk];
          } else if(i==nr-1) {
            dRdi = R[ijk] - R[ijk-1];
            dZdi = Z[ijk] - Z[ijk-1];
          } else {
            dRdi = (R[ijk+1] - R[ijk-1])/2.;
            dZdi = (Z[ijk+1] - Z[ijk-1])/2.;
          };
          int km = (k==0 ? ntheta-1 : k-1);
          int kp = (k==ntheta-1 ? 0 : k+1);
          dRdk = (R[i + j*nr + kp*nr*nphi] - R[i + j*nr + km*nr*nphi])/2.;
          dZdk = (Z[i + j*nr + kp*nr*nphi] - Z[i + j*nr + km*nr*nphi])/2.;

          // flux through dtheta x dr
          dPsids[i+j*nr] += BPhi[ijk]*(dRdi*dZdk - dRdk*dZdi);
        }
      }
    }

    // calculate toroidal flux
    double fluxt = 0.;
    for(int i=0; i<nr; i++) {
      dpsi_t[i] = 0.;

      // calculate mean toroidal flux
      for(int j=0; j<nphi; j++)
        dpsi_t[i] += dPsids[i+j*nr];
        dpsi_t[i] /= nphi;
        fluxt += dpsi_t[i];
        psi_t[i] = fluxt;

      // calculate variance in toroidal flux
      double fluxt_var = 0.;
      for(int j=0; j<nphi; j++) {
        double fluxt_j = 0.;
        for(int ii=0; ii<=i; ii++)
          fluxt_j += dPsids[ii+j*nr];

          fluxt_var += (fluxt_j - psi_t[i])*(fluxt_j - psi_t[i]);
        }
      fluxt_var /= nphi;

      std::cerr << "Toroidal flux at surface " << i << " = "
          << psi_t[i] << " +/- " << sqrt(fluxt_var) << std::endl;
    }


    // Poloidal flux
    std::cerr << "Calculating poloidal flux" << std::endl;
    double* dpsids = new double[nr*ntheta];
    double* dpsix  = new double[nr*ntheta];
    double* dAds = new double[nr*ntheta];
    double* dVds = new double[nr*ntheta];

    for(int i=0; i<nr; i++) {
      for(int k=0; k<ntheta; k++) {
        dpsids[i+k*nr] = 0.;
        dpsix[i+k*nr] = 0.;
        dAds[i+k*nr] = 0.;
        dVds[i+k*nr] = 0.0;

        // Do integral over phi
        for(int j=0; j<nphi; j++) {
          double dRdi, dZdi, dRdj, dZdj, dRdk, dZdk, dphidj;
          int ijk = i + j*nr + k*nr*nphi;

          if(i==0) {
            dRdi = R[ijk + 1] - R[ijk];
            dZdi = Z[ijk + 1] - Z[ijk];
          } else if(i==nr-1) {
            dRdi = R[ijk] - R[ijk - 1];
            dZdi = Z[ijk] - Z[ijk - 1];
          } else {
            dRdi = (R[ijk+1] - R[ijk-1])/2.;
            dZdi = (Z[ijk+1] - Z[ijk-1])/2.;
          };
          int jm = (j==0 ? nphi-1 : j-1);
          int jp = (j==nphi-1 ? 0 : j+1);
          dRdj = (R[i + jp*nr + k*nr*nphi] - R[i + jm*nr + k*nr*nphi])/2.;
          dZdj = (Z[i + jp*nr + k*nr*nphi] - Z[i + jm*nr + k*nr*nphi])/2.;
          dphidj = 2.*M_PI / nphi;

          int km = (k==0 ? ntheta-1 : k-1);
          int kp = (k==ntheta-1 ? 0 : k+1);
          dRdk = (R[i + j*nr + kp*nr*nphi] - R[i + j*nr + km*nr*nphi])/2.;
          dZdk = (Z[i + j*nr + kp*nr*nphi] - Z[i + j*nr + km*nr*nphi])/2.;

          // B . dphi x dr (poloidal flux)
          dpsids[i + k*nr] -= R[ijk]*dphidj*
          ( BZ[ijk] * dRdi - BR[ijk] * dZdi) - 
            BPhi[ijk] * (dRdi*dZdj - dRdj*dZdi);
        
          // B . dtheta x dphi (normal flux)
          dpsix[i + k*nr] -= R[ijk]*dphidj*
          ( BR[ijk] * dZdk - BZ[ijk] * dRdk) -
          BPhi[ijk] * (dRdj*dZdk - dRdk*dZdj);

          dAds[i + k*nr] += Jac[ijk] / sqrt(dRdk*dRdk + dZdk*dZdk);
          dVds[i + k*nr] += Jac[ijk];
        }
      }
    }

    // calculate poloidal flux
    double fluxp = 0.;
    for(int i=0; i<nr; i++) {
      dpsi_p[i] = 0.;

      // calculate mean poloidal flux
      for(int k=0; k<ntheta; k++)
        dpsi_p[i] += dpsids[i+k*nr];
        dpsi_p[i] /= ntheta;
        fluxp += dpsi_p[i];   
        psi_p[i] = fluxp;

      // calculate variance in poloidal flux
      double fluxp_var = 0.;
      for(int k=0; k<ntheta; k++) {
        double fluxp_k = 0.;
        for(int ii=0; ii<=i; ii++)
          fluxp_k += dpsids[ii+k*nr];
          
          fluxp_var += (fluxp_k - psi_p[i])*(fluxp_k - psi_p[i]);
        }
        fluxp_var /= ntheta;

        std::cerr << "Poloidal flux at surface " << i << " = "
          << psi_p[i] << " +/- " << sqrt(fluxp_var) << std::endl;
    }

    // Sanity check: calculate normal flux through surfaces
    // (should be zero)
    for(int i=0; i<nr; i++) {
      double flux_norm = 0;
      for(int k=0; k<ntheta; k++) {
        flux_norm += dpsix[i + k*nr];
      }
      std::cerr << "Net flux through surface " << i << ": "
          << flux_norm << std::endl;
    }

    // Calculate q
    for(int i=0; i<nr; i++) {
      q[i] = dpsi_t[i] / dpsi_p[i];
      std::cerr << "q at surface " << i << ": "
          << q[i] << std::endl;
    }


    // Calculate volume inside each surface
    double volume = 0;
    for(int i=0; i<nr; i++) {
      for(int k=0; k<ntheta; k++) {
        volume += dVds[i + k*nr];
      }
  
      std::cerr << "Net volume inside surface " << i << ": "
          << volume << std::endl;
    }


    double ion_mass = 2.;  // TODO: generalize this
    int version = 4;


    // convert phi to degrees
    for(int i=0; i<nphi; i++) { 
      phi[i] = 180.*phi[i]/M_PI;
      //    std::cerr << "phi[" << i << "] = " << phi[i] << std::endl;
    }

    // Create netcdf file
    std::cerr << "Writing netcdf file" << std::endl;

    int ncid;
    result = nc_create("neo_input.nc", NC_CLOBBER, &ncid);
    if(result != 0) {
      std::cerr << "Error opening netcdf file\n" << nc_strerror(result)
          << std::endl;
    }

    //  std::cerr << " attributes..." << std::endl;  
    nc_put_att_int(ncid, NC_GLOBAL, "version", NC_SHORT, 1, &version);
    nc_put_att_double(ncid, NC_GLOBAL, "ion_mass", NC_FLOAT, 1, &ion_mass);
    nc_put_att_double(ncid, NC_GLOBAL, "psi_0", NC_FLOAT, 1, &psi0);
    nc_put_att_double(ncid, NC_GLOBAL, "psi_1", NC_FLOAT, 1, &psi1);

    //  std::cerr << " defining variables..." << std::endl;  
    int nr_dimid, np_dimid, nt_dimid, npsi_dimid, single_value_dimid;
    nc_def_dim(ncid, "npsi", nr,   &npsi_dimid);
    nc_def_dim(ncid, "nr", nr,     &nr_dimid);
    nc_def_dim(ncid, "np", ntheta, &np_dimid);
    nc_def_dim(ncid, "nt", nphi,   &nt_dimid);
    nc_def_dim(ncid, "single_value", 1,   &single_value_dimid);

    int phi_id;
    nc_def_var(ncid, "Phi", NC_FLOAT, 1, &nt_dimid, &phi_id);

    int axis_r_id;
    int axis_z_id;
    nc_def_var(ncid, "axis_R", NC_FLOAT, 1, &nt_dimid, &axis_r_id); // ver 4
    nc_def_var(ncid, "axis_Z", NC_FLOAT, 1, &nt_dimid, &axis_z_id); // ver 4
    
    int q_id, psi_id, psi_t_id, dpsi_t_id, psi_p_id, dpsi_p_id;
    nc_def_var(ncid, "q", NC_FLOAT, 1, &nr_dimid, &q_id);
    nc_def_var(ncid, "psi", NC_FLOAT, 1, &nr_dimid, &psi_id);
    nc_def_var(ncid, "Psi_t", NC_FLOAT, 1, &nr_dimid, &psi_t_id);    // ver 3
    nc_def_var(ncid, "dPsi_t", NC_FLOAT, 1, &nr_dimid, &dpsi_t_id);  // ver 3
    nc_def_var(ncid, "Psi_p", NC_FLOAT, 1, &nr_dimid, &psi_p_id);    // ver 4
    nc_def_var(ncid, "dPsi_p", NC_FLOAT, 1, &nr_dimid, &dpsi_p_id);  // ver 4  
    int ne_id, te_id, ni_id, ti_id, psi0_id, temax_id,bx_fa_id,by_fa_id,bz_fa_id,bmag_fa_id,b2_fa_id;
    int jx_fa_id,jy_fa_id,jz_fa_id;
    int JpdotB_fluxavg_id,JpdotB_dndpsi_fluxavg_id,JpdotB_dtedpsi_fluxavg_id,JpdotB_dtidpsi_fluxavg_id;
    int JpdotB_L31_fluxavg_id,JpdotB_L32_fluxavg_id,JpdotB_L34_fluxavg_id,JpdotB_alpha_fluxavg_id;
    nc_def_var(ncid, "psi0", NC_FLOAT, 1, &npsi_dimid, &psi0_id);
    nc_def_var(ncid, "ne0",  NC_FLOAT, 1, &npsi_dimid, &ne_id);
    nc_def_var(ncid, "Te0",  NC_FLOAT, 1, &npsi_dimid, &te_id);
    nc_def_var(ncid, "ni0",  NC_FLOAT, 1, &npsi_dimid, &ni_id);
    nc_def_var(ncid, "Ti0",  NC_FLOAT, 1, &npsi_dimid, &ti_id);
    nc_def_var(ncid, "bx_fa",  NC_FLOAT, 1, &npsi_dimid, &bx_fa_id);
    nc_def_var(ncid, "by_fa",  NC_FLOAT, 1, &npsi_dimid, &by_fa_id);
    nc_def_var(ncid, "bz_fa",  NC_FLOAT, 1, &npsi_dimid, &bz_fa_id);
    nc_def_var(ncid, "bmag_fa",  NC_FLOAT, 1, &npsi_dimid, &bmag_fa_id);
    nc_def_var(ncid, "b2_fa",  NC_FLOAT, 1, &npsi_dimid, &b2_fa_id);  
    nc_def_var(ncid, "jx_fa",  NC_FLOAT, 1, &npsi_dimid, &jx_fa_id);
    nc_def_var(ncid, "jy_fa",  NC_FLOAT, 1, &npsi_dimid, &jy_fa_id);
    nc_def_var(ncid, "jz_fa",  NC_FLOAT, 1, &npsi_dimid, &jz_fa_id);
    if (ibootstrap==1 || ibootstrap ==2 ||ibootstrap ==3){
    nc_def_var(ncid, "JpdotB_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_fluxavg_id);
    nc_def_var(ncid, "JpdotB_dndpsi_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_dndpsi_fluxavg_id);
    nc_def_var(ncid, "JpdotB_dtedpsi_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_dtedpsi_fluxavg_id);
    nc_def_var(ncid, "JpdotB_dtidpsi_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_dtidpsi_fluxavg_id);
    }
    if (ibootstrap ==3){
    nc_def_var(ncid, "JpdotB_L31_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_L31_fluxavg_id);
    nc_def_var(ncid, "JpdotB_L32_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_L32_fluxavg_id);
    nc_def_var(ncid, "JpdotB_L34_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_L34_fluxavg_id);
    nc_def_var(ncid, "JpdotB_alpha_fluxavg",  NC_FLOAT, 1, &npsi_dimid, &JpdotB_alpha_fluxavg_id);
    }
    nc_def_var(ncid, "temax",  NC_FLOAT, 0, &single_value_dimid, &temax_id); 

    int r_id, z_id, jac_id, br_id, bphi_id, bz_id,vecpot_r_id,vecpot_phi_id,vecpot_z_id;
    int dims[3];
    dims[2] = nr_dimid;
    dims[1] = nt_dimid;
    dims[0] = np_dimid;
    nc_def_var(ncid, "R", NC_FLOAT, 3, dims, &r_id);
    nc_def_var(ncid, "Z", NC_FLOAT, 3, dims, &z_id);
    nc_def_var(ncid, "Jac", NC_FLOAT, 3, dims, &jac_id);    // added in version 3
    nc_def_var(ncid, "B_R",   NC_FLOAT, 3, dims, &br_id);   // added in version 3
    nc_def_var(ncid, "B_Phi", NC_FLOAT, 3, dims, &bphi_id); // added in version 3
    nc_def_var(ncid, "B_Z",   NC_FLOAT, 3, dims, &bz_id);   // added in version 3
    if(ivecpot==1){  
      nc_def_var(ncid, "VecPot_R",   NC_FLOAT, 3, dims, &vecpot_r_id);   // added in version 3
      nc_def_var(ncid, "VecPot_Phi", NC_FLOAT, 3, dims, &vecpot_phi_id); // added in version 3
      nc_def_var(ncid, "VecPot_Z",   NC_FLOAT, 3, dims, &vecpot_z_id);   // added in version 3
    }
    
    int te3d_id, ti3d_id, ne3d_id, ni3d_id, JpdotB3d_id, p3d_id;
    nc_def_var(ncid, "te_3d", NC_FLOAT, 3, dims, &te3d_id);   // added in version 3
    nc_def_var(ncid, "ti_3d", NC_FLOAT, 3, dims, &ti3d_id); // added in version 3
    nc_def_var(ncid, "ne_3d", NC_FLOAT, 3, dims, &ne3d_id);   // added in version 3
    nc_def_var(ncid, "ni_3d", NC_FLOAT, 3, dims, &ni3d_id); // added in version 3
    nc_def_var(ncid, "p_3d", NC_FLOAT, 3, dims, &p3d_id);
    if (ibootstrap==1 || ibootstrap ==2 ||ibootstrap ==3){
    nc_def_var(ncid, "JpdotB_3d", NC_FLOAT, 3, dims, &JpdotB3d_id); // added in version 3
    }
    nc_enddef(ncid);

    //  std::cerr << " writing variables..." << std::endl;  
    nc_put_var_double(ncid, phi_id, phi);
    nc_put_var_double(ncid, axis_r_id, axis_3d[0]);
    nc_put_var_double(ncid, axis_z_id, axis_3d[2]);
    nc_put_var_double(ncid, q_id, q);
    nc_put_var_double(ncid, psi_id, psi_surf);
    nc_put_var_double(ncid, psi0_id, psi_surf);
    nc_put_var_double(ncid, ne_id, ne);
    nc_put_var_double(ncid, te_id, te);
    nc_put_var_double(ncid, ni_id, ni);
    nc_put_var_double(ncid, ti_id, ti);
    nc_put_var_double(ncid, bx_fa_id, bx_fa);
    nc_put_var_double(ncid, by_fa_id, by_fa);
    nc_put_var_double(ncid, bz_fa_id, bz_fa);
    nc_put_var_double(ncid, bmag_fa_id, bmag_fa);
    nc_put_var_double(ncid, b2_fa_id, b2_fa);
    nc_put_var_double(ncid, jx_fa_id, jx_fa);
    nc_put_var_double(ncid, jy_fa_id, jy_fa);
    nc_put_var_double(ncid, jz_fa_id, jz_fa);
    if (ibootstrap==1 || ibootstrap ==2 ||ibootstrap ==3){
      nc_put_var_double(ncid, JpdotB_fluxavg_id, JpdotB_fluxavg);
      nc_put_var_double(ncid, JpdotB_dndpsi_fluxavg_id, JpdotB_dndpsi_fluxavg);
      nc_put_var_double(ncid, JpdotB_dtedpsi_fluxavg_id, JpdotB_dtedpsi_fluxavg);
      nc_put_var_double(ncid, JpdotB_dtidpsi_fluxavg_id, JpdotB_dtidpsi_fluxavg);
    }
    if (ibootstrap ==3){
      nc_put_var_double(ncid, JpdotB_L31_fluxavg_id, JpdotB_L31_fluxavg);
      nc_put_var_double(ncid, JpdotB_L32_fluxavg_id, JpdotB_L32_fluxavg);
      nc_put_var_double(ncid, JpdotB_L34_fluxavg_id, JpdotB_L34_fluxavg);
      nc_put_var_double(ncid, JpdotB_alpha_fluxavg_id, JpdotB_alpha_fluxavg);
    }
    nc_put_var_double(ncid, temax_id, &te_max);
    nc_put_var_double(ncid, r_id, R);
    nc_put_var_double(ncid, z_id, Z);
    nc_put_var_double(ncid, jac_id, Jac);
    nc_put_var_double(ncid, br_id, BR);
    nc_put_var_double(ncid, bphi_id, BPhi);
    nc_put_var_double(ncid, bz_id, BZ);

    if(ivecpot==1){
      nc_put_var_double(ncid, vecpot_r_id, VecPot_R);
      nc_put_var_double(ncid, vecpot_phi_id, VecPot_Phi);
      nc_put_var_double(ncid, vecpot_z_id, VecPot_Z);
    }
    nc_put_var_double(ncid, te3d_id, Te_3d);
    nc_put_var_double(ncid, ti3d_id, Ti_3d);
    nc_put_var_double(ncid, ne3d_id, Ne_3d);
    nc_put_var_double(ncid, ni3d_id, Ni_3d);
    nc_put_var_double(ncid, p3d_id, P_3d);
    if (ibootstrap==1 || ibootstrap ==2 ||ibootstrap ==3){
      nc_put_var_double(ncid, JpdotB3d_id, JpdotB_3d);
    }

    nc_put_var_double(ncid, psi_t_id, psi_t);
    nc_put_var_double(ncid, dpsi_t_id, dpsi_t);
    nc_put_var_double(ncid, psi_p_id, psi_p);
    nc_put_var_double(ncid, dpsi_p_id, dpsi_p);

    nc_close(ncid);

    std::cerr << "Done writing netcdf file." << std::endl;

  // ========================================================
  // BOOTSTRAP CALCULATIONS
  // ========================================================

  // Derived Grids and Gradients
  double* dtebydpsit = new double[nr];
  double tnorm = 1000.0;

  std::vector<double> psi_norm_new(nr, 0.0);

  if (ibootstrap == 2) {
      // Mode 2: x is Te, but psi_norm (used internally) is (temax - Te)/temax
      for (int i = 0; i < nr; i++) {
          psi_norm_new[i] = te[i];
      }
  } else if (ibootstrap == 3) {
      // Mode 3: x is (temax - Te)/temax
      for (int i = 0; i < nr; i++) {
          psi_norm_new[i] = (te_max - te[i]) / te_max;
      }
  } else if (ibootstrap == 1) {
      // Mode 1: x is the linear psi grid
      for (int i = 0; i < nr; i++) {
          psi_norm_new[i] = (psi_end - psi_start) * ((double)i / (double)nr) + psi_start;
      }
  }

  for(int i = 0; i < nr; i++) {
      if (ibootstrap == 1 || ibootstrap == 2) {
          if (i == 0) 
              dtebydpsit[i] = (te[i+1] - te[i]) / (psi_t[i+1] - psi_t[i]) / tnorm;
          else if (i == nr - 1)
              dtebydpsit[i] = (te[i] - te[i-1]) / (psi_t[i] - psi_t[i-1]) / tnorm;
          else
              dtebydpsit[i] = (te[i+1] - te[i-1]) / (psi_t[i+1] - psi_t[i-1]) / tnorm;
      } 
      else if (ibootstrap == 3) {
          // Find d(psi_norm)/d(psi_t)
          if (i == 0)
              dtebydpsit[i] = (psi_norm_new[i+1] - psi_norm_new[i]) / (psi_t[i+1] - psi_t[i]);
          else if (i == nr - 1)
              dtebydpsit[i] = (psi_norm_new[i] - psi_norm_new[i-1]) / (psi_t[i] - psi_t[i-1]);
          else
              dtebydpsit[i] = (psi_norm_new[i+1] - psi_norm_new[i-1]) / (psi_t[i+1] - psi_t[i-1]);
      }
  }

  // ========================================================
  // Magnetic Extrema
  // ========================================================
  double* Bmax_fpsi = new double[nr];
  double* Bmin_fpsi = new double[nr];
  double* Bmag_3d = new double[nr * nphi * ntheta];

  for (int i = 0; i < nr; i++) {
      Bmax_fpsi[i] = -1e30;
      Bmin_fpsi[i] = 1e30;
      for (int j = 0; j < nphi; j++) {
          for (int k = 0; k < ntheta; k++) {
              int idx = i + j*nr + k*nr*nphi;
              
              double b_val = std::sqrt(std::pow(BR[i + j*nr + k*nr*nphi], 2) + 
                                      std::pow(BPhi[i + j*nr + k*nr*nphi], 2) + 
                                      std::pow(BZ[i + j*nr + k*nr*nphi], 2));
              Bmag_3d[idx] = b_val;
              if (b_val > Bmax_fpsi[i]) Bmax_fpsi[i] = b_val;
              if (b_val < Bmin_fpsi[i]) Bmin_fpsi[i] = b_val;
          }
      }
  }


  // ========================================================
  // Execute Landreman Integrals
  // ========================================================
  CurrentIntegrals landreman = calculate_current_integrals(
      nr, nphi, ntheta, nfp, Nzp, R, Z, BR, BPhi, BZ, q
  );

  // ========================================================
  // Execute f_trap
  // ========================================================

  std::vector<double> ftrap_vec = calc_ftrap(
      nr, nphi, ntheta, lambda_count, Bmax_fpsi, b2_fa, Bmag_3d, Jac, psi_norm_new.data(), tol_trap
  );

  // ========================================================
  // Zeff & Collisionalities
  // ========================================================
  double Ze, Zi, Zeff;
  calculate_Zcharge(ne, ni, nr, Zcharge_e_input, Ze, Zi, Zeff);
  std::cerr << " Z_effective = " << Zeff << "\n";

  // Prep Bmag_inversefa (1/<B>)
  double* Bmag_inv = new double[nr];
  for(int i=0; i<nr; i++) Bmag_inv[i] = 1.0 / std::max(bmag_fa[i], 1e-30);

  CollisionResults coll = calculate_collision_frequencies(
      nr, nphi, ntheta, ne, ni, te, ti, Ze, Zi, Bmax_fpsi, Bmin_fpsi, q, 
      landreman.GplusiI_fa.data(), Bmag_inv, R,UseBfore, Nzp
  );

  // ========================================================
  // Sauter/Redl Coefficients
  // ========================================================

  NeoCoeffs bootstrap_coeffs;


  if (UseSauter0Redl1 == 0) {
      
      SauterResults sr = calculate_sauter_terms(nr, Zeff, ftrap_vec, coll.nu_e_star, coll.nu_i_star);
      bootstrap_coeffs.L31 = sr.L31;
      bootstrap_coeffs.L32 = sr.L32;
      bootstrap_coeffs.L34 = sr.L34;
      bootstrap_coeffs.alpha_nu_i = sr.alpha_nu_i;
      
  } else if (UseSauter0Redl1 == 1) {
      RedlResults rr = calculate_redl_terms(nr, Zeff, ftrap_vec,  coll.nu_e_star, coll.nu_i_star);
      bootstrap_coeffs.L31 = rr.L31;
      bootstrap_coeffs.L32 = rr.L32;
      bootstrap_coeffs.L34 = rr.L34;
      bootstrap_coeffs.alpha_nu_i = rr.alpha_nu_i;
  }



  // ========================================================
  // File Output
  // ========================================================
  std::string fname;
  if (ibootstrap == 1) fname = "ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G";
  else if (ibootstrap == 2) fname = "ProfileJBSCoeff_Psi_L31_32_34_alpha_B2_dtedpsit_G";
  else if (ibootstrap == 3) fname = "ProfileJBSCoeff_Tenorm_L31_32_34_alpha_B2_dtedpsit_G_ft_qR_e_temax";

  if (ibootstrap >= 1 && ibootstrap <= 3) {
      std::ofstream outfile(fname);
      if (outfile.is_open()) {
          // Headers
          if (ibootstrap == 1) {
              outfile << std::setw(9) << "Te" << " " << std::setw(9) << "L31" << " " << std::setw(9) << "L32" << " "
                      << std::setw(9) << "L34" << " " << std::setw(9) << "alpha" << " " << std::setw(9) << "1/<B^2>" << " "
                      << std::setw(14) << "dtebydpsit*2*pi" << " " << std::setw(12) << "Gbar/(i-N)" << "\n";
          } else if (ibootstrap == 2) {
              outfile << std::setw(9) << "Psi" << " " << std::setw(9) << "L31" << " " << std::setw(9) << "L32" << " "
                      << std::setw(9) << "L34" << " " << std::setw(9) << "alpha" << " " << std::setw(9) << "1/<B^2>" << " "
                      << std::setw(14) << "dtebydpsit*2*pi" << " " << std::setw(12) << "Gbar/(i-N)" << "\n";
          } else if (ibootstrap == 3) {
              outfile << std::setw(9) << "Psi" << " " << std::setw(9) << "L31" << " " << std::setw(9) << "L32" << " "
                      << std::setw(9) << "L34" << " " << std::setw(9) << "alpha" << " " << std::setw(9) << "1/<B^2>" << " "
                      << std::setw(37) << "dtnormdpsit*2*pi (multiply by -temax)" << " " << std::setw(12) << "Gbar/(i-N)" << " "
                      << std::setw(9) << "ftrap" << " " << std::setw(9) << "qR" << " " << std::setw(9) << "epsilon" << " "
                      << std::setw(9) << "temax>" << "\n";
          }

          
          outfile << std::fixed << std::setprecision(8);
          for (int i = 1; i < nr; i++) {
              
              outfile << psi_norm_new[i] << " " 
                      << bootstrap_coeffs.L31[i] << " " 
                      << bootstrap_coeffs.L32[i] << " " 
                      << bootstrap_coeffs.L34[i] << " " 
                      << bootstrap_coeffs.alpha_nu_i[i] << " " 
                      << (1.0 / b2_fa[i]) << " " 
                      << (dtebydpsit[i] * 2.0 * M_PI) << " " 
                      << landreman.Gbar_by_iminusN[i];

              if (ibootstrap == 3) {
                  outfile << " " << ftrap_vec[i] << " " 
                          << coll.qR[i] << " " 
                          << coll.epsilon[i] << " " 
                          << (te_max / 1000.0);
              }
              outfile << "\n";
          }
          outfile.close();
          std::cerr << "Successfully wrote output to " << fname << std::endl;
      } else {
          std::cerr << "Error: Could not open file " << fname << " for writing." << std::endl;
      }
  }


    std::cerr << "Freeing memory." << std::endl;

}
// Cleanup local dynamic memory


  delete[] jac;
  delete[] theta;
  delete[] phi;
  delete[] q;
  delete[] psi_surf;
  delete[] psi_t;
  delete[] dpsi_t;
  delete[] psi_p;
  delete[] dpsi_p;
  delete[] ne;
  delete[] te;
  delete[] ni;
  delete[] ti;
  delete[] b[0];
  delete[] b[1];
  delete[] b[2];
  delete[] b;
  delete[] b_surf;
  delete[] currden[0];
  delete[] currden[1];
  delete[] currden[2];
  delete[] path_plane;
  delete[] path_surf;
  delete[] path[0];
  delete[] path[1];
  delete[] path[2];
  delete[] path;
  delete[] axis_3d[0];
  delete[] axis_3d[1];
  delete[] axis_3d[2];

  delete[] telec;
  delete[] tion;
  delete[] nelec;
  delete[] nion;
  if (ibootstrap==1 || ibootstrap ==2 ||ibootstrap ==3){
    delete[] jdotb;
  }

  // Deallocate fields
  src0->deallocate_search_hint(&h);
  delete_sources();

  std::cerr << "Done." << std::endl;
  MPI_Finalize();
  return 0;
}


bool create_source(const int type, const int argc, const std::string argv[]) 
{
  fio_source* src;
  fio_option_list fopt;
  int result;
  double slice_time = 0.;
  int timeslice = 0;

  std::cerr << "Opening source\n";

  switch(type) {
  case(FIO_M3DC1_SOURCE):
    std::cerr << " Type: M3DC1\n";
    src = new m3dc1_source();
    if(argc>=1) {
      result = src->open(argv[0].c_str());
    } else {
      result = src->open("C1.h5");
    }
    if(result != FIO_SUCCESS) {
      std::cerr << "Error opening file" << std::endl;
      delete(src);
      return false;
    };
    // set options for fields obtained from this source
    src->get_field_options(&fopt);
    fopt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
    if(argc<2) fopt.set_option(FIO_TIMESLICE,-1);  // default to timeslice=-1
    break;

  case(FIO_GEQDSK_SOURCE):
    std::cerr << " Type: GEQDSK\n";
    src = new geqdsk_source();
    if(argc>=1) {
      if((result = src->open(argv[0].c_str())) != FIO_SUCCESS) {
	std::cerr << "Error opening file" << std::endl;
	delete(src);
	return false;
      };
    } else {
      std::cerr << "Filename must be provided for geqdsk files." << std::endl;
      delete(src);
      return false;
    }
    src->get_field_options(&fopt);
    break;

  case(FIO_GPEC_SOURCE):
    std::cerr << " Type: GPEC\n";
    src = new gpec_source();
    if(argc>=1) {
      if((result = src->open(argv[0].c_str())) != FIO_SUCCESS) {
	std::cerr << "Error opening file" << std::endl;
	delete(src);
	return false;
      };
    } else {
      std::cerr << "Directory must be provided for gpec files." << std::endl;
      delete(src);
      return false;
    }
    src->get_field_options(&fopt);
    break;

  default:
    return false;
  }

  if(argc>=2) {
    timeslice = atoi(argv[1].c_str());
    std::cerr << " Time slice: " << timeslice << std::endl;
    fopt.set_option(FIO_TIMESLICE, timeslice);
    if(src->get_slice_time(timeslice, &slice_time)==FIO_SUCCESS) {
      std::cerr << "  Slice time = " << slice_time << std::endl;
    }
  }
  if(argc>=3) {
    double scale = atof(argv[2].c_str());
    std::cerr << " Linear scale factor: " << scale << std::endl;
    fopt.set_option(FIO_LINEAR_SCALE, scale);
  }
  if(argc>=4) {
    double phase = atof(argv[3].c_str())*M_PI/180.;
    std::cerr << " Phase: " << phase << std::endl;
    fopt.set_option(FIO_PHASE, phase);
  }

  // Allocate search hint
  fio_hint hint;
  result = src->allocate_search_hint(&hint);
  if(result != FIO_SUCCESS) {
    std::cerr << "Warning: couldn't allocate search hint" << std::endl;
  }



  if(sources.size()==0) {
    // If this is the first source being read, use this source to 
    // determine equilibrium fields and magnetic axis.

    fio_series *magaxis[2];
    fio_series *magaxis_psi, *lcfs_psi;

    // Find magnetic axis
    if((result = src->get_series(FIO_MAGAXIS_R, &magaxis[0])) != FIO_SUCCESS) {
      std::cerr << "Couldn't load MAGAXIS_R" << std::endl;
      delete(src);
      return result;
    }
    if((result = src->get_series(FIO_MAGAXIS_Z, &magaxis[1])) != FIO_SUCCESS) {
      std::cerr << "Couldn't load MAGAXIS_Z" << std::endl;
      delete(src);
      return result;
    }
    if((result = magaxis[0]->eval(slice_time, &(axis[0]))) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate MAGAXIS_R" << std::endl;
	delete(src);
	return result;
    }
    if((result = magaxis[1]->eval(slice_time, &(axis[2]))) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate MAGAXIS_Z" << std::endl;
	delete(src);
	return result;
      }
    delete(magaxis[0]);
    delete(magaxis[1]);


    std::cerr << "Magnetic axis found at (" 
	      << axis[0] << ", " << axis[2] << ")." << std::endl;


    if((result = src->get_series(FIO_MAGAXIS_PSI, &magaxis_psi)) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't load MAGAXIS_PSI" << std::endl;
	delete(src);
	return result;
      }
    if((result = src->get_series(FIO_LCFS_PSI, &lcfs_psi)) != FIO_SUCCESS) 
      {
	std::cerr << "Couldn't load LCFS_PSI" << std::endl;
	delete(src);
	return result;
      }
    if((result = magaxis_psi->eval(slice_time, &psi0)) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate MAGAXIS_PSI" << std::endl;
	delete(src);
	return result;
      }
    if((result = lcfs_psi->eval(slice_time, &psi1)) != FIO_SUCCESS)
      {
	std::cerr << "Couldn't evaluate LCFS_PSI" << std::endl;
	delete(src);
	return result;
      }
    delete(magaxis_psi);
    delete(lcfs_psi);
    
    std::cerr << "psi0 = " << psi0 << std::endl;
    std::cerr << "psi1 = " << psi1 << std::endl;


    fopt.set_option(FIO_PART, FIO_TOTAL);
    
  } else {
    fopt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
  }
  
  



  // Read fields
  fio_field* field;
  if (ibootstrap==1 || ibootstrap ==2 ||ibootstrap ==3){
   result = src->get_field(FIO_JBS, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field " << std::endl;
    delete(src);
    return result;
   }
   JpdotB.add_field(field, FIO_ADD, 1., hint);

   result = src->get_field(FIO_JBS_dndpsi, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field dndpsi " << std::endl;
    delete(src);
    return result;
   }
   JpdotB_dndpsi.add_field(field, FIO_ADD, 1., hint);

   result = src->get_field(FIO_JBS_dtedpsi, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field dtedpsi" << std::endl;
    delete(src);
    return result;
   }
   JpdotB_dtedpsi.add_field(field, FIO_ADD, 1., hint);

   result = src->get_field(FIO_JBS_dtidpsi, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field dtidpsi" << std::endl;
    delete(src);
    return result;
   }
   JpdotB_dtidpsi.add_field(field, FIO_ADD, 1., hint);

  } 

  if (ibootstrap ==3){
   result = src->get_field(FIO_JBS_L31, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field L31 " << std::endl;
    delete(src);
    return result;
   }
   JpdotB_L31.add_field(field, FIO_ADD, 1., hint);

   result = src->get_field(FIO_JBS_L32, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field L32 " << std::endl;
    delete(src);
    return result;
   }
   JpdotB_L32.add_field(field, FIO_ADD, 1., hint);

   result = src->get_field(FIO_JBS_L34, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field L34" << std::endl;
    delete(src);
    return result;
   }
   JpdotB_L34.add_field(field, FIO_ADD, 1., hint);

   result = src->get_field(FIO_JBS_alpha, &field, &fopt);
   if(result != FIO_SUCCESS) {
    std::cerr << "Error opening bootstrap current field alpha" << std::endl;
    delete(src);
    return result;
   }
   JpdotB_alpha.add_field(field, FIO_ADD, 1., hint);

  } 

  if(ibootstrap !=1 ||ibootstrap!=2 ||ibootstrap!=3) {
    std::cerr << "With bootstrap model, ibootstrap=" 
          << ibootstrap
          << std::endl;
  }

  result = src->get_field(FIO_MAGNETIC_FIELD, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening magnetic field field" << std::endl;
    delete(src);
    return result;
  }
  mag.add_field(field, FIO_ADD, 1., hint);

  result = src->get_field(FIO_CURRENT_DENSITY, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening current density field" << std::endl;
    delete(src);
    return result;
  }
  current_density.add_field(field, FIO_ADD, 1., hint);

  result = src->get_field(FIO_VECTOR_POTENTIAL, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening vector potential field" << std::endl;
    delete(src);
    return result;
  }
  vector_potential.add_field(field, FIO_ADD, 1., hint);
  
  result = src->get_field(FIO_POLOIDAL_FLUX_NORM, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening psi norm field" << std::endl;
    delete(src);
    return result;
  }
  psin.add_field(field, FIO_ADD, 1., hint);


  // open fields
  fopt.set_option(FIO_SPECIES, FIO_ELECTRON);
  result = src->get_field(FIO_TEMPERATURE, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron temperature field" << std::endl;
    delete(src);
    return result;
  }
  electron_temperature.add_field(field, FIO_ADD, 1., hint);

  result = src->get_field(FIO_DENSITY, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening electron density field" << std::endl;
    delete(src);
    return result;
  }
  electron_density.add_field(field, FIO_ADD, 1., hint);

  fopt.set_option(FIO_SPECIES, FIO_MAIN_ION);
  result = src->get_field(FIO_DENSITY, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening ion density field" << std::endl;
    delete(src);
    return result;
  }
  ion_density.add_field(field, FIO_ADD, 1., hint);

  result = src->get_field(FIO_TEMPERATURE, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening ion temperature field" << std::endl;
    delete(src);
    return result;
  }
  ion_temperature.add_field(field, FIO_ADD, 1., hint);  

  result = src->get_field(FIO_TOTAL_PRESSURE, &field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening pressure field" << std::endl;
    delete(src);
    return result;
  }
  pressure.add_field(field, FIO_ADD, 1., hint);

  // Add source to list
  sources.push_back(src);
    

  return FIO_SUCCESS;
}




int process_command_line(int argc, char* argv[])
{
  const int max_args = 4;
  const int num_opts = 18;
  std::string arg_list[num_opts] = 
    { "-dl","-dR0", "-bootstrap", "-vecpot","-m3dc1", "-max_step", "-nphi", "-nphi", "-npsi", "-nr",
      "-ntheta", "-psi_end", "-psi_start", "-R0", "-te_start", "-te_end", "-tol", "-hybrid"};
  std::string opt = "";
  std::string arg[max_args];
  int args = 0;
  bool is_opt;
  bool processed = true;
  int result;

  for(int i=1; i<argc; i++) {
    // determine if current cl arg is an option
    is_opt = false;
    for(int j=0; j<num_opts; j++) {
      if(argv[i]==arg_list[j]) {
	is_opt = true;
	break;
      }
    }
    
    if(is_opt) {     // if so, process current option
      if(!processed)
	if((result = process_line(opt, args, arg)) != FIO_SUCCESS) 
	  return result;

      opt = argv[i];
      args = 0;
      processed = false;
      
    } else {         // otherwise, add argument
      if(args >= max_args)
	std::cerr << "Too many arguments for " << opt << std::endl;
      else 
	arg[args++] = argv[i];
    }
  }

  if(!processed)
    if((result = process_line(opt, args, arg)) != FIO_SUCCESS)
      return result;

  return FIO_SUCCESS;

}

int process_line(const std::string& opt, const int argc, const std::string argv[])
{
  bool argc_err = false;


  if(opt=="-dl") {
    if(argc==1) dl_pol = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dR0") {
    if(argc==1) dR0 = atof(argv[0].c_str());
    else argc_err = true;    
  } else if(opt=="-bootstrap") {
    if(argc==1) ibootstrap = atof(argv[0].c_str());
    else ibootstrap = 0; 
    if(ibootstrap !=0 && ibootstrap !=1 && ibootstrap !=2 && ibootstrap !=3){
    std::cerr << "Error: bootstrap option either 0 or 1 or 2 or 3 \n bootstrap="
              << ibootstrap
              << std::endl;
    return FIO_UNSUPPORTED;
    }
  } else if(opt=="-vecpot") {
    if(argc==1) ivecpot = atof(argv[0].c_str());
    else argc_err = true;    
  } else if(opt=="-m3dc1") {
    return create_source(FIO_M3DC1_SOURCE, argc, argv);
  } else if(opt=="-max_step") {
    if(argc==1) max_step = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-nphi") {
    if(argc==1) nphi = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-nr") {
    if(argc==1) nr = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-ntheta") {
    if(argc==1) ntheta = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-psi_end") {
    if(argc==1) psi_end = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-psi_start") {
    if(argc==1) psi_start = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-R0") {
    if(argc==1) R0 = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-te_end") {
    if(argc==1) te_end = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-te_start") {
    if(argc==1) te_start = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-tol") {
    if(argc==1) tol = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-hybrid") {
    if(argc==1) hybrid = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-npsi") {
    std::cerr << "Error: option -npsi is deprecated\n"
	      << " profiles are now only calculated at surface locations"
	      << std::endl;
    return FIO_UNSUPPORTED;
  } else if(opt=="-pert_prof") {
    std::cerr << "Error: option -pert_prof is deprecated\n"
	      << " profiles are now always calculated on perturbed surfaces"
	      << std::endl;
    return FIO_UNSUPPORTED;
  } else {
    std::cerr << "Unrecognized option " << opt << std::endl;
    return FIO_UNSUPPORTED;
  }

  if(argc_err) {
    std::cerr << "Incorrect number of arguments for option " 
	      << opt << std::endl;
    return FIO_UNSUPPORTED;
  }

  return FIO_SUCCESS;
}


void print_usage()
{
  std::cerr << "write_neo_input"
  	  << " -dl <dl>"
	    << " -dR0 <dR0>"
      << " -bootstrap <0/1/2/3>"
      << " -vecpot <0/1>"
	    << " -m3dc1 <m3dc1_source> <time> <scale> <phase>"
	    << " -nphi <nphi>"
	    << " -nr <nr>"
	    << " -ntheta <ntheta>"
	    << " -psi_end <psi_end>"
	    << " -psi_start <psi_start>"	    
      << " -te_end <te_end>"
	    << " -te_start <te_start>"	
	    << std::endl;

  std::cerr
    << "<dl>:           poloidal step size when tracing isosurface\n"
    << "<dR0>:          offset to major radius\n"
    << "<bootstrap 0/1/2/3>:flag to output <j.B> if the bootstrap model is on (1) in M3D-C1\n"
    << "<vecpot 0/1>:   flag to output Vector Potential\n"
    << "<m3dc1_source>: filename of M3D-C1 source file\n"
    << "<nphi>:         number of toroidal points per surface\n"
    << "<nr>:           number of surfaces\n"
    << "<ntheta>:       number of poloidal points per surface\n"
    << "<phase>:        phase factor to apply to fields\n"
    << "<psi_end>:      psi_norm of outermost surface\n"
    << "<psi_start>:    psi_norm of innermost surface\n"
    << "<scale>:        scale factor to apply to linear perturbation\n"
    << "<time>:         timeslice of fields to read\n"
    << "<tol>:          tolerance for finding Te isourface (in eV)\n"
    << "<hybrid>:       hybrid <0/1>\n";
}

void delete_sources()
{
  for(auto it = sources.cbegin(); it != sources.end(); it++) {
    (*it)->close();
    delete(*it);
  }

  sources.clear();
}

CurrentIntegrals calculate_current_integrals(
    int nr, int nphi, int ntheta, int nfp, int Nzp,
    const double* R, const double* Z, 
    const double* BR, const double* BPhi, const double* BZ,
    const double* q) 
{
    CurrentIntegrals res;
    res.Izp_sumj.assign(nr, 0.0);
    res.Izt_sumk.assign(nr, 0.0);
    res.GplusiI_fa.assign(nr, 0.0);
    res.Gbar_by_iminusN.assign(nr, 0.0);

    const double dphidj = 2.0 * M_PI / (nphi * nfp);

    for (int i = 0; i < nr; ++i) {
        // Initialize sums for this specific radial surface
        double Izp_accumulator = 0.0;
        double Izt_accumulator = 0.0;

        for (int j = 0; j < nphi; ++j) {
            for (int k = 0; k < ntheta; ++k) {
                int ijk = i + j * nr + k * nr * nphi;

                // Finite Difference along phi (j) ---
                double dRdj, dZdj;
                if (j == 0) {
                    dRdj = R[i + (j + 1) * nr + k * nr * nphi] - R[ijk];
                    dZdj = Z[i + (j + 1) * nr + k * nr * nphi] - Z[ijk];
                } else if (j == nphi - 1) {
                    dRdj = R[ijk] - R[i + (j - 1) * nr + k * nr * nphi];
                    dZdj = Z[ijk] - Z[i + (j - 1) * nr + k * nr * nphi];
                } else {
                    dRdj = 0.5 * (R[i + (j + 1) * nr + k * nr * nphi] - R[i + (j - 1) * nr + k * nr * nphi]);
                    dZdj = 0.5 * (Z[i + (j + 1) * nr + k * nr * nphi] - Z[i + (j - 1) * nr + k * nr * nphi]);
                }

                // Finite Difference along theta (k) ---
                double dRdk, dZdk;
                if (k == 0) {
                    dRdk = R[i + j * nr + (k + 1) * nr * nphi] - R[ijk];
                    dZdk = Z[i + j * nr + (k + 1) * nr * nphi] - Z[ijk];
                } else if (k == ntheta - 1) {
                    dRdk = R[ijk] - R[i + j * nr + (k - 1) * nr * nphi];
                    dZdk = Z[ijk] - Z[i + j * nr + (k - 1) * nr * nphi];
                } else {
                    dRdk = 0.5 * (R[i + j * nr + (k + 1) * nr * nphi] - R[i + j * nr + (k - 1) * nr * nphi]);
                    dZdk = 0.5 * (Z[i + j * nr + (k + 1) * nr * nphi] - Z[i + j * nr + (k - 1) * nr * nphi]);
                }

                // Izp term summed over k
                double term_zp = (R[ijk] * BPhi[ijk] * dphidj) + (BR[ijk] * dRdj) + (BZ[ijk] * dZdj);
                Izp_accumulator += term_zp;

                // Izt term Izt_sum[i] += Br*dRdk + Bz*dZdk
                double term_zt = (BR[ijk] * dRdk + BZ[ijk] * dZdk);
                Izt_accumulator += term_zt;
            }
        }

        //  Izt_sumk = Izt_sum / nphi / (2*np.pi)
        res.Izt_sumk[i] = Izt_accumulator / (double)nphi / (2.0 * M_PI);
        
        // Izp_sumj = Izp_sum / ntheta / (2*np.pi)
        res.Izp_sumj[i] = Izp_accumulator / (double)ntheta / (2.0 * M_PI);

        res.GplusiI_fa[i] = res.Izp_sumj[i] + (res.Izt_sumk[i] / q[i]);
        res.Gbar_by_iminusN[i] = (res.Izp_sumj[i] + Nzp * res.Izt_sumk[i]) / (1.0 / q[i] - Nzp);
    }

    return res;
}



// --- Trapped Fraction Calculation ---
std::vector<double> calc_ftrap(
    int nr, int nphi, int ntheta, int lambda_count, 
    const double* Bmax_fpsi, const double* B2_fa, const double* Bmag, 
    const double* Jac, const double* psi_eval, double tol) 
{
    std::vector<double> ftrap_out(nr, 0.0);
    
    // 5-pt Gauss-Legendre nodes/weights on [-1,1]
    const double abcissae[5] = {-0.90618, -0.53847, 0.0, 0.53847, 0.9061};
    const double weight[5]   = {0.23693, 0.47863, 0.56889, 0.47863, 0.23693};

    for (int i = 0; i < nr; ++i) {
   
        //  dl = (1.0/Bmax - 0)/lambda_count; lam0 = 0 + dl/2
        double dl = (1.0 / Bmax_fpsi[i]) / (double)lambda_count;
        double a_python = dl / 2.0;
        double b_python = a_python + dl * (lambda_count - 1);
        
        double current_psi = psi_eval[i];
        double integral_sum = 0.0;

        //  Gauss-Legendre Loop
        for (int n = 0; n < 5; ++n) {
            
            double lam_gauss = 0.5 * (b_python - a_python) * abcissae[n] + 0.5 * (b_python + a_python);
            
            double dV = 0.0;
            double acc = 0.0;

            // Flux surface integral 
            for (int j = 0; j < nphi; ++j) {
                for (int k = 0; k < ntheta; ++k) {
                    int idx = i + j * nr + k * nr * nphi;
                    
                    // if (psin <= current_psi + tol) and (psin >= current_psi - tol)
                    if (psi_eval[i] <= (current_psi + tol) && psi_eval[i] >= (current_psi - tol)) {
                        double val = lam_gauss * Bmag[idx];
                        if (val <= 1.0) {
                            acc += std::sqrt(1.0 - val) * Jac[idx];
                            dV += Jac[idx];
                        }
                    }
                }
            }

            if (dV > 0) {
                double denom = acc / dV;
                double f_gauss = lam_gauss / denom;
                integral_sum += weight[n] * f_gauss;
            }
        }


        double integral_gauss = 0.5 * (b_python - a_python) * integral_sum;
        ftrap_out[i] = 1.0 - 0.75 * B2_fa[i] * integral_gauss;
    }
    return ftrap_out;
}

// --- Zcharge and Zeff Calculation ---
void calculate_Zcharge(const double* ne_fa, const double* ni_fa, int nr, double Zinput, 
                       double& Z_e, double& Z_i, double& Zeff_out) 
{
    double Zs = 1.0;
    double Z_alpha = 1.0;
    double ne_tot = 0.0;
    double ni_tot = 0.0;

    for (int i = 0; i < nr; ++i) {
        ne_tot += ne_fa[i];
        ni_tot += ni_fa[i];
    }

    double Z_bar = ne_tot / ni_tot;
    Zeff_out = (ni_tot * Zs * Zs) / ne_tot;

    Z_e = Zinput;
    Z_i = std::pow((std::pow(Z_alpha, 2) * Z_bar * Zeff_out), 0.25);
}


CollisionResults calculate_collision_frequencies(
    int nr, int nphi, int ntheta, 
    const double* ne_fa, const double* ni_fa, 
    const double* Te_fa, const double* Ti_fa,
    double Ze, double Zi, const double* Bmax, const double* Bmin,
    const double* q, const double* GplusiI, const double* Binv_fa,
    const double* R, 
    int UseBfore, int Nzp) 
{
    CollisionResults res;
    res.nu_e_star.resize(nr); 
    res.nu_i_star.resize(nr);
    res.epsilon.resize(nr);  
    res.qR.resize(nr);

    for (int i = 0; i < nr; ++i) {
        if (UseBfore == 1) {
            // Magnetic-based epsilon and qR
            res.epsilon[i] = (Bmax[i] - Bmin[i]) / (Bmax[i] + Bmin[i]);
            res.qR[i] = GplusiI[i] * Binv_fa[i] / (1.0 / q[i] - Nzp);
        } else {
            // Geometric-based epsilon and qR (R_plus/R_minus logic)
            double r_max = -1e30;
            double r_min = 1e30;

            for (int j = 0; j < nphi; ++j) {
                for (int k = 0; k < ntheta; ++k) {
                    int idx = i + j * nr + k * nr * nphi;
                    if (R[idx] > r_max) r_max = R[idx];
                    if (R[idx] < r_min) r_min = R[idx];
                }
            }

            double Rmin_geom = 0.5 * (r_max - r_min);
            double Rc_geom   = 0.5 * (r_max + r_min);
            
            res.epsilon[i] = Rmin_geom / Rc_geom;
            res.qR[i] = q[i] * Rc_geom;
        }

        // ln_lambda calculations
        double ln_lam_e = 31.3 - std::log(std::sqrt(ne_fa[i]) / Te_fa[i]);
        double ln_lam_i = 30.0 - std::log((std::pow(Zi, 3) * std::sqrt(ni_fa[i])) / std::pow(Ti_fa[i], 1.5));

     
        // Formula: 6.921e-18 * Ze * ne * ln_lam / (Te^2 * eps^1.5) * qR
        res.nu_e_star[i] = std::abs(6.921e-18 * Ze * ne_fa[i] * ln_lam_e / 
                           (std::pow(Te_fa[i], 2) * std::pow(res.epsilon[i], 1.5)) * res.qR[i]);

  
        // Formula: 4.9e-18 * ni * Zi^4 * ln_lam / (Ti^2 * eps^1.5) * qR
        res.nu_i_star[i] = std::abs(4.9e-18 * ni_fa[i] * std::pow(Zi, 4) * ln_lam_i / 
                           (std::pow(Ti_fa[i], 2) * std::pow(res.epsilon[i], 1.5)) * res.qR[i]);
    }
    return res;
}

RedlResults calculate_redl_terms(int nr, double Zeff_phys, 
                                 const std::vector<double>& ftrap, 
                                 const std::vector<double>& nue, 
                                 const std::vector<double>& nui) 
{
    RedlResults res;
    res.L31.resize(nr); res.L32.resize(nr); res.L34.resize(nr); res.alpha_nu_i.resize(nr);

    // Guard Zc to prevent Zc^1.2 - 0.71 from hitting zero
    
    double Zc = std::max(Zeff_phys, 1.0); 
    double Z_denom = std::pow(Zc, 1.2) - 0.71;
    if (std::abs(Z_denom) < 1e-5) Z_denom = 1e-5; // Prevent singularity

    for (int i = 0; i < nr; ++i) {
        double f = ftrap[i];
        double ne = nue[i];
        double ni = nui[i];

        //  f_t31
        double ft31 = f / (1.0 + 0.67 * (1.0 - 0.7 * f) * std::sqrt(ne) / (0.56 + 0.44 * Zc) +
                      (0.52 + 0.086 * std::sqrt(ne)) * (1.0 + 0.87 * f) * ne / (1.0 + 1.13 * std::sqrt(std::max(Zc - 1.0, 0.0))));

        //  L31 and L34 
        double Zfac = 1.0 / Z_denom;
        res.L31[i] = (1.0 + 0.15 * Zfac) * ft31 - 0.22 * Zfac * std::pow(ft31, 2) + 
                      0.01 * Zfac * std::pow(ft31, 3) + 0.06 * Zfac * std::pow(ft31, 4);
        res.L34[i] = res.L31[i];

        //  f_t32ee and f_t32ei
        double ft32ee = f / (1.0 + 0.23 * (1.0 - 0.96 * f) * std::sqrt(ne) / std::sqrt(Zc) + 
                        0.13 * (1.0 - 0.38 * f) * ne / (Zc * Zc) * (std::sqrt(1.0 + 2.0 * std::sqrt(std::max(Zc - 1.0, 0.0))) + 
                         std::pow(f, 2) * std::sqrt((0.075 + 0.25 * std::pow(Zc - 1.0, 2)) * ne)));

        double ft32ei = f / (1.0 + 0.87 * (1.0 + 0.39 * f) * std::sqrt(ne) / (1.0 + 2.95 * std::pow(Zc - 1.0, 2)) + 
                        1.53 * (1.0 - 0.37 * f) * ne * (2.0 + 0.375 * (Zc - 1.0)));

        //  F32 terms
        double t_ee = ft32ee;
        double F32ee = (0.1 + 0.6 * Zc) / (Zc * (0.77 + 0.63 * (1.0 + std::pow(Zc - 1.0, 1.1)))) * (t_ee - std::pow(t_ee, 4)) +
                       0.7 / (1.0 + 0.2 * Zc) * (std::pow(t_ee, 2) - std::pow(t_ee, 4) - 1.2 * (std::pow(t_ee, 3) - std::pow(t_ee, 4))) +
                       1.3 / (1.0 + 0.5 * Zc) * std::pow(t_ee, 4);

        double t_ei = ft32ei;
        double F32ei = -(0.4 + 1.93 * Zc) / (Zc * (0.8 + 0.6 * Zc)) * (t_ei - std::pow(t_ei, 4)) +
                       5.5 / (1.5 + 2.0 * Zc) * (std::pow(t_ei, 2) - std::pow(t_ei, 4) - 0.8 * (std::pow(t_ei, 3) - std::pow(t_ei, 4))) -
                       1.3 / (1.0 + 0.5 * Zc) * std::pow(t_ei, 4);

        res.L32[i] = F32ee + F32ei;

        // alpha_nu_i
        double a0 = -(0.62 + 0.055 * (Zc - 1.0)) / (0.53 + 0.17 * (Zc - 1.0)) * (1.0 - f) / (1.0 - (0.31 - 0.065 * (Zc - 1.0)) * f - 0.25 * f * f);
        res.alpha_nu_i[i] = ((a0 + 0.7 * Zc * std::sqrt(f) * std::sqrt(ni)) / (1.0 + 0.18 * std::sqrt(ni)) - 0.002 * ni * ni * std::pow(f, 6)) / (1.0 + 0.004 * ni * ni * std::pow(f, 6));
        
        // Final check: if any result is NaN due to ne/ni at axis, force to 0
        if (std::isnan(res.L32[i])) res.L32[i] = 0.0;
        if (std::isnan(res.L31[i])) res.L31[i] = 0.0;
    }
    return res;
}

SauterResults calculate_sauter_terms(int nr, double Zeff, 
                                     const std::vector<double>& ftrap, 
                                     const std::vector<double>& nue, 
                                     const std::vector<double>& nui) 
{

    SauterResults res;
    res.L31.resize(nr); res.L32.resize(nr); res.L34.resize(nr); res.alpha_nu_i.resize(nr);

    for (int i = 0; i < nr; ++i) {
        double f = ftrap[i];
        double ne = nue[i];
        double ni = nui[i];
        double Zc = std::max(Zeff, 1e-12); // Prevent division by zero

        // Effective trapped fractions (Sauter specific fits)
        double ft31 = f / (1.0 + (1.0 - 0.1 * f) * std::sqrt(ne) + 0.5 * (1.0 - f) * ne / Zc);
        double ft32ee = f / (1.0 + 0.26 * (1.0 - f) * std::sqrt(ne) + 0.18 * (1.0 - 0.37 * f) * ne / std::sqrt(Zc));
        double ft32ei = f / (1.0 + (1.0 + 0.6 * f) * std::sqrt(ne) + 0.85 * (1.0 - 0.37 * f) * ne * (1.0 + Zc));
        double ft34 = f / (1.0 + (1.0 - 0.1 * f) * std::sqrt(ne) + 0.5 * (1.0 - 0.5 * f) * ne / Zc);

        // Alpha (Ion collisionality coefficient)
        double alpha_0 = -1.17 * (1.0 - f) / (1.0 - 0.22 * f - 0.19 * f * f);
        res.alpha_nu_i[i] = ((alpha_0 + 0.25 * (1.0 - f * f) * std::sqrt(ni)) / (1.0 + 0.5 * std::sqrt(ni)) 
                            + 0.315 * std::pow(ni, 2) * std::pow(f, 6)) 
                            / (1.0 + 0.15 * std::pow(ni, 2) * std::pow(f, 6));

        // L31
        res.L31[i] = (1.0 + 1.4 / (Zc + 1.0)) * ft31 
                   - 1.9 / (Zc + 1.0) * std::pow(ft31, 2) 
                   + 0.3 / (Zc + 1.0) * std::pow(ft31, 3) 
                   + 0.2 / (Zc + 1.0) * std::pow(ft31, 4);

        // L34
        res.L34[i] = (1.0 + 1.4 / (Zc + 1.0)) * ft34 
                   - 1.9 / (Zc + 1.0) * std::pow(ft34, 2) 
                   + 0.3 / (Zc + 1.0) * std::pow(ft34, 3) 
                   + 0.2 / (Zc + 1.0) * std::pow(ft34, 4);

        // F32ee and F32ei
        double F32ee = (0.05 + 0.62 * Zc) / (Zc * (1.0 + 0.44 * Zc)) * (ft32ee - std::pow(ft32ee, 4)) 
                     + 1.0 / (1.0 + 0.22 * Zc) * (std::pow(ft32ee, 2) - std::pow(ft32ee, 4) - 1.2 * (std::pow(ft32ee, 3) - std::pow(ft32ee, 4))) 
                     + 1.2 / (1.0 + 0.5 * Zc) * std::pow(ft32ee, 4);

        double F32ei = -(0.56 + 1.93 * Zc) / (Zc * (1.0 + 0.44 * Zc)) * (ft32ei - std::pow(ft32ei, 4)) 
                     + 4.95 / (1.0 + 2.48 * Zc) * (std::pow(ft32ei, 2) - std::pow(ft32ei, 4) - 0.55 * (std::pow(ft32ei, 3) - std::pow(ft32ei, 4))) 
                     - 1.2 / (1.0 + 0.5 * Zc) * std::pow(ft32ei, 4);

        res.L32[i] = F32ee + F32ei;
    }
    return res;
}