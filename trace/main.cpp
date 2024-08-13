#include "trace_integrator.h"

#include <fusion_io_source.h>
#include <m3dc1_source.h>
#include <geqdsk_source.h>
#include <gpec_source.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <stdlib.h>

trace_integrator tracer;
int transits = 100;
int steps_per_transit = 100;
int surfaces = 11;
int tpts = 1;
double dR = 0.;
double dZ = 0.;
double dR0 = 0.;
double dZ0 = 0.;
double Phi = 0.;
double Phi0 = 0.;
double R0 = 0.;
double Z0 = 0.;
double angle = 0.;
bool qout = true;
bool pout = true;
bool reverse = false;
double ds = 0.;
int nplanes = 1;

bool R0_set = false;
bool Z0_set = false;

void print_help();
void print_parameters();
bool process_command_line(int argc, char* argv[]);
bool process_line(const std::string& opt, const int argc, 
		  const std::string argv[]);
void delete_sources();

void partition();
void gather(const double* in, double* out);

int size, rank;
int *offset, *local_surfaces;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!process_command_line(argc, argv)) {
    print_help();
    return 1;
  }

  tracer.tpts = tpts;
  tracer.plane = angle;
  tracer.nplanes = nplanes;
  tracer.set_reverse(reverse);

  if(tracer.sources.size() == 0) {
    std::cerr << "No sources specified.  Returning." << std::endl;
    print_help();
    return 1;
  }

  int t;
  if(tracer.sources[0].source->get_int_parameter(FIO_GEOMETRY,&t) == FIO_SUCCESS) {
    tracer.toroidal = (t==FIO_CYLINDRICAL);
  } else {
    std::cerr << "Warning: can't determine geometry.  Assuming toroidal."
	      << std::endl;
    tracer.toroidal = true;
  }
  if(tracer.toroidal) {
    std::cerr << "Geometry: toroidal" << std::endl;
  } else {
    std::cerr << "Geometry: straight cylinder" << std::endl;
  }


  double p;
  if(tracer.sources[0].source->get_real_parameter(FIO_PERIOD,&p) == FIO_SUCCESS) {
    tracer.period = p;
  } else {
    std::cerr << "Warning: can't determine period.  Assuming 2*pi"
	      << std::endl;
    tracer.period = 2.*M_PI;
  }
  std::cerr << "Period: " << tracer.period << std::endl;


  double R_axis, Z_axis;
  tracer.center(&R_axis, &Z_axis);

  std::cerr << "Magnetic axis is at ( " << R_axis << ", "
	    << Z_axis << ")" << std::endl;

  if(!R0_set) R0 = R_axis;
  if(!Z0_set) Z0 = Z_axis;

  if(dR0==0.) dR0 = dR;

  print_parameters();

  std::fstream plotfile;
  char outfile[256];

  if(pout) {
    plotfile.open("gplot", std::fstream::out | std::fstream::trunc);
    plotfile << "set size ratio -1" << '\n'
	     << "set xlabel 'R'" << '\n'
	     << "set ylabel 'Z'" << '\n'
	     << "set pointsize 0.1" << '\n'
	     << "set mapping cartesian" << '\n'
	     << "plot ";

    for(int j=0; j<surfaces; j++) {
      sprintf(outfile, "out%d", j);     
      if(j>0) plotfile << ",\\\n";
      plotfile << "'" << outfile << "'" << " using 2:3 title ''";
    }
    plotfile << std::endl;
    plotfile.close();

    plotfile.open("gplotpsi", std::fstream::out | std::fstream::trunc);
    plotfile << "set size ratio 1" << '\n'
	      << "set xlabel 'theta'" << '\n'
	      << "set ylabel 'psi'" << '\n'
	      << "set pointsize 0.1" << '\n'
	      << "set mapping cartesian" << '\n'
	      << "plot [-180:180] ";

    for(int j=0; j<surfaces; j++) {
      sprintf(outfile, "out%d", j);     
      if(j>0) plotfile << ",\\\n";
      plotfile << "'" << outfile << "'" << " using 4:5 title ''";
    }
    plotfile << std::endl;
    plotfile.close();

    plotfile.open("gplot3D", std::fstream::out | std::fstream::trunc);
    plotfile << "set size ratio 1" << '\n'
	      << "set xlabel 'X'" << '\n'
	      << "set ylabel 'Y'" << '\n'
	      << "set zlabel 'Z'" << '\n'
	      << "set pointsize 0.1" << '\n'
	      << "set mapping cylindrical" << '\n'
	      << "set angles degrees" << '\n'
	      << "splot ";

    for(int j=0; j<surfaces; j++) {
      sprintf(outfile, "out%d", j);     
      if(j>0) plotfile << ",\\\n";
      if(nplanes > 12) {
	plotfile << "'" << outfile << "'" << " u 1:3:2 w lines title ''";
      } else {
	plotfile << "'" << outfile << "'" << " u 1:3:2 title ''";
      }
    }
    plotfile << std::endl;
    plotfile.close();
  }

  offset = new int[size];
  local_surfaces = new int[size];

  partition();
  
  double *my_q_min, *my_q_max, *my_q_mean, *my_r, *my_phase;
  my_r = new double[local_surfaces[rank]];
  my_q_min = new double[local_surfaces[rank]];
  my_q_max = new double[local_surfaces[rank]];
  my_q_mean = new double[local_surfaces[rank]];
  my_phase = new double[local_surfaces[rank]];

  time_t t1 = time(0);
  for(int j=0; j<local_surfaces[rank]; j++) {
    int s = rank + j*size;
    double R = R0 + dR0 + dR*s;
    double Z = Z0 + dZ0 + dZ*s;
    double *rr, *zz;
    bool result;

    int n;

    if(pout)  {
      sprintf(outfile, "out%d", s);
      tracer.open_file(outfile);
    }

    if(ds==0.) {
      n = 1;
      rr = new double[1];
      zz = new double[1];
      rr[0] = R;
      zz[0] = Z;
    } else { 
      if(!tracer.get_surface(R, Phi, Z, ds, &rr, &zz, &n))
	std::cerr << "Warning: not on closed surface" << std::endl;

      if(n==0) {
	std::cerr << "Error: no points found on surface" << std::endl;
	break;
      }
    }

    my_r[j] = sqrt((R - R_axis)*(R - R_axis) + (Z - Z_axis)*(Z - Z_axis));
    my_phase[j] = Phi;

    std::cerr << "Surface " << s+1 << " of " << surfaces
	      << " (" << n << " pts at this surface)" << std::endl;
    std::cerr << "(" << R << "," << Phi <<  ", "  << Z << ") ... ";

    // Cycle through surface points
    trace_integrator::integrator_data data;
    my_q_mean[j] = 0.;
    my_q_max[j] = 0.;
    my_q_min[j] = 0.;
    int tt = 0;
    result = false;
    for(int k=0; k<n; k++) {   
      Phi = Phi0;

      std::cerr << atan2(zz[k] - Z_axis, rr[k] - R_axis) << "\t";
      //      dfile << rr[k] << '\t' <<  zz[k] << '\t';
      tracer.set_pos(rr[k],Phi,zz[k]);

      // perform integration
      result = tracer.integrate(transits, steps_per_transit, &data);
      if(!result) {
	double RR, PP, ZZ;
	tracer.get_pos(&RR, &PP, &ZZ);
	std::cerr << "lost at (" << RR << ", " << ZZ <<  ")" << std::endl;
      }

      // write connection length
      //      if(result)
      //	dfile << 0. << std::endl;
      //      else
      //      dfile << data.distance << '\t' << log10(data.distance) << std::endl;

      if(k==0 || data.q < my_q_min[j]) my_q_min[j] = data.q;
      if(k==0 || data.q > my_q_max[j]) my_q_max[j] = data.q;
				   
      my_q_mean[j] += data.q*data.toroidal_transits;
      tt += data.toroidal_transits;
    }
    my_q_mean[j] /= (double)tt;
    
    std::cerr << tt << " tor. transits, q = " << my_q_mean[j] << std::endl;

    tracer.close_file();

    delete[] rr;
    delete[] zz;
  }
  time_t t2 = time(0);
  std::cerr << "Computation completed in " << t2 - t1 << " s." << std::endl;
  
  // assemble and write data
  // ~~~~~~~~~~~~~~~~~~~~~~~

  // allocate global data
  double *total_r, *total_q_min, *total_q_max, *total_q_mean, *total_phase;

  if(rank==0) {
    total_r = new double[surfaces];

    if(qout) {
      total_q_min = new double[surfaces];
      total_q_max = new double[surfaces];
      total_q_mean = new double[surfaces];
      total_phase = new double[surfaces];
    }
  }

  gather(my_r, total_r);
  if(qout) {
    gather(my_q_min, total_q_min);
    gather(my_q_max, total_q_max);
    gather(my_q_mean, total_q_mean);
    gather(my_phase, total_phase);
  }
  /*
  MPI_Gatherv(my_r, local_surfaces[rank], MPI_DOUBLE, 
	      total_r, local_surfaces, offset, MPI_DOUBLE, 
	      0, MPI_COMM_WORLD);
  if(qout) {
    MPI_Gatherv(my_q_min, local_surfaces[rank], MPI_DOUBLE, 
		total_q_min, local_surfaces, offset, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
    MPI_Gatherv(my_q_max, local_surfaces[rank], MPI_DOUBLE, 
		total_q_max, local_surfaces, offset, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
    MPI_Gatherv(my_q_mean, local_surfaces[rank], MPI_DOUBLE, 
		total_q_mean, local_surfaces, offset, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
    MPI_Gatherv(my_phase, local_surfaces[rank], MPI_DOUBLE, 
		total_phase, local_surfaces, offset, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
  }
  */

  if(rank==0) {
    // write data
    if(qout) {
      std::fstream qfile;     
      qfile.open("q.out", std::fstream::out | std::fstream::trunc);
      for(int i=0; i<surfaces; i++) {
	qfile << total_r[i] << '\t' 
	      << total_q_mean[i] << '\t'
	      << total_q_min[i] << '\t'
	      << total_q_max[i] << '\n';
      }
      qfile.close();
    }
    
    // delete global data
    delete[] total_r;
    if(qout) {
      delete[] total_q_min;
      delete[] total_q_max;
      delete[] total_q_mean;
      delete[] total_phase;
    }
  }

  // delete local data
  delete[] my_r;
  delete[] my_q_min;
  delete[] my_q_max;
  delete[] my_q_mean;
  delete[] my_phase;

  delete[] offset;
  delete[] local_surfaces;

  std::cerr << "Deleting sources." << std::endl;
  delete_sources();

  std::cerr << "===============================\n" 
	    << "Poincare computation complete.\n"
	    << "To view Poincare plot, use 'gnuplot gplot'" << std::endl;
  MPI_Finalize();

  return 0;
}

void delete_sources()
{
  trace_source_list::iterator i;
  i = tracer.sources.begin();

  while(i != tracer.sources.end()) {
    i->free();
    i++;
  }
}

bool process_command_line(int argc, char* argv[])
{
  const int max_args = 4;
  const int num_opts = 21;
  std::string arg_list[num_opts] = 
    { "-gpec", "-geqdsk", "-m3dc1", "-diiid-i",
      "-dR", "-dZ", "-dR0", "-dZ0","-R0","-Z0",
      "-ds", "-p", "-t", "-s", "-a",
      "-pout", "-qout", "-phi0", "-n", 
      "-reverse", "-tavg" };
  std::string opt = "";
  std::string arg[max_args];
  int args = 0;
  bool is_opt;
  bool processed = true;

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
	if(!process_line(opt, args, arg)) return false;

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
    if(!process_line(opt, args, arg)) return false;

  return true;

}

void partition()
{
  for(int i=0; i<size; i++)
    local_surfaces[i] = 0.;

  int r = 0;
  for(int i=0; i<surfaces; i++) {
    local_surfaces[r]++;
    r = (r+1) % size;
  }

  offset[0] = 0.;
  for(int i=1; i<size; i++)
    offset[i] = offset[i-1]+local_surfaces[i-1];
}

void gather(const double* in, double* out)
{
  double* buf = new double[surfaces];

  MPI_Gatherv((void*)in, local_surfaces[rank], MPI_DOUBLE, 
	      buf, local_surfaces, offset, MPI_DOUBLE, 
	      0, MPI_COMM_WORLD); 

  if(rank == 0) {
    int n=0;
    int r=0;
    for(int i=0; i<surfaces; i++) {
      out[i] = buf[offset[r]+n];
      r++;
      if(r==size) {
	n++;
	r = 0;
      }
    }
  }

  delete[] buf;
}

bool create_source(const int type, const int argc, const std::string argv[]) 
{
  trace_source src;
  fio_option_list fopt;
  fio_series* magaxis[2];
  int result;
  double slice_time = 0.;

  switch(type) {
  case(FIO_M3DC1_SOURCE):
    src.source = new m3dc1_source();
    if(argc>=1) {
      result = src.source->open(argv[0].c_str());
    } else {
      result = src.source->open("C1.h5");
    }
    if(result != FIO_SUCCESS) {
      std::cerr << "Error opening file" << std::endl;
      src.free();
      return false;
    };
    // set options for fields obtained from this source
    src.source->get_field_options(&fopt);
    fopt.set_option(FIO_PART, FIO_PERTURBED_ONLY);
    if(argc<2) fopt.set_option(FIO_TIMESLICE,-1);  // default to timeslice=-1
    break;

  case(FIO_GEQDSK_SOURCE):
    src.source = new geqdsk_source();
    if(argc>=1) {
      if((result = src.source->open(argv[0].c_str())) != FIO_SUCCESS) {
	std::cerr << "Error opening file" << std::endl;
	src.free();
	return false;
      };
    } else {
      std::cerr << "Filename must be provided for geqdsk files." << std::endl;
      src.free();
      return false;
    }
    src.source->get_field_options(&fopt);
    break;

  case(FIO_GPEC_SOURCE):
    src.source = new gpec_source();
    if(argc>=1) {
      if((result = src.source->open(argv[0].c_str())) != FIO_SUCCESS) {
	std::cerr << "Error opening file" << std::endl;
	src.free();
	return false;
      };
    } else {
      std::cerr << "Directory must be provided for gpec files." << std::endl;
      src.free();
      return false;
    }
    src.source->get_field_options(&fopt);
    break;

  default:
    return false;
  }

  if(argc>=2) {
    int timeslice = atoi(argv[1].c_str());
    std::cerr << "Time slice " << timeslice << std::endl;
    fopt.set_option(FIO_TIMESLICE, timeslice);
    if(src.source->get_slice_time(timeslice, &slice_time)==FIO_SUCCESS) {
      std::cerr << "Slice time = " << slice_time << std::endl;
    }
  }
  if(argc>=3) fopt.set_option(FIO_LINEAR_SCALE, atof(argv[2].c_str()));
  if(argc>=4) fopt.set_option(FIO_PHASE, atof(argv[3].c_str())*M_PI/180.);

  // Read magnetic field
  result = src.source->get_field(FIO_MAGNETIC_FIELD, &src.field, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error reading magnetic field" << std::endl;
    src.free();
    return false;
  };

  // Read Psi_Norm field
  fopt.set_option(FIO_PART, FIO_EQUILIBRIUM_ONLY);
  result = src.source->get_field(FIO_POLOIDAL_FLUX_NORM, &src.psi_norm, &fopt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Warning: couldn't open psi_norm field" << std::endl;
  };
  
  // Allocate search hint
  result = src.source->allocate_search_hint(&src.hint);
  if(result != FIO_SUCCESS) {
    std::cerr << "Warning: couldn't allocate search hint" << std::endl;
  }

  // Find magnetic axis
  if((result = src.source->get_series(FIO_MAGAXIS_R, &magaxis[0])) != FIO_SUCCESS) {
    std::cerr << "Couldn't load MAGAXIS_R" << std::endl;
    src.free();
    return false;
  }
  if((result = src.source->get_series(FIO_MAGAXIS_Z, &magaxis[1])) != FIO_SUCCESS) {
    std::cerr << "Couldn't load MAGAXIS_Z" << std::endl;
    src.free();
    return false;
  }
  if((result = magaxis[0]->eval(slice_time, &src.magaxis[0])) != FIO_SUCCESS)
    {
      std::cerr << "Couldn't evaluate MAGAXIS_R" << std::endl;
      src.free();
      return false;
    }
  if((result = magaxis[1]->eval(slice_time, &src.magaxis[1])) != FIO_SUCCESS)
    {
      std::cerr << "Couldn't evaluate MAGAXIS_Z" << std::endl;
      src.free();
      return false;
    }
  delete(magaxis[0]);
  delete(magaxis[1]);


  // Add source to list
  tracer.sources.push_back(src);
    
  return true;
}

bool process_line(const std::string& opt, const int argc, const std::string argv[])
{
  bool argc_err = false;

  if(opt=="-m3dc1") {
    return create_source(FIO_M3DC1_SOURCE, argc, argv);
  } else if(opt=="-geqdsk") {
    return create_source(FIO_GEQDSK_SOURCE, argc, argv);
  } else if(opt=="-gpec") {
    return create_source(FIO_GPEC_SOURCE, argc, argv);
  } else if(opt=="-diiid-i") {
    /*
    coil_source* s = new coil_source();
    double current;
    int n;
    double phase;
    if(argc>=1) current = atof(argv[0].c_str()); else current = 1e3;
    if(argc>=2) n = atoi(argv[1].c_str());       else n = 1;
    if(argc>=3) phase = atof(argv[2].c_str());   else phase = 0.;
    diiid_icoils(s, current, n, 0., phase);
    tracer.sources.push_back(s);
    */
  } else if(opt=="-dR") {
    if(argc==1) dR = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dZ") {
    if(argc==1) dZ = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dR0") {
    if(argc==1) dR0 = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dZ0") {
    if(argc==1) dZ0 = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-R0") {
    if(argc==1) {
      R0 = atof(argv[0].c_str());
      R0_set = true;
    } else argc_err = true;
  } else if(opt=="-Z0") {
    if(argc==1) {
      Z0 = atof(argv[0].c_str());
      Z0_set = true;
    }
    else argc_err = true;
  } else if(opt=="-ds") {
    if(argc==1) ds = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-p") {
    if(argc==1) surfaces = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-t") {
    if(argc==1) transits = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-tavg") {
    if(argc==1) tpts = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-s") {
    if(argc==1) steps_per_transit = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-a") {
    if(argc==1) angle = atof(argv[0].c_str())*M_PI/180.;
    else argc_err = true;
  } else if(opt=="-phi0") {
    if(argc==1) {
      Phi0 = atof(argv[0].c_str())*M_PI/180.;
    } else argc_err = true;
  } else if(opt=="-n") {
    if(argc==1) nplanes = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-pout") {
    if(argc==1) pout = (atoi(argv[0].c_str()) == 1);
    else argc_err = true;
  } else if(opt=="-qout") {
    if(argc==1) qout = (atoi(argv[0].c_str()) == 1);
    else argc_err = true;
  } else if(opt=="-reverse") {
    if(argc==0) reverse = true;
    else argc_err = true;
  } else {
    std::cerr << "Unrecognized option " << opt << std::endl;
    return false;
  }

  if(argc_err) {
    std::cerr << "Incorrect number of arguments for option " 
	      << opt << std::endl;
    return false;
  }

  return true;
}


void print_help()
{
  std::cout 
    << "Usage:\n"
    << "\ttrace <sources> -dR <dR> -dZ <dZ> -dR0 <dR0> -dZ0 <dZ0> /\n"
    << "\t   -p <pts> -t <trans> -s <steps> -a <angle> -phi0 <phi0> \n"
    << "\t   -pplot <pplot> -qplot <qplot> -n <nplanes> -tavg <tavg> \n\n"
    << " <angle>   The toroidal angle of the plane in degrees (default = 0)\n"
    << " <dR>      R-spacing of seed points (default = major radius/(2*pts))\n"
    << " <dR0>     R-distance of first seed point from axis (default = dR)\n"
    << " <dZ>      Z-spacing of seed points (default = 0)\n"
    << " <dZ0>     Z-distance of first seed point from axis (default = dZ)\n"
    << " <phi0>    Toroidal angle of all seed points in degrees\n"
    << " <pts>     Number of seed points (default = 11)\n"
    << " <steps>   Integration steps per toroidal transit (default = 100)\n"
    << " <pplot>   If 1, Generate poincare plot data (default = 1)\n"
    << " <qplot>   If 1, Generate q-profile data (default = 1)\n"
    << " <trans>   Toroidal transits per seed point (default = 100)\n"
    << " <tavg>    Number of toroidal planes for toroidal averaging (default = 1)\n"
    << " <nplanes> Number of toroidal planes for output (default = 1)\n"
    << " <sources> May be one or more of the following:\n"
    << "\n  -m3dc1 <c1_file> <ts> <factor> <shift>\n"
    << "   * Loads field information from M3D-C1 data file\n"
    << "   <c1_file> filename of M3D-C1 hdf5 file (default = C1.h5)\n" 
    << "   <ts>      timeslice of M3D-C1 hdf5 file (default = -1)\n"
    << "   <factor>  factor by which to multiply field (default = 1)\n"
    << "   <shift>   toroidal phase shift of field, in degrees (default = 0)\n"
    << "\n  -geqdsk <gfile>\n"
    << "   * Loads field information from EFIT g-file\n"
    << "   <gfile>   filename of EFIT g-file\n"
    << "\n  -diiid-i <curr> <n> <phase>\n"
    << "   * Calculates field from DIII-D I-coil set\n"
    << "   <curr>    current (in Amps)\n"
    << "   <n>       toroidal modenumber\n"
    << "   <phase>   phase difference between top and bottom coil sets (deg)\n"
    << "\nExample:\n"
    << "\ttrace -m3dc1 C1.h5 -m3dc1 C1.h5 1 1e4\n\n"
    << " Loads field from time slices -1 and 1 from M3D-C1 file 'C1.h5'.\n"
    << " The field from time slice 1 is multiplied by 1e4.\n"
    << std::endl;
}

void print_parameters()
{
  std::cerr << "dR0 = " << dR0 << '\n'
	    << "dR = " << dR << '\n'
	    << "dZ0 = " << dZ0 << '\n'
	    << "dZ = " << dZ << '\n'
	    << "Number of surfaces = " << surfaces << '\n'
	    << "Toroidal transits per seed = " << transits << '\n'
	    << "Integrator steps per transit = " << steps_per_transit 
	    << std::endl;
}
