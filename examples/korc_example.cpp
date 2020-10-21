#include <m3dc1_source.h>
#include <m3dc1_field.h>
#include <fusion_io.h>
#include <options.h>
#include <compound_field.h>

#include <iostream>
#include <random>
#include <cmath>
#include <array>
#include <vector>
#include <chrono>

void write_time(const std::string &name, const std::chrono::nanoseconds time);

int main() {
    int result;
    fio_source *src;
    fio_field *magnetic_field;
    fio_option_list opt;

// Open an m3dc1 source
//    result = fio_open_source(&src, FIO_M3DC1_SOURCE, "C1.h5");
    result = fio_open_source(&src, FIO_M3DC1_SOURCE, "/Users/nferraro/data/Teng/it27_nl/C1.h5");
    if(result != FIO_SUCCESS) {
        std::cerr << "Error opening file" << std::endl;
        delete(src);
        return result;
    };

// set options for fields obtained from this source
    src->get_field_options(&opt);
    opt.set_option(FIO_TIMESLICE, 0);

// open fields
    result = src->get_field(FIO_MAGNETIC_FIELD, &magnetic_field, &opt);
    if(result != FIO_SUCCESS) {
        std::cerr << "Error opening magnetic field" << std::endl;
        delete(src);
        return result;
    };
  
    size_t npts = 1000;
    size_t nsteps = 1;
    std::vector<std::array<double, 3> > x_array(npts);
    double b[3];

    size_t hint_size = src->sizeof_search_hint();
    std::cout << "Size of search hint: " << hint_size << std::endl;
    void* hint = calloc(npts, hint_size);
  
    std::mt19937_64 engine;
    std::uniform_real_distribution<double> r_dist(-0.1, 0.1);
    std::uniform_real_distribution<double> z_dist(-0.1, 0.1);
    std::uniform_real_distribution<double> phi_dist(0, 2.0*M_PI);

    std::cout << "Starting Random" << std::endl;
    for(std::array<double, 3> &x : x_array) {
        x[0] = r_dist(engine);
        x[1] = phi_dist(engine);
        x[2] = z_dist(engine);
    }
    std::cout << "Random Finished" << std::endl;
    
    std::cout << "Starting Interpolation" << std::endl;
    const std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    for(int i=0; i<nsteps; i++) {
      for(int p=0; p<npts; p++) {
	void* h = ((char*)hint) + hint_size*p;
	int h0 = *((int*)h);
        result = magnetic_field->eval(x_array[p].data(), b, h);
	int h1 = *((int*)h);
	if(i > 0 && h0 != h1) 
	  std::cerr << "Bad guess. " << h0 << " " << h1 << std::endl;
      }
    }
    const std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Interpolation Finished" << std::endl;

    const auto total_time = end_time - start_time;
    const std::chrono::nanoseconds total_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds> (total_time);
    std::cout << "Interpolation Finished" << std::endl;
    
    std::cout << std::endl << "Timing:" << std::endl;
    write_time("  Interpolation time : ", total_time_ns);
    std::cout << std::endl;
    
    fio_close_field(&magnetic_field);
    fio_close_source(&src);
    free(hint);

    return 0;
}

void write_time(const std::string &name, const std::chrono::nanoseconds time) {
    if (time.count() < 1000) {
        std::cout << name << time.count()               << " ns" << std::endl;
    } else if (time.count() < 1000000) {
        std::cout << name << time.count()/1000.0        << " μs" << std::endl;
    } else if (time.count() < 1000000000) {
        std::cout << name << time.count()/1000000.0     << " ms" << std::endl;
    } else if (time.count() < 60000000000) {
        std::cout << name << time.count()/1000000000.0  << " s" << std::endl;
    } else if (time.count() < 3600000000000) {
        std::cout << name << time.count()/60000000000.0 << " min" << std::endl;
    } else {
        std::cout << name << time.count()/3600000000000 << " h" << std::endl;
    }
}

