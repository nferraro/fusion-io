#include "fusion_io.h"
#include <iostream>

int mars_field::load(const fio_option_list* opt)
{
  opt->get_option(FIO_TIMESLICE, &time);
  opt->get_option(FIO_LINEAR_SCALE, &linfac);
  opt->get_option(FIO_PART, &part);
  opt->get_option(FIO_PHASE, &phase);

  if(linfac != 1. && part==FIO_EQUILIBRIUM_ONLY) {
    std::cerr << "Linear scale linfac is ignored for equilibrium data." 
	      << std::endl;
    linfac = 1.;
  }
  

  return FIO_SUCCESS;
}


int mars_magnetic_field::load(const fio_option_list* opt)
{
  mars_field::load(opt);

  int result;

  result = source->load_eigdata("BPLASMA.OUT",&(source->bplasma),false);
  if(result != FIO_SUCCESS) return result;

  return FIO_SUCCESS;
}


int mars_magnetic_field::eval(const double* x, double* v)
{  
  return FIO_SUCCESS;
}


int mars_fluid_velocity::load(const fio_option_list* opt)
{
  mars_field::load(opt);

  return FIO_SUCCESS;
}

int mars_fluid_velocity::eval(const double* x, double* v)
{  
  return FIO_SUCCESS;
}
