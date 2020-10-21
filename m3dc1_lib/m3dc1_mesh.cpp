#include "m3dc1_mesh.h"

#include <iostream>
#include <math.h>

m3dc1_mesh::m3dc1_mesh(int n)
{
  a = new double[n];
  b = new double[n];
  c = new double[n];
  co = new double[n];
  sn = new double[n];
  x = new double[n];
  z = new double[n];
  bound = new int[n];
  region = new int[n];
  nelms = n;
  memory_depth = 0;
  last_elm = -1;
  next_elm = 0;
  hits = 0;
  misses = 0;
  period = 2.*M_PI;
  neighbor = 0;
  nneighbors = 0;
  nplanes = 1;
}

m3dc1_mesh::~m3dc1_mesh()
{
  /*
  int evals = hits+misses;
  std::cerr << "hits = " << hits
	    << " (" << 100.*(double)hits/(double)evals << "%)\n"
	    << "misses = " << misses
	    << " (" << 100.*(double)misses/(double)evals << "%)\n"
	    << "hits+misses = " << evals
	    << std::endl;
  */

  delete[] a;
  delete[] b;
  delete[] c;
  delete[] co;
  delete[] sn;
  delete[] x;
  delete[] z;
  delete[] bound;
  delete[] region;
  if(neighbor != 0) {
    for(int i=0; i<nelms; i++)
      delete[] neighbor[i];
    delete[] neighbor;
  }
  if(nneighbors != 0) delete[] nneighbors;
  clear_memory();
}

void m3dc1_mesh::extent(double* X0, double* X1,
			double* Phi0, double* Phi1,
			double* Z0, double* Z1) const
{
  if(nelms==0) {
    *X0 = *X1 = *Phi0 = *Phi1 = *Z0 = *Z1 = 0.;
    return;
  }

  double R[3], Z[3];

  *X0 = x[0];
  *X1 = x[0];
  *Z0 = z[0];
  *Z1 = z[0];

  for(int i=0; i<nelms; i++) {
    R[0] = x[i];
    R[1] = x[i] + (a[i]+b[i])*co[i];
    R[2] = x[i] + b[i]*co[i] - c[i]*sn[i];
    Z[0] = z[i];
    Z[1] = z[i] + (a[i]+b[i])*sn[i];
    Z[2] = z[i] + c[i]*co[i] + b[i]*sn[i];

    for(int j=0; j<3; j++) {
      if(R[j] < *X0) *X0 = R[j];
      if(R[j] > *X1) *X1 = R[j];
      if(Z[j] < *Z0) *Z0 = Z[j];
      if(Z[j] > *Z1) *Z1 = Z[j];
    }
  }

  *Phi0 = 0.;
  *Phi1 = period;
}

void m3dc1_mesh::clear_memory()
{
  if(next_elm) {
    for(int d=0; d<memory_depth; d++)
      delete[] next_elm[d];
    
    delete[] next_elm;

    next_elm = 0;
  }
  last_elm = 0;
  memory_depth = 0;
}

bool m3dc1_mesh::set_memory_depth(const int d)
{
  if(memory_depth) clear_memory();

  memory_depth = d;
  if(memory_depth==0) return true;

  next_elm = new int*[memory_depth];
  for(int d=0; d<memory_depth; d++) {
    next_elm[d] = new int[nelms];
    for(int i=0; i<nelms; i++) 
      next_elm[d][i] = -1;
  }

  return true;
}

int m3dc1_mesh::in_element_threadsafe(double X, double Phi, double Z,
				      double* xi, double* zi, double* eta,
				      int guess)
{
  // if a guess is provided, test it
  if(guess >= 0) {
    if(is_in_element(guess,X,Phi,Z,xi,zi,eta)) {
      hits++;
      return guess;
    }

    // Test neighbors
    for(int n=0; n<nneighbors[guess]; n++) {
      int m = neighbor[guess][n];
      if(is_in_element(m,X,Phi,Z,xi,zi,eta)) {
	hits++;
	return m;
      }
    }

    // Test neighbors' neighbors
    for(int n=0; n<nneighbors[guess]; n++) {
      int m = neighbor[guess][n];
      for(int nn=0; nn<nneighbors[m]; nn++) {
	int l = neighbor[m][nn];
	if(is_in_element(l,X,Phi,Z,xi,zi,eta)) {
	  hits++;
	  return l;
	}
      }
    }
  }

  misses++;

  // Now, search randomly.
  for(int e=0; e<nelms; e++) {
    if(is_in_element(e,X,Phi,Z,xi,zi,eta))
      return e;
  }

  // failed to find an elm containing coordinates.
  return -1;

}

int m3dc1_mesh::in_element_memory(double X, double Phi, double Z,
				  double* xi, double* zi, double* eta,
				  int guess)
{
  // if a guess is provided, test it
  if(guess >= 0) {
    if(is_in_element(guess,X,Phi,Z,xi,zi,eta)) {
      hits++;
      last_elm = guess;
      return guess;
    }

    // Test neighbors
    for(int n=0; n<nneighbors[guess]; n++) {
      int m = neighbor[guess][n];
      if(is_in_element(m,X,Phi,Z,xi,zi,eta)) {
	hits++;
	last_elm = m;
	return m;
      }
    }

    // Test neighbors' neighbors
    for(int n=0; n<nneighbors[guess]; n++) {
      int m = neighbor[guess][n];
      for(int nn=0; nn<nneighbors[m]; nn++) {
	int l = neighbor[m][nn];
	if(is_in_element(l,X,Phi,Z,xi,zi,eta)) {
	  hits++;
	  last_elm = l;
	  return l;
	}
      }
    }
  }

  if(last_elm >= 0) {
    // first, check last elm
    if(is_in_element(last_elm,X,Phi,Z,xi,zi,eta)) {
      hits++;
      return last_elm;
    }

    // now, check each element in memory
    for(int m=0; m<memory_depth; m++) {
      int e = next_elm[m][last_elm];
      if(e==-1) break;
      if(is_in_element(e,X,Phi,Z,xi,zi,eta)) {
	hits++;
	last_elm = e;
	return e;
      }
    }
  }

  misses++;

  // Now, search randomly.
  for(int e=0; e<nelms; e++) {
    if(is_in_element(e,X,Phi,Z,xi,zi,eta)) {     
      // append this to memory
      if(memory_depth>0 && last_elm >= 0) {
	for(int m=memory_depth-1; m>0; m--)
	  next_elm[m][last_elm] = next_elm[m-1][last_elm];
	next_elm[0][last_elm] = e;
      }
      last_elm = e;
      return e;
    }
  }

  // failed to find an elm containing coordinates.
  return -1;

}

int m3dc1_mesh::in_element(double X, double Phi, double Z, 
			   double* xi, double* zi, double* eta,
			   int guess)
{
  if(memory_depth>0) {
    return in_element_memory(X, Phi, Z, xi, zi, eta, guess);
  } else {
    return in_element_threadsafe(X, Phi, Z, xi, zi, eta, guess); 
  }
}

int m3dc1_mesh::shared_nodes(const int i, const int j)
{
  int match = 0;
  double R[3], Z[3];
  const double t = (a[i] + b[i] + c[i])*TOL;

  R[0] = x[i];
  R[1] = x[i] + (a[i]+b[i])*co[i];
  R[2] = x[i] + b[i]*co[i] - c[i]*sn[i];
  Z[0] = z[i];
  Z[1] = z[i] + (a[i]+b[i])*sn[i];
  Z[2] = z[i] + c[i]*co[i] + b[i]*sn[i];

  for(int k=0; k<3; k++) {
    if(fabs(R[k] - x[j]) < t) {
      if(fabs(Z[k] - z[j]) < t) match++;
    } 
    if(fabs(R[k] - (x[j] + (a[j]+b[j])*co[j])) < t) {
      if(fabs(Z[k] - (z[j] + (a[j]+b[j])*sn[j])) < t) match++;
    } 
    if(fabs(R[k] - (x[j] + b[j]*co[j] - c[j]*sn[j])) < t) {
      if(fabs(Z[k] - (z[j] + c[j]*co[j] + b[j]*sn[j])) < t) match++;
    }
  }
  
  return match;
}

bool m3dc1_mesh::elements_are_neighbors(const int i, const int j)
{
  return (shared_nodes(i, j)>=2);
}

/*
bool m3dc1_3d_mesh::elements_are_neighbors(const int i, const int j)
{
  return (shared_nodes(i, j)>=3);
}

int m3dc1_3d_mesh::shared_nodes(const int i, const int j)
{
  const double t = d[i]*TOL;

  if((fabs(phi[i] - phi[j]) < t) ||
     (fabs(phi[i] - phi[j] - period) < t) ||
     (fabs(phi[i] - phi[j] + period) < t)) {
     if(i==105) std::cerr << "Hit! " << m3dc1_mesh::shared_nodes(i,j)
			  << std::endl;
    return m3dc1_mesh::shared_nodes(i,j)*2;

  } else if((fabs(phi[i] + d[i] - phi[j]) < t) ||
	    (fabs(phi[i] + d[i] - phi[j] - period) < t) || 
	    (fabs(phi[i] + d[i] - phi[j] + period) < t)) {
    return m3dc1_mesh::shared_nodes(i,j);
    
  } else if((fabs(phi[i] - (phi[j] + d[j])) < t) || 
	    (fabs(phi[i] - (phi[j] + d[j]) - period) < t) || 
	    (fabs(phi[i] - (phi[j] + d[j]) + period) < t)) {
    return m3dc1_mesh::shared_nodes(i,j);
    
  } else {
    return 0;
  }

  return m3dc1_mesh::shared_nodes(i,j);
}
*/

void m3dc1_mesh::find_neighbors()
{
  std::cerr << "Calculating M3D-C1 mesh connectivity..." << std::endl;
  nneighbors = new int[nelms];
  neighbor = new int*[nelms];

  for(int i=0; i<nelms; i++) {
    nneighbors[i] = 0;
    neighbor[i] = new int[max_neighbors()];
  }

  for(int i=0; i<nelms/nplanes; i++) {
    for(int j=i+1; j<nelms/nplanes; j++) {
      if(elements_are_neighbors(i, j)) {
	if(nneighbors[i] >= max_neighbors()) {
	  std::cerr << "Error: element " << i << " has too many neighbors!"
		    << std::endl;
	} else {
	  neighbor[i][nneighbors[i]] = j;
	  nneighbors[i]++;
	}

	if(nneighbors[j] >= max_neighbors()) {
	  std::cerr << "Error: element " << j << " has too many neighbors!"
		    << std::endl;
	} else {
	  neighbor[j][nneighbors[j]] = i;
	  nneighbors[j]++;
	}
      }
    }

    if(nneighbors[i]==0) {
      std::cerr << "Error: Element " << i << " has 0 neighbors!" 
		<< std::endl;
    }
  }
  std::cerr << "Done calculating M3D-C1 mesh connectivity..." << std::endl;
}

void m3dc1_3d_mesh::find_neighbors()
{
  int k;

  m3dc1_mesh::find_neighbors();
  
  // copy the neighbors found in the first plane to all the other planes
  for(int p=1; p < nplanes; p++) {
    for(int i=0; i<nelms/nplanes; i++) {
      k = i + nelms*p/nplanes;
      nneighbors[k] = nneighbors[i];
      for(int j=0; j<nneighbors[i]; j++) {
	neighbor[k][j] = neighbor[i][j] + nelms*p/nplanes;
      }
    }
  }

  // add the toroidally adjacent elements
  for(int i=0; i<nelms; i++) {
    if((i < nelms/nplanes) && (phi[i] != 0.)) {
      std::cerr << "Element " << i << " wrong plane " << phi[i] << std::endl;
    }

    k = i - nelms/nplanes;
    if(k < 0) k += nelms;
    neighbor[i][nneighbors[i]++] = k;
    
    k = i + nelms/nplanes;
    if(k >= nelms) k -= nelms;
    neighbor[i][nneighbors[i]++] = k;
  }
}


m3dc1_3d_mesh::m3dc1_3d_mesh(int n)
  : m3dc1_mesh(n)
{
  d = new double[n];
  phi = new double[n];
}

m3dc1_3d_mesh::~m3dc1_3d_mesh()
{
  delete[] d;
  delete[] phi;
}
