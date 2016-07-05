#include "hamiltonian.h"

#include <iostream>

using namespace std;
using Eigen::MatrixXcf;
using Eigen::SelfAdjointEigenSolver;


void Hamiltonian::test()
{
  cout << C3_vals[N_k-1] << endl;
}


Hamiltonian::Hamiltonian(double length, double mass_pi, unsigned int num_momenta, Model model,
			 unsigned int num_bare, unsigned int num_channel)
{
  // Set private variables
  L = length;
  m_pi = mass_pi;
  N_k = num_momenta;
  mod = model;
  N_b = num_bare;
  N_c = num_channel;
  
  // Calculate all finite factors
  C3_vals.resize(N_k);
  for (unsigned int n=0; n < N_k; ++n) {
    C3_vals[n] = C3(n);
    cout << n << ": " << C3(n) << ", " << C3_vals[n] << endl;
  }
}

unsigned int Hamiltonian::C3(int n) const {

  int num_ways = 0;
  int max = sqrt(n);

  // There is definitely a faster way to do this.
  for (int n1 = -max; n1 <= max; ++n1) {
    for (int n2 = -max; n2 <= max; ++n2) {
      for (int n3 = -max; n3 <= max; ++n3) {
	if (n1*n1 + n2*n2 + n3*n3 == n) {
	  ++num_ways;
	}
      }
    }
  }

  return num_ways;
}



