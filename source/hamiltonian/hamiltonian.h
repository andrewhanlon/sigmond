#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>

/*******************************************************************************
 * This file contains the declaration for the Hamiltonian classes.
 * The Hamiltonians are constructed following Phys. Rev. C 90, 055206 (2014)
 * There is an abstract base class from which the concrete Hamiltonians are
 * derived
 *
 * Only one type of concrete class exists thus far: 1b-1c

 ******************************************************************************/

enum Model {A, B, C};

class Hamiltonian {

public:

  Hamiltonian() {}
  
  Hamiltonian(double length, double mass_pi, unsigned int num_momenta, Model model,
	      unsigned int num_bare, unsigned int num_channel);

  ~Hamiltonian() {}

  void test();

private:

  double L, m_pi;
  unsigned int N_k;
  unsigned int N_b, N_c;
  Model mod;

  const double finite_factor = std::pow(2.*M_PI/L, 3./2.) / sqrt(4.*M_PI);

  Eigen::MatrixXcf H0, HI;
  std::vector<double> spectrum;
  std::vector<double> params;
  std::vector<double> C3_vals;

  double get_momenta(int n) const
  { return std::sqrt(n) * 2. * M_PI / L; }

  unsigned int C3(int n) const;
};



#endif // HAMILTONIAN_H
