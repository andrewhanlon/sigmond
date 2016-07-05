#include "hamiltonian_fcn.h"
#include "hamiltonian.h"

#include <cassert>

using namespace std;

Hamiltonian1b1cFcn::Hamiltonian1b1cFcn(const vector<double>& meas, const vector<double>& pos,
				       const vector<double>& mvar) {

  theMeasurements = meas;
  thePositions = pos;
  theMVariances = mvar;
  theErrorDef = 1.;
  
}

double Hamiltonian1b1cFcn::operator()(const vector<double>& par) const {

  assert(par.size() == 5);
}
