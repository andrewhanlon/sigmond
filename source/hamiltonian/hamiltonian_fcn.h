#ifndef HAMILTONIAN_FCN_H
#define HAMILTONIAN_FCN_H

#include "Minuit2/FCNBase.h"

#include <vector>

class HamiltonianFcn : public ROOT::Minuit2::FCNBase {

public:

  HamiltonianFcn(const std::vector<double>& meas, const std::vector<double>& pos,
		 const std::vector<double>& mvar) :
    theMeasurements(meas), thePositions(pos), theMVariances(mvar),
    theErrorDef(1.) {}

  ~HamiltonianFcn() {}

  virtual double up() const { return theErrorDef; }
  virtual double operator()(const std::vector<double>&) const = 0;

  std::vector<double> measurements() const { return theMeasurements; }
  std::vector<double> positions() const { return thePositions; }
  std::vector<double> variances() const { return theMVariances; }

  void setErrorDef(double def) { theErrorDef = def; }

  virtual void initializeHamiltonian(

protected:

  std::vector<double> theMeasurements;
  std::vector<double> thePositions;
  std::vector<double> theMVariances;
  double theErrorDef;
  
};

class Hamiltonian1b1cFcn : public HamiltonianFcn {

public:
  
  Hamiltonian1b1cFcn(const std::vector<double>& meas, const std::vector<double>& pos,
		     const std::vector<double>& mvar);

  virtual double operator()(const std::vector<double>&) const;

private:

  Hamiltonian1b1c Ham;
};

#endif //HAMILTONIAN_FCN_H
