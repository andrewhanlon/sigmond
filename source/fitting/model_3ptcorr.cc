#include "model_3ptcorr.h"
#include <cmath>
#include <string>
#include "task_utils.h"
using namespace std;


// ******************************************************************************

void create_3ptcorr_model(const string& modeltype, uint in_Tperiod,
                          ThreePointCorrelatorModel* &mptr)
{
  if (modeltype == "TwoStateRatio") {
    mptr = new TwoStateRatio(in_Tperiod);
  }
  else if (modeltype == "TwoStateNumerator") {
    mptr = new TwoStateNumerator(in_Tperiod);
  }
  else if (modeltype == "TwoStateNumeratorSymmetric") {
    mptr = new TwoStateNumeratorSymmetric(in_Tperiod);
  }
  else if (modeltype == "TwoStateSimple") {
    mptr = new TwoStateSimple(in_Tperiod);
  }
  else if (modeltype == "SummationLinear") {
    mptr = new SummationLinear(in_Tperiod);
  }
  else {
    mptr = 0;
    throw(std::invalid_argument(string("Invalid Model in RealThreePointCorrelatorFit: ")+modeltype));}
}

// ******************************************************************************

void ThreePointCorrelatorModel::setupInfos(XMLHandler& xmlm, vector<MCObsInfo>& fitparam_info, int taskcount) const
{
  try {
    fitparam_info.resize(getNumberOfParams());
    int param_count = 0;
    for (vector<string>::const_iterator param_name_it=param_names.begin(); param_name_it!=param_names.end(); ++param_name_it, ++param_count) {
      XMLHandler xmlparam(xmlm, *param_name_it);
      string name; int index;
      xmlreadchild(xmlparam, "Name", name);
      if (name.empty()) throw(std::invalid_argument("Must provide name for parameter " + *param_name_it));
      index = taskcount;
      xmlreadifchild(xmlparam, "IDIndex", index);
      fitparam_info[param_count] = MCObsInfo(name,index);
    }

    for (uint k = 0; k < getNumberOfParams(); k++)
    for (uint l = k + 1; l < getNumberOfParams(); l++)
      if (fitparam_info[k] == fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));
  }
  catch(const std::exception& errmsg) {
    throw(std::invalid_argument(string(model_name)+" -- "+string(errmsg.what())));
  }
}

void ThreePointCorrelatorModel::setupInfos(map<string,MCObsInfo> model_params, vector<MCObsInfo>& fitparam_info) const
{
  try {
    fitparam_info.resize(getNumberOfParams());
    int param_count = 0;
    for (vector<string>::const_iterator param_name_it=param_names.begin(); param_name_it!=param_names.end(); ++param_name_it, ++param_count) {
      fitparam_info[param_count] = model_params.at(*param_name_it);
    }

    for (uint k = 0; k < getNumberOfParams(); k++)
    for (uint l = k + 1; l < getNumberOfParams(); l++)
      if (fitparam_info[k] == fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));
  }
  catch(const std::exception& errmsg) {
    throw(std::invalid_argument(string(model_name)+" -- "+string(errmsg.what())));
  }
}



void ThreePointCorrelatorModel::output_tag(XMLHandler& xmlout) const
{
  xmlout.set_root("Model", model_name);
}




// ******************************************************************************


      //       f(tsep, tins) = (B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep))
      //                        / (1 + A*exp(-dE*tsep))
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           B3 = fitparams[3].
      //            A = fitparams[4].
      //           dE = fitparams[5]


void TwoStateRatio::evaluate(const vector<double>& fitparams, double tsep, double tins,
                             double& value) const
{
  eval_func(fitparams[0], fitparams[1], fitparams[2], fitparams[3], fitparams[4], fitparams[5],
            tsep, tins, value);
}


void TwoStateRatio::evalGradient(const vector<double>& fitparams, double tsep, double tins,
                                 vector<double>& grad) const
{
  eval_grad(fitparams[0], fitparams[1], fitparams[2], fitparams[3], fitparams[4], fitparams[5],
            tsep, tins, grad[0], grad[1], grad[2], grad[3], grad[4], grad[5]);
}


void TwoStateRatio::guessInitialParamValues(
                   const vector<double>& data, const map<uint,set<uint> >& tvals,
                   vector<double>& fitparams) const
{
  throw(std::invalid_argument("TwoStateRatio -- Error: no automatic initial guesses"));
}


      //       f(tsep, tins) = (B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep))
      //                        / (1 + A*exp(-dE*tsep))

 
void TwoStateRatio::eval_func(double B0, double B1, double B2, double B3, double A, double dE,
                              double tsep, double tins, double& funcval) const
{
  funcval = (B0 + B1*exp(-dE*(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep)) / (1. + A*exp(-dE*tsep));
}


void TwoStateRatio::eval_grad(double B0, double B1, double B2, double B3, double A, double dE,
                              double tsep, double tins, double& dB0val, double& dB1val,
                              double& dB2val, double& dB3val, double& dAval, double& ddEval) const
{
  double dE_tsep_tins = exp(-dE*(tsep - tins));
  double dE_tsep = exp(-dE*tsep);
  double dE_tins = exp(-dE*tins);
  double funcval;
  eval_func(dE, B0, B1, B2, B3, A, tsep, tins, funcval);
  double inv_denom = 1./(1. + A*exp(-dE*tsep));

  ddEval = inv_denom*(funcval*tsep*A*dE_tsep - B1*(tsep - tins)*dE_tsep_tins - B2*tins*dE_tins - B3*tsep*dE_tsep);
  dB0val = inv_denom;
  dB1val = inv_denom*dE_tsep_tins;
  dB2val = inv_denom*dE_tins;
  dB3val = inv_denom*dE_tsep;
  dAval = inv_denom*funcval*tsep*A*dE_tsep;
}


// ******************************************************************************


      //       f(tsep, tins) = B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep)
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           B3 = fitparams[3].
      //           dE = fitparams[4]


void TwoStateNumerator::evaluate(const vector<double>& fitparams, double tsep, double tins,
                             double& value) const
{
  eval_func(fitparams[0], fitparams[1], fitparams[2], fitparams[3], fitparams[4],
            tsep, tins, value);
}


void TwoStateNumerator::evalGradient(const vector<double>& fitparams, double tsep, double tins,
                                 vector<double>& grad) const
{
  eval_grad(fitparams[0], fitparams[1], fitparams[2], fitparams[3], fitparams[4],
            tsep, tins, grad[0], grad[1], grad[2], grad[3], grad[4]);
}


void TwoStateNumerator::guessInitialParamValues(
                   const vector<double>& data, const map<uint,set<uint> >& tvals,
                   vector<double>& fitparams) const
{
  throw(std::invalid_argument("TwoStateNumerator -- Error: no automatic initial guesses"));
}


      //       f(tsep, tins) = B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep)

 
void TwoStateNumerator::eval_func(double B0, double B1, double B2, double B3, double dE,
                              double tsep, double tins, double& funcval) const
{
  funcval = B0 + B1*exp(-dE*(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep);
}


void TwoStateNumerator::eval_grad(double B0, double B1, double B2, double B3, double dE,
                                  double tsep, double tins, double& dB0val, double& dB1val,
                                  double& dB2val, double& dB3val, double& ddEval) const
{
  double dE_tsep_tins = exp(-dE*(tsep - tins));
  double dE_tsep = exp(-dE*tsep);
  double dE_tins = exp(-dE*tins);

  ddEval = -B1*(tsep - tins)*dE_tsep_tins - B2*tins*dE_tins - B3*tsep*dE_tsep;
  dB0val = 1.;
  dB1val = dE_tsep_tins;
  dB2val = dE_tins;
  dB3val = dE_tsep;
}


// ******************************************************************************


      //       f(tsep, tins) = B0 + B1*[exp(-dE(tsep - tins)) + exp(-dE*tins)] + B2*exp(-dE*tsep)
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           dE = fitparams[4]


void TwoStateNumeratorSymmetric::evaluate(const vector<double>& fitparams, double tsep, double tins,
                             double& value) const
{
  eval_func(fitparams[0], fitparams[1], fitparams[2], fitparams[3], 
            tsep, tins, value);
}


void TwoStateNumeratorSymmetric::evalGradient(const vector<double>& fitparams, double tsep, double tins,
                                 vector<double>& grad) const
{
  eval_grad(fitparams[0], fitparams[1], fitparams[2], fitparams[3],
            tsep, tins, grad[0], grad[1], grad[2], grad[3]);
}


void TwoStateNumeratorSymmetric::guessInitialParamValues(
                   const vector<double>& data, const map<uint,set<uint> >& tvals,
                   vector<double>& fitparams) const
{
  throw(std::invalid_argument("TwoStateNumeratorSymmetric -- Error: no automatic initial guesses"));
}


      //       f(tsep, tins) = B0 + B1*[exp(-dE(tsep - tins)) + exp(-dE*tins)] + B2*exp(-dE*tsep)

 
void TwoStateNumeratorSymmetric::eval_func(double B0, double B1, double B2, double dE,
                              double tsep, double tins, double& funcval) const
{
  funcval = B0 + B1*(exp(-dE*(tsep - tins)) + exp(-dE*tins)) + B2*exp(-dE*tsep);
}


void TwoStateNumeratorSymmetric::eval_grad(double B0, double B1, double B2, double dE,
                                  double tsep, double tins, double& dB0val, double& dB1val,
                                  double& dB2val, double& ddEval) const
{
  double dE_tsep_tins = exp(-dE*(tsep - tins));
  double dE_tsep = exp(-dE*tsep);
  double dE_tins = exp(-dE*tins);

  ddEval = -B1*((tsep - tins)*dE_tsep_tins + tins*dE_tins) - B2*tsep*dE_tsep;
  dB0val = 1.;
  dB1val = dE_tsep_tins + dE_tins;
  dB2val = dE_tsep;
}

// ******************************************************************************


      //       f(tsep, tins) = B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins)
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           dE = fitparams[3]


void TwoStateSimple::evaluate(const vector<double>& fitparams, double tsep, double tins,
                             double& value) const
{
  eval_func(fitparams[0], fitparams[1], fitparams[2], fitparams[3],
            tsep, tins, value);
}


void TwoStateSimple::evalGradient(const vector<double>& fitparams, double tsep, double tins,
                                 vector<double>& grad) const
{
  eval_grad(fitparams[0], fitparams[1], fitparams[2], fitparams[3],
            tsep, tins, grad[0], grad[1], grad[2], grad[3]);
}


void TwoStateSimple::guessInitialParamValues(
                   const vector<double>& data, const map<uint,set<uint> >& tvals,
                   vector<double>& fitparams) const
{
  throw(std::invalid_argument("TwoStateSimple -- Error: no automatic initial guesses"));
}


      //       f(tsep, tins) = B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins)

 
void TwoStateSimple::eval_func(double B0, double B1, double B2, double dE,
                              double tsep, double tins, double& funcval) const
{
  funcval = B0 + B1*exp(-dE*(tsep - tins)) + B2*exp(-dE*tins);
}


void TwoStateSimple::eval_grad(double B0, double B1, double B2, double dE,
                               double tsep, double tins, double& dB0val, double& dB1val,
                               double& dB2val, double& ddEval) const
{
  double dE_tsep_tins = exp(-dE*(tsep - tins));
  double dE_tins = exp(-dE*tins);

  ddEval = -B1*(tsep - tins)*dE_tsep_tins - B2*tins*dE_tins;
  dB0val = 1.;
  dB1val = dE_tsep_tins;
  dB2val = dE_tins;
}

// ******************************************************************************


      //       f(tsep) = A + B0*tsep
      //
      // where 
      //           B0 = fitparams[0].
      //            A = fitparams[1].


void SummationLinear::evaluate(const vector<double>& fitparams, double tsep, double tins,
                             double& value) const
{
  eval_func(fitparams[0], fitparams[1], tsep, value);
}


void SummationLinear::evalGradient(const vector<double>& fitparams, double tsep, double tins,
                                 vector<double>& grad) const
{
  eval_grad(fitparams[0], fitparams[1], tsep, grad[0], grad[1]);
}


void SummationLinear::guessInitialParamValues(
                   const vector<double>& data, const map<uint,set<uint> >& tvals,
                   vector<double>& fitparams) const
{
  throw(std::invalid_argument("SummationLinear -- Error: no automatic initial guesses"));
}


      //       f(tsep) = A + B0*tsep

 
void SummationLinear::eval_func(double B0, double A, double tsep, double& funcval) const
{
  funcval = A + B0*tsep;
}


void SummationLinear::eval_grad(double B0, double A, double tsep, double& dB0val, double& dAval) const
{
  dAval = 1.;
  dB0val = tsep;
}



 // ***********************************************************************************
