#include "model_disp.h"
#include <cmath>
#include <string>
#include "task_utils.h"
using namespace std;


// ******************************************************************************


void create_disp_model(const string& modeltype, uint in_Xextent, uint in_Yextent, uint in_Zextent,
                       DispersionModel* &mptr)
{
  if (modeltype == "Anisotropy") {
    mptr = new Anisotropy(in_Xextent, in_Yextent, in_Zextent);
  }
  else if (modeltype == "Continuum") {
    mptr = new Continuum(in_Xextent, in_Yextent, in_Zextent);
  }
  else if (modeltype == "FourthOrderMomentum") {
    mptr = new FourthOrderMomentum(in_Xextent, in_Yextent, in_Zextent);
  }
  else {
    mptr = 0;
    throw(std::invalid_argument(string("Invalid Model in DispersionFit: ")+modeltype));
  }
}

// ******************************************************************************

void DispersionModel::setupInfos(XMLHandler& xmlm, vector<MCObsInfo>& fitparam_info, int taskcount) const
{
  try {
    fitparam_info.resize(m_nparams);
    int param_count = 0;
    for (auto& param_name: param_names) {
      XMLHandler xmlparam(xmlm, param_name);
      string name; int index;
      xmlreadchild(xmlparam, "Name", name);
      if (name.empty()) throw(std::invalid_argument("Must provide name for parameter " + param_name));
      index = taskcount;
      xmlreadifchild(xmlparam, "IDIndex", index);
      fitparam_info[param_count++] = MCObsInfo(name, index);
    }

    for (uint k = 0; k < m_nparams; k++) {
      for (uint l = k+1; l < m_nparams; l++) {
        if (fitparam_info[k] == fitparam_info[l]) {
          throw(std::invalid_argument("Fit parameter infos must all differ"));
        }
      }
    }
  }
  catch(const std::exception& errmsg) {
    throw(std::invalid_argument(string(model_name)+" -- "+string(errmsg.what())));
  }
}

void DispersionModel::setupInfos(map<string,MCObsInfo> model_params, vector<MCObsInfo>& fitparam_info) const
{
  try {
    fitparam_info.resize(m_nparams);
    int param_count = 0;
    for (auto& param_name: param_names) {
      fitparam_info[param_count++] = model_params.at(param_name);
    }

    for (uint k = 0; k < m_nparams; k++) {
      for (uint l = k+1; l < m_nparams; l++) {
        if (fitparam_info[k] == fitparam_info[l]) {
          throw(std::invalid_argument("Fit parameter infos must all differ"));
        }
      }
    }
  }
  catch(const std::exception& errmsg) {
    throw(std::invalid_argument(string(model_name)+" -- "+string(errmsg.what())));
  }
}


void DispersionModel::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model",model_name);
}



// ******************************************************************************


      // Fitting function is 
      //
      //          (a_t E)^2 = restmass_sq + p^2 / xi^2
      //
      // where 
      //           restmass_sq = fitparams[0]
      //                    xi = fitparams[1].
      //


void Anisotropy::evaluate(const vector<double>& fitparams, double psq, double& value) const
{
  eval_func(fitparams[0], fitparams[1], psq, value);
}


void Anisotropy::evalGradient(const vector<double>& fitparams, double psq, vector<double>& grad) const
{
 eval_grad(fitparams[1], psq, grad[0], grad[1]);
}


void Anisotropy::guessInitialParamValues(
                   const vector<double>& data, const vector<double>& psqvals,
                   vector<double>& fitparams) const
{
  uint p = 0;
  for (uint k = 1; k < psqvals.size(); ++k) {
    if (psqvals[k] != psqvals[0]) {
      p = k;
      break;
    }
  }
  if (p == 0) throw(std::invalid_argument("Could not guess initial parameter values"));
  
  fitparams[1] = sqrt((psqvals[p] - psqvals[0])/(data[p] - data[0]));

  bool flag = false;
  for (uint k = 0; k < psqvals.size(); ++k) {
    if (psqvals[k] == 0.) {
      flag = true;
      fitparams[0] = data[k];
      break;
    }
  }

  if (!flag) {
    fitparams[0] = data[0] - psqvals[0]/(fitparams[1]*fitparams[1]);
  }
}


      //          (a_t E)^2 = restmass_sq + p^2 / xi^2

 
void Anisotropy::eval_func(double msq, double xi, double psq, double& funcval) const
{
  funcval = msq + psq / pow(xi, 2.);
}


void Anisotropy::eval_grad(double xi, double psq, double& dmsqval, double& dxival) const
{
  dmsqval = 1.;
  dxival = -2.*psq/pow(xi, 3.);
}


// ******************************************************************************


      // Fitting function is 
      //
      //          (a_t E)^2 = restmass_sq + c * p^2
      //
      // where 
      //           restmass_sq = fitparams[0]
      //                     c = fitparams[1].
      //


void Continuum::evaluate(const vector<double>& fitparams, double psq, double& value) const
{
  eval_func(fitparams[0], fitparams[1], psq, value);
}


void Continuum::evalGradient(const vector<double>& fitparams, double psq, vector<double>& grad) const
{
 eval_grad(psq, grad[0], grad[1]);
}


void Continuum::guessInitialParamValues(
                   const vector<double>& data, const vector<double>& psqvals,
                   vector<double>& fitparams) const
{
  uint p = 0;
  for (uint k = 1; k < psqvals.size(); ++k) {
    if (psqvals[k] != psqvals[0]) {
      p = k;
      break;
    }
  }
  if (p == 0) throw(std::invalid_argument("Could not guess initial parameter values"));
  
  fitparams[1] = (data[p] - data[0])/(psqvals[p] - psqvals[0]);

  bool flag = false;
  for (uint k = 0; k < psqvals.size(); ++k) {
    if (psqvals[k] == 0.) {
      flag = true;
      fitparams[0] = data[k];
      break;
    }
  }

  if (!flag) {
    fitparams[0] = data[0] - fitparams[1]*psqvals[0];
  }
}


      //          (a_t E)^2 = restmass_sq + c * p^2

 
void Continuum::eval_func(double msq, double c, double psq, double& funcval) const
{
  funcval = msq + c*psq;
}


void Continuum::eval_grad(double psq, double& dmsqval, double& dcval) const
{
  dmsqval = 1.;
  dcval = psq;
}


// ******************************************************************************


      // Fitting function is 
      //
      //      (a_t E)^2 = restmass_sq + c1 * p^2 - c2 * [(restmass_sq)^2 + 2 * p^4 + restmass_sq * p^2]
      //
      // where 
      //           restmass_sq = fitparams[0]
      //                    c1 = fitparams[1].
      //                    c2 = fitparams[2].


void FourthOrderMomentum::evaluate(const vector<double>& fitparams, double psq, double& value) const
{
  eval_func(fitparams[0], fitparams[1], fitparams[2], psq, value);
}


void FourthOrderMomentum::evalGradient(const vector<double>& fitparams, double psq, vector<double>& grad) const
{
 eval_grad(fitparams[0], fitparams[2], psq, grad[0], grad[1], grad[2]);
}


void FourthOrderMomentum::guessInitialParamValues(
                   const vector<double>& data, const vector<double>& psqvals,
                   vector<double>& fitparams) const
{
  uint p = 0;
  for (uint k = 1; k < psqvals.size(); ++k) {
    if (psqvals[k] != psqvals[0]) {
      p = k;
      break;
    }
  }
  if (p == 0) throw(std::invalid_argument("Could not guess initial parameter values"));
  
  fitparams[1] = (data[p] - data[0])/(psqvals[p] - psqvals[0]);
  fitparams[2] = 0.1*fitparams[1];

  bool flag = false;
  for (uint k = 0; k < psqvals.size(); ++k) {
    if (psqvals[k] == 0.) {
      flag = true;
      fitparams[0] = data[k];
      break;
    }
  }

  if (!flag) {
    fitparams[0] = data[0] - fitparams[1]*psqvals[0];
  }
}


      //      (a_t E)^2 = restmass_sq + c1 * p^2 - c2 * [(restmass_sq)^2 + 2 * p^4 + restmass_sq * p^2]

 
void FourthOrderMomentum::eval_func(double msq, double c1, double c2, double psq, double& funcval) const
{
  funcval = msq + c1*psq - c2*(msq*msq + 2.*psq*psq + msq*psq);
}


void FourthOrderMomentum::eval_grad(double msq, double c2, double psq, double& dmsqval, double& dc1val, double& dc2val) const
{
  dmsqval = 1. - c2*(2.*msq + psq);
  dc1val = psq;
  dc2val = msq*msq + 2.*psq*psq + msq*psq;
}



 // ***********************************************************************************
