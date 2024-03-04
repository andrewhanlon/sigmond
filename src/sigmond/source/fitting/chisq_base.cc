#include "chisq_base.h"
#include "task_utils.h"
using namespace std;

// *************************************************************
/*#ChiSquare::ChiSquare(const ChiSquare& cs)
#{
#  m_obs = cs.m_obs;
#  m_nobs = cs.m_nobs;
#  m_nparams = cs.m_nparams;
#  m_obs_info = cs.m_obs_info;
#  m_fitparam_info = cs.m_fitparam_info;
#  m_priors = cs.m_priors;
#  m_means = cs.m_means;
#  m_inv_cov_cholesky = cs.m_inv_cov_cholesky;
#}*/

void ChiSquare::allocate_obs_memory()
{
 m_obs_info.resize(m_nobs);
 m_means.resize(m_nobs);
 m_inv_cov_cholesky.resize(m_nobs);
}



void ChiSquare::output(XMLHandler& xmlout) const
{
 do_output(xmlout);
}


std::string ChiSquare::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}


std::string ChiSquare::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}


void ChiSquare::setObsMean()
{
 try{
    for (uint k=0;k<m_nobs;++k)
       m_means[k]=m_obs->getCurrentSamplingValue(m_obs_info[k]);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Could not setObsMean: ")
      +string(errmsg.what())));}
}


void ChiSquare::setObsMeanCov()
{
 try{
    for (uint k=0;k<m_nobs;++k)
       m_means[k]=m_obs->getCurrentSamplingValue(m_obs_info[k]);
    RealSymmetricMatrix cov(m_nobs);
    for (uint k=0;k<m_nobs;++k)
    for (uint j=0;j<=k;++j)
       if ((j==k)||(m_obs->isCorrelated()))
          cov(j,k)=m_obs->getCovariance(m_obs_info[j],m_obs_info[k]);
       else
          cov(j,k)=0.;
    CholeskyDecomposer CHD;
    CHD.getCholeskyOfInverse(cov,m_inv_cov_cholesky);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Could not setObsMeanCov: ")
      +string(errmsg.what())));}
}


void ChiSquare::setObsMeanCov(RVector& coveigvals)
{
 try{
    for (uint k=0;k<m_nobs;++k)
       m_means[k]=m_obs->getCurrentSamplingValue(m_obs_info[k]);
    RealSymmetricMatrix cov(m_nobs);
    for (uint k=0;k<m_nobs;++k)
    for (uint j=0;j<=k;++j)
       if ((j==k)||(m_obs->isCorrelated()))
          cov(j,k)=m_obs->getCovariance(m_obs_info[j],m_obs_info[k]);
       else
          cov(j,k)=0.;
    Diagonalizer Dc;
    Dc.getEigenvalues(cov,coveigvals);
    CholeskyDecomposer CHD;
    CHD.getCholeskyOfInverse(cov,m_inv_cov_cholesky);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Could not setObsMeanCov: ")
             +string(errmsg.what())));}
}


void ChiSquare::guessInitialFitParamValues(vector<double>& fitparams)
{
 setObsMeanCov();
 guessInitialParamValues(m_means,fitparams);
}

void ChiSquare::addPriors(map<string,Prior> in_priors)
{
  map<string,Prior>::iterator prior_it;
  for (uint param_i = 0; param_i < m_nparams; ++param_i) {
    string param_name = getParameterName(param_i);
    prior_it = in_priors.find(param_name);
    if (prior_it != in_priors.end()) {
      m_priors.insert(pair<uint,Prior>(param_i, prior_it->second));
      m_npriors++;
    }
  }
  int dof = m_nobs-m_nparams+m_npriors;
  if (dof < 1) throw(std::invalid_argument("Degrees of Freedom must be greater than zero"));
}

void ChiSquare::evalResiduals(const vector<double>& fitparams,
                              vector<double>& residuals) const
{
 evalModelPoints(fitparams,residuals);
 for (uint k=0;k<m_nobs;++k)
    residuals[k]-=m_means[k];
 for (int i=m_nobs-1;i>=0;--i){
    double tmp=0.0;
    for (int j=0;j<=i;++j)
       tmp+=m_inv_cov_cholesky(i,j)*residuals[j];
    residuals[i]=tmp;}
 int i=m_nobs;
 for (map<uint,Prior>::const_iterator prior_it=m_priors.begin(); prior_it!=m_priors.end(); ++prior_it,++i){
   //  residuals[i]=(fitparams[prior_it->first] - prior_it->second.mean()) / (prior_it->second.error());
    residuals[i]=prior_it->second.evalPriorResidual(fitparams[prior_it->first]);
 }
}



void ChiSquare::evalResGradients(const vector<double>& fitparams,
                                 RMatrix& gradients) const
{
 evalGradients(fitparams,gradients);
 for (uint p=0;p<m_nparams;++p)
 for (int i=m_nobs-1;i>=0;--i){
    double tmp=0.0;
    for (int j=0;j<=i;++j)
       tmp+=m_inv_cov_cholesky(i,j)*gradients(j,p);
    gradients(i,p)=tmp;}
 int i=m_nobs;
 for (map<uint,Prior>::const_iterator prior_it=m_priors.begin(); prior_it!=m_priors.end(); ++prior_it,++i){
   //  gradients(i,prior_it->first)=1./(prior_it->second.error());
    gradients(i,prior_it->first)=prior_it->second.evalPriorGradient(fitparams[prior_it->first]);
 }
}


double ChiSquare::evalChiSquare(const vector<double>& residuals) const
{
 double tmp=0.0;
 for (uint k=0;k<m_nobs+m_npriors;++k)
    tmp+=residuals[k]*residuals[k];
 return tmp;
}



// *************************************************************
