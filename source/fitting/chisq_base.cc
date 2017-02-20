#include "chisq_base.h"
#include "task_utils.h"
using namespace std;

// *************************************************************


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
       cov(j,k)=m_obs->getCovariance(m_obs_info[j],m_obs_info[k]);
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
       cov(j,k)=m_obs->getCovariance(m_obs_info[j],m_obs_info[k]);
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
}


double ChiSquare::evalChiSquare(const vector<double>& residuals) const
{
 double tmp=0.0;
 for (uint k=0;k<m_nobs;++k)
    tmp+=residuals[k]*residuals[k];
 return tmp;
}



// *************************************************************
