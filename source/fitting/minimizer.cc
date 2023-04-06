#include "minimizer.h"
#include <iostream>
using namespace std;


// ***************************************************************************


ChiSquareMinimizerInfo::ChiSquareMinimizerInfo(XMLHandler& xmlin)
{
 XMLHandler xmlr(xmlin,"MinimizerInfo");
 m_method=defaultmethod;
 string reply;
 if (xmlreadifchild(xmlr,"Method",reply)){
    if (reply=="LMDer") m_method='L';
#ifndef NO_MINUIT
    else if (reply=="Minuit2") m_method='M';
    else if (reply=="Minuit2NoGradient") m_method='F';
#endif
    else if (reply=="NL2Sol") m_method='N';
    else throw(std::invalid_argument("Invalid <Method> tag in ChiSquareMinimizerInfo"));}
 m_param_reltol=1e-6;
 double tol;
 if (xmlreadifchild(xmlr,"ParameterRelTol",tol)) 
    setParameterRelativeTolerance(tol);
 m_chisq_reltol=1e-4;
 if (xmlreadifchild(xmlr,"ChiSquareRelTol",tol)) 
    setChiSquareRelativeTolerance(tol);
 m_max_its=1024;
 int its;
 if (xmlreadifchild(xmlr,"MaximumIterations",its)) 
    setMaximumIterations(its);
 m_verbosity='L';
 if (xmlreadifchild(xmlr,"Verbosity",reply)){
    if (reply=="Low") m_verbosity='L';
    else if (reply=="Medium") m_verbosity='M';
    else if (reply=="High") m_verbosity='H';
    else throw(std::invalid_argument("Invalid <Verbosity> tag in ChiSquareMinimizerInfo"));}
}



ChiSquareMinimizerInfo::ChiSquareMinimizerInfo(char method, 
                                               double param_reltol, double chisq_reltol, 
                                               int max_its, char verbosity)
{
 setMethod(method);
 setParameterRelativeTolerance(param_reltol);
 setChiSquareRelativeTolerance(chisq_reltol);
 setMaximumIterations(max_its);
 setVerbosity(verbosity);
}


ChiSquareMinimizerInfo::ChiSquareMinimizerInfo(const ChiSquareMinimizerInfo& info)
   :   m_method(info.m_method),  m_param_reltol(info.m_param_reltol),
       m_chisq_reltol(info.m_chisq_reltol),  m_max_its(info.m_max_its),
       m_verbosity(info.m_verbosity)
{}


ChiSquareMinimizerInfo& ChiSquareMinimizerInfo::operator=(const ChiSquareMinimizerInfo& info)
{
 m_method=info.m_method;      
 m_param_reltol=info.m_param_reltol;
 m_chisq_reltol=info.m_chisq_reltol;
 m_max_its=info.m_max_its;
 m_verbosity=info.m_verbosity;   
 return *this;
}


void ChiSquareMinimizerInfo::setMethod(char method)
{
 if ((method=='L')||(method=='N')
#ifndef NO_MINUIT
    ||(method=='M')||(method=='F')
#endif
    ){
    m_method=method;
    return;}
 m_method=defaultmethod;
 throw(std::invalid_argument("Invalid method character in ChiSquareMinimizerInfo::setMethod"));
}


void ChiSquareMinimizerInfo::setMinuit2()
{
#ifndef NO_MINUIT
 m_method='M';
#else
 throw(std::invalid_argument("Minuit2 library not available in ChiSquareMinimizerInfo::setMethod"));
#endif
}

void ChiSquareMinimizerInfo::setMinuit2NoGradient()
{
#ifndef NO_MINUIT
 m_method='F';
#else
 throw(std::invalid_argument("Minuit2 library not available in ChiSquareMinimizerInfo::setMethod"));
#endif
}

void ChiSquareMinimizerInfo::setParameterRelativeTolerance(double rtol)
{
 if (rtol>0.0){
    m_param_reltol=rtol; return;}
 throw(std::invalid_argument("Invalid input in ChiSquareMinimizerInfo::setParameterRelativeTolerance"));
}


void ChiSquareMinimizerInfo::setChiSquareRelativeTolerance(double rtol)
{
 if (rtol>0.0){
    m_chisq_reltol=rtol; return;}
 throw(std::invalid_argument("Invalid input in ChiSquareMinimizerInfo::setChiSquareRelativeTolerance"));
}


void ChiSquareMinimizerInfo::setMaximumIterations(unsigned int maxit)
{
 if (maxit>10){
    m_max_its=maxit; return;}
 throw(std::invalid_argument("Invalid input in ChiSquareMinimizerInfo::setMaximumIterations"));
}


void ChiSquareMinimizerInfo::setVerbosity(char verbosity)
{
 if ((verbosity=='L')||(verbosity=='M')||(verbosity=='H')){
    m_verbosity=verbosity; return;}
 throw(std::invalid_argument("Invalid verbosity character in ChiSquareMinimizerInfo::setVerbosity"));
}



void ChiSquareMinimizerInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("MinimizerInfo");
 if (m_method=='M') xmlout.put_child("Method","Minuit2");
 else if (m_method=='F') xmlout.put_child("Method","Minuit2NoGradient");
 else if (m_method=='L') xmlout.put_child("Method","LMDer");
 else if (m_method=='N') xmlout.put_child("Method","NL2Sol");
 xmlout.put_child("ParameterRelTol",make_string(m_param_reltol));
 xmlout.put_child("ChiSquareRelTol",make_string(m_chisq_reltol));
 xmlout.put_child("MaximumIterations",make_string(m_max_its));
 if (m_verbosity=='L') xmlout.put_child("Verbosity","Low");
 else if (m_verbosity=='M') xmlout.put_child("Verbosity","Medium");
 else if (m_verbosity=='H') xmlout.put_child("Verbosity","High");
}


string ChiSquareMinimizerInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}


string ChiSquareMinimizerInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}


// ***************************************************************************


ChiSquareMinimizer::ChiSquareMinimizer(ChiSquare &in_chisq)
   : m_chisq(&in_chisq),  m_lmder(0), m_nl2sol(0)
#ifndef NO_MINUIT
     , m_minuit2(0), m_minuit2ng(0)
#endif
{
 alloc_method();
}


ChiSquareMinimizer::ChiSquareMinimizer(ChiSquare &in_chisq, const ChiSquareMinimizerInfo& info)
   : m_chisq(&in_chisq),  m_lmder(0), m_nl2sol(0), m_info(info)
#ifndef NO_MINUIT
     , m_minuit2(0), m_minuit2ng(0)
#endif
{
 alloc_method();
}



ChiSquareMinimizer::~ChiSquareMinimizer()
{
 dealloc_method();
}


void ChiSquareMinimizer::dealloc_method()
{
 delete m_lmder; m_lmder=0;
 delete m_nl2sol; m_nl2sol=0;
#ifndef NO_MINUIT
 delete m_minuit2; m_minuit2=0;
 delete m_minuit2ng; m_minuit2ng=0;
#endif
}


void ChiSquareMinimizer::alloc_method()
{
 if (m_info.m_method=='L')
    m_lmder=new LMDerMinimizer(*m_chisq);
 else if (m_info.m_method=='N')
    m_nl2sol=new NL2SolMinimizer(*m_chisq);
#ifndef NO_MINUIT
 else if (m_info.m_method=='M')
    m_minuit2=new Minuit2ChiSquare(*m_chisq);
 else if (m_info.m_method=='F')
    m_minuit2ng=new Minuit2NoGradChiSquare(*m_chisq);
#endif
}



void ChiSquareMinimizer::reset(const ChiSquareMinimizerInfo& info)
{
 if (m_info.m_method==info.m_method) return;
 dealloc_method();
 m_info=info;
 alloc_method();
}



bool ChiSquareMinimizer::find_minimum(const vector<double>& starting_params,
                                      double& chisq_min, 
                                      vector<double>& params_at_minimum,
                                      XMLHandler& xmlout, char verbosity)
{
 if (m_info.m_method=='L')
    return find_minimum_lmder(starting_params,chisq_min,params_at_minimum,
                              xmlout,verbosity);
#ifndef NO_MINUIT
 else if (m_info.m_method=='M')
    return find_minimum_minuit2(starting_params,chisq_min,params_at_minimum,
                                xmlout,verbosity);
 else if (m_info.m_method=='F')
    return find_minimum_minuit2ng(starting_params,chisq_min,params_at_minimum,
                                  xmlout,verbosity);
#endif
 else if (m_info.m_method=='N')
    return find_minimum_nl2sol(starting_params,chisq_min,params_at_minimum,
                               xmlout,verbosity);
 return false;
}



bool ChiSquareMinimizer::findMinimum(const vector<double>& starting_params,
                                     double& chisq_min, 
                                     vector<double>& params_at_minimum,
                                     XMLHandler& xmlout)
{
 return find_minimum(starting_params,chisq_min,params_at_minimum,
                     xmlout,m_info.m_verbosity);
}

bool ChiSquareMinimizer::findMinimum(double& chisq_min, 
                                     vector<double>& params_at_minimum,
                                     XMLHandler& xmlout)
{
 vector<double> starting_params(m_chisq->getNumberOfParams());
 m_chisq->guessInitialFitParamValues(starting_params);
 return findMinimum(starting_params,chisq_min,params_at_minimum,xmlout);
}


bool ChiSquareMinimizer::findMinimum(const vector<double>& starting_params,
                                     double& chisq_min, 
                                     vector<double>& params_at_minimum)
{
 XMLHandler xmlout;
 return find_minimum(starting_params,chisq_min,params_at_minimum,xmlout,'L');
}


bool ChiSquareMinimizer::findMinimum(double& chisq_min, 
                                     vector<double>& params_at_minimum)
{
 vector<double> starting_params(m_chisq->getNumberOfParams());
 m_chisq->guessInitialFitParamValues(starting_params);
 XMLHandler xmlout;
 return find_minimum(starting_params,chisq_min,params_at_minimum,xmlout,'L');
}



#ifndef NO_MINUIT

bool ChiSquareMinimizer::find_minimum_minuit2(const vector<double>& starting_params,
                                              double& chisq_min, 
                                              vector<double>& params_at_minimum,
                                              XMLHandler& xmlout, char verbosity)
{
 xmlout.clear();
 uint nparam=m_chisq->getNumberOfParams();
 if (starting_params.size()!=nparam)
    throw(std::invalid_argument("Invalid starting parameters"));

 std::vector<double> unc(nparam);
 for (uint p=0;p<nparam;++p)
    unc[p]=0.01*starting_params[p];   // set up initial uncertainties

 unsigned int strategylevel=2;  // 0 = low, 1 = med, 2 = high quality
                                // lower level means faster, higher means
                                // more reliable minimization

 ROOT::Minuit2::MnMinimize M(*m_minuit2, starting_params, unc, strategylevel); 

       //  Now do the actual minimization!!
 ROOT::Minuit2::FunctionMinimum csmin = M(m_info.m_max_its,m_info.m_chisq_reltol);

 if (verbosity!='L'){
    ostringstream outlog;
    outlog<<"Minuit2 Minimization Result:\n"<<csmin;
    xmlformat("Minuit2Log",outlog.str(),xmlout);}

 if (csmin.IsValid()){
    chisq_min=csmin.Fval();
    params_at_minimum.resize(nparam);
    for (uint p=0;p<nparam;++p)
        params_at_minimum[p]=csmin.UserParameters().Value(p);}
 else{
    chisq_min=-1.0;
    params_at_minimum.clear();}

 return csmin.IsValid();
}

bool ChiSquareMinimizer::find_minimum_minuit2ng(const vector<double>& starting_params,
                                               double& chisq_min, 
                                               vector<double>& params_at_minimum,
                                               XMLHandler& xmlout, char verbosity)
{
 xmlout.clear();
 uint nparam=m_chisq->getNumberOfParams();
 if (starting_params.size()!=nparam)
    throw(std::invalid_argument("Invalid starting parameters"));

 std::vector<double> unc(nparam);
 for (uint p=0;p<nparam;++p)
    unc[p]=0.01*starting_params[p];   // set up initial uncertainties

 unsigned int strategylevel=2;  // 0 = low, 1 = med, 2 = high quality
                                // lower level means faster, higher means
                                // more reliable minimization

 ROOT::Minuit2::MnMinimize M(*m_minuit2ng, starting_params, unc, strategylevel); 

       //  Now do the actual minimization!!
 ROOT::Minuit2::FunctionMinimum csmin = M(m_info.m_max_its,m_info.m_chisq_reltol);

 if (verbosity!='L'){
    ostringstream outlog;
    outlog<<"Minuit2 Minimization Result:\n"<<csmin;
    xmlformat("Minuit2Log",outlog.str(),xmlout);}

// if (!(csmin.IsValid())){
//    M.SetPrecision(1e-12);
//    csmin = M(m_info.m_max_its,m_info.m_chisq_reltol);}

 if (csmin.IsValid()){
    chisq_min=csmin.Fval();
    params_at_minimum.resize(nparam);
    for (uint p=0;p<nparam;++p)
        params_at_minimum[p]=csmin.UserParameters().Value(p);}
 else{
    chisq_min=-1.0;
    params_at_minimum.clear();}

 return csmin.IsValid();
}

#endif


bool ChiSquareMinimizer::find_minimum_lmder(const vector<double>& starting_params,
                                            double& chisq_min,
                                            vector<double>& params_at_minimum,
                                            XMLHandler& xmlout, char verbosity)
{
 xmlout.clear();
 ostringstream outlog;
 m_lmder->setInitialFitParamValues(starting_params);
 int flag=m_lmder->chisq_fit(m_info.m_param_reltol,m_info.m_chisq_reltol,verbosity, 
                             m_info.m_max_its,chisq_min,outlog);
 if (verbosity!='L')
    xmlformat("LMDerLog",outlog.str(),xmlout);

 if (flag<=2){
    params_at_minimum=m_lmder->getFitParams();
    return true;}

 params_at_minimum.clear();
 chisq_min=-1.0;
 return false;
}



bool ChiSquareMinimizer::find_minimum_nl2sol(const vector<double>& starting_params,
                                             double& chisq_min, 
                                             vector<double>& params_at_minimum,
                                             XMLHandler& xmlout, char verbosity)
{
 xmlout.clear();
 ostringstream outlog;
 m_nl2sol->setInitialFitParamValues(starting_params);
 int flag=m_nl2sol->chisq_fit(m_info.m_param_reltol,m_info.m_chisq_reltol,m_info.m_verbosity, 
                             m_info.m_max_its,chisq_min,outlog);

 if (verbosity!='L')
    xmlformat("NL2SolLog",outlog.str(),xmlout);

 if (flag<=2){
    params_at_minimum=m_nl2sol->getFitParams();
    return true;}

 params_at_minimum.clear();
 chisq_min=-1.0;
 return false;
}


void ChiSquareMinimizer::xmlformat(const string& roottag, const string& inlogstr,
                                   XMLHandler& xmlout)
{
 string logstr(inlogstr);
 list<string> outlines;
 size_t pos=logstr.find_first_of("\n");
 size_t maxlength=0;
 while (pos!=string::npos){
    if (pos>maxlength) maxlength=pos;
    outlines.push_back(logstr.substr(0,pos));
    logstr.erase(0,pos+1);
    pos=logstr.find_first_of("\n");}
 xmlout.set_root(roottag);
 for (list<string>::iterator it=outlines.begin();it!=outlines.end();++it){
    it->resize(maxlength,' ');
    *it+=" :";
    xmlout.put_child("o",":  "+*it);}
}

// ********************************************************************



#ifndef NO_MINUIT

void Minuit2ChiSquare::guessInitialFitParamValues(vector<double>& params)
{
 params.resize(m_nparams);
 m_chisq->guessInitialFitParamValues(params);
}


double Minuit2ChiSquare::operator()(const vector<double>& params) const
{
 m_chisq->evalResiduals(params,m_residuals);
 return m_chisq->evalChiSquare(m_residuals);
}

double Minuit2ChiSquare::evalChiSquare(const vector<double>& params) const
{
 m_chisq->evalResiduals(params,m_residuals);
 return m_chisq->evalChiSquare(m_residuals); //add priors?
}

vector<double> Minuit2ChiSquare::Gradient(const vector<double>& params) const
{
 m_chisq->evalResiduals(params,m_residuals);
 m_chisq->evalResGradients(params,m_gradients);
 vector<double> grad(m_nparams);
 for (uint p=0;p<m_nparams;++p){
    double tmp=0.0;
    for (uint k=0;k<m_nobs;++k)
       tmp+=m_residuals[k]*m_gradients(k,p);
    grad[p]=2.0*tmp;}
 return grad;
}

void Minuit2NoGradChiSquare::guessInitialFitParamValues(vector<double>& params)
{
 params.resize(m_nparams);
 m_chisq->guessInitialFitParamValues(params);
}


double Minuit2NoGradChiSquare::operator()(const vector<double>& params) const
{
 m_chisq->evalResiduals(params,m_residuals);
 return m_chisq->evalChiSquare(m_residuals);
}

double Minuit2NoGradChiSquare::evalChiSquare(const vector<double>& params) const
{
 m_chisq->evalResiduals(params,m_residuals);
 return m_chisq->evalChiSquare(m_residuals);
}


#endif

// ********************************************************************



namespace MinPack {

int lmder(ChiSquare *p, int m, int n, vector<double>& params, 
          vector<double>& residuals, RMatrix& gradients, 
          int ldfjac, double ftol, double xtol, double gtol, int maxfev, 
          int mode, double factor, int nprint, int &nfev, 
          int &njev, vector<double>& vdiag, vector<double>& vwa1, 
          vector<double>& vwa2, vector<double>& vwa3, vector<double>& vqtf,
          vector<double> vwa4, vector<int> vipvt, ostringstream& outlog);

double enormsq(int n, double *x);

int nl2sol(ChiSquare *M, int n, int p, vector<double>& params,
           vector<double>& residuals, RMatrix& gradients,
           vector<int>& viv, vector<double>& vv, ostringstream& outlog);

int dfault(int *iv, double *v);
};


int LMDerMinimizer::chisq_fit(double paramreltol, double chisqreltol, char verbosity, 
                              int max_its, double& chisq, ostringstream& outlog)
{
 int nparams=m_chisq->getNumberOfParams();
 int nobs=m_chisq->getNumberOfObservables();

 double gtol = 0.0;
 double factor=100.0;
 int mode = 1;
 int nprint=(verbosity=='L')?0:((verbosity=='M')?10:1);
 int nfev,njev;

 int info=MinPack::lmder(m_chisq,nobs,nparams,m_fitparams,m_residuals,
                         m_gradients,nobs,chisqreltol,paramreltol,gtol,max_its,
                         mode,factor,nprint,nfev,njev,vdiag,vwa1,
                         vwa2,vwa3,vqtf,vwa4,vipvt,outlog);

    // return codes

 if (info==3) info=0;             // both x and f convergence
 else if (info==2) info=1;        // x convergence only
 else if (info==1) info=2;        // f convergence only
 else if (info==5) info=3;        // max its exceeded
 else if ((info==4)||(info==6)||(info==7)||(info==8)) info=4; // false
 else info=5;                    // bad input or miscellaneous

 chisq = m_chisq->evalChiSquare(m_residuals);
 return info;
}



void LMDerMinimizer::guessInitialFitParamValues()
{
 m_chisq->guessInitialFitParamValues(m_fitparams);
}



void LMDerMinimizer::setInitialFitParamValues(const std::vector<double>& start_params)
{
 if (start_params.size()!=m_nparams)
    throw(std::invalid_argument("Invalid starting parameters"));
 for (uint k=0;k<m_nparams;++k)
    m_fitparams[k]=start_params[k];
}


// ***************************************************************


int NL2SolMinimizer::chisq_fit(double paramreltol, double chisqreltol, char verbosity, 
                               int max_its, double& chisq, ostringstream& outlog)
{
 int p=m_chisq->getNumberOfParams();
 int n=m_chisq->getNumberOfObservables();

      //  set up nl2sol parameters

 MinPack::dfault(&iv[0],&v[0]);
 iv[16]=4*max_its;  // max func evals
 iv[17]=max_its;    // max iterations
 iv[18]=(verbosity=='L')?0:((verbosity=='H')?1:10);      // output flag
 v[31]=chisqreltol;         // relative func tolerance
 v[32]=paramreltol;         // relative solution tolerance

      //  perform minimization

 int info=MinPack::nl2sol(m_chisq,n,p,m_fitparams,m_residuals,
                          m_gradients,iv,v,outlog);

     //  return code
 
 if (iv[0]==5) info=0;                       // both x and f convergence
 else if (iv[0]==3) info=1;                  // x convergence only
 else if (iv[0]==4) info=2;                  // f convergence only
 else if ((iv[0]==9)||(iv[0]==10)) info=3;   // max its exceeded
 else if ((iv[0]==7)||(iv[0]==8)) info=4;    // false or singular convergence
 else info=5;                                // bad input or miscellaneous

 chisq=2.0*v[9]; 

 return info;
}

void NL2SolMinimizer::guessInitialFitParamValues()
{
 m_chisq->guessInitialFitParamValues(m_fitparams);
}


void NL2SolMinimizer::setInitialFitParamValues(const std::vector<double>& start_params)
{
 if (start_params.size()!=m_nparams)
    throw(std::invalid_argument("Invalid starting parameters"));
 for (uint k=0;k<m_nparams;++k)
    m_fitparams[k]=start_params[k];
}


