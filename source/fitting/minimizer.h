#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "chisq_base.h"
#include "minpack.h"

#ifndef NO_MINUIT
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnPrint.h"
#endif

class LMDerMinimizer;
class NL2SolMinimizer;
class Minuit2ChiSquare;
class ChiSquareMinimizerInfo;


// ********************************************************************************
// *                                                                              *
// *   "ChiSquareMinimizer" is the all-important class for performing the         *
// *   correlated least-squares fitting in SigMonD.  A model function with a      *
// *   vector of parameters is assumed to describe a vector of observables        *
// *   obs[j], and the best fit values of the parameters are determined by        *
// *   minimizing the correlated chi-square function:                             *
// *                                                                              *
// *    chi-square = sum_j (model[j]-obs[j]) * inv_cov(j,k) * (model[k]-obs[k])   *
// *                                                                              *
// *   where inv_cov = the inverse of the covariance matrix of the observables.   *
// *   A Cholesky decomposition of inv_cov can be done:                           *
// *              inv_cov = transpose(L) * L,   L = lower triangular              *
// *   The fit parameters determine the model[j] values, and the fit parameters   *
// *   are adjusted until chi-square is a minimum.                                *
// *                                                                              *
// *   Three different methods for performing the minimization are available:     *
// *                                                                              *
// *        Minuit2          method = 'M'                                         *
// *        Netlib's LMDer   method = 'L'                                         *
// *        Netlib's NL2Sol  method = 'N'                                         *
// *                                                                              *
// *   Minuit2 is the most sophisticated of the methods and the most modern       *
// *   software, but there is limited control over tolerances and it is slower,   *
// *   although more reliable.  The other two methods are from old Fortran code   *
// *   and available through NetLib.  The updated C++ code for LMDer and NL2Sol   *
// *   are included in SigMonD, whereas Minuit2 is used as a library.             *
// *                                                                              *
// *   The key routine of the class is "findMinimum", which returns "true" if     *
// *   the minimum is found to within the tolerance requested, or "false"         *
// *   otherwise:                                                                 *
// *                                                                              *
// *         ChiSquare chi(...)  <- specifies the model and observations          *
// *         ChiSquareMinimizer CSM(chi);                                         *
// *                                                                              *
// *         vector<double> init_params(...);                                     *
// *         double chisq;                                                        *
// *         vector<double> params_at_min(...);                                   *
// *         XMLHandler xmlout;                                                   *
// *                                                                              *
// *               // version using specified starting parameter values           *
// *         bool flag=CSM.findMinimum(init_params,chisq,params_at_min,xmlout);   *
// *                                                                              *
// *               // version with initial values specified by ChiSquare          *
// *         bool flag=CSM.findMinimum(chisq,params_at_min,xmlout);               *
// *                                                                              *
// *               // version with starting parameter values, no output           *
// *         bool flag=CSM.findMinimum(init_params,chisq,params_at_min);          *
// *                                                                              *
// *               // version with no output, ChiSquare for init values           *
// *         bool flag=CSM.findMinimum(chisq,params_at_min);                      *
// *                                                                              *
// *                                                                              *
// *   A requested relative tolerance on the solution and the chi-square value    *
// *   can be used.  The maximum number of iterations can be set, although        *
// *   default values are provided for all fitting tolerances, etc.  Set          *
// *   "verbosity" to low 'L' for no output, and to medium 'M' or 'H' for some    *
// *   logging.  These details are passed into the minimizer through an object    *
// *   of class "ChiSquareMinimizerInfo" which can be created with XML of the     *
// *   form:                                                                      *
// *                                                                              *
// *      <MinimizerInfo>                                                         *
// *         <Method>Minuit2</Method>    LMDer, NL2Sol other options              *
// *         <ParameterRelTol>1e-6</ParameterRelTol>                              *
// *         <ChiSquareRelTol>1e-4</ChiSquareRelTol>                              *
// *         <MaximumIterations>1024</MaximumIterations>                          *
// *         <Verbosity>Low</Verbosity>   Medium, High other options              *
// *      </MinimizerInfo>                                                        *
// *                                                                              *
// ********************************************************************************


class ChiSquareMinimizerInfo
{

    char m_method;      // 'L' = lmder, 'N' = nl2sol, 'M' = minuit2
    double m_param_reltol;
    double m_chisq_reltol;
    uint m_max_its;
    char m_verbosity;   // 'L' = low, 'M' = medium,  'H' = high

#ifndef NO_MINUIT
    static const char defaultmethod='M';
#else
    static const char defaultmethod='L';
#endif

 public:

    ChiSquareMinimizerInfo(XMLHandler& xmlin);
    ChiSquareMinimizerInfo(char method=defaultmethod, 
                           double param_reltol=1e-6, double chisq_reltol=1e-4, 
                           int max_its=1024, char verbosity='L');
    ChiSquareMinimizerInfo(const ChiSquareMinimizerInfo& info);
    ChiSquareMinimizerInfo& operator=(const ChiSquareMinimizerInfo& info);
    ~ChiSquareMinimizerInfo() {}

    void setMethod(char method);
    void setLMDer() {m_method='L';}
    void setNL2Sol()  {m_method='N';}
    void setMinuit2();
    void setChiSquareRelativeTolerance(double rtol);
    void setParameterRelativeTolerance(double rtol);
    void setMaximumIterations(unsigned int maxit);
    void setVerbosity(char verbosity);
    void setLowVerbosity() {m_verbosity='L';}
    void setMediumVerbosity() {m_verbosity='M';}
    void setHighVerbosity() {m_verbosity='H';}

    bool usingLMDer() const {return (m_method=='L');}
    bool usingNL2Sol() const {return (m_method=='N');}
    bool usingMinuit2() const {return (m_method=='M');}
    unsigned int getMaximumIterations() const {return m_max_its;}
    double getParameterRelativeTolerance() const {return m_param_reltol;}
    double getChiSquareRelativeTolerance() const {return m_chisq_reltol;}
    bool isLowVerbosity() const {return (m_verbosity=='L');}
    bool isMediumVerbosity() const {return (m_verbosity=='M');}
    bool isHighVerbosity() const {return (m_verbosity=='H');}
    bool isNotLowVerbosity() const {return (m_verbosity!='L');}

    void output(XMLHandler& xmlout) const;
    std::string output(int indent=0) const;  // XML output 
    std::string str() const;  // XML output

 friend class ChiSquareMinimizer;
   
};


// ******************************************************************************


class ChiSquareMinimizer
{

    ChiSquare *m_chisq;
    LMDerMinimizer *m_lmder;
    NL2SolMinimizer *m_nl2sol;
    ChiSquareMinimizerInfo m_info;

#ifndef NO_MINUIT
    Minuit2ChiSquare *m_minuit2;
#endif

#ifndef NO_CXX11
    ChiSquareMinimizer() = delete;
    ChiSquareMinimizer(const ChiSquareMinimizer&) = delete;
    ChiSquareMinimizer& operator=(const ChiSquareMinimizer&) = delete;
#else
    ChiSquareMinimizer();
    ChiSquareMinimizer(const ChiSquareMinimizer&);
    ChiSquareMinimizer& operator=(const ChiSquareMinimizer&);
#endif

 public:


    ChiSquareMinimizer(ChiSquare &in_chisq);
    ChiSquareMinimizer(ChiSquare &in_chisq, const ChiSquareMinimizerInfo& info);
    ~ChiSquareMinimizer();
    void reset(const ChiSquareMinimizerInfo& info);

    ChiSquareMinimizerInfo getInfo() const {return m_info;}
    ChiSquare& getChiSquare() {return *m_chisq;}

    bool findMinimum(const std::vector<double>& starting_params,
                     double& chisq_min, std::vector<double>& params_at_minimum,
                     XMLHandler& xmlout);
                 
    bool findMinimum(double& chisq_min, std::vector<double>& params_at_minimum,
                     XMLHandler& xmlout);

    bool findMinimum(const std::vector<double>& starting_params,
                     double& chisq_min, std::vector<double>& params_at_minimum);
                 
    bool findMinimum(double& chisq_min, std::vector<double>& params_at_minimum);

 private:

    void dealloc_method();
    void alloc_method();
    void set_method(char method);

    bool find_minimum(const std::vector<double>& starting_params,
                      double& chisq_min, std::vector<double>& params_at_minimum,
                      XMLHandler& xmlout, char verbosity);
    bool find_minimum_minuit2(const std::vector<double>& starting_params,
                              double& chisq_min, std::vector<double>& params_at_minimum,
                              XMLHandler& xmlout, char verbosity);
    bool find_minimum_lmder(const std::vector<double>& starting_params,
                            double& chisq_min, std::vector<double>& params_at_minimum,
                            XMLHandler& xmlout, char verbosity);
    bool find_minimum_nl2sol(const std::vector<double>& starting_params,
                             double& chisq_min, std::vector<double>& params_at_minimum,
                             XMLHandler& xmlout, char verbosity);

    void xmlformat(const std::string& roottag, const std::string& logstr,
                   XMLHandler& xmlout);
};


// ******************************************************************************

#ifndef NO_MINUIT
class Minuit2ChiSquare : public ROOT::Minuit2::FCNGradientBase
{

    ChiSquare *m_chisq;
    uint m_nobs, m_nparams;
    mutable std::vector<double> m_residuals;
    mutable RMatrix m_gradients;

    Minuit2ChiSquare(ChiSquare &in_chisq)
         : m_chisq(&in_chisq), m_nobs(m_chisq->getNumberOfObervables()),
           m_nparams(m_chisq->getNumberOfParams()),
           m_residuals(m_nobs), m_gradients(m_nobs,m_nparams) {}

    void guessInitialFitParamValues(std::vector<double>& params);

 public:

    double operator()(const std::vector<double>& params) const;

    double evalChiSquare(const std::vector<double>& params) const;

    std::vector<double> Gradient(const std::vector<double>& params) const;
    
    double Up() const {return 1.0;}

    friend class ChiSquareMinimizer;
};
#endif




class LMDerMinimizer
{

    ChiSquare *m_chisq;
    uint m_nobs, m_nparams;
    std::vector<double> m_fitparams, m_residuals;
    RMatrix m_gradients;
    std::vector<double> vdiag,vwa1,vwa2,vwa3,vqtf,vwa4;
    std::vector<int> vipvt;

    LMDerMinimizer(ChiSquare &in_chisq)
         : m_chisq(&in_chisq), m_nobs(m_chisq->getNumberOfObervables()),
           m_nparams(m_chisq->getNumberOfParams()),
           m_fitparams(m_nparams), m_residuals(m_nobs),
           m_gradients(m_nobs,m_nparams), vdiag(m_nparams), vwa1(m_nparams), 
           vwa2(m_nparams), vwa3(m_nparams), vqtf(m_nparams), vwa4(m_nobs),
           vipvt(m_nparams) {}

    void guessInitialFitParamValues();
    void setInitialFitParamValues(const std::vector<double>& start_params);

    int chisq_fit(double paramreltol, double chisqreltol, char verbosity, 
                  int max_its, double& chisq, std::ostringstream& outlog);

    const std::vector<double>& getFitParams() const {return m_fitparams;}

    friend class ChiSquareMinimizer;

};


class NL2SolMinimizer
{

    ChiSquare *m_chisq;
    uint m_nobs, m_nparams;
    std::vector<double> m_fitparams, m_residuals;
    RMatrix m_gradients;
    std::vector<double> v;
    std::vector<int> iv;
    

    NL2SolMinimizer(ChiSquare &in_chisq)
         : m_chisq(&in_chisq), m_nobs(m_chisq->getNumberOfObervables()),
           m_nparams(m_chisq->getNumberOfParams()),
           m_fitparams(m_nparams), m_residuals(m_nobs),
           m_gradients(m_nobs,m_nparams),
           v(94+2*m_nobs+(m_nparams*(3*m_nparams+33))/2),
           iv(60+m_nparams) {}

    void guessInitialFitParamValues();
    void setInitialFitParamValues(const std::vector<double>& start_params);

    int chisq_fit(double paramreltol, double chisqreltol, char verbosity, 
                  int max_its, double& chisq, std::ostringstream& outlog);

    const std::vector<double>& getFitParams() const {return m_fitparams;}

    friend class ChiSquareMinimizer;

};


// *************************************************************
#endif
