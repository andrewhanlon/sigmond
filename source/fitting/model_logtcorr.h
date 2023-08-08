#ifndef MODEL_LOGTCORR_H
#define MODEL_LOGTCORR_H
#include <vector>
#include "mcobs_info.h"
#include "model_tcorr.h"
#include "grace_plot.h"
#include "mc_estimate.h"

// ********************************************************************************
// *                                                                              *
// *   Real temporal correlators are key observables for extracting physical      *
// *   information in lattice QCD computations.  To extract the physics, a        *
// *   model fit function which describes the correlator is needed.  Here, an     *
// *   abstract base class for such a model is defined, and then different fit    *
// *   functions correspond to different classes derived from the base.           *
// *                                                                              *
// *   The important base class "LogTemporalCorrelatorModel" is defined in this   *
// *   file.  A variety of classes derived from "LogTemporalCorrelatorModel" are  *
// *   also defined.  Each class is used for a different model fit function.      *
// *   The constructor of the base class has the form                             *
// *                                                                              *
// *    LogTemporalCorrelatorModel(in_nparams,in_Tperiod,in_efftype);             *
// *                                                                              *
// *   This sets the number of parameters, the temporal extent of the lattice,    *
// *   and the integer code for the effective mass plot type.                     *
// *   The key member subroutines of the base class are:                          *
// *                                                                              *
// *     virtual void setupInfos(XMLHandler& xmlin,                               *
// *                             vector<MCObsInfo>& fitparam_info,                *
// *                             int taskcount);                                  *
// *                                                                              *
// *     virtual void evaluate(const ModelParams& p, double& value);              *
// *                                                                              *
// *     virtual void evalGradient(const ModelParams& p,                          *
// *                               std:vector<double>& grad);                     *
// *                                                                              *
// *     virtual void guessInitialParamValues(const std::vector<double>& data,    *
// *                                          ModelParams& p);                    *
// *                                                                              *
// *     virtual void output_tag(XMLHandler& xmlout) const = 0;                   *
// *                                                                              *
// *                                                                              *
// *     virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,    *
// *                 const std::vector<MCEstimate>& fitparams, uint fit_tmin,     *
// *                 uint fit_tmax, bool show_approach, uint meff_timestep,       *
// *                 double chisq_dof, double qual,                               *
// *                 TCorrFitInfo& fitinfo) const = 0;                            *
// *                                                                              *
// *                                                                              *
// *   The first member "setupInfos" sets up the MCObsInfo for the fit parameters.*
// *   The second member "evaluate" above evaluates the model function for a      *
// *   given set of parameter values, the third "evalGradient" evaluates the      *
// *   derivatives wrt these parameters at a given set of parameter values, and   *
// *   the fourth "guessInitialParamValues" is used to set initial values for     *
// *   the parameters given some known "data" when starting a fit.   The member   *
// *   "setFitInfo" sets up the data structure "TCorrFitInfo" which contains      *
// *   the necessary information to produce an effective energy plot with         *
// *   the fit information (see below).                                           *
// *                                                                              *
// *   Classes derived from "LogTemporalCorrelatorModel" below all have a         *
// *   constructor of the form                                                    *
// *                                                                              *
// *         DerivedClass(in_Tperiod);                                            *
// *                                                                              *
// *    The classes derived below from "LogTemporalCorrelatorModel" are           *
// *                                                                              *
// *        LogTimeForwardSingleExponential                                       *
// *        LogTimeForwardTwoExponential                                          *
// *                                                                              *
// *                                                                              *
// *   A useful routine for dynamically allocating an object of base class        *
// *  "LogTemporalCorrelatorModel" given a model type specified by a string is    *
// *   also defined in this file:                                                 *
// *                                                                              *
// *      void create_tcorr_model(const std::string& modeltype,                   *
// *                              uint in_Tperiod,                                *
// *                              LogTemporalCorrelatorModel* &mptr);             *
// *                                                                              *
// *   Polymorphism is used here.  An object of a derived class is actually       *
// *   created, but it can be accessed through the pointer to the base class.     *
// *                                                                              *
// ********************************************************************************



class LogTemporalCorrelatorModel
{

 protected:     // derived classes have access to the protected members

    uint m_nparams;  // number of fit parameters
    std::vector<std::string> param_names;
    uint T_period;   // temporal extent of lattice in number of sites
    uint m_effmasstype;   // effective mass type for plotting


 private:
          // disallow copying

#ifndef NO_CXX11
    LogTemporalCorrelatorModel() = delete;
    LogTemporalCorrelatorModel(const LogTemporalCorrelatorModel&) = delete;
    LogTemporalCorrelatorModel& operator=(const LogTemporalCorrelatorModel&) = delete;
#else
    LogTemporalCorrelatorModel();
    LogTemporalCorrelatorModel(const LogTemporalCorrelatorModel&);
    LogTemporalCorrelatorModel& operator=(const LogTemporalCorrelatorModel&);
#endif

 protected:

    LogTemporalCorrelatorModel(uint in_nparams, uint in_Tperiod, uint in_efftype) 
                : m_nparams(in_nparams), T_period(in_Tperiod), m_effmasstype(in_efftype) {}

 public:

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info,
                            int taskcount) const = 0;

    virtual void evaluate(const std::vector<double>& fitparams, double tval, 
                          double& value) const = 0;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const = 0;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals,
                                         std::vector<double>& fitparam) const = 0;    
    virtual void output_tag(XMLHandler& xmlout) const = 0;

    virtual ~LogTemporalCorrelatorModel(){}
                                         
    std::string getParameterName(uint param_index) const
     {return param_names[param_index];}

    uint getNumberOfParams() const
     {return m_nparams;}

    uint getEffMassType() const
     {return m_effmasstype;}


    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const = 0;


 protected:

    void simpleSetFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                          const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                          uint fit_tmax, double chisq_dof, double qual,
                          TCorrFitInfo& fitinfo) const;

    void approachSetFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, uint meff_step, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo, bool added_constant=false) const;
};

// *****************************************************************************

     //   A useful routine for dynamically allocating an object of base
     //   class "LogTemporalCorrelatorModel" given a model type specified by
     //   a string.  Polymorphism is used here.  An object of a derived
     //   class is actually created, but it can be accessed through the
     //   pointer to the base class.

void create_logtcorr_model(const std::string& modeltype, uint in_Tperiod,
                           LogTemporalCorrelatorModel* &mptr);





// *****************************************************************************
// *
// *
// *        Currently supported models:    names+ID index should be unique in each run
// *
// *                       log(A) - m*t
// *         <Model>
// *             <Type>LogTimeForwardSingleExponential</Type>          
// *             <Energy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *             <LogAmplitude>
// *                <Name>Amp</Name><IDIndex>0</IDIndex>
// *             </LogAmplitude>
// *         </Model>
// *
// *
// *                       log(A) - mt + log(1 + B * exp(-D^2*t) ) 
// *         <Model>
// *             <Type>LogTimeForwardTwoExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <LogFirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </LogFirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *         </Model>
// *

// ******************************************************************************


      // Fitting function is single exponential time-forward only:
      //
      //       f(t) = log(A) - m*t
      //
      // where 
      //             m = fitparams[0]
      //        log(A) = fitparams[1].
      //
      // For initial guess, need corr[tmin], corr[tmin+1]


class LogTimeForwardSingleExponential :  public LogTemporalCorrelatorModel 
{

#ifndef NO_CXX11
    LogTimeForwardSingleExponential() = delete;
    LogTimeForwardSingleExponential(const LogTimeForwardSingleExponential&) = delete;
    LogTimeForwardSingleExponential& operator=(const LogTimeForwardSingleExponential&) = delete;
#else
    LogTimeForwardSingleExponential();
    LogTimeForwardSingleExponential(const LogTimeForwardSingleExponential&);
    LogTimeForwardSingleExponential& operator=(const LogTimeForwardSingleExponential&);
#endif

 public:

    LogTimeForwardSingleExponential(uint in_Tperiod) 
          : LogTemporalCorrelatorModel(2,in_Tperiod,0) {}   // nparams = 2, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) const;

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual void output_tag(XMLHandler& xmlout) const;

    virtual ~LogTimeForwardSingleExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

    void eval_func(double logA, double m, double t, double& funcval) const;

    void eval_grad(double logA, double m, double t, double& dlogAval, double& dmval) const;

};


// ******************************************************************************


      // Fitting function is two exponential time-forward only:
      //
      //    f(t) = log(A)  - m*t + log[ 1 + B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //     log(A) = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //
      // For initial guess, need 4 corr values


class LogTimeForwardTwoExponential :  public LogTemporalCorrelatorModel 
{

#ifndef NO_CXX11
    LogTimeForwardTwoExponential() = delete;
    LogTimeForwardTwoExponential(const LogTimeForwardTwoExponential&) = delete;
    LogTimeForwardTwoExponential& operator=(const LogTimeForwardTwoExponential&) = delete;
#else
    LogTimeForwardTwoExponential();
    LogTimeForwardTwoExponential(const LogTimeForwardTwoExponential&);
    LogTimeForwardTwoExponential& operator=(const LogTimeForwardTwoExponential&);
#endif

 public:

    LogTimeForwardTwoExponential(uint in_Tperiod) 
          : LogTemporalCorrelatorModel(4,in_Tperiod,0) {}   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) const;

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual void output_tag(XMLHandler& xmlout) const;

    virtual ~LogTimeForwardTwoExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

    void eval_func(double logA, double m, double B, double DD,
                   double t, double& funcval) const;

    void eval_grad(double logA, double m, double B, double DD,
                   double t, double& dlogAval, double& dmval,
                   double& dBval, double& dDDval) const;

};


// ******************************************************************************
#endif
