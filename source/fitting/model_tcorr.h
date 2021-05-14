#ifndef MODEL_TCORR_H
#define MODEL_TCORR_H
#include <vector>
#include <algorithm>
#include "mcobs_info.h"
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
// *   The important base class "TemporalCorrelatorModel" is defined in this      *
// *   file.  A variety of classes derived from "TemporalCorrelatorModel" are     *
// *   also defined.  Each class is used for a different model fit function.      *
// *   The constructor of the base class has the form                             *
// *                                                                              *
// *    TemporalCorrelatorModel(in_nparams,in_Tperiod,in_efftype);                *
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
// *                                                                              *
// *   The first member "setupInfos" sets up the MCObsInfo for the fit parameters.*
// *   The second member "evaluate" above evaluates the model function for a      *
// *   given set of parameter values, the third "evalGradient" evaluates the      *
// *   derivatives wrt these parameters at a given set of parameter values, and   *
// *   the fourth "guessInitialParamValues" is used to set initial values for     *
// *   the parameters given some known "data" when starting a fit.                *
// *                                                                              *
// *   Classes derived from "TemporalCorrelatorModel" below all have a            *
// *   constructor of the form                                                    *
// *                                                                              *
// *         DerivedClass(in_Tperiod);                                            *
// *                                                                              *
// *    The classes derived below from "TemporalCorrelatorModel" are              *
// *                                                                              *
// *        TimeForwardSingleExponential                                          *
// *        TimeSymSingleExponential                                              *
// *        TimeForwardSingleExponentialPlusConstant                              *
// *        TimeSymSingleExponentialPlusConstant                                  *
// *        TimeForwardTwoExponential                                             *
// *        TimeSymTwoExponential                                                 *
// *        TimeForwardTwoExponentialPlusConstant                                 *
// *        TimeSymTwoExponentialPlusConstant                                     *
// *        TimeForwardThreeExponential                                           *
// *        TimeSymThreeExponential                                               *
// *        TimeForwardThreeExponentialPlusConstant                               *
// *        TimeSymThreeExponentialPlusConstant                                   *
// *        TimeForwardGeomSeriesExponential                                      *
// *        TimeSymGeomSeriesExponential                                          *
// *                                                                              *
// *                                                                              *
// *   A useful routine for dynamically allocating an object of base class        *
// *  "TemporalCorrelatorModel" given a model type specified by a string is       *
// *   also defined in this file:                                                 *
// *                                                                              *
// *      void create_tcorr_model(const std::string& modeltype,                   *
// *                              uint in_Tperiod,                                *
// *                              TemporalCorrelatorModel* &mptr);                *
// *                                                                              *
// *   Polymorphism is used here.  An object of a derived class is actually       *
// *   created, but it can be accessed through the pointer to the base class.     *
// *                                                                              *
// ********************************************************************************



class TemporalCorrelatorModel
{

 protected:     // derived classes have access to the protected members

    std::string model_name;
    std::vector<std::string> param_names;

    uint m_nparams;  // number of fit parameters
    uint T_period;   // temporal extent of lattice in number of sites
    uint m_effmasstype;   // effective mass type for plotting


 private:
          // disallow copying

    TemporalCorrelatorModel() = delete;
    TemporalCorrelatorModel(const TemporalCorrelatorModel&) = delete;
    TemporalCorrelatorModel& operator=(const TemporalCorrelatorModel&) = delete;

 protected:

    TemporalCorrelatorModel(uint in_nparams, uint in_Tperiod, uint in_efftype) 
                : m_nparams(in_nparams), T_period(in_Tperiod), m_effmasstype(in_efftype) {}

 public:

    void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) const;

    void setupInfos(std::map<std::string,MCObsInfo> model_params, std::vector<MCObsInfo>& fitparam_info) const;

    virtual void evaluate(const std::vector<double>& fitparams, double tval, 
                          double& value) const = 0;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const = 0;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals,
                                         std::vector<double>& fitparam) const = 0;    

    void output_tag(XMLHandler& xmlout) const;

    virtual ~TemporalCorrelatorModel(){}

    std::string getModelName() const
     {return model_name;}

    const std::vector<std::string>& getParameterNames() const
    {return param_names;}

    std::string getParameterName(uint param_index) const
     {return param_names[param_index];}

    uint getNumberOfParams() const
     {return m_nparams;}

    uint getEffMassType() const
     {return m_effmasstype;}

    uint getEnergyIndex(uint energy_level=0) const;

    uint getNumberOfEnergies() const;

    bool nonTrivialApproach() const
    {return (std::find(param_names.begin(), param_names.end(), "SqrtGapToSecondEnergy") != param_names.end());}

    bool addedConstant() const
    {return (std::find(param_names.begin(), param_names.end(), "AddedConstant") != param_names.end());}

    std::vector<XYPoint> getEffEnergyApproach(
        const std::vector<MCEstimate>& fitparams,
        uint fit_tmin,
        uint fit_tmax,
        uint meff_step) const;

};

// *****************************************************************************

     //   A useful routine for dynamically allocating an object of base
     //   class "TemporalCorrelatorModel" given a model type specified by
     //   a string.  Polymorphism is used here.  An object of a derived
     //   class is actually created, but it can be accessed through the
     //   pointer to the base class.

void create_tcorr_model(const std::string& modeltype, uint in_Tperiod,
                        TemporalCorrelatorModel* &mptr);




// *****************************************************************************
// *
// *
// *        Currently supported models:    names+ID index should be unique in each run
// *
// *                       A * exp(-m*t)
// *         <Model>
// *             <Type>TimeForwardSingleExponential</Type>          
// *             <Energy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *             <Amplitude>
// *                <Name>Amp</Name><IDIndex>0</IDIndex>
// *             </Amplitude>
// *         </Model>
// *
// *                       A * (exp(-m*t) + exp(-m*(T-t)) )
// *         <Model>
// *             <Type>TimeSymSingleExponential</Type>             
// *             <Energy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *             <Amplitude>
// *                <Name>Amp</Name><IDIndex>0</IDIndex>
// *             </Amplitude>
// *         </Model>
// *
// *                       A * exp(-m*t) + C0
// *         <Model>
// *             <Type>TimeForwardSingleExponentialPlusConstant</Type>          
// *             <Energy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *             <Amplitude>
// *                <Name>Amp</Name><IDIndex>0</IDIndex>
// *             </Amplitude>
// *             <AddedConstant>
// *                <Name>Cnst</Name><IDIndex>0</IDIndex>
// *             </AddedConstant>
// *         </Model>
// *
// *                       A * (exp(-m*t) + exp(-m*(T-t)) ) + C0
// *         <Model>
// *             <Type>TimeSymSingleExponentialPlusConstant</Type>  
// *             <Energy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *             <Amplitude>
// *                <Name>Amp</Name><IDIndex>0</IDIndex>
// *             </Amplitude>
// *             <AddedConstant>
// *                <Name>Cnst</Name><IDIndex>0</IDIndex>
// *             </AddedConstant>
// *         </Model>
// *
// *                       A * exp(-m*t)*(1 + B * exp(-D^2*t) ) 
// *         <Model>
// *             <Type>TimeForwardTwoExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *         </Model>
// *
// *                       A * (exp(-m*t)*(1 + B * exp(-D^2*t) )  
// *                          + exp(-m*(T-t))*(1 + B * exp(-D^2*(T-t)) ) )
// *         <Model>
// *             <Type>TimeSymTwoExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *         </Model>
// *
// *                       A * (exp(-m*t)*(1 + B * exp(-D^2*t) )  + C0
// *         <Model>
// *             <Type>TimeForwardTwoExponentialPlusConstant</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *             <AddedConstant>
// *                <Name>Cnst</Name><IDIndex>0</IDIndex>
// *             </AddedConstant>
// *         </Model>
// *
// *                       A * (exp(-m*t)*(1 + B * exp(-D^2*t) )  
// *                          + exp(-m*(T-t))*(1 + B * exp(-D^2*(T-t)) ) ) + C0
// *         <Model>
// *             <Type>TimeSymTwoExponentialPlusConstant</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *             <AddedConstant>
// *                <Name>Cnst</Name><IDIndex>0</IDIndex>
// *             </AddedConstant>
// *         </Model>
// *
// *                       A * exp(-m*t)*(1 + B * exp(-D^2*t) + F * exp(-G^2*t) ) 
// *         <Model>
// *             <Type>TimeForwardThreeExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *             <SqrtGapToThirdEnergy>
// *                <Name>pionprimeprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToThirdEnergy>
// *             <ThirdAmplitudeRatio>
// *                <Name>Amp2</Name><IDIndex>0</IDIndex>
// *             </ThirdAmplitudeRatio>
// *         </Model>
// *
// *                       A * (exp(-m*t)*(1 + B * exp(-D^2*t) + F * exp(-G^2*t) )  
// *                          + exp(-m*(T-t))*(1 + B * exp(-D^2*(T-t)) + F * exp(-G^2*(T-t)) ) )
// *         <Model>
// *             <Type>TimeSymThreeExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *             <SqrtGapToThirdEnergy>
// *                <Name>pionprimeprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToThirdEnergy>
// *             <ThirdAmplitudeRatio>
// *                <Name>Amp2</Name><IDIndex>0</IDIndex>
// *             </ThirdAmplitudeRatio>
// *         </Model>
// *
// *                       A * (exp(-m*t)*(1 + B * exp(-D^2*t) + F * exp(-G^2*t) )  + C0
// *         <Model>
// *             <Type>TimeForwardThreeExponentialPlusConstant</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *             <SqrtGapToThirdEnergy>
// *                <Name>pionprimeprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToThirdEnergy>
// *             <ThirdAmplitudeRatio>
// *                <Name>Amp2</Name><IDIndex>0</IDIndex>
// *             </ThirdAmplitudeRatio>
// *             <AddedConstant>
// *                <Name>Cnst</Name><IDIndex>0</IDIndex>
// *             </AddedConstant>
// *         </Model>
// *
// *                       A * (exp(-m*t)*(1 + B * exp(-D^2*t) + F * exp(-G^2*t) )  
// *                          + exp(-m*(T-t))*(1 + B * exp(-D^2*(T-t)) + F * exp(-G^2*t) ) ) + C0
// *         <Model>
// *             <Type>TimeSymThreeExponentialPlusConstant</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *             <SqrtGapToThirdEnergy>
// *                <Name>pionprimeprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToThirdEnergy>
// *             <ThirdAmplitudeRatio>
// *                <Name>Amp2</Name><IDIndex>0</IDIndex>
// *             </ThirdAmplitudeRatio>
// *             <AddedConstant>
// *                <Name>Cnst</Name><IDIndex>0</IDIndex>
// *             </AddedConstant>
// *         </Model>
// *
// *
// *                       A * exp(-m*t)/(1 - B * exp(-D^2*t) )
// *         <Model>
// *             <Type>TimeForwardGeomSeriesExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
// *             <SqrtGapToSecondEnergy>
// *                <Name>pionprime</Name><IDIndex>0</IDIndex>
// *             </SqrtGapToSecondEnergy>
// *             <SecondAmplitudeRatio>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </SecondAmplitudeRatio>
// *         </Model>
// *
// *                       A * (exp(-m*t)/(1 - B * exp(-D^2*t) )  
// *                          + exp(-m*(T-t))/(1 - B * exp(-D^2*(T-t)) ) )
// *         <Model>
// *             <Type>TimeSymGeomSeriesExponential</Type>
// *             <FirstEnergy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </FirstEnergy>
// *             <FirstAmplitude>
// *                <Name>Amp0</Name><IDIndex>0</IDIndex>
// *             </FirstAmplitude>
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
      //       f(t) = A * exp( -m*t ) 
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].
      //
      // For initial guess, need corr[tmin], corr[tmin+1]


class TimeForwardSingleExponential :  public TemporalCorrelatorModel 
{

    TimeForwardSingleExponential() = delete;
    TimeForwardSingleExponential(const TimeForwardSingleExponential&) = delete;
    TimeForwardSingleExponential& operator=(const TimeForwardSingleExponential&) = delete;

 public:

    TimeForwardSingleExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(2,in_Tperiod,0)    // nparams = 2, efftype = 0
    {
        model_name = "TimeForwardSingleExponential";
        param_names = {
            "Energy",
            "Amplitude"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardSingleExponential(){}

 private:

    void eval_func(double A, double m, double t, double& funcval) const;

    void eval_grad(double A, double m, double t, double& dAval, double& dmval) const;

    static void get_exp_guess(const std::vector<uint>& tvals, 
                              const std::vector<double>& corrvals,
                              double& energy0, double& amp0);

/*  static void get_exp_guess(int tval, double corrt, int tnext, double corrtnext, 
                              double& A, double& m); */

    friend class TimeSymSingleExponential;
    friend class TimeForwardTwoExponential;
    friend class TimeForwardSingleExponentialPlusConstant;
    friend class TimeSymTwoExponential;
    friend class TimeForwardGeomSeriesExponential;
    friend class TimeSymGeomSeriesExponential;
    friend class TimeForwardTwoExponentialPlusConstant;
    friend class TimeForwardThreeExponential;
    friend class TimeForwardThreeExponentialPlusConstant;

};


// ******************************************************************************


      // Fitting function is single exponential time-symmetric
      // (forwards and backwards):
      //
      //       f(t) = A * {exp( -m*t ) +  exp( -m*(T_period-t) )}
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].
      //
      // For initial guess, need corr[tmin], corr[tmin+1]


class TimeSymSingleExponential :  public TemporalCorrelatorModel 
{

    TimeSymSingleExponential() = delete;
    TimeSymSingleExponential(const TimeSymSingleExponential&) = delete;
    TimeSymSingleExponential& operator=(const TimeSymSingleExponential&) = delete;

 public:

    TimeSymSingleExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(2,in_Tperiod,1)    // nparams = 2, efftype = 1
    {
        model_name = "TimeSymSingleExponential";
        param_names = {
            "Energy",
            "Amplitude"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymSingleExponential(){}

 private:

    void eval_func(double A, double m, double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double t, int Nt, double& dAval, double& dmval) const;

};


// ******************************************************************************


      // The fitting function is a single exponential time-forward
      // with an added constant:
      //
      //    f(t) = A * exp(-m*t)  + c0
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1]
      //          c0 = fitparams[2]
      //
      // For initial guess, need corr[tmin], corr[tmin+1], corr[tmin+2]


class TimeForwardSingleExponentialPlusConstant :  public TemporalCorrelatorModel 
{

    TimeForwardSingleExponentialPlusConstant() = delete;
    TimeForwardSingleExponentialPlusConstant(const TimeForwardSingleExponentialPlusConstant&) = delete;
    TimeForwardSingleExponentialPlusConstant& operator=(const TimeForwardSingleExponentialPlusConstant&) = delete;

 public:

    TimeForwardSingleExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(3,in_Tperiod,2)    // nparams = 3, efftype = 2
    {
        model_name = "TimeForwardSingleExponentialPlusConstant";
        param_names = {
            "Energy",
            "Amplitude",
            "AddedConstant"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardSingleExponentialPlusConstant(){}

 private:

    void eval_func(double A, double m, double c0, double t, double& funcval) const;

    void eval_grad(double A, double m, double t, double& dAval, double& dmval,
                   double& dc0val) const;

    static void take_diff(const std::vector<uint>& tvals, 
               const std::vector<double>& corrvals,
               std::vector<uint>& tdf, std::vector<double>& corrdiff, uint& step);

    static void get_exp_plus_const_guess(
               const std::vector<uint>& tvals, const std::vector<double>& corrvals,
               double& energy0, double& amp0, double& c0);

/*  static void get_exp_plus_const_guess(
                int tval, double corrt, int tp1, double corrtp1, int tp2, double corrtp2,
                double& A, double& m, double& c0); */

    friend class TimeSymSingleExponentialPlusConstant;
    friend class TimeForwardTwoExponentialPlusConstant;
    friend class TimeSymTwoExponentialPlusConstant;
};


// ******************************************************************************


      // The fitting function is a single exponential time-symmetric
      // (forwards and backwards) with an added constant:
      //
      //       f(t) = A * {exp( -m*t ) +  exp( -m*(T_period-t) )} + c0
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1]
      //          c0 = fitparams[2]
      //
      // For initial guess, need corr[tmin], corr[tmin+1], corr[tmin+2]


class TimeSymSingleExponentialPlusConstant :  public TemporalCorrelatorModel 
{

    TimeSymSingleExponentialPlusConstant() = delete;
    TimeSymSingleExponentialPlusConstant(const TimeSymSingleExponentialPlusConstant&) = delete;
    TimeSymSingleExponentialPlusConstant& operator=(const TimeSymSingleExponentialPlusConstant&) = delete;

 public:

    TimeSymSingleExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(3,in_Tperiod,3)    // nparams = 3, efftype = 3
    {
        model_name = "TimeSymSingleExponentialPlusConstant";
        param_names = {
            "Energy",
            "Amplitude",
            "AddedConstant"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymSingleExponentialPlusConstant(){}

 private:

    void eval_func(double A, double m, double c0, double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double t, int Nt, double& dAval, double& dmval,
                   double& dc0val) const;


};


// ******************************************************************************


      // Fitting function is two exponential time-forward only:
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //
      // For initial guess, need 4 corr values


class TimeForwardTwoExponential :  public TemporalCorrelatorModel 
{

    TimeForwardTwoExponential() = delete;
    TimeForwardTwoExponential(const TimeForwardTwoExponential&) = delete;
    TimeForwardTwoExponential& operator=(const TimeForwardTwoExponential&) = delete;

 public:

    TimeForwardTwoExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,0)    // nparams = 4, efftype = 0
    {
        model_name = "TimeForwardTwoExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardTwoExponential(){}

 private:

    void eval_func(double A, double m, double B, double DD,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD,
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval) const;

    static void get_two_exp_guess(const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                                  double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                                  double tasymfrac=0.33);

/*  static void get_two_exp_guess(
              int tfar, double ffar, int tfarnext, double ffarnext, 
              int tnear, double fnear, int tnearnext, double fnearnext,
              double& A, double& m, double& B, double& DD);
    static void get_two_exp_guess(
                     const std::vector<double>& data, const std::vector<uint>& tvals,
                     std::vector<double>& fitparams); */

    friend class TimeSymTwoExponential;
    friend class TimeForwardGeomSeriesExponential;
    friend class TimeSymGeomSeriesExponential;
    friend class TimeForwardThreeExponential;
};


// ******************************************************************************



      // The fitting function is a sum of two exponentials, time-symmetric:
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //          + exp(-m*(T_period-t)) * [ 1 + B*exp(-Delta^2*(T_period-t)) ] }
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //
      // For initial guess, need 4 corr values


class TimeSymTwoExponential :  public TemporalCorrelatorModel 
{

    TimeSymTwoExponential() = delete;
    TimeSymTwoExponential(const TimeSymTwoExponential&) = delete;
    TimeSymTwoExponential& operator=(const TimeSymTwoExponential&) = delete;

 public:

    TimeSymTwoExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,1)    // nparams = 4, efftype = 1
    {
        model_name = "TimeSymTwoExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymTwoExponential(){}

 private:

    void eval_func(double A, double m, double B, double DD,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD,
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval) const;

};


// ******************************************************************************


      // The fitting function is a sum of two exponentials, time-forward + constant
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ] + c0
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         c0 = fitparams[4]
      //
      // For initial guess, need 5 corr values


class TimeForwardTwoExponentialPlusConstant :  public TemporalCorrelatorModel 
{

    TimeForwardTwoExponentialPlusConstant() = delete;
    TimeForwardTwoExponentialPlusConstant(const TimeForwardTwoExponentialPlusConstant&) = delete;
    TimeForwardTwoExponentialPlusConstant& operator=(const TimeForwardTwoExponentialPlusConstant&) = delete;

 public:

    TimeForwardTwoExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(5,in_Tperiod,2)    // nparams = 5, efftype = 2
    {
        model_name = "TimeForwardTwoExponentialPlusConstant";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "AddedConstant"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardTwoExponentialPlusConstant(){}

 private:

    void eval_func(double A, double m, double B, double DD, double c0,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD,
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dc0val) const;

    static void get_two_exp_plus_const_guess(
                       const std::vector<uint>& tvals, 
                       const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double& c0, double tasymfrac=0.33);

/*  static void get_two_exp_plus_const_guess(
              int tfar, double ffar, int tfarp1, double ffarp1, int tfarp2, double ffarp2,
              int tnear, double fnear, int tnearnext, double fnearnext,
              double& A, double& m, double& B, double& DD, double& c0);

    static void get_two_exp_plus_const_guess(
                        const std::vector<double>& data, const std::vector<uint>& tvals,
                        std::vector<double>& fitparams);  */

    friend class TimeSymTwoExponentialPlusConstant;
    friend class TimeForwardThreeExponentialPlusConstant;
};


// ******************************************************************************


      // The fitting function is a sum of two exponentials, time-symmetric + constant
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //          + exp(-m*(T_period-t)) * [ 1 + B*exp(-Delta^2*(T_period-t)) ] }  + c0
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         c0 = fitparams[4]
      //
      // For initial guess, need 5 corr values


class TimeSymTwoExponentialPlusConstant :  public TemporalCorrelatorModel 
{

    TimeSymTwoExponentialPlusConstant() = delete;
    TimeSymTwoExponentialPlusConstant(const TimeSymTwoExponentialPlusConstant&) = delete;
    TimeSymTwoExponentialPlusConstant& operator=(const TimeSymTwoExponentialPlusConstant&) = delete;

 public:

    TimeSymTwoExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(5,in_Tperiod,3)    // nparams = 5, efftype = 3
    {
        model_name = "TimeSymTwoExponentialPlusConstant";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "AddedConstant"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymTwoExponentialPlusConstant(){}

 private:

    void eval_func(double A, double m, double B, double DD, double c0,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, 
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dc0val) const;

};

// ******************************************************************************


      // Fitting function is three exponential time-forward only:
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B*exp(-Delta^2*t) + F*exp(-GG^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         GG = fitparams[4]
      //          F = fitparams[5]
      //
      // For initial guess, need 6 corr values


class TimeForwardThreeExponential :  public TemporalCorrelatorModel 
{

    TimeForwardThreeExponential() = delete;
    TimeForwardThreeExponential(const TimeForwardThreeExponential&) = delete;
    TimeForwardThreeExponential& operator=(const TimeForwardThreeExponential&) = delete;

 public:

    TimeForwardThreeExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(6,in_Tperiod,0)    // nparams = 6, efftype = 0
    {
        model_name = "TimeForwardThreeExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardThreeExponential(){}

 private:

    void eval_func(double A, double m, double B, double DD, double F, double GG,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double F, double GG,
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dFval, double& dGGval) const;

    static void get_three_exp_guess(const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                                  double& energy0, double& amp0, double& gapsqrt1, double& gapamp1,
                                  double& gapsqrt2, double& gapamp2, double tasymfrac=0.33);


    friend class TimeSymThreeExponential;
};


// ******************************************************************************



      // The fitting function is a sum of three exponentials, time-symmetric:
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) + F*exp(-GG^2*t) ]
      //          + exp(-m*(T_period-t)) * [ 1 + B*exp(-Delta^2*(T_period-t)) + F*exp(-GG^2*(T_period-t)) ] }
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         GG = fitparams[4]
      //          F = fitparams[5]
      //
      // For initial guess, need 6 corr values


class TimeSymThreeExponential :  public TemporalCorrelatorModel 
{

    TimeSymThreeExponential() = delete;
    TimeSymThreeExponential(const TimeSymThreeExponential&) = delete;
    TimeSymThreeExponential& operator=(const TimeSymThreeExponential&) = delete;

 public:

    TimeSymThreeExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(6,in_Tperiod,1)    // nparams = 6, efftype = 1
    {
        model_name = "TimeSymThreeExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymThreeExponential(){}

 private:

    void eval_func(double A, double m, double B, double DD, double F, double GG,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double F, double GG,
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dFval, double& dGGval) const;

};


// ******************************************************************************


      // The fitting function is a sum of three exponentials, time-forward + constant
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B*exp(-Delta^2*t) + F*exp(-GG^2*t) ] + c0
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         GG = fitparams[4]
      //          F = fitparams[5]
      //         c0 = fitparams[6]
      //
      // For initial guess, need 7 corr values


class TimeForwardThreeExponentialPlusConstant :  public TemporalCorrelatorModel 
{

    TimeForwardThreeExponentialPlusConstant() = delete;
    TimeForwardThreeExponentialPlusConstant(const TimeForwardThreeExponentialPlusConstant&) = delete;
    TimeForwardThreeExponentialPlusConstant& operator=(const TimeForwardThreeExponentialPlusConstant&) = delete;

 public:

    TimeForwardThreeExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(7,in_Tperiod,2)    // nparams = 7, efftype = 2
    {
        model_name = "TimeForwardThreeExponentialPlusConstant";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio",
            "AddedConstant"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardThreeExponentialPlusConstant(){}

 private:

    void eval_func(double A, double m, double B, double DD, double F, double GG, double c0,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double F, double GG,
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dFval, double& dGGval,
                   double& dc0val) const;

    static void get_three_exp_plus_const_guess(
                       const std::vector<uint>& tvals, 
                       const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt1, double& gapamp1,
                       double& gapsqrt2, double& gapamp2, double& c0, double tasymfrac=0.33);

    friend class TimeSymThreeExponentialPlusConstant;
};


// ******************************************************************************


      // The fitting function is a sum of three exponentials, time-symmetric + constant
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) + F*exp(-GG^2*t) ]
      //          + exp(-m*(T_period-t)) * [ 1 + B*exp(-Delta^2*(T_period-t)) + F*exp(-GG^2*(T_period-t)) ] }  + c0
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         GG = fitparams[4]
      //          F = fitparams[5]
      //         c0 = fitparams[6]
      //
      // For initial guess, need 7 corr values


class TimeSymThreeExponentialPlusConstant :  public TemporalCorrelatorModel 
{

    TimeSymThreeExponentialPlusConstant() = delete;
    TimeSymThreeExponentialPlusConstant(const TimeSymThreeExponentialPlusConstant&) = delete;
    TimeSymThreeExponentialPlusConstant& operator=(const TimeSymThreeExponentialPlusConstant&) = delete;

 public:

    TimeSymThreeExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(7,in_Tperiod,3)    // nparams = 7, efftype = 3
    {
        model_name = "TimeSymThreeExponentialPlusConstant";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio",
            "AddedConstant"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymThreeExponentialPlusConstant(){}

 private:

    void eval_func(double A, double m, double B, double DD, double F, double GG, double c0,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double F, double GG,
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dFval, double& dGGval,
                   double& dc0val) const;

};

// ******************************************************************************


      // Fitting function is sum of a geometric series of exponentials time-forward only:
      //
      //    f(t) = A * exp(-m*t) / [ 1 - B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //
      // For initial guess, need 4 corr values


class TimeForwardGeomSeriesExponential :  public TemporalCorrelatorModel 
{

    TimeForwardGeomSeriesExponential() = delete;
    TimeForwardGeomSeriesExponential(const TimeForwardGeomSeriesExponential&) = delete;
    TimeForwardGeomSeriesExponential& operator=(const TimeForwardGeomSeriesExponential&) = delete;

 public:

    TimeForwardGeomSeriesExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,0)    // nparams = 4, efftype = 0
    {
        model_name = "TimeForwardGeomSeriesExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardGeomSeriesExponential(){}

 private:

    void eval_func(double A, double m, double B, double DD,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD,
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval) const;

    friend class TimeSymGeomSeriesExponential;

};


// ******************************************************************************



      // The fitting function is a sum of a geometric series of exponentials, time-symmetric:
      //
      //    f(t) = A * { exp(-m*t) / [ 1 - B*exp(-Delta^2*t) ]
      //          + exp(-m*(T_period-t)) / [ 1 - B*exp(-Delta^2*(T_period-t)) ] }
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //
      // For initial guess, need 4 corr values


class TimeSymGeomSeriesExponential :  public TemporalCorrelatorModel 
{

    TimeSymGeomSeriesExponential() = delete;
    TimeSymGeomSeriesExponential(const TimeSymGeomSeriesExponential&) = delete;
    TimeSymGeomSeriesExponential& operator=(const TimeSymGeomSeriesExponential&) = delete;

 public:

    TimeSymGeomSeriesExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,1)    // nparams = 4, efftype = 1
    {
        model_name = "TimeSymGeomSeriesExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymGeomSeriesExponential(){}

 private:

    void eval_func(double A, double m, double B, double DD,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD,
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval) const;

};


// ******************************************************************************
#endif
