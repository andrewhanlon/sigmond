#ifndef MODEL_TCORR_H
#define MODEL_TCORR_H
#include <vector>
#include "mcobs_info.h"
#include "grace_plot.h"
#include "mc_estimate.h"
#include "prior.h"

class TCorrFitInfo;

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
// *        TimeForwardGeomSeriesExponential                                      *
// *        TimeSymGeomSeriesExponential                                          *
// *        TimeForwardGMO                                                        *
// *        TimeForwardTwoIndExp                                                  *
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
    std::map<uint,Prior> m_priors;
    std::vector<uint> m_prior_type;


 private:
          // disallow copying

#ifndef NO_CXX11
    TemporalCorrelatorModel() = delete;
    TemporalCorrelatorModel(const TemporalCorrelatorModel&) = delete;
    TemporalCorrelatorModel& operator=(const TemporalCorrelatorModel&) = delete;
#else
    TemporalCorrelatorModel();
    TemporalCorrelatorModel(const TemporalCorrelatorModel&);
    TemporalCorrelatorModel& operator=(const TemporalCorrelatorModel&);
#endif

 protected:

    TemporalCorrelatorModel(uint in_nparams, uint in_Tperiod, uint in_efftype) 
                : m_nparams(in_nparams), T_period(in_Tperiod), m_effmasstype(in_efftype) {
                    m_prior_type.resize(m_nparams);
                    fill(m_prior_type.begin(), m_prior_type.end(), 0);
                }

 public:

    void simpleSetupInfo(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) const;

    void setupInfos(std::map<std::string,MCObsInfo> model_params, std::vector<MCObsInfo>& fitparam_info) const;

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) = 0;
    // virtual void setupInfos(std::map<std::string,MCObsInfo> model_params, 
    //                             std::vector<MCObsInfo>& fitparam_info) const = 0;

    virtual void evaluate(const std::vector<double>& fitparams, double tval, 
                          double& value) const = 0;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const = 0;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals,
                                         std::vector<double>& fitparam) const = 0;    

    void output_tag(XMLHandler& xmlout) const;

    virtual ~TemporalCorrelatorModel(){}

    std::string getParameterName(uint param_index) const
     {return param_names[param_index];}

    uint getNumberOfParams() const
     {return m_nparams;}

    uint getEffMassType() const
     {return m_effmasstype;}

    void setupPriors( std::map<uint,Prior> in_priors){
        m_priors=in_priors;
        for(std::map<uint,Prior>::iterator prior_it=m_priors.begin(); prior_it!=m_priors.end(); ++prior_it) {
            if(m_prior_type[prior_it->first]>0) prior_it->second.setType(m_prior_type[prior_it->first]);
        }
    }

    void initializeParametersWithPriors( std::vector<double>& fitparams ) const;

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const = 0;

    double eval(const std::vector<double>& fitparams, double tval){
        double value;
        evaluate(fitparams, tval, value);
        return value;
    }

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
     //   class "TemporalCorrelatorModel" given a model type specified by
     //   a string.  Polymorphism is used here.  An object of a derived
     //   class is actually created, but it can be accessed through the
     //   pointer to the base class.

void create_tcorr_model(const std::string& modeltype, uint in_Tperiod,
                        TemporalCorrelatorModel* &mptr);




// *****************************************************************************

     //   This struct is used for plotting fit results on an effective
     //   energy plot.   For fits more complicated then single exponential,
     //   such as two-exponential fits, "meff_approach" contains points
     //   to plot that show how the asymptotic energies are reached.
     //   It is also used by the pivoting classes in computing Z-factors.
     //   The "energy_key" and "amplitude_key" are needed for extracting
     //   the resamplings of the energy and amplitude when needed.


class TCorrFitInfo
{
 public:
    uint tmin, tmax;
    double energy_mean, energy_err;
    MCObsInfo energy_key, amplitude_key;
    std::vector<XYPoint> meff_approach;
    double chisq_dof, quality;
    //prior param deviations?
};




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
// *                 AL * AS^(1/3) exp( -( mL+mS/3-2*mN/3-2*mX/3 )*t ) 
// *                                    /  AN^(2/3) / AX^(2/3)
// *         <Model>
// *             <Type>TimeForwardGMO</Type> (Only created for a simultaneous fit)
// *             <LambdaEnergy><Name>EnergyL</Name><IDIndex>0</IDIndex></LambdaEnergy>
// *             <LambdaAmplitude><Name>AmpL</Name><IDIndex>0</IDIndex></LambdaAmplitude>
// *             <SigmaEnergy><Name>EnergyS</Name><IDIndex>0</IDIndex></SigmaEnergy>
// *             <SigmaAmplitude><Name>AmpS</Name><IDIndex>0</IDIndex></SigmaAmplitude>
// *             <NucleonEnergy><Name>EnergyN</Name><IDIndex>0</IDIndex></NucleonEnergy>
// *             <NucleonAmplitude><Name>AmpN</Name><IDIndex>0</IDIndex></NucleonAmplitude>
// *             <XiEnergy><Name>EnergyX</Name><IDIndex>0</IDIndex></XiEnergy>
// *             <XiAmplitude><Name>AmpX</Name><IDIndex>0</IDIndex></XiAmplitude>
// *             <InitialEnergies>
// *                 <Lambda>0.363420379631</Lambda> (Good initial guesses of the fit parameters.
// *                 <Sigma>0.383032259747</Sigma>     Recommend using individual fit results.)
// *                 <Nucleon>0.314290654908</Nucleon>
// *                 <Xi>0.415428561009</Xi>
// *             </InitialEnergies>
// *             <InitialAmplitudes>
// *                 <Lambda>0.0499660266931</Lambda>
// *                 <Sigma>0.0730132545349</Sigma>
// *                 <Nucleon>0.0543509286276</Nucleon>
// *                 <Xi>0.0919270433468</Xi>
// *             </InitialAmplitudes>
// *         </Model>
// *
// *
// *                       A * exp(-m*t) + A1 * exp(-m1*t)
// *         <Model>
// *             <Type>TimeForwardTwoIndExp</Type>          
// *             <Energy>
// *                <Name>pion</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *             <Amplitude>
// *                <Name>Amp</Name><IDIndex>0</IDIndex>
// *             </Amplitude>
// *             <Energy1>
// *                <Name>pion1</Name><IDIndex>0</IDIndex>
// *             </Energy1>
// *             <Amplitude1>
// *                <Name>Amp1</Name><IDIndex>0</IDIndex>
// *             </Amplitude1>
// *         </Model>

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

#ifndef NO_CXX11
    TimeForwardSingleExponential() = delete;
    TimeForwardSingleExponential(const TimeForwardSingleExponential&) = delete;
    TimeForwardSingleExponential& operator=(const TimeForwardSingleExponential&) = delete;
#else
    TimeForwardSingleExponential();
    TimeForwardSingleExponential(const TimeForwardSingleExponential&);
    TimeForwardSingleExponential& operator=(const TimeForwardSingleExponential&);
#endif

 public:

    TimeForwardSingleExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(2,in_Tperiod,0) {
	    model_name = "TimeForwardSingleExponential";
        param_names = {
            "Energy",
            "Amplitude"
        };	  
	  
    }   // nparams = 2, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardSingleExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

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
    friend class LogTimeForwardSingleExponential;
    friend class TimeForwardTwoIndExp;

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

#ifndef NO_CXX11
    TimeSymSingleExponential() = delete;
    TimeSymSingleExponential(const TimeSymSingleExponential&) = delete;
    TimeSymSingleExponential& operator=(const TimeSymSingleExponential&) = delete;
#else
    TimeSymSingleExponential();
    TimeSymSingleExponential(const TimeSymSingleExponential&);
    TimeSymSingleExponential& operator=(const TimeSymSingleExponential&);
#endif

 public:

    TimeSymSingleExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(2,in_Tperiod,1) {
            model_name = "TimeSymSingleExponential";
          }   // nparams = 2, efftype = 1

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymSingleExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

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

#ifndef NO_CXX11
    TimeForwardSingleExponentialPlusConstant() = delete;
    TimeForwardSingleExponentialPlusConstant(const TimeForwardSingleExponentialPlusConstant&) = delete;
    TimeForwardSingleExponentialPlusConstant& operator=(const TimeForwardSingleExponentialPlusConstant&) = delete;
#else
    TimeForwardSingleExponentialPlusConstant();
    TimeForwardSingleExponentialPlusConstant(const TimeForwardSingleExponentialPlusConstant&);
    TimeForwardSingleExponentialPlusConstant& operator=(const TimeForwardSingleExponentialPlusConstant&);
#endif

 public:

    TimeForwardSingleExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(3,in_Tperiod,2) {
            model_name = "TimeForwardSingleExponentialPlusConstant";
          }   // nparams = 3, efftype = 2

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardSingleExponentialPlusConstant(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

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

#ifndef NO_CXX11
    TimeSymSingleExponentialPlusConstant() = delete;
    TimeSymSingleExponentialPlusConstant(const TimeSymSingleExponentialPlusConstant&) = delete;
    TimeSymSingleExponentialPlusConstant& operator=(const TimeSymSingleExponentialPlusConstant&) = delete;
#else
    TimeSymSingleExponentialPlusConstant();
    TimeSymSingleExponentialPlusConstant(const TimeSymSingleExponentialPlusConstant&);
    TimeSymSingleExponentialPlusConstant& operator=(const TimeSymSingleExponentialPlusConstant&);
#endif

 public:

    TimeSymSingleExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(3,in_Tperiod,3) {
            model_name = "TimeSymSingleExponentialPlusConstant";
          }   // nparams = 3, efftype = 3

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymSingleExponentialPlusConstant(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

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

#ifndef NO_CXX11
    TimeForwardTwoExponential() = delete;
    TimeForwardTwoExponential(const TimeForwardTwoExponential&) = delete;
    TimeForwardTwoExponential& operator=(const TimeForwardTwoExponential&) = delete;
#else
    TimeForwardTwoExponential();
    TimeForwardTwoExponential(const TimeForwardTwoExponential&);
    TimeForwardTwoExponential& operator=(const TimeForwardTwoExponential&);
#endif

 public:

    TimeForwardTwoExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,0) {
        model_name = "TimeForwardTwoExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio"
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardTwoExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

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
    friend class LogTimeForwardTwoExponential;
    friend class TimeForwardDoubleExpRatio1;
    friend class TimeForwardTwoExponentialForCons;
    friend class TimeForwardDoubleExpRatio2;
    friend class TimeForwardTwoIndExp;
    friend class TimeForwardThreeExponential;
    friend class TwoExpConspiracy;
    friend class DegTwoExpConspiracy;
    friend class DegThreeExpConspiracy;
    friend class TimeForwardThreeIndExponential;
    friend class TimeForwardFourExponential;
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

#ifndef NO_CXX11
    TimeSymTwoExponential() = delete;
    TimeSymTwoExponential(const TimeSymTwoExponential&) = delete;
    TimeSymTwoExponential& operator=(const TimeSymTwoExponential&) = delete;
#else
    TimeSymTwoExponential();
    TimeSymTwoExponential(const TimeSymTwoExponential&);
    TimeSymTwoExponential& operator=(const TimeSymTwoExponential&);
#endif

 public:

    TimeSymTwoExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,1) {
            model_name = "TimeSymTwoExponential";
          }   // nparams = 4, efftype = 1

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymTwoExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

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

#ifndef NO_CXX11
    TimeForwardTwoExponentialPlusConstant() = delete;
    TimeForwardTwoExponentialPlusConstant(const TimeForwardTwoExponentialPlusConstant&) = delete;
    TimeForwardTwoExponentialPlusConstant& operator=(const TimeForwardTwoExponentialPlusConstant&) = delete;
#else
    TimeForwardTwoExponentialPlusConstant();
    TimeForwardTwoExponentialPlusConstant(const TimeForwardTwoExponentialPlusConstant&);
    TimeForwardTwoExponentialPlusConstant& operator=(const TimeForwardTwoExponentialPlusConstant&);
#endif

 public:

    TimeForwardTwoExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(5,in_Tperiod,2) {
            model_name = "TimeForwardTwoExponentialPlusConstant";
          }   // nparams = 5, efftype = 2

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardTwoExponentialPlusConstant(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

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

#ifndef NO_CXX11
    TimeSymTwoExponentialPlusConstant() = delete;
    TimeSymTwoExponentialPlusConstant(const TimeSymTwoExponentialPlusConstant&) = delete;
    TimeSymTwoExponentialPlusConstant& operator=(const TimeSymTwoExponentialPlusConstant&) = delete;
#else
    TimeSymTwoExponentialPlusConstant();
    TimeSymTwoExponentialPlusConstant(const TimeSymTwoExponentialPlusConstant&);
    TimeSymTwoExponentialPlusConstant& operator=(const TimeSymTwoExponentialPlusConstant&);
#endif

 public:

    TimeSymTwoExponentialPlusConstant(uint in_Tperiod) 
          : TemporalCorrelatorModel(5,in_Tperiod,3) {
            model_name = "TimeSymTwoExponentialPlusConstant";
          }   // nparams = 5, efftype = 3

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymTwoExponentialPlusConstant(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double A, double m, double B, double DD, double c0,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, 
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dc0val) const;

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


class TimeForwardTwoExponentialForCons :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TimeForwardTwoExponentialForCons() = delete;
    TimeForwardTwoExponentialForCons(const TimeForwardTwoExponentialForCons&) = delete;
    TimeForwardTwoExponentialForCons& operator=(const TimeForwardTwoExponentialForCons&) = delete;
#else
    TimeForwardTwoExponentialForCons();
    TimeForwardTwoExponentialForCons(const TimeForwardTwoExponentialForCons&);
    TimeForwardTwoExponentialForCons& operator=(const TimeForwardTwoExponentialForCons&);
#endif

 public:

    TimeForwardTwoExponentialForCons(uint in_Tperiod) 
          : TemporalCorrelatorModel(5,in_Tperiod,0) {
        model_name = "TimeForwardTwoExponentialForCons";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergyShift",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio"
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardTwoExponentialForCons(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double m, double A, double shift, double DD, double B, 
                   double t, double& funcval) const;

    void eval_grad(double m, double A, double shift, double DD, double B, 
                   double t, double& dmval, double& dAval, double& dshiftval,
                   double& dDDval, double& dBval) const;
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

#ifndef NO_CXX11
    TimeForwardGeomSeriesExponential() = delete;
    TimeForwardGeomSeriesExponential(const TimeForwardGeomSeriesExponential&) = delete;
    TimeForwardGeomSeriesExponential& operator=(const TimeForwardGeomSeriesExponential&) = delete;
#else
    TimeForwardGeomSeriesExponential();
    TimeForwardGeomSeriesExponential(const TimeForwardGeomSeriesExponential&);
    TimeForwardGeomSeriesExponential& operator=(const TimeForwardGeomSeriesExponential&);
#endif

 public:

    TimeForwardGeomSeriesExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,0) {
            model_name = "TimeForwardGeomSeriesExponential";
          }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardGeomSeriesExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

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

#ifndef NO_CXX11
    TimeSymGeomSeriesExponential() = delete;
    TimeSymGeomSeriesExponential(const TimeSymGeomSeriesExponential&) = delete;
    TimeSymGeomSeriesExponential& operator=(const TimeSymGeomSeriesExponential&) = delete;
#else
    TimeSymGeomSeriesExponential();
    TimeSymGeomSeriesExponential(const TimeSymGeomSeriesExponential&);
    TimeSymGeomSeriesExponential& operator=(const TimeSymGeomSeriesExponential&);
#endif

 public:

    TimeSymGeomSeriesExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,1) {
            model_name = "TimeSymGeomSeriesExponential";
          }   // nparams = 4, efftype = 1

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeSymGeomSeriesExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double A, double m, double B, double DD,
                   double t, int Nt, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD,
                   double t, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval) const;

};


// ******************************************************************************


      // Fitting function is GMO constructed correlator:
      //
      //       f(t) = AL * AS^(1/3) exp( -( mL+mS/3-2*mN/3-2*mX/3 )*t ) /  AN^(2/3) / AX^(2/3)
      //
      // where 
      //           mL = fitparams[0]
      //           AL = fitparams[1]
      //           mS = fitparams[2]
      //           AS = fitparams[3]
      //           mN = fitparams[4]
      //           AN = fitparams[5]
      //           mX = fitparams[6]
      //           AX = fitparams[7].
      //
      // For initial guess, D200 results have been hardcoded in model_corr.cc
      //   Need to change this.


class TimeForwardGMO :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TimeForwardGMO() = delete;
    TimeForwardGMO(const TimeForwardGMO&) = delete;
    TimeForwardGMO& operator=(const TimeForwardGMO&) = delete;
#else
    TimeForwardGMO();
    TimeForwardGMO(const TimeForwardGMO&);
    TimeForwardGMO& operator=(const TimeForwardGMO&);
#endif

    double m_mN_init, m_mL_init, m_mS_init, m_mX_init;
    double m_AN_init, m_AL_init, m_AS_init, m_AX_init;
    
 public:

    TimeForwardGMO(uint in_Tperiod) 
          : TemporalCorrelatorModel(8,in_Tperiod,0) {
            model_name = "TimeForwardGMO";
          }   // nparams = 8, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardGMO(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

    void eval_func(double AL, double mL, double AS, double mS, double AN, double mN, 
                               double AX, double mX, double t, double& funcval) const;

    void eval_grad(double AL, double mL, double AS, double mS, double AN, double mN, 
                               double AX, double mX, double t, double& dALval, double& dmLval, 
                               double& dASval, double& dmSval, double& dANval, double& dmNval, 
                               double& dAXval, double& dmXval) const;

};

class TimeForwardDoubleExpRatio1 :  public TemporalCorrelatorModel 
{

    std::vector<double> m_init_params;

#ifndef NO_CXX11
    TimeForwardDoubleExpRatio1() = delete;
    TimeForwardDoubleExpRatio1(const TimeForwardDoubleExpRatio1&) = delete;
    TimeForwardDoubleExpRatio1& operator=(const TimeForwardDoubleExpRatio1&) = delete;
#else
    TimeForwardDoubleExpRatio1();
    TimeForwardDoubleExpRatio1(const TimeForwardDoubleExpRatio1&);
    TimeForwardDoubleExpRatio1& operator=(const TimeForwardDoubleExpRatio1&);
#endif

 public:

    TimeForwardDoubleExpRatio1(uint in_Tperiod) 
          : TemporalCorrelatorModel(6,in_Tperiod,0) {
        model_name = "TimeForwardDoubleExpRatio1";
        param_names = {
            "Energy",
            "Amplitude",
            "NumGap",
            "NumGapAmp",
            "SHGap",
            "SHGapAmp"
        };
      
    }   // nparams = 8, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardDoubleExpRatio1(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

    void eval_func(double A, double m, double AN, double mN, double ASH1, double mSH1, double tf, double& funcval) const;

    void eval_grad(double A, double m, double AN, double mN, double ASH1, double mSH1, 
                                         double tf, double& dAval, double& dmval,
                                         double& dANval, double& dmNval,double& dASH1val, double& dmSH1val) const;
                                         
     
    void simpleSetFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                          const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                          uint fit_tmax, double chisq_dof, double qual,
                          TCorrFitInfo& fitinfo) const;

    void approachSetFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, uint meff_step, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo, bool added_constant=false) const;

};


class TimeForwardDoubleExpRatio2 :  public TemporalCorrelatorModel 
{

    std::vector<double> m_init_params;

#ifndef NO_CXX11
    TimeForwardDoubleExpRatio2() = delete;
    TimeForwardDoubleExpRatio2(const TimeForwardDoubleExpRatio2&) = delete;
    TimeForwardDoubleExpRatio2& operator=(const TimeForwardDoubleExpRatio2&) = delete;
#else
    TimeForwardDoubleExpRatio2();
    TimeForwardDoubleExpRatio2(const TimeForwardDoubleExpRatio2&);
    TimeForwardDoubleExpRatio2& operator=(const TimeForwardDoubleExpRatio2&);
#endif

 public:

    TimeForwardDoubleExpRatio2(uint in_Tperiod) 
          : TemporalCorrelatorModel(8,in_Tperiod,0) {
        model_name = "TimeForwardDoubleExpRatio2";
        param_names = {
            "Energy",
            "Amplitude",
            "NumGap",
            "NumGapAmp",
            "SH1Gap",
            "SH1GapAmp",
            "SH2Gap",
            "SH2GapAmp"
        };
      
    }   // nparams = 8, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardDoubleExpRatio2(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

    void eval_func(double A, double m, double AN, double mN, double ASH1, double mSH1, 
                                          double ASH2, double mSH2, double tf, double& funcval) const;

    void eval_grad(double A, double m, double AN, double mN, double ASH1, double mSH1, 
                                          double ASH2, double mSH2, double tf, double& dAval, double& dmval,
                                         double& dANval, double& dmNval,double& dASH1val, double& dmSH1val,
                                         double& dASH2val, double& dmSH2val) const;
                                         
     
    void simpleSetFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                          const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                          uint fit_tmax, double chisq_dof, double qual,
                          TCorrFitInfo& fitinfo) const;

    void approachSetFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, uint meff_step, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo, bool added_constant=false) const;

};

// ******************************************************************************


      // Fitting function is independent double exponential time-forward only:
      //
      //       f(t) = A * exp( -m*t ) + A1 * exp( -m1*t )
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1]
      //           m1 = fitparams[0]
      //           A1 = fitparams[1].
      //
      // For initial guess, need corr[tmin], corr[tmin+1]


class TimeForwardTwoIndExp :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TimeForwardTwoIndExp() = delete;
    TimeForwardTwoIndExp(const TimeForwardTwoIndExp&) = delete;
    TimeForwardTwoIndExp& operator=(const TimeForwardTwoIndExp&) = delete;
#else
    TimeForwardTwoIndExp();
    TimeForwardTwoIndExp(const TimeForwardTwoIndExp&);
    TimeForwardTwoIndExp& operator=(const TimeForwardTwoIndExp&);
#endif

 public:

    TimeForwardTwoIndExp(uint in_Tperiod) 
          : TemporalCorrelatorModel(4,in_Tperiod,0) {
            model_name = "TimeForwardTwoIndExp";
          }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardTwoIndExp(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    static void setup(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount);

    void eval_func(double A, double m, double A1, double m1, double t, double& funcval) const;

    void eval_grad(double A, double m, double A1, double m1, double t, 
                double& dAval, double& dmval, double& dA1val, double& dm1val) const;


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


class TimeForwardThreeExponential :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TimeForwardThreeExponential() = delete;
    TimeForwardThreeExponential(const TimeForwardThreeExponential&) = delete;
    TimeForwardThreeExponential& operator=(const TimeForwardThreeExponential&) = delete;
#else
    TimeForwardThreeExponential();
    TimeForwardThreeExponential(const TimeForwardThreeExponential&);
    TimeForwardThreeExponential& operator=(const TimeForwardThreeExponential&);
#endif

 public:

    TimeForwardThreeExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(6,in_Tperiod,0) {
        model_name = "TimeForwardThreeExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio"
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardThreeExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double A, double m, double B, double DD, double C, double DDD,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double C, double DDD,
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval,
                   double& dCval, double& dDDDval) const;
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


class TwoExpConspiracy :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TwoExpConspiracy() = delete;
    TwoExpConspiracy(const TwoExpConspiracy&) = delete;
    TwoExpConspiracy& operator=(const TwoExpConspiracy&) = delete;
#else
    TwoExpConspiracy();
    TwoExpConspiracy(const TwoExpConspiracy&);
    TwoExpConspiracy& operator=(const TwoExpConspiracy&);
#endif

 public:

    TwoExpConspiracy(uint in_Tperiod) 
          : TemporalCorrelatorModel(10,in_Tperiod,0) {
        model_name = "TwoExpConspiracy";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SqrtGapToThirdEnergy",
            "SecondAmplitudeRatio",
            "ThirdAmplitudeRatio",
            "FourthAmplitudeRatio",
            "delta2",
            "delta3",
            "delta4",
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TwoExpConspiracy(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double m, double A, double DD1, double DD2, double B, double C, double D, 
                    double d2, double d3, double d4, double t, double& funcval) const;

    void eval_grad(double m, double A, double DD1, double DD2, double B, double C, double D, 
                    double d2, double d3, double d4, double t, 
                    double& dmval, double& dAval, double& dDD1val, double& dDD2val, double& dBval, 
                    double& dCval, double& dDval, double& dd2val, double& dd3val, double& dd4val) const;
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


class DegTwoExpConspiracy :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    DegTwoExpConspiracy() = delete;
    DegTwoExpConspiracy(const DegTwoExpConspiracy&) = delete;
    DegTwoExpConspiracy& operator=(const DegTwoExpConspiracy&) = delete;
#else
    DegTwoExpConspiracy();
    DegTwoExpConspiracy(const DegTwoExpConspiracy&);
    DegTwoExpConspiracy& operator=(const DegTwoExpConspiracy&);
#endif

 public:

    DegTwoExpConspiracy(uint in_Tperiod) 
          : TemporalCorrelatorModel(7,in_Tperiod,0) {
        model_name = "DegTwoExpConspiracy";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "delta1",
            "ThirdAmplitudeRatio",
            "delta2",
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~DegTwoExpConspiracy(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double A, double m, double B, double DD, double C, double d1, double d2,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double C, double d1, double d2, 
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval,
                   double& dCval, double& dd1, double& dd2) const;
};

// ******************************************************************************


class DegThreeExpConspiracy :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    DegThreeExpConspiracy() = delete;
    DegThreeExpConspiracy(const DegThreeExpConspiracy&) = delete;
    DegThreeExpConspiracy& operator=(const DegThreeExpConspiracy&) = delete;
#else
    DegThreeExpConspiracy();
    DegThreeExpConspiracy(const DegThreeExpConspiracy&);
    DegThreeExpConspiracy& operator=(const DegThreeExpConspiracy&);
#endif

 public:

    DegThreeExpConspiracy(uint in_Tperiod) 
          : TemporalCorrelatorModel(14,in_Tperiod,0) {
        model_name = "DegThreeExpConspiracy";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "delta2",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio",
            "delta3",
            "FourthAmplitudeRatio",
            "delta4",
            "FifthAmplitudeRatio",
            "delta5",
            "SixthAmplitudeRatio",
            "delta6",
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~DegThreeExpConspiracy(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double m, double A, 
                   double DD, double B, double d2, 
                   double DDD, double C, double d3, 
                   double D, double d4,
                   double E, double d5,
                   double F, double d6,
                   double t, double& funcval) const;

    void eval_grad( double m, double A,
                   double DD, double B, double d2, 
                   double DDD, double C, double d3, 
                   double D, double d4,
                   double E, double d5,
                   double F, double d6, 
                   double t,
                   double dmval, double dAval,
                   double dDDval, double dBval, double dd2val, 
                   double dDDDval, double dCval, double dd3val, 
                   double dDval, double dd4val,
                   double dEval, double dd5val,
                   double dFval , double dd6val 
                   ) const;
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


class TimeForwardThreeIndExponential :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TimeForwardThreeIndExponential() = delete;
    TimeForwardThreeIndExponential(const TimeForwardThreeIndExponential&) = delete;
    TimeForwardThreeIndExponential& operator=(const TimeForwardThreeIndExponential&) = delete;
#else
    TimeForwardThreeIndExponential();
    TimeForwardThreeIndExponential(const TimeForwardThreeIndExponential&);
    TimeForwardThreeIndExponential& operator=(const TimeForwardThreeIndExponential&);
#endif

 public:

    TimeForwardThreeIndExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(6,in_Tperiod,0) {
        model_name = "TimeForwardThreeIndExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SingleHadronEnergy",
            "ThirdAmplitudeRatio"
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardThreeIndExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double A, double m, double B, double DD, double C, double mSH,
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double C, double mSH, 
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval,
                   double& dCval, double& dmSHval) const;
};

// ******************************************************************************

class TimeForwardFourExponential :  public TemporalCorrelatorModel 
{

#ifndef NO_CXX11
    TimeForwardFourExponential() = delete;
    TimeForwardFourExponential(const TimeForwardFourExponential&) = delete;
    TimeForwardFourExponential& operator=(const TimeForwardFourExponential&) = delete;
#else
    TimeForwardFourExponential();
    TimeForwardFourExponential(const TimeForwardFourExponential&);
    TimeForwardFourExponential& operator=(const TimeForwardFourExponential&);
#endif

 public:

    TimeForwardFourExponential(uint in_Tperiod) 
          : TemporalCorrelatorModel(7,in_Tperiod,0) {
        model_name = "TimeForwardFourExponential";
        param_names = {
            "FirstEnergy",
            "FirstAmplitude",
            "SqrtGapToSecondEnergy",
            "SecondAmplitudeRatio",
            "SqrtGapToThirdEnergy",
            "ThirdAmplitudeRatio",
            "FourthAmplitudeRatio",
        };
      
    }   // nparams = 4, efftype = 0

    virtual void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount);

    virtual void evaluate(const std::vector<double>& fitparams, double tval, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tval, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<uint>& tvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~TimeForwardFourExponential(){}

    virtual void setFitInfo(const std::vector<MCObsInfo>& fitparams_info,
                            const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                            uint fit_tmax, bool show_approach,
                            uint meff_timestep, double chisq_dof, double qual,
                            TCorrFitInfo& fitinfo) const;

 private:

    void eval_func(double A, double m, double B, double DD, double C, double DDD, double E, 
                   double t, double& funcval) const;

    void eval_grad(double A, double m, double B, double DD, double C, double DDD, double E, 
                   double t, double& dAval, double& dmval,
                   double& dBval, double& dDDval,
                   double& dCval, double& dDDDval, double& dEval) const;
};


// ******************************************************************************
#endif