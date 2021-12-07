#ifndef MODEL_3PTCORR_H
#define MODEL_3PTCORR_H
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
// *   The important base class "ThreePointCorrelatorModel" is defined in this      *
// *   file.  A variety of classes derived from "ThreePointCorrelatorModel" are     *
// *   also defined.  Each class is used for a different model fit function.      *
// *   The constructor of the base class has the form                             *
// *                                                                              *
// *    ThreePointCorrelatorModel(in_nparams,in_Tperiod,in_efftype);                *
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
// *   Classes derived from "ThreePointCorrelatorModel" below all have a            *
// *   constructor of the form                                                    *
// *                                                                              *
// *         DerivedClass(in_Tperiod);                                            *
// *                                                                              *
// *    The classes derived below from "ThreePointCorrelatorModel" are              *
// *                                                                              *
// *        SingleExponential                                          *
// *                                                                              *
// *                                                                              *
// *   A useful routine for dynamically allocating an object of base class        *
// *  "ThreePointCorrelatorModel" given a model type specified by a string is       *
// *   also defined in this file:                                                 *
// *                                                                              *
// *      void create_tcorr_model(const std::string& modeltype,                   *
// *                              uint in_Tperiod,                                *
// *                              ThreePointCorrelatorModel* &mptr);                *
// *                                                                              *
// *   Polymorphism is used here.  An object of a derived class is actually       *
// *   created, but it can be accessed through the pointer to the base class.     *
// *                                                                              *
// ********************************************************************************

enum ObsType { Ratio, SummationRatio };

class ThreePointCorrelatorModel
{

 protected:     // derived classes have access to the protected members

    std::string model_name;
    std::vector<std::string> param_names;
    ObsType obs_type;

    uint T_period;   // temporal extent of lattice in number of sites


 private:
          // disallow copying

    ThreePointCorrelatorModel() = delete;
    ThreePointCorrelatorModel(const ThreePointCorrelatorModel&) = delete;
    ThreePointCorrelatorModel& operator=(const ThreePointCorrelatorModel&) = delete;

 protected:

    ThreePointCorrelatorModel(uint in_Tperiod) 
                : T_period(in_Tperiod) {}

 public:

    void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) const;

    void setupInfos(std::map<std::string,MCObsInfo> model_params, std::vector<MCObsInfo>& fitparam_info) const;

    virtual void evaluate(const std::vector<double>& fitparams, double tsep, double tins, double& value) const = 0;

    virtual void evalGradient(const std::vector<double>& fitparams, double tsep, double tins,
                              std::vector<double>& grad) const = 0;

    virtual void guessInitialParamValues(const std::vector<double>& data,
                                         const std::map<uint,std::set<uint> >& tvals,
                                         std::vector<double>& fitparam) const = 0;    

    void output_tag(XMLHandler& xmlout) const;

    virtual ~ThreePointCorrelatorModel(){}

    ObsType getObsType() const
     {return obs_type;}

    std::string getModelName() const
     {return model_name;}

    const std::vector<std::string>& getParameterNames() const
    {return param_names;}

    std::string getParameterName(uint param_index) const
     {return param_names[param_index];}

    uint getNumberOfParams() const
     {return param_names.size();}

};

// *****************************************************************************

     //   A useful routine for dynamically allocating an object of base
     //   class "ThreePointCorrelatorModel" given a model type specified by
     //   a string.  Polymorphism is used here.  An object of a derived
     //   class is actually created, but it can be accessed through the
     //   pointer to the base class.

void create_3ptcorr_model(const std::string& modeltype, uint in_Tperiod,
                          ThreePointCorrelatorModel* &mptr);




// *****************************************************************************
// *
// *
// *        Currently supported models:    names+ID index should be unique in each run
// *
// *                (B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep))
// *                   / (1 + A*exp(-dE*tsep))
// *         <Model>
// *             <Type>TwoStateRatio</Type>          
// *             <B0>
// *                <Name>b0</Name><IDIndex>0</IDIndex>
// *             </B0>
// *             <B1>
// *                <Name>b1</Name><IDIndex>0</IDIndex>
// *             </B1>
// *             <B2>
// *                <Name>b2</Name><IDIndex>0</IDIndex>
// *             </B2>
// *             <B3>
// *                <Name>b3</Name><IDIndex>0</IDIndex>
// *             </B3>
// *             <A>
// *                <Name>a</Name><IDIndex>0</IDIndex>
// *             </A>
// *             <EnergyShift>
// *                <Name>e_shift</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *         </Model>
// *
// *
// *                B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep)
// *         <Model>
// *             <Type>TwoStateNumerator</Type>          
// *             <B0>
// *                <Name>b0</Name><IDIndex>0</IDIndex>
// *             </B0>
// *             <B1>
// *                <Name>b1</Name><IDIndex>0</IDIndex>
// *             </B1>
// *             <B2>
// *                <Name>b2</Name><IDIndex>0</IDIndex>
// *             </B2>
// *             <B3>
// *                <Name>b3</Name><IDIndex>0</IDIndex>
// *             </B3>
// *             <EnergyShift>
// *                <Name>e_shift</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *         </Model>
// *
// *
// *                B0 + B1*[exp(-dE(tsep - tins)) + exp(-dE*tins)] + B2*exp(-dE*tsep)
// *         <Model>
// *             <Type>TwoStateNumeratorSymmetric</Type>          
// *             <B0>
// *                <Name>b0</Name><IDIndex>0</IDIndex>
// *             </B0>
// *             <B1>
// *                <Name>b1</Name><IDIndex>0</IDIndex>
// *             </B1>
// *             <B2>
// *                <Name>b2</Name><IDIndex>0</IDIndex>
// *             </B2>
// *             <EnergyShift>
// *                <Name>e_shift</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *         </Model>
// *
// *
// *                B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins)
// *         <Model>
// *             <Type>TwoStateSimple</Type>          
// *             <B0>
// *                <Name>b0</Name><IDIndex>0</IDIndex>
// *             </B0>
// *             <B1>
// *                <Name>b1</Name><IDIndex>0</IDIndex>
// *             </B1>
// *             <B2>
// *                <Name>b2</Name><IDIndex>0</IDIndex>
// *             </B2>
// *             <EnergyShift>
// *                <Name>e_shift</Name><IDIndex>0</IDIndex>
// *             </Energy>
// *         </Model>
// *
// *
// *                A + B0*tsep
// *         <Model>
// *             <Type>SummationLinear</Type>          
// *             <B0>
// *                <Name>b0</Name><IDIndex>0</IDIndex>
// *             </B0>
// *             <A>
// *                <Name>a</Name><IDIndex>0</IDIndex>
// *             </A>
// *         </Model>
// *
// *
// ******************************************************************************


      // Fitting function is single exponential time-forward only:
      //
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
      //


class TwoStateRatio :  public ThreePointCorrelatorModel 
{

    TwoStateRatio() = delete;
    TwoStateRatio(const TwoStateRatio&) = delete;
    TwoStateRatio& operator=(const TwoStateRatio&) = delete;

 public:

    TwoStateRatio(uint in_Tperiod) 
          : ThreePointCorrelatorModel(in_Tperiod)
    {
        model_name = "TwoStateRatio";
        param_names = {
            "B0",
            "B1",
            "B2",
            "B3",
            "A",
            "EnergyShift"
        };
        obs_type = Ratio;
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tsep, double tins, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tsep, double tins,
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data,
                                         const std::map<uint,std::set<uint> >& tvals,
                                         std::vector<double>& fitparam) const;    

    virtual ~TwoStateRatio(){}

 private:

    void eval_func(double B0, double B1, double B2, double B3, double A, double dE,
                   double tsep, double tins, double& funcval) const;

    void eval_grad(double B0, double B1, double B2, double B3, double A, double dE,
                   double tsep, double tins, double& dB0val, double& dB1val, double& dB2val,
                   double& dB3val, double& dAval, double& ddEval) const;
};


      // Fitting function is single exponential time-forward only:
      //
      //       f(tsep, tins) = B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins) + B3*exp(-dE*tsep)
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           B3 = fitparams[3].
      //           dE = fitparams[4]
      //


class TwoStateNumerator :  public ThreePointCorrelatorModel 
{

    TwoStateNumerator() = delete;
    TwoStateNumerator(const TwoStateNumerator&) = delete;
    TwoStateNumerator& operator=(const TwoStateNumerator&) = delete;

 public:

    TwoStateNumerator(uint in_Tperiod) 
          : ThreePointCorrelatorModel(in_Tperiod)
    {
        model_name = "TwoStateNumerator";
        param_names = {
            "B0",
            "B1",
            "B2",
            "B3",
            "EnergyShift"
        };
        obs_type = Ratio;
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tsep, double tins, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tsep, double tins,
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data,
                                         const std::map<uint,std::set<uint> >& tvals,
                                         std::vector<double>& fitparam) const;    

    virtual ~TwoStateNumerator(){}

 private:

    void eval_func(double B0, double B1, double B2, double B3, double dE,
                   double tsep, double tins, double& funcval) const;

    void eval_grad(double B0, double B1, double B2, double B3, double dE,
                   double tsep, double tins, double& dB0val, double& dB1val,
                   double& dB2val, double& dB3val, double& ddEval) const;
};


      // Fitting function is single exponential time-forward only:
      //
      //       f(tsep, tins) = B0 + B1*[exp(-dE(tsep - tins)) + exp(-dE*tins)] + B3*exp(-dE*tsep)
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           dE = fitparams[4]
      //


class TwoStateNumeratorSymmetric :  public ThreePointCorrelatorModel 
{

    TwoStateNumeratorSymmetric() = delete;
    TwoStateNumeratorSymmetric(const TwoStateNumeratorSymmetric&) = delete;
    TwoStateNumeratorSymmetric& operator=(const TwoStateNumeratorSymmetric&) = delete;

 public:

    TwoStateNumeratorSymmetric(uint in_Tperiod) 
          : ThreePointCorrelatorModel(in_Tperiod)
    {
        model_name = "TwoStateNumeratorSymmetric";
        param_names = {
            "B0",
            "B1",
            "B2",
            "EnergyShift"
        };
        obs_type = Ratio;
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tsep, double tins, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tsep, double tins,
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data,
                                         const std::map<uint,std::set<uint> >& tvals,
                                         std::vector<double>& fitparam) const;    

    virtual ~TwoStateNumeratorSymmetric(){}

 private:

    void eval_func(double B0, double B1, double B2, double dE,
                   double tsep, double tins, double& funcval) const;

    void eval_grad(double B0, double B1, double B2, double dE,
                   double tsep, double tins, double& dB0val, double& dB1val,
                   double& dB2val, double& ddEval) const;
};

      // Fitting function is single exponential time-forward only:
      //
      //       f(tsep, tins) = B0 + B1*exp(-dE(tsep - tins)) + B2*exp(-dE*tins)
      //
      // where 
      //           B0 = fitparams[0].
      //           B1 = fitparams[1].
      //           B2 = fitparams[2].
      //           dE = fitparams[3]
      //


class TwoStateSimple :  public ThreePointCorrelatorModel 
{

    TwoStateSimple() = delete;
    TwoStateSimple(const TwoStateSimple&) = delete;
    TwoStateSimple& operator=(const TwoStateSimple&) = delete;

 public:

    TwoStateSimple(uint in_Tperiod) 
          : ThreePointCorrelatorModel(in_Tperiod)
    {
        model_name = "TwoStateSimple";
        param_names = {
            "B0",
            "B1",
            "B2",
            "EnergyShift"
        };
        obs_type = Ratio;
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tsep, double tins, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tsep, double tins,
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data,
                                         const std::map<uint,std::set<uint> >& tvals,
                                         std::vector<double>& fitparam) const;    

    virtual ~TwoStateSimple(){}

 private:

    void eval_func(double B0, double B1, double B2, double dE,
                   double tsep, double tins, double& funcval) const;

    void eval_grad(double B0, double B1, double B2, double dE,
                   double tsep, double tins, double& dB0val,
                   double& dB1val, double& dB2val, double& ddEval) const;
};

      // Fitting function is single exponential time-forward only:
      //
      //       f(tsep) = A + B0*tsep
      //
      // where 
      //           B0 = fitparams[0].
      //            A = fitparams[1].
      //
      // For initial guess, need corr[tmin], corr[tmin+1]


class SummationLinear :  public ThreePointCorrelatorModel 
{

    SummationLinear() = delete;
    SummationLinear(const SummationLinear&) = delete;
    SummationLinear& operator=(const SummationLinear&) = delete;

 public:

    SummationLinear(uint in_Tperiod) 
          : ThreePointCorrelatorModel(in_Tperiod)
    {
        model_name = "SummationLinear";
        param_names = {
            "B0",
            "A"
        };
        obs_type = SummationRatio;
    }

    virtual void evaluate(const std::vector<double>& fitparams, double tsep, double tins, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double tsep, double tins,
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data,
                                         const std::map<uint,std::set<uint> >& tvals,
                                         std::vector<double>& fitparam) const;    

    virtual ~SummationLinear(){}

 private:

    void eval_func(double B0, double A, double tsep, double& funcval) const;

    void eval_grad(double B0, double A, double tsep, double& dB0val, double& dAval) const;
};

// ******************************************************************************
#endif
