#ifndef MODEL_DISP_H
#define MODEL_DISP_H
#include <vector>
#include <algorithm>
#include "mcobs_info.h"
#include "grace_plot.h"
#include "mc_estimate.h"

// ********************************************************************************
// *                                                                              *
// *   The important base class "DispersionModel" is defined in this              *
// *   file.  A variety of classes derived from "DispersionModel" are             *
// *   also defined.  Each class is used for a different model fit function.      *
// *   The constructor of the base class has the form                             *
// *                                                                              *
// *    DispersionModel(in_nparams, in_Xextent, in_Yextent, in_Zextent);          *
// *                                                                              *
// *   This sets the number of parameters and the spatial extents of the lattice. *
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
// *   Classes derived from "DispersionModel" below all have a                    *
// *   constructor of the form                                                    *
// *                                                                              *
// *         DerivedClass(in_Xextent, in_Yextent, in_Zextent);                    *
// *                                                                              *
// *    The classes derived below from "DispersionModel" are                      *
// *                                                                              *
// *        Anisotropy                                                            *
// *        Continuum                                                             *
// *        FourthOrderMomentum                                                   *
// *                                                                              *
// *                                                                              *
// *   A useful routine for dynamically allocating an object of base class        *
// *  "DispersionModel" given a model type specified by a string is               *
// *   also defined in this file:                                                 *
// *                                                                              *
// *      void create_disp_model(const std::string& modeltype,                    *
// *                             uint in_Xextent, uint in_Yextent, uint in_Zextent, *
// *                              DispersionModel* &mptr);                        *
// *                                                                              *
// *   Polymorphism is used here.  An object of a derived class is actually       *
// *   created, but it can be accessed through the pointer to the base class.     *
// *                                                                              *
// ********************************************************************************



class DispersionModel
{

 protected:     // derived classes have access to the protected members

    std::string model_name;
    std::vector<std::string> param_names;

    uint m_nparams;  // number of fit parameters

    uint Xextent; 
    uint Yextent; 
    uint Zextent; 


 private:
          // disallow copying

    DispersionModel() = delete;
    DispersionModel(const DispersionModel&) = delete;
    DispersionModel& operator=(const DispersionModel&) = delete;

 protected:

    DispersionModel(uint in_nparams, uint in_Xextent, uint in_Yextent, uint in_Zextent) 
                : m_nparams(in_nparams), Xextent(in_Xextent), Yextent(in_Yextent), Zextent(in_Zextent) {}

 public:

    void setupInfos(XMLHandler& xmlin, std::vector<MCObsInfo>& fitparam_info, int taskcount) const;

    void setupInfos(std::map<std::string,MCObsInfo> model_params, std::vector<MCObsInfo>& fitparam_info) const;

    virtual void evaluate(const std::vector<double>& fitparams, double psq, 
                          double& value) const = 0;

    virtual void evalGradient(const std::vector<double>& fitparams, double psq, 
                              std::vector<double>& grad) const = 0;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<double>& psqvals,
                                         std::vector<double>& fitparam) const = 0;    

    void output_tag(XMLHandler& xmlout) const;

    virtual ~DispersionModel(){}

    std::string getModelName() const
     {return model_name;}

    const std::vector<std::string>& getParameterNames() const
    {return param_names;}

    std::string getParameterName(uint param_index) const
     {return param_names[param_index];}

    uint getNumberOfParams() const
     {return m_nparams;}

};

// *****************************************************************************

     //   A useful routine for dynamically allocating an object of base
     //   class "DispersionModel" given a model type specified by
     //   a string.  Polymorphism is used here.  An object of a derived
     //   class is actually created, but it can be accessed through the
     //   pointer to the base class.

void create_disp_model(const std::string& modeltype, uint in_Xextent, uint in_Yextent, uint in_Zextent,
                       DispersionModel* &mptr);




// *****************************************************************************
// *
// *          let p^2 = (2*Pi)^2 * sum_i (n_i / Ni)^2  (n_i is integer momentum in ith dir
// *                                                    Ni is number of sites in ith dir)
// *
// *        Currently supported models:    names+ID index should be unique in each run
// *
// *                  (a_t E)^2 = restmass_sq + p^2 / xi^2
// *
// *         <Model>
// *             <Type>Anisotropy</Type>          
// *             <RestMassSq>
// *                <Name>pion_mass</Name><IDIndex>0</IDIndex>
// *             </RestMassSq>
// *             <Anisotropy>
// *                <Name>aniso</Name><IDIndex>0</IDIndex>
// *             </Anistropy>
// *         </Model>
// *
// *
// *                  (a_t E)^2 = restmass_sq + c * p^2
// *
// *         <Model>
// *             <Type>Continuum</Type>          
// *             <RestMassSq>
// *                <Name>pion_mass</Name><IDIndex>0</IDIndex>
// *             </RestMassSq>
// *             <MomentumCoefficient>
// *                <Name>mom_coeff</Name><IDIndex>0</IDIndex>
// *             </MomentumCoefficient>
// *         </Model>
// *
// *                  (a_t E)^2 = restmass_sq + c1 * p^2 - c2 * [(restmass_sq)^2 + 2 * p^4 + restmass_sq * p^2]
// *         <Model>
// *             <Type>FourthOrderMomentum</Type>          
// *             <RestMassSq>
// *                <Name>pion_mass</Name><IDIndex>0</IDIndex>
// *             </RestMassSq>
// *             <FirstMomentumCoefficient>
// *                <Name>mom_coeff_1</Name><IDIndex>0</IDIndex>
// *             </FirstMomentumCoefficient>
// *             <SecondMomentumCoefficient>
// *                <Name>mom_coeff_2</Name><IDIndex>0</IDIndex>
// *             </SecondMomentumCoefficient>
// *         </Model>
// *
// *

// ******************************************************************************


      // Fitting function is 
      //
      //          (a_t E)^2 = restmass_sq + p^2 / xi^2
      //
      // where 
      //           restmass_sq = fitparams[0]
      //                    xi = fitparams[1].


class Anisotropy :  public DispersionModel 
{

    Anisotropy() = delete;
    Anisotropy(const Anisotropy&) = delete;
    Anisotropy& operator=(const Anisotropy&) = delete;

 public:

    Anisotropy(uint in_Xextent, uint in_Yextent, uint in_Zextent) 
          : DispersionModel(2, in_Xextent, in_Yextent, in_Zextent)    // nparams = 2
    {
        model_name = "Anisotropy";
        param_names = {
            "RestMassSq",
            "Anisotropy"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double psq, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double psq, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<double>& psqvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~Anisotropy() {}

 private:

    void eval_func(double msq, double xi, double psq, double& funcval) const;

    void eval_grad(double xi, double psq, double& dmsqval, double& dxival) const;

};


// ******************************************************************************


      // Fitting function is 
      //
      //          (a_t E)^2 = restmass_sq + c * p^2
      //
      // where 
      //           restmass_sq = fitparams[0]
      //                     c = fitparams[1].


class Continuum :  public DispersionModel 
{

    Continuum() = delete;
    Continuum(const Continuum&) = delete;
    Continuum& operator=(const Continuum&) = delete;

 public:

    Continuum(uint in_Xextent, uint in_Yextent, uint in_Zextent) 
          : DispersionModel(2, in_Xextent, in_Yextent, in_Zextent)    // nparams = 2
    {
        model_name = "Continuum";
        param_names = {
            "RestMassSq",
            "MomentumCoefficient"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double psq, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double psq, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<double>& psqvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~Continuum(){}

 private:

    void eval_func(double msq, double c, double psq, double& funcval) const;

    void eval_grad(double psq, double& dmsqval, double& dcval) const;

};


// ******************************************************************************


      // Fitting function is 
      //
      //      (a_t E)^2 = restmass_sq + c1 * p^2 - c2 * [(restmass_sq)^2 + 2 * p^4 + restmass_sq * p^2]
      //
      // where 
      //           restmass_sq = fitparams[0]
      //                    c1 = fitparams[1].
      //                    c2 = fitparams[2].


class FourthOrderMomentum :  public DispersionModel 
{

    FourthOrderMomentum() = delete;
    FourthOrderMomentum(const FourthOrderMomentum&) = delete;
    FourthOrderMomentum& operator=(const FourthOrderMomentum&) = delete;

 public:

    FourthOrderMomentum(uint in_Xextent, uint in_Yextent, uint in_Zextent) 
          : DispersionModel(3, in_Xextent, in_Yextent, in_Zextent)    // nparams = 3
    {
        model_name = "FourthOrderMomentum";
        param_names = {
            "RestMassSq",
            "FirstMomentumCoefficient",
            "SecondMomentumCoefficient"
        };
    }

    virtual void evaluate(const std::vector<double>& fitparams, double psq, double& value) const;

    virtual void evalGradient(const std::vector<double>& fitparams, double psq, 
                              std::vector<double>& grad) const;

    virtual void guessInitialParamValues(const std::vector<double>& data, const std::vector<double>& psqvals, 
                                         std::vector<double>& fitparam) const;    

    virtual ~FourthOrderMomentum() {}

 private:

    void eval_func(double msq, double c1, double c2, double psq, double& funcval) const;

    void eval_grad(double msq, double c2, double psq, double& dmsqval, double& dc1val, double& dc2val) const;

};


// ******************************************************************************
#endif
