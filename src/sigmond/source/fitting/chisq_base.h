#ifndef CHISQ_BASE_H
#define CHISQ_BASE_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"


// ********************************************************************************
// *                                                                              *
// *   "ChiSquare" is the all-important base class for correlated least-squares   *
// *   fitting in SigMonD.  It stores information about the observables and the   *
// *   fit parameters.  Using a pointer to an MCObsHandler, it can compute the    *
// *   estimates of the observables, as well as their covariances.  The base      *
// *   class also carries out the Cholesky decomposition of the inverse of the    *
// *   covariance matrix.  An object of the base class stores the mean values     *
// *   and the Cholesky lower triangular matrix.  The base class cannot be        *
// *   constructed on its own since it contains several purely virtual member     *
// *   functions.                                                                 *
// *                                                                              *
// *   Observable mean values are computed using the current resampling mode in   *
// *   the associated MCObsHandler, while the covariances are computed using the  *
// *   current **covariance matrix** resampling mode in the MCObsHandler.         *
// *   Hence, a different mode can be used when computing the covariance versus   *
// *   the observable means.                                                      *
// *                                                                              *
// *   For the fitting, the "cost" function is a correlated chi-square:           *
// *                                                                              *
// *    chi-square = sum_j (model[j]-obs[j]) * inv_cov(j,k) * (model[k]-obs[k])   *
// *                                                                              *
// *   where inv_cov = the inverse of the covariance matrix given by              *
// *                                                                              *
// *     cov(j,k) = cov( obs[j], obs[k] )                                         *
// *                                                                              *
// *   A very important assumption here is that the model is INDEPENDENT of the   *
// *   observations.   The model can only depend on various parameters to         *
// *   predict the observations, and the functions which evaluate these           *
// *   predictions cannot use any of the observations.  (This could be changed    *
// *   by using cov(k) = cov( obs[j]-model[j], obs[k]-model[j] ), but this        *
// *   requires recomputation every time the fit parameters are adjusted during   *
// *   the fit, which is undesirable.)                                            *
// *                                                                              *
// *   A Cholesky decomposition of inv_cov can be done:                           *
// *              inv_cov = transpose(L) * L,   L = lower triangular              *
// *   The fit parameters determine the model[j] values, and the fit parameters   *
// *   are adjusted until the chi-square is a minimum.                            *
// *                                                                              *
// *                                                                              *
// *   A class "DerivedFit" derived from "ChiSquare" must have a constructor      *
// *   of the form                                                                *
// *                                                                              *
// *         DerivedFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);      *
// *                                                                              *
// *   This constructor must                                                      *
// *     --  call the base constructor with the initializer :  ChiSquare(OH)      *
// *     --  first, set "m_nobs" the number of observables                        *
// *     --  second, set "m_nparams" the number of parameters                     *
// *     --  then call "allocate_obs_memory()" to set up the needed memory        *
// *     --  use information in "xmlin" to initialize the data members            *
// *              std::vector<MCObsInfo> m_obs_info;                              *
// *              std::vector<MCObsInfo> m_fitparam_info;                         *
// *                                                                              *
// *   The derived class must define the members below.  NOTE: memory is          *
// *   already allocated, so the members below should not resize the              *
// *   "fitparams", etc. vectors and matrices.                                    *
// *                                                                              *
// *      virtual void evalModelPoints(const vector<double>& fitparams,           *
// *                                   vector<double>& modelpoints) const;        *
// *                                                                              *
// *           // gradients(i,p) is the derivative of the i-th                    *
// *           // observable with respect to the p-th parameter.                  *
// *                                                                              *
// *      virtual void evalGradients(const vector<double>& fitparams,             *
// *                                 RMatrix& gradients) const;                   *
// *                                                                              *
// *      virtual void guessInitialParamValues(const RVector& datapoints,         *
// *                                           vector<double>& fitparams) const;  *
// *                                                                              *
// *      virtual void do_output(XMLHandler& xmlout) const;                       *
// *                                                                              *
// *                                                                              *
// *                                                                              *
// *   Objects of classes derived from "ChiSquare" will generally be accessed     *
// *   through either a pointer or a reference to the base class to allow         *
// *   the needed polymorphism for the different fittings.  The important         *
// *   routines are described below:                                              *
// *                                                                              *
// *     XMLHandler xmlin;                                                        *
// *     MCObsHandler OH;                                                         *
// *     ChiSquare& CHSQ=new DerivedFit(xmlin,OH,taskcount);                      *
// *                                                                              *
// *     CHSQ.getNumberOfParams();       // returns number of parameters          *
// *                                                                              *
// *     CHSQ.getNumberOfObervables();   // number of observables                 *
// *                                                                              *
// *     for (OH.begin();!OH.end();++OH){    // iterate over resamplings          *
// *                                                                              *
// *        CHSQ.setObsMean();     // evaluate and store observable               *
// *                               // means using current resampling mode         *
// *                                                                              *
// *        CHSQ.setObsMeanCov();  // evaluate and store observable               *
// *                               // means using current resampling mode and     *
// *                               // also compute and store covariance using     *
// *                               // current covariance matrix resampling mode   *
// *                                                                              *
// *        vector<double> fitparams(nparams);                                    *
// *        CHSQ.guessInitialFitParamValues(fitparams); // set initial            *
// *                                                    // parameter values       *
// *                                                                              *
// *        vector<double> residual(nobs);                                        *
// *        CHSQ.evalResiduals(fitparams,residuals);                              *
// *             // updates the residual values using "fitparams", the fit        *
// *             // parameter values, that were recently changed;                 *
// *             // obtain i-th residual from res[i].  Note that these            *
// *             // residuals do not correspond to the observables, since         *
// *             // a Cholesky "rotation" has been done to convert the            *
// *             // correlated-chisquare into an ordinary sum of squares          *
// *                                                                              *
// *        RMatrix gradients(nobs,nparams);                                      *
// *        CHSQ.evalResGradients(fitparams,gradients);                           *
// *             // updates the residual gradients using the "fitparams"          *
// *             // fit parameter values,  grad(i,p) is the derivative of the     *
// *             // i-th observable with respect to the p-th parameter.           *
// *                                                                              *
// *        double chisq=CHSQ.evalChiSquare(residuals);  }                        *
// *             // get the value of the chi-square sum of                        *
// *             // residual squares                                              *
// *                                                                              *
// *                                                                              *
// ********************************************************************************

 

class ChiSquare
{

 protected:     // derived classes have access to the protected members

    MCObsHandler *m_obs;
    uint m_nobs;
    uint m_nparams;
    std::vector<MCObsInfo> m_obs_info; 
    std::vector<MCObsInfo> m_fitparam_info;
    RVector m_means;
    LowerTriangularMatrix<double> m_inv_cov_cholesky;


    ChiSquare(MCObsHandler& OH) : m_obs(&OH) {}

    virtual ~ChiSquare(){}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const = 0;

           // gradients(i,p) is the derivative of the i-th
           // observable with respect to the p-th parameter.

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const = 0;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const = 0;

    virtual void do_output(XMLHandler& xmlout) const = 0;

    void allocate_obs_memory();


 private:

#ifndef NO_CXX11
    ChiSquare() = delete;
    ChiSquare(const ChiSquare&) = delete;
    ChiSquare& operator=(const ChiSquare&) = delete;
#else
    ChiSquare();
    ChiSquare(const ChiSquare&);
    ChiSquare& operator=(const ChiSquare&);
#endif


 public:

    uint getNumberOfParams() const
     {return m_nparams;}

    uint getNumberOfObervables() const
     {return m_nobs;}

    const std::vector<MCObsInfo>& getObsInfos() const 
     {return m_obs_info;}

    const std::vector<MCObsInfo>& getFitParamInfos() const 
     {return m_fitparam_info;}

    MCObsHandler* getMCObsHandlerPtr() { return m_obs;}

    SamplingMode getObsMeansSamplingMode() const 
     {return m_obs->getCurrentSamplingMode();}

    SamplingMode getCovMatSamplingMode() const 
     {return m_obs->getCovMatCurrentSamplingMode();}

    void output(XMLHandler& xmlout) const;

    std::string output(int indent=0) const;

    std::string str() const;

          // calculates and stores observable means with current sampling mode

    void setObsMean(); 

          // calculates and stores observable means with current sampling mode
          // and also computes and store covariances using current covariance
          // matrix resampling mode

    void setObsMeanCov();

          // this version computes and returns the eigenvalues
          // of the covariance matrix

    void setObsMeanCov(RVector& coveigvals);


    void guessInitialFitParamValues(std::vector<double>& fitparams);

    void evalResiduals(const std::vector<double>& fitparams,
                       std::vector<double>& residuals) const;

    void evalResGradients(const std::vector<double>& fitparams,
                          RMatrix& gradients) const;

    double evalChiSquare(const std::vector<double>& residuals) const;
   

};


// *****************************************************************
#endif
