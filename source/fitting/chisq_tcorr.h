#ifndef CHISQ_TCORR_H
#define CHISQ_TCORR_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "model_tcorr.h"

// *********************************************************************
// *                                                                   *
// *    The class "RealTemporalCorrelatorFit", derived from the        *
// *    base class "ChiSquare", is defined here.  It evaluates the     *
// *    chi^2 value associated with fits to a real-valued temporal     *
// *    correlator.  The input XML for constructing such an object     *
// *    is shown below:                                                *
// *                                                                   *
// *       <TemporalCorrelatorFit>                                     *
// *         <Operator>.... </Operator>                                *
// *         <SubtractVEV/>             (as appropriate)               *
// *         <MinimumTimeSeparation>3</MinimumTimeSeparation>          *
// *         <MaximumTimeSeparation>12</MaximumTimeSeparation>         *
// *         <ExcludeTimes>4 8</ExcludeTimes>  (optional)              *
// *         <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>          *
// *         <Model>...</Model>   (see "model_tcorr.h")                *
// *       </TemporalCorrelatorFit>                                    *
// *                                                                   *
// *    "LargeTimeNoiseCutoff" will lower the maximum time             *
// *    separation in the fit once the ratio of the error in the       *
// *    correlator over the correlator is smaller than                 *
// *    "LargeTimeNoiseCutoff".                                        *
// *                                                                   *
// *********************************************************************



class RealTemporalCorrelatorFit :  public ChiSquare
{
    std::vector<uint> m_tvalues;
    uint T_period;
    OperatorInfo m_op;
    bool m_subt_vev;
    double m_noisecutoff;
    TemporalCorrelatorModel *m_model_ptr;

 public:

    RealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~RealTemporalCorrelatorFit();

    uint getTmin() const {return m_tvalues.front();}

    uint getTmax() const {return m_tvalues.back();}
    
    uint getNumberOfParams() const {return m_model_ptr->getNumberOfParams();}
    
    const std::vector<uint>& getTvalues() const {return m_tvalues;}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;
    
    virtual void plot(const XMLHandler& xmlin, const int taskcount, const double qual, const double goodness, 
                      const uint lat_time_extent, const std::vector<MCEstimate>& bestfit_params, XMLHandler& xmlout) const;

    friend class TaskHandler;

};



// *********************************************************************
// *                                                                   *
// *    The class "TwoRealTemporalCorrelatorFit", derived from the     *
// *    base class "ChiSquare", is defined here.  It evaluates the     *
// *    chi^2 value associated with a fit to two real-valued temporal  *
// *    correlators.  The input XML for constructing such an object    *
// *    is shown below:                                                *
// *                                                                   *
// *       <TwoTemporalCorrelatorFit>                                  *
// *         <CorrelatorOne>                                           *
// *           <Operator>.... </Operator>                              *
// *           <SubtractVEV/>             (as appropriate)             *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>        *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>       *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)            *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>        *
// *           <Model>...</Model>                                      *
// *         </CorrelatorOne>                                          *
// *         <CorrelatorTwo>                                           *
// *           <Operator>.... </Operator>                              *
// *           <SubtractVEV/>             (as appropriate)             *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>        *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>       *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)            *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>        *
// *           <Model>...</Model>                                      *
// *         </CorrelatorTwo>                                          *
// *         <EnergyRatio>                                             *
// *            <Name>pion</Name><IDIndex>0</IDIndex>                  *
// *         </EnergyRatio>                                            *
// *       </TwoTemporalCorrelatorFit>                                 *
// *                                                                   *
// *                                                                   *
// *********************************************************************



class TwoRealTemporalCorrelatorFit :  public ChiSquare
{
    std::vector<uint> m_tvalues1, m_tvalues2;
    uint T_period;
    OperatorInfo m_op1, m_op2;
    bool m_subt_vev1, m_subt_vev2;
    double m_noisecutoff1, m_noisecutoff2;
    TemporalCorrelatorModel *m_model1_ptr, *m_model2_ptr;
    MCObsInfo m_energyratio;

 public:

    TwoRealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~TwoRealTemporalCorrelatorFit();

    uint getTmin1() const {return m_tvalues1.front();}

    uint getTmax1() const {return m_tvalues1.back();}

    uint getTmin2() const {return m_tvalues2.front();}

    uint getTmax2() const {return m_tvalues2.back();}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;

    friend class TaskHandler;

};

class NSimRealTemporalCorrelatorFit :  public ChiSquare
{
    std::vector<RealTemporalCorrelatorFit*> m_fits;
//     std::vector<std::vector<uint>> m_tvalues;
//     uint T_period;
//     std::vector<OperatorInfo> m_op;
//     std::vector<bool> m_subt_vev;
//     std::vector<double> m_noisecutoff;
//     std::vector<TemporalCorrelatorModel> *m_model_ptrs;
//     MCObsInfo m_energyratio;

 public:

    NSimRealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~NSimRealTemporalCorrelatorFit();

//     uint getTmin(uint i) const {return m_tvalues[i].front();}

//     uint getTmax(uint i) const {return m_tvalues[i].back();}

    uint getFitNum() const {return m_fits.size();}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;
    
    virtual void plot(const uint index, const XMLHandler& xmlin, const int taskcount, 
                        const double qual, const double goodness, const uint lat_time_extent, 
                        const std::vector<MCEstimate>& bestfit_params, XMLHandler& xmlout) const;

    friend class TaskHandler;

};


// ************************************************************************
#endif
