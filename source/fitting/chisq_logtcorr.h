#ifndef CHISQ_LOGTCORR_H
#define CHISQ_LOGTCORR_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "chisq_tcorr.h"

// *********************************************************************
// *                                                                   *
// *    The class "LogRealTemporalCorrelatorFit", derived from the     *
// *    base class "ChiSquare", is defined here.  It evaluates the     *
// *    chi^2 value associated with fits to a real-valued temporal     *
// *    correlator.  The input XML for constructing such an object     *
// *    is shown below:                                                *
// *                                                                   *
// *       <LogTemporalCorrelatorFit>                                  *
// *         <Operator>.... </Operator>                                *
// *         <SubtractVEV/>             (as appropriate)               *
// *         <MinimumTimeSeparation>3</MinimumTimeSeparation>          *
// *         <MaximumTimeSeparation>12</MaximumTimeSeparation>         *
// *         <ExcludeTimes>4 8</ExcludeTimes>  (optional)              *
// *         <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>          *
// *         <Model>                                                   *
// *           <LogAmplitude>                                          *
// *             <Name>log_amp</Name><IDIndex>0</IDIndex>              *
// *           </LogAmplitude>                                         *
// *           <Energy>                                                *
// *             <Name>energy</Name><IDIndex>0</IDIndex>               *
// *           </Energy>                                               *
// *         </Model>                                                  *
// *       </LogTemporalCorrelatorFit>                                 *
// *                                                                   *
// *    "LargeTimeNoiseCutoff" will lower the maximum time             *
// *    separation in the fit once the ratio of the error in the       *
// *    correlator over the correlator is smaller than                 *
// *    "LargeTimeNoiseCutoff".                                        *
// *                                                                   *
// *    The model used, assuming C(t) = A e^{-Et}, will be             *
// *    log C(t) = log A - E t. However, to keep this fit linear,      *
// *     the log A term will itself be a fit parameter (not A).        *
// *                                                                   *
// *********************************************************************



class LogRealTemporalCorrelatorFit :  public ChiSquare
{
    std::vector<uint> m_tvalues;
    uint T_period;
    OperatorInfo m_op;
    bool m_subt_vev;
    double m_noisecutoff;
    TimeForwardSingleExponential *m_single_exp;

 public:

    LogRealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~LogRealTemporalCorrelatorFit();

    uint getTmin() const {return m_tvalues.front();}

    uint getTmax() const {return m_tvalues.back();}
    
    const std::vector<uint>& getTvalues() const {return m_tvalues;}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;


    friend class TaskHandler;

};


// ************************************************************************
#endif
