#ifndef CHISQ_TCORR_H
#define CHISQ_TCORR_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "model_tcorr.h"

// *****************************************************************


// *       <TemporalCorrelatorFit>
// *         <Operator>.... </Operator>
// *         <SubtractVEV/>             (as appropriate)
// *         <MinimumTimeSeparation>3</MinimumTimeSeparation>
// *         <MaximumTimeSeparation>12</MaximumTimeSeparation>
// *         <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>
// *         <Model>...</Model>   (see "model_tcorr.h")
// *       </TemporalCorrelatorFit>
// *
// *    "LargeTimeNoiseCutoff" will lower the maximum time
// *    separation in the fit once the ratio of the error in the 
// *    correlator over the correlator is smaller than
// *    "LargeTimeNoiseCutoff".



class RealTemporalCorrelatorFit :  public ChiSquare
{
    uint m_tmin, m_tmax, T_period;
    OperatorInfo m_op;
    bool m_subt_vev;
    double m_noisecutoff;
    TemporalCorrelatorModel *m_model_ptr;

 public:

    RealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~RealTemporalCorrelatorFit();

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;


    friend class TaskHandler;

};



// *****************************************************************

/*
// *       <TwoTemporalCorrelatorFit>
// *         <CorrelatorRef>
// *           <Operator>.... </Operator>
// *           <SubtractVEV/>             (as appropriate)
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>
// *           <Model>...</Model>
// *         </CorrelatorRef>
// *         <Correlator>
// *           <Operator>.... </Operator>
// *           <SubtractVEV/>             (as appropriate)
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>
// *           <Model>...</Model>
// *         </Correlator>
// *       </TwoTemporalCorrelatorFit>



class TwoRealTemporalCorrelatorFit :  public ChiSquare
{
    uint m_tmin, m_tmax, mref_tmin, mref_tmax, T_period;
    OperatorInfo m_op, mref_op;
    bool m_subt_vev, mref_subt_vev;
    double m_noisecutoff, mref_noisecutoff;
    TemporalCorrelatorModel *m_model_ptr, *mref_model_ptr;

 public:

    TwoRealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH);

    virtual ~TwoRealTemporalCorrelatorFit(){}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;

    friend class TaskHandler;

};

*/
// *************************************************************
#endif
