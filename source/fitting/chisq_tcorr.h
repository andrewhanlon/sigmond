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
    

 public:
 
    TemporalCorrelatorModel *m_model_ptr;

    RealTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~RealTemporalCorrelatorFit();

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


// ************************************************************************

class RealMultiTemporalCorrelatorFit :  public ChiSquare
{
    std::vector<uint> m_tvalues;
    std::vector<uint> m_tvalues_tot;
    std::vector<MCObsInfo> m_obs_info_tot; 
    uint T_period;
    OperatorInfo m_op;
    bool m_subt_vev;
    double m_noisecutoff;
    uint m_max_level;
    double m_initial_gap;
    double m_repeating_gaps;
    std::vector<bool> m_is_prior_param;
    std::vector<double> m_priors;
    std::vector<double> m_prior_range;

 public:
 
//     const XMLHandler xml_input;
    TemporalCorrelatorModel *m_model_ptr;

    RealMultiTemporalCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~RealMultiTemporalCorrelatorFit();

    uint getTmin() const {return m_tvalues.front();}

    uint getTmax() const {return m_tvalues.back();}
    
    const std::vector<uint>& getTvalues() const {return m_tvalues;}
    const std::vector<uint>& getTvalues_tot() const {return m_tvalues_tot;}

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;
    
    
//     void evalResiduals(const std::vector<double>& fitparams,
//                        std::vector<double>& residuals) const;

//     void evalResGradients(const std::vector<double>& fitparams,
//                           RMatrix& gradients) const;
    
    void pop_front_tlist();
    void reset_tlist( const std::vector<uint>& m_tvalues_tot, const std::vector<MCObsInfo>& m_obs_info_tot );
    
    uint getMaxLevel() const {return m_max_level;}
    double getInitialGap() const {return m_initial_gap;}
    double getRepeatingGaps() const {return m_repeating_gaps;}
    const std::vector<bool>& getPriorParams() const {return m_is_prior_param;}
    const std::vector<double>& getPriors() const {return m_priors;}
    const std::vector<double>& getPriorRanges() const {return m_prior_range;}

    friend class TaskHandler;

};
#endif
