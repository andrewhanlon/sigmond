#ifndef CHISQ_3PTCORR_H
#define CHISQ_3PTCORR_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "model_3ptcorr.h"

// *********************************************************************
// *                                                                   *
// *    The class "RealThreePointCorrelatorFit", derived from the      *
// *    base class "ChiSquare", is defined here.  It evaluates the     *
// *    chi^2 value associated with fits to a real-valued temporal     *
// *    correlator.  The input XML for constructing such an object     *
// *    is shown below:                                                *
// *                                                                   *
// *       <ThreePointCorrelatorFit>                                   *
// *         <Correlator>.... </Correlator>                            *
// *         <SubtractVEV/>             (as appropriate)               *
// *         <Figure out Time sep and ins stuff>
// *         <Model>...</Model>   (see "model_tcorr.h")                *
// *         <Summation/>      (optional)                              *
// *         <Priors>...</Priors>   (see "prior.h")                    *
// *       </ThreePointCorrelatorFit>                                  *
// *                                                                   *
// *    "LargeTimeNoiseCutoff" will lower the maximum time             *
// *    separation in the fit once the ratio of the error in the       *
// *    correlator over the correlator is smaller than                 *
// *    "LargeTimeNoiseCutoff".                                        *
// *                                                                   *
// *********************************************************************



class RealThreePointCorrelatorFit :  public ChiSquare
{
    CorrelatorInfo m_corr;
    OperatorInfo m_snk_op_2pt;
    OperatorInfo m_src_op_2pt;
    CorrelatorInfo m_rat_corr;
    bool m_subt_vev;
    ComplexArg m_arg;
    std::map<uint,std::set<uint> > m_tvalues;
    uint T_period;
    ThreePointCorrelatorModel *m_model_ptr;

 public:

    //RealThreePointCorrelatorFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    RealThreePointCorrelatorFit(
        MCObsHandler& OH, CorrelatorInfo in_corr,
        const OperatorInfo& in_snk_op_2pt, const OperatorInfo& in_src_op_2pt,
        CorrelatorInfo in_rat_corr,
        bool subtractvev, ComplexArg in_arg, std::string model_name,
        const std::map<std::string,MCObsInfo>& model_params, std::map<uint,std::set<uint> > in_fit_times);

    RealThreePointCorrelatorFit(
        MCObsHandler& OH, CorrelatorInfo in_corr, CorrelatorInfo in_rat_corr,
        bool subtractvev, ComplexArg in_arg, std::string model_name,
        const std::map<std::string,MCObsInfo>& model_params, std::map<uint,std::set<uint> > in_fit_times);

    RealThreePointCorrelatorFit(const RealThreePointCorrelatorFit& rtpcf);

    virtual ~RealThreePointCorrelatorFit();

    CorrelatorInfo getCorrelatorInfo() const {return m_corr;}

    OperatorInfo getTwoPointSinkOperatorInfo() const {return m_snk_op_2pt;}

    OperatorInfo getTwoPointSourceOperatorInfo() const {return m_src_op_2pt;}

    ComplexArg getComplexArg() const {return m_arg;}

    bool subtractVEV() const {return m_subt_vev;}

    const std::map<uint,std::set<uint> >& getTvalues() const {return m_tvalues;}

    ObsType getObsType() const {return m_model_ptr->getObsType();}

    virtual std::string getParameterName(uint param_index) const;

    double evalModelPoint(const std::vector<double>& fitparams, const double tsep, const double tins) const;

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;

 private:

    void setup();


    friend class TaskHandler;

};



// ************************************************************************
#endif
