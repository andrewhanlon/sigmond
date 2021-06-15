#ifndef CHISQ_DISP_H
#define CHISQ_DISP_H

#include <algorithm>
#include "xml_handler.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "model_disp.h"

// *********************************************************************
// *                                                                   *
// *    The class "DispersionFit", derived from the                    *
// *    base class "ChiSquare", is defined here.  It evaluates the     *
// *    chi^2 value associated with fits to a dispersion relation      *
// *    The input XML for constructing such an object                  *
// *    is shown below:                                                *
// *                                                                   *
// *       <DispersionFit>                                             *
// *         <Energy>                                                  *
// *           <Name>pion</Name><IDIndex>0</IDIndex>                   *
// *           <Momentum>0 0 0</Momentum>                              *
// *         </Energy>                                                 *
// *         <Energy>                                                  *
// *           <Name>pion</Name><IDIndex>1</IDIndex>                   *
// *           <Momentum>0 0 1</Momentum>                              *
// *         </Energy>                                                 *
// *         <Energy>... </Energy>                                     *
// *         <Model>...</Model>   (see "model_disp.h")                 *
// *         <Priors>...</Priors>   (see "prior.h")                    *
// *       </DispersionFit>                                            *
// *                                                                   *
// *********************************************************************



class DispersionFit :  public ChiSquare
{
    uint Xextent, Yextent, Zextent;
    std::vector<double> m_momsq;
    DispersionModel *m_model_ptr;

 public:

    DispersionFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    DispersionFit(MCObsHandler& OH, const std::string& model_name, 
        const std::map<MCObsInfo,double>& energies,
        const std::map<std::string,MCObsInfo>& model_params);

    virtual ~DispersionFit();

    const std::vector<double>& getMomSqs() const { return m_momsq; }

    double getMinMomSq() const { return *std::min_element(m_momsq.begin(), m_momsq.end()); }

    double getMaxMomSq() const { return *std::max_element(m_momsq.begin(), m_momsq.end()); }

    std::string getParameterName(uint param_index) const;

    double evalDispersion(const double msq, const double psq) const;

    double evalModelPoint(const std::vector<double>& fitparams,
                                const double psq) const;

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;

    MCObsInfo getMomSqObsParamInfo(const double psq) const;

 private:

    void setup(const std::map<MCObsInfo,double>& energies);

    friend class TaskHandler;

};




// ************************************************************************
#endif
