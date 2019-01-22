#ifndef CHISQ_CORRINTRATIO_H
#define CHISQ_CORRINTRATIO_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"

// **********************************************************************
// *                                                                    *
// *      For extracting the `interaction energy' or energy difference           *
// *      between an interacting energy level and the nearest                    *
// *      non-interacting state, the ratio:                                      *
// *                               C_int(t)                                      *
// *               R(t) =  -------------------------                             *
// *                         prod_i C_non-int[i](t)                              *
// *      can be fit to the ansatz:  R(t) = A exp(-DeltaE*t), where:             *
// *        DeltaE          =  E_int - E_non-int,                                *
// *        C_int(t)        =  (diagonal) correlator for interacting level       *
// *        C_non-int[i](t) =  N correlators coresponding to closest             *
// *                           non-interacting level                             *
// *                                                                             *
// *      Example:                                                               *
// *      if the closest non-interacting level to C_int is pi(0)K(1)K(1):        *
// *        C_non-int[0](t) = pion correlator for P^2=0                          *
// *        C_non-int[1](t) = kaon correlator for P^2=1                          *
// *        C_non-int[2](t) = kaon correlator for P^2=1                          *
// *                                                                             *
// *      Note: the resultant correlator will have (if specified) all            *
// *            VEVs subtracted already.  Information about the ratio
// *            of correlators must be specified in an MCObsInfo object
// *            (which must be real and nonsimple).                             *



// *    The class "CorrelatorInteractionRatioFit", derived from the       *
// *    base class "ChiSquare", is defined here.  It evaluates the      *
// *    chi^2 value associated with fits to the dispersion relation     *
// *    of an individual free particle.  The observables are the        *
// *    free particle energies for various three-momenta squared,       *
// *    and there are two fit parameters: the rest mass of the particle *
// *    times a_t, and the lattice anisotropy xi = a_s/a_t, where a_s   *
// *    is the spatial lattice spacing and a_t is the temporal lattice  *
// *    spacing. The input XML for constructing such an object          *
// *    is shown below:                                                 *
// *                                                                    *
// *       <CorrelatorInteractionRatioFit>                                *
// *         <SpatialExtentNumSites>32</SpatialExtentNumSites>          *
// *         <Energy>                                                   *
// *           <Name>pion</Name><IDIndex>0</IDIndex>                    *
// *           <IntMomSquared>0</IntMomSquared>                         *
// *         </Energy>                                                  *
// *         <Energy>                                                   *
// *           <Name>pion</Name><IDIndex>1</IDIndex>                    *
// *           <IntMomSquared>1</IntMomSquared>                         *
// *         </Energy>                                                  *
// *         <Energy>... </Energy>                                      *
// *         <Anisotropy>                                               *
// *           <Name>PionXi</Name><IDIndex>0</IDIndex>                  *
// *         </Anisotropy>                                              *
// *         <RestMassSquared>                                          *
// *           <Name>PionRestMassSq</Name><IDIndex>0</IDIndex>          *
// *         </RestMassSquared>                                         *
// *       </CorrelatorInteractionRatioFit>                               *
// *                                                                    *
// *    The model used for the observables is                           *
// *                                                                    *
// *      (a_t E)^2 = restmass_sq +  (2*Pi/Ns)^2 * nsq / xi^2           *
// *                                                                    *
// *    where "restmass_sq" and "xi" are the two model parameters, and  *
// *    "Ns" is extent of the lattice in terms of number of sites       *
// *    in each of the three spatial directions, and "nsq" is the       *
// *    integer square of the three momentum.  Recall that the          *
// *    (a_s P)^2 = (2*Pi/Ns)^2 * nsq.                                  *
// *                                                                    *
// *                    xi = fitparams[0]                               *
// *           restmass_sq = fitparams[1]                               *
// *                                                                    *
// **********************************************************************



class CorrelatorInteractionRatioFit :  public ChiSquare
{
    std::vector<uint> m_tvalues;
    uint T_period;
    MCObsInfo m_corrratioinfo;
    double m_noisecutoff;
    TemporalCorrelatorModel *m_model_ptr;

 public:

    CorrelatorInteractionRatioFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~CorrelatorInteractionRatioFit();

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    uint getTmin() const {return m_tvalues.front();}

    uint getTmax() const {return m_tvalues.back();}

    virtual void do_output(XMLHandler& xmlout) const;

    friend class TaskHandler;

};


// **************************************************************************
#endif
