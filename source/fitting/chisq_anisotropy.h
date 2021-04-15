#ifndef CHISQ_ANISOTROPY_H
#define CHISQ_ANISOTROPY_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "chisq_base.h"

// **********************************************************************
// *                                                                    *
// *    The class "AnisotropyFromDispersionFit", derived from the       *
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
// *       <AnisotropyFromDispersionFit>                                *
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
// *       </AnisotropyFromDispersionFit>                               *
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



class AnisotropyFromDispersionFit :  public ChiSquare
{
    static const std::vector<std::string> param_names;

    uint m_lat_spatial_extent;
    double m_momsq_quantum;
    std::vector<uint> m_imomsq;

 public:

    AnisotropyFromDispersionFit(XMLHandler& xmlin, MCObsHandler& OH, int taskcount);

    virtual ~AnisotropyFromDispersionFit();

    virtual void evalModelPoints(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;

    virtual void evalGradients(const std::vector<double>& fitparams,
                               RMatrix& gradients) const;

    virtual void guessInitialParamValues(const RVector& datapoints,
                                         std::vector<double>& fitparams) const;

    virtual void do_output(XMLHandler& xmlout) const;

    MCObsInfo getAnisotropyKey() const {return m_fitparam_info[0];}

    MCObsInfo getRestMassSquaredKey() const {return m_fitparam_info[1];}

    friend class TaskHandler;

};


// **************************************************************************
#endif
