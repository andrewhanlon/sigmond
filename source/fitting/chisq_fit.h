#ifndef CHISQ_FIT_H
#define CHISQ_FIT_H

#include <iostream>
#include <string>
#include "xml_handler.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "minimizer.h"
#include "chisq_tcorr.h"


// *****************************************************************
// *                                                               *
// *   The routine below performs the actual chi-square fitting.   *
// *   It carries out the minimization on the full sample, then    *
// *   loops over all resamplings (jackknife or bootstrap) to      *
// *   compute the uncertainties on the fit parameters.            *
// *   The best fit parameters are returned in "bestfit_params"    *
// *   and the chi-square value per degrees of freedom is          *
// *   returned in "chisq_dof" with the fit quality in "fitqual".  *
// *   The fit log is returned in XML format in "xmlout".          *
// *   An exception is thrown in any failure occurs.               *
// *                                                               *
// *****************************************************************

void doChiSquareFitting(ChiSquare& chisq_ref, 
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual, 
                        std::vector<MCEstimate>& bestfit_params,
                        XMLHandler& xmlout);

void doMultiSeriesFitting(XMLHandler& fit_xml, const int taskcount, RealMultiTemporalCorrelatorFit& chisq_ref, 
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual, 
                        std::vector<MCEstimate>& bestfit_params,
                        XMLHandler& xmlout, uint& final_tmin);

   // ************************************************************
   //
   //  The incomplete gamma function ratios are defined below:
   //
   //       double Qgamma(double s, double x);    s>0, x>=0
   //       double Pgamma(double s, double x);
   //
   //  If an error occurs, an exception is thrown.
   //
   //  Recall that the upper incomplete gamma function
   //  is defined by
   //  
   //    Gamma(s,x) = int( t^(s-1) exp(-t), t = x .. infinity)
   //  
   //  and the lower incomplete gamma function is given by
   //
   //    gamma(s,x) = int( t^(s-1) exp(-t), t = 0 .. x)
   //
   //  and the ordinary gamma function is
   //
   //    Gamma(s) = Gamma(s,0) =  int( t^(s-1) exp(-t), t = 0 .. infinity)
   //
   //  The incomplete gamma function ratios P and Q are ratios of the
   //  above functions:
   //
   //      Q(s,x) = Gamma(s,x)/Gamma(s) = 1 - P(s,x)
   //
   //      P(s,x) = gamma(s,x)/Gamma(s)
   //
   //  Note that  Q(s,0) = 1 and Q(a,infinity) = 0 and decreases 
   //  monotonically, whereas P(s,0) = 0  and  P(s,infinity) = 1 and 
   //  is monotonically increasing.  These functions get steeper
   //  as "s" decreases.
   //
   // ****************************************************************

double Qgamma(double s, double x);
double Pgamma(double s, double x);

   

   //  Returns the chi-square quality of fit, given by
   //
   //       Qgamma( dof/2,  chisquare/2 )


double getChiSquareFitQuality(unsigned int dof, double chisquare);


// ****************************************************************
#endif
