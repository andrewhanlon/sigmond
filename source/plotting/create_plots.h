#ifndef CREATE_PLOTS_H
#define CREATE_PLOTS_H

#include "matrix.h"
#include "grace_plot.h"
#include "histogram.h"
#include "mcobs_info.h"
#include "scalar_defs.h"
#include "model_tcorr.h"
#include "model_disp.h"
#include "plot_info.h"
#include "chisq_tcorr.h"
#include "chisq_disp.h"

void createEffEnergyPlotWithFit(EffEnergyWithFitPlotInfo plot_info,
                                RealTemporalCorrelatorFit& rtc,
                                FitResult& fit_result,
                                MCObsHandler* m_obs,
                                XMLHandler& xmlout);

void createDataFitRatioPlot(DataFitRatioPlotInfo plot_info,
                            std::vector<RealTemporalCorrelatorFit>& rtcs,
                            MCObsHandler* m_obs,
                            XMLHandler& xmlout);

void createTminPlot(TminFitPlotInfo plot_info,
                    std::vector<std::vector<RealTemporalCorrelatorFit> >& rtcs,
                    std::vector<std::vector<FitResult> >& fit_results,
                    MCObsHandler* m_obs,
                    XMLHandler& xmlout);

void createDispersionFitPlot(DispersionFitPlotInfo plot_info,
                             DispersionFit& disp_fit,
                             FitResult& fit_result,
                             MCObsHandler* m_obs,
                             XMLHandler& xmlout);

// *************************************************************

void createMCValuesPlot(const Vector<double>& mcvalues, 
                        const std::string& observable_name,
                        double mean_value, double std_dev,
                        const std::string& filename,
                        const std::string& symbol="circle", 
                        const std::string& symbolcolor="blue",
                        double rescale=1.0,
                        bool drawtoscreen=false);


void createMCBootstrapPlot(const Vector<double>& mcbootvalues, 
                           const std::string& observable_name,
                           double mean_value, double low, double upp, 
                           const std::string& filename,
                           const std::string& symbol="circle", 
                           const std::string& symbolcolor="blue",
                           double rescale=1.0,
                           bool drawtoscreen=false);


void createMCJackknifePlot(const Vector<double>& mcjackvalues, 
                           const std::string& observable_name,
                           double mean_value, 
                           const std::string& filename,
                           const std::string& symbol="circle", 
                           const std::string& symbolcolor="blue",
                           double rescale=1.0,
                           bool drawtoscreen=false);


void createMCHistogramPlot(const Histogram& histo, 
                           const std::string& observable_name,
                           double mean_value, double std_dev,
                           const std::string& filename, 
                           const std::string& barcolor="cyan", 
                           double rescale=1.0,
                           bool drawtoscreen=false);


void createMCBootstrapHistogramPlot(const Histogram& histo, 
                                    const std::string& observable_name,
                                    double mean_value, double low, double upp, 
                                    const std::string& filename, 
                                    const std::string& barcolor="cyan", 
                                    double rescale=1.0,
                                    bool drawtoscreen=false);


void createMCJackknifeHistogramPlot(const Histogram& histo, 
                                    const std::string& observable_name,
                                    double mean_value, 
                                    const std::string& filename, 
                                    const std::string& barcolor="cyan", 
                                    double rescale=1.0,
                                    bool drawtoscreen=false);


void createCorrelatorPlot(const std::vector<XYDYPoint>& corrvals,    
                          const ComplexArg& arg,
                          const std::string& correlator_name,
                          const std::string& filename, 
                          const std::string& symbol="circle", 
                          const std::string& symbolcolor="blue",
                          double rescale=1.0,
                          bool drawtoscreen=false);


void createCorrelatorPlot(const std::vector<XYDYPoint>& corrvals,    
                          const ComplexArg& arg,
                          const std::string& correlator_name,
                          const std::string& filename, 
                          double verticalmin, double verticalmax,
                          const std::string& symbol="circle", 
                          const std::string& symbolcolor="blue",
                          double rescale=1.0,
                          bool drawtoscreen=false);


void createEffEnergyPlot(const std::vector<XYDYPoint>& meffvals,    
                         const ComplexArg& arg,
                         const std::string& correlator_name,
                         const std::string& filename, 
                         const std::string& symbol="circle", 
                         const std::string& symbolcolor="blue",
                         bool drawtoscreen=false);


void createTminPlot(const std::vector<std::vector<XYDYDYPoint> >& goodfits,
                    const std::vector<std::vector<XYDYDYPoint> >& badfits,
                    const std::string& plotlabel,
                    const std::string& filename, 
                    uint energy_level,
                    const std::vector<std::string>& labels,
                    const std::vector<std::string>& symbols,
                    const std::vector<std::string>& colors,
                    const XYDYDYPoint& chosen_fit,
                    bool drawtoscreen=false);

             // goodnesstype='Q','X','N'  fit quality Q, chi-square per dof, none

void createEffEnergyPlotWithFit(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                double energy_mean, double energy_err,
                                uint tmin, uint tmax,
                                std::vector<XYPoint> meff_approach,
                                char goodnesstype, double goodness,
                                const std::string& correlator_name,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen=false);


void createEffEnergyPlotWithFitAndEnergyRatio(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                double energy_mean, double energy_err,
                                uint tmin, uint tmax,
                                std::vector<XYPoint> meff_approach,
                                double energy_ratio, double energy_ratio_err,
                                char goodnesstype, double goodness,
                                const std::string& correlator_name,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen=false);


void createDispersionFitPlot(const std::vector<XYDYPoint>& energy_sqs,
                             const std::vector<XYPoint>& lowerfit,
                             const std::vector<XYPoint>& upperfit,
                             const std::vector<double> param_means,
                             const std::vector<double> param_errs,
                             char goodnesstype, double goodness,
                             const std::string& particle_name,
                             const std::string& filename, 
                             const std::string& symbol, 
                             const std::string& symbolcolor,
                             bool drawtoscreen=false);


void createCorrMatrixZMagSquaresPlot(const std::vector<XYDYPoint>& zmag_sqs,
                                     const std::string& observable_name,
                                     const std::string& filename,
                                     const std::string& barcolor, 
                                     bool drawtoscreen=false);

void createDataFitRatioPlot(const std::vector<std::vector<XYDYPoint> >& ratiovals,
                            const std::string& plotlabel,
                            const std::string& filename, 
                            const std::vector<std::string>& labels, 
                            const std::vector<std::string>& symbols, 
                            const std::vector<std::string>& symbolcolors,
                            bool drawtoscreen=false);


std::string getMCObsStandardName(const MCObsInfo& obs);

std::string getCorrelatorStandardName(const CorrelatorInfo& corr);

std::string getOpStandardName(const OperatorInfo& qcd_op);

    //  If "snkname" or "srcname" is "standard", use the standard name,
    //  else use the label given.
std::string getCorrelatorName(const CorrelatorInfo& corr, 
               const std::string& snkname, const std::string& srcname);

// **************************************************************
#endif
