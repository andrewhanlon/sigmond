#ifndef CREATE_PLOTS_H
#define CREATE_PLOTS_H

#include "matrix.h"
#include "grace_plot.h"
#include "histogram.h"
#include "mcobs_info.h"
#include "scalar_defs.h"
#include "model_tcorr.h"

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


void createEffEnergyPlot(const std::vector<XYDYPoint>& meffvals,    
                         const ComplexArg& arg,
                         const std::string& correlator_name,
                         const std::string& filename, 
                         const std::string& symbol="circle", 
                         const std::string& symbolcolor="blue",
                         bool drawtoscreen=false);


             // goodnesstype='Q','X','N'  fit quality Q, chi-square per dof, none

void createEffEnergyPlotWithFit(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                const TCorrFitInfo& fitinfo,
                                char goodnesstype, double goodness,
                                const std::string& correlator_name,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen=false);


void createEffEnergyPlotWithFitAndEnergyRatio(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                const TCorrFitInfo& fitinfo,
                                double energy_ratio, double energy_ratio_err,
                                char goodnesstype, double goodness,
                                const std::string& correlator_name,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen=false);


void createEnergyDispersionPlot(const std::vector<XYDYPoint>& energy_sqs,
                                double anisotropy_mean, double anisotropy_err,
                                char goodnesstype, double goodness,
                                const std::string& particle_name,
                                const std::vector<XYPoint>& lowerfit,
                                const std::vector<XYPoint>& upperfit,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen=false);


void createCorrMatrixZMagSquaresPlot(const std::vector<XYDYPoint>& zmag_sqs,
                                     const std::string& observable_name,
                                     const std::string& filename,
                                     const std::string& barcolor, 
                                     bool drawtoscreen=false);


std::string getMCObsStandardName(const MCObsInfo& obs);

std::string getCorrelatorStandardName(const CorrelatorInfo& corr);

std::string getOpStandardName(const OperatorInfo& qcd_op);


// **************************************************************
#endif
