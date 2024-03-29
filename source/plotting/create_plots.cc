#include "create_plots.h"
#include "xml_handler.h"
#include "mc_estimate.h"

using namespace std;

// *************************************************************


void createMCValuesPlot(const Vector<double>& mcvalues, const string& observable_name,
                        double in_mean_value, double in_std_dev,
                        const string& filename, const string& symbol, 
                        const string& symbolcolor, double rescale, bool drawtoscreen)
{
 GracePlot P("Markov Chain Index","Value");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double std_dev=rescale*in_std_dev;

 P.addXYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<mcvalues.size();ind++)
    P.addXYDataPoint(double(ind),rescale*mcvalues[ind]);

 P.addXYDataSet("none","none","solid",symbolcolor);
 P.addXYDataPoint(0,mean_value); P.addXYDataPoint(mcvalues.size(),mean_value);
 P.addXYDataSet("none","none","dash",symbolcolor);
 P.addXYDataPoint(0,mean_value+std_dev); P.addXYDataPoint(mcvalues.size(),mean_value+std_dev);
 P.addXYDataSet("none","none","dash",symbolcolor);
 P.addXYDataPoint(0,mean_value-std_dev); P.addXYDataPoint(mcvalues.size(),mean_value-std_dev);
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createMCBootstrapPlot(const Vector<double>& bootvals, const string& observable_name,
                           double in_mean_value, double in_low, double in_upp, 
                           const string& filename, const string& symbol, 
                           const string& symbolcolor, double rescale, bool drawtoscreen)
{
 GracePlot P("Bootstrap Index","Value");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double low=rescale*in_low;
 double upp=rescale*in_upp;

 P.addXYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<bootvals.size();ind++)
    P.addXYDataPoint(double(ind),rescale*bootvals[ind]);

 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(0,mean_value); P.addXYDataPoint(bootvals.size(),mean_value);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0,low); P.addXYDataPoint(bootvals.size(),low);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0,upp); P.addXYDataPoint(bootvals.size(),upp);
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createMCJackknifePlot(const Vector<double>& jackvals, const string& observable_name,
                           double in_mean_value, const string& filename, const string& symbol, 
                           const string& symbolcolor, double rescale, bool drawtoscreen)
{
 GracePlot P("Jackknife Index","Value");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;

 P.addXYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<jackvals.size();ind++)
    P.addXYDataPoint(double(ind),rescale*jackvals[ind]);

 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(0,mean_value); P.addXYDataPoint(jackvals.size(),mean_value);
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createMCHistogramPlot(const Histogram& histo, const string& observable_name,
                           double in_mean_value, double in_std_dev,
                           const string& filename, const string& barcolor, 
                           double rescale, bool drawtoscreen)
{
 GracePlot P("MC Value","Number");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double std_dev=rescale*in_std_dev;

 P.addBarDataSet(barcolor,"black",rescale*histo.getBarWidth());
 for (uint ind=0;ind<histo.getNumberOfBars();ind++)
    P.addBarDataPoint(rescale*histo.getBarMiddleLocation(ind),histo.getBarHeight(ind));

 double hmax=double(histo.getMaxHeight())*1.3;
 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(mean_value,0.0); P.addXYDataPoint(mean_value,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(mean_value+std_dev,0.0); P.addXYDataPoint(mean_value+std_dev,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(mean_value-std_dev,0.0); P.addXYDataPoint(mean_value-std_dev,hmax);
 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}


void createMCBootstrapHistogramPlot(const Histogram& histo, const string& observable_name,
                                    double in_mean_value, double in_low, double in_upp, 
                                    const string& filename, const string& barcolor, 
                                    double rescale, bool drawtoscreen)
{
 GracePlot P("Bootstrap Value","Number");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double low=rescale*in_low;
 double upp=rescale*in_upp; 

 P.addBarDataSet(barcolor,"black",rescale*histo.getBarWidth());
 for (uint ind=0;ind<histo.getNumberOfBars();ind++)
    P.addBarDataPoint(rescale*histo.getBarMiddleLocation(ind),histo.getBarHeight(ind));

 double hmax=double(histo.getMaxHeight())*1.3;
 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(mean_value,0.0); P.addXYDataPoint(mean_value,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(upp,0.0); P.addXYDataPoint(upp,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(low,0.0); P.addXYDataPoint(low,hmax);
 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createMCJackknifeHistogramPlot(const Histogram& histo, const string& observable_name,
                                    double in_mean_value, 
                                    const string& filename, const string& barcolor, 
                                    double rescale, bool drawtoscreen)
{
 GracePlot P("Jackknife Value","Number");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;

 P.addBarDataSet(barcolor,"black",rescale*histo.getBarWidth());
 for (uint ind=0;ind<histo.getNumberOfBars();ind++)
    P.addBarDataPoint(rescale*histo.getBarMiddleLocation(ind),histo.getBarHeight(ind));

 double hmax=double(histo.getMaxHeight())*1.3;
 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(mean_value,0.0); P.addXYDataPoint(mean_value,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}


void createCorrelatorPlot(const std::vector<XYDYPoint>& corrvals,
                          const ComplexArg& arg,
                          const std::string& correlator_name,
                          const std::string& filename, 
                          const std::string& symbol, 
                          const std::string& symbolcolor,
                          double rescale, bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t",prefix+" C(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 XYDYPoint cval;
 for (uint ind=0;ind<corrvals.size();ind++){
    cval=corrvals[ind]; cval.yval*=rescale; cval.yerr*=rescale;
    P.addXYDYDataPoint(cval);
    if (corrvals[ind].xval>tmax) tmax=corrvals[ind].xval;}
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0.0,0.0); P.addXYDataPoint(tmax+5.0,0.0);

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createCorrelatorPlot(const std::vector<XYDYPoint>& corrvals,
                          const ComplexArg& arg,
                          const std::string& correlator_name,
                          const std::string& filename, 
                          double verticalmin, double verticalmax,
                          const std::string& symbol, 
                          const std::string& symbolcolor,
                          double rescale, bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t",prefix+" C(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 XYDYPoint cval;
 for (uint ind=0;ind<corrvals.size();ind++){
    cval=corrvals[ind]; cval.yval*=rescale; cval.yerr*=rescale;
    P.addXYDYDataPoint(cval);
    if (corrvals[ind].xval>tmax) tmax=corrvals[ind].xval;}
 P.autoScale(0.02,0.02,0.2,0.2);
 P.setVerticalLimits(verticalmin,verticalmax);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0.0,0.0); P.addXYDataPoint(tmax+5.0,0.0);

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createEffEnergyPlot(const std::vector<XYDYPoint>& meffvals,
                         const ComplexArg& arg,
                         const std::string& correlator_name,
                         const std::string& filename, 
                         const std::string& symbol, 
                         const std::string& symbolcolor,
                         bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t","a\\st\\NE\\s\\f{0}eff\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 for (uint ind=0;ind<meffvals.size();ind++){
    P.addXYDYDataPoint(meffvals[ind]);
    if (meffvals[ind].xval>tmax) tmax=meffvals[ind].xval;}
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createEffEnergyPlotWithFit(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                const TCorrFitInfo& fitinfo,
                                char goodnesstype, double goodness,
                                const std::string& correlator_name,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t","a\\st\\NE\\s\\f{0}eff\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 for (uint ind=0;ind<meffvals.size();ind++){
    P.addXYDYDataPoint(meffvals[ind]);
    if (meffvals[ind].xval>tmax) tmax=meffvals[ind].xval;}

 double fitupper=fitinfo.energy_mean+fitinfo.energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fitinfo.tmin,fitupper);
 P.addXYDataPoint(fitinfo.tmax,fitupper);
 double fitlower=fitinfo.energy_mean-fitinfo.energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fitinfo.tmin,fitlower);
 P.addXYDataPoint(fitinfo.tmax,fitlower);

 if (!(fitinfo.meff_approach.empty())){
    P.addXYDataSet("none","open","dash",symbolcolor);
    P.addXYDataPoints(fitinfo.meff_approach);}

 SimpleMCEstimate fitres(fitinfo.energy_mean,fitinfo.energy_err);
 string fitenergy("\\f{1}a\\st\\NE\\f{}\\sfit\\N = ");
 fitenergy+=fitres.str(2);
 P.addText(fitenergy,0.90,0.85,true,1.7,"black","top-right");

 if (goodnesstype=='Q'){
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.80,true,0,"black","top-right");}
 else if (goodnesstype=='X'){
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.80,true,0,"black","top-right");}

 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createTMinPlot(const std::vector<XYDYDYPoint>& goodcorrelatedfits,
                    const std::vector<XYDYDYPoint>& gooduncorrelatedfits,
                    const std::vector<XYDYDYPoint>& badcorrelatedfits,
                    const std::vector<XYDYDYPoint>& baduncorrelatedfits,
                    const std::string& observable_name,
                    const std::string& filename, 
                    const std::string& symbol,
                    const std::string& goodfitcolor,
                    const std::string& badfitcolor,
                    bool correlatedfit_hollow, bool uncorrelatedfit_hollow,
                    const XYDYDYPoint& chosen_fit,
                    bool drawtoscreen, bool tmax)
{
 string xlabel = "t\\s\\f{0}min\\f{}";
 if( tmax ) xlabel = "t\\s\\f{0}max\\f{}";
 GracePlot P(xlabel,"a\\st\\NE\\s\\f{0}fit\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);
 double xmax=-10.0;
 double xmin=1e32;
 string correlatedfill=(correlatedfit_hollow) ? "open" : "solid";
 P.addXYDYDYDataSet(symbol,correlatedfill,"none",goodfitcolor);
 for (uint ind=0;ind<goodcorrelatedfits.size();ind++){
    double x=goodcorrelatedfits[ind].xval;
    if (x<xmin) xmin=x;
    if (x>xmax) xmax=x;
    P.addXYDYDYDataPoint(goodcorrelatedfits[ind]);}
 P.addXYDYDYDataSet(symbol,correlatedfill,"none",badfitcolor);
 for (uint ind=0;ind<badcorrelatedfits.size();ind++){
    double x=badcorrelatedfits[ind].xval;
    if (x<xmin) xmin=x;
    if (x>xmax) xmax=x;
    P.addXYDYDYDataPoint(badcorrelatedfits[ind]);}
 string uncorrelatedfill=(uncorrelatedfit_hollow) ? "open" : "solid";
 P.addXYDYDYDataSet(symbol,uncorrelatedfill,"none",goodfitcolor);
 for (uint ind=0;ind<gooduncorrelatedfits.size();ind++){
    double x=gooduncorrelatedfits[ind].xval;
    if (x<xmin) xmin=x;
    if (x>xmax) xmax=x;
    P.addXYDYDYDataPoint(gooduncorrelatedfits[ind]);}
 P.addXYDYDYDataSet(symbol,uncorrelatedfill,"none",badfitcolor);
 for (uint ind=0;ind<baduncorrelatedfits.size();ind++){
    double x=baduncorrelatedfits[ind].xval;
    if (x<xmin) xmin=x;
    if (x>xmax) xmax=x;
    P.addXYDYDYDataPoint(baduncorrelatedfits[ind]);}
 if (chosen_fit.xval>0.0){
    double chosenupper=chosen_fit.yval+chosen_fit.yuperr;
    P.addXYDataSet("none","open","dash","black");
    P.addXYDataPoint(xmin,chosenupper);
    P.addXYDataPoint(xmax,chosenupper);
    double chosen=chosen_fit.yval;
    P.addXYDataSet("none","solid","solid","black");
    P.addXYDataPoint(xmin,chosen);
    P.addXYDataPoint(xmax,chosen);
    double chosenlower=chosen_fit.yval-chosen_fit.ydnerr;
    P.addXYDataSet("none","open","dash","black");
    P.addXYDataPoint(xmin,chosenlower);
    P.addXYDataPoint(xmax,chosenlower);}
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}


void createEffEnergyPlotWithFitAndEnergyRatio(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                const TCorrFitInfo& fitinfo,
                                double energy_ratio, double energy_ratio_err,
                                char goodnesstype, double goodness,
                                const std::string& correlator_name,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t","a\\st\\NE\\s\\f{0}eff\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 for (uint ind=0;ind<meffvals.size();ind++){
    P.addXYDYDataPoint(meffvals[ind]);
    if (meffvals[ind].xval>tmax) tmax=meffvals[ind].xval;}

 double fitupper=fitinfo.energy_mean+fitinfo.energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fitinfo.tmin,fitupper);
 P.addXYDataPoint(fitinfo.tmax,fitupper);
 double fitlower=fitinfo.energy_mean-fitinfo.energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fitinfo.tmin,fitlower);
 P.addXYDataPoint(fitinfo.tmax,fitlower);

 if (!(fitinfo.meff_approach.empty())){
    P.addXYDataSet("none","open","dash",symbolcolor);
    P.addXYDataPoints(fitinfo.meff_approach);}

 SimpleMCEstimate fitres(fitinfo.energy_mean,fitinfo.energy_err);
 string fitenergy("\\f{1}a\\st\\NE\\f{}\\sfit\\N = ");
 fitenergy+=fitres.str(2);
 P.addText(fitenergy,0.90,0.9,true,1.7,"black","top-right");

 SimpleMCEstimate ratio(energy_ratio,energy_ratio_err);
 string ratiostr("\\f{1}E\\f{}\\sfit\\N\\f{1}/E\\f{}\\sref\\N = ");
 ratiostr+=ratio.str(2);
 P.addText(ratiostr,0.25,0.18,true,1.7,"black","bottom-left");

 if (goodnesstype=='Q'){
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.85,true,0,"black","top-right");}
 else if (goodnesstype=='X'){
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.85,true,0,"black","top-right");}

 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createAnisoEnergyDispersionPlot(const std::vector<XYDYPoint>& energy_sqs,
                                double anisotropy_mean, double anisotropy_err,
                                char goodnesstype, double goodness,
                                const std::string& particle_name,
                                const std::vector<XYPoint>& lowerfit,
                                const std::vector<XYPoint>& upperfit,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen)
{
 GracePlot P("n\\m{0}\\S2\\N\\M{0}\\s\\f{3}P\\N","a\\st\\N\\S2\\NE\\s\\f{0}fit\\N\\S2\\N\\f{}");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.25,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<energy_sqs.size();ind++){
    P.addXYDYDataPoint(energy_sqs[ind]);} 

 P.addXYDataSet("none","none","solid",symbolcolor);
 for (uint ind=0;ind<lowerfit.size();ind++)
    P.addXYDataPoint(lowerfit[ind]);
 P.addXYDataSet("none","none","solid",symbolcolor);
 for (uint ind=0;ind<upperfit.size();ind++)
    P.addXYDataPoint(upperfit[ind]);

 SimpleMCEstimate xires(anisotropy_mean,anisotropy_err);
 string xifit("\\f{1}a\\ss\\N/a\\st\\N\\f{} = ");
 xifit+=xires.str(2);
 P.addText(xifit,0.90,0.23,true,1.7,"black","bottom-right");

 if (goodnesstype=='Q'){
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.19,true,0,"black","bottom-right");}
 else if (goodnesstype=='X'){
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.19,true,0,"black","bottom-right");}

 P.autoScale(0.02,0.02,0.2,0.2);
 if (!particle_name.empty())
    P.addText(particle_name,0.30,0.9,true,1.5,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
} 


void createCoeffEnergyDispersionPlot(const std::vector<XYDYPoint>& energy_sqs,
                                double coeff_mean, double coeff_err,
                                char goodnesstype, double goodness,
                                const std::string& particle_name,
                                const std::vector<XYPoint>& lowerfit,
                                const std::vector<XYPoint>& upperfit,
                                const std::string& filename, 
                                const std::string& symbol, 
                                const std::string& symbolcolor,
                                bool drawtoscreen)
{
 GracePlot P("n\\m{0}\\S2\\N\\M{0}\\s\\f{3}P\\N","a\\st\\N\\S2\\NE\\s\\f{0}fit\\N\\S2\\N\\f{}");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.25,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<energy_sqs.size();ind++){
    P.addXYDYDataPoint(energy_sqs[ind]);} 

 P.addXYDataSet("none","none","solid",symbolcolor);
 for (uint ind=0;ind<lowerfit.size();ind++)
    P.addXYDataPoint(lowerfit[ind]);
 P.addXYDataSet("none","none","solid",symbolcolor);
 for (uint ind=0;ind<upperfit.size();ind++)
    P.addXYDataPoint(upperfit[ind]);

 SimpleMCEstimate cres(coeff_mean,coeff_err);
 string cfit("\\f{1}c\\f{} = ");
 cfit+=cres.str(2);
 P.addText(cfit,0.90,0.23,true,1.7,"black","bottom-right");

 if (goodnesstype=='Q'){
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.19,true,0,"black","bottom-right");}
 else if (goodnesstype=='X'){
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.19,true,0,"black","bottom-right");}

 P.autoScale(0.02,0.02,0.2,0.2);
 if (!particle_name.empty())
    P.addText(particle_name,0.30,0.9,true,1.5,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
} 



void createCorrMatrixZMagSquaresPlot(const vector<XYDYPoint>& zmag_sqs,
                                     const string& observable_name,
                                     const string& filename, 
                                     const string& barcolor, 
                                     bool drawtoscreen)
{ 
 if (zmag_sqs.empty()) return;
 GracePlot P("Level number \\f{1}n","\\x\\cj\\C\\f{}\\h{-0.2}Z\\S(n)\\h{0.2}\\N\\x\\cj\\C\\h{-0.2}\\S2\\N");
 P.setFonts("times-roman","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.25,0.95,0.15,0.95);

 double xmin=zmag_sqs[0].xval;
 double xmax=xmin;
 for (uint ind=0;ind<zmag_sqs.size();ind++){
    if (zmag_sqs[ind].xval>xmax) xmax=zmag_sqs[ind].xval;
    if (zmag_sqs[ind].xval<xmin) xmin=zmag_sqs[ind].xval;}
 double barwidth=(xmax-xmin)/zmag_sqs.size();

 P.addBarDataSet(barcolor,"black",barwidth);
 for (uint ind=0;ind<zmag_sqs.size();ind++){
    P.addBarDataPoint(zmag_sqs[ind].xval,zmag_sqs[ind].yval);}

 P.addXYDYDataSet("none","solid","none","black");
 for (uint ind=0;ind<zmag_sqs.size();ind++){
    P.addXYDYDataPoint(zmag_sqs[ind]);}

 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.35,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen(); 
}



        // *****************************************************


string getMomentumName(int xmom, int ymom, int zmom)
{
 return make_string(xmom*xmom+ymom*ymom+zmom*zmom);
}


string getHadronName(const string& flavor, int xmom, int ymom, int zmom,
                     const string& irrep, const string& sptype, uint spid)
{
 string flav;
 if (flavor=="pion") flav="\\xp\\f{1}";
 else if (flavor=="glueball") flav="G";
 else if (flavor=="eta") flav="\\xh\\f{1}";
 else if (flavor=="phi") flav="\\xf\\f{1}";
 else if (flavor=="kaon")  flav="K";              
 else if (flavor=="kbar")  flav="\\oK\\O";            
 else if (flavor=="nucleon") flav="N";
 else if (flavor=="delta") flav="\\xD\\f{1}";
 else if (flavor=="sigma") flav="\\xS\\f{1}";
 else if (flavor=="lambda") flav="\\xL\\f{1}";
 else if (flavor=="xi") flav="\\xX\\f{1}";
 else if (flavor=="omega") flav="\\xW\\f{1}";
 else throw(std::invalid_argument("Unsupported hadron flavor"));

 string subscript(sptype+make_string(spid));
 if (irrep.length()>subscript.length())
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\s"+subscript+"\\v{}\\z{}\\M{1}\\S"
          +irrep+"\\N\\f{}";
 else
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\S"+irrep+"\\v{}\\z{}\\M{1}\\s"
          +subscript+"\\N\\f{}";
}


string getTetraquarkName(const string& flavor, int xmom, int ymom, int zmom,
                         const string& irrep, const string& sptype, 
                         uint spid, uint colortype)
{
 string flav;
 if (flavor=="isosinglet_eta_eta") flav="tquuuu1";
 else if (flavor=="isotriplet_eta_pion") flav="tquudu3";
 else if (flavor=="isosinglet_pion_pion") flav="tqdudu1";
 else if (flavor=="isotriplet_pion_pion") flav="tqdudu3";
 else if (flavor=="isoquintet_pion_pion") flav="tqdudu5";
 else if (flavor=="isodoublet_kaon_eta") flav="tqsuuu2";
 else if (flavor=="isodoublet_kaon_pion") flav="tqsudu2";
 else if (flavor=="isoquartet_kaon_pion") flav="tqsudu4";
 else if (flavor=="isotriplet_phi_pion") flav="tqssdu3";
 else if (flavor=="isosinglet_eta_phi") flav="tquuss1";
 else if (flavor=="isodoublet_kaon_phi") flav="tqsuss2";
 else if (flavor=="isosinglet_phi_phi") flav="tqssss1";
 else throw(std::invalid_argument("Unsupported tetraquark flavor"));

 if (colortype==1)
   flav+="p";
 else
   flav+="m";
 string subscript(sptype+make_string(spid));
 if (irrep.length()>subscript.length())
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\s"+subscript+"\\v{}\\z{}\\M{1}\\S"
          +irrep+"\\N\\f{}";//+make_string(colortype);
 else
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\S"+irrep+"\\v{}\\z{}\\M{1}\\s"
          +subscript+"\\N\\f{}";//+make_string(colortype);
}


string getOpStandardName(const OperatorInfo& qcd_op)
{
 if (qcd_op.isBasicLapH()){
    BasicLapHOperatorInfo qcdop(qcd_op.getBasicLapH());
    if ((qcdop.isMeson())||(qcdop.isBaryon())||(qcdop.isGlueball())){
       return getHadronName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum(),
                    qcdop.getLGIrrep(),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1));}
    else if ((qcdop.isMesonMeson())||(qcdop.isMesonBaryon())){
       string had1=getHadronName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(1),qcdop.getYMomentum(1),qcdop.getZMomentum(1),
                    qcdop.getLGIrrep(1),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1));
       string had2=getHadronName(qcdop.getFlavor(2),
                    qcdop.getXMomentum(2),qcdop.getYMomentum(2),qcdop.getZMomentum(2),
                    qcdop.getLGIrrep(2),qcdop.getSpatialType(2),
                    qcdop.getSpatialIdNumber(2));
       string iso=qcdop.getIsospin();
       if (iso=="singlet") iso="I=0";
       else if (iso=="doublet") iso="2I=1";
       else if (iso=="triplet") iso="I=1";
       else if (iso=="quartet") iso="2I=3";
       else if (iso=="quintet") iso="I=2";
       else if (iso=="sextet") iso="2I=5";
       else throw(std::invalid_argument("Unsupported total isospin in getOpStandardName"));
       return string("\\f{0}[")+had1+" "+had2+"\\f{0}]\\f{1}("
              +getMomentumName(qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum())
              +")\\S\\m{2}\\f{1}"+qcdop.getLGIrrep()+"\\N\\M{2}\\s"+iso+"\\f{}\\N";}
    else if (qcdop.isTetraquark()){
       return getTetraquarkName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum(),
                    qcdop.getLGIrrep(),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1), qcdop.getTetraquarkColorType(1));}
    else if (qcdop.isMesonMesonMeson()){
       string had1=getHadronName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(1),qcdop.getYMomentum(1),qcdop.getZMomentum(1),
                    qcdop.getLGIrrep(1),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1));
       string had2=getHadronName(qcdop.getFlavor(2),
                    qcdop.getXMomentum(2),qcdop.getYMomentum(2),qcdop.getZMomentum(2),
                    qcdop.getLGIrrep(2),qcdop.getSpatialType(2),
                    qcdop.getSpatialIdNumber(2));
       string had3=getHadronName(qcdop.getFlavor(3),
                    qcdop.getXMomentum(3),qcdop.getYMomentum(3),qcdop.getZMomentum(3),
                    qcdop.getLGIrrep(3),qcdop.getSpatialType(3),
                    qcdop.getSpatialIdNumber(3));
       string iso=qcdop.getIsospin();
       if (iso=="singlet") iso="I=0";
       else if (iso=="doublet") iso="2I=1";
       else if (iso=="triplet") iso="I=1";
       else if (iso=="quartet") iso="2I=3";
       else if (iso=="quintet") iso="I=2";
       else if (iso=="sextet") iso="2I=5";
       else if (iso=="septet") iso="I=3";
       else throw(std::invalid_argument("Unsupported total isospin in getOpStandardName"));
       return string("\\f{0}[")+had1+" "+had2+" "+had3+"\\f{0}]\\f{1}("
              +getMomentumName(qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum())
              +")\\S\\m{2}\\f{1}"+qcdop.getLGIrrep()+"\\N\\M{2}\\s"+iso+"\\f{}\\N";}}
 else if (qcd_op.isGenIrrep()){
    GenIrrepOperatorInfo qcdop(qcd_op.getGenIrrep());
    //return qcdop.getIDName()+string(" Level ")+make_string(qcdop.getIDIndex());}
    return qcdop.getIDName()+string(" ")+make_string(qcdop.getIDIndex());}
 else throw(std::invalid_argument("Unsupported operator type in getOpStandardName"));
 return string("");
}


string getMCObsStandardName(const MCObsInfo& obs)
{
 string label;
 if (obs.isRealPart()) label="Re";
 else label="Im";

 if ((obs.isCorrelatorAtTime())||(obs.isHermitianCorrelatorAtTime())){
    string tstr=string("(")+make_string(obs.getCorrelatorTimeIndex())+"), ";
    OperatorInfo src(obs.getCorrelatorSourceInfo());
    OperatorInfo snk(obs.getCorrelatorSinkInfo());
    if (src==snk) label+=" \\f{1}C\\sAA\\N"+tstr;
    else label+=" \\f{1}C\\sAB\\N"+tstr;
    label+="  \\m{3}A="+getOpStandardName(snk);
    if (src!=snk){
       label+=", ";
       string add("\\f{1}B="); add+=getOpStandardName(src);
       if (label.length()+add.length()<150) label+=add;
       else label+="\\M{3}\\V{-1.7}"+add;}
    label+="\\f{}";}

 else if (obs.isVEV()){
    OperatorInfo qcdop(obs.getVEVInfo());
    return string("\\x\\ca\\C ")+getOpStandardName(qcdop)+" \\x\\cq\\C";}

 return label;
}

string getCorrelatorStandardName(const CorrelatorInfo& corr)
{
 string label;
 string tstr=string("(t), ");
 OperatorInfo src(corr.getSource());
 OperatorInfo snk(corr.getSink());
 if (src==snk) label+=" \\f{1}C\\sAA\\N"+tstr;
 else label+=" \\f{1}C\\sAB\\N"+tstr;
 label+="  \\m{3}A="+getOpStandardName(snk);
 if (src!=snk){
    label+=", ";
    string add("\\f{1}B="); add+=getOpStandardName(src);
    if (label.length()+add.length()<150) label+=add;
    else label+="\\M{3}\\V{-1.7}"+add;}
 label+="\\f{}";
 return label;
}

    //  If "snkname" or "srcname" is "standard", use the standard name,
    //  else use the label given.

string getCorrelatorName(const CorrelatorInfo& corr, 
                         const std::string& snkname, const std::string& srcname)
{
 string label;
 string tstr=string("(t), ");
 OperatorInfo src(corr.getSource());
 OperatorInfo snk(corr.getSink());
 if (src==snk) label+=" \\f{1}C\\sAA\\N"+tstr;
 else label+=" \\f{1}C\\sAB\\N"+tstr;
 if (snkname=="standard")
    label+="  \\m{3}A="+getOpStandardName(snk);
 else
    label+="  \\m{3}A="+snkname;
 if (src!=snk){
    label+=", ";
    string add("\\f{1}B=");
    if (srcname=="standard")
       add+=getOpStandardName(src);
    else
       add+=srcname;
    if (label.length()+add.length()<150) label+=add;
    else label+="\\M{3}\\V{-1.7}"+add;}
 label+="\\f{}";
 return label;
}


// *************************************************************
