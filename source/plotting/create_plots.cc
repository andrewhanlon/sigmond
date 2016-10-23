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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
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
 if (drawtoscreen) P.drawToScreen();
}




void createEffEnergyPlotWithFit(const std::vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                const TempCorrFitInfo& fitinfo,
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
 string fitenergy("\\f{1}a\\st\\N\\E\\f{}\\sfit\\N = ");
 fitenergy+=fitres.str(2);
 P.addText(fitenergy,0.90,0.9,true,1.7,"black","top-right");

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
 if (drawtoscreen) P.drawToScreen();
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
              +")\\S\\m{2}\\f{1}"+qcdop.getLGIrrep()+"\\N\\M{2}\\s"+iso+"\\f{}\\N";}}
 else if (qcd_op.isGenIrrep()){
    GenIrrepOperatorInfo qcdop(qcd_op.getGenIrrep());
    return qcdop.getIDName()+string(" Level ")+make_string(qcdop.getIDIndex());}
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


// *************************************************************
