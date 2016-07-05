#include "task_handler.h"
#include "histogram.h"
#include "create_plots.h"
#include "task_utils.h"
#include "diag_corr_set.h"

using namespace std;

// *******************************************************************************
// *                                                                             *
// *    XML format for plotting:                                                 *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>MCValues</Type>                                                 *
// *       <MCObservable> ... </MCObservable>                                    *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <ObsName> ... </ObsName>             (optional: none default)         *
// *       <SymbolColor> ... </SymbolColor>     (optional: blue default)         *
// *       <SymbolType> ... </SymbolType>       (optional: circle default)       *
// *    </Task>                                                                  *
// *         ... <ObsName>standard</ObsName> creates name for standard ops       *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>MCBootstraps</Type>                                             *
// *       <MCObservable> ... </MCObservable>                                    *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <ObsName> ... </ObsName>             (optional: none default)         *
// *       <SymbolColor> ... </SymbolColor>     (optional: blue default)         *
// *       <SymbolType> ... </SymbolType>       (optional: circle default)       *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>MCJackknives</Type>                                             *
// *       <MCObservable> ... </MCObservable>                                    *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <ObsName> ... </ObsName>             (optional: none default)         *
// *       <SymbolColor> ... </SymbolColor>     (optional: blue default)         *
// *       <SymbolType> ... </SymbolType>       (optional: circle default)       *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>MCHistogram</Type>                                              *
// *       <NumberOfBins>25</NumberOfBins>    <!-- optional: default is 40       *
// *       <MCObservable>... </MCObservable>                                     *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <ObsName> ... </ObsName>             (optional: none default)         *
// *       <BarColor> ... </BarColor>           (optional: cyan default)         *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>MCBootstrapHistogram</Type>                                     *
// *       <NumberOfBins>25</NumberOfBins>     (optional: default is 40)         *
// *       <MCObservable>... </MCObservable>                                     *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <ObsName> ... </ObsName>             (optional: none default)         *
// *       <BarColor> ... </BarColor>           (optional: cyan default)         *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>MCJackknifeHistogram</Type>                                     *
// *       <NumberOfBins>25</NumberOfBins>     (optional: default is 40)         *
// *       <MCObservable>... </MCObservable>                                     *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <ObsName> ... </ObsName>             (optional: none default)         *
// *       <BarColor> ... </BarColor>           (optional: cyan default)         *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>TemporalCorrelator</Type>                                       *
// *       <Correlator>... </Correlator>                                         *
// *       <Arg>Re</Arg>                                                         *
// *       <HermitianMatrix/>   (optional)                                       *
// *       <SubtractVEV/>   (optional)                                           *
// *       <SamplingMode>Bootstrap</SamplingMode>  (optional: Jackknife default) *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <CorrName> ... </CorrName>           (optional: none default)         *
// *       <SymbolColor> ... </SymbolColor>     (optional: blue default)         *
// *       <SymbolType> ... </SymbolType>       (optional: circle default)       *
// *    </Task>                                                                  *
// *         ... <CorrName>standard</CorrName> creates name for standard ops     *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>EffectiveEnergy</Type>                                          *
// *       <EffEnergyType>TimeForward</EffEnergyType> (opt: TimeForward default) *
// *            or <EffEnergyType>TimeSymmetric</EffEnergyType>                  *
// *       <TimeStep>3</TimeStep>  (optional: 1 default)                         *
// *       <Correlator>... </Correlator>                                         *
// *       <Arg>Re</Arg>                                                         *
// *       <HermitianMatrix/>   (optional)                                       *
// *       <SubtractVEV/>   (optional)                                           *
// *       <SamplingMode>Bootstrap</SamplingMode>  (optional: Jackknife default) *
// *       <PlotFile> ... </PlotFile>                                            *
// *       <CorrName> ... </CorrName>           (optional: none default)         *
// *       <SymbolColor> ... </SymbolColor>     (optional: blue default)         *
// *       <SymbolType> ... </SymbolType>       (optional: circle default)       *
// *       <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                      *
// *    </Task>                                                                  *
// *                                                                             *
// *        TimeStep => solves for energy using C(t+step)/C(t)                   *
// *        EffEnergyType =>  TimeForward: C(t) = A*exp(-m*t),                   *
// *                          TimeSymmetric: C(t) = A*(exp(-m*t)+exp(-m*(T-t)))  *
// *            TimeForwardPlusConst: C(t) = A*exp(-m*t) + B0,                   *
// *            TimeSymmetricPlusConst: C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0  *
// *        EffEnergyIdName: string with no space, 26 characters or less,        *
// *                  used for ID purposes only, not needed unless want to       *
// *                  reference in later task                                    *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoPlot</Action>                                                 *
// *       <Type>EffectiveEnergies</Type>                                        *
// *       <EffEnergyType>TimeForward</EffEnergyType> (opt: TimeForward default) *
// *            or <EffEnergyType>TimeSymmetric</EffEnergyType>                  *
// *       <TimeStep>3</TimeStep>  (optional: 1 default)                         *
// *       <DiagonalCorrelatorSet>... </DiagonalCorrelatorSet>                   *
// *       <SamplingMode>Bootstrap</SamplingMode>  (optional: Jackknife default) *
// *       <PlotFileStub> ... </PlotFileStub>                                    *
// *       <SymbolColor> ... </SymbolColor>     (optional: blue default)         *
// *       <SymbolType> ... </SymbolType>       (optional: circle default)       *
// *       <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                      *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *******************************************************************************


void TaskHandler::doPlot(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 string plottype;
 xmlreadchild(xmltask,"Type",plottype,"DoPlot");

 if (plottype=="MCValues"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    string plotfile;
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    string color("blue"),obsname,symboltype("circle");
    xmlreadifchild(xmltask,"SymbolColor",color);
    xmlreadifchild(xmltask,"ObsName",obsname);
    if (obsname=="standard") obsname=getMCObsStandardName(obs);
    xmlreadifchild(xmltask,"SymbolType",symboltype);
    const Vector<double>& bins=m_obs->getBins(obs);
    double mean=m_obs->getFullSampleValue(obs);
    double stddev=m_obs->getStandardDeviation(obs);
    createMCValuesPlot(bins,obsname,mean,stddev,plotfile,symboltype,color);

    xmlout.set_root("DoPlot");
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("Type","MCValues"); 
    xmlout.put_sibling("PlotFile",plotfile); 
    XMLHandler xmlm("Mean",make_string( m_obs->getFullSampleValue(obs)));
    xmlout.put_sibling(xmlm);
    xmlm.set_root("StandardDeviation",make_string( m_obs->getStandardDeviation(obs)));
    xmlout.put_sibling(xmlm);
    for (uint jacksize=1;jacksize<=8;jacksize*=2){
       XMLHandler xmlj("JackKnifeError");
       xmlj.put_child("KnifeSize",make_string(jacksize));
       xmlj.put_child("Value",make_string(m_obs->getJackKnifeError(obs,jacksize)));
       xmlout.put_sibling(xmlj);}
    for (uint markovtime=1;markovtime<=4;markovtime++){
       XMLHandler xmla("AutoCorrelation");
       xmla.put_child("MarkovTime",make_string(markovtime));
       xmla.put_child("Value",make_string(m_obs->getAutoCorrelation(obs,markovtime)));
       xmlout.put_sibling(xmla);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with MCValues type encountered an error: ")
             +string(errmsg.what())).c_str()));}
    }

 else if (plottype=="MCBootstraps"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    string plotfile; 
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    string color("blue"),obsname,symboltype("circle");
    xmlreadifchild(xmltask,"SymbolColor",color);
    xmlreadifchild(xmltask,"ObsName",obsname);
    if (obsname=="standard") obsname=getMCObsStandardName(obs);
    xmlreadifchild(xmltask,"SymbolType",symboltype);
    Vector<double> bootvals=m_obs->getBootstrapSamplingValues(obs);
    MCEstimate res=m_obs->getBootstrapEstimate(obs);
    double mean=res.getFullEstimate();
    double upp=res.getUpperConfLimit();
    double low=res.getLowerConfLimit();
    createMCBootstrapPlot(bootvals,obsname,mean,low,upp,plotfile,symboltype,color);
    xmlout.set_root("DoPlot"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("Type","MCBootstraps"); 
    xmlout.put_sibling("PlotFile",plotfile); 
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);}  
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with MCBootstraps type encountered an error: ")
              +string(errmsg.what())).c_str()));}
    }

 else if (plottype=="MCJackknives"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    string plotfile; 
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    string color("blue"),obsname,symboltype("circle");
    xmlreadifchild(xmltask,"SymbolColor",color);
    xmlreadifchild(xmltask,"ObsName",obsname);
    if (obsname=="standard") obsname=getMCObsStandardName(obs);
    xmlreadifchild(xmltask,"SymbolType",symboltype);
    Vector<double> jackvals=m_obs->getJackknifeSamplingValues(obs);
    MCEstimate res=m_obs->getJackknifeEstimate(obs);
    double mean=res.getFullEstimate();
    double stddev=res.getSymmetricError();
    createMCJackknifePlot(jackvals,obsname,mean,stddev,plotfile,symboltype,color);
    xmlout.set_root("DoPlot"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("Type","MCJackknives"); 
    xmlout.put_sibling("PlotFile",plotfile); 
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);}  
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with MCJackknives type encountered an error: ")
            +string(errmsg.what())).c_str()));}
    }

 else if (plottype=="MCHistogram"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    string plotfile;
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    int histobins=40;
    if (xmlreadif(xmltask,"NumberOfBins",histobins,"DoPlot")){
       if (histobins<8) histobins=8;
       if (histobins>256) histobins=256;}
    string barcolor("cyan"),obsname;
    xmlreadifchild(xmltask,"BarColor",barcolor);
    xmlreadifchild(xmltask,"ObsName",obsname);
    if (obsname=="standard") obsname=getMCObsStandardName(obs);

    const Vector<double>& bins=m_obs->getBins(obs);
    Histogram mch(bins,histobins);    
    double mean=m_obs->getFullSampleValue(obs);
    double stddev=m_obs->getStandardDeviation(obs);
    createMCHistogramPlot(mch,obsname,mean,stddev,plotfile,barcolor);

    xmlout.set_root("DoPlot");
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("Type","MCHistogram"); 
    xmlout.put_sibling("PlotFile",plotfile); 
    XMLHandler xmlm("Mean",make_string( m_obs->getFullSampleValue(obs)));
    xmlout.put_sibling(xmlm);
    xmlm.set_root("StandardDeviation",make_string( m_obs->getStandardDeviation(obs)));
    xmlout.put_sibling(xmlm);
    for (uint jacksize=1;jacksize<=8;jacksize*=2){
       XMLHandler xmlj("JackKnifeError");
       xmlj.put_child("KnifeSize",make_string(jacksize));
       xmlj.put_child("Value",make_string(m_obs->getJackKnifeError(obs,jacksize)));
       xmlout.put_sibling(xmlj);}
    for (uint markovtime=1;markovtime<=4;markovtime++){
       XMLHandler xmla("AutoCorrelation");
       xmla.put_child("MarkovTime",make_string(markovtime));
       xmla.put_child("Value",make_string(m_obs->getAutoCorrelation(obs,markovtime)));
       xmlout.put_sibling(xmla);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with MCHistogram type encountered an error: ")
            +string(errmsg.what())).c_str()));}
    }


 else if (plottype=="MCBootstrapHistogram"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    string plotfile;
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    int histobins=40;
    if (xmlreadif(xmltask,"NumberOfBins",histobins,"DoPlot")){
       if (histobins<8) histobins=8;
       if (histobins>256) histobins=256;}
    string barcolor("cyan"),obsname;
    xmlreadifchild(xmltask,"BarColor",barcolor);
    xmlreadifchild(xmltask,"ObsName",obsname);
    if (obsname=="standard") obsname=getMCObsStandardName(obs);

    Vector<double> bootvals=m_obs->getBootstrapSamplingValues(obs);
    Histogram bsh(bootvals,histobins);    
    MCEstimate res=m_obs->getBootstrapEstimate(obs);
    double mean=res.getFullEstimate();
    double upp=res.getUpperConfLimit();
    double low=res.getLowerConfLimit();
    createMCBootstrapHistogramPlot(bsh,obsname,mean,low,upp,plotfile,barcolor);

    xmlout.set_root("DoPlot");
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("Type","MCBootstrapHistogram"); 
    xmlout.put_sibling("PlotFile",plotfile); 
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);}  
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with MCBootstrapHistogram type encountered an error: ")
             +string(errmsg.what())).c_str()));}
    }


 else if (plottype=="MCJackknifeHistogram"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    string plotfile;
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    int histobins=40;
    if (xmlreadif(xmltask,"NumberOfBins",histobins,"DoPlot")){
       if (histobins<8) histobins=8;
       if (histobins>256) histobins=256;}
    string barcolor("cyan"),obsname;
    xmlreadifchild(xmltask,"BarColor",barcolor);
    xmlreadifchild(xmltask,"ObsName",obsname);
    if (obsname=="standard") obsname=getMCObsStandardName(obs);

    Vector<double> jackvals=m_obs->getJackknifeSamplingValues(obs);
    Histogram bsh(jackvals,histobins);    
    MCEstimate res=m_obs->getJackknifeEstimate(obs);
    double mean=res.getFullEstimate();
    double stddev=res.getSymmetricError();
    createMCJackknifeHistogramPlot(bsh,obsname,mean,stddev,plotfile,barcolor);

    xmlout.set_root("DoPlot");
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("Type","MCJackknifeHistogram"); 
    xmlout.put_sibling("PlotFile",plotfile); 
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);}  
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with MCJackknifeHistogram type encountered an error: ")
           +string(errmsg.what())).c_str()));}
    }


 else if (plottype=="TemporalCorrelator"){
    try{
    XMLHandler xmlc(xmltask,"Correlator");
    CorrelatorInfo corr(xmlc);
    ComplexArg arg;
    read_arg_type(xmltask,arg);
    bool hermitian=(xml_tag_count(xmltask,"HermitianMatrix")==1);
    bool subvev=(xml_tag_count(xmltask,"SubtractVEV")==1);
    SamplingMode mode=Jackknife;
    string modestr;
    if (xmlreadifchild(xmltask,"SamplingMode",modestr)){
       if (modestr=="Bootstrap") mode=Bootstrap;
       else if (modestr=="Jackknife") mode=Jackknife;
       else throw(std::invalid_argument("Bad sampling mode"));} 
    map<int,MCEstimate> results;
    getCorrelatorEstimates(m_obs,corr,hermitian,subvev,arg,mode,results);
    if (results.empty()) throw(std::invalid_argument("No correlator estimates could be obtained"));

    vector<XYDYPoint> corrvals(results.size());
    uint k=0;
    for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
       corrvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                            (rt->second).getSymmetricError());}
    string plotfile;
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    string color("blue"),corrname,symboltype("circle");
    xmlreadifchild(xmltask,"SymbolColor",color);
    xmlreadifchild(xmltask,"CorrName",corrname);
    if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
    xmlreadifchild(xmltask,"SymbolType",symboltype);
    createCorrelatorPlot(corrvals,arg,corrname,plotfile,symboltype,color);

    xmlout.set_root("DoPlot");
    xmlout.put_child("Type","TemporalCorrelator"); 
    xmlout.seek_first_child();
    xmlout.put_sibling("PlotFile",plotfile);
    XMLHandler xmlt;
    corr.output(xmlt);
    xmlout.put_sibling(xmlt);
    if (arg==RealPart) xmlout.put_sibling("Arg","RealPart");
    else xmlout.put_sibling("Arg","ImaginaryPart");
    if (hermitian) xmlout.put_sibling("HermitianMatrix");
    if (subvev) xmlout.put_sibling("SubtractVEV");
    if (mode==Jackknife) xmlout.put_sibling("SamplingMode","Jackknife");
    else xmlout.put_sibling("SamplingMode","Bootstrap");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with TemporalCorrelator type encountered an error: ")
             +string(errmsg.what())).c_str()));}
    }


 else if (plottype=="EffectiveEnergy"){
    try{
    XMLHandler xmlc(xmltask,"Correlator");
    CorrelatorInfo corr(xmlc);
    ComplexArg arg;
    read_arg_type(xmltask,arg);
    bool hermitian=(xml_tag_count(xmltask,"HermitianMatrix")==1);
    bool subvev=(xml_tag_count(xmltask,"SubtractVEV")==1);
    SamplingMode mode=Jackknife;
    string instr;
    if (xmlreadifchild(xmltask,"SamplingMode",instr)){
       if (instr=="Bootstrap") mode=Bootstrap;
       else if (instr=="Jackknife") mode=Jackknife;
       else throw(std::invalid_argument("Bad sampling mode"));} 
    uint efftype=0;
    if (xmlreadifchild(xmltask,"EffEnergyType",instr)){
       if (instr=="TimeSymmetric") efftype=1;
       else if (instr=="TimeForward") efftype=0;
       else if (instr=="TimeSymmetricPlusConst") efftype=3;
       else if (instr=="TimeForwardPlusConst") efftype=2;
       else throw(std::invalid_argument("Bad effective energy type"));} 
    uint step=1;
    if (xmlreadifchild(xmltask,"TimeStep",step)){
       if ((step<1)||(step>getLatticeTimeExtent()/4))
          throw(std::invalid_argument("Bad effective energy time step"));}
    map<int,MCEstimate> results;
    getEffectiveEnergy(m_obs,corr,hermitian,subvev,arg,mode,step,efftype,results);
    if (results.empty()) throw(std::invalid_argument("No effective energy estimates could be obtained"));
    double maxerror=0.0;
    if (xmlreadifchild(xmltask,"MaxErrorToPlot",maxerror)){
       map<int,MCEstimate> raw(results);
       results.clear();
       for (map<int,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getSymmetricError()<std::abs(maxerror)) results.insert(*it);}

    vector<XYDYPoint> meffvals(results.size());
    uint k=0;
    for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
       meffvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                            (rt->second).getSymmetricError());}
    string plotfile;
    xmlreadifchild(xmltask,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()) throw(std::invalid_argument("No plot file name"));
    string color("blue"),corrname,symboltype("circle");
    xmlreadifchild(xmltask,"SymbolColor",color);
    xmlreadifchild(xmltask,"CorrName",corrname);
    if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
    xmlreadifchild(xmltask,"SymbolType",symboltype);
    createEffEnergyPlot(meffvals,arg,corrname,plotfile,symboltype,color);

    xmlout.set_root("DoPlot");
    xmlout.put_child("Type","EffectiveEnergy"); 
    xmlout.seek_first_child();
    xmlout.put_sibling("PlotFile",plotfile);
    xmlout.put_sibling("TimeStep",make_string(step));
    if (efftype==0) xmlout.put_sibling("EffEnergyType","TimeForward");
    else if (efftype==1) xmlout.put_child("EffEnergyType","TimeSymmetric");
    else if (efftype==2) xmlout.put_child("EffEnergyType","TimeForwardPlusConst");
    else if (efftype==3) xmlout.put_child("EffEnergyType","TimeSymmetricPlusConst");
    XMLHandler xmlt;
    corr.output(xmlt);
    xmlout.put_sibling(xmlt);
    if (arg==RealPart) xmlout.put_sibling("Arg","RealPart");
    else xmlout.put_sibling("Arg","ImaginaryPart");
    if (hermitian) xmlout.put_sibling("HermitianMatrix");
    if (subvev) xmlout.put_sibling("SubtractVEV");
    if (mode==Jackknife) xmlout.put_sibling("SamplingMode","Jackknife");
    else xmlout.put_sibling("SamplingMode","Bootstrap");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with EffectiveEnergy type encountered an error: ")
          +string(errmsg.what())).c_str()));}
    }



 else if (plottype=="EffectiveEnergies"){
    try{
    ArgsHandler xmlc(xmltask);
    DiagonalCorrelatorSet corrset;
    xmlc.getItem<DiagonalCorrelatorSet>("DiagonalCorrelatorSet",corrset);
    LogHelper xmllog("PlotEffectiveEnergies");
    bool subvev=corrset.isVEVsubtracted();
    string instr("Jackknife");
    xmlc.getOptionalString("SamplingMode",instr);
    SamplingMode mode=Jackknife;
    if (instr=="Bootstrap") mode=Bootstrap;
    else if (instr=="Jackknife") mode=Jackknife;
    else throw(std::invalid_argument("Bad sampling mode"));
    instr="TimeForward";
    xmlc.getOptionalString("EffEnergyType",instr);
    uint efftype=0;
    if (instr=="TimeSymmetric") efftype=1;
    else if (instr=="TimeForward") efftype=0;
    else if (instr=="TimeSymmetricPlusConst") efftype=3;
    else if (instr=="TimeForwardPlusConst") efftype=2;
    else throw(std::invalid_argument("Bad effective energy type"));
    uint step=1;
    xmlc.getOptionalUInt("TimeStep",step);
    if ((step<1)||(step>getLatticeTimeExtent()/4))
       throw(std::invalid_argument("Bad effective energy time step"));
    string plotfilestub(xmlc.getString("PlotFileStub"));
    string color("blue"),symboltype("circle");
    xmlc.getOptionalString("SymbolColor",color);
    xmlc.getOptionalString("SymbolType",symboltype);
    double maxerror=0.0;
    xmlc.getOptionalReal("MaxErrorToPlot",maxerror);
    xmllog.putEcho(xmlc);
    uint nplots=corrset.getNumberOfCorrelators();
    xmllog.putUInt("NumberOfPlots",nplots);
    bool herm=true;
    for (uint kp=0;kp<nplots;kp++){
       LogHelper xmlkp("EffEnergyPlot");
       xmlkp.putUInt("Index",kp);
       map<int,MCEstimate> results;
       getEffectiveEnergy(m_obs,corrset.getCorrelatorInfo(kp),herm,subvev,RealPart,mode,step,efftype,results);
       if (results.empty()){
          xmlkp.putString("Error","Could not make plot");
          xmllog.put(xmlkp);
          continue;}  // skip this plot
       if (maxerror>0.0){
          map<int,MCEstimate> raw(results);
          results.clear();
          for (map<int,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
             if ((it->second).getSymmetricError()<std::abs(maxerror)) results.insert(*it);}
       vector<XYDYPoint> meffvals(results.size());
       uint k=0;
       for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
          meffvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                               (rt->second).getSymmetricError());}
       string plotfile(plotfilestub+"_"+make_string(kp)+".agr");
       string corrname("Corr");
       try{corrname=getCorrelatorStandardName(corrset.getCorrelatorInfo(kp));}
       catch(const std::exception& xp){}
       createEffEnergyPlot(meffvals,RealPart,corrname,plotfile,symboltype,color);
       xmlkp.putString("PlotStatus","Success");
       xmlkp.putString("PlotFile",plotfile);
       xmllog.put(xmlkp);}
    xmllog.output(xmlout);}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoPlot with EffectiveEnergies type encountered an error: ")
          +string(errmsg.what())).c_str()));}  
    }




 else{
    throw(std::invalid_argument("DoPlot encountered unsupported type: "));}

}


// ***************************************************************************************
 
