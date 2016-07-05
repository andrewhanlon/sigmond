#include "task_handler.h"
#include "histogram.h"
#include "task_utils.h"

using namespace std;

// *******************************************************************************
// *                                                                             *
// *    XML format for printing:                                                 *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>MCValues</Type>                                                 *
// *       <MCObservable> ... </MCObservable>     (must be simple)               *
// *       <Verbose/> or <ShowBins/>   (optional)                                *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>MCBootstraps</Type>                                             *
// *       <MCObservable> ... </MCObservable>                                    *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>MCHistogram</Type>                                              *
// *       <NumberOfBins>25</NumberOfBins>     (optional: default is 40)         *
// *       <MCObservable>... </MCObservable>                                     *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>MCJackknives</Type>                                             *
// *       <MCObservable> ... </MCObservable>                                    *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>MCBootstrapHistogram</Type>                                     *
// *       <NumberOfBins>25</NumberOfBins>     (optional: default is 40)         *
// *       <MCObservable>... </MCObservable>                                     *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>MCJackknifeHistogram</Type>                                     *
// *       <NumberOfBins>25</NumberOfBins>     (optional: default is 40)         *
// *       <MCObservable>... </MCObservable>                                     *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>TemporalCorrelator</Type>                                       *
// *       <Correlator>... </Correlator>                                         *
// *       <Arg>Re</Arg>                                                         *
// *       <HermitianMatrix/>   (optional)                                       *
// *       <SubtractVEV/>   (optional)                                           *
// *       <SamplingMode>Bootstrap</SamplingMode>  (optional: Jackknife default) *
// *    </Task>                                                                  *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>PrintXML</Action>                                               *
// *       <Type>EffectiveEnergy</Type>                                          *
// *       <EffEnergyType>TimeForward</EffEnergyType> (opt: TimeForward default) *
// *            or <EffEnergyType>TimeSymmetric</EffEnergyType>                  *
// *       <TimeStep>3</TimeStep>  (optional: 1 default)                         *
// *       <EffEnergyIdName>PionA1um</EffEnergyIdName> (optional)                *
// *       <Correlator>... </Correlator>                                         *
// *       <Arg>Re</Arg>                                                         *
// *       <HermitianMatrix/>   (optional)                                       *
// *       <SubtractVEV/>   (optional)                                           *
// *       <SamplingMode>Bootstrap</SamplingMode>  (optional: Jackknife default) *
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
// *******************************************************************************


void TaskHandler::printXML(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 string printtype;
 xmlread(xmltask,"Type",printtype,"PrintXML");

 if (printtype=="MCValues"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    const Vector<double>& bins=m_obs->getBins(obs);
    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
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
    if ((xml_child_tag_count(xmltask,"Verbose")>0)||
        (xml_child_tag_count(xmltask,"ShowBins")>0)){
       for (unsigned int k=0;k<bins.size();k++){
          XMLHandler xmlb("Bin");
          xmlb.put_child("BinCount",make_string(k));
          xmlb.put_child("Value",make_string(bins[k]));
          xmlout.put_sibling(xmlb);}}}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with MCValues type encountered an error: ")
            +string(errmsg.what())).c_str()));}
    }

 else if (printtype=="MCBootstraps"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    Vector<double> bootvals=m_obs->getBootstrapSamplingValues(obs);
    MCEstimate res=m_obs->getBootstrapEstimate(obs);
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);
    for (unsigned int k=0;k<bootvals.size();k++){
       XMLHandler xmlb("Bootstrap");
       xmlb.put_child("ResamplingCount",make_string(k));
       xmlb.put_child("Value",make_string( bootvals[k]));
       xmlout.put_sibling(xmlb);}}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with MCBootstraps type encountered an error: ")
             +string(errmsg.what())).c_str()));}
    }

 else if (printtype=="MCJackknives"){
    try{
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    Vector<double> jackvals=m_obs->getJackknifeSamplingValues(obs);
    MCEstimate res=m_obs->getJackknifeEstimate(obs);
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);
    for (unsigned int k=0;k<jackvals.size();k++){
       XMLHandler xmlb("Jackknife");
       xmlb.put_child("ResamplingCount",make_string(k));
       xmlb.put_child("Value",make_string( jackvals[k]));
       xmlout.put_sibling(xmlb);}}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with MCJackknives type encountered an error: ")
         +string(errmsg.what())).c_str()));}
    }

 else if (printtype=="MCHistogram"){
    try{
    int histobins=40;
    if (xmlreadif(xmltask,"NumberOfBins",histobins,"printXML")){
       if (histobins<8) histobins=8;
       if (histobins>256) histobins=256;}
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    const Vector<double>& bins=m_obs->getBins(obs);
    Histogram mch(bins,histobins);    
    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    MCEstimate res=m_obs->getJackknifeEstimate(obs);
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);
    xmlout.put_sibling("Histogram");
    xmlout.put_child("NumberOfBars",make_string(mch.getNumberOfBars()));
    xmlout.seek_first_child();
    xmlout.put_sibling("BarWidth",make_string(mch.getBarWidth()));
    for (uint jbar=0;jbar<mch.getNumberOfBars();jbar++){
       XMLHandler xmlb("Bar");
       xmlb.put_child("BarCount",make_string(jbar));
       xmlb.put_child("BarRangeLowerLimit",make_string(mch.getBarRangeLowerLimit(jbar)));
       xmlb.put_child("BarRangeUpperLimit",make_string(mch.getBarRangeUpperLimit(jbar)));
       xmlb.put_child("BarHeight",make_string(mch.getBarHeight(jbar)));
       xmlout.put_sibling(xmlb);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with MCHistogram type encountered an error: ")
           +string(errmsg.what())).c_str()));}
    }


 else if (printtype=="MCBootstrapHistogram"){
    try{
    int histobins=40;
    if (xmlreadif(xmltask,"NumberOfBins",histobins,"printXML")){
       if (histobins<8) histobins=8;
       if (histobins>256) histobins=256;}
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    Vector<double> results=m_obs->getBootstrapSamplingValues(obs);
    Histogram mch(results,histobins);
    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    MCEstimate res=m_obs->getBootstrapEstimate(obs);
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);
    xmlout.put_sibling("BootstrapHistogram");
    xmlout.put_child("NumberOfBars",make_string(mch.getNumberOfBars()));
    xmlout.seek_first_child();
    xmlout.put_sibling("BarWidth",make_string(mch.getBarWidth()));
    for (uint jbar=0;jbar<mch.getNumberOfBars();jbar++){
       XMLHandler xmlb("Bar");
       xmlb.put_child("BarCount",make_string(jbar));
       xmlb.put_child("BarRangeLowerLimit",make_string(mch.getBarRangeLowerLimit(jbar)));
       xmlb.put_child("BarRangeUpperLimit",make_string(mch.getBarRangeUpperLimit(jbar)));
       xmlb.put_child("BarHeight",make_string(mch.getBarHeight(jbar)));
       xmlout.put_sibling(xmlb);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with MCBootstrapHistogram type encountered an error: ")
           +string(errmsg.what())).c_str()));}
    }

 else if (printtype=="MCJackknifeHistogram"){
    try{
    int histobins=40;
    if (xmlreadif(xmltask,"NumberOfBins",histobins,"printXML")){
       if (histobins<8) histobins=8;
       if (histobins>256) histobins=256;}
    XMLHandler xmlo(xmltask,"MCObservable");
    MCObsInfo obs(xmlo);
    Vector<double> results=m_obs->getJackknifeSamplingValues(obs);
    Histogram mch(results,histobins);
    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    obs.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    MCEstimate res=m_obs->getJackknifeEstimate(obs);
    res.output(xmlt);  
    xmlout.put_sibling(xmlt);
    xmlout.put_sibling("JackknifeHistogram");
    xmlout.put_child("NumberOfBars",make_string(mch.getNumberOfBars()));
    xmlout.seek_first_child();
    xmlout.put_sibling("BarWidth",make_string(mch.getBarWidth()));
    for (uint jbar=0;jbar<mch.getNumberOfBars();jbar++){
       XMLHandler xmlb("Bar");
       xmlb.put_child("BarCount",make_string(jbar));
       xmlb.put_child("BarRangeLowerLimit",make_string(mch.getBarRangeLowerLimit(jbar)));
       xmlb.put_child("BarRangeUpperLimit",make_string(mch.getBarRangeUpperLimit(jbar)));
       xmlb.put_child("BarHeight",make_string(mch.getBarHeight(jbar)));
       xmlout.put_sibling(xmlb);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with MCBootstrapHistogram type encountered an error: ")
          +string(errmsg.what())).c_str()));}
    }

 else if (printtype=="TemporalCorrelator"){
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

    xmlout.set_root("PrintXML"); 
    XMLHandler xmlt;
    corr.output(xmlt);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    if (arg==RealPart) xmlout.put_sibling("Arg","RealPart");
    else xmlout.put_sibling("Arg","ImaginaryPart");
    if (hermitian) xmlout.put_sibling("HermitianMatrix");
    if (subvev) xmlout.put_sibling("SubtractVEV");
    if (mode==Jackknife) xmlout.put_sibling("SamplingMode","Jackknife");
    else xmlout.put_sibling("SamplingMode","Bootstrap");
    for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++){
       XMLHandler xmlr("Estimate");
       xmlr.put_child("TimeSeparation",make_string(rt->first));
       xmlr.put_child("MeanValue",make_string((rt->second).getFullEstimate()));
       xmlr.put_child("SymmError",make_string((rt->second).getSymmetricError()));
       xmlout.put_sibling(xmlr);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with TemporalCorrelator type encountered an error: ")
           +string(errmsg.what())).c_str()));}
    }

 else if (printtype=="EffectiveEnergy"){
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

    xmlout.set_root("PrintXML"); 
    xmlout.put_child("EffectiveEnergy");
    xmlout.seek_first_child();
    xmlout.put_child("TimeStep",make_string(step));
    if (efftype==0) xmlout.put_child("EffEnergyType","TimeForward");
    else if (efftype==1) xmlout.put_child("EffEnergyType","TimeSymmetric");
    else if (efftype==2) xmlout.put_child("EffEnergyType","TimeForwardPlusConst");
    else if (efftype==3) xmlout.put_child("EffEnergyType","TimeSymmetricPlusConst");
    XMLHandler xmlt;
    corr.output(xmlt);
    xmlout.put_child(xmlt);
    if (arg==RealPart) xmlout.put_child("Arg","RealPart");
    else xmlout.put_child("Arg","ImaginaryPart");
    if (hermitian) xmlout.put_child("HermitianMatrix");
    if (subvev) xmlout.put_child("SubtractVEV");
    if (mode==Jackknife) xmlout.put_child("SamplingMode","Jackknife");
    else xmlout.put_child("SamplingMode","Bootstrap");
    xmlout.seek_root();
    xmlout.seek_first_child();
    for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++){
       XMLHandler xmlr("Estimate");
       xmlr.put_child("TimeSeparation",make_string(rt->first));
       xmlr.put_child("MeanValue",make_string((rt->second).getFullEstimate()));
       xmlr.put_child("SymmError",make_string((rt->second).getSymmetricError()));
       xmlout.put_sibling(xmlr);}
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("PrintXML with EffectiveEnergy type encountered an error: ")
              +string(errmsg.what())).c_str()));}
    }

 else{
    throw(std::invalid_argument("PrintXML encountered unsupported type: "));}


}


// ***************************************************************************************
 
