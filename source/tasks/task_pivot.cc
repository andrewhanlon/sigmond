#include "task_handler.h"
#include "single_pivot.h"

using namespace std;

// ***********************************************************************************
// *                                                                                 *
// *   Two tasks are defined in this file:                                           *
// *                                                                                 *
// *      (a)  <Action>InsertIntoPivot</Action>                                      *
// *                                                                                 *
// *      (b)  <Action>GetFromPivot</Action>                                         *
// *                                                                                 *
// *     <Task>                                                                      *
// *       <Action>InsertIntoPivot</Action>                                          *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *                (short form below)                                               *
// *        <EnergyFitCommonName>En</EnergyFitCommonName>                            *
// *                (or specify individually)                                        *
// *        <EnergyFit>                                                              *
// *           <Level>0</Level>                                                      *
// *           <Name>A</Name><IDIndex>0</IDIndex> (name of fit energy observable)    *
// *        </EnergyFit>                                                             *  
// *            ... one for each level                                               *
// *                                                                                 *
// *                (short form below)                                               *
// *        <RotatedAmplitudeCommonName>Amp</RotatedAmplitudeCommonName>             *
// *                (or specify individually)                                        *
// *        <RotatedAmplitude>                                                       *
// *           <Level>0</Level>                                                      *
// *           <Name>A</Name><IDIndex>0</IDIndex>                                    *
// *        </RotatedAmplitude>                                                      *
// *            ... one for each level                                               *
// *     </Task>                                                                     *
// *                                                                                 *
// *                                                                                 *
// *     <Task>                                                                      *
// *       <Action>GetFromPivot</Action>                                             *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *        <ReferenceEnergy>  (optional: gives energies as a ratio over the ref.)   *
// *          <Name>kaon</Name><IDIndex>0<IDIndex>                                   *
// *        </ReferenceEnergy>                                                       *
// *     </Task>                                                                     *
// *                                                                                 *
// *                                                                                 *
// ***********************************************************************************



void TaskHandler::insertIntoPivot(XMLHandler& xml_task, XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("InsertIntoPivot");

 map<uint,MCObsInfo> energyfits;
 string encommon;
 xmltask.getOptionalString("EnergyFitCommonName",encommon);
 if (!encommon.empty()){
    xmlout.putString("EnergyFitCommonName",encommon);}
 else{
    list<XMLHandler> xmlens=xml_task.find_among_children("EnergyFit"); 
    for (list<XMLHandler>::iterator it=xmlens.begin();it!=xmlens.end();it++){
       ArgsHandler xmle(*it);
       uint level=xmle.getUInt("Level");
       string name(xmle.getString("Name"));
       uint index=taskcount;
       xmle.getOptionalUInt("IDIndex",index);
       MCObsInfo energykey(name,index);
       energyfits.insert(make_pair(level,energykey));
       xmlout.putEcho(xmle);}}

 map<uint,MCObsInfo> ampfits;
 string ampcommon;
 xmltask.getOptionalString("RotatedAmplitudeCommonName",ampcommon);
 if (!ampcommon.empty()){
    xmlout.putString("RotatedAmplitudeCommonName",ampcommon);}
 else{
    list<XMLHandler> xmlamps=xml_task.find_among_children("RotatedAmplitude"); 
    for (list<XMLHandler>::iterator it=xmlamps.begin();it!=xmlamps.end();it++){
       ArgsHandler xmla(*it);
       uint level=xmla.getUInt("Level");
       string name(xmla.getString("Name"));
       uint index=taskcount;
       xmla.getOptionalUInt("IDIndex",index);
       MCObsInfo ampkey(name,index);
       ampfits.insert(make_pair(level,ampkey));
       xmlout.putEcho(xmla);}}

 string rotatetype(xmltask.getString("Type"));
 
 if (rotatetype=="SinglePivot"){
    ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
    LogHelper xmllog;
    bool pkeep;
    SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
                             *this,xmlpiv,xmllog,pkeep);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}

    if (!encommon.empty()){
       MCObsInfo encommonkey(encommon,0);
       for (uint level=0;level<pivoter->getNumberOfLevels();level++){
          encommonkey.resetObsIndex(level);
          pivoter->insertEnergyFitInfo(level,encommonkey);}}
    else{
       for (map<uint,MCObsInfo>::iterator it=energyfits.begin();it!=energyfits.end();it++)
          pivoter->insertEnergyFitInfo(it->first,it->second);}

    
    if (!ampcommon.empty()){
       MCObsInfo ampcommonkey(ampcommon,0);
       for (uint level=0;level<pivoter->getNumberOfLevels();level++){
          ampcommonkey.resetObsIndex(level);
          pivoter->insertAmplitudeFitInfo(level,ampcommonkey);}}
    else{
       for (map<uint,MCObsInfo>::iterator it=ampfits.begin();it!=ampfits.end();it++)
          pivoter->insertAmplitudeFitInfo(it->first,it->second);}

    xmlout.putItem(xmllog);

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}



void TaskHandler::getFromPivot(XMLHandler& xml_task, XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("GetFromPivot");

 uint refcount=xml_task.count("ReferenceEnergy");
 MCObsInfo* refkey=0;
 if (refcount==1){
    XMLHandler xmlref(xml_task,"ReferenceEnergy");
    string refname; int refindex;
    xmlreadchild(xmlref,"Name",refname);
    if (refname.empty()) throw(std::invalid_argument("Must provide name for reference energy"));
    refindex=taskcount;
    xmlreadifchild(xmlref,"IDIndex",refindex);
    refkey = new MCObsInfo(refname,refindex);}
 
 string rotatetype(xmltask.getString("Type"));
 
 if (rotatetype=="SinglePivot"){
    ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
    LogHelper xmllog;
    bool pkeep;
    SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
                             *this,xmlpiv,xmllog,pkeep);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}

    uint nlevels=pivoter->getNumberOfLevels();
    if (pivoter->allEnergyFitInfoAvailable()){
       XMLHandler xmles("Energies");
       if (refkey!=0){
          XMLHandler xmlre("ReferenceEnergy");
          XMLHandler xmlrei;
          refkey->output(xmlrei);
          xmlre.put_child(xmlrei);
          MCEstimate refenergy=m_obs->getEstimate(*refkey);
          XMLHandler xmlree;
          refenergy.output(xmlree);
          xmlre.put_child(xmlree);
          xmles.put_child(xmlre);
       }
       for (uint level=0;level<nlevels;level++){
          MCObsInfo energyfitkey=pivoter->getEnergyKey(level);
          XMLHandler xmle("Energy");
          xmle.put_child("Level",make_string(level));
          XMLHandler xmlei;
          energyfitkey.output(xmlei);
          xmle.put_child(xmlei);
          if (refkey==0){
             MCEstimate energy=m_obs->getEstimate(energyfitkey);
             XMLHandler xmlee;
             energy.output(xmlee);
             xmle.put_child(xmlee);}
          else{
             MCObsInfo enratio(string("TempEnergyRatioGwiqb"),level);
             for (m_obs->setSamplingBegin();!m_obs->isSamplingEnd();m_obs->setSamplingNext()){
                double ratiovalue=m_obs->getCurrentSamplingValue(energyfitkey)
                                 /m_obs->getCurrentSamplingValue(*refkey);
                m_obs->putCurrentSamplingValue(enratio,ratiovalue);}
             MCEstimate ratioest=m_obs->getEstimate(enratio);
             XMLHandler xmlee;
             ratioest.output(xmlee);
             xmle.put_child(xmlee);}
          xmles.put_child(xmle);}
       xmlout.put(xmles);}
    else {
       xmlout.putString("Energies", "Not all energy info available");}

    if (pivoter->allAmplitudeFitInfoAvailable()){
       XMLHandler xmlas("Amplitudes");
       for (uint level=0;level<nlevels;level++){
          MCObsInfo Zrotfitkey=pivoter->getAmplitudeKey(level);
          XMLHandler xmla("Amplitude");
          xmla.put_child("Level",make_string(level));
          XMLHandler xmlai;
          Zrotfitkey.output(xmlai);
          xmla.put_child(xmlai);
          MCEstimate amplitude=m_obs->getEstimate(Zrotfitkey);
          XMLHandler xmlae;
          amplitude.output(xmlae);
          xmla.put_child(xmlae);
          xmlas.put_child(xmla);}
       xmlout.put(xmlas);}
    else {
       xmlout.putString("Amplitudes", "Not all amplitude info available");}

    xmlout.putItem(xmllog);

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 if (refkey!=0) delete refkey;
 xmlout.output(xml_out);
}


// ***************************************************************************************
 
