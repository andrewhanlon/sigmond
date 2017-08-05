#include "task_handler.h"
#include "single_pivot.h"

using namespace std;

// ***********************************************************************************
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
          LogHelper xmle("EnergyLevel");
          xmle.putUInt("Level",level);
          XMLHandler xmlei;
          energyfitkey.output(xmlei);
          xmle.put(xmlei);
          MCEstimate energy=m_obs->getEstimate(energyfitkey);
          LogHelper result("Energy"); 
          result.putReal("MeanValue",energy.getFullEstimate());
          result.putReal("StandardDeviation",energy.getSymmetricError());
          xmle.put(result);
          if (refkey!=0){
             MCObsInfo enratio(string("TempEnergyRatioGwiqb"),level);
             for (m_obs->setSamplingBegin();!m_obs->isSamplingEnd();m_obs->setSamplingNext()){
                double ratiovalue=m_obs->getCurrentSamplingValue(energyfitkey)
                                 /m_obs->getCurrentSamplingValue(*refkey);
                m_obs->putCurrentSamplingValue(enratio,ratiovalue);}
             MCEstimate ratioest=m_obs->getEstimate(enratio);
             LogHelper ratio_result("EnergyRatio"); 
             ratio_result.putReal("MeanValue",ratioest.getFullEstimate());
             ratio_result.putReal("StandardDeviation",ratioest.getSymmetricError());
             xmle.put(ratio_result);}
          XMLHandler xmlev;
          xmle.output(xmlev);
          xmles.put_child(xmlev);}
       xmlout.put(xmles);}
    else {
       xmlout.putString("Energies", "Not all energy info available");}

    if (pivoter->allAmplitudeFitInfoAvailable()){
       XMLHandler xmlas("Amplitudes");
       for (uint level=0;level<nlevels;level++){
          MCObsInfo Zrotfitkey=pivoter->getAmplitudeKey(level);
          LogHelper xmla("AmplitudeLevel");
          xmla.putUInt("Level",level);
          XMLHandler xmlai;
          Zrotfitkey.output(xmlai);
          xmla.put(xmlai);
          MCEstimate amplitude=m_obs->getEstimate(Zrotfitkey);
          LogHelper result("Amplitude");
          result.putReal("MeanValue",amplitude.getFullEstimate());
          result.putReal("StandardDeviation",amplitude.getSymmetricError());
          xmla.put(result);
          XMLHandler xmlev;
          xmla.output(xmlev);
          xmlas.put_child(xmlev);}
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
 
