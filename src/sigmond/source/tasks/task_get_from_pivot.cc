#include "task_handler.h"
#include "single_pivot.h"
#include "rolling_pivot.h"
#include "task_utils.h"
#include "pivoter.h"

using namespace std;

// ***********************************************************************************
// *                                                                                 *
// *                                                                                 *
// *     <Task>                                                                      *
// *       <Action>GetFromPivot</Action>                                             *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *        <EnergyName>reordered_energy</EnergyName>                                *
// *        <AmplitudeName>reordered_amplitude</AmplitudeName>                       *
// *     </Task>                                                                     *
// *                                                                                 *
// *                                                                                 *
// ***********************************************************************************



void TaskHandler::getFromPivot(XMLHandler& xml_task, XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("GetFromPivot");

 string rotatetype(xmltask.getString("Type"));
 string energy_name(xmltask.getString("EnergyName"));
 string amplitude_name(xmltask.getString("AmplitudeName"));
 if (rotatetype=="SinglePivot" || rotatetype=="RollingPivot"){
    ArgsHandler xmlpiv(xmltask,rotatetype+"Initiate");
    LogHelper xmllog;
    bool pkeep;
    Pivot pivoter;
      
    pivoter.setType(rotatetype);
    pivoter.initiatePivot(*this,xmlpiv,xmllog,pkeep);
    xmlout.putItem(xmllog);
    pivoter.checkInitiate(xmlout,xml_out);
    uint nlevels=pivoter.getNumberOfLevels();
    if (pivoter.allEnergyFitInfoAvailable()){
       XMLHandler xmles("Energies");
       for (uint level=0;level<nlevels;level++){
          MCObsInfo energyfitkey=pivoter.getEnergyKey(level);
          MCObsInfo new_energyfitkey(energy_name, level);
          doCopyBySamplings(*m_obs, energyfitkey, new_energyfitkey);
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
          XMLHandler xmlev;
          xmle.output(xmlev);
          xmles.put_child(xmlev);}
       xmlout.put(xmles);}
    else {
       xmlout.putString("Energies", "Not all energy info available");}

    if (pivoter.allAmplitudeFitInfoAvailable()){
       XMLHandler xmlas("Amplitudes");
       for (uint level=0;level<nlevels;level++){
          MCObsInfo Zrotfitkey=pivoter.getAmplitudeKey(level);
          MCObsInfo new_Zrotfitkey(amplitude_name, level);
          doCopyBySamplings(*m_obs, Zrotfitkey, new_Zrotfitkey);
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


       // delete pivoter if not put into persistent memory
    pivoter.deletePivoter(pkeep);}
    
//  if (rotatetype=="SinglePivot"){
//     ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
//     LogHelper xmllog;
//     bool pkeep;
//     SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
//                              *this,xmlpiv,xmllog,pkeep);
//     if (pivoter==0){
//        xmlout.output(xml_out);
//        throw(std::runtime_error("Could not initiate Single Pivot"));}

//     uint nlevels=pivoter->getNumberOfLevels();
//     if (pivoter->allEnergyFitInfoAvailable()){
//        XMLHandler xmles("Energies");
//        for (uint level=0;level<nlevels;level++){
//           MCObsInfo energyfitkey=pivoter->getEnergyKey(level);
//           MCObsInfo new_energyfitkey(energy_name, level);
//           doCopyBySamplings(*m_obs, energyfitkey, new_energyfitkey);
//           LogHelper xmle("EnergyLevel");
//           xmle.putUInt("Level",level);
//           XMLHandler xmlei;
//           energyfitkey.output(xmlei);
//           xmle.put(xmlei);
//           MCEstimate energy=m_obs->getEstimate(energyfitkey);
//           LogHelper result("Energy"); 
//           result.putReal("MeanValue",energy.getFullEstimate());
//           result.putReal("StandardDeviation",energy.getSymmetricError());
//           xmle.put(result);
//           XMLHandler xmlev;
//           xmle.output(xmlev);
//           xmles.put_child(xmlev);}
//        xmlout.put(xmles);}
//     else {
//        xmlout.putString("Energies", "Not all energy info available");}

//     if (pivoter->allAmplitudeFitInfoAvailable()){
//        XMLHandler xmlas("Amplitudes");
//        for (uint level=0;level<nlevels;level++){
//           MCObsInfo Zrotfitkey=pivoter->getAmplitudeKey(level);
//           MCObsInfo new_Zrotfitkey(amplitude_name, level);
//           doCopyBySamplings(*m_obs, Zrotfitkey, new_Zrotfitkey);
//           LogHelper xmla("AmplitudeLevel");
//           xmla.putUInt("Level",level);
//           XMLHandler xmlai;
//           Zrotfitkey.output(xmlai);
//           xmla.put(xmlai);
//           MCEstimate amplitude=m_obs->getEstimate(Zrotfitkey);
//           LogHelper result("Amplitude");
//           result.putReal("MeanValue",amplitude.getFullEstimate());
//           result.putReal("StandardDeviation",amplitude.getSymmetricError());
//           xmla.put(result);
//           XMLHandler xmlev;
//           xmla.output(xmlev);
//           xmlas.put_child(xmlev);}
//        xmlout.put(xmlas);}
//     else {
//        xmlout.putString("Amplitudes", "Not all amplitude info available");}

//     xmlout.putItem(xmllog);

//        // delete pivoter if not put into persistent memory
//     if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}


// ***************************************************************************************
 
