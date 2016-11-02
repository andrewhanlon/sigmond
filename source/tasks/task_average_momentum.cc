#include <string>
#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoAverageMomentum</Action>                                   *
// *     <CorrelatorsToAverage>                                               *
// *       <CorrelatorMatrixInfo>                                             *
// *         <BLOperator>...</BLOperator>                                     *
// *         <BLOperator>...</BLOperator>                                     *
// *              ...                                                         *
// *         <HermitianMatrix/>    (optional)                                 *
// *         <SubtractVEV/>        (optional)                                 *
// *       </CorrelatorMatrixInfo>                                            *
// *       <CorrelatorMatrixInfo>                                             *
// *            ...                                                           *
// *       </CorrelatorMatrixInfo>                                            *
// *          ...                                                             *
// *     </CorrelatorsToAverage>                                              *
// *     <TimeToCompare>3</TimeToCompare>                                     *
// *     <OutFileStub> ... </OutFileStub>                                     *
// *    </Task>                                                               *
// *                                                                          *
// ****************************************************************************

void average_correlator(MCObsHandler *m_obs, CorrelatorAtTimeInfo& corrt, vector<CorrelatorAtTimeInfo>& toAverage,
                        uint timeToCompare, string filename, XMLHandler& xmlout)
{
 double coef = 1./(1.+toAverage.size());
 // Get comparison values
 MCObsInfo obskeyRe(corrt,RealPart);
 MCObsInfo obskeyIm(corrt,ImaginaryPart);
 if ((!m_obs->queryBins(obskeyRe))&&(!m_obs->queryBins(obskeyIm))) return;
 double re_mean=0.0; double re_err=1.0;
 double im_mean=0.0; double im_err=1.0;
 if (m_obs->queryBins(obskeyRe)){
    MCEstimate est=m_obs->getJackknifeEstimate(obskeyRe);
    re_mean=std::abs(est.getFullEstimate());
    re_err=est.getSymmetricError();}
 if (m_obs->queryBins(obskeyIm)){
    MCEstimate est=m_obs->getJackknifeEstimate(obskeyIm);
    im_mean=std::abs(est.getFullEstimate());
    im_err=est.getSymmetricError();}

 cout << "first (Re): " << re_mean << "+/-" << re_err << endl;
 cout << "first (Im): " << im_mean << "+/-" << im_err << endl;

 for (vector<CorrelatorAtTimeInfo>::iterator corr_it=toAverage.begin();
      corr_it!=toAverage.end(); ++corr_it) {
    MCObsInfo obskeyCompRe(*corr_it,RealPart);
    MCObsInfo obskeyCompIm(*corr_it,ImaginaryPart);
    if ((!m_obs->queryBins(obskeyCompRe))&&(!m_obs->queryBins(obskeyCompIm))) return;
    if (m_obs->queryBins(obskeyCompRe)){
      MCEstimate est=m_obs->getJackknifeEstimate(obskeyCompRe);
      re_mean=std::abs(est.getFullEstimate());
      re_err=est.getSymmetricError();}
    if (m_obs->queryBins(obskeyCompIm)){
      MCEstimate est=m_obs->getJackknifeEstimate(obskeyCompIm);
      im_mean=std::abs(est.getFullEstimate());
      im_err=est.getSymmetricError();}

    cout << "Compare (Re): " << re_mean << "+/-" << re_err << endl;
    cout << "Compare (Im): " << im_mean << "+/-" << im_err << endl;
 }

}


void TaskHandler::doAverageMomentum(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 try{
   xmlout.set_root("DoAverageMomentum");
   string outfilestub;
   xmlread(xmltask,"OutFileStub",outfilestub,"DoAverageMomentum");
   uint timeToCompare;
   xmlread(xmltask,"TimeToCompare",timeToCompare,"DoAverageMomentum");
   XMLHandler xmlf(xmltask,"CorrelatorsToAverage");
   list<string> tagnames;
   tagnames.push_back("CorrelatorMatrixInfo");
   list<XMLHandler> corrMatxml=xmlf.find_among_children(tagnames);
   vector<CorrelatorMatrixInfo> corrMatInfos;
   list<XMLHandler>::iterator ct=corrMatxml.begin();
   CorrelatorMatrixInfo first_corrMat(*ct);
   for (++ct; ct!=corrMatxml.end();++ct) {
     corrMatInfos.push_back(CorrelatorMatrixInfo(*ct));
   }

   int file_count=0;
   for (first_corrMat.begin(); !first_corrMat.end(); ++first_corrMat) {
     CorrelatorInfo corr=first_corrMat.getCurrentCorrelatorInfo();
     vector<CorrelatorAtTimeInfo> toAverage;
     for (vector<CorrelatorMatrixInfo>::iterator corrMat_it=corrMatInfos.begin();
          corrMat_it!=corrMatInfos.end(); ++corrMat_it) {
       int n = 0;
       for (corrMat_it->begin(); !corrMat_it->end(); ++(*corrMat_it)) {
         CorrelatorInfo corr_compare = corrMat_it->getCurrentCorrelatorInfo();
         if (corr.rotationallyEquivalent(corr_compare)) {
           n++;
           CorrelatorAtTimeInfo corrt_compare(corr_compare,timeToCompare,corrMat_it->isHermitian(),corrMat_it->isVEVSubtracted());
           toAverage.push_back(corrt_compare);
         }
       }
       if (n!=1) cout << "warning" << endl;
     }
  
     string filename = outfilestub + "." + to_string(file_count++);
     CorrelatorAtTimeInfo corrt(corr,timeToCompare,first_corrMat.isHermitian(),first_corrMat.isVEVSubtracted());
     average_correlator(m_obs,corrt,toAverage,timeToCompare,filename,xmlout);
   }

 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
