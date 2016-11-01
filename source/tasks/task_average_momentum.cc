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

void print_obs(MCObsHandler* m_obs, const MCObsInfo& obskey, XMLHandler& xmlout)
{
 MCObsInfo obskeyRe(obskey); obskeyRe.setToRealPart();
 MCObsInfo obskeyIm(obskey); obskeyIm.setToImaginaryPart();
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
 
 XMLHandler xmlc; obskey.output(xmlc,false);
 xmlout.put_sibling(xmlc);
}

bool equivalent_correlators(CorrelatorInfo& corr1, CorrelatorInfo& corr2)
{
 return false;
}

void average_correlator(MCObsHandler *m_obs, CorrelatorMatrixInfo& first_corrMat, vector<CorrelatorMatrixInfo>& corrMatInfos,
                        uint timeToCompare, string filename, XMLHandler& xmlout)
{
 CorrelatorInfo corr=first_corrMat.getCurrentCorrelatorInfo();
 CorrelatorAtTimeInfo corrt(corr,timeToCompare,first_corrMat.isHermitian(),first_corrMat.isVEVSubtracted());
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

 //cout << obskeyRe.output() << endl << re_mean << "+/-" << re_err << endl << endl;
 //cout << obskeyIm.output() << endl << im_mean << "+/-" << im_err << endl << endl;

 for (vector<CorrelatorMatrixInfo>::iterator corrMat_it=corrMatInfos.begin();
      corrMat_it!=corrMatInfos.end(); ++corrMat_it) {
   for (corrMat_it->begin(); !corrMat_it->end(); ++(*corrMat_it)) {
     CorrelatorInfo corr_compare = corrMat_it->getCurrentCorrelatorInfo();
     if (equivalent_correlators(corr,corr_compare)) {
       cout << endl << endl << corr.output() << endl << "equals" << endl << corr_compare.output() << endl << endl;
     }
   }
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

   // Check that they all have the same momentum
   int mom_sqr = -1;
   /*
   for (first_corrMat.begin(); first_corrMat.end(); ++first_corrMat) {
     const CorrelatorInfo& corr=first_corrMat.getCurrentCorrelatorInfo();
     const OperatorInfo& sourceInfo = corr.getSource();
     const OperatorInfo& sinkInfo = corr.getSink();

   }
   */

   int file_count=0;
   for (first_corrMat.begin(); !first_corrMat.end(); ++first_corrMat) {
     string filename = outfilestub + "." + to_string(file_count++);
     average_correlator(m_obs,first_corrMat,corrMatInfos,timeToCompare,filename,xmlout);
   }

 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
