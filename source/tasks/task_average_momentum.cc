#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoAverageMomentum</Action>                                   *
// *     <CorrelatorsToAverage>                                               *
// *       <CorrelatorMatrix>                                                 *
// *         <BLOperator>...</BLOperator>                                     *
// *         <BLOperator>...</BLOperator>                                     *
// *              ...                                                         *
// *         <HermitianMatrix/>    (optional)                                 *
// *         <SubtractVEV/>        (optional)                                 *
// *       </CorrelatorMatrix>                                                *
// *       <CorrelatorMatrix>                                                 *
// *            ...                                                           *
// *       </CorrelatorMatrix>                                                *
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
   tagnames.push_back("CorrelatorMatrix");
   list<XMLHandler> opxml=xmlf.find_among_children(tagnames);
   vector<CorrelatorMatrixInfo> m_corrinfos;
   for (list<XMLHandler>::iterator ot=opxml.begin(); ot!=opxml.end();++ot) {
     m_corrinfos.push_back(CorrelatorMatrixInfo(*ot));
   }

   // Check that they all have the same momentum, and see if there are any missing

   for (vector<CorrelatorMatrixInfo>::iterator corrmat_it=m_corrinfos.begin();
        corrmat_it!=m_corrinfos.end(); ++corrmat_it) {
     for (corrmat_it->begin(); corrmat_it->end(); ++corrmat_it) {
       const CorrelatorInfo& corr=corrmat_it->getCurrentCorrelatorInfo();
       CorrelatorAtTimeInfo corrt(corr,timeToCompare,corrmat_it->isHermitian(),corrmat_it->isVEVSubtracted());
       MCObsInfo obskey(corrt,RealPart);
       print_obs(m_obs,obskey,xmlout);
     }
   }
 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
