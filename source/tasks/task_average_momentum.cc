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

void TaskHandler::doAverageMomentum(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 try{
   string outfilestub;
   xmlread(xmltask,"OutFileStub",outfilestub,"DoAverageMomentum");
   uint timeToCompare;
   xmlread(xmltask,"TimeToCompare",timeToCompare,"DoAverageMomentum");
   XMLHandler xmlf(xmltask,"CorrelatorsToAverage");
   list<string> tagnames;
   tagnames.push_back("CorrelatorMatrix");
   list<XMLHandler> opxml=xmlf.find_among_children(tagnames);
   set<CorrelatorMatrixInfo> m_corrinfos;
   for (list<XMLHandler>::iterator ot=opxml.begin(); ot!=opxml.end();++ot) {
     m_corrinfos.insert(CorrelatorMatrixInfo(*ot));
   }

   // Check that they all have the same momentum, and see if there are any missing

   for (set<CorrelatorMatrixInfo>::iterator corrm_it=m_corrinfos.begin();
        corrm_it!=m_corrinfos.end(); ++corrm_it) {
     //for (CorrelatorMatrixInfo& corr=corrm_it->begin();
     //     corr!=corrm_it->end(); ++corr) {
     }
   }
 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
