#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoAverageMomentum</Action>                                     *
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
// *     <OutFileStub> ... </OutFileStub>                                     *
// *    </Task>                                                               *
// *                                                                          *
// ****************************************************************************

void TaskHandler::doAverageMomentum(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 try{
   string outfilestub;
   xmlread(xmltask,"OutFileStub",outfilestub,"DoAverageMomentum");
   XMLHandler xmlf(xmltask,"CorrelatorsToAverage");
   list<string> tagnames;
   tagnames.push_back("CorrelatorMatrix");
   list<XMLHandler> opxml=xmlf.find_among_children(tagnames);
   set<CorrelatorMatrixInfo> m_corrinfos;
   for (list<XMLHandler>::iterator ot=opxml.begin(); ot!=opxml.end();++ot) {
     m_corrinfos.insert(CorrelatorMatrixInfo(*ot));
   }
 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
