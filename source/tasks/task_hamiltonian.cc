#include "task_handler.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    XML format for Finite Hamiltonian:                                    *
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoHamiltonian</Action>                                       *
// *     <Energies>                                                           *
// *      <SamplingMode>Jackknife</SamplingMode>  (or Bootstrap or Current)   *
// *      <MCObservable>...</MCObservable>                                    *
// *      <MCObservable>...</MCObservable>                                    *
// *       ...                                                                *
// *     </Energies>                                                          *
// *    </Task>                                                               *
// *                                                                          *
// *                                                                          *
// ****************************************************************************

void TaskHandler::doHamiltonian(XMLHandler& xmltask, XMLHandler& xmlout, int task_count)
{
  xmlout.set_root("DoHamiltonian");

  try {
    
    XMLHandler xmle(xmltask,"Energies");
    string smode;
    xmlreadchild(xmle,"SamplingMode",smode,"TaskHandler");
    SamplingMode mode;
    if (smode=="Bootstrap") mode=Bootstrap;
    else if (smode=="Jackknife") mode=Jackknife;
    else if (smode=="Current") mode=m_obs->getCurrentSamplingMode();
    else throw(string("Invalide Sampling Mode"));
    
    list<XMLHandler> xmlh=xmle.find("MCObservable");
    set<MCObsInfo> obskeys;
    for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++) {
      obskeys.insert(MCObsInfo(*tt));
    }

    for (set<MCObsInfo>::iterator m_it=obskeys.begin(); m_it!=obskeys.end(); m_it++) {
      clog << m_it->output() << endl;

      clog << m_obs->getFullSampleValue(*m_it) << endl;
      Vector<double> Vsamps;
      if (mode==Jackknife) {
        Vsamps = m_obs->getJackknifeSamplingValues(*m_it);
      }
      else {
        Vsamps = m_obs->getBootstrapSamplingValues(*m_it);
      }

      clog << m_obs->getCurrentSamplingCovariance(*m_it,*m_it) << endl;
    }
  }
  catch(const string& errmsg) {
    xmlout.put_child("Error", string("DoHamiltonian encountered an error: ") + errmsg);
  }
}
