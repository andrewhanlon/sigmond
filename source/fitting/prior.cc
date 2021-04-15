#include "prior.h"
using namespace std;

Prior::Prior(XMLHandler& xmlin, MCObsHandler& OH)    : m_obs(&OH)
{
 m_resampled = (xmlin.count_among_children("Name")>0);
 if (m_resampled) {
    string obs_name; 
    uint obs_id;
    xmlreadchild(xmlin,"Name",obs_name);
    xmlreadchild(xmlin,"IDIndex",obs_id);
    m_prior = MCObsInfo(obs_name,obs_id);}
 else {
    xmlreadchild(xmlin,"Mean",m_mean);
    xmlreadchild(xmlin,"Error",m_error);
 }
}
