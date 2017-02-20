#include "testing.h"
#include "matrix.h"
#include "task_utils.h"
using namespace std;


void doareorder(const std::vector<double>& energies, uint nops)
{
 cout << "NumberOfOperators = "<<nops<<endl;
 for (uint level=0;level<energies.size();level++)
    cout << "energy["<<level<<"] = "<<energies[level]<<endl;

 list<pair<double,uint> > energylevels;
 for (uint level=0;level<energies.size();level++){
    energylevels.push_back(make_pair(energies[level],level));}
 energylevels.sort(level_compare);

 RMatrix *data=new RMatrix(nops,energies.size());
 for (uint level=0;level<energies.size();level++)
 for (uint k=0;k<nops;k++)
    (*data)(k,level)=level*10000+k;

 RMatrix *buf1=new RMatrix(data->size(0),data->size(1));
 uint level=0;
 for (list<pair<double,uint> >::const_iterator it=energylevels.begin();
      it!=energylevels.end();it++,level++){
    for (uint k=0;k<buf1->size(0);k++)
       (*buf1)(k,level)=(*data)(k,it->second);}
 delete data;

 for (uint level=0;level<energies.size();level++)
 for (uint k=0;k<nops;k++)
    cout << "buf1("<<k<<","<<level<<") = "<<(*buf1)(k,level)<<endl;

 delete buf1;

 map<uint,MCObsInfo> m_energykeys;
 for (uint level=0;level<energies.size();level++)
    m_energykeys.insert(make_pair(level,MCObsInfo("Tester",level)));

 map<uint,MCObsInfo> temp(m_energykeys);
 m_energykeys.clear();
 level=0;
 for (list<pair<double,uint> >::const_iterator it=energylevels.begin();
      it!=energylevels.end();it++,level++){
    m_energykeys.insert(make_pair(level,temp.at(it->second)));}

 for (map<uint,MCObsInfo>::iterator it=m_energykeys.begin();it!=m_energykeys.end();it++)
    cout << "level = "<<it->first<<" MCObsInfo = "<<it->second.output()<<endl;

}


void testReorder(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestReorder")==0)
 return;
 xml_in.seek_unique("TestReorder");
 cout << "Beginning test of reordering"<<endl<<endl;
 vector<double> energies;
 list<XMLHandler> xmlens=xml_in.find("Energy");
 for (list<XMLHandler>::iterator it=xmlens.begin();it!=xmlens.end();it++){
    double energy;
    xmlread(*it,"Energy",energy,"TestReorder");
    energies.push_back(energy);}
 uint nops;
 xmlread(xml_in,"NumberOfOperators",nops,"TestReorder");
 doareorder(energies,nops);
}


// ***********************************************
