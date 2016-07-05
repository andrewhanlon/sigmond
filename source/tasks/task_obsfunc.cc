#include "task_handler.h"

using namespace std;

// *******************************************************************************
// *                                                                             *
// *    XML format for Functions of Observables (ratios, linear                  *
// *             superposition, ... )                                            *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>Ratio</Type>                                                    *
// *       <Result>                                                              *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </Result>                                                             *
// *       <Numerator><MCObservable> ... </MCObservable></Numerator>             *
// *       <Denominator><MCObservable> ... </MCObservable></Denominator>         *
// *       <Mode>Jackknife</Mode> (optional)                                     *
// *                      (or Bootstrap or Current [default] or Bins )           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>LinearSuperposition</Type>                                      *
// *       <Result>                                                              *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </Result>                                                             *
// *       <Summand>                                                             *
// *            <MCObservable> ... </MCObservable>                               *
// *            <Coefficient>3.2</Coefficient>   (must be real)                  *
// *       </Summand>                                                            *
// *       <Summand>                                                             *
// *            <MCObservable> ... </MCObservable>                               *
// *            <Coefficient>-5.7</Coefficient>  (must be real)                  *
// *       </Summand>                                                            *
// *           ....                                                              *
// *       <Mode>Jackknife</Mode> (optional)                                     *
// *                      (or Bootstrap or Current [default] or Bins )           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *******************************************************************************


void TaskHandler::doObsFunction(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 string functype;
 xmlread(xmltask,"Type",functype,"DoObsFunction");

 if (functype=="Ratio"){
    xmlout.set_root("DoObsFunction"); 
    xmlout.put_child("Type","Ratio");
    try{
    XMLHandler xmlnum(xmltask,"Numerator");
    XMLHandler xmlt1,xmlt2;
    MCObsInfo obsnum(xmlnum);
    xmlt1.set_root("Numerator");
    obsnum.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    XMLHandler xmlden(xmltask,"Denominator");
    MCObsInfo obsden(xmlden);
    xmlt1.set_root("Denominator");
    obsden.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    string datamode="Current";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="Bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="Current"){
       if (m_obs->isJackknifeMode()){
          mcode='J'; datamode="Jackknife";}
       else{
          mcode='B'; datamode="Bootstrap";}}
    else throw(std::invalid_argument("Invalid Sampling Mode"));
    xmlout.put_child("Mode",datamode);

    XMLHandler xmlres(xmltask,"Result");
    string name; int index;
    xmlreadchild(xmlres,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for Ratio result"));
    index=taskcount;
    xmlreadifchild(xmlres,"IDIndex",index);
    MCObsInfo resinfo(name,index,mcode=='D');
    xmlt1.set_root("ResultInfo");
    resinfo.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    if (mcode=='D'){
       const Vector<double>& numerbins=m_obs->getBins(obsnum);
       const Vector<double>& denombins=m_obs->getBins(obsden);
       int nbins=numerbins.size();
       Vector<double> ratiovalues(nbins);
       for (int bin=0;bin<nbins;bin++)
          ratiovalues[bin]=numerbins[bin]/denombins[bin];
       m_obs->putBins(resinfo,ratiovalues);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       for (m_obs->setSamplingBegin();!m_obs->isSamplingEnd();m_obs->setSamplingNext()){
          double ratiovalue=m_obs->getCurrentSamplingValue(obsnum)
                           /m_obs->getCurrentSamplingValue(obsden);
          m_obs->putCurrentSamplingValue(resinfo,ratiovalue);}
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoObsFunction with type Ratio encountered an error: ")
                +string(errmsg.what())).c_str()));}
    }

 else if (functype=="LinearSuperposition"){
    xmlout.set_root("DoObsFunction"); 
    xmlout.put_child("Type","LinearSuperposition");
    try{
    list<XMLHandler> xmlsums=xmltask.find_among_children("Summand"); 
    vector<MCObsInfo> suminfos;
    vector<double> sumcoefs;
    XMLHandler xmlt1,xmlt2;
    for (list<XMLHandler>::iterator it=xmlsums.begin();it!=xmlsums.end();it++){
       MCObsInfo obskey(*it);
       double coef=1.0;
       xmlreadifchild(*it,"Coefficient",coef);
       suminfos.push_back(obskey);
       sumcoefs.push_back(coef);
       xmlt1.set_root("Summand");
       obskey.output(xmlt2);
       xmlt1.put_child(xmlt2);
       xmlt1.put_child("Coefficient",make_string(coef));
       xmlout.put_child(xmlt1);}
    if ((suminfos.size()!=sumcoefs.size())||(suminfos.size()!=xmlsums.size()))
       throw(std::invalid_argument("Had problem reading linear superposition specs in DoObsFunction"));
    if (suminfos.size()==0){
       xmlout.put_child("Warning","No summands");
       return;}

    string datamode="Current";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="Bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="Current"){
       if (m_obs->isJackknifeMode()){
          mcode='J'; datamode="Jackknife";}
       else{
          mcode='B'; datamode="Bootstrap";}}
    else throw(std::invalid_argument("Invalid Sampling Mode"));
    xmlout.put_child("Mode",datamode);

    XMLHandler xmlres(xmltask,"Result");
    string name; int index;
    xmlreadchild(xmlres,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for LinearSuperposition result"));
    index=taskcount;
    xmlreadifchild(xmlres,"IDIndex",index);
    MCObsInfo resinfo(name,index,mcode=='D');
    xmlt1.set_root("ResultInfo");
    resinfo.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    if (mcode=='D'){
       int nsummands=suminfos.size();
       vector<const Vector<double>* > bins(nsummands);
       for (int k=0;k<nsummands;k++)
          bins[k]=&(m_obs->getBins(suminfos[k]));
       int nbins=bins[0]->size();
       Vector<double> result(nbins);
       for (int bin=0;bin<nbins;bin++){
          double temp=0.0;
          for (int k=0;k<nsummands;k++)
             temp+=sumcoefs[k]*(*(bins[k]))[bin];
          result[bin]=temp;}
       m_obs->putBins(resinfo,result);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       int nsummands=suminfos.size();
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       for (m_obs->setSamplingBegin();!m_obs->isSamplingEnd();m_obs->setSamplingNext()){
          double result=0.0;
          for (int k=0;k<nsummands;k++)
             result+=sumcoefs[k]*m_obs->getCurrentSamplingValue(suminfos[k]);
          m_obs->putCurrentSamplingValue(resinfo,result);}
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument((string("DoObsFunction with type LinearSuperposition encountered an error: ")
             +string(errmsg.what())).c_str()));} }

 else{
    throw(std::invalid_argument("DoObsFunction encountered unsupported function: "));}

}


// ***************************************************************************************
 
