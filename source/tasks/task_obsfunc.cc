#include "task_handler.h"
#include "task_utils.h"

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
// *      For a correlator matrix, compute C(t)-C(t+1)                           *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>CorrelatorMatrixTimeDifferences</Type>                          *
// *       <NewOperatorOrderedList>                                              *
// *         <Operator>...</Operator>                                            *
// *            ...                                                              *
// *       </NewOperatorOrderedList>                                             *
// *       <OriginalOperatorOrderedList>                                         *
// *         <Operator>...</Operator>                                            *
// *            ...                                                              *
// *       </OriginalOperatorOrderedList>                                        *
// *       <MinimumTimeSeparation>3</MinimumTimeSeparation>                      *
// *       <MaximumTimeSeparation>12</MaximumTimeSeparation>                     *
// *       <HermitianMatrix/>  (if hermitian)                                    *
// *       <WriteToBinFile>filename</WriteToBinFile>                             *
// *       <FileMode>overwrite</FileMode>   (optional)                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      Combine several correlation matrices to make a resultant matrix:       *
// *          R[i,j] = sum_k  d[k][i]*d[k][j] * C[k][i,j]                        *
// *      where k labels the different matrices, i and j are the row and column, *
// *      and the superposition coefficients are specified in d[k][i].           *
// *      The linear superposition coefficients are given at the operator        *
// *      level.  Each correlator matrix to combine must be given in an          *
// *      <OperatorOrderedList> tag, which includes <Item>... tags.  Each        *
// *      item in each list must have an <Operator> tag, and an optional         *
// *      <Coefficient> tag; if absent, a coefficient equal to 1 is assumed.     *
// *      All coefficients must be REAL.                                         *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>CorrelatorMatrixSuperposition</Type>                            *
// *       <ResultOperatorOrderedList>                                           *
// *         <Operator>...</Operator>                                            *
// *            ...                                                              *
// *       </ResultOperatorOrderedList>                                          *
// *       <OperatorOrderedList>                                                 *
// *         <Item><Operator>...</Operator><Coeffient>-1.0</Coefficient></Item>  *
// *            ...                                                              *
// *       </OperatorOrderedList>                                                *
// *       <OperatorOrderedList>                                                 *
// *             ...                                                             *
// *       </OperatorOrderedList>                                                *
// *             ...                                                             *
// *       <MinimumTimeSeparation>3</MinimumTimeSeparation>                      *
// *       <MaximumTimeSeparation>12</MaximumTimeSeparation>                     *
// *       <HermitianMatrix/>  (if hermitian)                                    *
// *       <WriteToBinFile>filename</WriteToBinFile>                             *
// *       <FileMode>overwrite</FileMode>   (optional)                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      For boosting an energy specified in the <FrameEnergy> tag              *
// *      to a frame new frame using the momentum specified by the               *
// *      <IntMomSquared> tag. If the <BoostToCM/> tag is present                *
// *      it is assumed that the energy is in a lab frame specified by the       *
// *      <IntMomSquared> tag and the resulting energy is boosted to the         *
// *      center of momentum frame.                                              *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>BoostEnergy</Type>                                              *
// *       <BoostToCM/>     (optional)                                           *
// *       <Result>                                                              *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </Result>                                                             *
// *       <IntMomSquared>4</IntMomSquared>                                      *
// *       <SpatialExtentNumSites>32</SpatialExtentNumSites>                     *
// *       <FrameEnergy><MCObservable> ... </MCObservable></FrameEnergy>         *
// *       <Anisotropy><MCObservable> ...           (optional)                   *
// *                           </MCObservable></Anisotropy>                      *
// *       <Mode>Jackknife</Mode> (optional)                                     *
// *                      (or Bootstrap or Current [default] or Bins )           *
// *       <ReferenceEnergy>   (optional)                                        *
// *          <MCObservable> ... </MCObservable>                                 *
// *       </ReferenceEnergy>                                                    *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      For extracting the `interaction energy' or energy difference           *
// *      between an interacting energy level and the nearest                    *
// *      non-interacting state, the ratio:                                      *
// *                               C_int(t)                                      *
// *               R(t) =  -------------------------                             *
// *                         prod_i C_non-int[i](t)                              *
// *      can be fit to the ansatz:  R(t) = A exp(-DeltaE*t), where:             *
// *        DeltaE          =  E_int - E_non-int,                                *
// *        C_int(t)        =  (diagonal) correlator for interacting level       *
// *        C_non-int[i](t) =  N correlators coresponding to closest             *
// *                           non-interacting level                             *
// *                                                                             *
// *      Example:                                                               *
// *      if the closest non-interacting level to C_int is pi(0)K(1)K(1):        *
// *        C_non-int[0](t) = pion correlator for P^2=0                          *
// *        C_non-int[1](t) = kaon correlator for P^2=1                          *
// *        C_non-int[2](t) = kaon correlator for P^2=1                          *
// *                                                                             *
// *      Note: the resultant correlator will have (if specified) all            *
// *            VEVs subtracted in this task.                                    *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>CorrelatorInteractionRatio</Type>                               *
// *       <Result>                                                              *
// *          <Operator>...</Operator>                                           *
// *       </Result>                                                             *
// *       <InteractingOperator>                                                 *
// *          <Operator>...</Operator>                                           *
// *          <SubtractVEV />    (optional)                                      *
// *       </InteractingOperator>                                                *
// *       <NonInteractingOperator>                                              *
// *          <Operator>...</Operator>                                           *
// *          <SubtractVEV />    (optional)                                      *
// *       </NonInteractingOperator>                                             *
// *          ......                                                             *
// *       <NonInteractingOperator>                                              *
// *          <Operator>...</Operator>                                           *
// *          <SubtractVEV />    (optional)                                      *
// *       </NonInteractingOperator>                                             *
// *       <MinimumTimeSeparation>0</MinimumTimeSeparation>                      *
// *       <MaximumTimeSeparation>15</MaximumTimeSeparation>                     *
// *       <WriteToBinFile>filename</WriteToBinFile>   (optional)                *
// *       <FileMode>overwrite</FileMode>   (optional)                           *
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
       doRatioByBins(*m_obs,obsnum,obsden,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       doRatioBySamplings(*m_obs,obsnum,obsden,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type Ratio encountered an error: ")
                +string(errmsg.what())));}
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
       doLinearSuperpositionByBins(*m_obs,suminfos,sumcoefs,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       doLinearSuperpositionBySamplings(*m_obs,suminfos,sumcoefs,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type LinearSuperposition encountered an error: ")
             +string(errmsg.what())));} }

 else if (functype=="BoostEnergy"){
   xmlout.set_root("DoObsFunction");
   xmlout.put_child("Type","BoostEnergy");
   try{
     XMLHandler xmlrest(xmltask,"FrameEnergy");
     XMLHandler xmlt1,xmlt2;
     MCObsInfo obsframe(xmlrest);
     xmlt1.set_root("FrameEnergy");
     obsframe.output(xmlt2);
     xmlt1.put_child(xmlt2);
     xmlout.put_child(xmlt1);

     uint refcount=xmltask.count("ReferenceEnergy");
     MCObsInfo* refkey=0;
     if (refcount==1){
       XMLHandler xmlref(xmltask,"ReferenceEnergy");
       string refname; int refindex;
       xmlreadchild(xmlref,"Name",refname);
       if (refname.empty()) throw(std::invalid_argument("Must provide name for reference energy"));
       refindex=taskcount;
       xmlreadifchild(xmlref,"IDIndex",refindex);
       refkey = new MCObsInfo(refname,refindex);
       xmlt1.set_root("ReferenceEnergy");
       refkey->output(xmlt2);
       xmlt1.put_child(xmlt2);
       MCEstimate refenergy=m_obs->getEstimate(*refkey);
       XMLHandler xmlre;
       refenergy.output(xmlre);
       xmlt1.put_child(xmlre);
       xmlout.put_child(xmlt1);}

     uint aniscount=xmltask.count("Anisotropy");
     MCObsInfo* obsxi=0;
     if (aniscount==1){
       XMLHandler xmlxi(xmltask,"Anisotropy");
       obsxi = new MCObsInfo(xmlxi);
       xmlt1.set_root("Anisotropy");
       obsxi->output(xmlt2);
       xmlt1.put_child(xmlt2);
       xmlout.put_child(xmlt1);}

     int psq=-1;
     xmlreadifchild(xmltask,"IntMomSquared",psq);
     if (psq<0) throw(std::invalid_argument("Must provide positive Integer Momentum Squared for Boost"));

     uint m_lat_spatial_extent;
     xmlreadifchild(xmltask,"SpatialExtentNumSites",m_lat_spatial_extent);
     if (m_lat_spatial_extent<1)
       throw(std::invalid_argument("Lattice spatial extent must be a positive integer"));
     double m_momsq_quantum=6.2831853071795864770/double(m_lat_spatial_extent);
     m_momsq_quantum*=m_momsq_quantum;
     double psqfactor=psq*m_momsq_quantum;

     string datamode="Current";
     xmlreadifchild(xmltask,"Mode",datamode);
     char mcode;
     if (datamode=="Bootstrap") mcode='B';
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
     if (name.empty()) throw(std::invalid_argument("Must provide name for Boost result"));
     index=taskcount;
     xmlreadifchild(xmlres,"IDIndex",index);
     MCObsInfo resinfo(name,index,mcode=='D');
     xmlt1.set_root("ResultInfo");
     resinfo.output(xmlt2);
     xmlt1.put_child(xmlt2);
     xmlout.put_child(xmlt1);

     SamplingMode origmode=m_obs->getCurrentSamplingMode();
     if (mcode=='J') m_obs->setToJackknifeMode();
     else m_obs->setToBootstrapMode();

     ArgsHandler xmlcm(xmltask);
     bool boostcm=false;
     xmlcm.getOptionalBool("BoostToCM", boostcm);
     if (boostcm) psqfactor = -psqfactor;

     if (aniscount==1)
       doBoostBySamplings(*m_obs,obsframe,*obsxi,psqfactor,resinfo);
     else
       doBoostBySamplings(*m_obs,obsframe,psqfactor,resinfo);

     if (refcount==1)
       doRatioBySamplings(*m_obs,resinfo,*refkey,resinfo);

     if (obsxi!=0) delete obsxi;
     if (refkey!=0) delete refkey;

     MCEstimate est=m_obs->getEstimate(resinfo);
     est.output(xmlt1);
     xmlout.put_child(xmlt1);
     m_obs->setSamplingMode(origmode);}
   catch(const std::exception& errmsg){
     xmlout.clear();
     throw(std::invalid_argument(string("DoObsFunction with type BoostEnergy encountered an error: ")
				 +string(errmsg.what())));}
 }


 else if (functype=="CorrelatorMatrixTimeDifferences"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","CorrelatorMatrixTimeDifferences");
    try{
    list<string> tagnames;
    tagnames.push_back("Operator");
    tagnames.push_back("OperatorString");
    tagnames.push_back("BLOperator");
    tagnames.push_back("BLOperatorString");
    tagnames.push_back("GIOperator");
    tagnames.push_back("GIOperatorString");
    XMLHandler xmlnew(xmltask,"NewOperatorOrderedList");
    list<XMLHandler> newopxml=xmlnew.find_among_children(tagnames);
    XMLHandler xmlorig(xmltask,"OriginalOperatorOrderedList");
    list<XMLHandler> origopxml=xmlorig.find_among_children(tagnames);
    list<OperatorInfo> newops,origops;
    for (list<XMLHandler>::iterator ot=newopxml.begin();ot!=newopxml.end();++ot)
       newops.push_back(OperatorInfo(*ot));
    for (list<XMLHandler>::iterator ot=origopxml.begin();ot!=origopxml.end();++ot)
       origops.push_back(OperatorInfo(*ot));
    if (newops.size()!=origops.size())
       throw(std::runtime_error("Mismatch in number of original and new operators in CorrelatorMatrixTimeDifferences"));
    uint tmin,tmax;
    xmlreadchild(xmltask,"MinimumTimeSeparation",tmin);
    xmlreadchild(xmltask,"MaximumTimeSeparation",tmax);
    bool herm=(xmltask.count("HermitianMatrix")>0) ? true: false;
    string filename;
    bool writetofile=xmlreadifchild(xmltask,"WriteToBinFile",filename);
    if (filename.empty()) writetofile=false;
    bool overwrite = false;  // protect mode
    if ((writetofile)&&(xml_tag_count(xmltask,"FileMode")==1)){
       string fmode;
       xmlread(xmltask,"FileMode",fmode,"FileListInfo");
       fmode=tidyString(fmode);
       if (fmode=="overwrite") overwrite=true;}

    xmlout.put_child("MinimumTimeSeparation",make_string(tmin));
    xmlout.put_child("MaximumTimeSeparation",make_string(tmax));
    if (herm) xmlout.put_child("HermitianMatrix");
    XMLHandler xmlo("OriginalOperatorOrderedList");
    for (list<OperatorInfo>::const_iterator it=origops.begin();it!=origops.end();it++){
       XMLHandler xmloo; it->output(xmloo); xmlo.put_child(xmloo);}
    xmlout.put_child(xmlo);
    xmlo.set_root("NewOperatorOrderedList");
    for (list<OperatorInfo>::const_iterator it=newops.begin();it!=newops.end();it++){
       XMLHandler xmloo; it->output(xmloo); xmlo.put_child(xmloo);}
    xmlout.put_child(xmlo);
    list<OperatorInfo>::const_iterator oldrow,oldcol,newrow,newcol;
    newcol=newops.begin();
    uint count=0;
    set<MCObsInfo> obskeys;
    for (oldcol=origops.begin();oldcol!=origops.end();oldcol++,newcol++){
       newrow=(herm?newcol:newops.begin());
       for (oldrow=(herm?oldcol:origops.begin());oldrow!=origops.end();oldrow++,newrow++){
          CorrelatorAtTimeInfo origcorr(*oldrow,*oldcol,0,herm,false);
          CorrelatorAtTimeInfo newcorr(*newrow,*newcol,0,herm,false);
          for (uint t=tmin;t<tmax;t++){
             origcorr.resetTimeSeparation(t);
             const RVector& bins1=m_obs->getBins(MCObsInfo(origcorr));
             origcorr.resetTimeSeparation(t+1);
             const RVector& bins2=m_obs->getBins(MCObsInfo(origcorr));
             RVector newbins(bins1);
             newbins-=bins2;
             newcorr.resetTimeSeparation(t);
             MCObsInfo newkey(newcorr);
             m_obs->putBins(newkey,newbins);
             if (writetofile) obskeys.insert(newkey);
             origcorr.resetTimeSeparation(t);
             m_obs->eraseData(MCObsInfo(origcorr));
             count++;
#ifdef COMPLEXNUMBERS
             origcorr.resetTimeSeparation(t);
             const RVector& ibins1=m_obs->getBins(MCObsInfo(origcorr,ImaginaryPart));
             origcorr.resetTimeSeparation(t+1);
             const RVector& ibins2=m_obs->getBins(MCObsInfo(origcorr,ImaginaryPart));
             newbins=ibins1;
             newbins-=ibins2;
             newkey.setToImaginaryPart();
             m_obs->putBins(newkey,newbins);
             if (writetofile) obskeys.insert(newkey);
             origcorr.resetTimeSeparation(t);
             m_obs->eraseData(MCObsInfo(origcorr,ImaginaryPart));
             count++;
#endif
             }}}
       xmlout.put_child("NumberOfRealObservablesProcessed",make_string(count));
       if (writetofile){
          XMLHandler xmlf;
          m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
          xmlout.put_child("WriteBinsToFile",filename);}
       xmlout.put_child("Status","Done");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type CorrelatorMatrixTimeDifferences encountered an error: ")
             +string(errmsg.what())));} }


 else if (functype=="CorrelatorMatrixSuperposition"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","CorrelatorMatrixSuperposition");
    try{
    list<string> tagnames;
    tagnames.push_back("Operator");
    tagnames.push_back("OperatorString");
    tagnames.push_back("BLOperator");
    tagnames.push_back("BLOperatorString");
    tagnames.push_back("GIOperator");
    tagnames.push_back("GIOperatorString");
    XMLHandler xmlres(xmltask,"ResultOperatorOrderedList");
    list<XMLHandler> opxml=xmlres.find_among_children(tagnames);
    vector<OperatorInfo> resultops;
    for (list<XMLHandler>::iterator ot=opxml.begin();ot!=opxml.end();++ot)
       resultops.push_back(OperatorInfo(*ot));
    uint nops=resultops.size();
    list<vector<pair<OperatorInfo,double> > > superposition;
    list<XMLHandler> corxml=xmltask.find_among_children("OperatorOrderedList");
    for (list<XMLHandler>::iterator ct=corxml.begin();ct!=corxml.end();++ct){
       vector<pair<OperatorInfo,double> > newcor;
       list<XMLHandler> itemxml=ct->find_among_children("Item");
       if (itemxml.size()!=nops) throw(std::invalid_argument("Mismatch in number of operators in CorrelatorMatrixSuperposition"));
       for (list<XMLHandler>::iterator it=itemxml.begin();it!=itemxml.end();++it){
          OperatorInfo opinfo(*it);
          double coef=1.0;
          xmlreadifchild(*it,"Coefficient",coef);
          newcor.push_back(make_pair(opinfo,coef));}
       superposition.push_back(newcor);}
    uint nterms=superposition.size();
    if (nterms<1) throw(std::invalid_argument("Zero terms in CorrelatorMatrixSuperposition"));
    uint tmin,tmax;
    xmlreadchild(xmltask,"MinimumTimeSeparation",tmin);
    xmlreadchild(xmltask,"MaximumTimeSeparation",tmax);
    bool herm=(xmltask.count("HermitianMatrix")>0) ? true: false;
    string filename;
    bool writetofile=xmlreadifchild(xmltask,"WriteToBinFile",filename);
    if (filename.empty()) writetofile=false;
    bool overwrite = false;  // protect mode
    if ((writetofile)&&(xml_tag_count(xmltask,"FileMode")==1)){
       string fmode;
       xmlread(xmltask,"FileMode",fmode,"FileListInfo");
       fmode=tidyString(fmode);
       if (fmode=="overwrite") overwrite=true;}

    xmlout.put_child("MinimumTimeSeparation",make_string(tmin));
    xmlout.put_child("MaximumTimeSeparation",make_string(tmax));
    if (herm) xmlout.put_child("HermitianMatrix");
    XMLHandler xmlo("ResultOperatorOrderedList");
    for (vector<OperatorInfo>::const_iterator it=resultops.begin();it!=resultops.end();it++){
       XMLHandler xmloo; it->output(xmloo); xmlo.put_child(xmloo);}
    xmlout.put_child(xmlo);
    for (list<vector<pair<OperatorInfo,double> > >::const_iterator
         st=superposition.begin();st!=superposition.end();st++){
       xmlo.set_root("OperatorOrderedList");
       for (vector<pair<OperatorInfo,double> >::const_iterator it=st->begin();it!=st->end();it++){
          XMLHandler xmloi("Item");
          XMLHandler xmloo; it->first.output(xmloo); xmloi.put_child(xmloo);
          xmloi.put_child("Coefficient",make_string(it->second));
          xmlo.put_child(xmloi);}
       xmlout.put_child(xmlo);}

    uint count=0;
    set<MCObsInfo> obskeys;
    for (uint row=0;row<nops;row++){
       for (uint col=(herm?row:0);col<nops;col++){
          CorrelatorAtTimeInfo resultcorr(resultops[row],resultops[col],0,herm,false);
          vector<CorrelatorAtTimeInfo> corrterms; vector<double> wts;
          for (list<vector<pair<OperatorInfo,double> > >::const_iterator
               st=superposition.begin();st!=superposition.end();st++){
             corrterms.push_back(CorrelatorAtTimeInfo((*st)[row].first,(*st)[col].first,0,herm,false));
             wts.push_back((*st)[row].second*(*st)[col].second);}
          for (uint t=tmin;t<=tmax;t++){
             resultcorr.resetTimeSeparation(t);
             MCObsInfo newkey(resultcorr);
             try{
             for (uint k=0;k<nterms;k++){
                corrterms[k].resetTimeSeparation(t);}
             const RVector& bins1=m_obs->getBins(MCObsInfo(corrterms[0]));
             RVector newbins(bins1);
             newbins*=wts[0];
             for (uint k=1;k<nterms;k++){
                const RVector& binsk=m_obs->getBins(MCObsInfo(corrterms[k]));
                RVector addbins(binsk);
                addbins*=wts[k];
                newbins+=addbins;}
             m_obs->putBins(newkey,newbins);
             if (writetofile) obskeys.insert(newkey);
             for (uint k=0;k<nterms;k++){
                m_obs->eraseData(MCObsInfo(corrterms[k]));}  // erase original data
             count++;}
             catch(const std::exception& xp){
                XMLHandler xmlk; newkey.output(xmlk);
                XMLHandler xmlerr("FailureToCompute");
                xmlerr.put_child(xmlk);
                xmlout.put_child(xmlerr);}
#ifdef COMPLEXNUMBERS
             try{
             const RVector& ibins1=m_obs->getBins(MCObsInfo(corrterms[0],ImaginaryPart));
             RVector newbins(ibins1);
             newbins*=wts[0];
             for (uint k=1;k<nterms;k++){
                const RVector& ibinsk=m_obs->getBins(MCObsInfo(corrterms[k],ImaginaryPart));
                RVector addbins(ibinsk);
                addbins*=wts[k];
                newbins+=addbins;}
             newkey.setToImaginaryPart();
             m_obs->putBins(newkey,newbins);
             if (writetofile) obskeys.insert(newkey);
             for (uint k=0;k<nterms;k++){
                m_obs->eraseData(MCObsInfo(corrterms[k],ImaginaryPart));}  // erase original data
             count++;}
             catch(const std::exception& xp){
                XMLHandler xmlk; newkey.output(xmlk);
                XMLHandler xmlerr("FailureToCompute");
                xmlerr.put_child(xmlk);
                xmlout.put_child(xmlerr);}
#endif
             }}}
       xmlout.put_child("NumberOfRealObservablesProcessed",make_string(count));
       if (writetofile){
          XMLHandler xmlf;
          m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
          xmlout.put_child("WriteBinsToFile",filename);
          xmlout.put_child(xmlf);}
       xmlout.put_child("Status","Done");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type CorrelatorMatrixSuperposition encountered an error: ")
             +string(errmsg.what())));} }

 /**************************************************************************************************************************************/


 else if (functype=="CorrelatorInteractionRatio"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","CorrelatorInteractionRatio");
    try{
    list<string> tagnames;
    tagnames.push_back("Operator");
    tagnames.push_back("OperatorString");
    tagnames.push_back("BLOperator");
    tagnames.push_back("BLOperatorString");
    tagnames.push_back("GIOperator");
    tagnames.push_back("GIOperatorString");
    XMLHandler xmlres(xmltask,"Result");
    OperatorInfo resultop(xmlres);

    XMLHandler xmlint(xmltask,"InteractingOperator");
    bool numvev=(xmlint.count("SubtractVEV")>0) ? true: false;
    pair<OperatorInfo,bool> numerator=make_pair(OperatorInfo(xmlint),numvev);
    vector<pair<OperatorInfo,bool> > denominator;
    list<XMLHandler> denomxml=xmltask.find_among_children("NonInteractingOperator");
    for (list<XMLHandler>::iterator it=denomxml.begin();it!=denomxml.end();++it){
       OperatorInfo opinfo(*it);
       bool subvev=(it->count("SubtractVEV")>0) ? true: false;
       denominator.push_back(make_pair(opinfo,subvev));}
    uint nterms=denominator.size();
    if (nterms<1) throw(std::invalid_argument("Zero NonInteractingOperators found"));
    uint tmin,tmax;
    xmlreadchild(xmltask,"MinimumTimeSeparation",tmin);
    xmlreadchild(xmltask,"MaximumTimeSeparation",tmax);
    // bool herm=(xmltask.count("HermitianMatrix")>0) ? true: false;
    bool herm=true;
    string filename;
    bool writetofile=xmlreadifchild(xmltask,"WriteToBinFile",filename);
    if (filename.empty()) writetofile=false;
    bool overwrite = false;  // protect mode
    if ((writetofile)&&(xml_tag_count(xmltask,"FileMode")==1)){
       string fmode;
       xmlread(xmltask,"FileMode",fmode,"FileListInfo");
       fmode=tidyString(fmode);
       if (fmode=="overwrite") overwrite=true;}

    xmlout.put_child("MinimumTimeSeparation",make_string(tmin));
    xmlout.put_child("MaximumTimeSeparation",make_string(tmax));
    // if (herm) xmlout.put_child("HermitianMatrix");
    XMLHandler xmlo("Result");
    resultop.output(xmlo);
    xmlout.put_child(xmlo);
    xmlo.set_root("NonInteractingOperators");
    for (vector<pair<OperatorInfo,bool> >::const_iterator
	   it=denominator.begin();it!=denominator.end();it++){
      XMLHandler xmloo; it->first.output(xmloo); xmlo.put_child(xmloo);}
    xmlout.put_child(xmlo);

    uint count=0;
    set<MCObsInfo> obskeys;
    CorrelatorAtTimeInfo resultcorr(resultop,resultop,0,herm,false);
    CorrelatorAtTimeInfo origcorrInfo(numerator.first,numerator.first,0,herm,false);
    // each tuple: <CorrelatorAtTimeInfo,vevFlag,vevOp>
    tuple<CorrelatorAtTimeInfo,bool,OperatorInfo> origcorr=make_tuple(origcorrInfo,numerator.second,numerator.first);
    vector<tuple<CorrelatorAtTimeInfo,bool,OperatorInfo> > denomcorrs;
    for (vector<pair<OperatorInfo,bool> >::const_iterator
           st=denominator.begin();st!=denominator.end();st++){
      bool subVEV=(*st).second;
      denomcorrs.push_back(make_tuple(CorrelatorAtTimeInfo((*st).first,(*st).first,0,herm,false),subVEV,(*st).first));}
    for (uint t=tmin;t<=tmax;t++){
      resultcorr.resetTimeSeparation(t);
      get<0>(origcorr).resetTimeSeparation(t);
      for (uint k=0;k<nterms;k++) get<0>(denomcorrs[k]).resetTimeSeparation(t);
      MCObsInfo newkey(resultcorr);
#ifdef COMPLEXNUMBERS
      try{
        const RVector& bins1=m_obs->getBins(MCObsInfo(get<0>(origcorr),RealPart));
        const RVector& ibins1=m_obs->getBins(MCObsInfo(get<0>(origcorr),ImaginaryPart));
        RVector newbins(bins1);
        RVector inewbins(ibins1);
        if(get<1>(origcorr)){
          const RVector& bins1VEV=m_obs->getBins(MCObsInfo(get<2>(origcorr),RealPart));
          const RVector& ibins1VEV=m_obs->getBins(MCObsInfo(get<2>(origcorr),ImaginaryPart));
          RVector origVEVbins(bins1VEV);
          RVector iorigVEVbins(ibins1VEV);
          origVEVbins*=origVEVbins;
          iorigVEVbins*=iorigVEVbins;
          newbins-=origVEVbins;
          newbins-=iorigVEVbins;}
        // loop over terms in denominator and divide
        for (uint k=0;k<nterms;k++){
	      const RVector& binsk=m_obs->getBins(MCObsInfo(get<0>(denomcorrs[k]),RealPart));
          const RVector& ibinsk=m_obs->getBins(MCObsInfo(get<2>(denomcorrs[k]),ImaginaryPart));
          RVector denombins(binsk);
          RVector idenombins(ibinsk);
          if(get<1>(denomcorrs[k])){
            const RVector& binskVEV=m_obs->getBins(MCObsInfo(get<2>(denomcorrs[k]),RealPart));
            const RVector& ibinskVEV=m_obs->getBins(MCObsInfo(get<2>(denomcorrs[k]),ImaginaryPart));
            RVector denomVEVbins(binskVEV);
            RVector idenomVEVbins(ibinskVEV);
            denomVEVbins*=denomVEVbins;
            idenomVEVbins*=idenomVEVbins;
            denombins-=denomVEVbins;
            denombins-=idenomVEVbins;}
          // now divide newbins by denombins
          for (uint j=0;j<newbins.size();j++){
            double norm=denombins[j]*denombins[j]+idenombins[j]*idenombins[j];
            double res=newbins[j]*denombins[j]+inewbins[j]*idenombins[j];
            double ires=inewbins[j]*denombins[j]-newbins[j]*idenombins[j];
            res/=norm; ires/=norm;
            newbins[j]=res;
            inewbins[j]=ires;}}
        m_obs->putBins(newkey,newbins);
        if (writetofile) obskeys.insert(newkey);
        newkey.setToImaginaryPart();
        m_obs->putBins(newkey,inewbins);
        if (writetofile) obskeys.insert(newkey);
        for (uint k=0;k<nterms;k++){
          m_obs->eraseData(MCObsInfo(get<0>(denomcorrs[k]),RealPart));  // erase original data
          m_obs->eraseData(MCObsInfo(get<0>(denomcorrs[k]),ImaginaryPart));}
        m_obs->eraseData(MCObsInfo(get<0>(origcorr),RealPart));
        m_obs->eraseData(MCObsInfo(get<0>(origcorr),ImaginaryPart));
        count+=2;}
      catch(const std::exception& xp){
        XMLHandler xmlk; newkey.output(xmlk);
        XMLHandler xmlerr("FailureToCompute");
        xmlerr.put_child(xmlk);
        xmlout.put_child(xmlerr);}
#else
      try{
        const RVector& bins1=m_obs->getBins(MCObsInfo(get<0>(origcorr),RealPart));
        RVector newbins(bins1);
        if (get<1>(origcorr)){
          const RVector& bins1VEV=m_obs->getBins(MCObsInfo(get<2>(origcorr),RealPart));
          RVector origVEVbins(bins1VEV);
          origVEVbins*=origVEVbins;
          newbins-=origVEVbins;}
        for (uint k=0;k<nterms;k++){
	      const RVector& binsk=m_obs->getBins(MCObsInfo(get<0>(denomcorrs[k]),RealPart));
          RVector denombins(binsk);
          if(get<1>(denomcorrs[k])){
            const RVector& binskVEV=m_obs->getBins(MCObsInfo(get<2>(denomcorrs[k]),RealPart));
            RVector denomVEVbins(binskVEV);
            denomVEVbins*=denomVEVbins;
            denombins-=denomVEVbins;}
          newbins/=denombins;}
        m_obs->putBins(newkey,newbins);
        if (writetofile) obskeys.insert(newkey);
        for (uint k=0;k<nterms;k++)
          m_obs->eraseData(MCObsInfo(get<0>(denomcorrs[k])));
        m_obs->eraseData(MCObsInfo(get<0>(origcorr))); // erase original data
        count++;}
      catch(const std::exception& xp){
      	XMLHandler xmlk; newkey.output(xmlk);
      	XMLHandler xmlerr("FailureToCompute");
      	xmlerr.put_child(xmlk);
      	xmlout.put_child(xmlerr);}
#endif
    }
    // erase VEV data when no longer needed
    if(get<1>(origcorr)){
      m_obs->eraseData(MCObsInfo(get<2>(origcorr),RealPart));
#ifdef COMPLEXNUMBERS
      m_obs->eraseData(MCObsInfo(get<2>(origcorr),ImaginaryPart));
#endif
    }
    for (uint k=0;k<nterms;k++){
      if(get<1>(denomcorrs[k])){
        m_obs->eraseData(MCObsInfo(get<2>(denomcorrs[k]),RealPart));
#ifdef COMPLEXNUMBERS
        m_obs->eraseData(MCObsInfo(get<2>(denomcorrs[k]),ImaginaryPart));
#endif
      }}
       xmlout.put_child("NumberOfRealObservablesProcessed",make_string(count));
       if (writetofile){
          XMLHandler xmlf;
          m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
          xmlout.put_child("WriteBinsToFile",filename);
          xmlout.put_child(xmlf);}
       xmlout.put_child("Status","Done");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type CorrelatorInteractionRatio encountered an error: ")
             +string(errmsg.what())));} }

 else{
    throw(std::invalid_argument("DoObsFunction encountered unsupported function: "));}

}


// ***************************************************************************************
