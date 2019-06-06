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
// *      Notes: 1) the resultant correlator will have (if specified) all        *
// *             VEVs subtracted in this task.                                   *
// *             2) correlator ratio is a non-simple observable so is only       *
// *             available as resamplings.                                       *
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
// *          <SubtractVEV/>    (optional)                                       *
// *       </InteractingOperator>                                                *
// *       <NonInteractingOperator>                                              *
// *          <Operator>...</Operator>                                           *
// *          <SubtractVEV/>    (optional)                                       *
// *       </NonInteractingOperator>                                             *
// *          ......                                                             *
// *       <NonInteractingOperator>                                              *
// *          <Operator>...</Operator>                                           *
// *          <SubtractVEV/>    (optional)                                       *
// *       </NonInteractingOperator>                                             *
// *       <Reweight/>                                                           *
// *       <MinimumTimeSeparation>0</MinimumTimeSeparation>                      *
// *       <MaximumTimeSeparation>15</MaximumTimeSeparation>                     *
// *       <WriteToSamplingsFile>filename</WriteToSamplingsFile>   (optional)    *
// *       <SamplingMode>Bootstrap</SamplingMode>    (or Jackknife)              *
// *       <FileMode>overwrite</FileMode>   (optional)                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      You may want to obtain the total energy after having extracted         *
// *      the energy difference from a ratio fit. You can use the following      *
// *      task for this.                                                         *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>ReconstructEnergy</Type>                                        *
// *       <SpatialExtentNumSites>32</SpatialExtentNumSites>                     *
// *       <Anisotropy>       (optional)                                         *
// *          <Name>aniso_fit_name</Name>                                        *
// *          <IDIndex>0</IDIndex>                                               *
// *       </Anisotropy>                                                         *
// *       <Result>                                                              *
// *          <Name>result-name</Name>                                           *
// *          <IDIndex>0</IDIndex>                                               *
// *       </Result>                                                             *
// *       <EnergyDifferenceFit>                                                 *
// *          <Name>diff_fit_name</Name>                                         *
// *          <IDIndex>0</IDIndex>                                               *
// *       </EnergyDifferenceFit>                                                *
// *       <ScatteringParticleEnergyFit>                                         *
// *          <IntMomSquared>4</IntMomSquared>                                   *
// *          <Name>scatting_part_atrest_energy_fit_name</Name>                  *
// *          <IDIndex>0</IDIndex>                                               *
// *       </ScatteringParticleEnergyFit>                                        *
// *       <ScatteringParticleEnergyFit>                                         *
// *          <IntMomSquared>2</IntMomSquared>                                   *
// *          <Name>scatting_part_atrest_energy_fit_name</Name>                  *
// *          <IDIndex>0</IDIndex>                                               *
// *       </ScatteringParticleEnergyFit>                                        *
// *          ...                                                                *
// *       <Mode>Jackknife</Mode> (optional)                                     *
// *                      (or Bootstrap or Current [default] or Bins )           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      You may want to obtain the the amplitude of the original correlator    *
// *      after obtaining a fit to the amplitude of the ratio correlator. You    *
// *      can use the following task for this.                                   *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>ReconstructAmplitude</Type>                                     *
// *       <Result>                                                              *
// *          <Name>result-name</Name>                                           *
// *          <IDIndex>0</IDIndex>                                               *
// *       </Result>                                                             *
// *       <EnergyDifferenceAmplitudeFit>                                        *
// *          <Name>diff_fit_name</Name>                                         *
// *          <IDIndex>0</IDIndex>                                               *
// *       </EnergyDifferenceAmplitudeFit>                                       *
// *       <ScatteringParticleAmplitudeFit>                                      *
// *          <Name>scatting_part_atrest_amp_fit_name</Name>                     *
// *          <IDIndex>0</IDIndex>                                               *
// *       </ScatteringParticleAmplitudeFit>                                     *
// *       <ScatteringParticleAmplitudeFit>                                      *
// *          <Name>scatting_part_atrest_amp_fit_name</Name>                     *
// *          <IDIndex>0</IDIndex>                                               *
// *       </ScatteringParticleAmplitudeFit>                                     *
// *          ...                                                                *
// *       <Mode>Jackknife</Mode> (optional)                                     *
// *                      (or Bootstrap or Current [default] or Bins )           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>Copy</Type>                                                     *
// *       <From>                                                                *
// *            <MCObservable> ... </MCObservable>                               *
// *       </From>                                                               *
// *       <To>                                                                  *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </To>                                                                 *
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
          CorrelatorAtTimeInfo origcorr(*oldrow,*oldcol,0,herm,false,false);
          CorrelatorAtTimeInfo newcorr(*newrow,*newcol,0,herm,false,false);
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
          CorrelatorAtTimeInfo resultcorr(resultops[row],resultops[col],0,herm,false,false);
          vector<CorrelatorAtTimeInfo> corrterms; vector<double> wts;
          for (list<vector<pair<OperatorInfo,double> > >::const_iterator
               st=superposition.begin();st!=superposition.end();st++){
             corrterms.push_back(CorrelatorAtTimeInfo((*st)[row].first,(*st)[col].first,0,herm,false,false));
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
     bool reweight=(xmltask.count("Reweight")>0) ? true: false;
     // bool herm=(xmltask.count("HermitianMatrix")>0) ? true: false;
     bool herm=true;
     string filename;
     bool writetofile=xmlreadifchild(xmltask,"WriteToSamplingsFile",filename);
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
     SamplingMode mode=Jackknife;
     string instr;
     if (xmlreadifchild(xmltask,"SamplingMode",instr)){
       if (instr=="Bootstrap") mode=Bootstrap;
       else if (instr=="Jackknife") mode=Jackknife;
       else throw(std::invalid_argument("Bad sampling mode"));}
     if (mode==Bootstrap){
       xmlout.put_child("SamplingMode","Bootstrap");
       m_obs->setToBootstrapMode();}
     else{
       xmlout.put_child("SamplingMode","Jackknife");
       m_obs->setToJackknifeMode();}

     XMLHandler xmlo, xmlp;
     xmlo.set_root("Result");
     resultop.output(xmlp);
     xmlo.put_child(xmlp);
     xmlout.put_child(xmlo);
     xmlo.set_root("InteractingOperator");
     numerator.first.output(xmlp);
     xmlo.put_child(xmlp);
     xmlout.put_child(xmlo);
     xmlo.set_root("NonInteractingOperators");
     for (vector<pair<OperatorInfo,bool> >::const_iterator
            it=denominator.begin();it!=denominator.end();it++){
       XMLHandler xmloo; it->first.output(xmloo); xmlo.put_child(xmloo);}
     xmlout.put_child(xmlo);

     uint count=0;
     set<MCObsInfo> obskeys;
     CorrelatorAtTimeInfo resultcorr(resultop,resultop,0,herm,false,reweight);
     CorrelatorAtTimeInfo origcorr(numerator.first,numerator.first,0,herm,numerator.second,reweight);
     vector<CorrelatorAtTimeInfo> denomcorrs;
     for (vector<pair<OperatorInfo,bool> >::const_iterator
            st=denominator.begin();st!=denominator.end();st++){
       denomcorrs.push_back(CorrelatorAtTimeInfo(st->first,st->first,0,herm,st->second,reweight));}
     for (uint t=tmin;t<=tmax;t++){
       resultcorr.resetTimeSeparation(t);
       origcorr.resetTimeSeparation(t);
       for (uint k=0;k<nterms;k++) denomcorrs[k].resetTimeSeparation(t);
       MCObsInfo newkey(resultcorr,RealPart);
       try{
         for (m_obs->begin();!m_obs->end();++(*m_obs)){
           double newval=m_obs->getCurrentSamplingValue(MCObsInfo(origcorr,RealPart));
           for (uint k=0;k<nterms;k++){
             double denomval=m_obs->getCurrentSamplingValue(MCObsInfo(denomcorrs[k],RealPart));
             newval/=denomval;}
           m_obs->putCurrentSamplingValue(newkey,newval,true);}

         if (writetofile) obskeys.insert(newkey);
         count++;}
       catch(const std::exception& xp){
         XMLHandler xmlk; newkey.output(xmlk);
         XMLHandler xmlerr("FailureToCompute");
         xmlerr.put_child(xmlk);
         xmlout.put_child(xmlerr);}
     }

     xmlout.put_child("NumberOfRealObservablesProcessed",make_string(count));
     if (writetofile){
       XMLHandler xmlf;
       m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,overwrite);
       xmlout.put_child("WriteSamplingsToFile",filename);
       xmlout.put_child(xmlf);}
     xmlout.put_child("Status","Done");}
   catch(const std::exception& errmsg){
     xmlout.clear();
     throw(std::invalid_argument(string("DoObsFunction with type CorrelatorInteractionRatio encountered an error: ")
                                 +string(errmsg.what())));} }

 else if (functype=="ReconstructEnergy"){
   xmlout.set_root("DoObsFunction");
   xmlout.put_child("Type","ReconstructEnergy");
   try{
     XMLHandler xml_diff_fit(xmltask,"EnergyDifferenceFit");
     ArgsHandler arg_diff_fit(xml_diff_fit);
     string diff_fit_name = arg_diff_fit.getString("Name");
     uint diff_fit_index=0;
     arg_diff_fit.getOptionalUInt("IDIndex",diff_fit_index);
     MCObsInfo diff_fit(diff_fit_name,diff_fit_index);

     XMLHandler xmlt1,xmlt2;
     xmlt1.set_root("EnergyDifference");
     diff_fit.output(xmlt2);
     xmlt1.put_child(xmlt2);
     xmlout.put_child(xmlt1);

     uint aniscount=xmltask.count("Anisotropy");
     MCObsInfo* obsxi=0;
     if (aniscount==1){
       XMLHandler xmlxi(xmltask,"Anisotropy");
       ArgsHandler arg_aniso(xmlxi);
       string aniso_name = arg_aniso.getString("Name");
       uint aniso_index=0;
       arg_aniso.getOptionalUInt("IDIndex",aniso_index);
       obsxi = new MCObsInfo(aniso_name,aniso_index);
       xmlt1.set_root("Anisotropy");
       obsxi->output(xmlt2);
       xmlt1.put_child(xmlt2);
       xmlout.put_child(xmlt1);}

     uint m_lat_spatial_extent;
     xmlreadifchild(xmltask,"SpatialExtentNumSites",m_lat_spatial_extent);
     if (m_lat_spatial_extent<1)
       throw(std::invalid_argument("Lattice spatial extent must be a positive integer"));

     list<XMLHandler> xmlscats=xmltask.find_among_children("ScatteringParticleEnergyFit");
     list<pair<MCObsInfo,double> > scattering_particles;
     for (list<XMLHandler>::iterator it=xmlscats.begin();it!=xmlscats.end();it++){
        ArgsHandler xmlscat(*it);
        uint psq=xmlscat.getUInt("IntMomSquared");

        double m_momsq_quantum=6.2831853071795864770/double(m_lat_spatial_extent);
        m_momsq_quantum*=m_momsq_quantum;
        double psqfactor=psq*m_momsq_quantum;

        string name(xmlscat.getString("Name"));
        uint index=0;
        xmlscat.getOptionalUInt("IDIndex",index);
        MCObsInfo scat_info(name,index);

        xmlt1.set_root("ScatteringParticleEnergyFit");
        xmlt1.put_child("IntMomSquared", to_string(psq));
        scat_info.output(xmlt2);
        xmlt1.put_child(xmlt2);
        xmlout.put_child(xmlt1);

        scattering_particles.push_back(make_pair(scat_info,psqfactor));}

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

     SamplingMode origmode=m_obs->getCurrentSamplingMode();
     if (mcode=='J') m_obs->setToJackknifeMode();
     else m_obs->setToBootstrapMode();

     XMLHandler xmlres(xmltask,"Result");
     string name; int index;
     xmlreadchild(xmlres,"Name",name);
     if (name.empty()) throw(std::invalid_argument("Must provide name for Energy result"));
     index=taskcount;
     xmlreadifchild(xmlres,"IDIndex",index);
     MCObsInfo resinfo(name,index,mcode=='D');
     xmlt1.set_root("ResultInfo");
     resinfo.output(xmlt2);
     xmlt1.put_child(xmlt2);
     xmlout.put_child(xmlt1);

     if (aniscount==1)
       doReconstructEnergyBySamplings(*m_obs,diff_fit,*obsxi,scattering_particles,resinfo);
     else
       doReconstructEnergyBySamplings(*m_obs,diff_fit,scattering_particles,resinfo);

     if (obsxi!=0) delete obsxi;

     MCEstimate est=m_obs->getEstimate(resinfo);
     est.output(xmlt1);
     xmlout.put_child(xmlt1);
     m_obs->setSamplingMode(origmode);}
   catch(const std::exception& errmsg){
     xmlout.clear();
     throw(std::invalid_argument(string("DoObsFunction with type ReconstructEnergy encountered an error: ")
                                 +string(errmsg.what())));} }

 else if (functype=="ReconstructAmplitude"){
   xmlout.set_root("DoObsFunction");
   xmlout.put_child("Type","ReconstructAmplitude");
   try{
     XMLHandler xml_diff_fit(xmltask,"EnergyDifferenceAmplitudeFit");
     ArgsHandler arg_diff_fit(xml_diff_fit);
     string diff_fit_name = arg_diff_fit.getString("Name");
     uint diff_fit_index=0;
     arg_diff_fit.getOptionalUInt("IDIndex",diff_fit_index);
     MCObsInfo diff_fit(diff_fit_name,diff_fit_index);

     XMLHandler xmlt1,xmlt2;
     xmlt1.set_root("EnergyDifferenceAmplitude");
     diff_fit.output(xmlt2);
     xmlt1.put_child(xmlt2);
     xmlout.put_child(xmlt1);

     list<XMLHandler> xmlscatamps=xmltask.find_among_children("ScatteringParticleAmplitudeFit");
     list<MCObsInfo> scattering_particles_amps;
     for (list<XMLHandler>::iterator it=xmlscatamps.begin();it!=xmlscatamps.end();it++){
        ArgsHandler xmlscatamp(*it);
        string name(xmlscatamp.getString("Name"));
        uint index=0;
        xmlscatamp.getOptionalUInt("IDIndex",index);
        MCObsInfo scat_info(name,index);

        xmlt1.set_root("ScatteringParticleAmplitudeFit");
        scat_info.output(xmlt2);
        xmlt1.put_child(xmlt2);
        xmlout.put_child(xmlt1);

        scattering_particles_amps.push_back(scat_info);}

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
     if (name.empty()) throw(std::invalid_argument("Must provide name for Amplitude result"));
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

     doReconstructAmplitudeBySamplings(*m_obs,diff_fit,scattering_particles_amps,resinfo);

     MCEstimate est=m_obs->getEstimate(resinfo);
     est.output(xmlt1);
     xmlout.put_child(xmlt1);
     m_obs->setSamplingMode(origmode);}
   catch(const std::exception& errmsg){
     xmlout.clear();
     throw(std::invalid_argument(string("DoObsFunction with type ReconstructAmplitude encountered an error: ")
                                 +string(errmsg.what())));} }

 else if (functype=="Copy"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","Copy");
    try{
    XMLHandler xmlfrom(xmltask,"From");
    XMLHandler xmlt1,xmlt2;
    MCObsInfo obsfrom(xmlfrom);
    xmlt1.set_root("From");
    obsfrom.output(xmlt2);
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

    XMLHandler xmlto(xmltask,"To");
    string name; int index;
    xmlreadchild(xmlto,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for Copy result"));
    index=taskcount;
    xmlreadifchild(xmlto,"IDIndex",index);
    MCObsInfo toinfo(name,index,mcode=='D');
    xmlt1.set_root("ToInfo");
    toinfo.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    if (mcode=='D'){
       doCopyByBins(*m_obs,obsfrom,toinfo);
       MCEstimate est=m_obs->getEstimate(toinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       doCopyBySamplings(*m_obs,obsfrom,toinfo);
       MCEstimate est=m_obs->getEstimate(toinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type Copy encountered an error: ")
                                   +string(errmsg.what())));}
    }

 else{
   throw(std::invalid_argument("DoObsFunction encountered unsupported function: "));}

}


// ***************************************************************************************
