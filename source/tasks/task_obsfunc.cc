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
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife or bins )                   *
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
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife or bins )                   *
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
// *       <Mode>samplings</Mode> ( or bins )                                    *
// *       <WriteToFile>                                                         *
// *          <FileName>name</FileName>                                          *
// *          <FileType>bins</FileType> (or samplings)                           *
// *          <FileFormat>fstr</FileFormat> (or hdf5: default if absent)         *
// *          <WriteMode>overwrite</WriteMode> (protect, update, overwrite)      *
// *       </WriteToFile>                                                        *
// *    </Task>                                                                  *
// *             <Mode> and <FileType> must be same if both present;             *
// *             at least one of these tags must be present                      *
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
// *       <Mode>samplings</Mode> ( or bins )                                    *
// *       <WriteToFile>                                                         *
// *          <FileName>name</FileName>                                          *
// *          <FileType>bins</FileType> (or samplings)                           *
// *          <FileFormat>fstr</FileFormat> (or hdf5: default if absent)         *
// *          <WriteMode>overwrite</WriteMode> (protect, update, overwrite)      *
// *       </WriteToFile>                                                        *
// *       <IgnoreMissing/>  (optional to ignore missing correlators)            *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *                                                                             *
// *      Performs the following transformation on a correlation matrix C:       *
// *                                                                             *
// *         Ctrans[i,j] = T^dagger[i,k] C[k,m] T[m,j]                           *
// *                                                                             *
// *      If C is NxN, then the transformation matrix T must be NXM, where       *
// *      M <= N.  In terms of operators, this means the transformed operators   *
// *      W[k] are given in terms of the original operator O[m] by               *
// *                Wbar[j] = Obar[m] T[m,j]                                     *
// *                   W[i] = T^dag[i,k] O[k]                                    *
// *      The input XML must have the form                                       *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>TransformCorrelatorMatrix</Type>                                *
// *         <TransformedOperator>   (creation operator)                         *
// *            <OpName><Operator>...</Operator></OpName>                        *
// *            <OpTerm>                                                         *
// *                 <Operator>...</Operator>                                    *
// *                 <Coefficient>(0.43,-0.121)</Coefficient>                    *
// *            </OpTerm>                                                        *
// *             ... other opterms                                               *
// *         </TransformedOperator>                                              *
// *            ...   other transformed operators                                *
// *       <MinimumTimeSeparation>3</MinimumTimeSeparation>                      *
// *       <MaximumTimeSeparation>12</MaximumTimeSeparation>                     *
// *       <HermitianMatrix/>  (if hermitian)                                    *
// *       <SubtractVEV/> (optional)                                             *
// *       <Mode>samplings</Mode> ( or bins )                                    *
// *       <WriteToFile>                                                         *
// *          <FileName>name</FileName>                                          *
// *          <FileType>bins</FileType> (or samplings)                           *
// *          <FileFormat>fstr</FileFormat> (or hdf5: default if absent)         *
// *          <WriteMode>overwrite</WriteMode> (protect, update, overwrite)      *
// *          <SeparateVEVWrite/> (optional for FileType=samplings: default no)  *
// *       </WriteToFile>                                                        *
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
// *       <Anisotropy><MCObservable> ... </MCObservable></Anisotropy>           *
// *       <ReferenceEnergy>   (optional)                                        *
// *          <MCObservable> ... </MCObservable>                                 *
// *       </ReferenceEnergy>                                                    *
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife or bins )                   *
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
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife )                           *
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
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife )                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      You may want to obtain the energy difference after having extracted    *
// *      the absolute energy from a non-ratio fit. You can use the following    *
// *      task for this.                                                         *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>EnergyDifference</Type>                                         *
// *       <SpatialExtentNumSites>32</SpatialExtentNumSites>                     *
// *       <Anisotropy>       (optional)                                         *
// *          <Name>aniso_fit_name</Name>                                        *
// *          <IDIndex>0</IDIndex>                                               *
// *       </Anisotropy>                                                         *
// *       <Result>                                                              *
// *          <Name>result-name</Name>                                           *
// *          <IDIndex>0</IDIndex>                                               *
// *       </Result>                                                             *
// *       <EnergyFit>                                                           *
// *          <Name>fit_name</Name>                                              *
// *          <IDIndex>0</IDIndex>                                               *
// *       </EnergyFit>                                                          *
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
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife )                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *      For performing an exponential of a function                            *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>Exp</Type>                                                      *
// *       <Result>                                                              *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </Result>                                                             *
// *       <InObservable><MCObservable> ... </MCObservable></InObservable>       *
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife or bins )                   *
// *    </Task>                                                                  *
// *                                                                             *
// *      For performing the log of a function                                   *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>Log</Type>                                                      *
// *       <Result>                                                              *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </Result>                                                             *
// *       <InObservable><MCObservable> ... </MCObservable></InObservable>       *
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife or bins )                   *
// *    </Task>                                                                  *
// *                                                                             *
// *      For copying data from one observable to another                        *
// *    <Task>                                                                   *
// *     <Action>DoObsFunction</Action>                                          *
// *       <Type>Copy</Type>                                                     *
// *       <Result>                                                              *
// *          <Name>result-name</Name><IDIndex>0</IDIndex>                       *
// *       </Result>                                                             *
// *       <InObservable><MCObservable> ... </MCObservable></InObservable>       *
// *       <Mode>samplings</Mode> (default: current sampling method)             *
// *                      (or Bootstrap or Jackknife or bins )                   *
// *    </Task>                                                                  *
// *                                                                             *
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

    string datamode="samplings";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="samplings"){
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

    string datamode="samplings";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="samplings"){
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

     string datamode="samplings";
     xmlreadifchild(xmltask,"Mode",datamode);
     char mcode;
     if (datamode=="Bootstrap") mcode='B';
     else if (datamode=="Jackknife") mcode='J';
     else if (datamode=="samplings"){
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
    string datamode="none";
    xmlreadifchild(xmltask,"Mode",datamode);
    if ((datamode!="bins")&&(datamode!="samplings")&&(datamode!="none"))
       throw(std::invalid_argument("Invalid Sampling Mode in CorrelatorMatrixTimeDifferences"));
    xmlout.put_child("Mode",datamode);
    string ftype(datamode),filename;
    string fformat("default"); char ffmt='D';
    WriteMode wmode = Protect;  // protect mode
    bool writetofile = false;
    if (xmltask.count("WriteToFile")==1){
       XMLHandler xmlw(xmltask,"WriteToFile");
       xmlreadifchild(xmlw,"FileType",ftype);
       if ((ftype!="bins")&&(ftype!="samplings"))
          throw(std::invalid_argument("<FileType> must be bins or samplings in CorrelatorMatrixTimeDifferences"));
       if ((ftype!=datamode)&&(datamode!="none"))
          throw(std::invalid_argument("<Mode> and <FileType> must match in CorrelatorMatrixTimeDifferences"));
       xmlreadifchild(xmlw,"FileFormat",fformat);
       if (fformat=="fstr") ffmt='F';
       else if (fformat=="hdf5") ffmt='H';
       else if (fformat=="default") ffmt='D';
       else throw(std::invalid_argument("<FileFormat> must be ftr or hdf5 or default in CorrelatorMatrixTimeDifferences"));
       xmlreadchild(xmlw,"FileName",filename,"TaskHandler");
       if (xml_tag_count(xmltask,"WriteMode")==1){
          string fmode;
          xmlread(xmltask,"WriteMode",fmode,"CorrelatorMatrixTimeDifferences");
          fmode=tidyString(fmode);
          if (fmode=="overwrite") wmode=Overwrite;
          else if (fmode=="update") wmode=Update;}
       writetofile=true;}
    if (ftype=="none")
       throw(std::invalid_argument("At least one of <Mode> and <FileType> must appear in CorrelatorMatrixTimeDifferences"));
    xmlout.put_child("Mode",ftype);
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
    bool erase_orig=true;
    set<MCObsInfo> obskeys;
    if (ftype=="bins")
       doCorrelatorMatrixTimeDifferencesByBins(*m_obs,origops,newops,herm,tmin,tmax,
                                               obskeys,erase_orig);
    else
       doCorrelatorMatrixTimeDifferencesBySamplings(*m_obs,origops,newops,herm,tmin,tmax,
                                                    obskeys,erase_orig);
    xmlout.put_child("NumberOfRealObservablesProcessed",make_string(obskeys.size()));
    if (writetofile){
       XMLHandler xmlf;
       xmlout.put_child("FileType",ftype);
       if (ftype=="bins")
          m_obs->writeBinsToFile(obskeys,filename,xmlf,wmode,ffmt);
       else
          m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,wmode,ffmt);
       xmlout.put_child("WriteToFile",filename);}
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
    bool ignore_missing=(xmltask.count("IgnoreMissing")>0) ? true: false;
    bool herm=(xmltask.count("HermitianMatrix")>0) ? true: false;
    string datamode="none";
    xmlreadifchild(xmltask,"Mode",datamode);
    if ((datamode!="bins")&&(datamode!="samplings")&&(datamode!="none"))
       throw(std::invalid_argument("Invalid Sampling Mode in CorrelatorMatrixSuperposition"));
    string ftype(datamode),filename;
    string fformat("default"); char ffmt='D';
    WriteMode wmode = Protect;  // protect mode
    bool writetofile = false;
    if (xmltask.count("WriteToFile")==1){
       XMLHandler xmlw(xmltask,"WriteToFile");
       xmlreadifchild(xmlw,"FileType",ftype);
       if ((ftype!="bins")&&(ftype!="samplings"))
          throw(std::invalid_argument("<FileType> must be bins or samplings in CorrelatorMatrixSuperposition"));
       if ((ftype!=datamode)&&(datamode!="none"))
          throw(std::invalid_argument("<Mode> and <FileType> must match in CorrelatorMatrixSuperposition"));
       xmlreadifchild(xmlw,"FileFormat",fformat);
       if (fformat=="fstr") ffmt='F';
       else if (fformat=="hdf5") ffmt='H';
       else if (fformat=="default") ffmt='D';
       else throw(std::invalid_argument("<FileFormat> must be ftr or hdf5 or default in CorrelatorMatrixSuperposition"));
       xmlreadchild(xmlw,"FileName",filename,"TaskHandler");
       if (xml_tag_count(xmltask,"WriteMode")==1){
          string fmode;
          xmlread(xmltask,"WriteMode",fmode,"CorrelatorMatrixSuperposition");
          fmode=tidyString(fmode);
          if (fmode=="overwrite") wmode=Overwrite;
          else if (fmode=="update") wmode=Update;}
       writetofile=true;}
    if (ftype=="none")
       throw(std::invalid_argument("At least one of <Mode> and <FileType> must appear in CorrelatorMatrixSuperposition"));
    xmlout.put_child("Mode",ftype);
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
    set<MCObsInfo> obskeys;
    bool erase_orig=true;
    if (ftype=="bins"){
       doCorrelatorMatrixSuperpositionByBins(*m_obs,superposition,resultops,herm,
                                             tmin,tmax,obskeys,erase_orig,ignore_missing);}
    else{
       doCorrelatorMatrixSuperpositionBySamplings(*m_obs,superposition,resultops,herm,
                                                  tmin,tmax,obskeys,erase_orig,ignore_missing);}
    xmlout.put_child("NumberOfRealObservablesProcessed",make_string(obskeys.size()));
    for (set<MCObsInfo>::const_iterator kt=obskeys.begin();kt!=obskeys.end();++kt){
       xmlout.put_child("ProcessedKey",kt->str());}
    if (writetofile){
       XMLHandler xmlf;
       xmlout.put_child("FileType",ftype);
       if (ftype=="bins")
          m_obs->writeBinsToFile(obskeys,filename,xmlf,wmode,ffmt);
       else
          m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,wmode,ffmt);
       xmlout.put_child("WriteToFile",filename);}
    xmlout.put_child("Status","Done");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type CorrelatorMatrixSuperposition encountered an error: ")
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

     string datamode="samplings";
     xmlreadifchild(xmltask,"Mode",datamode);
     char mcode;
     if (datamode=="Bootstrap") mcode='B';
     else if (datamode=="Jackknife") mcode='J';
     else if (datamode=="samplings"){
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

     string datamode="samplings";
     xmlreadifchild(xmltask,"Mode",datamode);
     char mcode;
     if (datamode=="Bootstrap") mcode='B';
     else if (datamode=="Jackknife") mcode='J';
     else if (datamode=="samplings"){
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

 else if (functype=="EnergyDifference"){
   xmlout.set_root("DoObsFunction");
   xmlout.put_child("Type","EnergyDifference");
   try{
     XMLHandler xml_fit(xmltask,"EnergyFit");
     ArgsHandler arg_fit(xml_fit);
     string fit_name = arg_fit.getString("Name");
     uint fit_index=0;
     arg_fit.getOptionalUInt("IDIndex",fit_index);
     MCObsInfo fit(fit_name,fit_index);

     XMLHandler xmlt1,xmlt2;
     xmlt1.set_root("Energy");
     fit.output(xmlt2);
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

     string datamode="samplings";
     xmlreadifchild(xmltask,"Mode",datamode);
     char mcode;
     if (datamode=="Bootstrap") mcode='B';
     else if (datamode=="Jackknife") mcode='J';
     else if (datamode=="samplings"){
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
       doEnergyDifferenceBySamplings(*m_obs,fit,*obsxi,scattering_particles,resinfo);
     else
       doEnergyDifferenceBySamplings(*m_obs,fit,scattering_particles,resinfo);

     if (obsxi!=0) delete obsxi;

     MCEstimate est=m_obs->getEstimate(resinfo);
     est.output(xmlt1);
     xmlout.put_child(xmlt1);
     m_obs->setSamplingMode(origmode);}
   catch(const std::exception& errmsg){
     xmlout.clear();
     throw(std::invalid_argument(string("DoObsFunction with type ReconstructEnergy encountered an error: ")
                                 +string(errmsg.what())));} }




 else if (functype=="TransformCorrelatorMatrix"){
    xmlout.set_root("DoObsFunction");
    LogHelper logger("TransformCorrelatorMatrix");
    ArgsHandler gin(xmltask);
    try{
    map<OperatorInfo,list<pair<OperatorInfo,Scalar> > > transops;
    set<OperatorInfo> origops;
    list<ArgsHandler> glist=gin.getSubHandlers("TransformedOperator");
    for (list<ArgsHandler>::iterator gt=glist.begin();gt!=glist.end();++gt){
       ArgsHandler gop(*gt,"OpName");
       LogHelper glog("TransformedOperator");
       OperatorInfo transop;
       gop.getItem("TransformedOpName",transop);
       if (transops.find(transop)!=transops.end())
          throw(std::invalid_argument("Duplicate TransformedOperator in TransformCorrelatorMatrix"));
       list<pair<OperatorInfo,Scalar> >& todef=transops[transop];  // inserts new element into map
       list<ArgsHandler> gglist=gt->getSubHandlers("OpTerm");
       for (list<ArgsHandler>::iterator ggt=gglist.begin();ggt!=gglist.end();++ggt){
          OperatorInfo opterm;
          ggt->getItem("OrigOperator",opterm);
          origops.insert(opterm);
          Scalar coef;
          ggt->getScalar("Coefficient",coef);
          glog.putEcho(*ggt,"OpTerm");
          todef.push_back(make_pair(opterm,coef));}
       logger.put(glog);}
         // for the vectors "origopsvec" and "transopsvec" and the transformation matrix "T"
    vector<OperatorInfo> origopsvec(origops.begin(),origops.end());
    vector<OperatorInfo> transopsvec;
    uint norig=origopsvec.size();
    uint ntrans=transops.size();
    TransMatrix T(norig,ntrans,0.0);
    uint kk=0;
    for (map<OperatorInfo,list<pair<OperatorInfo,Scalar> > >::const_iterator 
         tt=transops.begin();tt!=transops.end();++tt,++kk){
       transopsvec.push_back(tt->first);
       for (list<pair<OperatorInfo,Scalar> >::const_iterator st=tt->second.begin();st!=tt->second.end();++st){
          vector<OperatorInfo>::iterator ft = std::find(origopsvec.begin(), origopsvec.end(),st->first);
          if (ft==origopsvec.end())
             throw(std::invalid_argument("Could not construct transformation matrix in TransformCorrelatorMatrix"));
          int index = std::distance(origopsvec.begin(), ft);
          T(index,kk)+=st->second;}}
    uint tmin=gin.getUInt("MinimumTimeSeparation");
    uint tmax=gin.getUInt("MaximumTimeSeparation");
    bool herm=false; gin.getOptionalBool("HermitianMatrix",herm);
    bool subvev=false; gin.getOptionalBool("SubtractVEV",subvev); 
    logger.putEcho(gin,"OtherInput");
    string datamode="none";
    gin.getOptionalString("Mode",datamode);
    if ((datamode!="bins")&&(datamode!="samplings")&&(datamode!="none"))
       throw(std::invalid_argument("Invalid Sampling Mode in TransformCorrelatorMatrix"));
    string ftype(datamode),filename;
    WriteMode wmode = Protect;  // protect mode
    string fformat("default"); char ffmt='D';
    bool writetofile=false;
    bool vsep=false;
    if (gin.queryTag("WriteToFile")){
       writetofile=true;
       ArgsHandler ggin(gin,"WriteToFile");
       filename=ggin.getName("FileName");
       ggin.getOptionalString("FileType",ftype);
       if ((ftype!="bins")&&(ftype!="samplings"))
          throw(std::invalid_argument("<FileType> must be bins or samplings in TransformCorrelatorMatrix"));
       if ((ftype!=datamode)&&(datamode!="none"))
          throw(std::invalid_argument("<Mode> and <FileType> must match in TransformCorrelatorMatrix"));
       ggin.getOptionalString("FileFormat",fformat);
       if (fformat=="fstr") ffmt='F';
       else if (fformat=="hdf5") ffmt='H';
       else if (fformat=="default") ffmt='D';
       else throw(std::invalid_argument("<FileFormat> must be ftr or hdf5 or default in TransformCorrelatorMatrix"));
       string fmode="protect"; ggin.getOptionalString("WriteMode",fmode);
       if (fmode=="overwrite") wmode=Overwrite;
       else if (fmode=="update") wmode=Update;
       if (subvev && ftype=="samplings"){
          ggin.getOptionalBool("SeparateVEVWrite",vsep);}
       logger.putEcho(ggin);}
    if (ftype=="none")
       throw(std::invalid_argument("At least one of <Mode> and <FileType> must appear in TransformCorrelatorMatrix"));
    XMLHandler xmlo; logger.output(xmlo);
    xmlout.put_child(xmlo);
    set<MCObsInfo> obskeys;
    bool eraseorig=true;
    if (ftype=="bins"){
       doTransformedCorrMatrixByBins(*m_obs,origopsvec,transopsvec,herm,subvev,tmin,tmax,
                                      T,eraseorig,obskeys);}
    else if (ftype=="samplings"){
       doTransformedCorrMatrixBySamplings(*m_obs,origopsvec,transopsvec,herm,subvev,tmin,tmax,
                                           T,eraseorig,obskeys,vsep);}
    xmlout.put_child("NumberOfRealObservablesProcessed",make_string(obskeys.size()));
    if (writetofile){
       XMLHandler xmlf;
       if (ftype=="bins")
          m_obs->writeBinsToFile(obskeys,filename,xmlf,wmode,ffmt);
       else
          m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,wmode,ffmt);
       xmlout.put_child(xmlf);}
    xmlout.put_child("Status","Done");}
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type TransformCorrelatorMatrix encountered an error: ")
             +string(errmsg.what())));} }


 else if (functype=="Exp"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","Exp");
    try{
    XMLHandler xmlinobs(xmltask,"InObservable");
    XMLHandler xmlt1,xmlt2;
    MCObsInfo obsin(xmlinobs);
    xmlt1.set_root("InObservable");
    obsin.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    string datamode="samplings";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="samplings"){
       if (m_obs->isJackknifeMode()){
          mcode='J'; datamode="Jackknife";}
       else{
          mcode='B'; datamode="Bootstrap";}}
    else throw(std::invalid_argument("Invalid Sampling Mode"));
    xmlout.put_child("Mode",datamode);

    XMLHandler xmlres(xmltask,"Result");
    string name; int index;
    xmlreadchild(xmlres,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for Exp result"));
    index=taskcount;
    xmlreadifchild(xmlres,"IDIndex",index);
    MCObsInfo resinfo(name,index,mcode=='D');
    xmlt1.set_root("ResultInfo");
    resinfo.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    if (mcode=='D'){
       doExpByBins(*m_obs,obsin,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       doExpBySamplings(*m_obs,obsin,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type Exp encountered an error: ")
                +string(errmsg.what())));}
    }

 else if (functype=="Log"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","Log");
    try{
    XMLHandler xmlinobs(xmltask,"InObservable");
    XMLHandler xmlt1,xmlt2;
    MCObsInfo obsin(xmlinobs);
    xmlt1.set_root("InObservable");
    obsin.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    string datamode="samplings";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="samplings"){
       if (m_obs->isJackknifeMode()){
          mcode='J'; datamode="Jackknife";}
       else{
          mcode='B'; datamode="Bootstrap";}}
    else throw(std::invalid_argument("Invalid Sampling Mode"));
    xmlout.put_child("Mode",datamode);

    XMLHandler xmlres(xmltask,"Result");
    string name; int index;
    xmlreadchild(xmlres,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for Log result"));
    index=taskcount;
    xmlreadifchild(xmlres,"IDIndex",index);
    MCObsInfo resinfo(name,index,mcode=='D');
    xmlt1.set_root("ResultInfo");
    resinfo.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    if (mcode=='D'){
       doLogByBins(*m_obs,obsin,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       doLogBySamplings(*m_obs,obsin,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);
       m_obs->setSamplingMode(origmode);} }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("DoObsFunction with type Log encountered an error: ")
                +string(errmsg.what())));}
    }


 else if (functype=="Copy"){
    xmlout.set_root("DoObsFunction");
    xmlout.put_child("Type","Copy");
    try{
    XMLHandler xmlinobs(xmltask,"InObservable");
    XMLHandler xmlt1,xmlt2;
    MCObsInfo obsin(xmlinobs);
    xmlt1.set_root("InObservable");
    obsin.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    string datamode="samplings";
    xmlreadifchild(xmltask,"Mode",datamode);
    char mcode;
    if (datamode=="bins") mcode='D';
    else if (datamode=="Bootstrap") mcode='B';
    else if (datamode=="Jackknife") mcode='J';
    else if (datamode=="samplings"){
       if (m_obs->isJackknifeMode()){
          mcode='J'; datamode="Jackknife";}
       else{
          mcode='B'; datamode="Bootstrap";}}
    else throw(std::invalid_argument("Invalid Sampling Mode"));
    xmlout.put_child("Mode",datamode);

    XMLHandler xmlres(xmltask,"Result");
    string name; int index;
    xmlreadchild(xmlres,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for Copy result"));
    index=taskcount;
    xmlreadifchild(xmlres,"IDIndex",index);
    MCObsInfo resinfo(name,index,mcode=='D');
    xmlt1.set_root("ResultInfo");
    resinfo.output(xmlt2);
    xmlt1.put_child(xmlt2);
    xmlout.put_child(xmlt1);

    if (mcode=='D'){
       doCopyByBins(*m_obs,obsin,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
       est.output(xmlt1);
       xmlout.put_child(xmlt1);}
    else{
       SamplingMode origmode=m_obs->getCurrentSamplingMode();
       if (mcode=='J') m_obs->setToJackknifeMode();
       else m_obs->setToBootstrapMode();
       doCopyBySamplings(*m_obs,obsin,resinfo);
       MCEstimate est=m_obs->getEstimate(resinfo);
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
 
