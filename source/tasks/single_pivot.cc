#include "single_pivot.h"

using namespace std;
using namespace LaphEnv;

 // ******************************************************************


SinglePivotOfCorrMat::SinglePivotOfCorrMat(TaskHandler& taskhandler, ArgsHandler& xml_in,
                                           LogHelper& xmlout)
                        : m_moh(taskhandler.getMCObsHandler()), m_cormat_info(0), 
                          m_rotated_info(0), m_Zmat(0), m_transmat(0), m_rotcorset(0)
{
 try{
    ArgsHandler xmlin(xml_in,"SinglePivotInitiate");
    if (xmlin.queryTag("ReadPivotFromFile")){
       ArgsHandler xmlf(xmlin,"ReadPivotFromFile");
       initiate_from_file(xmlf,xmlout);}
    else
       initiate_new(xmlin,xmlout);}
 catch(const std::exception& errmsg){
    clear();
    throw(std::invalid_argument((string("Constructing SinglePivotOfCorrMat failed: ")
          +string(errmsg.what())).c_str()));}
}



void SinglePivotOfCorrMat::initiate_new(ArgsHandler& xmlin, LogHelper& xmlout)
{
 xmlout.reset("InitiateNew");
 try{
    string rotate_name(xmlin.getName("RotatedCorrelatorName"));
    m_rotated_info=new OperatorInfo(rotate_name,OperatorInfo::Rotated);
    m_cormat_info=new CorrelatorMatrixInfo(
           xmlin.getItem<CorrelatorMatrixInfo>("CorrelatorMatrixInfo"));
    if (!m_cormat_info->isHermitian()){
       throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for rotation"));}
    if (m_cormat_info->getNumberOfOperators()<2){
       throw(std::invalid_argument("CorrelatorMatrix must have at least 2 operators in RotatedCorrelatorMatrix"));}
    xmlin.getUInt("NormTime",m_tauN);
    xmlin.getUInt("MetricTime",m_tau0);
    xmlin.getUInt("DiagonalizeTime",m_tauD);
    xmlin.getReal("MinimumInverseConditionNumber",m_min_inv_condnum);
    m_neg_eig_alarm=-5.0*m_min_inv_condnum;
    xmlin.getOptionalReal("NegativeEigenvalueAlarm",m_neg_eig_alarm);
    bool checkMetricErrors=false, checkCommonNullSpace=false;
    xmlin.getOptionalBool("CheckMetricErrors",checkMetricErrors);
    xmlin.getOptionalBool("CheckCommonMetricMatrixNullSpace",checkCommonNullSpace);
    if ((m_min_inv_condnum<0.0)||(m_tau0==m_tauD)||(m_tauN>m_tau0)||(m_tauN>m_tauD))
       throw(std::invalid_argument("Invalid parameters in SinglePivot"));
    xmlout.putEcho(xmlin);
    LogHelper xmlc;
    create_pivot(xmlc,checkMetricErrors,checkCommonNullSpace);
    xmlout.put(xmlc);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument((string("SinglePivot failed in doCorrMatrixRotation: ")
           +string(errmsg.what())).c_str()));}
 if (xmlin.queryTag("WritePivotToFile")){
    ArgsHandler xmlf(xmlin,"WritePivotToFile");
    string fname(xmlf.getName("PivotFileName"));
    bool overwrite=xmlf.getBool("Overwrite");
    write_to_file(fname,overwrite);
    xmlout.putEcho(xmlf);}
}


void SinglePivotOfCorrMat::initiate_from_file(ArgsHandler& xmlin, LogHelper& xmlout)
{
 string fname(xmlin.getName("PivotFileName"));
 IOMap<UIntKey,ArrayBuf> iom;
 string filetypeid("Sigmond--SinglePivotFile");
 string header;
 iom.openReadOnly(fname,filetypeid,header);
 if (!iom.isOpen())
     throw(std::invalid_argument("File could not be opened for input"));
 XMLHandler xmlh; xmlh.set_from_string(header);
 ArgsHandler xmlr(xmlh,"SinglePivotOfCorrMat");
 string rotate_name(xmlr.getName("RotatedCorrelatorName"));
 m_rotated_info=new OperatorInfo(rotate_name,OperatorInfo::Rotated);
 m_cormat_info=new CorrelatorMatrixInfo(
           xmlr.getItem<CorrelatorMatrixInfo>("CorrelatorMatrixInfo"));
 xmlr.getUInt("NormTime",m_tauN);
 xmlr.getUInt("MetricTime",m_tau0);
 xmlr.getUInt("DiagonalizeTime",m_tauD);
 xmlr.getReal("MinimumInverseConditionNumber",m_min_inv_condnum);
 xmlout.putEcho(xmlr);
 ArrayBuf buffer;
 iom.get(UIntKey(0),buffer);
 TransMatrix *matptr=new TransMatrix;
 array_to_matrix(buffer,*matptr);
 m_transmat=matptr;  // now pointer to const
 iom.get(UIntKey(1),buffer);
 TransMatrix *zptr=new TransMatrix;
 array_to_matrix(buffer,*zptr);
 m_Zmat=zptr;  // now pointer to const
 iom.close();
 uint nlevels=m_transmat->size(1);
 m_rotcorset=new DiagonalCorrelatorSet;
 for (uint opindex=0;opindex<nlevels;opindex++){
    m_rotcorset->addCorrelator(m_rotated_info->resetRotatedLevel(opindex));}
}



void SinglePivotOfCorrMat::write_to_file(const string& filename, bool overwrite)
{
 string fname=tidyString(filename);
 if (fname.empty()){
    throw(std::invalid_argument("Error in SinglePivotWriteToFile:: Empty file name"));}
 if ((fileExists(fname))&&(!overwrite)){
    throw(std::invalid_argument("Error in SingePivotWriteToFile:: File exists and cannot overwrite"));}
 XMLHandler xmlout("SinglePivotOfCorrMat");
 XMLHandler xmlt; m_cormat_info->output(xmlt,false);
 xmlout.put_child(xmlt);
 xmlout.seek_first_child();
 xmlout.put_sibling("RotatedCorrelatorName",m_rotated_info->getRotatedName());
 xmlout.put_sibling("NormTime",make_string(m_tauN));
 xmlout.put_sibling("MetricTime",make_string(m_tau0));
 xmlout.put_sibling("DiagonalizeTime",make_string(m_tauD));
 xmlout.put_sibling("MinimumInverseConditionNumber",make_string(m_min_inv_condnum));

 IOMap<UIntKey,ArrayBuf> iom;
 string filetypeid("Sigmond--SinglePivotFile");   
 iom.openNew(fname,filetypeid,xmlout.str(),false,'N',false,overwrite);
 if (!iom.isOpen())
     throw(std::invalid_argument("File could not be opened for output"));
 ArrayBuf buffer;
 matrix_to_array(*m_transmat,buffer);
 iom.put(UIntKey(0),buffer);
 matrix_to_array(*m_Zmat,buffer);
 iom.put(UIntKey(1),buffer);
 iom.close();
}



uint SinglePivotOfCorrMat::getNumberOfOperators() const
{
 return (m_cormat_info ? m_cormat_info->getNumberOfOperators() : 0);
}

uint SinglePivotOfCorrMat::getNumberOfLevels() const
{
 return (m_transmat ? m_transmat->size(1) : 0);
}


OperatorInfo SinglePivotOfCorrMat::getRotatedOperator() const
{
 if (m_rotated_info==0)
    throw(std::runtime_error("RotatedOperator not yet set in SinglePivotOfCorrMat"));
 return *m_rotated_info;
}


bool SinglePivotOfCorrMat::isVEVsubtracted() const
{
 return (m_cormat_info ? m_cormat_info->isVEVSubtracted() : false);
}


void SinglePivotOfCorrMat::create_pivot(LogHelper& xmlout, bool checkMetricErrors,
                                        bool checkCommonNullSpace)
{
 xmlout.reset("CreatePivot");
 HermMatrix corrN,corr0,corrD;
 VVector vev;
 bool subvev=m_cormat_info->isVEVSubtracted();
 m_moh->setSamplingBegin();   // rotate using full estimates
 vector<MCEstimate> corr0_diag,corrD_diag;
 bool save_memory=false;//true;
 try{
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,*m_cormat_info,m_tauN,corrN);
    if (subvev) getHermCorrelatorMatrixVEVs_CurrentSampling(m_moh,*m_cormat_info,vev);
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,*m_cormat_info,m_tau0,corr0);
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,*m_cormat_info,m_tauD,corrD);
    getDiagonalCorrelatorsAtTimeEstimates(m_moh,*m_cormat_info,m_tau0,corr0_diag);
    getDiagonalCorrelatorsAtTimeEstimates(m_moh,*m_cormat_info,m_tauD,corrD_diag);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument((string("get Correlator matrix failed in SinglePivot: ")
          +string(errmsg.what())).c_str()));}

    // output fractional errors in diagonal elements of C(tau0), C(tauD)
    // for informational purposes
 uint opnum=0;
 LogHelper xmlcd("DiagonalCorrelatorFractionalErrors");
 const set<OperatorInfo>& corrops=m_cormat_info->getOperators();
 for (set<OperatorInfo>::const_iterator itop=corrops.begin();itop!=corrops.end();itop++,opnum++){
    LogHelper xmlv("DiagonalCorrelator");
    XMLHandler xmlop;
    itop->output(xmlop,false);  // short form
    xmlv.put(xmlop);
    const MCEstimate& mc0=corr0_diag[opnum];
    double fracerr=mc0.getSymmetricError()/std::abs(mc0.getFullEstimate());
    xmlv.putReal("MetricTimeFractionalError",fracerr);
    const MCEstimate& mcD=corrD_diag[opnum];
    fracerr=mcD.getSymmetricError()/std::abs(mcD.getFullEstimate());
    xmlv.putReal("MatrixTimeFractionalError",fracerr);
    xmlcd.put(xmlv);}
 xmlout.put(xmlcd);

         //  jackknife the largest and smallest eigenvalues of 
         //  renormalized metric
 
 if (checkMetricErrors){
 try{
 SamplingMode currmode=m_moh->getCurrentSamplingMode();
 m_moh->setToJackknifeMode();
 uint nops=m_cormat_info->getNumberOfOperators();
 Diagonalizer DG;
 HermMatrix corrjN,corrj0;
 RVector lambda;
 vector<MCObsInfo> tempkeys;
 for (uint k=0;k<nops;k++)
    tempkeys.push_back(MCObsInfo("TempMetricEigevalue",k));
 for (m_moh->begin();!m_moh->end();++(*m_moh)){
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,*m_cormat_info,m_tauN,corrjN);
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,*m_cormat_info,m_tau0,corrj0);
    doRescaleByDiagonals(corrj0,corrjN);
    DG.getEigenvalues(corrj0,lambda);
    for (uint k=0;k<nops;k++)
       m_moh->putCurrentSamplingValue(tempkeys[k],lambda[k],true);}
 LogHelper xmllambda("RescaledMetricConditionVariation");
 for (uint k=0;k<nops;k++){
    MCEstimate lambdaest=m_moh->getEstimate(tempkeys[k]);
    m_moh->eraseData(tempkeys[k]);
    LogHelper xmlmetriceig("MetricEigenvalue");
    xmlmetriceig.putUInt("Level",k);
    XMLHandler xmltemp;
    lambdaest.output(xmltemp);
    xmlmetriceig.put(xmltemp);
    xmllambda.put(xmlmetriceig);}
 xmlout.put(xmllambda);
 m_moh->setSamplingMode(currmode);
 m_moh->setSamplingBegin();}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument((string("get Correlator matrix failed in SinglePivot: ")
          +string(errmsg.what())).c_str()));}}

 if (save_memory){
 try{
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_cormat_info,m_tauN);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_cormat_info,m_tau0);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_cormat_info,m_tauD);
    if (subvev) eraseHermCorrelatorMatrixVEVs(m_moh,*m_cormat_info);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument((string("get Correlator matrix failed in SinglePivot: ")
          +string(errmsg.what())).c_str()));}}


      // rescale corr0 and corrD using corrN......
 doRescaleByDiagonals(corr0,corrN);
 doRescaleByDiagonals(corrD,corrN);

      // set the metric
 DiagonalizerWithMetric DM(m_min_inv_condnum,m_neg_eig_alarm);
 DM.setExceptionsOff();
 LogHelper logmetric;
 int info=DM.setMetric(corr0,logmetric);
 xmlout.putItem(logmetric);
 if (info!=0) 
    throw(std::invalid_argument((string("setMetric encountered problem in SinglePivot: ")
                 +DiagonalizerWithMetric::getRotateMetricCode(info)+string("Log: \n\n")
                 +xmlout.output()).c_str()));

     // set the matrix
 LogHelper logmatrix;
 info=DM.setMatrix(corrD,logmatrix,checkCommonNullSpace); 
 xmlout.putItem(logmatrix);
 if ((info!=0)&&(info!=-5)) 
    throw(std::invalid_argument((string("setMatrix encountered problem in SinglePivot: ")
        +DiagonalizerWithMetric::getRotateMatrixCode(info)+string("Log: \n\n")
        +xmlout.output()).c_str()));

         //  set the rotation matrix and the Zmatrix
         //  (remember to include the rescaling)
 TransMatrix rotationMatrix,Zmat;
 DM.getEigenvectors(rotationMatrix);
 doRescaleTransformation(rotationMatrix,corrN);
 DM.getZMatrix(Zmat);
 doRescaleTransformation(Zmat,corrN);
 uint nlevels=rotationMatrix.size(1);

         // if there are nonzero VEVs, rephase rotated operators
         // so all VEVs are real and positive
 if (subvev){
    uint nops=rotationMatrix.size(0);
    doVectorRotation(vev,rotationMatrix);
    for (uint col=0;col<nlevels;col++){
#if defined COMPLEXNUMBERS
       complex<double> phase(vev[col]/std::abs(vev[col]));
#else
       double phase=(vev[col]>=0.0)?1.0:-1.0; 
#endif
       for (uint row=0;row<nops;row++){
          rotationMatrix(row,col)*=phase;
          Zmat(row,col)*=phase;}}}

 m_transmat=new TransMatrix(rotationMatrix);
 m_Zmat=new TransMatrix(Zmat);

 m_rotcorset=new DiagonalCorrelatorSet;
 for (uint opindex=0;opindex<nlevels;opindex++){
    m_rotcorset->addCorrelator(m_rotated_info->resetRotatedLevel(opindex));}
}


void SinglePivotOfCorrMat::clear()
{
 delete m_cormat_info;
 delete m_rotated_info;
 delete m_Zmat;
 delete m_transmat;
 delete m_rotcorset;
 m_cormat_info=0; 
 m_rotated_info=0;
 m_Zmat=0;
 m_transmat=0;
 m_rotcorset=0;
}


SinglePivotOfCorrMat::~SinglePivotOfCorrMat()
{
 clear();
}


SinglePivotOfCorrMat* SinglePivotOfCorrMat::initiateFromMemory(
                          TaskHandler& taskhandler, 
                          ArgsHandler& xml_in, LogHelper& xmlout)
{
 xmlout.reset("InitiateFromMemory");
 ArgsHandler xmlin(xml_in,"SinglePivotInitiate");
 if (!xmlin.queryTag("GetFromMemory"))
    return 0;
 try{
    ArgsHandler xmln(xmlin,"GetFromMemory");
    string idname(xmln.getName("IDName"));
    TaskHandlerData* ptr=taskhandler.get_task_data(idname);
    if (ptr){
       xmlout.putEcho(xmlin);
       SinglePivotOfCorrMat *pptr=dynamic_cast<SinglePivotOfCorrMat*>(ptr);
       return pptr;}
    else{
       xmlout.putString("Error","Id name of SinglePivotOfCorrMat not in memory");
       throw(std::invalid_argument("Id name of SinglePivotOfCorrMat not in memory"));}}
 catch(const std::exception& msg){
    throw;}
 return 0;
}



bool SinglePivotOfCorrMat::putInMemory(TaskHandler& taskhandler, ArgsHandler& xml_in,
                                       LogHelper& xmlout, SinglePivotOfCorrMat* pivot)
{
 xmlout.reset("PutIntoMemory");
 if (!pivot) return false;
 try{
    ArgsHandler xmlin(xml_in,"SinglePivotInitiate");
    string idname(xmlin.getName("AssignName"));
    taskhandler.insert_task_data(idname,pivot);
    xmlout.putString("AssignedName",idname);
    return true;}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument((string("Error trying to putInMemory SinglePivotOfCorrMat: ")
         +string(errmsg.what())).c_str()));}
}



SinglePivotOfCorrMat* SinglePivotOfCorrMat::initiateSinglePivot(
                   TaskHandler& taskhandler, ArgsHandler& xmlin,
                   LogHelper& xmlout, bool& keep_in_task_map)
{
 ArgsHandler xmlpiv(xmlin,"SinglePivotInitiate");
 LogHelper xmlt;
 SinglePivotOfCorrMat* pivoter=0;
 try{
    pivoter=SinglePivotOfCorrMat::initiateFromMemory(taskhandler,xmlpiv,xmlt);
    if (pivoter){
       keep_in_task_map=true;  // was already in memory
       xmlout.putEcho(xmlpiv);
       return pivoter;}}
 catch(const std::exception& errmsg){
    xmlout.putItem(xmlt); 
    throw(std::invalid_argument((string("Error in SinglePivotOfCorrMat::initiateFromMemory: ")
           +string(errmsg.what())).c_str()));}
 keep_in_task_map=true;
 try{
    pivoter=new SinglePivotOfCorrMat(taskhandler,xmlpiv,xmlt);
    xmlout.putEcho(xmlpiv);}
 catch(const std::exception& errmsg){
    xmlout.putItem(xmlt); 
    throw(std::invalid_argument((string("Error in SinglePivotOfCorrMat::initiating new: ")
          +string(errmsg.what())).c_str()));}
 LogHelper xmlp;
 keep_in_task_map=SinglePivotOfCorrMat::putInMemory(taskhandler,xmlpiv,xmlp,pivoter);
 if (keep_in_task_map){ xmlt.putItem(xmlp); xmlout.putItem(xmlt);}
 return pivoter;
}




void SinglePivotOfCorrMat::doRotation(uint tmin, uint tmax, LogHelper& xmllog)
{
 xmllog.reset("DoRotation");
 bool vevs=m_cormat_info->isVEVSubtracted();
 bool flag=true;
 if (vevs){
    try{
       do_vev_rotation();
       xmllog.putString("VEVRotation","Success");}
    catch(const std::exception& errmsg){
       flag=false;
       xmllog.putString("VEVRotation",string("Failure: ")+string(errmsg.what()));
       throw(std::invalid_argument((string("VEVRotation failed ")
               +string(errmsg.what())).c_str()));}}
 for (uint tval=tmin;tval<=tmax;tval++){
    bool diagonly=(tval<=m_tauD) ? true : false;
    LogHelper xmlc("CorrelatorRotation");
    xmlc.putUInt("TimeValue",tval);
    try{
       do_corr_rotation(tval,diagonly);
       xmlc.putString("Status","Success");}
    catch(const std::exception& errmsg){
       flag=false;
       xmlc.putString(string("Status"),string("Failure: ")+string(errmsg.what()));}
    if (!diagonly){
           // perform some tests of off-diagonal correlators for t>tauD
           // these should be zero, so check this; remove from memory
      uint nlevels=getNumberOfLevels();
      uint onesigma=0, twosigma=0, threesigma=0;
      uint foursigma=0, large=0; double maxviolation=0.0;
      OperatorInfo rowop(*m_rotated_info);
      OperatorInfo colop(*m_rotated_info);
      for (uint col=0;col<nlevels;col++){
         colop.resetRotatedLevel(col);
         for (uint row=0;row<col;row++){
            rowop.resetRotatedLevel(row);
            MCObsInfo obskey(rowop,colop,tval,true,RealPart,vevs); 
            MCEstimate est_re=m_moh->getEstimate(obskey,Bootstrap);
            m_moh->eraseData(obskey);
            double offratio=std::abs(est_re.getFullEstimate())/est_re.getSymmetricError();
            if (offratio<=1.0) onesigma++;
            else if (offratio<=2.0) twosigma++;
            else if (offratio<=3.0) threesigma++;
            else if (offratio<=4.0) foursigma++;
            else large++;
            if (offratio>maxviolation) maxviolation=offratio;
#ifdef COMPLEXNUMBERS
            obskey.setToImaginaryPart();
            MCEstimate est_im=m_moh->getEstimate(obskey,Bootstrap);
            m_moh->eraseData(obskey);
            offratio=std::abs(est_im.getFullEstimate())/est_im.getSymmetricError();
            if (offratio<=1.0) onesigma++;
            else if (offratio<=2.0) twosigma++;
            else if (offratio<=3.0) threesigma++;
            else if (offratio<=4.0) foursigma++;
            else large++;
            if (offratio>maxviolation) maxviolation=offratio;
#endif
            }}
      LogHelper xmloff("OffDiagonalChecks");
      LogHelper xmlck("MaximumDeviationFromZero");
      xmlck.putReal("RelativeToError",maxviolation);
      xmloff.put(xmlck);
      xmlck.reset("PercentDeviationsFromZero");
      uint ntotal=onesigma+twosigma+threesigma+foursigma+large;
      xmlck.putReal("GreaterThanOneSigma",
           100.0*double(twosigma+threesigma+foursigma+large)/double(ntotal));
      xmlck.putReal("GreaterThanTwoSigma",
           100.0*double(threesigma+foursigma+large)/double(ntotal));
      xmlck.putReal("GreaterThanThreeSigma",
           100.0*double(foursigma+large)/double(ntotal));
      xmlck.putReal("GreaterThanFourSigma",
           100.0*double(large)/double(ntotal));
      xmloff.put(xmlck);
      xmlc.put(xmloff);}
    xmllog.put(xmlc);}
 if (!flag) xmllog.putString("Warning","some rotations failed");
}


    //  Computes the VEVs of the rotated operators and puts
    //  them into memory.  The VEVs of the original operators are
    //  read from file and put into memory, but afterwards, the
    //  memory used by them is freed.

#ifdef COMPLEXNUMBERS


void SinglePivotOfCorrMat::do_vev_rotation()
{
 uint nops=getNumberOfOperators();
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_cormat_info->getOperators();
                // read original bins, arrange pointers in certain way
 vector<const Vector<double>* > binptrs(2*nops);  // pointers to original bins
 uint count=0;
 for (set<OperatorInfo>::const_iterator it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    binptrs[count++]=&(m_moh->getBins(vev_info));
    vev_info.setToImaginaryPart();
    binptrs[count++]=&(m_moh->getBins(vev_info));}

 vector<Vector<double> > VEVrotated(nlevels,nbins);  // only real parts
 CVector VEVbuffer(nops);

      // loop over bins
 for (uint bin=0;bin<nbins;bin++){
            // read this one bin into a vector
    count=0;
    for (uint k=0;k<nops;k++){
       double br=(*binptrs[count++])[bin];
       double bi=(*binptrs[count++])[bin];
       VEVbuffer[k]=complex<double>(br,bi);}
              // do the rotation
    doVectorRotation(VEVbuffer,*m_transmat);
              // store results in VEVrotated
    count=0;
    for (uint level=0;level<nlevels;level++){
       VEVrotated[count++][bin]=VEVbuffer[level].real();}
    VEVbuffer.resize(nops);}

             //  free up memory no longer needed (original bins)

 for (set<OperatorInfo>::const_iterator it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    m_moh->eraseData(vev_info);
    vev_info.setToImaginaryPart();
    m_moh->eraseData(vev_info);}

       // put rotated bins into memory (imaginary parts should be zero)

 count=0;
 for (uint level=0;level<nlevels;level++){
    m_rotated_info->resetRotatedLevel(level);
    MCObsInfo rotvev_info(*m_rotated_info,RealPart);
    m_moh->putBins(rotvev_info,VEVrotated[count++]);}

}



void SinglePivotOfCorrMat::do_corr_rotation(uint timeval, bool diagonly)
{
 uint nops=getNumberOfOperators();
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_cormat_info->getOperators();
                // read original bins, arrange pointers in certain way
 vector<const Vector<double>* > binptrs(nops*nops);  // pointers to original bins
 uint count=0;
 set<OperatorInfo>::const_iterator itrow,itcol;
 try{
 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,true,false);
       MCObsInfo obskey(corrt,RealPart);
       binptrs[count++]=&(m_moh->getBins(obskey));
       obskey.setToImaginaryPart();
       binptrs[count++]=&(m_moh->getBins(obskey));}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,true,false);
    binptrs[count++]=&(m_moh->getBins(MCObsInfo(corrt,RealPart)));}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument("Could not read bins in do_corr_rotation"));}

 uint ncompute=(diagonly ? nlevels : nlevels*nlevels);
 vector<Vector<double> > Crotated(ncompute,nbins);
 ComplexHermitianMatrix Cbuffer(nops);

      // loop over bins
 for (uint bin=0;bin<nbins;bin++){
            // read this one bin into a matrix
    count=0;
    for (uint col=0;col<nops;col++){
       for (uint row=0;row<col;row++){
          double br=(*binptrs[count++])[bin];
          double bi=(*binptrs[count++])[bin];
          Cbuffer.put(row,col,complex<double>(br,bi));}
       double br=(*binptrs[count++])[bin];
       Cbuffer.put(col,col,complex<double>(br,0.0));}
              // do the rotation
    if (diagonly){
       RVector diagbuf;
       doMatrixRotation(Cbuffer,*m_transmat,diagbuf);
                // store results in Crotated
       for (uint level=0;level<nlevels;level++)
          Crotated[level][bin]=diagbuf[level];}
    else{
       doMatrixRotation(Cbuffer,*m_transmat);
                // store results in Crotated
       count=0;
       for (uint col=0;col<nlevels;col++){
       for (uint row=0;row<col;row++){
             Crotated[count++][bin]=Cbuffer(row,col).real();
             Crotated[count++][bin]=Cbuffer(row,col).imag();}
          Crotated[count++][bin]=Cbuffer(col,col).real();}
       Cbuffer.resize(nops);}}

       //  free up memory for original bins

 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,true,false);
       MCObsInfo obskey(corrt,RealPart);
       m_moh->eraseData(obskey);
       obskey.setToImaginaryPart();
       m_moh->eraseData(obskey);}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,true,false);
    m_moh->eraseData(MCObsInfo(corrt,RealPart));}

       // put rotated bins into memory
 if (diagonly){
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetRotatedLevel(level);
       MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,true,RealPart,false);
       m_moh->putBins(obskey,Crotated[level]);}}
 else{
    OperatorInfo rowop(*m_rotated_info);
    OperatorInfo colop(*m_rotated_info);
    count=0;
    for (uint col=0;col<nlevels;col++){
       colop.resetRotatedLevel(col);
       for (uint row=0;row<col;row++){
          rowop.resetRotatedLevel(row);
          MCObsInfo obskey(rowop,colop,timeval,true,RealPart,false);
          m_moh->putBins(obskey,Crotated[count++]);
          obskey.setToImaginaryPart();
          m_moh->putBins(obskey,Crotated[count++]);}
       MCObsInfo obskey(colop,colop,timeval,true,RealPart,false);
       m_moh->putBins(obskey,Crotated[count++]);}}

}


#elif defined REALNUMBERS


void SinglePivotOfCorrMat::do_vev_rotation()
{
 uint nops=getNumberOfOperators();
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_cormat_info->getOperators();
                // read original bins, arrange pointers in certain way
 vector<const Vector<double>* > binptrs(nops);  // pointers to original bins
 uint count=0;
 for (set<OperatorInfo>::const_iterator it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    binptrs[count++]=&(m_moh->getBins(vev_info));}

 vector<Vector<double> > VEVrotated(nlevels,nbins);
 RVector VEVbuffer(nops);

      // loop over bins
 for (uint bin=0;bin<nbins;bin++){
            // read this one bin into a vector
    count=0;
    for (uint k=0;k<nops;k++){
       VEVbuffer[k]=(*binptrs[count++])[bin];}
              // do the rotation
    doVectorRotation(VEVbuffer,*m_transmat);
              // store results in VEVrotated
    count=0;
    for (uint level=0;level<nlevels;level++){
       VEVrotated[count++][bin]=VEVbuffer[level];}
    VEVbuffer.resize(nops);}

             //  free up memory no longer needed (original bins)

 for (set<OperatorInfo>::const_iterator it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    m_moh->eraseData(vev_info);}

       // put rotated bins into memory

 count=0;
 for (uint level=0;level<nlevels;level++){
    m_rotated_info->resetRotatedLevel(level);
    MCObsInfo rotvev_info(*m_rotated_info,RealPart);
    m_moh->putBins(rotvev_info,VEVrotated[count++]);}

}



void SinglePivotOfCorrMat::do_corr_rotation(uint timeval, bool diagonly)
{
 uint nops=getNumberOfOperators();
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_cormat_info->getOperators();
                // read original bins, arrange pointers in certain way
 vector<const Vector<double>* > binptrs((nops*(nops+1))/2);  // pointers to original bins
 uint count=0;
 set<OperatorInfo>::const_iterator itrow,itcol;
 try{
 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,true,false);
       MCObsInfo obskey(corrt,RealPart);
       binptrs[count++]=&(m_moh->getBins(obskey));}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,true,false);
    binptrs[count++]=&(m_moh->getBins(MCObsInfo(corrt,RealPart)));}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument("Could not read bins in do_corr_rotation"));}

 uint ncompute=(diagonly ? nlevels : (nlevels*(nlevels+1))/2);
 vector<Vector<double> > Crotated(ncompute,nbins);
 RealSymmetricMatrix Cbuffer(nops);

      // loop over bins
 for (uint bin=0;bin<nbins;bin++){
            // read this one bin into a matrix
    count=0;
    for (uint col=0;col<nops;col++)
       for (uint row=0;row<=col;row++){
          Cbuffer(row,col)=(*binptrs[count++])[bin];}
              // do the rotation
    if (diagonly){
       RVector diagbuf;
       doMatrixRotation(Cbuffer,*m_transmat,diagbuf);
                // store results in Crotated
       for (uint level=0;level<nlevels;level++)
          Crotated[level][bin]=diagbuf[level];}
    else{
       doMatrixRotation(Cbuffer,*m_transmat);
                // store results in Crotated
       count=0;
       for (uint col=0;col<nlevels;col++)
       for (uint row=0;row<=col;row++){
             Crotated[count++][bin]=Cbuffer(row,col);}
       Cbuffer.resize(nops);}}

       //  free up memory for original bins

 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,true,false);
       MCObsInfo obskey(corrt,RealPart);
       m_moh->eraseData(obskey);}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,true,false);
    m_moh->eraseData(MCObsInfo(corrt,RealPart));}

       // put rotated bins into memory
 if (diagonly){
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetRotatedLevel(level);
       MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,true,RealPart,false);
       m_moh->putBins(obskey,Crotated[level]);}}
 else{
    OperatorInfo rowop(*m_rotated_info);
    OperatorInfo colop(*m_rotated_info);
    count=0;
    for (uint col=0;col<nlevels;col++){
       colop.resetRotatedLevel(col);
       for (uint row=0;row<=col;row++){
          rowop.resetRotatedLevel(row);
          MCObsInfo obskey(rowop,colop,timeval,true,RealPart,false);
          m_moh->putBins(obskey,Crotated[count++]);}}}

}

#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif



void SinglePivotOfCorrMat::writeRotated(uint tmin, uint tmax, const string& corrfile,
                                        bool overwrite, LogHelper& xmlout)
{
 xmlout.reset("WriteRotated");
 uint nlevels=getNumberOfLevels();
 set<MCObsInfo> obskeys;
 bool vevs=m_cormat_info->isVEVSubtracted();
 if (vevs){
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetRotatedLevel(level);
       MCObsInfo obskey(*m_rotated_info,RealPart);
       obskeys.insert(obskey);}}
 for (uint level=0;level<nlevels;level++){
    m_rotated_info->resetRotatedLevel(level);
    for (uint timeval=tmin;timeval<=tmax;timeval++){
       MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,true,RealPart,false);
       obskeys.insert(obskey);}}
 XMLHandler xmlf;
 m_moh->writeBinsToFile(obskeys,corrfile,xmlf,overwrite);
 xmlout.put(xmlf);
}


 // ******************************************************************


  //  get |Z(opindex,level)|^2 for all operators for all levels

void SinglePivotOfCorrMat::computeZMagnitudesSquared(Matrix<MCEstimate>& ZMagSq)
{
 if (!m_rotcorset->allFitInfoAvailable())
    throw(std::runtime_error("Not all fit info available to compute ZMagSquares"));

 uint nops=getNumberOfOperators();
 uint nlevels=getNumberOfLevels();

   //  Method:
   //  each diagonal element of rotated correlator matrix is
   //  fit to  A * exp(-E_n t )    then   Zrotsq = std::abs(A)

 ZMagSq.resize(nops,nlevels);   //  final results put in here
 bool overwrite=true;
 MCObsInfo obskey("TempZMagSq",0);   // temporary key
 for (uint level=0;level<nlevels;level++){
    MCObsInfo Zrotfitkey(m_rotcorset->getAmplitudeKey(level));
    for (m_moh->begin();!m_moh->end();++(*m_moh)){   // loop over resamplings
       double Zrotsq=std::abs(m_moh->getCurrentSamplingValue(Zrotfitkey));
       for (uint opindex=0;opindex<nops;opindex++){
          obskey.resetObsIndex(opindex);
          m_moh->putCurrentSamplingValue(obskey,
                  sqr((*m_Zmat)(opindex,level))*Zrotsq,overwrite);}}
    for (uint opindex=0;opindex<nops;opindex++){
       obskey.resetObsIndex(opindex);
       ZMagSq(opindex,level)=m_moh->getEstimate(obskey);
       m_moh->eraseSamplings(obskey); }}

}

 // ******************************************************************