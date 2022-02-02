#include "single_pivot.h"

using namespace std;
using namespace LaphEnv;

 // ******************************************************************


SinglePivotOfCorrMat::SinglePivotOfCorrMat(TaskHandler& taskhandler, ArgsHandler& xml_in,
                                           LogHelper& xmlout)
                        : m_moh(taskhandler.getMCObsHandler()), m_cormat_info(0),
                          m_orig_cormat_info(0), m_rotated_info(0), m_Zmat(0),
                          m_transmat(0), m_imp_trans(0), m_vevs_avail(true)
{
 try{
    ArgsHandler xmlin(xml_in,"SinglePivotInitiate");
    if (xmlin.queryTag("ReadPivotFromFile")){
       ArgsHandler xmlf(xmlin,"ReadPivotFromFile");
       initiate_from_file(xmlf,xmlout);}
    else{
       initiate_new(xmlin,xmlout);}}
 catch(const std::exception& errmsg){
    clear();
    throw(std::invalid_argument(string("Constructing SinglePivotOfCorrMat failed: ")
          +string(errmsg.what())));}
}



void SinglePivotOfCorrMat::initiate_new(ArgsHandler& xmlin, LogHelper& xmlout)
{
 xmlout.reset("InitiateNew");
 try{
    ArgsHandler xmlg(xmlin,"RotatedCorrelator");
    m_rotated_info=new GenIrrepOperatorInfo(
           xmlg.getItem<GenIrrepOperatorInfo>("RotatedCorrelator"));
    m_cormat_info=new CorrelatorMatrixInfo(
           xmlin.getItem<CorrelatorMatrixInfo>("CorrelatorMatrixInfo"));
    if (!m_cormat_info->isHermitian()){
       throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for rotation"));}
    if (m_cormat_info->getNumberOfOperators()<2){
       throw(std::invalid_argument("CorrelatorMatrixInfo must have at least 2 operators"));}
    if (xmlin.queryTag("ImprovedOperators")){
       map<OperatorInfo,map<OperatorInfo,Scalar> > trans;
       read_trans(xmlin,trans);
       setup_improved_operators(trans);}
    else{
       m_orig_cormat_info=m_cormat_info;}
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
    bool set_imag_zero=false;
    xmlin.getOptionalBool("SetImaginaryPartsZero",set_imag_zero);
    list<CorrelatorInfo> set_to_zero;
    xmlin.getMultiItems("SetToZero",set_to_zero);
    xmlout.putEcho(xmlin);
    LogHelper xmlc;
    create_pivot(xmlc,checkMetricErrors,checkCommonNullSpace,set_to_zero,set_imag_zero);
    xmlout.put(xmlc);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("SinglePivot failed in doCorrMatrixRotation: ")
           +string(errmsg.what())));}
 if (xmlin.queryTag("WritePivotToFile")){
    ArgsHandler xmlf(xmlin,"WritePivotToFile");
    string fname(xmlf.getName("PivotFileName"));
    bool overwrite=xmlf.getBool("Overwrite");
    string fformat("default"); char ffmt='D';
    xmlf.getOptionalString("FileFormat",fformat);
    if (fformat=="fstr") ffmt='F';
    else if (fformat=="hdf5") ffmt='H';
    else if (fformat=="default") ffmt='D';
    else throw(std::invalid_argument("<FileFormat> must be ftr or hdf5 or default in doCorrMatrixRotation"));
    write_to_file(fname,overwrite,ffmt);
    xmlout.putEcho(xmlf);}
 if (xmlin.queryTag("PrintTransformationMatrix")){
    LogHelper xmltrans;
    print_trans(xmltrans);
    xmlout.put(xmltrans);}
}


void SinglePivotOfCorrMat::initiate_from_file(ArgsHandler& xmlin, LogHelper& xmlout)
{
 xmlout.reset("InitiatedFromFile");
 string fname(xmlin.getName("PivotFileName"));
 IOMap<UIntKey,ArrayBuf> iom;
#ifdef COMPLEXNUMBERS
 string filetypeid("Sigmond--SinglePivotFile-CN");
#else
 string filetypeid("Sigmond--SinglePivotFile-RN");
#endif
 string header;
 iom.openReadOnly(fname,filetypeid,header);
 if (!iom.isOpen())
     throw(std::invalid_argument("File could not be opened for input"));
 XMLHandler xmlh; xmlh.set_from_string(header);
 ArgsHandler xmlr(xmlh,"SinglePivotOfCorrMat");
 ArgsHandler xmlc(xmlr,"RotatedCorrelator");
 m_rotated_info=new GenIrrepOperatorInfo(
           xmlc.getItem<GenIrrepOperatorInfo>("RotatedCorrelator"));
 ArgsHandler xmlmd(xmlr,"MatrixDefinition");
 m_cormat_info=new CorrelatorMatrixInfo(
           xmlmd.getItem<CorrelatorMatrixInfo>("CorrelatorMatrixInfo"));
 if (xmlr.queryTag("OriginalMatrix")){
    ArgsHandler xmlo(xmlr,"OriginalMatrix");
    m_orig_cormat_info=new CorrelatorMatrixInfo(
           xmlo.getItem<CorrelatorMatrixInfo>("CorrelatorMatrixInfo"));}
 else{
    m_orig_cormat_info=m_cormat_info;}
 xmlr.getUInt("NormTime",m_tauN);
 xmlr.getUInt("MetricTime",m_tau0);
 xmlr.getUInt("DiagonalizeTime",m_tauD);
 xmlr.getReal("MinimumInverseConditionNumber",m_min_inv_condnum);
 xmlout.putString("FileName",fname);
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
}



void SinglePivotOfCorrMat::write_to_file(const string& filename, bool overwrite, char file_format)
{
 string fname=tidyString(filename);
 if (fname.empty()){
    throw(std::invalid_argument("Error in SinglePivotWriteToFile:: Empty file name"));}
 if ((fileExists(fname))&&(!overwrite)){
    throw(std::invalid_argument("Error in SingePivotWriteToFile:: File exists and cannot overwrite"));}
 XMLHandler xmlout("SinglePivotOfCorrMat");
 XMLHandler xmlt; m_cormat_info->output(xmlt,false);
 XMLHandler xmlmatdef("MatrixDefinition");
 xmlmatdef.put_child(xmlt);
 xmlout.put_child(xmlmatdef);
 if (m_cormat_info!=m_orig_cormat_info){
    m_orig_cormat_info->output(xmlt,false);
    XMLHandler xmlorig("OriginalMatrix");
    xmlorig.put_child(xmlt);
    xmlout.put_child(xmlorig);}
 XMLHandler xmlrot; m_rotated_info->output(xmlrot,false);
 xmlt.set_root("RotatedCorrelator"); xmlt.put_child(xmlrot);
 xmlout.put_child(xmlt);
 xmlout.seek_first_child();
 xmlout.put_sibling("NormTime",make_string(m_tauN));
 xmlout.put_sibling("MetricTime",make_string(m_tau0));
 xmlout.put_sibling("DiagonalizeTime",make_string(m_tauD));
 xmlout.put_sibling("MinimumInverseConditionNumber",make_string(m_min_inv_condnum));

 IOMap<UIntKey,ArrayBuf> iom;
#ifdef COMPLEXNUMBERS
 string filetypeid("Sigmond--SinglePivotFile-CN");
#else
 string filetypeid("Sigmond--SinglePivotFile-RN");
#endif
 iom.openNew(fname,filetypeid,xmlout.str(),false,'N',false,overwrite,file_format);
 if (!iom.isOpen())
     throw(std::invalid_argument("File could not be opened for output"));
 ArrayBuf buffer;
 matrix_to_array(*m_transmat,buffer);
 iom.put(UIntKey(0),buffer);
 matrix_to_array(*m_Zmat,buffer);
 iom.put(UIntKey(1),buffer);
 iom.close();
}



void SinglePivotOfCorrMat::print_trans(LogHelper& xmlout)
{
 const set<OperatorInfo>& origops=m_orig_cormat_info->getOperators();
 xmlout.reset("TransformationMatrix");
 LogHelper xmli("ImprovedOperators");
 for (uint level=0;level<m_transmat->size(1);level++){
    LogHelper xmllevel("ImprovedOperator");
    LogHelper xmlname("OpName");
    xmlname.putItem(m_rotated_info->resetIDIndex(level));
    xmllevel.put(xmlname);
    uint opindex=0;
    for (set<OperatorInfo>::const_iterator it=origops.begin();it!=origops.end();it++,opindex++){
       LogHelper term("OpTerm");
       term.putItem(*it);
       term.putString("Coefficient",make_string((*m_transmat)(opindex,level)));
       xmllevel.put(term);}
    xmli.put(xmllevel);}
 xmlout.put(xmli);
}


void SinglePivotOfCorrMat::read_trans(ArgsHandler& xmlin,
                    map<OperatorInfo,map<OperatorInfo,Scalar> >& trans)
{
 trans.clear();
 ArgsHandler xmli(xmlin,"ImprovedOperators");
 list<ArgsHandler> xmlo(xmli.getSubHandlers("ImprovedOperator"));
 for (list<ArgsHandler>::iterator it=xmlo.begin();it!=xmlo.end();it++){
    ArgsHandler xmlopname(*it,"OpName");
    OperatorInfo opimp;
    xmlopname.getItem("OpName",opimp);
    if (trans.find(opimp)!=trans.end())
       throw(std::invalid_argument("Repeated improved operator in SinglePivot"));
    it->insert(xmlopname);
    map<OperatorInfo,Scalar> opdef;
    list<ArgsHandler> xmlc(it->getSubHandlers("OpTerm"));
    for (list<ArgsHandler>::iterator ct=xmlc.begin();ct!=xmlc.end();ct++){
       OperatorInfo opterm;
       ct->getItem("OpTerm",opterm);
       if (opdef.find(opimp)!=opdef.end())
          throw(std::invalid_argument("Repeated operator term in ImprovedOperator in SinglePivot"));
       Scalar cf;
       ct->getScalar("Coefficient",cf);
       opdef.insert(make_pair(opterm,cf));
      it->insert(*ct);}
    trans.insert(make_pair(opimp,opdef));
    xmli.insert(*it);}
 xmlin.insert(xmli);
}


void SinglePivotOfCorrMat::setup_improved_operators(const std::map<OperatorInfo,
                std::map<OperatorInfo,Scalar> >& trans)
{
 set<OperatorInfo> addset=m_cormat_info->getOperators();
    //  from improved operators, get original set of operators
 set<OperatorInfo> origset;
 for (map<OperatorInfo,map<OperatorInfo,Scalar> >::const_iterator
     it=trans.begin();it!=trans.end();it++){
    if (addset.erase(it->first)==0)     // remove improved operators from "addset"
       throw(std::invalid_argument("Improved operator not in CorrelatorMatrixInfo"));
    for (map<OperatorInfo,Scalar>::const_iterator
       ct=(it->second).begin();ct!=(it->second).end();ct++)
          origset.insert(ct->first);}
   // operators remaining in "addset" must be original operators (not improved)
 origset.insert(addset.begin(),addset.end());

 m_orig_cormat_info=new CorrelatorMatrixInfo(origset,true,m_cormat_info->subtractVEV());
 map<OperatorInfo,uint> origind;
 uint count=0;
 for (set<OperatorInfo>::iterator ot=origset.begin();ot!=origset.end();ot++,count++)
    origind.insert(make_pair(*ot,count));
 uint norig=origset.size();
 uint nops=getNumberOfOperators();
 TransMatrix tt(norig,nops);
 count=0;
 const set<OperatorInfo>& opset=m_cormat_info->getOperators();
 for (set<OperatorInfo>::const_iterator it=opset.begin();it!=opset.end();it++,count++){
    map<OperatorInfo,map<OperatorInfo,Scalar> >::const_iterator kt=trans.find(*it);
    if (kt==trans.end()){
       map<OperatorInfo,uint>::const_iterator ut=origind.find(*it);
       if (ut==origind.end())
          throw(std::invalid_argument("Something went wrong with indexing in setup_improved_operators"));
       uint oind=ut->second;
       tt(oind,count)=1.0;}
    else{
       for (map<OperatorInfo,Scalar>::const_iterator
            ct=(kt->second).begin();ct!=(kt->second).end();ct++){
          Scalar cf=ct->second;
          const OperatorInfo& opt=ct->first;
          map<OperatorInfo,uint>::const_iterator ut=origind.find(opt);
          if (ut==origind.end())
             throw(std::invalid_argument("Something went wrong with indexing in setup_improved_operators"));
          uint oind=ut->second;
          tt(oind,count)=cf;}}}
 m_imp_trans=new TransMatrix(tt);
}




uint SinglePivotOfCorrMat::getNumberOfOperators() const
{
 return (m_cormat_info ? m_cormat_info->getNumberOfOperators() : 0);
}

uint SinglePivotOfCorrMat::getNumberOfLevels() const
{
 return (m_transmat ? m_transmat->size(1) : 0);
}


const std::set<OperatorInfo>& SinglePivotOfCorrMat::getOperators() const
{
 if (m_cormat_info==0)
    throw(std::runtime_error("CorrelatorMatrixInfo not yet set in SinglePivotOfCorrMat"));
 return m_cormat_info->getOperators();
}


const std::set<OperatorInfo>& SinglePivotOfCorrMat::getOriginalOperators() const
{
 if (m_orig_cormat_info==0)
    throw(std::runtime_error("CorrelatorMatrixInfo not yet set in SinglePivotOfCorrMat"));
 return m_orig_cormat_info->getOperators();
}


GenIrrepOperatorInfo SinglePivotOfCorrMat::getRotatedOperator() const
{
 if (m_rotated_info==0)
    throw(std::runtime_error("RotatedOperator not yet set in SinglePivotOfCorrMat"));
 return *m_rotated_info;
}


bool SinglePivotOfCorrMat::subtractVEV() const
{
 return (m_cormat_info ? m_cormat_info->subtractVEV() : false);
}


void SinglePivotOfCorrMat::create_pivot(LogHelper& xmlout, bool checkMetricErrors,
                                        bool checkCommonNullSpace, 
                                        const std::list<CorrelatorInfo>& set_to_zero,
                                        bool setImagPartsZero)
{
 xmlout.reset("CreatePivot");
 if (m_moh->isJackknifeMode()) xmlout.putString("ResamplingMode","Jackknife");
 else xmlout.putString("ResamplingMode","Bootstrap");
 HermMatrix corrN,corr0,corrD;
 VVector vev;
 bool subvev=m_cormat_info->subtractVEV();
 m_moh->setSamplingBegin();   // rotate using full estimates
 vector<MCEstimate> corr0_diag,corrD_diag;
 bool save_memory=false;
      // if some elements have been requested to be set to zero, determine their rows,columns
 list<pair<uint,uint> > elems_set_to_zero(getCorrelatorMatrixIndices(*m_cormat_info,set_to_zero,true));
 try{
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,m_cormat_info,m_tauN,corrN,m_orig_cormat_info,m_imp_trans);
    if (subvev){
       try{
       getHermCorrelatorMatrixVEVs_CurrentSampling(m_moh,m_cormat_info,vev,m_orig_cormat_info,m_imp_trans);}
       catch(const std::exception& xp){
          m_vevs_avail=false;}}
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,m_cormat_info,m_tau0,corr0,m_orig_cormat_info,m_imp_trans);
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,m_cormat_info,m_tauD,corrD,m_orig_cormat_info,m_imp_trans);
    if (m_orig_cormat_info==m_cormat_info){
       getDiagonalCorrelatorsAtTimeEstimates(m_moh,*m_cormat_info,m_tau0,corr0_diag);
       getDiagonalCorrelatorsAtTimeEstimates(m_moh,*m_cormat_info,m_tauD,corrD_diag);}
    if (setImagPartsZero){
       setImaginaryPartsToZero(corrN);
       setImaginaryPartsToZero(corr0);
       setImaginaryPartsToZero(corrD);}
    for (list<pair<uint,uint> >::iterator zt=elems_set_to_zero.begin();zt!=elems_set_to_zero.end();++zt){
       corrN.setToZero(zt->first,zt->second);
       corr0.setToZero(zt->first,zt->second);
       corrD.setToZero(zt->first,zt->second);}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("get Correlator matrix failed in SinglePivot: ")
          +string(errmsg.what())));}
    // output fractional errors in diagonal elements of C(tau0), C(tauD)
    // for informational purposes
 if (m_orig_cormat_info==m_cormat_info){
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
    xmlout.put(xmlcd);}

         //  resample the largest and smallest eigenvalues of
         //  renormalized metric

 if (checkMetricErrors){
 try{
 uint nops=m_cormat_info->getNumberOfOperators();
 Diagonalizer DG;
 HermMatrix corrjN,corrj0;
 RVector lambda;
 vector<MCObsInfo> tempkeys;
 for (uint k=0;k<nops;k++)
    tempkeys.push_back(MCObsInfo("TempMetricEigevalue",k));
 for (m_moh->begin();!m_moh->end();++(*m_moh)){
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,m_cormat_info,m_tauN,corrjN,m_orig_cormat_info,m_imp_trans);
    getHermCorrelatorMatrixAtTime_CurrentSampling(m_moh,m_cormat_info,m_tau0,corrj0,m_orig_cormat_info,m_imp_trans);
    if (setImagPartsZero){
       setImaginaryPartsToZero(corrj0);
       setImaginaryPartsToZero(corrjN);}
    for (list<pair<uint,uint> >::iterator zt=elems_set_to_zero.begin();zt!=elems_set_to_zero.end();++zt){
       corrjN.setToZero(zt->first,zt->second);
       corrj0.setToZero(zt->first,zt->second);}
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
 m_moh->setSamplingBegin();}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Check Metric Errors failed in SinglePivot: ")
          +string(errmsg.what())));}}

 if (save_memory){
 try{
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_cormat_info,m_tauN);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_cormat_info,m_tau0);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_cormat_info,m_tauD);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_orig_cormat_info,m_tauN);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_orig_cormat_info,m_tau0);
    eraseHermCorrelatorMatrixAtTime(m_moh,*m_orig_cormat_info,m_tauD);
    if (subvev){
       eraseHermCorrelatorMatrixVEVs(m_moh,*m_cormat_info);
       eraseHermCorrelatorMatrixVEVs(m_moh,*m_orig_cormat_info);}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("save memory failed in SinglePivot: ")
          +string(errmsg.what())));}}

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
    throw(std::invalid_argument(string("setMetric encountered problem in SinglePivot: ")
                 +DiagonalizerWithMetric::getRotateMetricCode(info)+string("Log: \n\n")
                 +xmlout.output()));

     // set the matrix
 LogHelper logmatrix;
 info=DM.setMatrix(corrD,logmatrix,checkCommonNullSpace);
 xmlout.putItem(logmatrix);
 if ((info!=0)&&(info!=-5))
    throw(std::invalid_argument(string("setMatrix encountered problem in SinglePivot: ")
        +DiagonalizerWithMetric::getRotateMatrixCode(info)+string("Log: \n\n")
        +xmlout.output()));

         //  set the rotation matrix and the Zmatrix
         //  (remember to include the rescaling)
 TransMatrix rotationMatrix,Zmat;
 DM.getEigenvectors(rotationMatrix);
 doRescaleTransformation(rotationMatrix,corrN);
 DM.getZMatrix(Zmat);
 uint nlevels=rotationMatrix.size(1);

         // if there are nonzero VEVs, rephase rotated operators
         // so all VEVs are real and positive
 if (subvev && m_vevs_avail){
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

 if (m_imp_trans!=0){
    doMatrixMultiply(*m_imp_trans,rotationMatrix);}
 m_transmat=new TransMatrix(rotationMatrix);
 m_Zmat=new TransMatrix(Zmat);
}


void SinglePivotOfCorrMat::clear()
{
 if (m_cormat_info==m_orig_cormat_info){
    delete m_cormat_info;}
 else{
    delete m_cormat_info;
    delete m_orig_cormat_info;}
 m_cormat_info=0;
 m_orig_cormat_info=0;
 delete m_rotated_info;
 delete m_Zmat;
 delete m_transmat;
 delete m_imp_trans;
 m_rotated_info=0;
 m_Zmat=0;
 m_transmat=0;
 m_imp_trans=0;
 m_ampkeys.clear();
 m_energykeys.clear();
 m_reorder.clear();
 m_vevs_avail=true;
}


SinglePivotOfCorrMat::~SinglePivotOfCorrMat()
{
 clear();
}


SinglePivotOfCorrMat* SinglePivotOfCorrMat::initiateFromMemory(
                          TaskHandler& taskhandler,
                          ArgsHandler& xml_in, LogHelper& xmlout)
{
 ArgsHandler xmlin(xml_in,"SinglePivotInitiate");
 if (!xmlin.queryTag("GetFromMemory"))
    return 0;
 try{
    ArgsHandler xmln(xmlin,"GetFromMemory");
    string idname(xmln.getName("IDName"));
    TaskHandlerData* ptr=taskhandler.get_task_data(idname);
    if (ptr){
       xmlout.reset("InitiateFromMemory");
       xmlout.putString("IDName",idname);
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
    throw(std::invalid_argument(string("Error trying to putInMemory SinglePivotOfCorrMat: ")
         +string(errmsg.what())));}
}



SinglePivotOfCorrMat* SinglePivotOfCorrMat::initiateSinglePivot(
                   TaskHandler& taskhandler, ArgsHandler& xmlin,
                   LogHelper& xmlout, bool& keep_in_task_map)
{
 ArgsHandler xmlpiv(xmlin,"SinglePivotInitiate");
 LogHelper xmlt;
 xmlout.reset("SinglePivot");
 SinglePivotOfCorrMat* pivoter=0;
 try{
    pivoter=SinglePivotOfCorrMat::initiateFromMemory(taskhandler,xmlpiv,xmlt);
    if (pivoter){
       xmlout.put(xmlt);
       keep_in_task_map=true;
       xmlout.putEcho(xmlpiv);
       return pivoter;}}
 catch(const std::exception& errmsg){
    xmlout.putItem(xmlt);
    throw(std::invalid_argument(string("Error in SinglePivotOfCorrMat::initiateFromMemory: ")
           +string(errmsg.what())));}
 keep_in_task_map=false;
 try{
    pivoter=new SinglePivotOfCorrMat(taskhandler,xmlpiv,xmlt);
    xmlout.put(xmlt);}
 catch(const std::exception& errmsg){
  xmlout.put(xmlt);
    throw(std::invalid_argument(string("Error in SinglePivotOfCorrMat::initiating: ")
          +string(errmsg.what())));}
 if (xmlpiv.queryTag("AssignName")){
    LogHelper xmlp;
    keep_in_task_map=SinglePivotOfCorrMat::putInMemory(taskhandler,xmlpiv,xmlp,pivoter);
    if (keep_in_task_map){ xmlout.put(xmlp);}}

 return pivoter;
}


   // mode 'B' => rotate by bins
   //      'S' => rotate by samplings and only the vev-subtracted correlators
   //      'U' => rotate by samplings the unsubtracted correlators and vevs
   //      'A' => all rotate by samplings the vevs and subtracted and unsubtracted correlators
 
void SinglePivotOfCorrMat::doRotation(uint tmin, uint tmax, char mode, LogHelper& xmllog)
{
 xmllog.reset("DoRotation");
 bool vevs=m_cormat_info->subtractVEV();
 bool flag=true;
 bool herm=true;
 if (vevs && m_vevs_avail && (mode!='S')){
    try{
    if (mode=='B')
       do_vev_rotation_by_bins();
    else if ((mode=='U')||(mode=='A'))
       do_vev_rotation_by_samplings();
    xmllog.putString("VEVRotation","Success");}
    catch(const std::exception& errmsg){
       flag=false;
       xmllog.putString("VEVRotation",string("Failure: ")+string(errmsg.what()));
       throw(std::invalid_argument(string("VEVRotation failed ")
               +string(errmsg.what())));}}
 for (uint tval=tmin;tval<=tmax;tval++){
    bool diagonly=(tval<=m_tauD) ? true : false;
    LogHelper xmlc("CorrelatorRotation");
    xmlc.putUInt("TimeValue",tval);
    try{
       if (mode=='B')
          do_corr_rotation_by_bins(tval,diagonly);
       else if ((mode=='S')||(!vevs))
          do_corr_rotation_by_samplings(tval,diagonly,mode);
       else if ((mode=='A')||(mode=='U')){
          do_corr_rotation_by_samplings(tval,diagonly,'S');
          do_corr_rotation_by_samplings(tval,diagonly,'U');}
       else
         throw(std::invalid_argument(string("Invalid mode in Pivot doRotation")));
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
      GenIrrepOperatorInfo rowop(*m_rotated_info);
      GenIrrepOperatorInfo colop(*m_rotated_info);
      for (uint col=0;col<nlevels;col++){
         colop.resetIDIndex(col);
         for (uint row=0;row<col;row++){
            rowop.resetIDIndex(row);
            MCObsInfo obskey(OperatorInfo(rowop),OperatorInfo(colop),tval,herm,RealPart,vevs);
            MCEstimate est_re=m_moh->getEstimate(obskey);
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
            MCEstimate est_im=m_moh->getEstimate(obskey);
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
      if (m_moh->isJackknifeMode()) xmloff.putString("ResamplingMode","Jackknife");
      else xmloff.putString("ResamplingMode","Bootstrap");
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


void SinglePivotOfCorrMat::do_vev_rotation_by_bins()
{
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();

                // read original bins, arrange pointers in certain way
 vector<const Vector<double>* > binptrs(2*nops);  // pointers to original bins
 uint count=0;
 for (set<OperatorInfo>::const_iterator it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    binptrs[count++]=&(m_moh->getBins(vev_info));
    vev_info.setToImaginaryPart();
    binptrs[count++]=&(m_moh->getBins(vev_info));}

 vector<Vector<double> > VEVrotated(nlevels,nbins);  // only real parts
 CVector VEVbuffer;

      // loop over bins
 for (uint bin=0;bin<nbins;bin++){
    VEVbuffer.resize(nops);
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
       VEVrotated[count++][bin]=VEVbuffer[level].real();}}

             //  free up memory no longer needed (original bins)

 for (set<OperatorInfo>::const_iterator it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    m_moh->eraseData(vev_info);
    vev_info.setToImaginaryPart();
    m_moh->eraseData(vev_info);}

       // put rotated bins into memory (imaginary parts should be zero)

 Vector<double> imVEVrotated(nbins,0.0);  // only real parts
 count=0;
 for (uint level=0;level<nlevels;level++){
    m_rotated_info->resetIDIndex(level);
    MCObsInfo rotvev_info(*m_rotated_info,RealPart);
    m_moh->putBins(rotvev_info,VEVrotated[count++]);
    rotvev_info.setToImaginaryPart();
    m_moh->putBins(rotvev_info,imVEVrotated);}

}


void SinglePivotOfCorrMat::do_vev_rotation_by_samplings()
{
 uint nlevels=getNumberOfLevels();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
 set<OperatorInfo>::const_iterator it;
 CVector Vbuffer;
 
      // loop over samplings
 for (m_moh->begin(); !m_moh->end(); m_moh->setSamplingNext()){
    Vbuffer.resize(nops);
            // read original VEVs for this sampling
    try{
       uint count=0;
       for (it=ops.begin();it!=ops.end();++it){
          MCObsInfo vev_info(*it,RealPart);
          double vr=m_moh->getCurrentSamplingValue(vev_info);
          vev_info.setToImaginaryPart();
          double vi=m_moh->getCurrentSamplingValue(vev_info);
          Vbuffer.put(count++,complex<double>(vr,vi));}}
    catch(const std::exception& errmsg){
       throw(std::invalid_argument("Could not obtain samplings in do_vev_rotation"));}
              // do the rotation
    doVectorRotation(Vbuffer,*m_transmat);
       // put rotated sampling into memory (imaginary parts should be zero)
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetIDIndex(level);
       MCObsInfo obskey(*m_rotated_info,RealPart);
       m_moh->putCurrentSamplingValue(obskey,Vbuffer[level].real());
       obskey.setToImaginaryPart();
       m_moh->putCurrentSamplingValue(obskey,0.0);}}

             //  free up memory no longer needed (original vevs)
 for (it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    m_moh->eraseData(vev_info);
    vev_info.setToImaginaryPart();
    m_moh->eraseData(vev_info);}
}


void SinglePivotOfCorrMat::do_corr_rotation_by_bins(uint timeval, bool diagonly)
{
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
 bool herm=true;
                // read original bins, arrange pointers in certain way
 vector<const Vector<double>* > binptrs(nops*nops);  // pointers to original bins
 uint count=0;
 set<OperatorInfo>::const_iterator itrow,itcol;
 try{
 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,herm,false);
       MCObsInfo obskey(corrt,RealPart);
       binptrs[count++]=&(m_moh->getBins(obskey));
       obskey.setToImaginaryPart();
       binptrs[count++]=&(m_moh->getBins(obskey));}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,herm,false);
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
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,herm,false);
       MCObsInfo obskey(corrt,RealPart);
       m_moh->eraseData(obskey);
       obskey.setToImaginaryPart();
       m_moh->eraseData(obskey);}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,herm,false);
    MCObsInfo obskey(corrt,RealPart);
    m_moh->eraseData(obskey);
    obskey.setToImaginaryPart();
    m_moh->eraseData(obskey);}

       // put rotated bins into memory
 if (diagonly){
    Vector<double> imCdiag(nbins,0.0);
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetIDIndex(level);
       MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,herm,RealPart,false);
       m_moh->putBins(obskey,Crotated[level]);
       obskey.setToImaginaryPart();
       m_moh->putBins(obskey,imCdiag);}}
 else{
    Vector<double> imCdiag(nbins,0.0);
    GenIrrepOperatorInfo rowop(*m_rotated_info);
    GenIrrepOperatorInfo colop(*m_rotated_info);
    count=0;
    for (uint col=0;col<nlevels;col++){
       colop.resetIDIndex(col);
       for (uint row=0;row<col;row++){
          rowop.resetIDIndex(row);
          MCObsInfo obskey(rowop,colop,timeval,herm,RealPart,false);
          m_moh->putBins(obskey,Crotated[count++]);
          obskey.setToImaginaryPart();
          m_moh->putBins(obskey,Crotated[count++]);}
       MCObsInfo obskey(colop,colop,timeval,herm,RealPart,false);
       m_moh->putBins(obskey,Crotated[count++]);
       obskey.setToImaginaryPart();
       m_moh->putBins(obskey,imCdiag);}}

}

void SinglePivotOfCorrMat::do_corr_rotation_by_samplings(uint timeval, bool diagonly, char mode)
{
 uint nlevels=getNumberOfLevels();
 bool subvev=(mode=='U') ? false : m_cormat_info->subtractVEV();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
 bool herm=true;

 set<OperatorInfo>::const_iterator itrow,itcol;
 ComplexHermitianMatrix Cbuffer;

      // loop over samplings
 for (m_moh->begin(); !m_moh->end(); m_moh->setSamplingNext()){
    Cbuffer.resize(nops);
            // read this one sampling into a matrix
    try{
    uint col=0;
    for (itcol=ops.begin();itcol!=ops.end();itcol++,col++){
      uint row=0;
      for (itrow=ops.begin();itrow!=itcol;itrow++,row++){
         CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,herm,subvev);
         MCObsInfo obskey(corrt,RealPart);
         double br=m_moh->getCurrentSamplingValue(obskey);
         obskey.setToImaginaryPart();
         double bi=m_moh->getCurrentSamplingValue(obskey);
         Cbuffer.put(row,col,complex<double>(br,bi));}
      CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,herm,subvev);
      double br=m_moh->getCurrentSamplingValue(MCObsInfo(corrt,RealPart));
      Cbuffer.put(col,col,complex<double>(br,0.0));}}
    catch(const std::exception& errmsg){
       throw(std::invalid_argument("Could not obtain samplings in do_corr_rotation"));}

              // do the rotation
    if (diagonly){
       RVector diagbuf;
       doMatrixRotation(Cbuffer,*m_transmat,diagbuf);
       for (uint level=0;level<nlevels;level++){
          m_rotated_info->resetIDIndex(level);
          MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,herm,RealPart,subvev);
          m_moh->putCurrentSamplingValue(obskey,diagbuf[level]);
          obskey.setToImaginaryPart();
          m_moh->putCurrentSamplingValue(obskey,0.0);}}
    else{
       doMatrixRotation(Cbuffer,*m_transmat);
       GenIrrepOperatorInfo rowop(*m_rotated_info);
       GenIrrepOperatorInfo colop(*m_rotated_info);
       for (uint col=0;col<nlevels;col++){
           colop.resetIDIndex(col);
           for (uint row=0;row<col;row++){
              rowop.resetIDIndex(row);
              MCObsInfo obskey(rowop,colop,timeval,herm,RealPart,subvev);
              m_moh->putCurrentSamplingValue(obskey,Cbuffer(row,col).real());
              obskey.setToImaginaryPart();
              m_moh->putCurrentSamplingValue(obskey,Cbuffer(row,col).imag());}
           MCObsInfo obskey(colop,colop,timeval,herm,RealPart,subvev);
           m_moh->putCurrentSamplingValue(obskey,Cbuffer(col,col).real());
           obskey.setToImaginaryPart();
           m_moh->putCurrentSamplingValue(obskey,0.0);}}}

       //  free up memory for original samplings

 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,herm,subvev);
       MCObsInfo obskey(corrt,RealPart);
       m_moh->eraseData(obskey);
       obskey.setToImaginaryPart();
       m_moh->eraseData(obskey);}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,herm,subvev);
    MCObsInfo obskey(corrt,RealPart);
    obskey.setToImaginaryPart();
    m_moh->eraseData(obskey);}
}


#elif defined REALNUMBERS


void SinglePivotOfCorrMat::do_vev_rotation_by_bins()
{
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
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
    m_rotated_info->resetIDIndex(level);
    MCObsInfo rotvev_info(*m_rotated_info,RealPart);
    m_moh->putBins(rotvev_info,VEVrotated[count++]);}

}


void SinglePivotOfCorrMat::do_vev_rotation_by_samplings()
{
 uint nlevels=getNumberOfLevels();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
 set<OperatorInfo>::const_iterator it;
 RVector Vbuffer;
 
      // loop over samplings
 for (m_moh->begin(); !m_moh->end(); m_moh->setSamplingNext()){
    Vbuffer.resize(nops);
            // read original VEVs for this sampling
    try{
       uint count=0;
       for (it=ops.begin();it!=ops.end();++it){
          MCObsInfo vev_info(*it,RealPart);
          double vr=m_moh->getCurrentSamplingValue(vev_info);
          Vbuffer[count++]=vr;}}
    catch(const std::exception& errmsg){
       throw(std::invalid_argument("Could not obtain samplings in do_vev_rotation"));}
              // do the rotation
    doVectorRotation(Vbuffer,*m_transmat);
       // put rotated sampling into memory
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetIDIndex(level);
       MCObsInfo obskey(*m_rotated_info,RealPart);
       m_moh->putCurrentSamplingValue(obskey,Vbuffer[level]);}}

             //  free up memory no longer needed (original vevs)
 for (it=ops.begin();it!=ops.end();it++){
    MCObsInfo vev_info(*it,RealPart);
    m_moh->eraseData(vev_info);}
}


void SinglePivotOfCorrMat::do_corr_rotation_by_bins(uint timeval, bool diagonly)
{
 uint nlevels=getNumberOfLevels();
 uint nbins=m_moh->getNumberOfBins();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
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
       m_rotated_info->resetIDIndex(level);
       MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,true,RealPart,false);
       m_moh->putBins(obskey,Crotated[level]);}}
 else{
    GenIrrepOperatorInfo rowop(*m_rotated_info);
    GenIrrepOperatorInfo colop(*m_rotated_info);
    count=0;
    for (uint col=0;col<nlevels;col++){
       colop.resetIDIndex(col);
       for (uint row=0;row<=col;row++){
          rowop.resetIDIndex(row);
          MCObsInfo obskey(OperatorInfo(rowop),OperatorInfo(colop),timeval,true,RealPart,false);
          m_moh->putBins(obskey,Crotated[count++]);}}}

}

void SinglePivotOfCorrMat::do_corr_rotation_by_samplings(uint timeval, bool diagonly, char mode)
{
 uint nlevels=getNumberOfLevels();
 bool subvev=(mode=='U') ? false : m_cormat_info->subtractVEV();
 const set<OperatorInfo>& ops=m_orig_cormat_info->getOperators();
 uint nops=ops.size();
 bool herm=true;

 set<OperatorInfo>::const_iterator itrow,itcol;
 RealSymmetricMatrix Rbuffer;

      // loop over samplings
 for (m_moh->begin(); !m_moh->end(); m_moh->setSamplingNext()){
    Rbuffer.resize(nops);
            // read this one sampling into a matrix
    try{
    uint col=0;
    for (itcol=ops.begin();itcol!=ops.end();itcol++,col++){
      uint row=0;
      for (itrow=ops.begin();itrow!=itcol;itrow++,row++){
         CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,herm,subvev);
         MCObsInfo obskey(corrt,RealPart);
         Rbuffer(row,col)=m_moh->getCurrentSamplingValue(obskey);}
      CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,herm,subvev);
      Rbuffer(col,col)=m_moh->getCurrentSamplingValue(MCObsInfo(corrt,RealPart));}}
    catch(const std::exception& errmsg){
       throw(std::invalid_argument("Could not obtain samplings in do_corr_rotation"));}

              // do the rotation
    if (diagonly){
       RVector diagbuf;
       doMatrixRotation(Rbuffer,*m_transmat,diagbuf);
       for (uint level=0;level<nlevels;level++){
          m_rotated_info->resetIDIndex(level);
          MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,herm,RealPart,subvev);
          m_moh->putCurrentSamplingValue(obskey,diagbuf[level]);}}
    else{
       doMatrixRotation(Rbuffer,*m_transmat);
       GenIrrepOperatorInfo rowop(*m_rotated_info);
       GenIrrepOperatorInfo colop(*m_rotated_info);
       for (uint col=0;col<nlevels;col++){
           colop.resetIDIndex(col);
           for (uint row=0;row<=col;row++){
              rowop.resetIDIndex(row);
              MCObsInfo obskey(OperatorInfo(rowop),OperatorInfo(colop),timeval,herm,RealPart,subvev);
              m_moh->putCurrentSamplingValue(obskey,Rbuffer(row,col));}}}}

       //  free up memory for original samplings

 for (itcol=ops.begin();itcol!=ops.end();itcol++){
    for (itrow=ops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,herm,subvev);
       MCObsInfo obskey(corrt,RealPart);
       m_moh->eraseData(obskey);}
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,herm,subvev);
    m_moh->eraseData(MCObsInfo(corrt,RealPart));}
}


#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif



void SinglePivotOfCorrMat::writeRotated(uint tmin, uint tmax, const string& corrfile,
                                        WriteMode wmode, LogHelper& xmlout, char mode,
                                        char file_format)
{
 xmlout.reset("WriteRotated");
 uint nlevels=getNumberOfLevels();
 set<MCObsInfo> obskeys;
 bool vevs=m_cormat_info->subtractVEV();
 if (vevs && (mode!='S')){
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetIDIndex(level);
       MCObsInfo obskey(*m_rotated_info,RealPart);
       obskeys.insert(obskey);
#ifdef COMPLEXNUMBERS
       obskey.setToImaginaryPart();
       obskeys.insert(obskey);
#endif
       }}
 if (mode!='S'){
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetIDIndex(level);
       for (uint timeval=tmin;timeval<=tmax;timeval++){
          MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,true,RealPart,false);
          obskeys.insert(obskey);
#ifdef COMPLEXNUMBERS
          obskey.setToImaginaryPart();
          obskeys.insert(obskey);
#endif
          }}}
 if (vevs && ((mode=='S')||(mode=='A'))){
    for (uint level=0;level<nlevels;level++){
       m_rotated_info->resetIDIndex(level);
       for (uint timeval=tmin;timeval<=tmax;timeval++){
          MCObsInfo obskey(*m_rotated_info,*m_rotated_info,timeval,true,RealPart,true);
          obskeys.insert(obskey);
#ifdef COMPLEXNUMBERS
          obskey.setToImaginaryPart();
          obskeys.insert(obskey);
#endif
          }}}
 XMLHandler xmlf;
 if (mode=='B')
    m_moh->writeBinsToFile(obskeys,corrfile,xmlf,wmode,file_format);
 else
    m_moh->writeSamplingValuesToFile(obskeys,corrfile,xmlf,wmode,file_format);
 xmlout.put(xmlf);
}



void SinglePivotOfCorrMat::insertAmplitudeFitInfo(uint level, const MCObsInfo& ampinfo)
{
 if (!m_reorder.empty())
    throw(std::invalid_argument("Cannot insert Amplitude info after reordering already done"));
 try{
    if (level<getNumberOfLevels())
       m_ampkeys.insert(make_pair(level,ampinfo));}
 catch(std::exception& xp){
    throw(std::invalid_argument(string("Could not insertAmplitudeFitInfo: ")+string(xp.what())));}
}


uint SinglePivotOfCorrMat::get_orig_level(uint level) const
{
 return (m_reorder.empty())?level:m_reorder[level];
}


MCObsInfo SinglePivotOfCorrMat::getAmplitudeKey(uint level) const
{
 if (level>=getNumberOfLevels())
    throw(std::invalid_argument("invalid level index in getAmplitudeKey"));
 std::map<uint,MCObsInfo>::const_iterator it=m_ampkeys.find(get_orig_level(level));
 if (it!=m_ampkeys.end()) return it->second;
 throw(std::invalid_argument("could not find AmplitudeKey"));
}



void SinglePivotOfCorrMat::insertEnergyFitInfo(uint level, const MCObsInfo& energyinfo)
{
 if (!m_reorder.empty())
    throw(std::invalid_argument("Cannot insert Energy info after reordering already done"));
 try{
    if (level<getNumberOfLevels())
       m_energykeys.insert(make_pair(level,energyinfo));}
 catch(std::exception& xp){
    throw(std::invalid_argument(string("Could not insertEnergyFitInfo: ")+string(xp.what())));}
}


MCObsInfo SinglePivotOfCorrMat::getEnergyKey(uint level) const
{
 if (level>=getNumberOfLevels())
    throw(std::invalid_argument("invalid level index in getEnergyKey"));
 std::map<uint,MCObsInfo>::const_iterator it=m_energykeys.find(get_orig_level(level));
 if (it!=m_energykeys.end()) return it->second;
 throw(std::invalid_argument("could not find EnergyKey"));
}


  // This sets m_reorder.  m_reorder[0] contains original level number of lowest energy
  // state, m_reorder[1] contains original level number of first excited state, etc.

void SinglePivotOfCorrMat::reorderLevelsByFitEnergy(LogHelper& xmllog)
{
 m_reorder.clear();
 if (!allEnergyFitInfoAvailable())
    throw(std::runtime_error("Not all Energy fit info available to reorderLevelsByFitEnergy"));
 if (!allAmplitudeFitInfoAvailable())
    throw(std::runtime_error("can reorderLevelsByFitEnergy only after inserting all amplitude infos"));

 try{
 uint nlevels=getNumberOfLevels();
 list<pair<double,uint> > energylevels;
 for (uint level=0;level<nlevels;level++){
    MCObsInfo energyfitkey=getEnergyKey(level);
    double energy=m_moh->getEstimate(energyfitkey).getFullEstimate();
    energylevels.push_back(make_pair(energy,level));}
 energylevels.sort(level_compare);
 list<pair<double,uint> >::const_iterator it;
 m_reorder.resize(nlevels);
 it=energylevels.begin();
 for (uint level=0;level<nlevels;level++,it++){
    m_reorder[level]=it->second;}
 xmllog.reset("ReorderEnergies");
 for (uint level=0;level<nlevels;level++){
    MCObsInfo energyfitkey=getEnergyKey(level);
    LogHelper xmllev("EnergyLevel");
    xmllev.putUInt("LevelIndex",level);
    xmllev.putUInt("OriginalIndex",m_reorder[level]);
    LogHelper result("Energy");
    MCEstimate energy=m_moh->getEstimate(energyfitkey);
    result.putReal("MeanValue",energy.getFullEstimate());
    result.putReal("StandardDeviation",energy.getSymmetricError());
    XMLHandler xmlinfo; energyfitkey.output(xmlinfo);
    xmllev.put(xmlinfo);
    xmllev.put(result);
    xmllog.put(xmllev);
    }


}
 catch(const std::exception& xp){
    throw(std::runtime_error(string("Not all energies known so cannot reorder: ")+xp.what()));}
}


void SinglePivotOfCorrMat::clearReordering()
{
 m_reorder.clear();
}


  //  get |Z(opindex,level)|^2 for all operators for all levels

void SinglePivotOfCorrMat::computeZMagnitudesSquared(Matrix<MCEstimate>& ZMagSq)
{
 if (!allAmplitudeFitInfoAvailable())
    throw(std::runtime_error("Not all Amplitude fit info available to compute ZMagSquares"));

 uint nops=getNumberOfOperators();
 uint nlevels=getNumberOfLevels();

   //  Method:
   //  each diagonal element of rotated correlator matrix is
   //  fit to  A * exp(-E_n t )    then   Zrotsq = std::abs(A)

 try{
 ZMagSq.resize(nops,nlevels);   //  final results put in here
 bool overwrite=true;
 MCObsInfo obskey("TempZMagSq",0);   // temporary key
 for (uint level=0;level<nlevels;level++){
    MCObsInfo Zrotfitkey=getAmplitudeKey(level);
    uint origlevel=get_orig_level(level);
    for (m_moh->begin();!m_moh->end();++(*m_moh)){   // loop over resamplings
       double Zrotsq=std::abs(m_moh->getCurrentSamplingValue(Zrotfitkey));
       for (uint opindex=0;opindex<nops;opindex++){
          obskey.resetObsIndex(opindex);
          m_moh->putCurrentSamplingValue(obskey,
                  sqr((*m_Zmat)(opindex,origlevel))*Zrotsq,overwrite);}}
    for (uint opindex=0;opindex<nops;opindex++){
       obskey.resetObsIndex(opindex);
       ZMagSq(opindex,level)=m_moh->getEstimate(obskey);
       m_moh->eraseSamplings(obskey); }}}
 catch(const std::exception& xp){
    throw(std::runtime_error(string("Not all fit amplitudes known so cannot Z factors: ")+xp.what()));}

}

 // ******************************************************************
