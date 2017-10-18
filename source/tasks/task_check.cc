#include "task_handler.h"
#include "correlator_matrix_info.h"
#include "histogram.h"


using namespace std;

// *******************************************************************************
// *                                                                             *
// *    XML format for correlator matrix checking: checks that correlators       *
// *    and VEVs, if required, are present, and checks each observable for       *
// *    "outliers", which might indicate corrupt data.  Outliers are identified  *
// *    as follows: the data is sorted, then the data range "a..b" for the middle*
// *    half of the data points is found.  Let the mid-range be "m=(a+b)/2"      *
// *    then outliers are points outside the range  "m-v .. m+v"  where          *
// *    "v = outlier_scale * (b-a)/2 ".                                          *
// *                                                                             *
// *    If <Type>TemporalCorrelatorMatrix</Type> is used, this task also         *
// *    checks to see if any vevs or diagonal correlators at minimum time        *
// *    separation are consistent with zero.                                     *
// *                                                                             *
// *    If <Type>TemporalCorrelatorMatrixIsHermitian</Type> is used, it checks   *
// *    that the correlator matrix is hermitian.                                 *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoChecks</Action>                                               *
// *       <Type>TemporalCorrelatorMatrix</Type>                                 *
// *       <CorrelatorMatrix>                                                    *
// *           <Operator>...</Operator>                                          *
// *           <Operator>...</Operator>                                          *
// *                ...                                                          *
// *           <HermitianMatrix\>    (optional)                                  *
// *           <SubtractVEV\>    (optional)                                      *
// *       </CorrelatorMatrix>                                                   *
// *       <MinTimeSep>3</MinTimeSep>                                            *
// *       <MaxTimeSep>25</MaxTimeSep>                                           *
// *       <Verbose/>                        (optional)                          *
// *       <OutlierScale>12.5</OutlierScale> (optional: default 9.5)             *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoChecks</Action>                                               *
// *       <Type>TemporalCorrelatorMatrixIsHermitian</Type>                      *
// *       <CorrelatorMatrix>                                                    *
// *           <Operator>...</Operator>                                          *
// *           <Operator>...</Operator>                                          *
// *                ...                                                          *
// *       </CorrelatorMatrix>                                                   *
// *       <MinTimeSep>3</MinTimeSep>                                            *
// *       <MaxTimeSep>25</MaxTimeSep>                                           *
// *       <Verbose/>                        (optional)                          *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *******************************************************************************



void do_obs_check(MCObsHandler *m_obs, const MCObsInfo& obskey, XMLHandler& xmlout,
                  double outlier_scale, bool verbose, bool& checkflag)
{
 if (!m_obs->queryBins(obskey)){
    if (verbose){
       XMLHandler xmlc; obskey.output(xmlc,false);
       XMLHandler xmlres; xmlres.set_root("Missing");
       xmlres.put_child(xmlc);
       xmlout.put_sibling(xmlres);}
    checkflag=false;}
 else{
    const Vector<double>& bins=m_obs->getBins(obskey);
    for (uint k=0;k<bins.size();k++)
       if (std::isnan(bins[k])){
          XMLHandler xmlc; obskey.output(xmlc,false);
          XMLHandler xmlres; xmlres.set_root("NaNValuesPresent");
          xmlres.put_child(xmlc);
          xmlout.put_sibling(xmlres);
          checkflag=false;
          m_obs->eraseData(obskey);
          return;}
    Vector<uint> outliers;
    getOutliers(bins,outliers,outlier_scale);
    if (outliers.size()>0){
       checkflag=false;
       if (verbose){
          XMLHandler xmlc; obskey.output(xmlc,false);
          XMLHandler xmlres; xmlres.set_root("Outliers");
          xmlres.put_child(xmlc);
          xmlres.seek_first_child();
          for (uint k=0;k<outliers.size();k++)
             xmlres.put_sibling("SerialIndex",make_string(outliers[k]));
          xmlout.put_sibling(xmlres);}}
    m_obs->eraseData(obskey);}
}



void do_snk_src_check(MCObsHandler *m_obs, const OperatorInfo& snk, const OperatorInfo& src,
                      uint mintimesep, uint maxtimesep,
                      XMLHandler& xmlout, bool verbose, uint& total, 
                      uint& onesig, uint& twosig, uint& foursig, 
                      uint& missing)
{
 CorrelatorAtTimeInfo corA(snk,src,0,false,false);
 CorrelatorAtTimeInfo corB(src,snk,0,false,false);
 for (uint t=mintimesep;t<=maxtimesep;t++){
    corA.resetTimeSeparation(t);
    corB.resetTimeSeparation(t);
    MCObsInfo obsAre(corA,RealPart);
    MCObsInfo obsAim(corA,ImaginaryPart);
    MCObsInfo obsBre(corB,RealPart);
    MCObsInfo obsBim(corB,ImaginaryPart);
    bool availA=m_obs->queryBins(obsAre)&&m_obs->queryBins(obsAim);
    bool availB=m_obs->queryBins(obsBre)&&m_obs->queryBins(obsBim);
    if ((!availA)&&(!availB)){
       missing++;
       if (verbose){
          XMLHandler xmlsnk; snk.output(xmlsnk,false);
          XMLHandler xmlsrc; src.output(xmlsrc,false);
          XMLHandler xmlres; xmlres.set_root("MissingCorrelation");
          xmlres.put_child(xmlsnk);
          xmlres.put_child(xmlsrc);
          xmlres.put_child("TimeSeparation",make_string(t));
          xmlout.put_sibling(xmlres);}}
    else if (availA && availB){
       const Vector<double>& bufAre=m_obs->getBins(obsAre);
       const Vector<double>& bufBre=m_obs->getBins(obsBre);
       uint nbins=m_obs->getNumberOfBins();
       Vector<double> buf(nbins);
       for (uint k=0;k<nbins;k++){
          buf[k]=bufAre[k]-bufBre[k];}
       MCObsInfo oredif("RealPartDiff",0,true);
       m_obs->putBins(oredif,buf);
       MCEstimate mredif=m_obs->getJackknifeEstimate(oredif);
       double df1=std::abs(mredif.getFullEstimate())/mredif.getSymmetricError();
       total++;
       if (df1>1.0) onesig++;
       if (df1>2.0) twosig++;
       if (df1>4.0) foursig++;
       m_obs->eraseData(oredif);
       const Vector<double>& bufAim=m_obs->getBins(obsAim);
       const Vector<double>& bufBim=m_obs->getBins(obsBim);
       for (uint k=0;k<nbins;k++){
          buf[k]=bufAim[k]+bufBim[k];}
       MCObsInfo oimsum("ImagPartSum",0,true);
       m_obs->putBins(oimsum,buf);
       MCEstimate mimsum=m_obs->getJackknifeEstimate(oimsum);
       double df2=std::abs(mimsum.getFullEstimate())/mimsum.getSymmetricError();
       total++;
       if (df2>1.0) onesig++;
       if (df2>2.0) twosig++;
       if (df2>4.0) foursig++;
       m_obs->eraseData(oimsum);
       if (((df1>4.0)||(df2>4.0))&&(verbose)){
          XMLHandler xmlsnk; snk.output(xmlsnk,false);
          XMLHandler xmlsrc; src.output(xmlsrc,false);
          XMLHandler xmlres; xmlres.set_root("Discrepancy");
          xmlres.put_child(xmlsnk);
          xmlres.put_child(xmlsrc);
          xmlres.put_child("TimeSeparation",make_string(t));
          xmlres.put_child("DiffRealPartRatio",make_string(df1));
          xmlres.put_child("SumImagPartRatio",make_string(df2));
          xmlout.put_sibling(xmlres);}}
    m_obs->eraseData(obsAre);
    m_obs->eraseData(obsAim);
    m_obs->eraseData(obsBre);
    m_obs->eraseData(obsBim);
    }
}



void is_obs_zero(MCObsHandler *m_obs, const MCObsInfo& obskey,
                 XMLHandler& xmlout, bool& checkflag)
{
 MCObsInfo obskey1(obskey); obskey1.setToRealPart();
 MCObsInfo obskey2(obskey); obskey2.setToImaginaryPart();
 if ((!m_obs->queryBins(obskey1))&&(!m_obs->queryBins(obskey2))) return;
 double re_mean=0.0; double re_err=1.0;
 double im_mean=0.0; double im_err=1.0;
 if (m_obs->queryBins(obskey1)){
    MCEstimate est=m_obs->getJackknifeEstimate(obskey1);
    re_mean=std::abs(est.getFullEstimate());
    re_err=est.getSymmetricError();}
 if (m_obs->queryBins(obskey2)){
    MCEstimate est=m_obs->getJackknifeEstimate(obskey2);
    im_mean=std::abs(est.getFullEstimate());
    im_err=est.getSymmetricError();}
 if ((re_mean<=4.0*re_err)&&(im_mean<=4.0*im_err)){
    checkflag=false;
    XMLHandler xmlc; obskey.output(xmlc,false);
    XMLHandler xmlres; xmlres.set_root("Zero");
    xmlres.put_child(xmlc);
    xmlout.put_sibling(xmlres);}
 m_obs->eraseData(obskey1);
 m_obs->eraseData(obskey2);
}




void TaskHandler::doChecks(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 string printtype;
 xmlread(xmltask,"Type",printtype,"DoChecks");

 if (printtype=="TemporalCorrelatorMatrix"){
    try{
    xmlout.set_root("DoChecks"); 
    bool checkflag=true;
    uint mintimesep,maxtimesep;
    xmlreadchild(xmltask,"MinTimeSep",mintimesep);
    xmlreadchild(xmltask,"MaxTimeSep",maxtimesep);
    XMLHandler xmlcm(xmltask,"CorrelatorMatrixInfo");
    CorrelatorMatrixInfo cormat(xmlcm);
    bool herm=cormat.isHermitian();
    bool vevs=cormat.subtractVEV();  
    bool verbose=(xml_tag_count(xmltask,"Verbose")>0);
    double outlier_scale=9.5;
    xmlreadifchild(xmltask,"OutlierScale",outlier_scale);
    XMLHandler xmlt; cormat.output(xmlt,false);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("MinTimeSep",make_string(mintimesep));
    xmlout.put_sibling("MaxTimeSep",make_string(maxtimesep));
    xmlout.put_sibling("OutlierScale",make_string(outlier_scale));
    vector<ComplexArg> args(2);
    args[0]=RealPart; args[1]=ImaginaryPart;
    for (cormat.begin();!cormat.end();++cormat){
       const CorrelatorInfo& cor=cormat.getCurrentCorrelatorInfo();
       CorrelatorAtTimeInfo cort(cor,0,herm,false);
       for (uint t=mintimesep;t<=maxtimesep;t++){
          cort.resetTimeSeparation(t);
          for (uint kk=0;kk<args.size();kk++){
             MCObsInfo obskey(cort,args[kk]);
             do_obs_check(m_obs,obskey,xmlout,outlier_scale,verbose,checkflag);}}
       if (cor.isSinkSourceSame()){
          cort.resetTimeSeparation(mintimesep);
          MCObsInfo obskey(cort,RealPart);
          is_obs_zero(m_obs,obskey,xmlout,checkflag);}}
    if (vevs){
       const set<OperatorInfo>& vevops=cormat.getOperators();
       for (set<OperatorInfo>::const_iterator it=vevops.begin();it!=vevops.end();it++){
          for (uint kk=0;kk<args.size();kk++){
             MCObsInfo obskey(*it,args[kk]);
             do_obs_check(m_obs,obskey,xmlout,outlier_scale,verbose,checkflag);}
             MCObsInfo obskey2(*it,RealPart);
             is_obs_zero(m_obs,obskey2,xmlout,checkflag);}}
    if (checkflag){
       xmlout.put_sibling("Status","All checks PASSED");}
    else{
       xmlout.put_sibling("Status","Some checks FAILED");}
    }

    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("doChecks with TemporalCorrelatorMatrix type encountered an error: ")
                +string(errmsg.what())));}
    }


 else if (printtype=="TemporalCorrelatorMatrixIsHermitian"){
    try{
    xmlout.set_root("DoChecks"); 
    uint total=0, onesig=0, twosig=0, foursig=0, missing=0;
    uint mintimesep,maxtimesep;
    xmlreadchild(xmltask,"MinTimeSep",mintimesep);
    xmlreadchild(xmltask,"MaxTimeSep",maxtimesep);
    XMLHandler xmlcm(xmltask,"CorrelatorMatrixInfo");
    CorrelatorMatrixInfo cormat(xmlcm);
    bool verbose=(xml_tag_count(xmltask,"Verbose")>0);
    XMLHandler xmlt; cormat.output(xmlt,false);
    xmlout.put_child(xmlt);
    xmlout.seek_first_child();
    xmlout.put_sibling("MinTimeSep",make_string(mintimesep));
    xmlout.put_sibling("MaxTimeSep",make_string(maxtimesep));
    xmlout.put_sibling("Type","CheckIfHermitian");
    const set<OperatorInfo>& corrops=cormat.getOperators();
    for (set<OperatorInfo>::const_iterator snk=corrops.begin();snk!=corrops.end();snk++)
    for (set<OperatorInfo>::const_iterator src=snk;src!=corrops.end();src++){
       do_snk_src_check(m_obs,*snk,*src,mintimesep,maxtimesep,xmlout,verbose,total, 
                        onesig,twosig,foursig,missing);}
    XMLHandler xmls("Status");
    xmls.put_child("NumberMissing",make_string(missing));
    xmls.put_child("TotalNumberTests",make_string(total));
    if (total>0){
        xmls.put_child("FractionOneSigma",make_string(double(onesig)/double(total)));
        xmls.put_child("FractionTwoSigma",make_string(double(twosig)/double(total)));
        xmls.put_child("FractionFourSigma",make_string(double(foursig)/double(total)));}
    xmlout.put_sibling(xmls);
    }
    catch(const std::exception& errmsg){
       xmlout.clear();
       throw(std::invalid_argument(string("doChecks with TemporalCorrelatorMatrixIsHermitian type encountered an error: ")
           +string(errmsg.what())));}
    }


 else{
    throw(std::invalid_argument("doChecks encountered unsupported type: "));}


}


// ***************************************************************************************
 
