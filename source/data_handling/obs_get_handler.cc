#include "obs_get_handler.h"
#include "filelist_info.h"
#include "correlator_matrix_info.h"

using namespace std;

// *************************************************************************


MCObsGetHandler::MCObsGetHandler(XMLHandler& xmlin, const MCBinsInfo& bins_info, 
                                 const MCSamplingInfo& samp_info)
                     : m_corrdh(0), m_vevdh(0), m_binsdh(0), m_sampsdh(0),
                       m_bins_info(bins_info), m_sampling_info(samp_info),
                       m_use_checksums(false)
{
 XMLHandler xmlr(xmlin);
 list<FileListInfo> corrinputfiles;
 list<FileListInfo> vevinputfiles;
 set<string> binfiles;
 set<string> sampfiles;
 xml_tag_assert(xmlr,"MCObservables");

    //  read the input file lists to find the data

 if (xmlr.query_unique_to_among_children("BLCorrelatorData")){
    XMLHandler xmlp(xmlr,"BLCorrelatorData");
    {list<XMLHandler> infxml;
    infxml=xmlp.find("FileListInfo");  
    for (list<XMLHandler>::iterator 
        it=infxml.begin();it!=infxml.end();++it)
       corrinputfiles.push_back(FileListInfo(*it));}}

 if (xmlr.query_unique_to_among_children("BLVEVData")){
    XMLHandler xmlp(xmlr,"BLVEVData");
    {list<XMLHandler> infxml;
    infxml=xmlp.find("FileListInfo");  
    for (list<XMLHandler>::iterator 
        it=infxml.begin();it!=infxml.end();++it)
       vevinputfiles.push_back(FileListInfo(*it));}}

 if (xmlr.query_unique_to_among_children("BinData")){
    XMLHandler xmlp(xmlr,"BinData");
    {list<XMLHandler> infxml;
    infxml=xmlp.find("FileName");  
    for (list<XMLHandler>::iterator 
        it=infxml.begin();it!=infxml.end();++it){
       ArgsHandler xmls(*it);
       string fname(xmls.getString("FileName"));
       binfiles.insert(fname);}}}

 if (xmlr.query_unique_to_among_children("SamplingData")){
    XMLHandler xmlp(xmlr,"SamplingData");
    {list<XMLHandler> infxml;
    infxml=xmlp.find("FileName");  
    for (list<XMLHandler>::iterator 
        it=infxml.begin();it!=infxml.end();++it){
       ArgsHandler xmls(*it);
       string fname(xmls.getString("FileName"));
       sampfiles.insert(fname);}}}

/*
 cout << "Correlator Input Files:"<<endl;
 for (list<FileListInfo>::iterator it=corrinputfiles.begin();it!=corrinputfiles.end();it++)
    cout << it->output()<<endl;
 cout << "VEV Input Files:"<<endl;
 for (list<FileListInfo>::iterator it=vevinputfiles.begin();it!=vevinputfiles.end();it++)
    cout << it->output()<<endl;
 cout << "Bin Input Files:"<<endl;
 for (set<string>::iterator it=binfiles.begin();it!=binfiles.end();it++)
    cout << *it<<endl;
 cout << "Sampling Input Files:"<<endl;
 for (set<string>::iterator it=sampfiles.begin();it!=sampfiles.end();it++)
    cout << *it<<endl;
*/

 bool nodata=(corrinputfiles.empty())&&(vevinputfiles.empty())
             &&(binfiles.empty())&&(sampfiles.empty());

 if (xml_child_tag_count(xmlr,"UseCheckSums")>1){
    m_use_checksums=true;}

    //  get and rebin the weights, if ensemble is weighted

 m_is_weighted=m_bins_info.getMCEnsembleInfo().isWeighted();
 if (m_is_weighted)
    get_the_weights();

 if (nodata) return;

    //  read requested observables (optionally given)

 try{

 set<CorrelatorInfo> corrSetNoSym;
 set<CorrelatorInfo> corrSetSym;
 set<OperatorInfo> vevSet;
 set<MCObsInfo> obsbinSet;
 set<MCObsInfo> obssampSet;

 int speccount=xml_tag_count(xmlr,"Specifications");
 if (speccount>1){
    throw(std::invalid_argument("Multiple <Specifications> tags not allowed"));}
 else if (speccount==1){

    bool allhermitian=false;
    XMLHandler xmlspec(xmlr,"Specifications");
    if (xml_tag_count(xmlspec,"AllHermitian")>=1){
       allhermitian=true;}
    set<CorrelatorInfo>& cset=(allhermitian)? corrSetSym : corrSetNoSym;

          // individual correlators
    setup_correlators(xmlspec,"Correlator",false,cset,vevSet);
          // individual correlators with VEV subtraction
    setup_correlators(xmlspec,"CorrelatorWithVEV",true,cset,vevSet);
          // individual VEVs
    setup_vevs(xmlspec,"VEV",vevSet);
          // Hermitian correlation matrices
    setup_correlator_matrices(xmlspec,"HermitianCorrelationMatrix",
                              true,false,corrSetSym,vevSet);
          // Hermitian correlation matrices with VEVs
    setup_correlator_matrices(xmlspec,"HermitianCorrelationMatrixWithVEVs",
                              true,true,corrSetSym,vevSet);
          // non-Hermitian correlation matrices
    setup_correlator_matrices(xmlspec,"CorrelationMatrix",
                              allhermitian,false,cset,vevSet);
          // non-Hermitian correlation matrices with VEVs
    setup_correlator_matrices(xmlspec,"CorrelationMatrixWithVEVs",
                              allhermitian,true,cset,vevSet);
          // sigmond observables in bin files
    setup_obsset(xmlspec,"ObsBins",obsbinSet);
          // sigmond observables in sampling files
    setup_obsset(xmlspec,"ObsSamplings",obssampSet);
    }

// bool cspec=(!corrSetNoSym.empty())||(!corrSetSym.empty());
// bool vspec=!vevSet.empty();
// if ((cspec)&&(corrinputfiles.empty()))
//    throw(std::invalid_argument("No correlator input files but correlators requested"));
// if ((vspec)&&(vevinputfiles.empty()))
//    throw(std::invalid_argument("No VEV input files but VEVs requested"));

// if ((!vevinputfiles.empty()) && (vspec || !cspec)){
 if (!vevinputfiles.empty()){
    m_vevdh=new LaphEnv::BLVEVDataHandler(
                vevinputfiles,vevSet,&(m_bins_info.getMCEnsembleInfo()),m_use_checksums);}
// if ((!corrinputfiles.empty()) && (cspec || !vspec)){
 if (!corrinputfiles.empty()){
    m_corrdh=new LaphEnv::BLCorrelatorDataHandler(
                corrinputfiles,corrSetNoSym,
                corrSetSym,&(m_bins_info.getMCEnsembleInfo()),m_use_checksums);}

 bool bspec=!obsbinSet.empty();
// if ((bspec)&&(binfiles.empty()))
//    throw(std::invalid_argument("No bin input files but observables requested"));
 if (!binfiles.empty()){
    m_binsdh=new BinsGetHandler(m_bins_info,binfiles,m_use_checksums);
    if ((m_binsdh!=0)&&(bspec))
       if (!(m_binsdh->keepKeys(obsbinSet)))
          throw(std::runtime_error("Requested observable not available in input bin files"));}

 bool sspec=!obssampSet.empty();
 if ((sspec)&&(sampfiles.empty()))
    throw(std::invalid_argument("No sampling input files but observables requested"));
 if (!sampfiles.empty()){
    m_sampsdh=new SamplingsGetHandler(m_bins_info,m_sampling_info, 
                                      sampfiles,m_use_checksums);
    if ((m_sampsdh!=0)&&(sspec))
       if (!(m_sampsdh->keepKeys(obssampSet)))
          throw(std::runtime_error("Requested observable not available in input sampling files"));}

// if ((m_corrdh==0)&&(m_vevdh==0)&&(m_binsdh==0)&&(m_sampsdh==0))
//    throw(std::invalid_argument("No input data/samplings: nothing to do"));

 }
 catch(const std::exception& msg){
    clear(); throw;}
}



MCObsGetHandler::~MCObsGetHandler()
{
 clear();
}


void MCObsGetHandler::connectBinsFile(const std::string& file_name)
{
 try{
    if (m_binsdh==0){
       set<string> binfiles; binfiles.insert(file_name);
       m_binsdh=new BinsGetHandler(m_bins_info,binfiles,m_use_checksums);}
    else{
       m_binsdh->addFile(file_name);}}
 catch(const std::exception& msg){
    clear(); throw;}
}


void MCObsGetHandler::connectBinsFile(const std::string& file_name, 
                                      const std::set<MCObsInfo>& keys_to_keep)
{
 try{
    if (m_binsdh==0){
       set<string> binfiles; binfiles.insert(file_name);
       m_binsdh=new BinsGetHandler(m_bins_info,binfiles,m_use_checksums);
       if (!(m_binsdh->keepKeys(keys_to_keep)))
          throw(std::runtime_error("Requested observable not available in input bin files"));}
    else{
       if (!(m_binsdh->addFile(file_name,keys_to_keep)))
          throw(std::runtime_error("Requested observable not available in input bin files"));}}
 catch(const std::exception& msg){
    clear(); throw;}
}


void MCObsGetHandler::disconnectBinsFile(const std::string& file_name)
{
 if (m_binsdh!=0)
    m_binsdh->removeFile(file_name);
}


void MCObsGetHandler::connectSamplingsFile(const std::string& file_name)
{
 try{
    if (m_sampsdh==0){
       set<string> sampfiles; sampfiles.insert(file_name);
       m_sampsdh=new SamplingsGetHandler(m_bins_info,m_sampling_info, 
                                         sampfiles,m_use_checksums);}
    else{
       m_sampsdh->addFile(file_name);}}
 catch(const std::exception& msg){
    clear(); throw;}
}

void MCObsGetHandler::connectSamplingsFile(const std::string& file_name, 
                                           const std::set<MCObsInfo>& keys_to_keep)
{
 try{
    if (m_sampsdh==0){
       set<string> sampfiles; sampfiles.insert(file_name);
       m_sampsdh=new SamplingsGetHandler(m_bins_info,m_sampling_info, 
                                         sampfiles,m_use_checksums);
       if (!(m_sampsdh->keepKeys(keys_to_keep)))
          throw(std::runtime_error("Requested observable not available in input sampling files"));}
    else{
       if (!(m_sampsdh->addFile(file_name,keys_to_keep)))
          throw(std::runtime_error("Requested observable not available in input sampling files"));}}
 catch(const std::exception& msg){
    clear(); throw;}
}

void MCObsGetHandler::disconnectSamplingsFile(const std::string& file_name)
{
 if (m_sampsdh!=0)
    m_sampsdh->removeFile(file_name);
}



void MCObsGetHandler::clear()
{
 delete m_corrdh; m_corrdh=0;
 delete m_vevdh; m_vevdh=0;
 delete m_binsdh; m_binsdh=0;
 delete m_sampsdh; m_sampsdh=0;
 m_wts.clear(); m_BLorig_wts.clear();
 m_is_weighted=false;
}

string MCObsGetHandler::getEnsembleId() const
{
 return (m_bins_info.getMCEnsembleInfo()).getId();
}


MCEnsembleInfo MCObsGetHandler::getEnsembleInfo() const
{
 return m_bins_info.getMCEnsembleInfo();
}


unsigned int MCObsGetHandler::getNumberOfMeasurements() const
{
 return m_bins_info.getNumberOfMeasurements();
}


unsigned int MCObsGetHandler::getNumberOfDefaultResamplings() const
{
 return m_sampling_info.getNumberOfReSamplings(m_bins_info);
}


SamplingMode MCObsGetHandler::getDefaultSamplingMode() const
{
 return m_sampling_info.getSamplingMode();
}


const MCBinsInfo& MCObsGetHandler::getBinsInfo() const
{
 return m_bins_info;
}


const MCSamplingInfo& MCObsGetHandler::getSamplingInfo() const
{
 return m_sampling_info;
}


bool MCObsGetHandler::useCheckSums() const
{
 return m_use_checksums;
}


void MCObsGetHandler::getBasicLapHData(const MCObsInfo& obsinfo, 
                                       int serial_index, Scalar& data)
{
 if (!((obsinfo.isBasicLapH())&&(obsinfo.isSimple())))
    throw(std::invalid_argument(string("getBasicLapHData failed for ")+obsinfo.str()));
 if (obsinfo.isHermitianCorrelatorAtTime()){
    if (m_corrdh==0) throw(std::invalid_argument(string("getBasicLapHData failed for ")+obsinfo.str()));
    m_corrdh->getSymData(obsinfo.getCorrelatorAtTimeInfo(),serial_index,data);}
 else if (obsinfo.isCorrelatorAtTime()){
    if (m_corrdh==0) throw(std::invalid_argument(string("getBasicLapHData failed for ")+obsinfo.str()));
    m_corrdh->getData(obsinfo.getCorrelatorAtTimeInfo(),serial_index,data);}
 else if (obsinfo.isVEV()){
    if (m_vevdh==0) throw(std::invalid_argument(string("getBasicLapHData failed for ")+obsinfo.str()));
    m_vevdh->getData(obsinfo.getVEVInfo(),serial_index,data);}
 else
    throw(std::invalid_argument(string("getBasicLapHData failed for ")+obsinfo.str()));
}


bool MCObsGetHandler::queryBasicLapHData(const MCObsInfo& obsinfo, int serial_index)
{
 if (!((obsinfo.isBasicLapH())&&(obsinfo.isSimple()))) return false;
 if (obsinfo.isHermitianCorrelatorAtTime()){
    if (m_corrdh==0) return false;
    return m_corrdh->querySymData(obsinfo.getCorrelatorAtTimeInfo(),serial_index);}
 else if (obsinfo.isCorrelatorAtTime()){
    if (m_corrdh==0) return false;
    return m_corrdh->queryData(obsinfo.getCorrelatorAtTimeInfo(),serial_index);}
 else if (obsinfo.isVEV()){
    if (m_vevdh==0) return false;
    return m_vevdh->queryData(obsinfo.getVEVInfo(),serial_index);}
 return false;
}




#ifdef COMPLEXNUMBERS

void MCObsGetHandler::getBins(const MCObsInfo& obsinfo, RVector& bins)
{
 if (obsinfo.isNonSimple())
    throw(std::invalid_argument(string("cannot getBins for non simple observable for ")+obsinfo.str()));
 if (obsinfo.isBasicLapH()){
    try{
    RVector bin_discard;
    RVector *b1, *b2;
    if (obsinfo.isRealPart()){
       b1=&bins; b2=&bin_discard;}
    else{
       b2=&bins; b1=&bin_discard;}
    if (obsinfo.isHermitianCorrelatorAtTime()){
       if (m_corrdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHCorrSymGetter p(obsinfo,m_corrdh);
       get_bl_data(p,*b1,*b2); return;}
    else if (obsinfo.isCorrelatorAtTime()){
       if (m_corrdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHCorrGetter p(obsinfo,m_corrdh);
       get_bl_data(p,*b1,*b2); return;}
    else if (obsinfo.isVEV()){
       if (m_vevdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHVEVGetter p(obsinfo,m_vevdh);
       get_bl_data(p,*b1,*b2); return;}}
    catch(const std::exception& xp){}}
 get_sig_bin_data(obsinfo,bins);
}

#else

void MCObsGetHandler::getBins(const MCObsInfo& obsinfo, RVector& bins)
{
 if (obsinfo.isNonSimple())
    throw(std::invalid_argument(string("cannot getBins for non simple observable for ")+obsinfo.str()));
 if (obsinfo.isBasicLapH()){
    try{
    if (obsinfo.isImaginaryPart())
       throw(std::invalid_argument(string("cannot getBins for observable for ")+obsinfo.str()));
    if (obsinfo.isHermitianCorrelatorAtTime()){
       if (m_corrdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHCorrSymGetter p(obsinfo,m_corrdh);
       get_bl_data(p,bins); return;}
    else if (obsinfo.isCorrelatorAtTime()){
       if (m_corrdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHCorrGetter p(obsinfo,m_corrdh);
       get_bl_data(p,bins); return;}
    else if (obsinfo.isVEV()){
       if (m_vevdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHVEVGetter p(obsinfo,m_vevdh);
       get_bl_data(p,bins); return;}}
    catch(const std::exception& xp){}}
 get_sig_bin_data(obsinfo,bins);
}

#endif

void MCObsGetHandler::getSamplings(const MCObsInfo& obsinfo, RVector& samplings)
{
 if (m_sampsdh==0)
       throw(std::invalid_argument(string("getSamplings fails due to unavailable sampling for ")+obsinfo.str()));
 m_sampsdh->getSymData(obsinfo,samplings);
}


bool MCObsGetHandler::getSamplingsMaybe(const MCObsInfo& obsinfo, RVector& samplings)
{
 samplings.clear();
 if (m_sampsdh==0) return false;
 return m_sampsdh->getSymDataMaybe(obsinfo,samplings);
}


bool MCObsGetHandler::query_bl_bins(const MCObsInfo& obsinfo)
{
 if (obsinfo.isHermitianCorrelatorAtTime()){
    if (m_corrdh==0) return false;
    BasicLapHCorrSymGetter p(obsinfo,m_corrdh);
    return query_bl_data(p);}
 else if (obsinfo.isCorrelatorAtTime()){
    if (m_corrdh==0) return false;
    BasicLapHCorrGetter p(obsinfo,m_corrdh);
    return query_bl_data(p);}
 else if (obsinfo.isVEV()){
    if (m_vevdh==0) return false;
    BasicLapHVEVGetter p(obsinfo,m_vevdh);
    return query_bl_data(p);}
 return false;
}


bool MCObsGetHandler::queryBins(const MCObsInfo& obsinfo)
{
 if (obsinfo.isNonSimple()) return false;
 if (obsinfo.isBasicLapH()){
#ifdef REALNUMBERS
    if (obsinfo.isImaginaryPart()) return false;
#endif
    if (query_bl_bins(obsinfo)) return true;}
 if (m_binsdh==0) return false;
 return m_binsdh->querySymData(obsinfo);
}


bool MCObsGetHandler::querySamplings(const MCObsInfo& obsinfo)
{
 if (m_sampsdh==0) return false;
 return m_sampsdh->querySymData(obsinfo);
}


#ifdef COMPLEXNUMBERS

void MCObsGetHandler::getBinsComplex(const MCObsInfo& obsinfo, RVector& bins_re, 
                                     RVector& bins_im)
{
 if (obsinfo.isNonSimple())
    throw(std::invalid_argument(string("cannot getBins for non simple observable for ")+obsinfo.str()));
 if (obsinfo.isBasicLapH()){
    try{
    if (obsinfo.isHermitianCorrelatorAtTime()){
       if (m_corrdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHCorrSymGetter p(obsinfo,m_corrdh);
       get_bl_data(p,bins_re,bins_im); return;}
    else if (obsinfo.isCorrelatorAtTime()){
       if (m_corrdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHCorrGetter p(obsinfo,m_corrdh);
       get_bl_data(p,bins_re,bins_im); return;}
    else if (obsinfo.isVEV()){
       if (m_vevdh==0) throw(std::invalid_argument(string("getData failed for ")+obsinfo.str()));
       BasicLapHVEVGetter p(obsinfo,m_vevdh);
       get_bl_data(p,bins_re,bins_im); return;}}
    catch(const std::exception& xp){}}
 if (obsinfo.isRealPart()){
    get_sig_bin_data(obsinfo,bins_re);
    MCObsInfo oim(obsinfo); oim.setToImaginaryPart();
    get_sig_bin_data(oim,bins_im);}
 else{
    get_sig_bin_data(obsinfo,bins_im);
    MCObsInfo ore(obsinfo); ore.setToRealPart();
    get_sig_bin_data(ore,bins_re);}
}


void MCObsGetHandler::getSamplingsComplex(const MCObsInfo& obsinfo, RVector& samp_re, 
                                          RVector& samp_im)
{
 if (obsinfo.isRealPart()){
    getSamplings(obsinfo,samp_re);
    MCObsInfo oim(obsinfo); oim.setToImaginaryPart();
    getSamplings(oim,samp_im);}
 else{
    getSamplings(obsinfo,samp_im);
    MCObsInfo ore(obsinfo); ore.setToRealPart();
    getSamplings(ore,samp_re);}
}

bool MCObsGetHandler::querySamplingsComplex(const MCObsInfo& obsinfo)
{
 MCObsInfo obs(obsinfo);
 obs.setToRealPart();
 if (!querySamplings(obs)) return false;
 obs.setToImaginaryPart();
 return querySamplings(obs);
}

bool MCObsGetHandler::queryBinsComplex(const MCObsInfo& obsinfo)
{
 MCObsInfo obs(obsinfo);
 obs.setToRealPart();
 if (!queryBins(obs)) return false;
 obs.setToImaginaryPart();
 return queryBins(obs);
}

#endif


void MCObsGetHandler::close()
{
 if (m_corrdh) m_corrdh->close();
 if (m_vevdh) m_vevdh->close();
 if (m_binsdh) m_binsdh->close();
 if (m_sampsdh) m_sampsdh->close();
}


void MCObsGetHandler::getFileMap(XMLHandler& xmlout) const
{
 if (m_corrdh){
    m_corrdh->getFileMap(xmlout);}
 else{
    xmlout.set_root("FileMap");}
 if (m_vevdh){
    XMLHandler xmlv;
    m_vevdh->getFileMap(xmlv);
    list<XMLHandler> ventries=xmlv.find("Entry");
    for (list<XMLHandler>::const_iterator 
         it=ventries.begin();it!=ventries.end();it++)
       xmlout.put_child(*it);
    }
 if (m_binsdh){
    XMLHandler xmlb;
    m_binsdh->getFileMap(xmlb);
    list<XMLHandler> bentries=xmlb.find("Entry");
    for (list<XMLHandler>::const_iterator 
         it=bentries.begin();it!=bentries.end();it++)
       xmlout.put_child(*it);
    }
 if (m_sampsdh){
    XMLHandler xmls;
    m_sampsdh->getFileMap(xmls);
    list<XMLHandler> sentries=xmls.find("Entry");
    for (list<XMLHandler>::const_iterator 
         it=sentries.begin();it!=sentries.end();it++)
       xmlout.put_child(*it);
    }
}


set<OperatorInfo> MCObsGetHandler::getVEVInfos() const
{
 if (m_vevdh){
    return m_vevdh->getOperatorSet();}
 return set<OperatorInfo>();
}


set<CorrelatorInfo> MCObsGetHandler::getCorrelatorInfos() const
{
 if (m_corrdh){
    return m_corrdh->getCorrelatorSet();}
 return set<CorrelatorInfo>();
}




   //   private members

           // read individual correlators with/without VEV subtraction

void MCObsGetHandler::setup_correlators(XMLHandler& xmlin, 
                          const string& tagname, bool vevs, 
                          set<CorrelatorInfo>& corrSet,
                          set<OperatorInfo>& vevSet)
{
 try{
    list<XMLHandler> corrxml;
    corrxml=xmlin.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       it=corrxml.begin();it!=corrxml.end();++it){
       CorrelatorInfo corr(*it);
       corrSet.insert(corr);
       if (vevs){
          vevSet.insert(corr.getSource());
          vevSet.insert(corr.getSink());}}}
 catch(const std::exception& xp){
    throw(std::invalid_argument(string("SetupCorrelator failed: ")+xp.what()));}
}


void MCObsGetHandler::setup_vevs(XMLHandler& xmls, 
                                const string& tagname,
                                set<OperatorInfo>& vevSet)
{
 try{
    list<XMLHandler> vevxml;
    vevxml=xmls.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       it=vevxml.begin();it!=vevxml.end();++it){
       OperatorInfo vev(*it);
       vevSet.insert(vev);}}
 catch(const std::exception& xp){
    throw(std::invalid_argument(string("SetupVEVs failed: ")+xp.what()));}
}


          // read correlation matrices

void MCObsGetHandler::setup_correlator_matrices(XMLHandler& xmlin,
                          const string& tagname,
                          bool hermitian, bool vevs, 
                          set<CorrelatorInfo>& corrSet,
                          set<OperatorInfo>& vevSet)
{
 try{
    list<XMLHandler> corrxml;
    corrxml=xmlin.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       it=corrxml.begin();it!=corrxml.end();++it){
       list<string> tagnames;
       tagnames.push_back("Operator");
       tagnames.push_back("OperatorString");
       tagnames.push_back("BLOperator");
       tagnames.push_back("BLOperatorString");
       list<XMLHandler> opxml=it->find_among_children(tagnames);
       set<OperatorInfo> opSet;
       for (list<XMLHandler>::iterator
          ot=opxml.begin();ot!=opxml.end();++ot)
             opSet.insert(OperatorInfo(*ot));
       for (set<OperatorInfo>::const_iterator
          ita=opSet.begin();ita!=opSet.end();ita++){
          if (vevs){ vevSet.insert(*ita);}
          for (set<OperatorInfo>::const_iterator
             itb=(hermitian)?ita:opSet.begin();itb!=opSet.end();itb++){
             CorrelatorInfo corr(*ita,*itb);
             corrSet.insert(corr);}}
       int namecount=it->count("AssignName");
       if (namecount>1) throw(std::invalid_argument("Multiple <AssignName> tags"));
       if (namecount==1){
          string name; xmlreadchild(*it,"AssignName",name);
          CorrelatorMatrixInfo corrkeep(opSet,hermitian,vevs);
          corrkeep.setName(name);}}}
 catch(const std::exception& xp){
    throw(std::invalid_argument(string("SetupCorrelatorMatrices failed: ")+xp.what()));}
}


                //  read a set of observables

void MCObsGetHandler::setup_obsset(XMLHandler& xmlin, const std::string& tagname, 
                                   std::set<MCObsInfo>& obsSet)
{
 try{
    list<XMLHandler> obsxml;
    obsxml=xmlin.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       ot=obsxml.begin();ot!=obsxml.end();++ot){
       list<XMLHandler> obsxml;
       obsxml=ot->find_among_children("MCObservable");
       for (list<XMLHandler>::iterator 
          it=obsxml.begin();it!=obsxml.end();++it){
          MCObsInfo obskey(*it);
          obsSet.insert(obskey);}}}
 catch(const std::exception& xp){
    throw(std::invalid_argument(string("SetupObsSet failed: ")+xp.what()));}
}


void MCObsGetHandler::read_bl_data(MCObsGetHandler::BasicLapHGetter& getter, Vector<Scalar>& result)
{
 uint nbins=m_bins_info.getNumberOfBins();
 uint rebin=m_bins_info.getRebinFactor();
 const set<unsigned int>& omit=m_bins_info.getOmissions();
 try{
   result.resize(nbins);
   if ((rebin==1)&&(omit.empty())){
      for (unsigned int k=0;k<nbins;k++){
         getter.getData(k,result[k]);}}
   else if ((rebin>1)&&(omit.empty())){
      if (!m_is_weighted){
         unsigned int count=0;
         double r=1.0/double(rebin);
         Scalar buffer,buffer2;
         for (unsigned int k=0;k<nbins;k++){
            getter.getData(count,buffer); ++count;
            for (unsigned int j=1;j<rebin;j++){
               getter.getData(count,buffer2); ++count;
               buffer+=buffer2;}
            result[k]=buffer*r;}}
      else{
         unsigned int count=0;
         Scalar buffer,buffer2;
         for (unsigned int k=0;k<nbins;k++){
            getter.getData(count,buffer);
            buffer*=m_BLorig_wts[count];
            double rd=m_BLorig_wts[count]; ++count;
            for (unsigned int j=1;j<rebin;j++){
               getter.getData(count,buffer2);
               buffer+=buffer2*m_BLorig_wts[count];
               rd+=m_BLorig_wts[count]; ++count;}
            result[k]=buffer/rd;}}}
   else if ((rebin==1)&&(!omit.empty())){
      set<unsigned int>::const_iterator om=omit.begin();
      unsigned int count=0;
      while ((om!=omit.end())&&(count==*om)){om++; ++count;}
      for (unsigned int k=0;k<nbins;k++){
         getter.getData(count,result[k]);
         ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}}
   else{
      if (!m_is_weighted){
         double r=1.0/double(rebin);
         Scalar buffer,buffer2;
         set<unsigned int>::const_iterator om=omit.begin();
         unsigned int count=0;
         while ((om!=omit.end())&&(count==*om)){om++; ++count;}
         for (unsigned int k=0;k<nbins;k++){
            getter.getData(count,buffer);
            ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
            for (unsigned int j=1;j<rebin;j++){
               getter.getData(count,buffer2);
               ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
               buffer+=buffer2;}
            result[k]=buffer*r;}}
      else{
         Scalar buffer,buffer2;
         set<unsigned int>::const_iterator om=omit.begin();
         unsigned int count=0;
         while ((om!=omit.end())&&(count==*om)){om++; ++count;}
         for (unsigned int k=0;k<nbins;k++){
            getter.getData(count,buffer);
            buffer*=m_BLorig_wts[count];
            double rd=m_BLorig_wts[count];
            ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
            for (unsigned int j=1;j<rebin;j++){
               getter.getData(count,buffer2);
               buffer2*=m_BLorig_wts[count];
               rd+=m_BLorig_wts[count];
               ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
               buffer+=buffer2;}
            result[k]=buffer/rd;}}}}
 catch(const std::exception& errmsg){
    //cout << errmsg <<endl;
    throw(std::invalid_argument(string("read_data failed  ")+string(errmsg.what())));}
}


bool MCObsGetHandler::query_bl_data(MCObsGetHandler::BasicLapHGetter& getter)
{
 uint nbins=m_bins_info.getNumberOfBins();
 uint rebin=m_bins_info.getRebinFactor();
 const set<unsigned int>& omit=m_bins_info.getOmissions();
 if ((rebin==1)&&(omit.empty())){
    for (unsigned int k=0;k<nbins;k++){
       if (!getter.queryData(k)) return false;}}
 else if ((rebin>1)&&(omit.empty())){
    unsigned int count=0;
    for (unsigned int k=0;k<nbins;k++){
       if (!getter.queryData(count++)) return false;
       for (unsigned int j=1;j<rebin;j++)
          if (!getter.queryData(count++)) return false;}}
 else if ((rebin==1)&&(!omit.empty())){
    set<unsigned int>::const_iterator om=omit.begin();
    unsigned int count=0;
    while ((om!=omit.end())&&(count==*om)){om++; ++count;}
    for (unsigned int k=0;k<nbins;k++){
       if (!getter.queryData(count)) return false;
       ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}}
 else{
    set<unsigned int>::const_iterator om=omit.begin();
    unsigned int count=0;
    while ((om!=omit.end())&&(count==*om)){om++; ++count;}
    for (unsigned int k=0;k<nbins;k++){
       if (!getter.queryData(count)) return false;
       ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
       for (unsigned int j=1;j<rebin;j++){
          if (!getter.queryData(count)) return false;
          ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}}}
 return true;
}


void MCObsGetHandler::get_sig_bin_data(const MCObsInfo& obsinfo, RVector& bins)
{
 if (m_binsdh==0)
    throw(std::invalid_argument(string("getBins fails due to unavailable bins for ")+obsinfo.str()));
 m_binsdh->getSymData(obsinfo,bins);
 uint nbins=m_bins_info.getNumberOfBins();
 if (bins.size()==nbins)
    return;
 uint rebin=bins.size()/nbins;
 if (!m_is_weighted){
   RVector temp(bins);
   bins.resize(nbins);
   unsigned int count=0;
   double r=1.0/double(rebin);
   for (unsigned int k=0;k<nbins;k++){
      double res=0.0;
      for (unsigned int j=0;j<rebin;j++){
         res+=temp[count]; ++count;}
      bins[k]=res*r;}}
 else{
   RVector temp(bins);
   bins.resize(nbins);
   vector<double> bwts;
   transform_weights(m_BLorig_wts,bwts,m_bins_info.getRebinFactor()/rebin,
                     m_bins_info.getOmissions());
   unsigned int count=0;
   for (unsigned int k=0;k<nbins;k++){
      double res=0.0;
      double rd=0.0;
      for (unsigned int j=0;j<rebin;j++){
         res+=bwts[count]*temp[count];
         rd+=bwts[count]; ++count;}
      bins[k]=res/rd;}}
}


#ifdef COMPLEXNUMBERS

    //  If numbers are stored as complex in the data files, then might as
    //  well read and store both the real and imaginary parts, since most
    //  likely, both will eventually be needed.


void MCObsGetHandler::get_bl_data(MCObsGetHandler::BasicLapHGetter& getter,
                                  RVector& results_re, RVector& results_im)
{
 try{
   Vector<Scalar> results;
   read_bl_data(getter,results);
   uint n=results.size();
   results_re.resize(n);
   results_im.resize(n);
   for (uint k=0;k<n;k++){
      results_re[k]=realpart(results[k]);
      results_im[k]=imaginarypart(results[k]);}}
 catch(const std::exception& errmsg){
    //cout << errmsg <<endl;
    throw(std::invalid_argument(string("get_data failed")+string(errmsg.what())));}
}


#else


void MCObsGetHandler::get_bl_data(MCObsGetHandler::BasicLapHGetter& getter,
                                  RVector& results)
{
 try{
   read_bl_data(getter,results);}
 catch(const std::exception& errmsg){
    //cout << errmsg <<endl;
    throw(std::invalid_argument(string("get_data failed: ")+string(errmsg.what())));}
}

#endif

// ****************************************************************

   //  This transforms the CLS weights by rebinning and omitting
   //  certain weights.  A rebinned weight is simply the sum of the
   //  original weights in the bin.

void MCObsGetHandler::transform_weights(const vector<double>& worig,
                                        vector<double>& wnew, uint rebin,
                                        const std::set<uint>& omit)
{
 uint nomit=omit.size();
 uint norig=worig.size();
 for (set<uint>::const_reverse_iterator rm=omit.rbegin();rm!=omit.rend();++rm){
    if (*rm>norig) --nomit;
    else break;}
 if (rebin==1){
    if (nomit==0){         // no rebinning, no omissions
       wnew=worig;
       return;}
    else{                  // no rebinning, omissions
       uint nbins=norig-nomit;
       set<unsigned int>::const_iterator om=omit.begin();
       unsigned int count=0;
       wnew.resize(nbins);
       while ((om!=omit.end())&&(count==*om)){om++; ++count;}
       for (unsigned int k=0;k<nbins;k++){
          wnew[k]=worig[count];
          ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}
       return;}}
 if (nomit==0){         // rebinning but no omissions
    uint nbins=norig/rebin;
    unsigned int count=0;
    wnew.resize(nbins);
    for (unsigned int k=0;k<nbins;k++){
       double r=worig[count]; ++count;
       for (unsigned int j=1;j<rebin;j++){
          r+=worig[count]; ++count;}
       wnew[k]=r;}
    return;}
 uint nbins=(norig-nomit)/rebin;   // rebinning and omissions
 set<unsigned int>::const_iterator om=omit.begin();
 unsigned int count=0;
 wnew.resize(nbins);
 while ((om!=omit.end())&&(count==*om)){om++; ++count;}
 for (unsigned int k=0;k<nbins;k++){
    double r=worig[count];
    ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
    for (unsigned int j=1;j<rebin;j++){
       r+=worig[count];
       ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}
     wnew[k]=r;}
}



void MCObsGetHandler::get_the_weights()
{
 m_bins_info.getMCEnsembleInfo().getWeights(m_BLorig_wts);
 transform_weights(m_BLorig_wts,m_wts,m_bins_info.getRebinFactor(),
                   m_bins_info.getOmissions());
}


// ***************************************************************************************
 
