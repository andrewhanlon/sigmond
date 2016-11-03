#include <string>
#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoAverageMomentum</Action>                                   *
// *     <CorrelatorsToAverage>                                               *
// *       <CorrelatorMatrixInfo>                                             *
// *         <BLOperator>...</BLOperator>                                     *
// *         <BLOperator>...</BLOperator>                                     *
// *              ...                                                         *
// *       </CorrelatorMatrixInfo>                                            *
// *       <CorrelatorMatrixInfo>                                             *
// *            ...                                                           *
// *       </CorrelatorMatrixInfo>                                            *
// *          ...                                                             *
// *     </CorrelatorsToAverage>                                              *
// *     <MinimumTimeSeparation>3</MinimumTimeSeparation>                     *
// *     <MaximumTimeSeparation>20</MaximumTimeSeparation>                    *
// *     <FileName>name_of_file</FileName>                                    *
// *     <FileMode>overwrite</FileMode>   (optional)                          *
// *    </Task>                                                               *
// *                                                                          *
// ****************************************************************************

void store_in_memory(MCObsHandler *m_obs, CorrelatorAtTimeInfo& corrt_result,
                     vector<CorrelatorAtTimeInfo>& toAverage, vector<double>& coefs,
                     uint minTime, uint maxTime, set<MCObsInfo>& obskeys, XMLHandler& xmlout)
{
 vector<ComplexArg> args(2);
 args[0]=RealPart; args[1]=ImaginaryPart;
 for (uint t=minTime; t<=maxTime; ++t) {
    for (uint kk=0; kk<args.size();kk++) {
      vector<const Vector<double>* > bins(toAverage.size());
      vector<const Vector<double>* >::iterator bins_it;
      bins_it = bins.begin();
      for (vector<CorrelatorAtTimeInfo>::iterator corrt=toAverage.begin();
           corrt!=toAverage.end(); ++corrt) {
        corrt->resetTimeSeparation(t);
        *bins_it=&(m_obs->getBins(MCObsInfo(*corrt,args[kk])));
        ++bins_it;
      }
      int nbins=bins[0]->size();
      Vector<double> result(nbins);
      for (int bin=0; bin<nbins; bin++) {
        double temp=0.0;
        bins_it = bins.begin();
        for (vector<double>::iterator coefs_it = coefs.begin(); coefs_it != coefs.end(); ++coefs_it) {
          temp+=(*coefs_it)*(*(*bins_it))[bin];
          ++bins_it;
        }
        result[bin] = temp;
      }
      corrt_result.resetTimeSeparation(t);
      MCObsInfo averaged_corrt(corrt_result, args[kk]);
      m_obs->putBins(averaged_corrt,result);
      obskeys.insert(averaged_corrt);
    }
 }
}

int relative_sign(double re_mean, double im_mean, double re_comp_mean, double im_comp_mean)
{
 bool real = (re_mean > 0.) ^ (re_comp_mean < 0.);
 bool imag = (im_mean > 0.) ^ (im_comp_mean < 0.);

 if (real != imag) {
   cerr << "Relative Sign Issue" << endl;
   exit(1);
 }

 if (real) return 1;
 else return -1;
}

void zero_warning(double re_mean, double im_mean, double re_err, double im_err)
{
 if ((abs(re_mean)<=4.0*re_err)&&(abs(im_mean)<=4.0*im_err)) {
   cout << "Zero Warning" << endl;
 }
}

vector<double> get_coefs(MCObsHandler *m_obs, CorrelatorAtTimeInfo& corrt, vector<CorrelatorAtTimeInfo>& toAverage,
                         XMLHandler& xmlout)
{
 double coef = 1./(1.+toAverage.size());
 vector<double> coefs;
 // Get comparison values
 MCObsInfo obskeyRe(corrt,RealPart);
 MCObsInfo obskeyIm(corrt,ImaginaryPart);
 if ((!m_obs->queryBins(obskeyRe))&&(!m_obs->queryBins(obskeyIm))) {
    cerr << "Data not found...weird...aborting" << endl;
    exit(1);}
 double re_mean=0.0; double re_err=1.0;
 double im_mean=0.0; double im_err=1.0;
 if (m_obs->queryBins(obskeyRe)){
    MCEstimate est=m_obs->getJackknifeEstimate(obskeyRe);
    re_mean=std::abs(est.getFullEstimate());
    re_err=est.getSymmetricError();}
 if (m_obs->queryBins(obskeyIm)){
    MCEstimate est=m_obs->getJackknifeEstimate(obskeyIm);
    im_mean=std::abs(est.getFullEstimate());
    im_err=est.getSymmetricError();}
 
 zero_warning(re_mean,im_mean,re_err,im_err);
 

 //cout << "first (Re): " << re_mean << "+/-" << re_err << endl;
 //cout << "first (Im): " << im_mean << "+/-" << im_err << endl;

 double re_comp_mean=0.0; double re_comp_err=1.0;
 double im_comp_mean=0.0; double im_comp_err=1.0;
 vector<CorrelatorAtTimeInfo>::iterator corr_it=toAverage.begin();
 while(corr_it!=toAverage.end()) {
    MCObsInfo obskeyCompRe(*corr_it,RealPart);
    MCObsInfo obskeyCompIm(*corr_it,ImaginaryPart);
    if ((!m_obs->queryBins(obskeyCompRe))&&(!m_obs->queryBins(obskeyCompIm))) {
       cout << "Data not found in same correlator" << endl;
       corr_it = toAverage.erase(corr_it);
       continue;}
    if (m_obs->queryBins(obskeyCompRe)){
      MCEstimate est=m_obs->getJackknifeEstimate(obskeyCompRe);
      re_comp_mean=std::abs(est.getFullEstimate());
      re_comp_err=est.getSymmetricError();}
    if (m_obs->queryBins(obskeyCompIm)){
      MCEstimate est=m_obs->getJackknifeEstimate(obskeyCompIm);
      im_comp_mean=std::abs(est.getFullEstimate());
      im_comp_err=est.getSymmetricError();}

    zero_warning(re_comp_mean,im_comp_mean,re_comp_err,im_comp_err);

    coefs.push_back(coef*relative_sign(re_mean,im_mean,re_comp_mean,im_comp_mean));
    ++corr_it;
 }

 coefs.push_back(coef);
 toAverage.push_back(corrt);

 return coefs;
}

void TaskHandler::doAverageMomentum(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 // This map is used, because BasicLapHOperatorInfo::getIsospin() returns the flavor
 // if numHadrons is one. But, the GenIrrepOperatorInfop requires the isospin
 map<string,string> isospinMap;
 isospinMap["eta"]="singlet";
 isospinMap["phi"]="singlet";
 isospinMap["kaon"]="doublet";
 isospinMap["kbar"]="doublet";
 isospinMap["pion"]="triplet";
 isospinMap["lambda"]="singlet";
 isospinMap["omega"]="singlet";
 isospinMap["nuclon"]="doublet";
 isospinMap["xi"]="doublet";
 isospinMap["sigma"]="triplet";
 isospinMap["delta"]="quartet";

 try{
   // Read in XML
   xmlout.set_root("DoAverageMomentum");
   string filename;
   xmlread(xmltask,"FileName",filename,"DoAverageMomentum");
   bool overwrite = false;  // protect mode
   if (xml_tag_count(xmltask,"FileMode")==1){
     string fmode;
     xmlread(xmltask,"FileMode",fmode,"FileListInfo");
     fmode=tidyString(fmode);
     if (fmode=="overwrite") overwrite=true;}
   uint minTime, maxTime;
   xmlread(xmltask,"MinimumTimeSeparation",minTime,"DoAverageMomentum");
   xmlread(xmltask,"MaximumTimeSeparation",maxTime,"DoAverageMomentum");
   XMLHandler xmlf(xmltask,"CorrelatorsToAverage");
   list<string> tagnames;
   tagnames.push_back("CorrelatorMatrixInfo");
   list<XMLHandler> corrMatxml=xmlf.find_among_children(tagnames);

   // Put first CorrelatorMatrixInfo in first_corrMat and build opIdMap
   list<XMLHandler>::iterator ct=corrMatxml.begin();
   CorrelatorMatrixInfo first_corrMat(*ct);
   map<BasicLapHOperatorInfo,string> opIdMap;
   {
     map<string,int> flavorMap;
     const set<OperatorInfo>& opList = first_corrMat.getOperators();
     for (set<OperatorInfo>::const_iterator op_it=opList.begin(); op_it!=opList.end(); ++op_it) {
       BasicLapHOperatorInfo blOp = op_it->getBasicLapH();
       string flavor = "";
       for (uint hadron=1; hadron<=blOp.getNumberOfHadrons(); hadron++) {
         flavor += blOp.getFlavor(hadron);
         if (hadron < blOp.getNumberOfHadrons()) flavor += "_";
       }
       if (flavorMap.find(flavor) == flavorMap.end()) {
         flavorMap.insert(pair<string,int>(flavor,1));
       }
       else {
         flavorMap[flavor]++;
       }
       opIdMap.insert(pair<BasicLapHOperatorInfo,string>(blOp,to_string(flavorMap[flavor])));
     }
   }

   // Put the rest of the CorrelatorMatrixInfo's into corrMatInfos
   vector<CorrelatorMatrixInfo> corrMatInfos;
   for (++ct; ct!=corrMatxml.end();++ct) {
     corrMatInfos.push_back(CorrelatorMatrixInfo(*ct));
   }

   // Begin looping over the correlator elements of first_corrMat 
   set<MCObsInfo> obskeys;
   for (first_corrMat.begin(); !first_corrMat.end(); ++first_corrMat) {
     CorrelatorInfo corr=first_corrMat.getCurrentCorrelatorInfo();
     CorrelatorAtTimeInfo corrt(corr,minTime,false,false);
     // check the data exists
     MCObsInfo obskeyRe(corrt,RealPart);
     MCObsInfo obskeyIm(corrt,ImaginaryPart);
     if ((!m_obs->queryBins(obskeyRe))||(!m_obs->queryBins(obskeyIm))) {
       cout << "Warning: data not found" << endl;
       continue;}
     vector<CorrelatorAtTimeInfo> toAverage;
     for (vector<CorrelatorMatrixInfo>::iterator corrMat_it=corrMatInfos.begin();
          corrMat_it!=corrMatInfos.end(); ++corrMat_it) {
       uint num_compare = 0;
       for (corrMat_it->begin(); !corrMat_it->end(); ++(*corrMat_it)) {
         CorrelatorInfo corr_compare = corrMat_it->getCurrentCorrelatorInfo();
         if (corr.rotationallyEquivalent(corr_compare)) {
           num_compare++;
           if (num_compare==1) {
             CorrelatorAtTimeInfo corrt_compare(corr_compare,minTime,false,false);
             toAverage.push_back(corrt_compare);
           }
         }
       }
       if (num_compare==0) {
         cerr << "Couldn't Find a matching correlator entry" << endl;
         exit(1);
       }
       if (num_compare>1) cout << "Not 1-to-1 Warning" << endl;
     }
  
     vector<double> coefs = get_coefs(m_obs,corrt,toAverage,xmlout);
     // Make CorrelatorAtTimeInfo from GenIrrepOperatorInfo's
     BasicLapHOperatorInfo srcOp = corrt.getSource().getBasicLapH();
     BasicLapHOperatorInfo snkOp = corrt.getSink().getBasicLapH();
     string srcIsospin = srcOp.getIsospin();
     string snkIsospin = snkOp.getIsospin();
     if (srcOp.getNumberOfHadrons()==1) srcIsospin=isospinMap[srcIsospin];
     if (snkOp.getNumberOfHadrons()==1) snkIsospin=isospinMap[snkIsospin];
     string srcOpString = "iso"+srcIsospin+" P=("+to_string(srcOp.getXMomentum())
                        + ","+to_string(srcOp.getYMomentum())+","+to_string(srcOp.getZMomentum())+") "
                        + srcOp.getLGIrrep() + "_" + to_string(srcOp.getLGIrrepRow()) + " ";
     string snkOpString = "iso"+snkIsospin+" P=("+to_string(snkOp.getXMomentum())
                        + ","+to_string(snkOp.getYMomentum())+","+to_string(snkOp.getZMomentum())+") "
                        + snkOp.getLGIrrep() + "_" + to_string(snkOp.getLGIrrepRow()) + " ";

     for (uint hadron=1; hadron<=srcOp.getNumberOfHadrons(); hadron++) {
       srcOpString += srcOp.getFlavor(hadron);
       if (hadron < srcOp.getNumberOfHadrons()) srcOpString += "_";
     }
     for (uint hadron=1; hadron<=snkOp.getNumberOfHadrons(); hadron++) {
       snkOpString += snkOp.getFlavor(hadron);
       if (hadron < snkOp.getNumberOfHadrons()) snkOpString += "_";
     }
     srcOpString += " " + opIdMap[srcOp];
     snkOpString += " " + opIdMap[snkOp];

     CorrelatorInfo corr_result(OperatorInfo(snkOpString,OperatorInfo::GenIrrep),OperatorInfo(srcOpString,OperatorInfo::GenIrrep));
     CorrelatorAtTimeInfo corrt_result(corr_result,minTime,false,false);
     store_in_memory(m_obs,corrt_result,toAverage,coefs,minTime,maxTime,obskeys,xmlout);
   }

   XMLHandler xmlff;
   m_obs->writeSamplingValuesToFile(obskeys,filename,xmlff,overwrite);
   xmlout.put_child(xmlff);
 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
