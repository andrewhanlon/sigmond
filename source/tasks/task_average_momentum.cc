#include <string>
#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoAverageMomentum</Action>                                   *
// *     <CoefficientsProvided/>       (optional)                             *
// *     <CorrelatorMatrix>                                                   *
// *       <CorrelatorMatrixInfo>                                             *
// *         <BLOperator>...</BLOperator>                                     *
// *         <BLOperator>...</BLOperator>                                     *
// *              ...                                                         *
// *       </CorrelatorMatrixInfo>                                            *
// *       <Coefficients>+ - - + ....</Coefficients>  (optional)              *
// *     </CorrelatorMatrix>                                                  *
// *     <CorrelatorMatrix>                                                   *
// *          ...                                                             *
// *          ...                                                             *
// *     </CorrelatorMatrix>                                                  *
// *     <MinimumTimeSeparation>3</MinimumTimeSeparation>                     *
// *     <MaximumTimeSeparation>20</MaximumTimeSeparation>                    *
// *     <FileName>name_of_file</FileName>                                    *
// *     <FileMode>overwrite</FileMode>   (optional)                          *
// *    </Task>                                                               *
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>CompareCorrelators</Action>                                  *
// *     <CorrelatorMatrixInfo>                                               *
// *       <Operator>...</Operator>                                           *
// *       <Operator>...</Operator>                                           *
// *            ...                                                           *
// *     </CorrelatorMatrixInfo>                                              *
// *     <CompareID>1</CompareID>                                             *
// *     <MinimumTimeSeparation>3</MinimumTimeSeparation>                     *
// *     <MaximumTimeSeparation>20</MaximumTimeSeparation>                    *
// *    </Task>                                                               *
// *                                                                          *
// ****************************************************************************

/*
int relative_sign(double re_mean, double im_mean, double re_comp_mean, double im_comp_mean)
{
 bool real = (re_mean > 0.) ^ (re_comp_mean < 0.);
 bool imag = (im_mean > 0.) ^ (im_comp_mean < 0.);

 if (real != imag) {
   if (abs(re_comp_mean) > abs(im_comp_mean)) {
     imag = real;
   }
   else {
     real = imag;
   }
 }

 if (real) return 1;
 else return -1;
}

bool zero_warning(double re_mean, double im_mean, double re_err, double im_err)
{
 return ((abs(re_mean)<=4.0*re_err)&&(abs(im_mean)<=4.0*im_err));
}

vector<double> get_coefs(MCObsHandler *m_obs, CorrelatorAtTimeInfo& corrt, vector<CorrelatorAtTimeInfo>& toAverage,
                         XMLHandler& xmlout)
{
 XMLHandler xml_ave, xml_corr_first;
 xml_ave.set_root("Average");
 xml_corr_first.set_root("Correlator");
 XMLHandler xml_corr_first_out;
 corrt.output(xml_corr_first_out);
 xml_corr_first.put_child(xml_corr_first_out);

 double coef = 1./(1.+toAverage.size());
 vector<double> coefs;
 // Get comparison values
 MCObsInfo obskeyRe(corrt,RealPart);
 MCObsInfo obskeyIm(corrt,ImaginaryPart);
 if ((!m_obs->queryBins(obskeyRe))||(!m_obs->queryBins(obskeyIm)))
    throw(string("Data was here just a minute ago..."));
 double re_mean=0.0; double re_err=1.0;
 double im_mean=0.0; double im_err=1.0;
 MCEstimate est;
 est=m_obs->getJackknifeEstimate(obskeyRe);
 re_mean=est.getFullEstimate();
 re_err=est.getSymmetricError();
 est=m_obs->getJackknifeEstimate(obskeyIm);
 im_mean=est.getFullEstimate();
 im_err=est.getSymmetricError();

 xml_corr_first.put_child("Value", "("+to_string(re_mean)+"+/-"+to_string(re_err)+","+to_string(im_mean)+"+/-"+to_string(im_err)+")");
 xml_corr_first.put_child("Coefficient", to_string(coef));
 
 if (zero_warning(re_mean,im_mean,re_err,im_err)) xml_corr_first.put_child("Zero");
 xml_ave.put_child(xml_corr_first);
 

 //cout << "first (Re): " << re_mean << "+/-" << re_err << endl;
 //cout << "first (Im): " << im_mean << "+/-" << im_err << endl;

 double re_comp_mean=0.0; double re_comp_err=1.0;
 double im_comp_mean=0.0; double im_comp_err=1.0;
 vector<CorrelatorAtTimeInfo>::iterator corr_it=toAverage.begin();
 while(corr_it!=toAverage.end()) {
    XMLHandler xml_corr;
    xml_corr.set_root("Correlator");
    XMLHandler xml_corr_out;
    corr_it->output(xml_corr_out);
    xml_corr.put_child(xml_corr_out);
    MCObsInfo obskeyCompRe(*corr_it,RealPart);
    MCObsInfo obskeyCompIm(*corr_it,ImaginaryPart);
    if ((!m_obs->queryBins(obskeyCompRe))||(!m_obs->queryBins(obskeyCompIm))) {
       XMLHandler xml_err;
       xml_err.set_root("Error");
       xml_err.put_child("Found",corrt.output());
       xml_err.put_child("NotFound",corr_it->output());
       xmlout.put_child(xml_err);
       throw(string("Data not found in same correlator"));}
    est=m_obs->getJackknifeEstimate(obskeyCompRe);
    re_comp_mean=est.getFullEstimate();
    re_comp_err=est.getSymmetricError();
    est=m_obs->getJackknifeEstimate(obskeyCompIm);
    im_comp_mean=est.getFullEstimate();
    im_comp_err=est.getSymmetricError();

    xml_corr.put_child("Value", "("+to_string(re_comp_mean)+"+/-"+to_string(re_comp_err)+","+to_string(im_comp_mean)+"+/-"+to_string(im_comp_err)+")");
    xml_corr.put_child("Coefficient",to_string(coef*relative_sign(re_mean,im_mean,re_comp_mean,im_comp_mean)));
    if (zero_warning(re_comp_mean,im_comp_mean,re_comp_err,im_comp_err)) xml_corr.put_child("Zero");
    xml_ave.put_child(xml_corr);

    coefs.push_back(coef*relative_sign(re_mean,im_mean,re_comp_mean,im_comp_mean));
    ++corr_it;
 }

 xmlout.put_child(xml_ave);

 coefs.push_back(coef);
 toAverage.push_back(corrt);

 return coefs;
}

void find_coefficients(CorrelatorMatrixInfo& first_corrMat, vector<CorrelatorMatrixInfo>& corrMatInfos,
                       map<OperatorInfo,char>& coefs)
{
  const set<OperatorInfo>& first_ops = first_corrMat.getOperators();
  for (set<OperatorInfo>::const_iterator first_op=first_ops.begin(); first_op!=first_ops.end(); ++first_op) {
    coefs.insert(pair<OperatorInfo,char>(*first_op,'+'));
  }


}
*/

void determine_coefficients(MCObsHandler* m_obs, map<OperatorInfo,int>& coefs, vector<CorrelatorMatrixInfo>& corrmats)
{
}

void store_in_memory(MCObsHandler *m_obs, CorrelatorAtTimeInfo& corrt_result,
                     vector<CorrelatorAtTimeInfo>& to_average, map<OperatorInfo,int>& coefs_map,
                     uint minTime, uint maxTime, set<MCObsInfo>& obskeys, XMLHandler& xmlout)
{
  XMLHandler xml_corr;
  xml_corr.set_root("Correlator");

  XMLHandler xml_ave;
  xml_ave.set_root("Average");
  vector<double> coefs;
  for (vector<CorrelatorAtTimeInfo>::iterator corrt=to_average.begin();
       corrt!=to_average.end(); ++corrt) {
    double coef = double(coefs_map[corrt->getSource()]*coefs_map[corrt->getSink()]); // to_average.size();
    coefs.push_back(coef);
    
    XMLHandler xml_corr_out;
    corrt->output(xml_corr_out);
    xml_ave.put_child(xml_corr_out);
    xml_ave.put_child("Coefficient",to_string(coef));
    
    MCObsInfo obskeyRe(*corrt,RealPart);
    MCObsInfo obskeyIm(*corrt,ImaginaryPart);
    if ((!m_obs->queryBins(obskeyRe))||(!m_obs->queryBins(obskeyIm)))
      throw(string("Data was here just a minute ago..."));
    double re_mean=0.0; double re_err=1.0;
    double im_mean=0.0; double im_err=1.0;
    MCEstimate est;
    est=m_obs->getEstimate(obskeyRe);
    re_mean=est.getAverageEstimate();
    re_err=est.getSymmetricError();
    est=m_obs->getEstimate(obskeyIm);
    im_mean=est.getAverageEstimate();
    im_err=est.getSymmetricError();
    xml_ave.put_child("Value", "("+to_string(re_mean)+"+/-"+to_string(re_err)+","+to_string(im_mean)+"+/-"+to_string(im_err)+")");
  }
  xml_corr.put_child(xml_ave);

  vector<ComplexArg> args(2);
  args[0]=RealPart; args[1]=ImaginaryPart;
  for (uint t=minTime; t<=maxTime; ++t) {
    for (uint kk=0; kk<args.size();kk++) {
      vector<const Vector<double>* > bins(to_average.size());
      vector<const Vector<double>* >::iterator bins_it;
      bins_it = bins.begin();
      for (vector<CorrelatorAtTimeInfo>::iterator corrt=to_average.begin();
           corrt!=to_average.end(); ++corrt) {
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
    if (t==minTime) {
      XMLHandler xml_result;
      xml_result.set_root("Result");
      XMLHandler xml_corr_out;
      corrt_result.output(xml_corr_out);
      xml_result.put_child(xml_corr_out);
      MCObsInfo obskeyRe(corrt_result,RealPart);
      MCObsInfo obskeyIm(corrt_result,ImaginaryPart);
      if ((!m_obs->queryBins(obskeyRe))||(!m_obs->queryBins(obskeyIm)))
        throw(string("Data was here just a minute ago..."));
      double re_mean=0.0; double re_err=1.0;
      double im_mean=0.0; double im_err=1.0;
      MCEstimate est;
      est=m_obs->getEstimate(obskeyRe);
      re_mean=est.getAverageEstimate();
      re_err=est.getSymmetricError();
      est=m_obs->getEstimate(obskeyIm);
      im_mean=est.getAverageEstimate();
      im_err=est.getSymmetricError();
      xml_result.put_child("Value", "("+to_string(re_mean)+"+/-"+to_string(re_err)+","+to_string(im_mean)+"+/-"+to_string(im_err)+")");
      xml_corr.put_child(xml_result);
   }
 }
 xmlout.put_child(xml_corr);
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
      if (fmode=="overwrite") overwrite=true;
    }
    uint minTime, maxTime;
    xmlread(xmltask,"MinimumTimeSeparation",minTime,"DoAverageMomentum");
    xmlread(xmltask,"MaximumTimeSeparation",maxTime,"DoAverageMomentum");
    list<XMLHandler> corrmatsXML=xmltask.find_among_children("CorrelatorMatrix");
    bool coefs_provided = false;
    if (xml_tag_count(xmltask,"CoefficientsProvided")==1) coefs_provided = true;
    
    // Insert the correlator matrices into corrmats
    vector<CorrelatorMatrixInfo> corrmats;
    map<OperatorInfo,int> coefs;
    list<string> tagnames;
    tagnames.push_back("Operator");
    tagnames.push_back("OperatorString");
    tagnames.push_back("BLOperator");
    tagnames.push_back("BLOperatorString");
    for (list<XMLHandler>::iterator corrmatXML=corrmatsXML.begin(); corrmatXML!=corrmatsXML.end(); ++corrmatXML) {
      XMLHandler corrmatXMLhandler(*corrmatXML,"CorrelatorMatrixInfo");
      if (corrmatXML==corrmatsXML.begin()) corrmatXMLhandler.put_child("HermitianMatrix");
      else {
        corrmatXMLhandler.seek_unique("HermitianMatrix");
        if (corrmatXMLhandler.good()) corrmatXMLhandler.erase_current_element();
        corrmatXMLhandler.seek_root();
      }

      CorrelatorMatrixInfo corrmat(corrmatXMLhandler);
      corrmats.push_back(corrmat);
      if (coefs_provided) {
        string coefs_input;
        xmlread(*corrmatXML,"Coefficients",coefs_input,"DoAverageMomentum");
        stringstream coefs_stream(coefs_input);

        list<XMLHandler> opsXML=corrmatXMLhandler.find_among_children(tagnames);
        for (list<XMLHandler>::iterator xml_it=opsXML.begin(); xml_it!=opsXML.end(); xml_it++) {
          string sign;
          if (getline(coefs_stream, sign, ' ')) {
            coefs.insert(pair<OperatorInfo,int>(OperatorInfo(*xml_it), stoi(sign)));
          }
        }
      }
    }

    if (!coefs_provided) determine_coefficients(m_obs,coefs,corrmats);
    
    set<MCObsInfo> obskeys;

    vector<CorrelatorMatrixInfo>::iterator corrmat_it=corrmats.begin();
    CorrelatorMatrixInfo first_corrmat = *corrmat_it;
    for (first_corrmat.begin(); !first_corrmat.end(); ++first_corrmat) {
      CorrelatorInfo corr=first_corrmat.getCurrentCorrelatorInfo();
      CorrelatorAtTimeInfo corrt(corr,minTime,true,false);
      vector<CorrelatorAtTimeInfo> to_average;
      to_average.push_back(corrt);
      for (corrmat_it++; corrmat_it!=corrmats.end(); ++corrmat_it) {
        CorrelatorMatrixInfo corrmat_comp = *corrmat_it;
        for (corrmat_comp.begin(); !corrmat_comp.end(); ++corrmat_comp) {
          CorrelatorInfo corr_compare = corrmat_comp.getCurrentCorrelatorInfo();
          CorrelatorAtTimeInfo corrt_compare(corr_compare,minTime,true,false);
          if (corr.rotationallyEquivalent(corr_compare)) to_average.push_back(corrt_compare);
        }
      }
      corrmat_it = corrmats.begin();
      
      // Make CorrelatorAtTimeInfo from GenIrrepOperatorInfo's
      BasicLapHOperatorInfo srcOp = corr.getSource().getBasicLapH();
      BasicLapHOperatorInfo snkOp = corr.getSink().getBasicLapH();
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
        srcOpString += srcOp.getFlavor(hadron).substr(0,2);
        if (srcOp.getNumberOfHadrons() > 1) srcOpString += "_" + srcOp.getLGIrrep(hadron);
        srcOpString += "_" + srcOp.getSpatialType(hadron) + "_" + to_string(srcOp.getSpatialIdNumber(hadron));
        if (hadron < srcOp.getNumberOfHadrons()) srcOpString += "_";
      }
      for (uint hadron=1; hadron<=snkOp.getNumberOfHadrons(); hadron++) {
        snkOpString += snkOp.getFlavor(hadron).substr(0,2);
        if (snkOp.getNumberOfHadrons() > 1) snkOpString += "_" + snkOp.getLGIrrep(hadron);
        snkOpString += "_" + snkOp.getSpatialType(hadron) + "_" + to_string(snkOp.getSpatialIdNumber(hadron));
        if (hadron < snkOp.getNumberOfHadrons()) snkOpString += "_";
      }

      CorrelatorInfo corr_result(OperatorInfo(snkOpString,OperatorInfo::GenIrrep),OperatorInfo(srcOpString,OperatorInfo::GenIrrep));
      CorrelatorAtTimeInfo corrt_result(corr_result,minTime,true,false);
      store_in_memory(m_obs,corrt_result,to_average,coefs,minTime,maxTime,obskeys,xmlout);
    }
    
    XMLHandler xmlf;
    m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
    xmlout.put_child(xmlf);
  }
  catch(const std::exception& errmsg){
    throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
         +string(errmsg.what())).c_str()));
  }
}

bool compare_bins(const Vector<double>* bins, const Vector<double>* bins_compare)
{
 double epsilon=1.5e-07;
 double diff;

 
 for (uint n=0; n<bins->size();++n) {
   //cout << 3.*(*bins)[n] << " = " << (*bins_compare)[n] << endl;
   diff = abs((*bins)[n]-(*bins_compare)[n]);
   if (diff > epsilon) return false;
 }
 return true;
}

bool compare_correlators(MCObsHandler *m_obs, CorrelatorInfo& corr, CorrelatorInfo& corr_compare, uint minTime, uint maxTime)
{
 vector<ComplexArg> args(2);
 args[0]=RealPart; args[1]=ImaginaryPart;
 for (uint t=minTime; t<=maxTime; ++t) {
   CorrelatorAtTimeInfo corrt(corr,t,true,false);
   CorrelatorAtTimeInfo corrt_compare(corr_compare,t,true,false);
   for (uint kk=0; kk<args.size(); kk++) {
     const Vector<double>* bins=&(m_obs->getBins(MCObsInfo(corrt,args[kk])));
     const Vector<double>* bins_compare=&(m_obs->getBins(MCObsInfo(corrt_compare,args[kk])));
     if (!compare_bins(bins,bins_compare)) return false;
   }
 }
 return true;
}

void TaskHandler::compareCorrelators(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 try {
   // Read in XML
   xmlout.set_root("CompareCorrelators");
   uint minTime, maxTime;
   xmlread(xmltask,"MinimumTimeSeparation",minTime,"CompareCorrelators");
   xmlread(xmltask,"MaximumTimeSeparation",maxTime,"CompareCorrelators");
   uint compareID;
   xmlread(xmltask,"CompareID",compareID,"CompareCorrelators");
   XMLHandler xmlcm(xmltask,"CorrelatorMatrixInfo");
   CorrelatorMatrixInfo cormat(xmlcm);
   for (cormat.begin(); !cormat.end(); ++cormat) {
     CorrelatorInfo corr=cormat.getCurrentCorrelatorInfo();
     string snk = "iso" + corr.getSink().getGenIrrep().short_output();
     string src = "iso" + corr.getSource().getGenIrrep().short_output();
     snk.replace(snk.find_last_of("0"),1,"1");
     src.replace(src.find_last_of("0"),1,"1");
     //cout << "snk: " << snk << endl;
     //cout << "src: " << src << endl;
     CorrelatorInfo corr_compare(OperatorInfo(snk,OperatorInfo::GenIrrep),OperatorInfo(src,OperatorInfo::GenIrrep));
     if (!compare_correlators(m_obs,corr,corr_compare,minTime,maxTime)) {
       cout << "Correlators not the same!" << endl;
       exit(1);
     }
   }
 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task CompareCorrelators: ")
        +string(errmsg.what())).c_str()));
 }
}
