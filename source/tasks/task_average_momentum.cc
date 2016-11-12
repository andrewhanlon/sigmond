#include <string>
#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// ****************************************************************************
// *                                                                          *
// *    <Task>                                                                *
// *     <Action>DoAverageMomentum</Action>                                   *
// *     <CoefficientsProvided/>       (optional)                             *
// *     <MultiCompare/>           (optional)                                 *
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


int relative_sign(double mean, double comp_mean, double err, double comp_err)
{
  if ((abs(mean)<=4.0*err)||(abs(comp_mean)<=4.0*comp_err)) return 0;
  
  bool sign = (mean > 0.) ^ (comp_mean < 0.);

  if (sign) return 1;
  else return -1;
}

void determine_coefficients(MCObsHandler* m_obs, map<OperatorInfo,int>& coefs, vector<CorrelatorMatrixInfo>& corrmats,
                            int minTime, bool multi_compare, XMLHandler& xmlout)
{
  map<OperatorInfo,uint> coef_info;
  // Set all the coefficients for operators in the first correlation matrix to one
  vector<CorrelatorMatrixInfo>::iterator corrmat_it=corrmats.begin();
  const set<OperatorInfo> first_ops=corrmat_it->getOperators();
  for (set<OperatorInfo>::const_iterator first_op=first_ops.begin(); first_op!=first_ops.end(); ++first_op) {
    coefs.insert(pair<OperatorInfo,int>(*first_op,1));
  }
  // Set all the coefficients for single-particle operators in the other
  // Correlation matrices to one, and all other operators to zero
  for (++corrmat_it; corrmat_it!=corrmats.end(); ++corrmat_it) {
    const set<OperatorInfo> ops=corrmat_it->getOperators();
    for (set<OperatorInfo>::const_iterator op=ops.begin(); op!=ops.end(); ++op) {
      if ((op->isBasicLapH())&&(op->getBasicLapH().getNumberOfHadrons() == 1)) {
        coefs.insert(pair<OperatorInfo,int>(*op,1));
      }
      else {
        coefs.insert(pair<OperatorInfo,int>(*op,0));
        coef_info.insert(pair<OperatorInfo,uint>(*op,0));
      }
    }
  }
  
  double re_mean=0.0; double re_err=1.0;
  double im_mean=0.0; double im_err=1.0;
  double re_comp_mean=0.0; double re_comp_err=1.0;
  double im_comp_mean=0.0; double im_comp_err=1.0;
  MCEstimate est;

  // Begin looping over the two-particle operators in the first correlator matrix.
  // We will then use these to form a correlator. We will then find the corresponding
  // correlators in the other correlation matrices. It is thus the two-particle
  // operator (the snk) in the other correlation matrix that will have its coefficient affected.
  for (set<OperatorInfo>::const_iterator snk_op=first_ops.begin(); snk_op!=first_ops.end(); ++snk_op) {
    if ((snk_op->isBasicLapH())&&(snk_op->getBasicLapH().getNumberOfHadrons() == 1)) continue;
    for (set<OperatorInfo>::const_iterator src_op=first_ops.begin(); src_op!=first_ops.end(); ++src_op) {
      if (snk_op==src_op) continue;
      if ((!multi_compare)&&(src_op->isBasicLapH())&&(src_op->getBasicLapH().getNumberOfHadrons() != 1)) continue;
      CorrelatorInfo first_corr(*snk_op,*src_op);
      CorrelatorAtTimeInfo first_corrt(first_corr,minTime,true,false);
      
      MCObsInfo obskeyRe(first_corrt,RealPart);
      MCObsInfo obskeyIm(first_corrt,ImaginaryPart);
      est=m_obs->getEstimate(obskeyRe);
      re_mean=est.getAverageEstimate();
      re_err=est.getSymmetricError();
      est=m_obs->getEstimate(obskeyIm);
      im_mean=est.getAverageEstimate();
      im_err=est.getSymmetricError();

      corrmat_it=corrmats.begin();
      for (++corrmat_it; corrmat_it!=corrmats.end(); ++corrmat_it) {
        CorrelatorMatrixInfo corrmat_comp = *corrmat_it;
        for (corrmat_comp.begin(); !corrmat_comp.end(); ++corrmat_comp) {
          CorrelatorInfo corr_comp = corrmat_comp.getCurrentCorrelatorInfo();
          CorrelatorAtTimeInfo corrt_comp(corr_comp,minTime,true,false);
          if (first_corr.rotationallyEquivalent(corr_comp)) {
            MCObsInfo obskeyReComp(corrt_comp,RealPart);
            MCObsInfo obskeyImComp(corrt_comp,ImaginaryPart);
            est=m_obs->getEstimate(obskeyReComp);
            re_comp_mean=est.getAverageEstimate();
            re_comp_err=est.getSymmetricError();
            est=m_obs->getEstimate(obskeyImComp);
            im_comp_mean=est.getAverageEstimate();
            im_comp_err=est.getSymmetricError();

            // We are only concerned with the sink op
            OperatorInfo snk_comp = corr_comp.getSink();
            int re_rel_sign = relative_sign(re_mean,re_comp_mean,re_err,re_comp_err);
            int im_rel_sign = relative_sign(im_mean,im_comp_mean,im_err,im_comp_err);
            coefs[snk_comp] += re_rel_sign + im_rel_sign;

            if (re_rel_sign == 1 || im_rel_sign == 1) coef_info[snk_comp] |= 0x1u;
            if (re_rel_sign == -1 || im_rel_sign == -1) coef_info[snk_comp] |= 0x2u;
          }
        }
      }
    }
  }

  XMLHandler coefs_out;
  coefs_out.set_root("Coefficients");

  for (map<OperatorInfo,uint>::iterator coef_it=coef_info.begin(); coef_it!=coef_info.end(); ++coef_it) {
    XMLHandler op_xml;
    op_xml.set_root("Operator");
    XMLHandler xml_corr_out;
    coef_it->first.output(xml_corr_out);
    op_xml.put_child(xml_corr_out);
    if (coefs[coef_it->first] == 0) {
      op_xml.put_child("Zero");
      coefs[coef_it->first] = 1;
    }
    if (coef_it->second == 0) {
      op_xml.put_child("Info", "Null");
    }
    else if (coef_it->second == 1) {
      op_xml.put_child("Info", "Increment");
    }
    else if (coef_it->second == 2) {
      op_xml.put_child("Info", "Decrement");
    }
    else if (coef_it->second == 3) {
      op_xml.put_child("Info", "Mixed");
    }
    coefs[coef_it->first] /= abs(coefs[coef_it->first]);
    op_xml.put_child("Coefficient",to_string(coefs[coef_it->first]));
    coefs_out.put_child(op_xml);
  }
  xmlout.put_child(coefs_out);
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
        if (!m_obs->queryBins(MCObsInfo(*corrt,args[kk]))) {
          xmlout.put_child(xml_corr);
          return;
        }
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
    bool multi_compare = false;
    if (xml_tag_count(xmltask,"MultiCompare")==1) multi_compare=true;
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

    if (!coefs_provided) determine_coefficients(m_obs,coefs,corrmats,minTime,multi_compare,xmlout);
    
    set<MCObsInfo> obskeys;

    vector<CorrelatorMatrixInfo>::iterator corrmat_it=corrmats.begin();
    CorrelatorMatrixInfo first_corrmat = *corrmat_it;
    set<OperatorInfo> result_operators;
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

      OperatorInfo snk_op(snkOpString,OperatorInfo::GenIrrep);
      OperatorInfo src_op(srcOpString,OperatorInfo::GenIrrep);
      result_operators.insert(snk_op);
      result_operators.insert(src_op);

      CorrelatorInfo corr_result(snk_op,src_op);
      CorrelatorAtTimeInfo corrt_result(corr_result,minTime,true,false);
      store_in_memory(m_obs,corrt_result,to_average,coefs,minTime,maxTime,obskeys,xmlout);
    }
    
    XMLHandler xmlf;
    m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
    xmlout.put_child(xmlf);

    CorrelatorMatrixInfo result_corrm(result_operators,true,false);

    XMLHandler result_corr_xml;
    result_corr_xml.set_root("Result");
    XMLHandler xml_corrm_out;
    result_corrm.output(xml_corrm_out);
    result_corr_xml.put_child(xml_corrm_out);
    xmlout.put_child(result_corr_xml);
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
