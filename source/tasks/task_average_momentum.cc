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
 re_mean=std::abs(est.getFullEstimate());
 re_err=est.getSymmetricError();
 est=m_obs->getJackknifeEstimate(obskeyIm);
 im_mean=std::abs(est.getFullEstimate());
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
    re_comp_mean=std::abs(est.getFullEstimate());
    re_comp_err=est.getSymmetricError();
    est=m_obs->getJackknifeEstimate(obskeyCompIm);
    im_comp_mean=std::abs(est.getFullEstimate());
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

void TaskHandler::doAverageMomentum(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
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

   // Put first CorrelatorMatrixInfo in first_corrMat, and
   // put the rest of the CorrelatorMatrixInfo's into corrMatInfos
   list<XMLHandler>::iterator ct=corrMatxml.begin();
   CorrelatorMatrixInfo first_corrMat(*ct);
   vector<CorrelatorMatrixInfo> corrMatInfos;
   for (++ct; ct!=corrMatxml.end();++ct) {
     corrMatInfos.push_back(CorrelatorMatrixInfo(*ct));
   }

   // Begin looping over the correlator elements of first_corrMat 
   set<MCObsInfo> obskeys;
   for (first_corrMat.begin(); !first_corrMat.end(); ++first_corrMat) {
     CorrelatorInfo corr=first_corrMat.getCurrentCorrelatorInfo();
     CorrelatorAtTimeInfo corrt(corr,minTime,true,false);
     // check the data exists
     MCObsInfo obskeyRe(corrt,RealPart);
     MCObsInfo obskeyIm(corrt,ImaginaryPart);
     if ((!m_obs->queryBins(obskeyRe))||(!m_obs->queryBins(obskeyIm))) continue;
     // build toAverage full of Correlators to average with corrt from the other Corrlator Matrices
     vector<CorrelatorAtTimeInfo> toAverage;
     for (vector<CorrelatorMatrixInfo>::iterator corrMat_it=corrMatInfos.begin();
          corrMat_it!=corrMatInfos.end(); ++corrMat_it) {
       uint num_compare = 0;
       for (corrMat_it->begin(); !corrMat_it->end(); ++(*corrMat_it)) {
         CorrelatorInfo corr_compare = corrMat_it->getCurrentCorrelatorInfo();
         if (corr.rotationallyEquivalent(corr_compare)) {
           num_compare++;
           if (num_compare==1) {
             toAverage.push_back(CorrelatorAtTimeInfo(corr_compare,minTime,true,false));
           }
         }
       }
       if (num_compare==0) throw(string("Could not find matching correlator"));
       if (num_compare>1) throw(string("Not a 1-to-1 match"));
     }
  
     vector<double> coefs = get_coefs(m_obs,corrt,toAverage,xmlout);
     store_in_memory(m_obs,corrt,toAverage,coefs,minTime,maxTime,obskeys,xmlout);
   }

   XMLHandler xmlff;
   m_obs->writeBinsToFile(obskeys,filename,xmlff,overwrite);
   xmlout.put_child(xmlff);
 }
 catch(const std::exception& errmsg){
   throw(std::invalid_argument((string("Invalid XML for task AverageMomentum: ")
        +string(errmsg.what())).c_str()));
 }
}
