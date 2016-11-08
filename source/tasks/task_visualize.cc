#include <iomanip>
#include "task_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

// *****************************************************************************
// *                                                                           *
// *    XML format for Visualizing Data                                        *
// *                                                                           *
// *    <Task>                                                                 *
// *     <Action>DoVisualization</Action>                                      *
// *     <Type>TemporalCorrelationMatrix</Type>                                *
// *     <CorrelatorMatrixInfo>                                                *
// *       <Operator>...</Operator>                                            *
// *       <Operator>...</Operator>                                            *
// *          ...                                                              *
// *       <HermitianMatrix/>    (optional)                                    *
// *       <SubtractVEV/>     (optional)                                       *
// *     </CorrelatorMatrixInfo>                                               *
// *     <MinTimeSep>3</MinTimeSep>                                            *
// *     <MaxTimeSep>25</MaxTimeSep>                                           *
// *     <Mode>Jackknife</Mode>  (optional)                                    *
// *                    (or Bootstrap or Current [default])                    *
// *     <NumberSize>6</NumberSize>                                            *
// *    </Task>                                                                *
// *                                                                           *
// *****************************************************************************

int getPrecision(double val, uint size_num)
{
  int precision = size_num - 1;
  if (val < 0.) --precision;
  if (abs(val) < 1.) {
    precision--;
  }
  else {
    precision -= (1 + (int)log10(abs(val)));
  }
  if (precision < 0) precision = 0;

  return precision;
}

void TaskHandler::doVisualization(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
  string visualType;
  xmlread(xmltask,"Type",visualType,"DoVisualization");

  if (visualType=="TemporalCorrelationMatrix") {
    xmlout.set_root("DoVisualization");
    xmlout.put_child("Type","TemporalCorrelationMatrix");
    try {
      XMLHandler xmlcm(xmltask,"CorrelatorMatrixInfo");
      CorrelatorMatrixInfo cormat(xmlcm);
      uint num_ops = cormat.getNumberOfOperators();
      const set<OperatorInfo>& ops = cormat.getOperators();
      uint size_num;
      uint tMin, tMax;
      xmlread(xmltask,"MinTimeSep",tMin,"DoVisualization");
      xmlread(xmltask,"MaxTimeSep",tMax,"DoVisualization");
      xmlread(xmltask,"NumberSize",size_num,"DoVisualization");
      string datamode="Current";
      xmlreadifchild(xmltask,"Mode",datamode);

      SamplingMode origmode=m_obs->getCurrentSamplingMode();
      if (datamode=="Bootstrap") m_obs->setToBootstrapMode();
      else if (datamode=="Jackknife") m_obs->setToJackknifeMode();
      xmlout.put_child("Mode",datamode);

      double re_mean,im_mean,re_err,im_err;
      for (uint t=tMin; t<=tMax; ++t) {
        XMLHandler xml_corr;
        xml_corr.set_root("CorrelatorMatrix");
        xml_corr.put_child("Time", to_string(t));
        stringstream corrmat_out;
        corrmat_out << setiosflags(ios::fixed);
        for (set<OperatorInfo>::const_iterator snk=ops.begin(); snk!=ops.end(); ++snk) {
          corrmat_out << endl << setfill('-')<<setw(num_ops*(2*size_num+7)+1)<<"-"<<endl << setfill(' ') << "|";
          for (uint i=0; i<num_ops;++i) corrmat_out << setw(2*size_num+7) << "|";
          corrmat_out << endl << "|";
          for (set<OperatorInfo>::const_iterator src=ops.begin(); src!=ops.end(); ++src) {
            CorrelatorAtTimeInfo corrt(*snk,*src,t,cormat.isHermitian(),cormat.isVEVSubtracted());
            MCObsInfo obskeyRe(corrt,RealPart);
            MCObsInfo obskeyIm(corrt,ImaginaryPart);
            if ((m_obs->queryBins(obskeyRe))&&(m_obs->queryBins(obskeyIm))) {
              string sign = "+";
              MCEstimate est;
              est=m_obs->getEstimate(obskeyRe);
              re_mean=est.getAverageEstimate();
              est=m_obs->getEstimate(obskeyIm);
              im_mean=est.getAverageEstimate();
              if (im_mean < 0.) {
                sign = "-";
                im_mean = abs(im_mean);
              }
              corrmat_out << " " << setprecision(getPrecision(re_mean,size_num)) << setw(size_num) << re_mean 
                          << " " << sign << " " << setprecision(getPrecision(im_mean,size_num)) << setw(size_num) << im_mean << "I" << " |";
            }
            else {
              corrmat_out << setw(size_num+6) << "No Data" << setw(size_num+1) << "|";
            }
          }
          corrmat_out << endl << "|";
          for (uint i=0; i<num_ops;++i) corrmat_out << setw(size_num+4) << "+/-" << setw(size_num+3) << "|";
          corrmat_out << endl << "|";
          for (set<OperatorInfo>::const_iterator src=ops.begin(); src!=ops.end(); ++src) {
            CorrelatorAtTimeInfo corrt(*snk,*src,t,cormat.isHermitian(),cormat.isVEVSubtracted());
            MCObsInfo obskeyRe(corrt,RealPart);
            MCObsInfo obskeyIm(corrt,ImaginaryPart);
            if ((m_obs->queryBins(obskeyRe))&&(m_obs->queryBins(obskeyIm))) {
              string sign = "+";
              MCEstimate est;
              est=m_obs->getEstimate(obskeyRe);
              re_err=est.getSymmetricError();
              est=m_obs->getEstimate(obskeyIm);
              im_err=est.getSymmetricError();
              corrmat_out << " " << setprecision(getPrecision(re_err,size_num)) << setw(size_num) << re_err 
                          << " + " << setprecision(getPrecision(im_err,size_num)) << setw(size_num) << im_err << "I" << " |";
            }
            else {
              corrmat_out << setw(size_num+6) << "No Data" << setw(size_num+1) << "|";
            }
          }
          corrmat_out << endl << "|";
          for (uint i=0; i<num_ops;++i) corrmat_out << setw(2*size_num+7) << "|";
        }
        corrmat_out << endl << setfill('-')<<setw(num_ops*(2*size_num+7)+1)<<"-"<<endl << setfill(' ') << "|" << endl;
        xml_corr.put_child_text_node_whitespace(corrmat_out.str());
        xmlout.put_child(xml_corr);
      }
      m_obs->setSamplingMode(origmode);
    }
    catch(const std::exception& errmsg){
      xmlout.clear();
      throw(std::invalid_argument((string("DoVisualization with type TemporalCorrelationMatrix encountered an error: ")
            +string(errmsg.what())).c_str()));
    }
  }
}
