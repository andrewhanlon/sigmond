#include "chisq_3ptcorr.h"
#include "task_utils.h"
#include <string>
#include <map>
using namespace std;


// *************************************************************

/*
RealThreePointCorrelatorFit::RealThreePointCorrelatorFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
 XMLHandler xmlf(xmlin,"TemporalCorrelatorFit");
 OperatorInfo COp(xmlf);
 m_op=COp;
 m_subt_vev=false;
 if (xmlf.count_among_children("SubtractVEV")>0) m_subt_vev=true;

 T_period=OH.getLatticeTimeExtent();
 uint tmin,tmax;
 xmlreadchild(xmlf,"MinimumTimeSeparation",tmin);
 xmlreadchild(xmlf,"MaximumTimeSeparation",tmax);
 if (tmin<0) throw(std::invalid_argument("Invalid MinimumTimeSeparation"));
 if (tmax<=tmin) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));
 vector<int> texclude;
 string exstr;
 if (xmlreadifchild(xmlf,"ExcludeTimes",exstr)){
    extract_from_string(exstr,texclude);}
 m_tvalues=form_tvalues(tmin,tmax,texclude);  // ascending order

 m_noisecutoff=0.0;
 xmlreadifchild(xmlf,"LargeTimeNoiseCutoff",m_noisecutoff);

 XMLHandler xmlm(xmlf,"Model");
 string modeltype;
 xmlreadchild(xmlm,"Type",modeltype);

 try{
    create_tcorr_model(modeltype,T_period,m_model_ptr);
    m_nparams=m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);}
 catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in RealThreePointCorrelatorFit: ")
                 +string(errmsg.what())));}
 
 // read in priors

 m_npriors=0;
 if (xmlf.count("Priors")>0){
    XMLHandler xmlp(xmlf,"Priors");
    for (uint param_i = 0; param_i < m_nparams; ++param_i) {
      string param_name = m_model_ptr->getParameterName(param_i);
      if (xmlp.count(param_name) > 0) {
        XMLHandler xmlpp(xmlp, param_name);
        m_priors.insert(pair<uint,Prior>(param_i, Prior(xmlpp, OH)));
        m_npriors++;
      }
    }
 }

 setup();

}
*/

RealThreePointCorrelatorFit::RealThreePointCorrelatorFit(
    MCObsHandler& OH, CorrelatorInfo in_corr,
    const OperatorInfo& in_snk_op_2pt, const OperatorInfo& in_src_op_2pt,
    CorrelatorInfo in_rat_corr,
    bool subtractvev, ComplexArg in_arg, string model_name,
    const map<string,MCObsInfo>& model_params, map<uint,set<uint> > fit_times)
  : ChiSquare(OH), m_corr(in_corr), m_snk_op_2pt(in_snk_op_2pt), m_src_op_2pt(in_src_op_2pt),
    m_rat_corr(in_rat_corr), m_subt_vev(subtractvev), m_arg(in_arg), m_tvalues(fit_times)
{
  T_period = OH.getLatticeTimeExtent();
  try {
    create_3ptcorr_model(model_name, T_period, m_model_ptr);
    m_nparams = m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(model_params, m_fitparam_info);
  }
  catch(const std::exception& errmsg) {
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in RealThreePointCorrelatorFit: ")
        +string(errmsg.what())));
  }

  m_npriors = 0;
 
  setup();
}

RealThreePointCorrelatorFit::RealThreePointCorrelatorFit(
    MCObsHandler& OH, CorrelatorInfo in_corr, CorrelatorInfo in_rat_corr,
    bool subtractvev, ComplexArg in_arg, string model_name,
    const map<string,MCObsInfo>& model_params, map<uint,set<uint> > fit_times)
  : ChiSquare(OH), m_corr(in_corr), m_snk_op_2pt(in_corr.getSink()), m_src_op_2pt(in_corr.getSource()),
    m_rat_corr(in_rat_corr), m_subt_vev(subtractvev), m_arg(in_arg), m_tvalues(fit_times)
{
  T_period = OH.getLatticeTimeExtent();
  try {
    create_3ptcorr_model(model_name, T_period, m_model_ptr);
    m_nparams = m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(model_params, m_fitparam_info);
  }
  catch(const std::exception& errmsg) {
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in RealThreePointCorrelatorFit: ")
        +string(errmsg.what())));
  }

  m_npriors = 0;
 
  setup();
}

RealThreePointCorrelatorFit::RealThreePointCorrelatorFit(const RealThreePointCorrelatorFit& rtpcf)
    : ChiSquare(rtpcf), m_corr(rtpcf.m_corr), m_rat_corr(rtpcf.m_rat_corr), m_subt_vev(rtpcf.m_subt_vev), m_arg(rtpcf.m_arg), m_tvalues(rtpcf.m_tvalues)
{
  T_period = rtpcf.T_period;
  try{
    create_3ptcorr_model(rtpcf.m_model_ptr->getModelName(), T_period, m_model_ptr);
  }
  catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in RealThreePointCorrelatorFit: ")
        +string(errmsg.what())));
  }
}

void RealThreePointCorrelatorFit::setup()
{
  set<MCObsInfo> re_obskeys, im_obskeys;
  bool hermitian = false;

  if (m_model_ptr->getObsType() == Ratio) {
    doCorrelatorRatioBySamplings(m_obs, m_corr, m_snk_op_2pt, m_src_op_2pt, m_rat_corr, hermitian, m_subt_vev, re_obskeys, im_obskeys);
  }
  else if (m_model_ptr->getObsType() == SummationRatio) {
    doCorrelatorRatioSummationBySamplings(m_obs, m_corr, m_snk_op_2pt, m_src_op_2pt, m_rat_corr, hermitian, m_subt_vev, re_obskeys, im_obskeys);
    for (auto&& [tsep, tinss]: m_tvalues) {
      tinss.clear();
      tinss.insert(0);
    }
  }



  vector<MCObsInfo> obskeys;
  m_nobs = 0;
  CorrelatorAtTimeInfo rat_corrt(m_rat_corr, 0, 0, hermitian, m_subt_vev);
  for (auto [tsep, tinss]: m_tvalues) {
    rat_corrt.resetTimeSeparation(tsep);
    for (auto tins: tinss) {
      rat_corrt.resetTimeInsertion(tins);
      if (!re_obskeys.contains(MCObsInfo(rat_corrt, RealPart))) {
        throw(std::invalid_argument(string("Data not available for tsep = ")+make_string(tsep)
                      +string(", tins = ")+make_string(tins)+string(" for ")+m_rat_corr.str()));
      }
      obskeys.push_back(MCObsInfo(rat_corrt, m_arg));
      m_nobs++;
    }
  }

  int dof = m_nobs-m_nparams+m_npriors;
  if (dof < 1) throw(std::invalid_argument("Degrees of Freedom must be greater than zero"));

  allocate_obs_memory();

  for (uint k = 0; k < m_nobs; k++) {
    m_obs_info[k] = obskeys[k];
  }
}


RealThreePointCorrelatorFit::~RealThreePointCorrelatorFit()
{
 delete m_model_ptr;
}

double RealThreePointCorrelatorFit::evalModelPoint(
                               const vector<double>& fitparams,
                               const double tsep,
                               const double tins) const
{
  double model_point;
  m_model_ptr->evaluate(fitparams, tsep, tins, model_point);
  return model_point;
}

void RealThreePointCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
  uint k = 0;
  for (auto [tsep, tinss]: m_tvalues) {
    for (auto tins: tinss) {
      m_model_ptr->evaluate(fitparams, double(tsep), double(tins), modelpoints[k++]);
    }
  }
}


void RealThreePointCorrelatorFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
  uint k = 0;
  uint nparam=m_model_ptr->getNumberOfParams();
  vector<double> grad(nparam);
  for (auto [tsep, tinss]: m_tvalues) {
    for (auto tins: tinss) {
      m_model_ptr->evalGradient(fitparams, double(tsep), double(tins), grad);
      for (int p = 0; p < int(nparam); p++) {
        gradients(k,p) = grad[p];
      }
      k++;
    }
  }
}


void RealThreePointCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
  vector<double> corr(datapoints.size());
  for (uint k = 0; k < datapoints.size(); k++) {
    corr[k] = datapoints[k];
  }
  m_model_ptr->guessInitialParamValues(corr, m_tvalues, fitparams);  
}

string RealThreePointCorrelatorFit::getParameterName(uint param_index) const
{
  return m_model_ptr->getParameterName(param_index);
}

void RealThreePointCorrelatorFit::do_output(XMLHandler& xmlout) const
{
  /*
  xmlout.set_root("TemporalCorrelatorFit");
  XMLHandler xmlop;
  m_op.output(xmlop);
  xmlout.put_child(xmlop);
  if (m_subt_vev) xmlout.put_child("SubtractVEV");
  xmlout.put_child("TimeSeparations",make_string(m_tvalues));
  XMLHandler xmlmodel;
  m_model_ptr->output_tag(xmlmodel);
  if (m_npriors > 0) {
    XMLHandler xmlprior;
    xmlprior.set_root("Priors");
    for (map<uint,Prior>::const_iterator prior_it=m_priors.begin(); prior_it!=m_priors.end(); ++prior_it) {
      string param_name=getParameterName(prior_it->first);
      xmlprior.put_child(param_name);
    }
    xmlout.put_child(xmlprior);
  }
  
  xmlout.put_child(xmlmodel); 
  */
}


// *********************************************************************
