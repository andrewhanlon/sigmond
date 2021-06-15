#include "chisq_disp.h"
#include "task_utils.h"
#include <string>
#include <map>
using namespace std;


// *************************************************************


DispersionFit::DispersionFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
  XMLHandler xmlf(xmlin, "DispersionFit");

  Xextent = OH.getLatticeXExtent();
  Yextent = OH.getLatticeYExtent();
  Zextent = OH.getLatticeZExtent();

  XMLHandler xmlm(xmlf, "Model");
  string modeltype;
  xmlreadchild(xmlm, "Type", modeltype);

  try{
    create_disp_model(modeltype, Xextent, Yextent, Zextent, m_model_ptr);
    m_nparams = m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);
  }
  catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in DispersionFit: ")
                 +string(errmsg.what())));
  }

  map<MCObsInfo,double> energies;
  list<XMLHandler> xmlobs = xmlf.find("Energy");
  for (list<XMLHandler>::iterator it = xmlobs.begin(); it != xmlobs.end(); it++) {
    string obsname;
    xmlreadchild(*it, "Name", obsname);
    uint index=taskcount;
    xmlreadifchild(*it, "IDIndex", index);
    MCObsInfo obskey(obsname, index);
    XMLHandler mom_xml(*it, "Momentum");
    vector<int> mom_ints;
    string mom_str = mom_xml.get_text_content();
    extract_from_string(mom_str, mom_ints);
    Momentum mom(mom_ints[0], mom_ints[1], mom_ints[2]);
    double psq = mom.getPsq(Xextent, Yextent, Zextent);
    energies.insert(pair<MCObsInfo,double>(obskey, psq));
  }
 
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

  setup(energies);
}

DispersionFit::DispersionFit(
    MCObsHandler& OH, const string& model_name, const map<MCObsInfo,double>& energies,
    const map<string,MCObsInfo>& model_params) : ChiSquare(OH)
{
  Xextent = OH.getLatticeXExtent();
  Yextent = OH.getLatticeYExtent();
  Zextent = OH.getLatticeZExtent();

  try{
    create_disp_model(model_name, Xextent, Yextent, Zextent, m_model_ptr);
    m_nparams = m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(model_params, m_fitparam_info);
  }
  catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in DispersionFit: ") +string(errmsg.what())));
  }

  m_npriors = 0;

  setup(energies);
}

void DispersionFit::setup(const map<MCObsInfo,double>& energies)
{
  m_nobs = energies.size();
  allocate_obs_memory();
  for (uint k = 0; auto& [energy_info, psq]: energies) {
    if (!m_obs->queryFullAndSamplings(energy_info)) {
      throw(std::invalid_argument(string("Observation: ")+energy_info.str() +string(" not available")));
    }
    m_momsq.push_back(psq);
    MCObsInfo esq_info(energy_info.getObsName() + "_sq", energy_info.getObsIndex());
    m_obs_info[k++] = esq_info;
    doSquareBySamplings(*m_obs, energy_info, esq_info);
  }

  int dof = m_nobs-m_nparams+m_npriors;
  if (dof < 1) throw(std::invalid_argument("Degrees of Freedom must be greater than zero"));
}

DispersionFit::~DispersionFit()
{
 delete m_model_ptr;
}

double DispersionFit::evalDispersion(const double msq, const double psq) const
{
  double disp_value;
  m_model_ptr->evalDispersion(msq, psq, disp_value);
  return disp_value;
}

double DispersionFit::evalModelPoint(const vector<double>& fitparams,
                                     const double psq) const
{
  double model_point;
  m_model_ptr->evaluate(fitparams, psq, model_point);
  return model_point;
}

void DispersionFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
  for (uint k = 0; k < m_nobs; ++k) {
    double psq = m_momsq[k];
    m_model_ptr->evaluate(fitparams, psq, modelpoints[k]);
  }
}


void DispersionFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
  uint nparam = m_model_ptr->getNumberOfParams();
  vector<double> grad(nparam);
  for (uint k = 0; k < m_nobs; ++k) {
    double psq = m_momsq[k];
    m_model_ptr->evalGradient(fitparams, psq, grad);
    for (uint p=0 ; p < nparam; ++p) {
      gradients(k,p) = grad[p];
    }
  }
}


void DispersionFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
  vector<double> energies(datapoints.size());
  for (uint k = 0; k< m_nobs; ++k) {
    energies[k] = datapoints[k];
  }
  m_model_ptr->guessInitialParamValues(energies, m_momsq, fitparams);  
}


string DispersionFit::getParameterName(uint param_index) const
{
 return m_model_ptr->getParameterName(param_index);
}

void DispersionFit::do_output(XMLHandler& xmlout) const
{
  xmlout.set_root("DispersionFit");
  for (uint k = 0; k < m_nobs; k++) {
    XMLHandler xmlp("EnergyObservable");
    XMLHandler xmln;
    m_obs_info[k].output(xmln); 
    xmlp.put_child(xmln); 
    xmlp.put_child("MomSquared", to_string(m_momsq[k]));
    xmlout.put_child(xmlp);
  }

  XMLHandler xmlmodel;
  m_model_ptr->output_tag(xmlmodel);
  if (m_npriors>0){
    XMLHandler xmlprior;
    xmlprior.set_root("Priors");
    for (auto& [param_i, prior]: m_priors) {
      string param_name = getParameterName(param_i);
      xmlprior.put_child(param_name);
    }
    xmlout.put_child(xmlprior);
  }
  xmlout.put_child(xmlmodel); 
}

MCObsInfo DispersionFit::getMomSqObsParamInfo(const double psq) const
{
  auto momsq_it = find(m_momsq.begin(), m_momsq.end(), psq);
  if (momsq_it != m_momsq.end()) {
    int index = distance(m_momsq.begin(), momsq_it);
    return m_obs_info.at(index);
  }
  else {
    throw(std::invalid_argument(string("No energy with that psq")));
  }
}

// *********************************************************************
