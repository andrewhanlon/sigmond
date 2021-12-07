#include "create_plots.h"
#include "xml_handler.h"
#include "mc_estimate.h"
#include "task_utils.h"
#include <stdexcept>

using namespace std;

// *************************************************************

void createEffEnergyWithFitPlot(EffEnergyWithFitPlotInfo plot_info, RealTemporalCorrelatorFit& rtc,
                                FitResult& fit_result, MCObsHandler* m_obs, XMLHandler& xmlout)
{
  if (plot_info.plotfile.empty()) {
    xmlout.put_child("Warning","No plot file but asked for plot!");
    return;
  }
  char goodtype = 'N';
  double goodness = fit_result.quality;
  if (plot_info.goodness == "qual") {
    goodtype = 'Q';
  }
  else if (plot_info.goodness == "chisq") {
    goodtype = 'X';
    goodness = fit_result.chisq_dof;
  }
  CorrelatorInfo corr(rtc.getOperatorInfo(), rtc.getOperatorInfo());
  string corrname = plot_info.corrname;
  if (corrname=="standard") corrname = getCorrelatorStandardName(corr);
  bool hermitian = true;
  bool subvev = rtc.subtractVEV();
  uint fit_tmin = rtc.getTmin();
  uint fit_tmax = rtc.getTmax();
  uint efftype = rtc.getEffMassType();
  double subt_const=0.0;
   // if (efftype>1){    // subtract fit constant
   //    efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
   //    subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
  SamplingMode mode=m_obs->getCurrentSamplingMode();

  map<double,MCEstimate> results;
  getEffectiveEnergy(m_obs, corr, hermitian, subvev, RealPart, mode, plot_info.timestep, efftype,
                     results, subt_const);
  if (results.empty()) {
    xmlout.put_child("PlotError","No effective energy estimates could be obtained");
    return;
  }
       // do some XML output
  xmlout.put_child("PlotFile", plot_info.plotfile);
  XMLHandler xmlef;
  xmlef.set_root("EffectiveEnergy");
  xmlef.put_child("TimeStep", make_string(plot_info.timestep));
  if (efftype==0) xmlef.put_child("EffEnergyType", "TimeForward");
  else if (efftype==1) xmlef.put_child("EffEnergyType", "TimeSymmetric");
  else if (efftype==2) xmlef.put_child("EffEnergyType", "TimeForwardPlusConst");
  else if (efftype==3) xmlef.put_child("EffEnergyType", "TimeSymmetricPlusConst");
  xmlef.seek_root();
  xmlef.seek_first_child();
  for (map<double,MCEstimate>::const_iterator rt = results.begin(); rt != results.end(); ++rt) {
    XMLHandler xmlr("Estimate");
    xmlr.put_child("TimeSeparation", make_string(rt->first));
    xmlr.put_child("MeanValue", make_string((rt->second).getFullEstimate()));
    xmlr.put_child("SymmError", make_string((rt->second).getSymmetricError()));
    xmlef.put_sibling(xmlr);
  }
  xmlout.put_child(xmlef);
  if (plot_info.maxerror > 0.) {
    map<double,MCEstimate> raw(results);
    results.clear();
    for (map<double,MCEstimate>::const_iterator it = raw.begin(); it != raw.end(); ++it) {
      if ((it->second).getSymmetricError() < abs(plot_info.maxerror)) results.insert(*it);
    }
  }

  vector<XYDYPoint> meffvals(results.size());
  uint k=0;
  for (map<double,MCEstimate>::const_iterator rt = results.begin(); rt != results.end(); ++rt, k++) {
    meffvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(), (rt->second).getSymmetricError());
  }

  vector<XYPoint> meff_approach = rtc.getEffEnergyApproach(fit_result.bestfit_params, plot_info.timestep);
  double energy_mean = fit_result.bestfit_params[0].getFullEstimate();
  double energy_err = fit_result.bestfit_params[0].getSymmetricError();

  if (plot_info.ref_energy.isVacuum()) {
    createEffEnergyWithFitPlot(meffvals, RealPart, energy_mean, energy_err, fit_tmin, fit_tmax,
                               meff_approach, goodtype, goodness, corrname, plot_info.plotfile,
                               plot_info.symboltype, plot_info.symbolcolor);
  }
  else {
    MCObsInfo enratio(string("TempEnergyRatioGwiqb"), plot_info.ref_energy.getObsIndex());
    doRatioBySamplings(*m_obs, rtc.getFitParamInfos()[0], plot_info.ref_energy, enratio);
    MCEstimate ratioest=m_obs->getEstimate(enratio);
    XMLHandler xmlrat("EnergyRatioFitResult");
    XMLHandler xmlrr;
    ratioest.output(xmlrr); xmlrat.put_child(xmlrr);
    xmlout.put_child(xmlrat);
    double energy_ratio = ratioest.getFullEstimate();
    double energy_ratio_err = ratioest.getSymmetricError();
    createEffEnergyWithFitAndEnergyRatioPlot(
        meffvals, RealPart, energy_mean, energy_err, fit_tmin, fit_tmax, meff_approach,
        energy_ratio, energy_ratio_err, goodtype, goodness, corrname, plot_info.plotfile,
        plot_info.symboltype, plot_info.symbolcolor);
    m_obs->eraseData(enratio);
  }
}

void createThreePointCorrelatorPlot(ThreePointCorrelatorPlotInfo plot_info,
                            const CorrelatorInfo& corr_3pt, bool subtract_vev,
                            MCObsHandler* m_obs, XMLHandler& xmlout)
{
  createThreePointCorrelatorPlot(plot_info, corr_3pt, corr_3pt.getSink(), corr_3pt.getSource(),
                                 subtract_vev, m_obs, xmlout);
}

void createThreePointCorrelatorPlot(ThreePointCorrelatorPlotInfo plot_info,
                            const CorrelatorInfo& corr_3pt,
                            const OperatorInfo& two_pt_snk_op, const OperatorInfo& two_pt_src_op,
                            bool subtract_vev, MCObsHandler* m_obs, XMLHandler& xmlout)
{
  if (plot_info.plotfile.empty()) {
    xmlout.put_child("Warning", "No plot file but asked for plot!");
    return;
  }

  bool hermitian = false;
  SamplingMode mode = m_obs->getCurrentSamplingMode();
  map<pair<double,double>,MCEstimate> real_results;
  map<pair<double,double>,MCEstimate> imag_results;
  getCorrelatorRatioEstimates(m_obs, corr_3pt, two_pt_snk_op, two_pt_src_op, hermitian,
                              subtract_vev, mode, real_results, imag_results);

  map<uint,set<uint> > timesavailable;
  getCorrelatorRatioAvailableTimes(m_obs, timesavailable, corr_3pt, two_pt_snk_op, two_pt_src_op,
                                   hermitian, true, true);
  vector<vector<XYDYPoint> > corrvals;
  for (auto [tsep, tinss]: timesavailable) {
    vector<XYDYPoint> tsep_corrvals;
    for (auto tins: tinss) {
      double x = tins - tsep/2.;
      MCEstimate est;
      if (plot_info.complex_arg == RealPart) {
        est = real_results.at(make_pair(double(tsep),double(tins)));
      }
      else {
        est = imag_results.at(make_pair(double(tsep),double(tins)));
      }
      tsep_corrvals.push_back(XYDYPoint(x, est.getFullEstimate(), est.getSymmetricError()));
    }
    corrvals.push_back(tsep_corrvals);
  }
  if (corrvals.empty()) {
    xmlout.put_child("PlotError","No estimates could be obtained");
    return;
  }

  createThreePointCorrelatorPlot(corrvals, plot_info.complex_arg, plot_info.plotlabel, plot_info.plotfile,
                         plot_info.labels, plot_info.symboltypes, plot_info.symbolcolors);
}

void createThreePointCorrelatorWithFitPlot(ThreePointCorrelatorWithFitPlotInfo plot_info,
                                RealThreePointCorrelatorFit& rtpc,
                                FitResult& fit_result,
                                MCObsHandler* m_obs,
                                XMLHandler& xmlout)
{
  if (plot_info.plotfile.empty()) {
    xmlout.put_child("Warning", "No plot file but asked for plot!");
    return;
  }

  bool hermitian = false;

  char goodtype = 'N';
  double goodness = fit_result.quality;
  if (plot_info.goodness == "qual") {
    goodtype = 'Q';
  }
  else if (plot_info.goodness == "chisq") {
    goodtype = 'X';
    goodness = fit_result.chisq_dof;
  }

  double b0_mean = fit_result.bestfit_params[0].getFullEstimate();
  double b0_err = fit_result.bestfit_params[0].getSymmetricError();

  SamplingMode mode = m_obs->getCurrentSamplingMode();

  if (rtpc.getObsType() == Ratio) {
    uint npoints = 400;
    map<uint,set<uint> > fit_times = rtpc.getTvalues();
    vector<vector<XYPoint> > lowerfits;
    vector<vector<XYPoint> > upperfits;

    vector<MCObsInfo> fitparam_infos = rtpc.getFitParamInfos();
    vector<double> fitparams(fitparam_infos.size());
    MCObsInfo func_info("func_val_temp", 10000);

    map<uint,set<uint> > timesavailable;
    getCorrelatorRatioAvailableTimes(m_obs, timesavailable, rtpc.getCorrelatorInfo(),
                                     rtpc.getTwoPointSinkOperatorInfo(), rtpc.getTwoPointSourceOperatorInfo(),
                                     hermitian, true, true);
    for (auto [tsep, _]: timesavailable) {
      if (find(plot_info.time_seps.begin(), plot_info.time_seps.end(), tsep) == plot_info.time_seps.end()) continue;
      if (fit_times.contains(tsep)) {
        set<uint> tinss = fit_times.at(tsep);
        vector<XYPoint> lowerfit(npoints+1);
        vector<XYPoint> upperfit(npoints+1);
        double lower_x = *(tinss.begin()) - tsep/2.;
        double upper_x = *(tinss.rbegin()) - tsep/2.;

        for (uint i = 0; i <= npoints; ++i) {
          double x = lower_x + i*(upper_x - lower_x)/npoints;
          double tins = x + tsep/2.;

          for (m_obs->setSamplingBegin(); !m_obs->isSamplingEnd(); m_obs->setSamplingNext()) {
            for (uint k = 0; auto param_info: fitparam_infos) {
              fitparams[k++] = m_obs->getCurrentSamplingValue(param_info);
            }
            double func_val = rtpc.evalModelPoint(fitparams, double(tsep), tins);
            m_obs->putCurrentSamplingValue(func_info, func_val);
          }

          MCEstimate func_est = m_obs->getEstimate(func_info);
          lowerfit[i].xval = x;
          lowerfit[i].yval = func_est.getFullEstimate() - func_est.getSymmetricError();
          upperfit[i].xval = x;
          upperfit[i].yval = func_est.getFullEstimate() + func_est.getSymmetricError();
        }

        reverse(lowerfit.begin(), lowerfit.end());
        lowerfits.push_back(lowerfit);
        upperfits.push_back(upperfit);
      }
      else {
        lowerfits.push_back(vector<XYPoint>());
        upperfits.push_back(vector<XYPoint>());
      }
    }

    map<pair<double,double>,MCEstimate> real_results;
    map<pair<double,double>,MCEstimate> imag_results;
    getCorrelatorRatioEstimates(m_obs, rtpc.getCorrelatorInfo(), rtpc.getTwoPointSinkOperatorInfo(),
                                rtpc.getTwoPointSourceOperatorInfo(), hermitian, rtpc.subtractVEV(),
                                mode, real_results, imag_results);

    vector<vector<XYDYPoint> > corrvals;
    for (auto [tsep, tinss]: timesavailable) {
      if (find(plot_info.time_seps.begin(), plot_info.time_seps.end(), tsep) == plot_info.time_seps.end()) continue;
      vector<XYDYPoint> tsep_corrvals;
      for (auto tins: tinss) {
        double x = tins - tsep/2.;
        MCEstimate est;
        if (rtpc.getComplexArg() == RealPart) {
          est = real_results.at(make_pair(double(tsep),double(tins)));
        }
        else {
          est = imag_results.at(make_pair(double(tsep),double(tins)));
        }
        tsep_corrvals.push_back(XYDYPoint(x, est.getFullEstimate(), est.getSymmetricError()));
      }
      corrvals.push_back(tsep_corrvals);
    }
    if (corrvals.empty()) {
      xmlout.put_child("PlotError","No estimates could be obtained");
      return;
    }

    createThreePointCorrelatorWithFitPlot(corrvals, rtpc.getComplexArg(), b0_mean, b0_err,
                                          lowerfits, upperfits, goodtype, goodness,
                                          plot_info.plotlabel, plot_info.plotfile, plot_info.labels,
                                          plot_info.symboltypes, plot_info.symbolcolors);
  }
  else {
    map<double,MCEstimate> real_results;
    map<double,MCEstimate> imag_results;
    getCorrelatorRatioSummationEstimates(m_obs, rtpc.getCorrelatorInfo(), rtpc.getTwoPointSinkOperatorInfo(),
                                rtpc.getTwoPointSourceOperatorInfo(), hermitian, rtpc.subtractVEV(),
                                mode, real_results, imag_results);

    map<uint,set<uint> > timesavailable;
    getCorrelatorRatioAvailableTimes(m_obs, timesavailable, rtpc.getCorrelatorInfo(),
                                     rtpc.getTwoPointSinkOperatorInfo(), rtpc.getTwoPointSourceOperatorInfo(),
                                     hermitian, true, true);
    vector<XYDYPoint> summvals;
    double min_tsep = 0;
    double max_tsep = 0;
    for (auto [tsep, _]: timesavailable) {
      if ((min_tsep == 0) || (min_tsep > tsep)) {
        min_tsep = double(tsep);
      }
      if ((max_tsep == 0) || (max_tsep < tsep)) {
        max_tsep = double(tsep);
      }

      MCEstimate est;
      if (rtpc.getComplexArg() == RealPart) {
        est = real_results.at(double(tsep));
      }
      else {
        est = imag_results.at(double(tsep));
      }
      summvals.push_back(XYDYPoint(double(tsep), est.getFullEstimate(), est.getSymmetricError()));
    }
    if (summvals.empty()) {
      xmlout.put_child("PlotError","No estimates could be obtained");
      return;
    }

    min_tsep -= .5;
    max_tsep += .5;

    uint npoints = 400;
    vector<XYPoint> lowerfit(npoints+1);
    vector<XYPoint> upperfit(npoints+1);

    vector<MCObsInfo> fitparam_infos = rtpc.getFitParamInfos();
    vector<double> fitparams(fitparam_infos.size());
    MCObsInfo func_info("func_val_temp", 10000);
    for (uint i = 0; i <= npoints; ++i) {
      double tsep = min_tsep + i*(max_tsep - min_tsep)/npoints;
      for (m_obs->setSamplingBegin(); !m_obs->isSamplingEnd(); m_obs->setSamplingNext()) {
        for (uint k = 0; auto param_info: fitparam_infos) {
          fitparams[k++] = m_obs->getCurrentSamplingValue(param_info);
        }
        double func_val = rtpc.evalModelPoint(fitparams, tsep, 0.);
        m_obs->putCurrentSamplingValue(func_info, func_val);
      }
      MCEstimate func_est = m_obs->getEstimate(func_info);
      lowerfit[i].xval = tsep;
      lowerfit[i].yval = func_est.getFullEstimate() - func_est.getSymmetricError();
      upperfit[i].xval = tsep;
      upperfit[i].yval = func_est.getFullEstimate() + func_est.getSymmetricError();
    }

    createThreePointCorrelatorWithFitPlot(summvals, rtpc.getComplexArg(), b0_mean, b0_err,
                               lowerfit, upperfit, goodtype, goodness, plot_info.plotlabel, plot_info.plotfile,
                               plot_info.symboltypes.at(0), plot_info.symbolcolors.at(0));
  }
}

void createDataFitRatioPlot(DataFitRatioPlotInfo plot_info, vector<RealTemporalCorrelatorFit>& rtcs,
                            MCObsHandler* m_obs, XMLHandler& xmlout)
{
  if (plot_info.plotfile.empty()) {
    xmlout.put_child("Warning", "No plot file but asked for plot!");
    return;
  }

  vector<vector<XYDYPoint> > all_ratiovals;
  all_ratiovals.reserve(rtcs.size());
  bool hermitian = true;

  for (vector<RealTemporalCorrelatorFit>::iterator rtc_it = rtcs.begin(); rtc_it != rtcs.end(); ++rtc_it) {
    CorrelatorInfo corr(rtc_it->getOperatorInfo(), rtc_it->getOperatorInfo());
    bool subvev = rtc_it->subtractVEV();

    set<uint> tseps;
    getCorrelatorAvailableTimes(m_obs, tseps, corr, hermitian, RealPart);
    map<uint,MCObsInfo> ratios;
    for (set<uint>::iterator t_it = tseps.begin(); t_it != tseps.end(); ++t_it) {
      ratios.insert(pair<uint,MCObsInfo>(*t_it, MCObsInfo("ratio_temp", *t_it)));
    }

    CorrelatorAtTimeInfo corr_t(corr, 0, hermitian, subvev);

    vector<MCObsInfo> fitparam_infos = rtc_it->getFitParamInfos();
    vector<double> fitparams;
    fitparams.reserve(fitparam_infos.size());
    for (m_obs->setSamplingBegin(); !m_obs->isSamplingEnd(); m_obs->setSamplingNext()) {
      for (vector<MCObsInfo>::iterator param_it = fitparam_infos.begin(); param_it != fitparam_infos.end(); ++param_it) {
        fitparams.push_back(m_obs->getCurrentSamplingValue(*param_it));
      }

      for (map<uint,MCObsInfo>::iterator ratio_it = ratios.begin(); ratio_it != ratios.end(); ++ratio_it) {
        double model_point = rtc_it->evalModelPoint(fitparams, ratio_it->first);
        corr_t.resetTimeSeparation(ratio_it->first);
        MCObsInfo corr_t_obs(corr_t, RealPart);
        double ratio = m_obs->getCurrentSamplingValue(corr_t_obs) / model_point;
        m_obs->putCurrentSamplingValue(ratio_it->second, ratio);
      }
      fitparams.clear();
    }

    map<double,MCEstimate> results;
    for (map<uint,MCObsInfo>::iterator ratio_it = ratios.begin(); ratio_it != ratios.end(); ++ratio_it) {
      results.insert(pair<double,MCEstimate>(ratio_it->first, m_obs->getEstimate(ratio_it->second)));
    }

    if (results.empty()) {
      xmlout.put_child("PlotError","No effective energy estimates could be obtained");
      return;
    }
         // do some XML output
    xmlout.put_child("PlotFile", plot_info.plotfile);
    /*
    XMLHandler xmlef;
    xmlef.set_root("EffectiveEnergy");
    xmlef.put_child("TimeStep", make_string(plot_info.timestep));
    if (efftype==0) xmlef.put_child("EffEnergyType", "TimeForward");
    else if (efftype==1) xmlef.put_child("EffEnergyType", "TimeSymmetric");
    else if (efftype==2) xmlef.put_child("EffEnergyType", "TimeForwardPlusConst");
    else if (efftype==3) xmlef.put_child("EffEnergyType", "TimeSymmetricPlusConst");
    xmlef.seek_root();
    xmlef.seek_first_child();
    for (map<double,MCEstimate>::const_iterator rt = results.begin(); rt != results.end(); ++rt) {
      XMLHandler xmlr("Estimate");
      xmlr.put_child("TimeSeparation", make_string(rt->first));
      xmlr.put_child("MeanValue", make_string((rt->second).getFullEstimate()));
      xmlr.put_child("SymmError", make_string((rt->second).getSymmetricError()));
      xmlef.put_sibling(xmlr);
    }
    xmlout.put_child(xmlef);
    */
    if (plot_info.maxerror > 0.) {
      map<double,MCEstimate> raw(results);
      results.clear();
      for (map<double,MCEstimate>::const_iterator it = raw.begin(); it != raw.end(); ++it) {
        if ((it->second).getSymmetricError() < abs(plot_info.maxerror)) results.insert(*it);
      }
    }

    vector<XYDYPoint> ratiovals(results.size());
    uint k=0;
    for (map<double,MCEstimate>::const_iterator rt = results.begin(); rt != results.end(); ++rt, k++) {
      ratiovals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(), (rt->second).getSymmetricError());
    }
    all_ratiovals.push_back(ratiovals);
  }

  createDataFitRatioPlot(all_ratiovals, plot_info.plotlabel, plot_info.plotfile,
                         plot_info.labels, plot_info.symboltypes, plot_info.symbolcolors);
}

void createTminPlot(TminFitPlotInfo plot_info,
                    vector<vector<RealTemporalCorrelatorFit> >& rtcs,
                    vector<vector<FitResult> >& fit_results,
                    MCObsHandler* m_obs,
                    XMLHandler& xmlout)
{
  if (plot_info.plotfile.empty()) {
    xmlout.put_child("Warning", "No plot file but asked for plot!");
    return;
  }

  if (rtcs.size() != fit_results.size()) {
    throw(invalid_argument("rtcs and fit_results must be same size"));
  }

  vector<vector<XYDYDYPoint> > all_goodfits;
  all_goodfits.reserve(rtcs.size());
  vector<vector<XYDYDYPoint> > all_badfits;
  all_badfits.reserve(rtcs.size());

  for (uint fit_type_i = 0; fit_type_i < rtcs.size(); ++fit_type_i) {
    if (rtcs[fit_type_i].size() != fit_results[fit_type_i].size()) {
      throw(invalid_argument("rtcs and fit_results must be same size"));
    }

    vector<XYDYDYPoint> goodfits;
    vector<XYDYDYPoint> badfits;
    for (uint fit_i = 0; fit_i < rtcs[fit_type_i].size(); ++fit_i) {
      uint tmin = rtcs[fit_type_i][fit_i].getTmin();
      MCObsInfo energy_info("full_energy_temp", 100000);
      rtcs[fit_type_i][fit_i].setEnergyInfo(plot_info.energy_level, energy_info);
      MCEstimate energy_est = m_obs->getEstimate(energy_info);

      double y = energy_est.getFullEstimate();
      double dyup, dydn;
      if (energy_est.isJackknifeMode()) {
        dyup = dydn = energy_est.getSymmetricError();
      }
      else {
        dyup = energy_est.getUpperConfLimit() - y;
        dydn = y - energy_est.getLowerConfLimit();
      }

      if ((abs(dyup) > plot_info.maxerror) || (abs(dydn) > plot_info.maxerror)) {
        continue;
      }

      if (fit_results[fit_type_i][fit_i].quality >= plot_info.quality_threshold) {
        goodfits.push_back(XYDYDYPoint(tmin, y, dyup, dydn));
      }
      else {
        badfits.push_back(XYDYDYPoint(tmin, y, dyup, dydn));
      }
    }

    all_goodfits.push_back(goodfits);
    all_badfits.push_back(badfits);
  }



  XYDYDYPoint chosen_fit(0, 0, 0, 0);
  if (!plot_info.chosen_fit.isVacuum()) {
    MCEstimate chosen_fit_estimate = m_obs->getEstimate(plot_info.chosen_fit);
    double y = chosen_fit_estimate.getFullEstimate();
    double dyup, dydn;
    if (chosen_fit_estimate.isJackknifeMode()) {
      dyup = dydn = chosen_fit_estimate.getSymmetricError();
    }
    else {
      dyup = chosen_fit_estimate.getUpperConfLimit() - y;
      dydn = y - chosen_fit_estimate.getLowerConfLimit();
    }
    chosen_fit = XYDYDYPoint(1 , y, dyup, dydn);
  }

  createTminPlot(all_goodfits, all_badfits, plot_info.plotlabel, plot_info.plotfile,
                 plot_info.energy_level, plot_info.labels, plot_info.symboltypes,
                 plot_info.symbolcolors, chosen_fit);
}


void createDispersionFitPlot(DispersionFitPlotInfo plot_info, DispersionFit& disp_fit,
                             FitResult& fit_result, MCObsHandler* m_obs, XMLHandler& xmlout,
                             bool show_fit)
{
  if (plot_info.plotfile.empty()) {
    xmlout.put_child("Warning","No plot file but asked for plot!");
    return;
  }
  char goodtype = 'N';
  double goodness = fit_result.quality;
  if (plot_info.goodness == "qual") {
    goodtype = 'Q';
  }
  else if (plot_info.goodness == "chisq") {
    goodtype = 'X';
    goodness = fit_result.chisq_dof;
  }

  const vector<MCObsInfo> obs_infos = disp_fit.getObsInfos();
  const vector<double> momsqs = disp_fit.getMomSqs();
  vector<XYDYPoint> energy_sqs(disp_fit.getNumberOfObervables());
  for (uint k=0; k < disp_fit.getNumberOfObervables(); ++k) {
    MCEstimate est = m_obs->getEstimate(obs_infos[k]);
    energy_sqs[k].xval = momsqs[k];
    energy_sqs[k].yval = est.getFullEstimate();
    energy_sqs[k].yerr = est.getSymmetricError();
  }

  uint npoints = 400;
  vector<XYPoint> lowerfit(npoints+1);
  vector<XYPoint> upperfit(npoints+1);
  vector<XYPoint> disp(npoints+1);
  double lower_psq = disp_fit.getMinMomSq();
  double upper_psq = disp_fit.getMaxMomSq();

  MCObsInfo energy_sq_info = disp_fit.getMomSqObsParamInfo(lower_psq);
  MCObsInfo mass_sq_info = MCObsInfo("temp_mass", 10000);
  doBoostBySamplings(*m_obs, energy_sq_info, -lower_psq, mass_sq_info);
  m_obs->setSamplingBegin();
  double msq = m_obs->getCurrentSamplingValue(mass_sq_info);

  vector<MCObsInfo> fitparam_infos = disp_fit.getFitParamInfos();
  vector<double> fitparams(fitparam_infos.size());
  MCObsInfo esq_info("esq_val_temp", 10000);
  for (uint i = 0; i <= npoints; ++i) {
    double psq = lower_psq + i*(upper_psq - lower_psq)/npoints;
    disp[i].xval = psq;
    disp[i].yval = disp_fit.evalDispersion(msq, psq);
    if (show_fit) {
      for (m_obs->setSamplingBegin(); !m_obs->isSamplingEnd(); m_obs->setSamplingNext()) {
        for (uint k = 0; auto param_info: fitparam_infos) {
          fitparams[k++] = m_obs->getCurrentSamplingValue(param_info);
        }
        double esq_val = disp_fit.evalModelPoint(fitparams, psq);
        m_obs->putCurrentSamplingValue(esq_info, esq_val);
      }
      MCEstimate esq_est = m_obs->getEstimate(esq_info);
      lowerfit[i].xval = psq;
      lowerfit[i].yval = esq_est.getFullEstimate() - esq_est.getSymmetricError();
      upperfit[i].xval = psq;
      upperfit[i].yval = esq_est.getFullEstimate() + esq_est.getSymmetricError();
    }
  }

  vector<double> param_means(fitparam_infos.size());
  vector<double> param_errs(fitparam_infos.size());
  if (show_fit) {
    for (uint k = 0; auto param_info: fitparam_infos) {
      MCEstimate param_est = m_obs->getEstimate(param_info);
      param_means[k] = param_est.getFullEstimate();
      param_errs[k] = param_est.getSymmetricError();
      k++;
    }
  }

  createDispersionFitPlot(energy_sqs, lowerfit, upperfit, disp, param_means, param_errs,
                          goodtype, goodness, plot_info.particle_name, plot_info.plotfile,
                          plot_info.symboltype, plot_info.symbolcolor, show_fit);
}


// *************************************************************

void createMCValuesPlot(const Vector<double>& mcvalues, const string& observable_name,
                        double in_mean_value, double in_std_dev,
                        const string& filename, const string& symbol, 
                        const string& symbolcolor, double rescale, bool drawtoscreen)
{
 GracePlot P("Markov Chain Index","Value");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double std_dev=rescale*in_std_dev;

 P.addXYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<mcvalues.size();ind++)
    P.addXYDataPoint(double(ind),rescale*mcvalues[ind]);

 P.addXYDataSet("none","none","solid",symbolcolor);
 P.addXYDataPoint(0,mean_value); P.addXYDataPoint(mcvalues.size(),mean_value);
 P.addXYDataSet("none","none","dash",symbolcolor);
 P.addXYDataPoint(0,mean_value+std_dev); P.addXYDataPoint(mcvalues.size(),mean_value+std_dev);
 P.addXYDataSet("none","none","dash",symbolcolor);
 P.addXYDataPoint(0,mean_value-std_dev); P.addXYDataPoint(mcvalues.size(),mean_value-std_dev);
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createMCBootstrapPlot(const Vector<double>& bootvals, const string& observable_name,
                           double in_mean_value, double in_low, double in_upp, 
                           const string& filename, const string& symbol, 
                           const string& symbolcolor, double rescale, bool drawtoscreen)
{
 GracePlot P("Bootstrap Index","Value");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double low=rescale*in_low;
 double upp=rescale*in_upp;

 P.addXYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<bootvals.size();ind++)
    P.addXYDataPoint(double(ind),rescale*bootvals[ind]);

 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(0,mean_value); P.addXYDataPoint(bootvals.size(),mean_value);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0,low); P.addXYDataPoint(bootvals.size(),low);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0,upp); P.addXYDataPoint(bootvals.size(),upp);
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createMCJackknifePlot(const Vector<double>& jackvals, const string& observable_name,
                           double in_mean_value, const string& filename, const string& symbol, 
                           const string& symbolcolor, double rescale, bool drawtoscreen)
{
 GracePlot P("Jackknife Index","Value");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;

 P.addXYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<jackvals.size();ind++)
    P.addXYDataPoint(double(ind),rescale*jackvals[ind]);

 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(0,mean_value); P.addXYDataPoint(jackvals.size(),mean_value);
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createMCHistogramPlot(const Histogram& histo, const string& observable_name,
                           double in_mean_value, double in_std_dev,
                           const string& filename, const string& barcolor, 
                           double rescale, bool drawtoscreen)
{
 GracePlot P("MC Value","Number");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double std_dev=rescale*in_std_dev;

 P.addBarDataSet(barcolor,"black",rescale*histo.getBarWidth());
 for (uint ind=0;ind<histo.getNumberOfBars();ind++)
    P.addBarDataPoint(rescale*histo.getBarMiddleLocation(ind),histo.getBarHeight(ind));

 double hmax=double(histo.getMaxHeight())*1.3;
 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(mean_value,0.0); P.addXYDataPoint(mean_value,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(mean_value+std_dev,0.0); P.addXYDataPoint(mean_value+std_dev,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(mean_value-std_dev,0.0); P.addXYDataPoint(mean_value-std_dev,hmax);
 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}


void createMCBootstrapHistogramPlot(const Histogram& histo, const string& observable_name,
                                    double in_mean_value, double in_low, double in_upp, 
                                    const string& filename, const string& barcolor, 
                                    double rescale, bool drawtoscreen)
{
 GracePlot P("Bootstrap Value","Number");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;
 double low=rescale*in_low;
 double upp=rescale*in_upp; 

 P.addBarDataSet(barcolor,"black",rescale*histo.getBarWidth());
 for (uint ind=0;ind<histo.getNumberOfBars();ind++)
    P.addBarDataPoint(rescale*histo.getBarMiddleLocation(ind),histo.getBarHeight(ind));

 double hmax=double(histo.getMaxHeight())*1.3;
 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(mean_value,0.0); P.addXYDataPoint(mean_value,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(upp,0.0); P.addXYDataPoint(upp,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(low,0.0); P.addXYDataPoint(low,hmax);
 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createMCJackknifeHistogramPlot(const Histogram& histo, const string& observable_name,
                                    double in_mean_value, 
                                    const string& filename, const string& barcolor, 
                                    double rescale, bool drawtoscreen)
{
 GracePlot P("Jackknife Value","Number");
 P.setFonts("times-roman","times-roman","times-roman","times-roman");
 P.setFontsizes(2.0,1.7,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);
 double mean_value=rescale*in_mean_value;

 P.addBarDataSet(barcolor,"black",rescale*histo.getBarWidth());
 for (uint ind=0;ind<histo.getNumberOfBars();ind++)
    P.addBarDataPoint(rescale*histo.getBarMiddleLocation(ind),histo.getBarHeight(ind));

 double hmax=double(histo.getMaxHeight())*1.3;
 P.addXYDataSet("none","none","solid","black");
 P.addXYDataPoint(mean_value,0.0); P.addXYDataPoint(mean_value,hmax);
 P.addXYDataSet("none","none","dash","black");
 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.25,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}


void createCorrelatorPlot(const vector<XYDYPoint>& corrvals,
                          const ComplexArg& arg,
                          const string& correlator_name,
                          const string& filename, 
                          const string& symbol, 
                          const string& symbolcolor,
                          double rescale, bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t",prefix+" C(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 XYDYPoint cval;
 for (uint ind=0;ind<corrvals.size();ind++){
    cval=corrvals[ind]; cval.yval*=rescale; cval.yerr*=rescale;
    P.addXYDYDataPoint(cval);
    if (corrvals[ind].xval>tmax) tmax=corrvals[ind].xval;}
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0.0,0.0); P.addXYDataPoint(tmax+5.0,0.0);

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createCorrelatorPlot(const vector<XYDYPoint>& corrvals,
                          const ComplexArg& arg,
                          const string& correlator_name,
                          const string& filename, 
                          double verticalmin, double verticalmax,
                          const string& symbol, 
                          const string& symbolcolor,
                          double rescale, bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t",prefix+" C(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.15,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 XYDYPoint cval;
 for (uint ind=0;ind<corrvals.size();ind++){
    cval=corrvals[ind]; cval.yval*=rescale; cval.yerr*=rescale;
    P.addXYDYDataPoint(cval);
    if (corrvals[ind].xval>tmax) tmax=corrvals[ind].xval;}
 P.autoScale(0.02,0.02,0.2,0.2);
 P.setVerticalLimits(verticalmin,verticalmax);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 P.addXYDataSet("none","none","dash","black");
 P.addXYDataPoint(0.0,0.0); P.addXYDataPoint(tmax+5.0,0.0);

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createEffEnergyPlot(const vector<XYDYPoint>& meffvals,
                         const ComplexArg& arg,
                         const string& correlator_name,
                         const string& filename, 
                         const string& symbol, 
                         const string& symbolcolor,
                         bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t","a\\st\\NE\\s\\f{0}eff\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 for (uint ind=0;ind<meffvals.size();ind++){
    P.addXYDYDataPoint(meffvals[ind]);
    if (meffvals[ind].xval>tmax) tmax=meffvals[ind].xval;}
 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}




void createEffEnergyWithFitPlot(const vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                double energy_mean, double energy_err,
                                uint fit_tmin, uint fit_tmax,
                                vector<XYPoint> meff_approach,
                                char goodnesstype, double goodness,
                                const string& correlator_name,
                                const string& filename, 
                                const string& symbol, 
                                const string& symbolcolor,
                                bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t","a\\st\\NE\\s\\f{0}eff\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 for (uint ind=0;ind<meffvals.size();ind++){
    P.addXYDYDataPoint(meffvals[ind]);}

 double fitupper=energy_mean+energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fit_tmin,fitupper);
 P.addXYDataPoint(fit_tmax,fitupper);
 double fitlower=energy_mean-energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fit_tmin,fitlower);
 P.addXYDataPoint(fit_tmax,fitlower);

 if (!(meff_approach.empty())){
    P.addXYDataSet("none","open","dash",symbolcolor);
    P.addXYDataPoints(meff_approach);}

 SimpleMCEstimate fitres(energy_mean,energy_err);
 string fitenergy("\\f{1}a\\st\\NE\\f{}\\sfit\\N = ");
 fitenergy+=fitres.str(2);
 P.addText(fitenergy,0.90,0.85,true,1.7,"black","top-right");

 if (goodnesstype=='Q'){
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.80,true,0,"black","top-right");}
 else if (goodnesstype=='X'){
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.80,true,0,"black","top-right");}

 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}

void createThreePointCorrelatorPlot(const vector<vector<XYDYPoint> >& corrvals,
                            const ComplexArg& arg,
                            const string& plotlabel,
                            const string& filename, 
                            const vector<string>& labels,
                            const vector<string>& symbols, 
                            const vector<string>& symbolcolors,
                            bool drawtoscreen)
{
  string prefix;
  if (arg==RealPart) prefix="\\f{0}Re\\f{}";
  else prefix="\\f{0}Im\\f{}";

  GracePlot P("t\\s\\f{0}ins\\N\\f{} - t\\s\\f{0}sep\\N\\f{}/2", prefix+"R(t\\s\\f{0}sep\\N\\f{}, t\\s\\f{0}ins\\N\\f{})");
  P.setFonts("times-italics", "times-italics", "times-roman", "times-roman");
  P.setFontsizes(2.0, 2.0, 1.5, 1.4);
  P.setView(0.2, 0.95, 0.15, 0.95);
  P.setLegend(0.75, 0.9);

  if ((corrvals.size() != labels.size()) || (corrvals.size() != symbols.size()) ||
      (corrvals.size() != symbolcolors.size())) {
    throw(invalid_argument("unequal number of fits"));
  }

  for (uint si = 0; si < corrvals.size(); ++si) {
    P.addXYDYDataSet(symbols[si], "solid", "none", symbolcolors[si], labels[si]);
    for (uint ind = 0; ind < corrvals[si].size(); ind++) {
      P.addXYDYDataPoint(corrvals[si][ind]);
    }
  }

  P.autoScale(0.02,0.02,0.2,0.2);
  if (!plotlabel.empty())
    P.addText(plotlabel,0.25,0.92,true,0,"black","top-left");

  if (!tidyString(filename).empty()) P.saveToFile(filename);
//  if (drawtoscreen) P.drawToScreen();
}

void createThreePointCorrelatorWithFitPlot(const vector<vector<XYDYPoint> >& corrvals,
                            const ComplexArg& arg,
                            const double b0_mean,
                            const double b0_err,
                            const vector<vector<XYPoint> >& lowerfits,
                            const vector<vector<XYPoint> >& upperfits,
                            char goodnesstype, double goodness,
                            const string& plotlabel,
                            const string& filename, 
                            const vector<string>& labels,
                            const vector<string>& symbols, 
                            const vector<string>& symbolcolors,
                            bool drawtoscreen)
{
  string prefix;
  if (arg==RealPart) prefix="\\f{0}Re\\f{}";
  else prefix="\\f{0}Im\\f{}";

  GracePlot P("t\\s\\f{0}ins\\N\\f{} - t\\s\\f{0}sep\\N\\f{}/2", prefix+"R(t\\s\\f{0}sep\\N\\f{}, t\\s\\f{0}ins\\N\\f{})");
  P.setFonts("times-italics", "times-italics", "times-roman", "times-roman");
  P.setFontsizes(2.0, 2.0, 1.5, 1.4);
  P.setView(0.2, 0.95, 0.15, 0.95);
  P.setLegend(0.75, 0.9);

  if ((corrvals.size() != labels.size()) || (corrvals.size() != symbols.size()) ||
      (corrvals.size() != symbolcolors.size()) || (corrvals.size() != upperfits.size()) ||
      (corrvals.size() != lowerfits.size())) {
    throw(invalid_argument("unequal number of fits"));
  }

  double min_x = -1., max_x = 1.;
  bool set_min_max = false;
  for (uint si = 0; si < corrvals.size(); ++si) {
    for (uint ind = 0; ind < corrvals[si].size(); ind++) {
      if (!set_min_max) {
        min_x = max_x = corrvals[si][ind].xval;
        set_min_max = true;
      }
      if (corrvals[si][ind].xval > max_x) max_x = corrvals[si][ind].xval;
      else if (corrvals[si][ind].xval < min_x) min_x = corrvals[si][ind].xval;
    }
  }

  double fitupper = b0_mean + b0_err;
  double fitlower = b0_mean - b0_err;
  P.addFillDataSet("black","dotted");
  P.addFillDataPoint(min_x, fitupper);
  P.addFillDataPoint(max_x, fitupper);
  P.addFillDataPoint(max_x, fitlower);
  P.addFillDataPoint(min_x, fitlower);

  for (uint si = 0; si < upperfits.size(); ++si) {
    if (upperfits[si].size() == 0) continue;

    P.addFillDataSet(symbolcolors[si], "dotted");
    P.addFillDataPoints(upperfits[si]);
    P.addFillDataPoints(lowerfits[si]);
  }

  for (uint si = 0; si < corrvals.size(); ++si) {
    P.addXYDYDataSet(symbols[si], "solid", "none", symbolcolors[si], labels[si]);
    for (uint ind = 0; ind < corrvals[si].size(); ind++) {
      P.addXYDYDataPoint(corrvals[si][ind]);
    }
  }

  P.addXYDataSet("none", "open", "solid", "black");
  P.addXYDataPoint(min_x, fitupper);
  P.addXYDataPoint(max_x, fitupper);
  P.addXYDataSet("none", "open", "solid", "black");
  P.addXYDataPoint(min_x, fitlower);
  P.addXYDataPoint(max_x, fitlower);


  SimpleMCEstimate fitres(b0_mean, b0_err);
  string fitb0("\\f{1}B\\s0\\N\\f{}\\sfit\\N = ");
  fitb0 += fitres.str(2);
  P.addText(fitb0, 0.25, 0.87, true, 1.7, "black", "top-left");

  if (goodnesstype=='Q') {
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss << goodness; qualstr += ss.str();
    P.addText(qualstr, 0.25, 0.82, true, 0, "black", "top-left");
  }
  else if (goodnesstype=='X') {
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss << goodness; qualstr += ss.str();
    P.addText(qualstr, 0.25, 0.82, true, 0, "black", "top-left");
  }

  P.autoScale(0.02, 0.02, 0.2, 0.2);
  if (!plotlabel.empty())
    P.addText(plotlabel, 0.25, 0.92, true, 0, "black", "top-left");

  if (!tidyString(filename).empty()) P.saveToFile(filename);
//  if (drawtoscreen) P.drawToScreen();
}

void createThreePointCorrelatorWithFitPlot(const vector<XYDYPoint>& summvals,
                            const ComplexArg& arg,
                            const double b0_mean,
                            const double b0_err,
                            const vector<XYPoint>& lowerfit,
                            const vector<XYPoint>& upperfit,
                            char goodnesstype, double goodness,
                            const string& plotlabel,
                            const string& filename, 
                            const string& symbol, 
                            const string& symbolcolor,
                            bool drawtoscreen)
{
  string prefix;
  if (arg==RealPart) prefix="\\f{0}Re\\f{}";
  else prefix="\\f{0}Im\\f{}";

  GracePlot P("t\\s\\f{0}sep\\N\\f{}", prefix+"R\\s\\f{0}summ\\N\\f{}(t\\s\\f{0}sep\\N\\f{})");
  P.setFonts("times-italics", "times-italics", "times-roman", "times-roman");
  P.setFontsizes(2.0, 2.0, 1.5, 1.4);
  P.setView(0.2, 0.95, 0.15, 0.95);
  P.setLegend(0.75, 0.9);


  double min_x = -1., max_x = 1.;
  bool set_min_max = false;
  P.addXYDYDataSet(symbol, "solid", "none", symbolcolor);
  for (uint ind = 0; ind < summvals.size(); ind++) {
    P.addXYDYDataPoint(summvals[ind]);
    if (!set_min_max) {
      min_x = max_x = summvals[ind].xval;
      set_min_max = true;
    }
    if (summvals[ind].xval > max_x) max_x = summvals[ind].xval;
    else if (summvals[ind].xval < min_x) min_x = summvals[ind].xval;
  }

  P.addXYDataSet("none", "none", "solid", symbolcolor);
  P.addXYDataPoints(upperfit);
  P.addXYDataSet("none", "none", "solid", symbolcolor);
  P.addXYDataPoints(lowerfit);

  SimpleMCEstimate fitres(b0_mean, b0_err);
  string fitb0("\\f{1}B\\s0\\N\\f{}\\sfit\\N = ");
  fitb0 += fitres.str(2);
  P.addText(fitb0, 0.25, 0.87, true, 1.7, "black", "top-left");

  if (goodnesstype=='Q') {
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss << goodness; qualstr += ss.str();
    P.addText(qualstr, 0.25, 0.82, true, 0, "black", "top-left");
  }
  else if (goodnesstype=='X') {
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss << goodness; qualstr += ss.str();
    P.addText(qualstr, 0.25, 0.82, true, 0, "black", "top-left");
  }

  P.autoScale(0.02, 0.02, 0.2, 0.2);
  if (!plotlabel.empty())
    P.addText(plotlabel, 0.25, 0.92, true, 0, "black", "top-left");

  if (!tidyString(filename).empty()) P.saveToFile(filename);
//  if (drawtoscreen) P.drawToScreen();
}


void createTminPlot(const vector<vector<XYDYDYPoint> >& goodfits,
                    const vector<vector<XYDYDYPoint> >& badfits,
                    const string& plotlabel,
                    const string& filename, 
                    uint energy_level,
                    const vector<string>& labels,
                    const vector<string>& symbols,
                    const vector<string>& colors,
                    const XYDYDYPoint& chosen_fit,
                    bool drawtoscreen)
{
  GracePlot P("t\\s\\f{0}min\\f{}", "a\\st\\NE\\S("+to_string(energy_level)+")\\N\\s\\f{0}fit\\N\\f{}(t)");
  P.setFonts("times-italics", "times-italics", "times-roman", "times-roman");
  P.setFontsizes(2.0, 2.0, 1.5, 1.4);
  P.setView(0.2, 0.95, 0.15, 0.95);
  P.setLegend(0.65, 0.9);
  double xmax = -10.0;
  double xmin = 1e32;

  if ((goodfits.size() != badfits.size()) || (goodfits.size() != symbols.size()) ||
      (goodfits.size() != colors.size())) {
    throw(invalid_argument("unequal number of fits"));
  }

  for (uint si = 0; si < goodfits.size(); ++si) {
    P.addXYDYDYDataSet(symbols[si], "solid", "none", colors[si], labels[si]);
    for (uint ind=0; ind < goodfits[si].size(); ++ind) {
      double x = goodfits[si][ind].xval;
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      P.addXYDYDYDataPoint(goodfits[si][ind]);
    }

    if (goodfits[si].size() == 0) {
      P.addXYDYDYDataSet(symbols[si], "open", "none", colors[si], labels[si]);
    }
    else {
      P.addXYDYDYDataSet(symbols[si], "open", "none", colors[si]);
    }
    for (uint ind=0; ind < badfits[si].size(); ++ind) {
      double x = badfits[si][ind].xval;
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      P.addXYDYDYDataPoint(badfits[si][ind]);
    }
  }

  if (chosen_fit.xval > 0.0) {
    double chosenupper = chosen_fit.yval + chosen_fit.yuperr;
    P.addXYDataSet("none" ,"open" ,"dash" , "black");
    P.addXYDataPoint(xmin, chosenupper);
    P.addXYDataPoint(xmax, chosenupper);
    double chosen = chosen_fit.yval;
    P.addXYDataSet("none", "solid", "solid", "black");
    P.addXYDataPoint(xmin, chosen);
    P.addXYDataPoint(xmax, chosen);
    double chosenlower = chosen_fit.yval - chosen_fit.ydnerr;
    P.addXYDataSet("none", "open", "dash", "black");
    P.addXYDataPoint(xmin, chosenlower);
    P.addXYDataPoint(xmax, chosenlower);
  }

  P.autoScale(0.02, 0.02, 0.2, 0.2);
  if (!plotlabel.empty())
    P.addText(plotlabel, 0.25, 0.92, true, 0, "black", "top-left");

  if (!tidyString(filename).empty()) P.saveToFile(filename);
//  if (drawtoscreen) P.drawToScreen();
}


void createEffEnergyWithFitAndEnergyRatioPlot(const vector<XYDYPoint>& meffvals,
                                const ComplexArg& arg,
                                double energy_mean, double energy_err,
                                uint fit_tmin, uint fit_tmax,
                                vector<XYPoint> meff_approach,
                                double energy_ratio, double energy_ratio_err,
                                char goodnesstype, double goodness,
                                const string& correlator_name,
                                const string& filename, 
                                const string& symbol, 
                                const string& symbolcolor,
                                bool drawtoscreen)
{
 string prefix;
 if (arg==RealPart) prefix="\\f{0}Re\\f{}";
 else prefix="\\f{0}Im\\f{}";

 GracePlot P("t","a\\st\\NE\\s\\f{0}eff\\N\\f{}(t)");
 P.setFonts("times-italics","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.2,0.95,0.15,0.95);

 P.addXYDYDataSet(symbol,"solid","none",symbolcolor);
 int tmax=0;
 for (uint ind=0;ind<meffvals.size();ind++){
    P.addXYDYDataPoint(meffvals[ind]);
    if (meffvals[ind].xval>tmax) tmax=meffvals[ind].xval;}

 double fitupper=energy_mean+energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fit_tmin,fitupper);
 P.addXYDataPoint(fit_tmax,fitupper);
 double fitlower=energy_mean-energy_err;
 P.addXYDataSet("none","open","solid",symbolcolor);
 P.addXYDataPoint(fit_tmin,fitlower);
 P.addXYDataPoint(fit_tmax,fitlower);

 if (!(meff_approach.empty())){
    P.addXYDataSet("none","open","dash",symbolcolor);
    P.addXYDataPoints(meff_approach);}

 SimpleMCEstimate fitres(energy_mean,energy_err);
 string fitenergy("\\f{1}a\\st\\NE\\f{}\\sfit\\N = ");
 fitenergy+=fitres.str(2);
 P.addText(fitenergy,0.90,0.9,true,1.7,"black","top-right");

 SimpleMCEstimate ratio(energy_ratio,energy_ratio_err);
 string ratiostr("\\f{1}E\\f{}\\sfit\\N\\f{1}/E\\f{}\\sref\\N = ");
 ratiostr+=ratio.str(2);
 P.addText(ratiostr,0.25,0.18,true,1.7,"black","bottom-left");

 if (goodnesstype=='Q'){
    string qualstr("\\f{1}Q\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.85,true,0,"black","top-right");}
 else if (goodnesstype=='X'){
    string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
    stringstream ss; ss.precision(2); ss.setf(ios::fixed);
    ss<<goodness; qualstr+=ss.str();
    P.addText(qualstr,0.90,0.85,true,0,"black","top-right");}

 P.autoScale(0.02,0.02,0.2,0.2);
 if (!correlator_name.empty())
    P.addText(prefix+correlator_name,0.25,0.92,true,0,"black","top-left");

 if (!tidyString(filename).empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen();
}



void createDispersionFitPlot(const vector<XYDYPoint>& energy_sqs,
                             const vector<XYPoint>& lowerfit,
                             const vector<XYPoint>& upperfit,
                             const vector<XYPoint>& disp,
                             const vector<double> param_means,
                             const vector<double> param_errs,
                             char goodnesstype, double goodness,
                             const string& particle_name,
                             const string& filename, 
                             const string& symbol, 
                             const string& symbolcolor,
                             bool show_fit,
                             bool drawtoscreen)
{
  GracePlot P("a\\st\\N\\S2\\Np\\S2\\N", "a\\st\\N\\S2\\NE\\s\\f{0}fit\\N\\S2\\N\\f{}");
  P.setFonts("times-italics", "times-italics", "times-roman", "times-roman");
  P.setFontsizes(2.0, 2.0, 1.5, 1.4);
  P.setView(0.25, 0.95, 0.15, 0.95);

  P.addXYDYDataSet(symbol, "solid", "none", symbolcolor);
  for (uint ind = 0; ind < energy_sqs.size(); ind++) {
    P.addXYDYDataPoint(energy_sqs[ind]);
  } 

  if (show_fit) {
    if (goodnesstype=='Q') {
      string qualstr("\\f{1}Q\\f{} = "); 
      stringstream ss; ss.precision(2); ss.setf(ios::fixed);
      ss << goodness; 
      qualstr += ss.str();
      P.addText(qualstr, 0.90, 0.19, true, 0, "black", "bottom-right");
    }
    else if (goodnesstype=='X') {
      string qualstr("\\xc\\S2\\N/\\f{0}dof\\f{} = "); 
      stringstream ss; ss.precision(2); ss.setf(ios::fixed);
      ss << goodness;
      qualstr += ss.str();
      P.addText(qualstr, 0.90, 0.27, true, 0, "black", "bottom-right");
    }

    for (uint pi = 0; pi < param_means.size(); ++pi) {
      SimpleMCEstimate param_est(param_means[pi], param_errs[pi]);
      string param_str("\\f{1}c\\s" + to_string(pi) + "\\N\\f{} = ");
      param_str += param_est.str(2);
      double y_pos = 0.23 - 0.04*pi;
      P.addText(param_str, 0.90, y_pos, true, 0, "black", "bottom-right");
    }

    P.addXYDataSet("none", "none", "solid", symbolcolor);
    P.addXYDataPoints(upperfit);
    P.addXYDataSet("none", "none", "solid", symbolcolor);
    P.addXYDataPoints(lowerfit);
  }

  P.addXYDataSet("none", "open", "dash", "black");
  P.addXYDataPoints(disp);

  P.autoScale(0.02, 0.02, 0.2, 0.2);
  if (!particle_name.empty()) {
    P.addText(particle_name, 0.30, 0.9, true, 1.5, "black", "top-left");
  }

  if (!tidyString(filename).empty()) P.saveToFile(filename);
//  if (drawtoscreen) P.drawToScreen();
} 




void createCorrMatrixZMagSquaresPlot(const vector<XYDYPoint>& zmag_sqs,
                                     const string& observable_name,
                                     const string& filename, 
                                     const string& barcolor, 
                                     bool drawtoscreen)
{ 
 if (zmag_sqs.empty()) return;
 GracePlot P("Level number \\f{1}n","\\x\\cj\\C\\f{}\\h{-0.2}Z\\S(n)\\h{0.2}\\N\\x\\cj\\C\\h{-0.2}\\S2\\N");
 P.setFonts("times-roman","times-italics","times-roman","times-roman");
 P.setFontsizes(2.0,2.0,1.5,1.4);
 P.setView(0.25,0.95,0.15,0.95);

 double xmin=zmag_sqs[0].xval;
 double xmax=xmin;
 for (uint ind=0;ind<zmag_sqs.size();ind++){
    if (zmag_sqs[ind].xval>xmax) xmax=zmag_sqs[ind].xval;
    if (zmag_sqs[ind].xval<xmin) xmin=zmag_sqs[ind].xval;}
 double barwidth=(xmax-xmin)/zmag_sqs.size();

 P.addBarDataSet(barcolor,"black",barwidth);
 for (uint ind=0;ind<zmag_sqs.size();ind++){
    P.addBarDataPoint(zmag_sqs[ind].xval,zmag_sqs[ind].yval);}

 P.addXYDYDataSet("none","solid","none","black");
 for (uint ind=0;ind<zmag_sqs.size();ind++){
    P.addXYDYDataPoint(zmag_sqs[ind]);}

 P.autoScale(0.02,0.02,0.25,0.0);
 if (!observable_name.empty())
    P.addText(observable_name,0.35,0.92,true,0,"black","top-left");
 if (!filename.empty()) P.saveToFile(filename);
// if (drawtoscreen) P.drawToScreen(); 
}


void createDataFitRatioPlot(const vector<vector<XYDYPoint> >& ratiovals,
                            const string& plotlabel,
                            const string& filename, 
                            const vector<string>& labels,
                            const vector<string>& symbols, 
                            const vector<string>& symbolcolors,
                            bool drawtoscreen)
{
  GracePlot P("t", "C\\S\\f{0}dat\\N\\f{}(t)/C\\S\\f{0}fit\\N\\f{}(t)");
  P.setFonts("times-italics", "times-italics", "times-roman", "times-roman");
  P.setFontsizes(2.0, 2.0, 1.5, 1.4);
  P.setView(0.2, 0.95, 0.15, 0.95);
  P.setLegend(0.65, 0.9);

  if ((ratiovals.size() != labels.size()) || (ratiovals.size() != symbols.size()) ||
      (ratiovals.size() != symbolcolors.size())) {
    throw(invalid_argument("unequal number of fits"));
  }

  for (uint si = 0; si < ratiovals.size(); ++si) {
    P.addXYDYDataSet(symbols[si], "solid", "none", symbolcolors[si], labels[si]);
    for (uint ind = 0; ind < ratiovals[si].size(); ind++) {
      P.addXYDYDataPoint(ratiovals[si][ind]);
    }
  }

  P.autoScale(0.02,0.02,0.2,0.2);
  if (!plotlabel.empty())
    P.addText(plotlabel,0.25,0.92,true,0,"black","top-left");

  if (!tidyString(filename).empty()) P.saveToFile(filename);
//  if (drawtoscreen) P.drawToScreen();
}



        // *****************************************************


string getMomentumName(int xmom, int ymom, int zmom)
{
 return make_string(xmom*xmom+ymom*ymom+zmom*zmom);
}


string getHadronName(const string& flavor, int xmom, int ymom, int zmom,
                     const string& irrep, const string& sptype, uint spid)
{
 string flav;
 if (flavor=="pion") flav="\\xp\\f{1}";
 else if (flavor=="glueball") flav="G";
 else if (flavor=="eta") flav="\\xh\\f{1}";
 else if (flavor=="phi") flav="\\xf\\f{1}";
 else if (flavor=="kaon")  flav="K";              
 else if (flavor=="kbar")  flav="\\oK\\O";            
 else if (flavor=="nucleon") flav="N";
 else if (flavor=="delta") flav="\\xD\\f{1}";
 else if (flavor=="sigma") flav="\\xS\\f{1}";
 else if (flavor=="lambda") flav="\\xL\\f{1}";
 else if (flavor=="xi") flav="\\xX\\f{1}";
 else if (flavor=="omega") flav="\\xW\\f{1}";
 else throw(invalid_argument("Unsupported hadron flavor"));

 string subscript(sptype+make_string(spid));
 if (irrep.length()>subscript.length())
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\s"+subscript+"\\v{}\\z{}\\M{1}\\S"
          +irrep+"\\N\\f{}";
 else
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\S"+irrep+"\\v{}\\z{}\\M{1}\\s"
          +subscript+"\\N\\f{}";
}


string getTetraquarkName(const string& flavor, int xmom, int ymom, int zmom,
                         const string& irrep, const string& sptype, 
                         uint spid, uint colortype)
{
 string flav;
 if (flavor=="isosinglet_eta_eta") flav="tquuuu1";
 else if (flavor=="isotriplet_eta_pion") flav="tquudu3";
 else if (flavor=="isosinglet_pion_pion") flav="tqdudu1";
 else if (flavor=="isotriplet_pion_pion") flav="tqdudu3";
 else if (flavor=="isoquintet_pion_pion") flav="tqdudu5";
 else if (flavor=="isodoublet_kaon_eta") flav="tqsuuu2";
 else if (flavor=="isodoublet_kaon_pion") flav="tqsudu2";
 else if (flavor=="isoquartet_kaon_pion") flav="tqsudu4";
 else if (flavor=="isotriplet_phi_pion") flav="tqssdu3";
 else if (flavor=="isosinglet_eta_phi") flav="tquuss1";
 else if (flavor=="isodoublet_kaon_phi") flav="tqsuss2";
 else if (flavor=="isosinglet_phi_phi") flav="tqssss1";
 else throw(invalid_argument("Unsupported tetraquark flavor"));

 if (colortype==1)
   flav+="p";
 else
   flav+="m";
 string subscript(sptype+make_string(spid));
 if (irrep.length()>subscript.length())
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\s"+subscript+"\\v{}\\z{}\\M{1}\\S"
          +irrep+"\\N\\f{}";//+make_string(colortype);
 else
    return flav+"("+getMomentumName(xmom,ymom,zmom)+")\\m{1}\\S"+irrep+"\\v{}\\z{}\\M{1}\\s"
          +subscript+"\\N\\f{}";//+make_string(colortype);
}


string getOpStandardName(const OperatorInfo& qcd_op)
{
 if (qcd_op.isBasicLapH()){
    BasicLapHOperatorInfo qcdop(qcd_op.getBasicLapH());
    if ((qcdop.isMeson())||(qcdop.isBaryon())||(qcdop.isGlueball())){
       return getHadronName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum(),
                    qcdop.getLGIrrep(),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1));}
    else if ((qcdop.isMesonMeson())||(qcdop.isMesonBaryon())){
       string had1=getHadronName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(1),qcdop.getYMomentum(1),qcdop.getZMomentum(1),
                    qcdop.getLGIrrep(1),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1));
       string had2=getHadronName(qcdop.getFlavor(2),
                    qcdop.getXMomentum(2),qcdop.getYMomentum(2),qcdop.getZMomentum(2),
                    qcdop.getLGIrrep(2),qcdop.getSpatialType(2),
                    qcdop.getSpatialIdNumber(2));
       string iso=qcdop.getIsospin();
       if (iso=="singlet") iso="I=0";
       else if (iso=="doublet") iso="2I=1";
       else if (iso=="triplet") iso="I=1";
       else if (iso=="quartet") iso="2I=3";
       else if (iso=="quintet") iso="I=2";
       else if (iso=="sextet") iso="2I=5";
       else throw(invalid_argument("Unsupported total isospin in getOpStandardName"));
       return string("\\f{0}[")+had1+" "+had2+"\\f{0}]\\f{1}("
              +getMomentumName(qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum())
              +")\\S\\m{2}\\f{1}"+qcdop.getLGIrrep()+"\\N\\M{2}\\s"+iso+"\\f{}\\N";}
    else if (qcdop.isTetraquark()){
       return getTetraquarkName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum(),
                    qcdop.getLGIrrep(),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1), qcdop.getTetraquarkColorType(1));}
    else if (qcdop.isMesonMesonMeson()){
       string had1=getHadronName(qcdop.getFlavor(1),
                    qcdop.getXMomentum(1),qcdop.getYMomentum(1),qcdop.getZMomentum(1),
                    qcdop.getLGIrrep(1),qcdop.getSpatialType(1),
                    qcdop.getSpatialIdNumber(1));
       string had2=getHadronName(qcdop.getFlavor(2),
                    qcdop.getXMomentum(2),qcdop.getYMomentum(2),qcdop.getZMomentum(2),
                    qcdop.getLGIrrep(2),qcdop.getSpatialType(2),
                    qcdop.getSpatialIdNumber(2));
       string had3=getHadronName(qcdop.getFlavor(3),
                    qcdop.getXMomentum(3),qcdop.getYMomentum(3),qcdop.getZMomentum(3),
                    qcdop.getLGIrrep(3),qcdop.getSpatialType(3),
                    qcdop.getSpatialIdNumber(3));
       string iso=qcdop.getIsospin();
       if (iso=="singlet") iso="I=0";
       else if (iso=="doublet") iso="2I=1";
       else if (iso=="triplet") iso="I=1";
       else if (iso=="quartet") iso="2I=3";
       else if (iso=="quintet") iso="I=2";
       else if (iso=="sextet") iso="2I=5";
       else if (iso=="septet") iso="I=3";
       else throw(invalid_argument("Unsupported total isospin in getOpStandardName"));
       return string("\\f{0}[")+had1+" "+had2+" "+had3+"\\f{0}]\\f{1}("
              +getMomentumName(qcdop.getXMomentum(),qcdop.getYMomentum(),qcdop.getZMomentum())
              +")\\S\\m{2}\\f{1}"+qcdop.getLGIrrep()+"\\N\\M{2}\\s"+iso+"\\f{}\\N";}}
 else if (qcd_op.isGenIrrep()){
    GenIrrepOperatorInfo qcdop(qcd_op.getGenIrrep());
    //return qcdop.getIDName()+string(" Level ")+make_string(qcdop.getIDIndex());}
    return qcdop.getIDName()+string(" ")+make_string(qcdop.getIDIndex());}
 else throw(invalid_argument("Unsupported operator type in getOpStandardName"));
 return string("");
}


string getMCObsStandardName(const MCObsInfo& obs)
{
 string label;
 if (obs.isRealPart()) label="Re";
 else label="Im";

 if ((obs.isCorrelatorAtTime())||(obs.isHermitianCorrelatorAtTime())){
    string tstr=string("(")+make_string(obs.getCorrelatorTimeIndex())+"), ";
    OperatorInfo src(obs.getCorrelatorSourceInfo());
    OperatorInfo snk(obs.getCorrelatorSinkInfo());
    if (src==snk) label+=" \\f{1}C\\sAA\\N"+tstr;
    else label+=" \\f{1}C\\sAB\\N"+tstr;
    label+="  \\m{3}A="+getOpStandardName(snk);
    if (src!=snk){
       label+=", ";
       string add("\\f{1}B="); add+=getOpStandardName(src);
       if (label.length()+add.length()<150) label+=add;
       else label+="\\M{3}\\V{-1.7}"+add;}
    label+="\\f{}";}

 else if (obs.isVEV()){
    OperatorInfo qcdop(obs.getVEVInfo());
    return string("\\x\\ca\\C ")+getOpStandardName(qcdop)+" \\x\\cq\\C";}

 return label;
}

string getCorrelatorStandardName(const CorrelatorInfo& corr)
{
 string label;
 string tstr=string("(t), ");
 OperatorInfo src(corr.getSource());
 OperatorInfo snk(corr.getSink());
 if (src==snk) label+=" \\f{1}C\\sAA\\N"+tstr;
 else label+=" \\f{1}C\\sAB\\N"+tstr;
 label+="  \\m{3}A="+getOpStandardName(snk);
 if (src!=snk){
    label+=", ";
    string add("\\f{1}B="); add+=getOpStandardName(src);
    if (label.length()+add.length()<150) label+=add;
    else label+="\\M{3}\\V{-1.7}"+add;}
 label+="\\f{}";
 return label;
}

    //  If "snkname" or "srcname" is "standard", use the standard name,
    //  else use the label given.

string getCorrelatorName(const CorrelatorInfo& corr, 
                         const string& snkname, const string& srcname)
{
 string label;
 string tstr=string("(t), ");
 OperatorInfo src(corr.getSource());
 OperatorInfo snk(corr.getSink());
 if (src==snk) label+=" \\f{1}C\\sAA\\N"+tstr;
 else label+=" \\f{1}C\\sAB\\N"+tstr;
 if (snkname=="standard")
    label+="  \\m{3}A="+getOpStandardName(snk);
 else
    label+="  \\m{3}A="+snkname;
 if (src!=snk){
    label+=", ";
    string add("\\f{1}B=");
    if (srcname=="standard")
       add+=getOpStandardName(src);
    else
       add+=srcname;
    if (label.length()+add.length()<150) label+=add;
    else label+="\\M{3}\\V{-1.7}"+add;}
 label+="\\f{}";
 return label;
}


// *************************************************************
