#include "task_handler.h"
#include "task_utils.h"
#include "correlator_matrix_info.h"
#include "single_pivot.h"
#include "rolling_pivot.h"
#include "create_plots.h"
#include "chisq_tcorr.h"
#include "chisq_fit.h"
#include "pivoter.h"
#include "prior.h"
#include <tuple>
#include <typeinfo>

using namespace std;

struct newlevel_info {
    uint level_insert;
    double new_level_insert_value;
    string corr_name;
    MCEstimate amp;
    MCEstimate energy;
    double new_level_chisqr;
    double bestfit_chisqr;
    double bestfit_chisqrdof;
    bool before;
    bool priored;
};
bool sortkey(const newlevel_info &one, const newlevel_info &two){
    return one.energy.getFullEstimate()<two.energy.getFullEstimate();
}

void TaskHandler::doGEVPCheck(XMLHandler& xml_task, XMLHandler& xml_out, int taskcount)
{ 
    xml_out.set_root("DoGEVPCheck");
    try{
        ArgsHandler xmltask(xml_task);
        uint tmin,tmax;
        xmltask.getUInt("MinTimeSep",tmin);
        xmltask.getUInt("MaxTimeSep",tmax);
        string rotatetype(xmltask.getString("Type")); 
        string color("blue"),symboltype("circle");
        string gevp_data_file(xmltask.getString("GEVPDataFile"));
        string gevp_inserts_datafile(xmltask.getString("GEVPInsertsDataFile"));
        
        uint nlevelinserts = 1;
        xmltask.getOptionalUInt("NLevelInserts",nlevelinserts);

        //get minimizer info
        ChiSquareMinimizerInfo mz_info;  // default minimizer info
        if (xml_task.count_among_children("MinimizerInfo")>0){
            ChiSquareMinimizerInfo mz_user(xml_task);
            mz_info=mz_user;}
        XMLHandler xmlmz;
        mz_info.output(xmlmz);
        xml_out.put_child(xmlmz);
        
        //set uncorrelated if present
        bool uncorrelated=(xml_task.count_among_children("Uncorrelated")>0) ? true: false;
        if (uncorrelated){
        m_obs->setToUnCorrelated();
        // xml_out.put_child("Uncorrelated");
        } else m_obs->setToCorrelated();

        //input rotated correlation matrix info
        Pivot pivoter;
        try{

            LogHelper xmllog;
            bool pkeep = false; //an input?
            ArgsHandler xmlpiv(xmltask,rotatetype+"Initiate");

            pivoter.setType(rotatetype);
            pivoter.initiatePivot(*this,xmlpiv,xmllog,pkeep);
            pivoter.checkInitiate(xmllog,xml_out);

        }catch(const exception& msg){
            throw(invalid_argument("Read pivoter failed."));}
            
        const set<OperatorInfo> opset = pivoter.getOperators();

        uint nops=pivoter.getNumberOfOperators();
        uint nlevels=pivoter.getNumberOfLevels();

        //get time for initial normalization
        string normtime; uint tauN;
        try{xmltask.getString("NormTime",normtime);}
        catch(const exception& msg){xml_out.put_child("Error",string(msg.what()));}
        if(normtime=="tN") tauN = pivoter.getTauN(); 
        else if(normtime=="tZ") tauN = pivoter.getTauZ(); 
        else{
            try{
                tauN = stoull(normtime);
            }catch(const exception& msg){
                tauN = pivoter.getTauZ();
            }
        }
        xml_out.put_child("NormTime",to_string(tauN));

        vector<newlevel_info> newlevels;

        MCObsInfo zmagkey("TempZMagSq",0);
        uint opindex=0;
        uint num_tvals = 100; double dt = double(tmax-tmin+1)/double(num_tvals);
        for (set<OperatorInfo>::const_iterator opit=opset.begin();opit!=opset.end();opit++,opindex++){
            XMLHandler check_xml("GEVPCheck");
            CorrelatorInfo corrinfo(*opit,*opit);
            CorrelatorAtTimeInfo corrNinfo(corrinfo,tauN);

            //use input time to initialize the normalization of the gevp reconstruction
            MCObsInfo normalizekey("TempNormalizeRecon",0);
            for (m_obs->begin();!m_obs->end();++(*m_obs)){   // loop over resamplings
                double energy_sum = 0.0;
                double corrNval = m_obs->getCurrentSamplingValue(corrNinfo);
                for (uint level=0;level<nlevels;level++){
                    zmagkey.resetObsIndex(level*nops+opindex);
                    double Zmag=m_obs->getCurrentSamplingValue(zmagkey);
                    double energy=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level));
                    energy_sum += Zmag*exp(-energy*tauN);
                }   
                m_obs->putCurrentSamplingValue(normalizekey, corrNval/energy_sum);
            }

            //cycle through tmin and fit to the means and use minimum chi^2 to determine range and normalization?

            //fit to higher energies
            // fit to N*(A0*exp(-E0t)+...)+AK*exp(-(EN+dk*dk)t)

            //generate input xml
            XMLHandler xmlf, xmlo;
            xmlf.set_root("TemporalCorrelatorFit");
            opit->output(xmlo);
            xmlf.put_child(xmlo);
            xmlf.put_child("MinimumTimeSeparation",to_string(tmin)); 
            xmlf.put_child("MaximumTimeSeparation",to_string(tmax)); 
            XMLHandler xmlm("Model");
            xmlm.put_child("Type","TimeForwardGEVPReconWithHigherState");
            
            check_xml.put_child(xmlo);

            vector<string> param_names = { //TimeForwardGEVPReconWithHigherState
                "NormalizeGEVP",
                "KAmplitude",
                "deltak" 
            };
            vector<string> new_param_names = { //TimeForwardGEVPReconWithTwoHigherStates
                "NormalizeGEVP",
                "jAmplitude",
                "deltaj",
                "KAmplitude",
                "deltak"
            };
            vector<string> param_names2 = { //TimeForwardHiddenStateSearch
                // "olAmplitude",
                "nlAmplitude",
                "nlDelta"
            };

            // fit, param_name, MCObsInfo
            map<string,map<string,MCObsInfo>> param_keys;

            for(uint i=0;i<param_names.size();i++){
                XMLHandler xmlp(param_names[i]);
                string param_name = param_names[i]+"_"+getOpStandardName(*opit);
                size_t found = param_name.find(" ");
                if (found!=string::npos) param_name.replace(found,1,"_");
                param_keys["TimeForwardGEVPReconWithHigherState"][param_names[i]] = MCObsInfo(param_name,0);
                xmlp.put_child("Name",param_name);
                xmlp.put_child("IDIndex",to_string(0));
                xmlm.put_child(xmlp);
            }
            m_obs->begin();
            XMLHandler xmlni("NormalizeInit");
            XMLHandler xmlni0("MCObsInfo");
            xmlni0.put_child("ObsName",normalizekey.getObsName());
            xmlni0.put_child("Index",to_string(normalizekey.getObsIndex()));
            xmlni.put_child(xmlni0);
            xmlni.put_child("Mean",to_string(m_obs->getCurrentSamplingValue(normalizekey)));
            xmlm.put_child(xmlni);
            for (uint level=0;level<nlevels;level++){
                XMLHandler xmla("Amplitude"+to_string(level));
                XMLHandler xmla0("MCObsInfo");
                xmla0.put_child("ObsName","TempZMagSq");
                xmla0.put_child("Index",to_string(level*nops+opindex));
                xmla.put_child(xmla0);
                zmagkey.resetObsIndex(level*nops+opindex);
                xmla.put_child("Mean",to_string(m_obs->getCurrentSamplingValue(zmagkey)));
                xmlm.put_child(xmla);

                MCObsInfo energy_key = pivoter.getEnergyKey(level);
                XMLHandler xmle("Energy"+to_string(level));
                XMLHandler xmle0("MCObsInfo");
                xmle0.put_child("ObsName",energy_key.getObsName());
                xmle0.put_child("Index",to_string(energy_key.getObsIndex()));
                xmle.put_child(xmle0);
                xmle.put_child("Mean",to_string(m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level))));
                xmlm.put_child(xmle);
            }
            xmlf.put_child(xmlm);

            //set up 2 state fit xml
            MCObsInfo AkInit("AkInit",0), dkInit("dkInit",0);
            XMLHandler xmlf2state(xmlf, XMLHandler::copy);
            XMLHandler xmlm2state(xmlf2state,"Model");
            XMLHandler xmlak("AkInit");
            XMLHandler xmlak0("MCObsInfo");
            xmlak0.put_child("ObsName","AkInit");
            xmlak0.put_child("Index",to_string(0));
            xmlak.put_child(xmlak0);
            xmlm2state.put_child(xmlak);
            XMLHandler xmldk("dkInit");
            XMLHandler xmldk0("MCObsInfo");
            xmldk0.put_child("ObsName","dkInit");
            xmldk0.put_child("Index",to_string(0));
            xmldk.put_child(xmldk0);
            xmlm2state.put_child(xmldk);
            zmagkey.resetObsIndex((nlevels-1)*nops+opindex);
            for (m_obs->begin();!m_obs->end();++(*m_obs)){
                double Ak = 0.5*m_obs->getCurrentSamplingValue(normalizekey)*m_obs->getCurrentSamplingValue(zmagkey);
                m_obs->putCurrentSamplingValue(AkInit, Ak);
                double dk = sqrt(m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(nlevels-1)));
                m_obs->putCurrentSamplingValue(dkInit, dk);
            }
            for(uint i=0;i<param_names.size();i++){
                xmlm2state.seek_unique(param_names[i]);
                xmlm2state.erase_current_element();
            }
            for(uint i=0;i<new_param_names.size();i++){
                XMLHandler xmlp(new_param_names[i]);
                string param_name = "two_"+new_param_names[i]+"_"+getOpStandardName(*opit);
                size_t found = param_name.find(" ");
                if (found!=string::npos) param_name.replace(found,1,"_");
                param_keys["TimeForwardGEVPReconWithTwoHigherStates"][new_param_names[i]] = MCObsInfo(param_name,0);
                xmlp.put_child("Name",param_name);
                xmlp.put_child("IDIndex",to_string(0));
                xmlm2state.put_child(xmlp);
            }
            xmlm2state.seek_unique("Type");
            xmlm2state.seek_next_node();       
            xmlm2state.set_text_content("TimeForwardGEVPReconWithTwoHigherStates"); 

            //fit higher state
            bool fit_success = false; bool two_states = false;
            double best_chisqdof = 200, good_chisqdof = 1.4; double best_chisq; uint best_tmin;
            for( uint tmin_vary=tmin; tmin_vary<tmin+8; tmin_vary++){ //make input
                string this_fit = "TimeForwardGEVPReconWithHigherState";
                for (map<string,MCObsInfo>::const_iterator rt=param_keys[this_fit].begin(); rt!=param_keys[this_fit].end();rt++) 
                     m_obs->eraseSamplings(rt->second);
                XMLHandler xmlfitlog("FitLog");
                xmlf.seek_unique("MinimumTimeSeparation");
                xmlf.seek_next_node(); 
                xmlf.set_text_content(make_string(tmin_vary)); 

                RealTemporalCorrelatorFit RTC(xmlf,*m_obs,taskcount);
                double chisq_dof,qual; vector<MCEstimate> bestfit_params;
                XMLHandler test("FitTest");
                test.put_child("Tmin",to_string(tmin_vary));
                try{
                    doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlfitlog);
                    test.put_child("ChiSqDOF",to_string(chisq_dof));

                    //put ak and dk into ak and dk init
                    if(chisq_dof<best_chisqdof){
                        best_chisqdof=chisq_dof; best_tmin=tmin_vary; fit_success=true; two_states=false;
                        for (m_obs->begin();!m_obs->end();++(*m_obs)){
                            double Ak = m_obs->getCurrentSamplingValue(param_keys[this_fit]["KAmplitude"]);
                            m_obs->putCurrentSamplingValue(AkInit, Ak);
                            double dk = m_obs->getCurrentSamplingValue(param_keys[this_fit]["deltak"]);
                            m_obs->putCurrentSamplingValue(dkInit, dk);
                        }
                    }
                    if(chisq_dof<=good_chisqdof){ check_xml.put_child(test); break; }
                }catch(const exception& msg){
                        test.put_child("Error","Fit fail: "+string(msg.what()));}
                check_xml.put_child(test);

            }
            // uint one_exp_best_tmin = best_tmin;
            for( uint tmin_vary=tmin; tmin_vary<tmin+7; tmin_vary++){
                string this_fit = "TimeForwardGEVPReconWithTwoHigherStates";
                for (map<string,MCObsInfo>::const_iterator rt=param_keys[this_fit].begin(); rt!=param_keys[this_fit].end();rt++) 
                     m_obs->eraseSamplings(rt->second);
                XMLHandler xmlfitlog2("FitLog");
                xmlf2state.seek_unique("MinimumTimeSeparation");
                xmlf2state.seek_next_node();       
                xmlf2state.set_text_content(make_string(tmin_vary)); 

                RealTemporalCorrelatorFit RTC2(xmlf2state,*m_obs,taskcount);
                double chisq_dof2,qual2; vector<MCEstimate> bestfit_params2;
                XMLHandler test2("FitTest2");
                test2.put_child("Tmin",to_string(tmin_vary));
                try{
                    doChiSquareFitting(RTC2,mz_info,chisq_dof2,qual2,bestfit_params2,xmlfitlog2);
                    // test2.put_child(xmlfitlog2);
                    test2.put_child("ChiSqDOF",to_string(chisq_dof2));
                    if(chisq_dof2<best_chisqdof){
                        best_chisqdof=chisq_dof2; best_tmin=tmin_vary; fit_success=true; two_states=true;
                        if(chisq_dof2<=good_chisqdof || tmin_vary==tmax){ check_xml.put_child(test2); break; }
                        for (m_obs->begin();!m_obs->end();++(*m_obs)){
                            double Ak = m_obs->getCurrentSamplingValue(param_keys[this_fit]["KAmplitude"]);
                            m_obs->putCurrentSamplingValue(AkInit, Ak);
                            double dk = m_obs->getCurrentSamplingValue(param_keys[this_fit]["deltak"]);
                            m_obs->putCurrentSamplingValue(dkInit, dk);
                        }
                    }
                }catch(const exception& msg){
                        // test2.put_child(xmlfitlog2);
                        test2.put_child("Error","Fit fail: "+string(msg.what()));}
                check_xml.put_child(test2);

                if(tmin_vary==tmax) break;
            }

            string this_fit = "TimeForwardGEVPReconWithHigherState";
            if(two_states) this_fit = "TimeForwardGEVPReconWithTwoHigherStates";
            if(fit_success){
                for (map<string,MCObsInfo>::const_iterator rt=param_keys[this_fit].begin(); rt!=param_keys[this_fit].end();rt++){ 
                     m_obs->eraseSamplings(rt->second);
                }
                XMLHandler xmlfitlog("HigherStateFitLog");
                XMLHandler bestfit_xml;
                if(two_states) bestfit_xml=xmlf2state;
                else bestfit_xml=xmlf;
                bestfit_xml.seek_unique("MinimumTimeSeparation");
                bestfit_xml.seek_next_node();       
                bestfit_xml.set_text_content(make_string(best_tmin)); 
                check_xml.put_child(bestfit_xml);
                RealTemporalCorrelatorFit RTC(bestfit_xml,*m_obs,taskcount);
                double chisq_dof,qual; vector<MCEstimate> bestfit_params;
                try{
                    doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlfitlog);
                } catch(const exception& msg){
                    check_xml.put_child(xmlfitlog);
                    xml_out.put_child(check_xml);
                    throw(invalid_argument("Failed fit reconstruction: "+string(msg.what())));
                }
                check_xml.put_child(xmlfitlog);

                double dof=double(RTC.getNumberOfObervables()-RTC.getNumberOfParams());
                best_chisq = chisq_dof*dof;
                check_xml.put_child("BestFitChiSqrDOF",to_string(chisq_dof));
                check_xml.put_child("BestFitChiSqr",to_string(best_chisq));
                check_xml.put_child("BestFitTmin",to_string(best_tmin));
                check_xml.put_child("BestFitTmax",to_string(tmax));

                vector<int> try_before = {0,1}; 
                vector<int> try_priors = {0,1}; 
                string hidden_fit = "TimeForwardHiddenStateSearch";
                double energy_low, energy_high, dE;
                energy_low = m_obs->getEstimate(pivoter.getEnergyKey(0)).getFullEstimate();
                energy_high = m_obs->getEstimate(pivoter.getEnergyKey(nlevels-1)).getFullEstimate();
                dE = (energy_high-energy_low)/double(nlevelinserts+1);


                for (uint il=1;il<=nlevelinserts;il++){
                    double new_level_insert = energy_low+double(il)*dE;
                    m_obs->begin();
                    uint level;
                    double below, above;
                    for (level=0;level<nlevels-2;level++){
                        below=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level));
                        above=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level+1));
                        if( (new_level_insert>below) && (new_level_insert<above) ) break;
                    }
                    double initial_dk = sqrt(new_level_insert-below);

                    for(uint ip=0;ip<try_priors.size();ip++){
                        for( uint ib=0;ib<try_before.size();ib++){
                            XMLHandler xmlf2(xmlf,XMLHandler::copy);
                            m_obs->begin();
                            XMLHandler xmlfitlog2("NewStateFitLog");
                            XMLHandler xmlm2(xmlf2,"Model");
                            bool fit_insert_success = false;

                            for (map<string,MCObsInfo>::const_iterator rt=param_keys[hidden_fit].begin(); rt!=param_keys[hidden_fit].end();rt++){
                                m_obs->eraseSamplings(rt->second);
                            } 
                            XMLHandler test("NewStateFit");
                            test.put_child("Level",to_string(level)); 
                            if(try_before[ib]){ 
                                xmlm2.put_child("Before", "true"); 
                                test.put_child("Before", "1");
                            } else {
                                xmlm2.put_child("Before", "false"); 
                            }
                            xmlm2.put_child("LevelInsertInitialize",to_string(initial_dk));

                            //delete old param input
                            for(uint i=0;i<param_names.size();i++){
                                xmlm2.seek_unique(param_names[i]);
                                xmlm2.erase_current_element();
                            } 
                            xmlm2.seek_root();

                            //set new constants input
                            m_obs->begin();
                            for (map<string,MCObsInfo>::const_iterator rt=param_keys[this_fit].begin(); rt!=param_keys[this_fit].end();rt++){
                                XMLHandler xmlp(rt->first);
                                XMLHandler xmlp0("MCObsInfo");
                                xmlp0.put_child("ObsName",rt->second.getObsName());
                                xmlp0.put_child("Index",to_string(rt->second.getObsIndex()));
                                xmlp.put_child(xmlp0);
                                xmlp.put_child("Mean",to_string(m_obs->getCurrentSamplingValue(rt->second)));
                                xmlm2.put_child(xmlp);
                            }
                            for(uint i=0;i<param_names2.size();i++){
                                XMLHandler xmlp(param_names2[i]);
                                string param_name = param_names2[i]+"_"+getOpStandardName(*opit);
                                size_t found = param_name.find(" ");
                                if (found!=string::npos) param_name.replace(found,1,"_");
                                param_keys[hidden_fit][param_names2[i]] = MCObsInfo(param_name,0);
                                xmlp.put_child("Name",param_name);
                                xmlp.put_child("IDIndex",to_string(0));
                                xmlm2.put_child(xmlp);
                            }
                            xmlm2.put_child("LevelInsertIndex",to_string(level));

                            xmlm2.seek_unique("Type");
                            xmlm2.seek_next_node();       
                            xmlm2.set_text_content("TimeForwardHiddenStateSearch"); 

                            if(try_priors[ip]){
                                XMLHandler priors("Priors");

                                XMLHandler nlDelta("nlDelta");
                                double level_below = m_obs->getEstimate(pivoter.getEnergyKey(level)).getFullEstimate();
                                double level_above = m_obs->getEstimate(pivoter.getEnergyKey(level+1)).getFullEstimate();
                                nlDelta.put_child("Mean",to_string(sqrt(new_level_insert-below)));
                                // nlDelta.put_child("Mean",to_string(sqrt((level_above-level_below)/2.0)));
                                double err = min(sqrt(new_level_insert-below),sqrt(above-new_level_insert));
                                nlDelta.put_child("Error",to_string(err)); //update later
                                // nlDelta.put_child("Error",to_string(sqrt(new_level_insert-below))); //update later
                                // nlDelta.put_child("Error",to_string(sqrt((level_above-level_below)/2.0))); //update later
                                priors.put_child(nlDelta);

                                XMLHandler nlamp("nlAmplitude");
                                if(try_before[ib]) zmagkey.resetObsIndex(level*nops+opindex);
                                else zmagkey.resetObsIndex((level+1)*nops+opindex);
                                double amp = m_obs->getEstimate(zmagkey).getFullEstimate();
                                nlamp.put_child("Mean",to_string(amp*0.5));
                                nlamp.put_child("Error",to_string(amp*0.5)); //update later
                                priors.put_child(nlamp);

                                xmlf2.put_child(priors);
                                test.put_child(priors);
                            }

                            double chisq_dof2,qual2; vector<MCEstimate> bestfit_params2;
                            try{ 
                                RealTemporalCorrelatorFit RTC2(xmlf2,*m_obs,taskcount);
                                doChiSquareFitting(RTC2,mz_info,chisq_dof2,qual2,bestfit_params2,xmlfitlog2);
                                MCObsInfo new_energy_key = MCObsInfo("NewEnergy",0);
                                test.put_child("AmplitudeValue",to_string(bestfit_params2[0].getFullEstimate()));
                                test.put_child("AmplitudeError",to_string(bestfit_params2[0].getSymmetricError()));
                                test.put_child("NewLevelValue",to_string(bestfit_params2[1].getFullEstimate()));
                                test.put_child("NewLevelError",to_string(bestfit_params2[1].getSymmetricError()));
                                for (m_obs->begin();!m_obs->end();++(*m_obs)){
                                    double energy=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level));
                                    double dk = m_obs->getCurrentSamplingValue(param_keys[hidden_fit]["nlDelta"]);
                                    m_obs->putCurrentSamplingValue(new_energy_key,energy+dk*dk);
                                }
                                MCEstimate new_energy = m_obs->getEstimate(new_energy_key);
                                test.put_child("NewLevelEnergy",to_string(new_energy.getFullEstimate()));
                                dof=double(RTC2.getNumberOfObervables()-RTC2.getNumberOfParams());
                                test.put_child("ChiSqr",to_string(chisq_dof2*dof));
                                if(bestfit_params2[0].getFullEstimate()>0.0) // if new level is not zero and above spectrum?
                                if( (chisq_dof2*dof<best_chisq*2.0) && (new_energy.getFullEstimate()<m_obs->getEstimate(pivoter.getEnergyKey(nlevels-1)).getFullEstimate()) )
                                    if(bestfit_params2[0].getFullEstimate()>bestfit_params2[0].getSymmetricError()){
                                        double ampsum = 0.0;
                                        for (uint level2=0;level2<nlevels;level2++){
                                            zmagkey.resetObsIndex(level2*nops+opindex);
                                            ampsum+=m_obs->getEstimate(zmagkey).getFullEstimate();
                                        }
                                        if(bestfit_params2[0].getFullEstimate()<2.0*ampsum){
                                            newlevels.push_back(newlevel_info{level,new_level_insert,make_string(opindex),bestfit_params2[0],new_energy,chisq_dof2*dof,
                                                                                best_chisq,chisq_dof, try_before[ib], try_priors[ip]});
                                            fit_insert_success = true;
                                        }
                                    }
                                m_obs->eraseSamplings(new_energy_key);
                                // test.put_child(xmlf2);
                                // test.put_child(xmlfitlog2); 
                                
                            }catch(const exception& msg){ 
                                // test.put_child(xmlf2);
                                test.put_child("Error",string(msg.what()));
                                // test.put_child(xmlfitlog2); 
                            }

                            if(fit_insert_success){
                                //get reconstruction estimates
                                map<double,MCEstimate> corrtrecon;
                                map<double,MCEstimate> efftrecon;
                                for (double tval=tmin;tval<=double(tmax);tval+=dt){ //finer t
                                    MCObsInfo obskey("TempCorrRecon",0);
                                    MCObsInfo obskey2("TempEffRecon",0);
                                    //off
                                    for (m_obs->begin();!m_obs->end();++(*m_obs)){   // loop over resamplings
                                        double energy_sum = 0.0;
                                        double eff_energy_sum = 0.0;
                                        // double Aol = m_obs->getCurrentSamplingValue(param_keys[param_names.size()]);
                                        double Ak = m_obs->getCurrentSamplingValue(param_keys[hidden_fit]["nlAmplitude"]);
                                        double dk = m_obs->getCurrentSamplingValue(param_keys[hidden_fit]["nlDelta"]);

                                        for (uint level2=0;level2<nlevels;level2++){
                                            zmagkey.resetObsIndex(level2*nops+opindex);
                                            double Zmag=m_obs->getCurrentSamplingValue(zmagkey);
                                            double energy=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level2));
                                            if( (level2==level) && (try_before[ib]) ) energy_sum+=0.5*Zmag*exp(-energy*tval)*(1.0+Ak*exp(-dk*dk*tval));
                                            else if( (level2==level) ) energy_sum+=Zmag*exp(-energy*tval)*(1.0+0.5*Ak*exp(-dk*dk*tval));
                                            else if( (level2==level+1) && (try_before[ib]) ) energy_sum+=Zmag*exp(-(energy)*tval); //+dk*dk
                                            else if( (level2==level+1)  ) energy_sum+=0.5*Zmag*exp(-(energy)*tval); //+dk*dk
                                            else energy_sum += Zmag*exp(-energy*tval);

                                            if( (level2==level) && (try_before[ib]) ){
                                                eff_energy_sum+=energy*0.5*Zmag*exp(-energy*tval);
                                                eff_energy_sum+=(energy+dk*dk)*0.25*Zmag*Ak*exp(-(energy+dk*dk)*tval);
                                            } else if (level2==level) {
                                                eff_energy_sum+=energy*Zmag*exp(-energy*tval);
                                                eff_energy_sum+=(energy+dk*dk)*0.5*Zmag*Ak*exp(-(energy+dk*dk)*tval);
                                            }
                                            else if( (level2==level+1) && (!try_before[ib]) ) eff_energy_sum+=energy*0.5*Zmag*exp(-(energy)*tval); //+dk*dk
                                            else eff_energy_sum += energy*Zmag*exp(-energy*tval);
                                            
                                        }   
                                        double fit_N = m_obs->getCurrentSamplingValue(param_keys[this_fit]["NormalizeGEVP"]);
                                        Ak = m_obs->getCurrentSamplingValue(param_keys[this_fit]["KAmplitude"]);
                                        dk = m_obs->getCurrentSamplingValue(param_keys[this_fit]["deltak"]);
                                        double Aj, dj, this_corrval, this_ndcorrval; 
                                        double energy=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(nlevels-1));
                                        if(two_states){
                                            Aj = m_obs->getCurrentSamplingValue(param_keys[this_fit]["jAmplitude"]);
                                            dj = m_obs->getCurrentSamplingValue(param_keys[this_fit]["deltaj"]);
                                            this_corrval = fit_N*energy_sum+Aj*exp(-(energy+dj*dj)*tval)+Ak*exp(-(energy+dj*dj+dk*dk)*tval);
                                            this_ndcorrval = fit_N*eff_energy_sum+Aj*(energy+dj*dj)*exp(-(energy+dj*dj)*tval)
                                                                +Ak*(energy+dj*dj+dk*dk)*exp(-(energy+dj*dj+dk*dk)*tval);
                                        }else{
                                            this_corrval = fit_N*energy_sum+Ak*exp(-(energy+dk*dk)*tval);
                                            this_ndcorrval = fit_N*eff_energy_sum+Ak*(energy+dk*dk)*exp(-(energy+dk*dk)*tval );
                                        }
                                        m_obs->putCurrentSamplingValue(obskey, this_corrval);
                                        m_obs->putCurrentSamplingValue(obskey2, this_ndcorrval/this_corrval);
                                    }
                                    //get estimates
                                    corrtrecon.insert( pair<double,MCEstimate>{tval, m_obs->getEstimate(obskey)} );
                                    efftrecon.insert( pair<double,MCEstimate>{tval, m_obs->getEstimate(obskey2)} );

                                    //erase samplings
                                    m_obs->eraseSamplings(obskey);
                                    m_obs->eraseSamplings(obskey2);
                                }
                                //plot corr
                                map<double,MCEstimate> results;
                                getCorrelatorEstimates(m_obs,corrinfo,true,false,RealPart,m_obs->getCurrentSamplingMode(),results);
                                uint type = 0;
                                string plotfilestub(xmltask.getString("PlotFileStub"));
                                vector<XYDYPoint> corrvals(results.size());
                                vector<XYPoint> corrreconvals(corrtrecon.size());
                                vector<XYPoint> corrreconvals_below(corrtrecon.size());
                                vector<XYPoint> corrreconvals_above(corrtrecon.size());
                                uint k=0;
                                for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
                                corrvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                                    (rt->second).getSymmetricError());}

                                k=0;
                                for (map<double,MCEstimate>::const_iterator rt=corrtrecon.begin();rt!=corrtrecon.end();rt++,k++){
                                    corrreconvals[k]=XYPoint(rt->first, (rt->second).getFullEstimate());
                                    corrreconvals_below[k]=XYPoint(rt->first, (rt->second).getFullEstimate()-(rt->second).getSymmetricError());
                                    corrreconvals_above[k]=XYPoint(rt->first, (rt->second).getFullEstimate()+(rt->second).getSymmetricError());}

                                string corrname("Corr");
                                try{corrname=getCorrelatorStandardName(corrinfo);}
                                catch(const exception& xp){}
                                string plotfile(plotfilestub+"_corr2_"+make_string(opindex)+"_levelInsert"+make_string(level)+".agr");
                                createCorrelatorPlotWithRecon(corrvals,corrreconvals,corrreconvals_below,corrreconvals_above,
                                                                RealPart,corrname,plotfile,symboltype,color,type);
                                test.put_child("CorrPlotFile",plotfile);

                                //plot energy
                                type = 1;
                                getEffectiveEnergy(m_obs,corrinfo,true,false,RealPart,m_obs->getCurrentSamplingMode(), 1, 0,results);
                                corrvals.resize(results.size());
                                corrreconvals.resize(efftrecon.size());
                                corrreconvals_below.resize(efftrecon.size());
                                corrreconvals_above.resize(efftrecon.size());
                                k=0;
                                for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
                                corrvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                                    (rt->second).getSymmetricError());}

                                k=0;
                                for (map<double,MCEstimate>::const_iterator rt=efftrecon.begin();rt!=efftrecon.end();rt++,k++){
                                    corrreconvals[k]=XYPoint(rt->first, (rt->second).getFullEstimate());
                                    corrreconvals_below[k]=XYPoint(rt->first, (rt->second).getFullEstimate()-(rt->second).getSymmetricError());
                                    corrreconvals_above[k]=XYPoint(rt->first, (rt->second).getFullEstimate()+(rt->second).getSymmetricError());}

                                plotfile = plotfilestub+"_eff2_"+make_string(opindex)+"_levelInsert"+make_string(level)+".agr";
                                createCorrelatorPlotWithRecon(corrvals,corrreconvals,corrreconvals_below,corrreconvals_above,
                                                                RealPart,corrname,plotfile,symboltype,color,type);
                                test.put_child("EffPlotFile",plotfile);
                            }

                            check_xml.put_child(test);
                        }
                    }
                }
            }
            
            // get reconstruction estimates
            map<double,MCEstimate> corrtrecon;
            map<double,MCEstimate> efftrecon;
            map<double,MCEstimate> corrtrecon2;
            map<double,MCEstimate> efftrecon2;
            for (double tval=tmin;tval<=double(tmax);tval+=dt){ //finer t
                MCObsInfo obskey("TempCorrRecon",0);
                MCObsInfo obskey2("TempEffRecon",0);
                MCObsInfo obskey3("TempCorrRecon",1);
                MCObsInfo obskey4("TempEffRecon",1);
                for (m_obs->begin();!m_obs->end();++(*m_obs)){   // loop over resamplings
                    double energy_sum = 0.0;
                    double eff_energy_sum = 0.0;
                    for (uint level=0;level<nlevels;level++){
                        zmagkey.resetObsIndex(level*nops+opindex);
                        double Zmag=m_obs->getCurrentSamplingValue(zmagkey);
                        double energy=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(level));
                        energy_sum += Zmag*exp(-energy*tval);
                        eff_energy_sum += energy*Zmag*exp(-energy*tval);
                    }   
                    double corrNval = m_obs->getCurrentSamplingValue(normalizekey);
                    m_obs->putCurrentSamplingValue(obskey, corrNval*energy_sum);
                    m_obs->putCurrentSamplingValue(obskey2, eff_energy_sum/energy_sum);

                    if(fit_success){
                        double fit_N = m_obs->getCurrentSamplingValue(param_keys[this_fit]["NormalizeGEVP"]);
                        double Ak = m_obs->getCurrentSamplingValue(param_keys[this_fit]["KAmplitude"]);
                        double dk = m_obs->getCurrentSamplingValue(param_keys[this_fit]["deltak"]);
                        double energy=m_obs->getCurrentSamplingValue(pivoter.getEnergyKey(nlevels-1));
                        double Aj, dj, this_corrval, this_ndcorrval; 
                        if(two_states){
                            Aj = m_obs->getCurrentSamplingValue(param_keys[this_fit]["jAmplitude"]);
                            dj = m_obs->getCurrentSamplingValue(param_keys[this_fit]["deltaj"]);
                            this_corrval = fit_N*energy_sum+Aj*exp(-(energy+dj*dj)*tval)+Ak*exp(-(energy+dj*dj+dk*dk)*tval);
                            this_ndcorrval = fit_N*eff_energy_sum+Aj*(energy+dj*dj)*exp(-(energy+dj*dj)*tval)
                                                +Ak*(energy+dj*dj+dk*dk)*exp(-(energy+dj*dj+dk*dk)*tval);
                        }else{
                            this_corrval = fit_N*energy_sum+Ak*exp(-(energy+dk*dk)*tval);
                            this_ndcorrval = fit_N*eff_energy_sum+Ak*(energy+dk*dk)*exp(-(energy+dk*dk)*tval );
                        }
                        m_obs->putCurrentSamplingValue(obskey3, this_corrval);
                        m_obs->putCurrentSamplingValue(obskey4, this_ndcorrval/this_corrval);
                    }
                }
                //get estimates
                corrtrecon.insert( pair<double,MCEstimate>{tval, m_obs->getEstimate(obskey)} );
                efftrecon.insert( pair<double,MCEstimate>{tval, m_obs->getEstimate(obskey2)} );
                if(fit_success){
                    corrtrecon2.insert( pair<double,MCEstimate>{tval, m_obs->getEstimate(obskey3)} );
                    efftrecon2.insert( pair<double,MCEstimate>{tval, m_obs->getEstimate(obskey4)} );
                }

                //erase samplings
                m_obs->eraseSamplings(obskey);
                m_obs->eraseSamplings(obskey2);
                if(fit_success){
                    m_obs->eraseSamplings(obskey3);
                    m_obs->eraseSamplings(obskey4);
                }
            }
            
            //erase samplings
            m_obs->eraseSamplings(normalizekey);

            if(fit_success)
                for (map<string,map<string,MCObsInfo>>::const_iterator rt1=param_keys.begin(); rt1!=param_keys.end();rt1++)
                    for (map<string,MCObsInfo>::const_iterator rt2=rt1->second.begin(); rt2!=rt1->second.end();rt2++)
                        m_obs->eraseSamplings(rt2->second);

            //get raw correlator
            map<double,MCEstimate> results;
            getCorrelatorEstimates(m_obs,corrinfo,true,false,RealPart,m_obs->getCurrentSamplingMode(),results);

            // XMLHandler xmloutfiles()//get opname

            //plot corr
            uint type = 0;
            string plotfilestub(xmltask.getString("PlotFileStub"));
            vector<XYDYPoint> corrvals(results.size());
            vector<XYPoint> corrreconvals(corrtrecon.size());
            vector<XYPoint> corrreconvals_below(corrtrecon.size());
            vector<XYPoint> corrreconvals_above(corrtrecon.size());
            uint k=0;
            for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
            corrvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                (rt->second).getSymmetricError());}

            k=0;
            for (map<double,MCEstimate>::const_iterator rt=corrtrecon.begin();rt!=corrtrecon.end();rt++,k++){
                corrreconvals[k]=XYPoint(rt->first, (rt->second).getFullEstimate());
                corrreconvals_below[k]=XYPoint(rt->first, (rt->second).getFullEstimate()-(rt->second).getSymmetricError());
                corrreconvals_above[k]=XYPoint(rt->first, (rt->second).getFullEstimate()+(rt->second).getSymmetricError());}

            string corrname("Corr");
            try{corrname=getCorrelatorStandardName(corrinfo);}
            catch(const exception& xp){}
            string plotfile(plotfilestub+"_corr_"+make_string(opindex)+".agr");
            createCorrelatorPlotWithRecon(corrvals,corrreconvals,corrreconvals_below,corrreconvals_above,
                                            RealPart,corrname,plotfile,symboltype,color,type);
            check_xml.put_child("CorrPlotFile",plotfile);

            //plot corr with extra state
            if(fit_success){
                corrreconvals.resize(corrtrecon2.size());
                corrreconvals_below.resize(corrtrecon2.size());
                corrreconvals_above.resize(corrtrecon2.size());

                k=0;
                for (map<double,MCEstimate>::const_iterator rt=corrtrecon2.begin();rt!=corrtrecon2.end();rt++,k++){
                    corrreconvals[k]=XYPoint(rt->first, (rt->second).getFullEstimate());
                    corrreconvals_below[k]=XYPoint(rt->first, (rt->second).getFullEstimate()-(rt->second).getSymmetricError());
                    corrreconvals_above[k]=XYPoint(rt->first, (rt->second).getFullEstimate()+(rt->second).getSymmetricError());}

                plotfile = plotfilestub+"_corr2_"+make_string(opindex)+".agr";
                createCorrelatorPlotWithRecon(corrvals,corrreconvals,corrreconvals_below,corrreconvals_above,
                                                RealPart,corrname,plotfile,symboltype,color,type);
                check_xml.put_child("CorrPlotFile2",plotfile);
            }

            //plot energy
            type = 1;
            getEffectiveEnergy(m_obs,corrinfo,true,false,RealPart,m_obs->getCurrentSamplingMode(), 1, 0,results);
            corrvals.resize(results.size());
            corrreconvals.resize(efftrecon.size());
            corrreconvals_below.resize(efftrecon.size());
            corrreconvals_above.resize(efftrecon.size());
            k=0;
            for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
            corrvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                (rt->second).getSymmetricError());}

            k=0;
            for (map<double,MCEstimate>::const_iterator rt=efftrecon.begin();rt!=efftrecon.end();rt++,k++){
                corrreconvals[k]=XYPoint(rt->first, (rt->second).getFullEstimate());
                corrreconvals_below[k]=XYPoint(rt->first, (rt->second).getFullEstimate()-(rt->second).getSymmetricError());
                corrreconvals_above[k]=XYPoint(rt->first, (rt->second).getFullEstimate()+(rt->second).getSymmetricError());}

            plotfile = plotfilestub+"_eff_"+make_string(opindex)+".agr";
            createCorrelatorPlotWithRecon(corrvals,corrreconvals,corrreconvals_below,corrreconvals_above,
                                            RealPart,corrname,plotfile,symboltype,color,type);
            check_xml.put_child("EffPlotFile",plotfile);

            //plot energy with extra state
            if(fit_success){
                corrreconvals.resize(efftrecon2.size());
                corrreconvals_below.resize(efftrecon2.size());
                corrreconvals_above.resize(efftrecon2.size());

                k=0;
                for (map<double,MCEstimate>::const_iterator rt=efftrecon2.begin();rt!=efftrecon2.end();rt++,k++){
                    corrreconvals[k]=XYPoint(rt->first, (rt->second).getFullEstimate());
                    corrreconvals_below[k]=XYPoint(rt->first, (rt->second).getFullEstimate()-(rt->second).getSymmetricError());
                    corrreconvals_above[k]=XYPoint(rt->first, (rt->second).getFullEstimate()+(rt->second).getSymmetricError());}

                plotfile = plotfilestub+"_eff2_"+make_string(opindex)+".agr";
                createCorrelatorPlotWithRecon(corrvals,corrreconvals,corrreconvals_below,corrreconvals_above,
                                                RealPart,corrname,plotfile,symboltype,color,type);
                check_xml.put_child("EffPlotFile2",plotfile);
            }
            xml_out.put_child(check_xml);
        }
        XMLHandler final_fits("FinalFits");
        sort(newlevels.begin(), newlevels.end(), sortkey);
        for (uint i=0;i<newlevels.size();i++){
            XMLHandler fit("Fit");
            fit.put_child("Level",to_string(newlevels[i].level_insert));
            fit.put_child("Corr",newlevels[i].corr_name);
            fit.put_child("Amplitude",to_string(newlevels[i].amp.getFullEstimate()));
            fit.put_child("AmplitudeErr",to_string(newlevels[i].amp.getSymmetricError()));
            fit.put_child("ConstrainStrength",to_string(newlevels[i].amp.getSymmetricError()/newlevels[i].amp.getFullEstimate()));
            fit.put_child("Energy",to_string(newlevels[i].energy.getFullEstimate()));
            fit.put_child("EnergyErr",to_string(newlevels[i].energy.getSymmetricError()));
            fit.put_child("BestFitChisqr",to_string(newlevels[i].bestfit_chisqr));
            fit.put_child("LevelInsertChisqr",to_string(newlevels[i].new_level_chisqr));
            fit.put_child("BestFitChisqrDOF",to_string(newlevels[i].bestfit_chisqrdof));
            fit.put_child("Before",make_string(newlevels[i].before));
            fit.put_child("Priored",make_string(newlevels[i].priored));
            final_fits.put_child(fit);//add chisqr
        }
        xml_out.put_child(final_fits);

        // m_obs->begin();
        ofstream f; f.open(gevp_data_file);
        f<<"level,energy,energyerr"<<endl;
        for (uint level=0;level<nlevels;level++){
            MCEstimate energy=m_obs->getEstimate(pivoter.getEnergyKey(level));
            f<<level<<","<<energy.getFullEstimate()<<","<<energy.getSymmetricError()<<endl;
        }
        f.close();f.open(gevp_inserts_datafile);
        f<<"corr,energy,energyerr,insert,prior"<<endl;
        for (uint i=0;i<newlevels.size();i++){
            f<<newlevels[i].corr_name<<","<<newlevels[i].energy.getFullEstimate()<<","<<newlevels[i].energy.getSymmetricError()<<","<<newlevels[i].new_level_insert_value<<","<<newlevels[i].priored<<endl;
        }
        f.close();

    }catch(const exception& msg){
            xml_out.put_child("Error","doGEVPCheck: "+string(msg.what()));}
}