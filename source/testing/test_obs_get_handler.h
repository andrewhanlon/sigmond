#ifndef TEST_OBS_GET_HANDLER_H
#define TEST_OBS_GET_HANDLER_H
#include "obs_get_handler.h"
#include "ensemble_info.h"
#include <map>
#include <set>
#include <string>


class VEVCorrect
{
    std::map<OperatorInfo,std::map<int,InScalar> > values;

 public:

    void addValue(const OperatorInfo& opinfo, int serind, const InScalar& invalue);
    
    bool getCorrect(const MCObsInfo& mcobs, int serind, double& result) const;

    bool getCorrectBins(const MCObsInfo& mcobs, const MCBinsInfo& bininfo, 
                        RVector& result) const;

};


class CorrCorrect
{
    std::map<CorrelatorAtTimeInfo,std::map<int,InScalar> > values;
    bool hermitian;

 public:

    CorrCorrect() : hermitian(false) {}

    void setHermitian() {hermitian=true;}

    void addValue(const CorrelatorAtTimeInfo& opinfo, int serind, const InScalar& invalue);
    
    bool getCorrect(const MCObsInfo& mcobs, int serind, double& result) const;

    bool getCorrectBins(const MCObsInfo& mcobs, const MCBinsInfo& bininfo, 
                        RVector& result) const;
};

bool compare_vectors(const RVector& v1, const RVector& v2, double factor);

bool compare_floats(const double& v1, const double& v2, double factor);


void make_fake_correlator_file(const MCEnsembleInfo& ens, const std::string& filename,
                               const CorrelatorInfo& corr,
                               int tmin, int tmax, double Aval, double Bval, 
                               double Cval, const std::set<int>& omit,
                               CorrCorrect& CC);

void make_fake_vev_file(const MCEnsembleInfo& ens, const std::string& filename,
                        const OperatorInfo& opinfo, double Aval, double Cval, 
                        const std::set<int>& omit, VEVCorrect& VC);



class BinsCorrect
{
    std::map<MCObsInfo,RVector > values;

 public:

    BinsCorrect() {}

    void addValue(const MCObsInfo& opinfo, const RVector& inbins);
    
    bool getCorrect(const MCObsInfo& mcobs, RVector& result) const;

};

void make_fake_bin_file(const MCObsInfo& obskey, const std::string& filename,
                        const MCBinsInfo& bininfo, double Aval, double Cval, 
                        BinsCorrect& BC);


double get_random_double();

class SamplingsCorrect
{
    std::map<MCObsInfo,RVector > values;

 public:

    SamplingsCorrect() {}

    void addValue(const MCObsInfo& opinfo, const RVector& insamp);
    
    bool getCorrect(const MCObsInfo& mcobs, RVector& result) const;

};

void make_fake_samplings_file(const MCObsInfo& obskey, const std::string& filename,
                              const MCBinsInfo& bininfo, const MCSamplingInfo& sampinfo, 
                              double Aval, double Cval, SamplingsCorrect& SC);

#endif
