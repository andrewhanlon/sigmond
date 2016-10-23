#include "mcobs_handler.h"
#include "test_obs_get_handler.h"
#include <cstdio>
#include <ctime>
#include <map>
#include <cstdlib>



 // **********************************

class BinMapper
{
    std::vector<std::list<int> > configs;

 public:

    BinMapper(int nconfigs, int rebin, const std::set<int>& omit);

    const std::list<int>& getBinConfigs(int bin) const {return configs.at(bin);}

};



bool getVEVBin(VEVCorrect& VC, const BinMapper& BM, 
               const MCObsInfo& mcobs, int bin, double& result);


bool getCorBin(CorrCorrect& CC, const BinMapper& BM, 
               const MCObsInfo& mcobs, int bin, double& result);


  // *************************************************************************
