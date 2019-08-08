#include <iostream>
#include <set>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
using namespace std;



typedef unsigned int   uint;

   //  This transforms the CLS weights by rebinning and omitting
   //  certain weights.  A rebinned weight is simply the sum of the
   //  original weights in the bin.

void transform_weights(const vector<double>& worig,
                                        vector<double>& wnew, uint rebin,
                                        const std::set<uint>& omit)
{
 uint nomit=omit.size();
 uint norig=worig.size();
 for (set<uint>::const_reverse_iterator rm=omit.rbegin();rm!=omit.rend();++rm){
    if (*rm>norig) --nomit;
    else break;}
 if (rebin==1){
    if (nomit==0){         // no rebinning, no omissions
       wnew=worig;
       return;}
    else{                  // no rebinning, omissions
       uint nbins=norig-nomit;
       set<unsigned int>::const_iterator om=omit.begin();
       unsigned int count=0;
       wnew.resize(nbins);
       while ((om!=omit.end())&&(count==*om)){om++; ++count;}
       for (unsigned int k=0;k<nbins;k++){
          wnew[k]=worig[count];
          ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}
       return;}}
 if (nomit==0){         // rebinning but no omissions
    uint nbins=norig/rebin;
    unsigned int count=0;
    wnew.resize(nbins);
    for (unsigned int k=0;k<nbins;k++){
       double r=worig[count++];
       for (unsigned int j=1;j<rebin;j++){
          r+=worig[count++];}
       wnew[k]=r;}
    return;}
 uint nbins=(norig-nomit)/rebin;   // rebinning and omissions
 set<unsigned int>::const_iterator om=omit.begin();
 unsigned int count=0;
 wnew.resize(nbins);
 while ((om!=omit.end())&&(count==*om)){om++; ++count;}
 for (unsigned int k=0;k<nbins;k++){
    double r=worig[count];
    ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}
    for (unsigned int j=1;j<rebin;j++){
       r+=worig[count];
       ++count; while ((om!=omit.end())&&(count==*om)){om++; ++count;}}
     wnew[k]=r;}
}



int main()
{

 uint n=25;
 vector<double> worig(n);
 for (uint k=0;k<n;++k)
    worig[k]=double(k);
 vector<double> wnew;
 set<uint> omit;
 omit.insert(54);
 omit.insert(24);
 omit.insert(8);
 omit.insert(13);
 omit.insert(44);
 uint rebin=3;

 transform_weights(worig,wnew,rebin,omit);
 for (uint k=0;k<wnew.size();++k)
    cout << "wnew["<<k<<"] = "<<wnew[k]<<endl;

 return 0;
}
