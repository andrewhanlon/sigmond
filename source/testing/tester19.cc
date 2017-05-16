#include <iostream>
#include <set>
#include <vector>
using namespace std;

typedef unsigned int   uint;


std::vector<uint> form_tvalues(uint tmin, uint tmax, const std::vector<int>& texclude)
{
 set<uint> tvals;  // values will automatically be sorted in set
 for (uint tt=tmin;tt<=tmax;tt++){
    tvals.insert(tt);}
 for (uint k=0;k<texclude.size();k++){
    tvals.erase(texclude[k]);}
 vector<uint> result(tvals.begin(),tvals.end());
 return result;
}


int main()
{

 uint tmin=3; 
 uint tmax=12;
 vector<int> texclude;
 texclude.push_back(18);
 texclude.push_back(12);
 texclude.push_back(5);

 vector<uint> res(form_tvalues(tmin,tmax,texclude));
 for (uint k=0;k<res.size();k++)
    cout << res[k]<<endl;

 return 0;
}
