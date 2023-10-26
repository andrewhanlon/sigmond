#include <string>
#include <iostream>
using namespace std;

extern "C" {

  // C only has pass by value

void sigmondReadSamplings(const char* fileName, const char* ensembleId,
                          double **results, int *nresults)
{
 if ((*results)!=0) delete *results;
 std::string fname(fileName);
 std::string ensId(ensembleId);
 *nresults=ensId.length();
 *results=new double[*nresults];
 for (int k=0;k<*nresults;++k) 
    *(*results+k)=double(k);
}

}  // extern C end


/*
void sigmondReadBins(const char* fileName, const char* ensembleId,
*/


int main(){

char fname[]="theFile";
char ens[]="theEnsemble";

double *results=0;
int nresults;

sigmondReadSamplings(fname,ens,&results,&nresults);

cout << "nresults = "<<nresults<<endl;
for (int k=0;k<nresults;++k)
   cout << "results["<<k<<"] = "<<*(results+k)<<endl;

delete results;

return 0;
}
