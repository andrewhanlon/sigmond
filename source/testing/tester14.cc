#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
using namespace std;



int main(){

unsigned int n=12;
std::vector<double> buffer(n,std::numeric_limits<double>::quiet_NaN());

buffer[0]=0.0;
buffer[1]=1.0;

for (unsigned int k=0;k<n;k++){
   if (!isnan(buffer[k]))
      cout << "buffer["<<k<<"] = "<<buffer[k]<<endl;}


return 0;
}
