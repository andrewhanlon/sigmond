#include "bootstrapper.h"
#include <iostream>
#include "xml_handler.h"
using namespace std;



void testBootstrapper(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestBootstrapper")==0)
 return;

 unsigned int nobjects=854;
 unsigned int nsamples=2400;
 unsigned long rngseed=871234;
 unsigned int skip=18;
 bool precompute=true;
 Bootstrapper B(nobjects,nsamples,rngseed,skip,precompute);

 nobjects=1344;
 nsamples=18920;
 rngseed=62101;
 skip=33;
 precompute=false;
 B.reset(nobjects,nsamples,rngseed,skip,precompute);

 Vector<unsigned int> counters(nobjects,0);

 cout << "number of objects = "<<B.getNumberOfObjects()<<endl;
 cout << "number of samples = "<<B.getNumberOfResamplings()<<endl;
 cout << "seed = "<<B.getRNGSeed()<<endl;
 cout << "skip = "<<B.getSkipValue()<<endl;

 Matrix<unsigned int> M(nsamples,nobjects);
 
 for (unsigned int k=0;k<B.getNumberOfResamplings();++k){
    const Vector<unsigned int>& sample=B.getResampling(k);
    cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
    for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
       if ((sample[kk]<0)||(sample[kk]>=nobjects)) cout << "OH WRONG!"<<endl;
       counters[sample[kk]]++;
       M(k,kk)=sample[kk];}}

 {int k=24;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}
 {int k=26;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}
 {int k=29;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}
 {int k=11;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}
 {int k=1;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}
 {int k=0;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}
 {int k=782;
 const Vector<unsigned int>& sample=B.getResampling(k);
 cout << "B current sample index = "<<B.getCurrentResamplingCount()<<endl;
 for (unsigned int kk=0;kk<B.getNumberOfObjects();++kk){
     if (M(k,kk)!=sample[kk]) cout << "FAIL"<<endl;}}

 cout << endl<<endl;
 double temp=0.0;
 for (unsigned int k=0;k<counters.size();++k)
    temp+=counters[k];

 for (unsigned int k=0;k<counters.size();++k)
    cout << "counter["<<k<<"] = "<<(double(counters[k])/temp*nobjects)<<endl;

}

