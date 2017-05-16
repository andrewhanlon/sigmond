#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <stdexcept>
using namespace std;

typedef unsigned int uint;


void eval_guess(int tval, double corrt, int tnext, double corrtnext, 
                double& A, double& m)
{
 double s=corrt/corrtnext;
 if ((s<=0.0)||(tval==tnext))
    throw(std::invalid_argument("SingleExponential -- could not compute a guess for exponential"));
 m=log(s)/(double(tnext)-double(tval));       // guess for m
 A=exp(m*double(tval))*corrt;                 // guess for A
}


void eval_guess(int tval, double corrt, int tp1, double corrtp1, int tp2, double corrtp2,
                double& A, double& m, double& c0)
{
 if (((tp2-tp1)!=(tp1-tval))||(tval==tp1))
    throw(std::invalid_argument("SingleExponentialPlusConst -- could not compute a guess for exponential"));
 double cor0=corrtp1-corrt;
 double cor1=corrtp2-corrtp1;
 double s=cor0/cor1;
 double d=double(tp1)-double(tval);
 if (s<=0.0)
    throw(std::invalid_argument("SingleExponentialPlusConst -- could not compute a guess for exponential"));
 m=log(s)/d;
 double rr=exp(-m*double(tval));
 A=cor0/(rr*(exp(-m*d)-1.0));
 c0=corrt-A*rr;  
}


     //  initial guess for a two-exponential fit 
     //     A * exp( -m*t ) * (1 + B * exp(-DD^2*t) )
     //
     //  choose tfar in large time region where second
     //  exponential is negligible, then choose tnear
     //  in small time region where second exponential
     //  can be exposed;  ffar = corr(tfar), ffarnext=corr(tfarnext)
     //  and fnear=corr(tnear), fnearnext=corr(tnearnext)

void eval_guess(
              int tfar, double ffar, int tfarnext, double ffarnext, 
              int tnear, double fnear, int tnearnext, double fnearnext,
              double& A, double& m, double& B, double& DD)
{
 double s=ffar/ffarnext;
 if ((s<=0.0)||(tfarnext==tfar))
    throw(std::invalid_argument("Two exponential -- could not compute an initial guess"));
 m=log(s)/(double(tfarnext)-double(tfar)); // guess for m
 A=exp(m*double(tfar))*ffar;               // guess for A
 double r1=exp(m*double(tnear))*fnear/A-1.0;
 double r2=exp(m*double(tnearnext))*fnearnext/A-1.0;
 s=r1/r2;
 if ((s<=1.0)||(tnearnext==tnear))
    throw(std::invalid_argument("Two exponential -- could not compute an initial guess"));
 double DDsq=log(s)/(double(tnearnext)-double(tnear)); // guess for DD
 DD=sqrt(DDsq);
 B=exp(DDsq*double(tnear))*r1;       // guess for B
}





     //  initial guess for a two-exponential fit 
     //     A * exp( -m*t ) * (1 + B * exp(-DD^2*t) ) + c0
     //
     //  choose tfar in large time region where second
     //  exponential is negligible, then choose tnear
     //  in small time region where second exponential
     //  can be exposed;  
     //   ffar = corr(tfar), ffarp1=corr(tfarp1), ffarp2=corr(tfarp2)
     //  and fnear=corr(tnear), fnearnext=corr(tnearnext)

void eval_guess(
              int tfar, double ffar, int tfarp1, double ffarp1, int tfarp2, double ffarp2,
              int tnear, double fnear, int tnearnext, double fnearnext,
              double& A, double& m, double& B, double& DD, double& c0)
{
 if (((tfarp2-tfarp1)!=(tfarp1-tfar))||(tfar==tfarp1))
    throw(std::invalid_argument("Two exponential + const -- could not compute an initial guess"));
 double cor0=ffarp1-ffar;
 double cor1=ffarp2-ffarp1;
 double s=cor0/cor1;
 double d=double(tfarp1)-double(tfar);
 if (s<=0.0)
    throw(std::invalid_argument("Two exponential + const -- could not compute an initial guess"));
 m=log(s)/d;                         // guess for m
 double rr=exp(-m*double(tfar));
 A=cor0/(rr*(exp(-m*d)-1.0));        // guess for A
 c0=ffar-A*rr;                       // guess for c0
 double r1=exp(m*double(tnear))*(fnear-c0)/A-1.0;   
 double r2=exp(m*double(tnearnext))*(fnearnext-c0)/A-1.0; 
 s=r1/r2;
 if ((s<=1.0)||(tnearnext==tnear))
    throw(std::invalid_argument("Two exponential + const -- could not compute an initial guess"));
 double DDsq=log(s)/(double(tnearnext)-double(tnear)); // guess for DD
 DD=sqrt(DDsq);
 B=exp(DDsq*double(tnear))*r1;       // guess for B
}




int main(){

cout.precision(15);

double c0=-4.21;
double A=34.7;
double m=0.234;
double B=0.562;
double DD=0.771;

double eA,em,eB,eDD,ec0;

cout << "Testing single exponential"<<endl<<endl;
vector<double> corr(30);
vector<int> tvals(30);
for (int t=0;t<30;t++){
   tvals[t]=t;
   corr[t]=A*exp(-m*double(t));}

for (uint k=0;k<15;k++)
for (uint j=1;j<8;j++){
   eA=0.0; em=0.0;
   try{
      eval_guess(tvals[k],corr[k],tvals[k+j],corr[k+j],eA,em);
      cout << "eA = "<<eA<<"  em = "<<em<<endl;}
   catch(const std::exception& xp){ cout << xp.what()<<endl;}}


cout << "Testing single exponential plus constant"<<endl<<endl;
for (int t=0;t<30;t++){
   tvals[t]=t;
   corr[t]=c0+A*exp(-m*double(t));}

for (uint k=0;k<15;k++)
for (uint j=1;j<8;j++)
for (uint l=1;l<4;l++){
   eA=0.0; em=0.0;
   try{
      eval_guess(tvals[k],corr[k],tvals[k+j],corr[k+j],tvals[k+j+l],corr[k+j+l],eA,em,ec0);
      cout << "eA = "<<eA<<"  em = "<<em<<" ec0 = "<<ec0<<endl;}
   catch(const std::exception& xp){ cout << xp.what()<<endl;}}

cout << "Testing two exponentials"<<endl<<endl;
for (int t=0;t<30;t++)
   corr[t]=A*exp(-m*double(t))*(1.0+B*exp(-DD*DD*double(t)));

eA=em=eDD=eB=0.0;
try{
   eval_guess(25,corr[25],26,corr[26],2,corr[2],3,corr[3],eA,em,eB,eDD);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}
eA=em=eDD=eB=0.0;
try{
   eval_guess(25,corr[25],27,corr[27],2,corr[2],4,corr[4],eA,em,eB,eDD);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}
eA=em=eDD=eB=0.0;
try{
   eval_guess(25,corr[25],28,corr[28],2,corr[2],5,corr[5],eA,em,eB,eDD);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}
eA=em=eDD=eB=0.0;
try{
   eval_guess(26,corr[26],26,corr[26],0,corr[0],3,corr[3],eA,em,eB,eDD);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}
eA=em=eDD=eB=0.0;
try{
   eval_guess(26,corr[26],27,corr[27],1,corr[1],3,corr[3],eA,em,eB,eDD);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}
eA=em=eDD=eB=0.0;
try{
   eval_guess(26,corr[26],28,corr[28],2,corr[2],4,corr[4],eA,em,eB,eDD);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}


cout << "Testing two exponentials + constant"<<endl<<endl;
for (int t=0;t<30;t++)
   corr[t]=c0+A*exp(-m*double(t))*(1.0+B*exp(-DD*DD*double(t)));

eA=em=eDD=eB=ec0=0.0;
try{
   eval_guess(25,corr[25],26,corr[26],27,corr[27],2,corr[2],3,corr[3],eA,em,eB,eDD,ec0);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<" ec0 = "<<ec0<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}

eA=em=eDD=eB=ec0=0.0;
try{
   eval_guess(25,corr[25],27,corr[27],29,corr[29],2,corr[2],4,corr[4],eA,em,eB,eDD,ec0);
   cout << "eA = "<<eA<< " em = "<<em<<" eB = "<<eB<<" eDD = "<<eDD<<" ec0 = "<<ec0<<endl;}
catch(const std::exception& xp){ cout << xp.what()<<endl;}

return 0;
}
