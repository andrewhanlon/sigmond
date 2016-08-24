#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;


     //  initial guess for a two-exponential fit 
     //     A * exp( -m*t ) * (1 + B * exp(-DD^2*t) )
     //
     //  choose tfar in large time region where second
     //  exponential is negligible, then choose tnear
     //  in small time region where second exponential
     //  can be exposed;  ffar = corr(tfar), ffarnext=corr(tfar+1)
     //  and fnear=corr(tnear), fnearnext=corr(tnear+1)

void eval_guess(
              int tfar, double ffar, double ffarnext, 
              int tnear, double fnear, double fnearnext,
              double& A, double& m, double& B, double& DD)
{
 double s=ffar/ffarnext;
 if (s<=0.0)
    throw(string("Two exponential -- could not compute an initial guess"));
 m=log(s);                         // guess for m
 A=exp(m*double(tfar))*ffar;       // guess for A
 double r1=exp(m*double(tnear))*fnear/A-1.0;
 double r2=exp(m*double(tnear+1))*fnearnext/A-1.0;
 s=r1/r2;
 if (s<=1.0)
    throw(string("Two exponential -- could not compute an initial guess"));
 double DDsq=log(s);                         // guess for DD
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
     //   ffar = corr(tfar), ffarp1=corr(tfar+1), ffarp2=corr(tfar+2)
     //  and fnear=corr(tnear), fnearnext=corr(tnear+1)

void eval_guess(
              int tfar, double ffar, double ffarp1, double ffarp2,
              int tnear, double fnear, double fnearnext,
              double& A, double& m, double& B, double& DD, double& c0)
{
 double cor0=ffarp1-ffar;
 double cor1=ffarp2-ffarp1;
 double s=cor0/cor1;
 if (s<=0.0)
    throw(string("Two exponential + const -- could not compute an initial guess"));
 m=log(s);                         // guess for m
 double rr=exp(-m*double(tfar));
 A=cor0/(rr*(exp(-m)-1.0));        // guess for A
 c0=ffar-A*rr;                     // guess for c0
 double r1=exp(m*double(tnear))*(fnear-c0)/A-1.0;   
 double r2=exp(m*double(tnear+1))*(fnearnext-c0)/A-1.0; 
 s=r1/r2;
 if (s<=1.0)
    throw(string("Two exponential + const -- could not compute an initial guess"));
 double DDsq=log(s);                         // guess for DD
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

vector<double> corr(30);
for (int t=0;t<30;t++)
   corr[t]=c0+A*exp(-m*double(t))*(1.0+B*exp(-DD*DD*double(t)));

double eA,em,eB,eDD,ec0;
try{eval_guess(25,corr[25],corr[26],corr[27],2,corr[2],corr[3],eA,em,eB,eDD,ec0);}
catch(const string& errmsg){ cout << errmsg<<endl;}

cout << "eA = "<<eA<<endl;
cout << "em = "<<em<<endl;
cout << "eB = "<<eB<<endl;
cout << "eDD = "<<eDD<<endl;
cout << "ec0 = "<<ec0<<endl;


return 0;
}
