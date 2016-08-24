#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;



void eval_func0(
              double A, double m, double B, double DD,
              double tf, double& funcval)
{
 funcval=A*(exp(-m*tf)/(1.0-B*exp(-DD*DD*tf)));
}


void eval_grad0(
              double A, double m, double B, double DD,
              double tf, double& dAval, double& dmval,
              double& dBval, double& dDDval)
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 double d1=1.0-B*r2;
 dAval=r1/d1;
 dmval=-tf*A*dAval;
 dBval=A*r1*r2/(d1*d1);
 dDDval=-2.0*tf*B*DD*dBval;
}


void eval_func(
                double A, double m, double B, double DD,
                double tf, int Nt, double& funcval)
{
 double tb=double(Nt)-tf;
 funcval=A*(exp(-m*tf)/(1.0-B*exp(-DD*DD*tf))+exp(-m*tb)/(1.0-B*exp(-DD*DD*tb)));
}


void eval_grad(
                double A, double m, double B, double DD,
                double tf, int Nt, double& dAval, double& dmval,
                double& dBval, double& dDDval)
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 double d1=1.0-B*r2;
 dAval=r1/d1;
 dmval=-tf*A*dAval;
 dBval=A*r1*r2/(d1*d1);
 dDDval=-2.0*tf*B*DD*dBval;
 double tb=double(Nt)-tf;
 r1=exp(-m*tb); r2=exp(-gap*tb);
 d1=1.0-B*r2;
 double s=r1/d1;
 dAval+=s;
 dmval-=tb*A*s;
 s=A*r1*r2/(d1*d1);
 dBval+=s;
 dDDval-=2.0*tb*B*DD*s;
}


int main(){

cout.precision(15);
cout << "#tester routine for geomseriesexponential:"<<endl;
cout << "Digits:=18:"<<endl;
double A=1.34;
double m=0.14;
double B=0.4;
double DD=0.27;
double tf=2.0;
int Nt=9;

cout << "func:=A*exp(-m*t)/(1-B*exp(-DD^2*t))+A*exp(-m*(Nt-t))/(1-B*exp(-DD^2*(Nt-t))):"<<endl;
cout << "dfunc_dA:=diff(func,A):"<<endl;
cout << "dfunc_dB:=diff(func,B):"<<endl;
cout << "dfunc_dm:=diff(func,m):"<<endl;
cout << "dfunc_dDD:=diff(func,DD):"<<endl;

cout << "A:="<<A<<":"<<endl;
cout << "B:="<<B<<":"<<endl;
cout << "m:="<<m<<":"<<endl;
cout << "DD:="<<DD<<":"<<endl;
cout << "t:="<<tf<<":"<<endl;
cout << "Nt:="<<Nt<<":"<<endl;

cout << "fval:=evalf(func):"<<endl;
double funcval;
eval_func(A,m,B,DD,tf,Nt,funcval);
cout << "fcheck:="<<funcval<<":"<<endl;
cout << "print(\"diff = \",fcheck-fval):"<<endl;

double dAval,dmval,dBval,dDDval;
eval_grad(A,m,B,DD,tf,Nt,dAval,dmval,dBval,dDDval);
cout << "res:=evalf(dfunc_dA):"<<endl;
cout << "chk:="<<dAval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;
cout << "res:=evalf(dfunc_dB):"<<endl;
cout << "chk:="<<dBval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;
cout << "res:=evalf(dfunc_dm):"<<endl;
cout << "chk:="<<dmval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;
cout << "res:=evalf(dfunc_dDD):"<<endl;
cout << "chk:="<<dDDval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;

cout << "A:='A':"<<endl;
cout << "B:='B':"<<endl;
cout << "m:='m':"<<endl;
cout << "DD:='DD':"<<endl;
cout << "t:='t':"<<endl;

cout << "func:=A*exp(-m*t)/(1-B*exp(-DD^2*t)):"<<endl;
cout << "dfunc_dA:=diff(func,A):"<<endl;
cout << "dfunc_dB:=diff(func,B):"<<endl;
cout << "dfunc_dm:=diff(func,m):"<<endl;
cout << "dfunc_dDD:=diff(func,DD):"<<endl;

cout << "A:="<<A<<":"<<endl;
cout << "B:="<<B<<":"<<endl;
cout << "m:="<<m<<":"<<endl;
cout << "DD:="<<DD<<":"<<endl;
cout << "t:="<<tf<<":"<<endl;

cout << "fval:=evalf(func):"<<endl;
eval_func0(A,m,B,DD,tf,funcval);
cout << "fcheck:="<<funcval<<":"<<endl;
cout << "print(\"diff = \",fcheck-fval):"<<endl;

eval_grad0(A,m,B,DD,tf,dAval,dmval,dBval,dDDval);
cout << "res:=evalf(dfunc_dA):"<<endl;
cout << "chk:="<<dAval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;
cout << "res:=evalf(dfunc_dB):"<<endl;
cout << "chk:="<<dBval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;
cout << "res:=evalf(dfunc_dm):"<<endl;
cout << "chk:="<<dmval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;
cout << "res:=evalf(dfunc_dDD):"<<endl;
cout << "chk:="<<dDDval<<":"<<endl;
cout << "print(\"diff = \",res-chk):"<<endl;


return 0;
}
