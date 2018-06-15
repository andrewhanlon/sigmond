#include <iostream>
#include <set>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
using namespace std;



typedef unsigned int   uint;

   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) 
   //
   //   Strategy:  fit a straight line to  ln(+/-corr(t)) vs t with least squares
   //   Discards correlator values that are opposite sign from small t corr.
   //   We assume that tvals[k] < tvals[k+1].

void get_exp_guess(const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                   double& energy0, double& amp0)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enought time vals in get_exp_guess"));
 for (uint k=0;k<(tvals.size()-1);k++)
    if (tvals[k]>=tvals[k+1]) throw(std::invalid_argument("increasing tvalues needed in get_exp_guess"));
 double sgn=(corrvals[0]>0.0)?1.0:-1.0;
 vector<double> x,y;
 for (uint k=0;k<corrvals.size();k++){
    if (sgn*corrvals[k]>0.0){
       x.push_back(tvals[k]);
       y.push_back(log(sgn*corrvals[k]));}
    else
       break;}
 uint npts=y.size();
 if (npts<2) 
    throw(std::invalid_argument("number of valid times too small for get_exp_guess"));
 double r0=npts;
 double rx=0.0,ry=0.0,rxx=0.0,rxy=0.0;    
 for (uint k=0;k<npts;k++){
    rx+=x[k]; ry+=y[k]; rxx+=x[k]*x[k]; rxy+=x[k]*y[k];}
 double den=r0*rxx-rx*rx;
 energy0=(rx*ry-r0*rxy)/den;
 amp0=sgn*exp((ry*rxx-rx*rxy)/den);
}


   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, gapsq, and gapamp, approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) * (1 + gapamp*exp(-gapsqrt*gapsqrt*t))
   //
   //   Key assumption: the second exponential must be negligible after "tasymfrac"
   //   of the way from the minimum to the maximum time separation.  Also,
   //   we assume that tvals[k]<tvals[k+1].
   //
   //   Strategy: - first, fit a straight line to  ln(+/-corr(t)) vs t with least squares
   //                from "tasymfrac" the way from tmin to tmax to obtain amp0, energy0;
   //             - form f = ln(+/-corr(t)) -ln(+/-amp0) + energy0 * t
   //             - fit straight line to ln(+/-f) with least squares to get gapamp, gapsqrt
   //             - if last part fails, just try gapamp=0.5, gapsqrt=0.5


void get_two_exp_guess(const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double tasymfrac=0.33)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_two_exp_guess"));
 if (tvals.size()<4)
    throw(std::invalid_argument("Need at least 4 points for get_two_exp_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<4)
    throw(std::invalid_argument("Insufficient time range for get_two_exp_guess"));
 double tasym=tmin+(tmax-tmin)*tasymfrac;
 vector<uint> tv; vector<double> cv;
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])>tasym){
       tv.push_back(tvals[k]);
       cv.push_back(corrvals[k]);}}
 get_exp_guess(tv,cv,energy0,amp0);
 tv.clear(); cv.clear();
 double sgn=(amp0>0.0)?1.0:-1.0;
 double A=log(sgn*amp0);
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])<tasym){
       double f=log(sgn*corrvals[k])-A+energy0*double(tvals[k]);
       tv.push_back(tvals[k]);
       cv.push_back(f);}}
  try{
     get_exp_guess(tv,cv,gapsqrt,gapamp);
     gapsqrt=sqrt(gapsqrt);}
  catch(const std::exception& xp){
     gapsqrt=0.5; gapamp=0.5;}
}

   //  Evaluates C(t)-C(t+step)
   //   tvals[k+1] - tvals[k] = step is required.

void take_diff(const std::vector<uint>& tvals, const std::vector<double>& corrvals,
               std::vector<uint>& tdf, std::vector<double>& corrdiff, uint& step)
{
 uint npts=tvals.size();
 if (corrvals.size()!=npts) 
    throw(std::invalid_argument("vector size mismatch in take_diff"));
 if (npts<3)
    throw(std::invalid_argument("Need at least 3 points for take_diff"));
 corrdiff.clear();
 tdf.clear();
 step=tvals[1]-tvals[0];
 for (uint k=0;k<(npts-1);k++){
    if ((tvals[k+1]-tvals[k])==step){
       tdf.push_back(tvals[k]); 
       corrdiff.push_back(corrvals[k]-corrvals[k+1]);}}
}



   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, c0  approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) + c0
   //
   //   Strategy: - compute D(t) = C(t)-C(t+step)
   //             - fit a straight line to  ln(+/-D(t)) vs t with least squares
   //             - extract c0 from first 1/4 of points
   //   Discards correlator values that are opposite sign from small t corr.
   //    tvals[k+1] - tvals[k] = step is required.

void get_exp_plus_const_guess(const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                              double& energy0, double& amp0, double& c0)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_exp_plus_const_guess"));
 vector<double> corrdiff; vector<uint> tdf;
 uint step;
 take_diff(tvals,corrvals,tdf,corrdiff,step);
 get_exp_guess(tdf,corrdiff,energy0,amp0);
 amp0/=(1.0-exp(-double(step)*energy0));
 c0=0.0;
 uint nc0=corrvals.size()/4;
 if (nc0<2) nc0=2;
 for (uint k=0;k<nc0;k++)
   c0+=corrvals[k]-amp0*exp(-energy0*tvals[k]);
 c0/=double(nc0);
}

   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, gapsq, gapamp, and c0, approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) * (1 + gapamp*exp(-gapsqrt*gapsqrt*t)) + c0
   //
   //   Key assumption: the second exponential must be negligible after "tasymfrac"
   //   of the way from the minimum to the maximum time separation.  Also,
   //    tvals[k+1] - tvals[k] = step is required.
   //
   //   Strategy: - for t=tasym ... tmax  compute D(t) = C(t)-C(t+step)
   //             - fit a straight line to  ln(+/-D(t)) vs t with least squares
   //             - extract c0 from first 1/4 of these points
   //             - form f = ln(+/-(corr(t)-c0)) -ln(+/-amp0) + energy0 * t
   //             - fit straight line to ln(+/-f) with least squares to get gapamp, gapsqrt
   //             - if last part fails, just try gapamp=0.5, gapsqrt=0.5

void get_two_exp_plus_const_guess(const std::vector<uint>& tvals, 
                       const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double& c0, double tasymfrac=0.33)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_two_exp_plus_const_guess"));
 if (tvals.size()<5)
    throw(std::invalid_argument("Need at least 5 points for get_two_exp_plus_const_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<5)
    throw(std::invalid_argument("Insufficient time range for get_two_exp_plus_const_guess"));
 double tasym=tmin+(tmax-tmin)*tasymfrac;
 vector<uint> tv; vector<double> cv;
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])>tasym){
       tv.push_back(tvals[k]);
       cv.push_back(corrvals[k]);}}
 get_exp_plus_const_guess(tv,cv,energy0,amp0,c0);
 tv.clear(); cv.clear();
 double sgn=(amp0>0.0)?1.0:-1.0;
 double A=log(sgn*amp0);
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])<tasym){
       double f=log(sgn*(corrvals[k]-c0))-A+energy0*double(tvals[k]);
       tv.push_back(tvals[k]);
       cv.push_back(f);}}
  try{
     get_exp_guess(tv,cv,gapsqrt,gapamp);
     if (gapsqrt>0.0) gapsqrt=sqrt(gapsqrt);
     else throw(std::runtime_error("invalid exponential decay"));}
  catch(const std::exception& xp){
     gapsqrt=0.3; gapamp=0.5;}
}



int main()
{

 vector<uint> tvals;
 vector<double> corrvals;
 double energy0,amp0,gapamp,gapsqrt;

 for (uint t=3;t<=26;t++){
    tvals.push_back(t);
    corrvals.push_back(-4.53*exp(-0.25*t));}

 get_exp_guess(tvals,corrvals,energy0,amp0);
 cout << "energy0 = "<<energy0<<endl;
 cout << "amp0 = "<<amp0<<endl<<endl;

 tvals.clear(); corrvals.clear();
 for (uint t=12;t<=26;t++){
    tvals.push_back(t);
    corrvals.push_back(-4.53*exp(-0.25*t)*(1.0+0.75*exp(-0.8*t)));}

 get_exp_guess(tvals,corrvals,energy0,amp0);
 cout << "energy0 = "<<energy0<<endl;
 cout << "amp0 = "<<amp0<<endl<<endl;

 tvals.clear(); corrvals.clear();
 for (uint t=3;t<=26;t++){
    tvals.push_back(t);
    corrvals.push_back( -4.53*exp(-0.25*t)*(1.0+0.75*exp(-0.8*t)) );}

 get_two_exp_guess(tvals,corrvals,energy0,amp0,gapsqrt,gapamp);
 cout << "energy0 = "<<energy0<<endl;
 cout << "amp0 = "<<amp0<<endl;
 cout << "gapsqrt = "<<gapsqrt<<endl;
 cout << "gapamp = "<<gapamp<<endl<<endl;

 tvals.clear(); corrvals.clear();
 for (uint t=3;t<=26;t+=2){
    tvals.push_back(t);
    corrvals.push_back( -4.53*exp(-0.25*t) - 3.4 );}

 double c0;
 get_exp_plus_const_guess(tvals,corrvals,energy0,amp0,c0);
 cout << "energy0 = "<<energy0<<endl;
 cout << "amp0 = "<<amp0<<endl;
 cout << "c0 = "<<c0<<endl<<endl;

 tvals.clear(); corrvals.clear();
 for (uint t=3;t<=26;t++){
    tvals.push_back(t);
    corrvals.push_back( -4.53*exp(-0.25*t)*(1.0-0.75*exp(-0.5*t)) - 3.4 );}

 get_two_exp_plus_const_guess(tvals,corrvals,energy0,amp0,gapsqrt,gapamp,c0);
 cout << "energy0 = "<<energy0<<endl;
 cout << "amp0 = "<<amp0<<endl;
 cout << "gapsqrt = "<<gapsqrt<<endl;
 cout << "gapamp = "<<gapamp<<endl;
 cout << "c0 = "<<c0<<endl<<endl;

 return 0;
}
