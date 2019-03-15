#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;

string double_to_string(double dval, unsigned int precision)
{
 const int fmax=32;  // first make the format string
 const int bmax=64;
 char format[fmax];
 char charbuffer[bmax];
 int prec=(precision>0)?precision:0;

 format[0]='%'; format[1]='.';
 int csize=snprintf(&format[2],fmax-6,"%i",prec);
 if (csize<0)
    throw(string(" error in converting double to string"));
 csize+=2;
 format[csize++]='g'; format[csize++]='\0';
 csize=snprintf(charbuffer,bmax,format,dval); // now make string from double
 if (csize<0){
    throw(string(" error in converting integer to charString"));}
 return string(charbuffer);
}

/*
double exp10(int power)
{
 if (power==0) return 1.0;
 double cf=1.0;
 int i;
 if (power<0)
    for (i=1;i<=-power;i++) cf*=0.1;
 else
    for (i=1;i<=power;i++) cf*=10.0;
 return cf;
}

int posint_exp10(int power)
{
 if (power<0){
    throw(string("posint_exp10 requires positive power"));
    return 0;}
 if (power==0) return 1;
 int cf=1;
 int i;
 for (i=1;i<=power;i++) cf*=10;
 return cf;
}
*/

string measure_to_string(unsigned int nerr_digits, double val, double err)
{
 ostringstream ch;
 long int exponent,ndec,errint,valint;
 double errv,valv;
 string res;

 if (nerr_digits<1) ndec=1;
 else ndec=nerr_digits;
 if (err<=0.0){
    return string("BAD_VALUE");}

 exponent=(int) floor(log10(err));
 exponent-=ndec-1;
     //   rescale so that error is now a number between 10^(nerr_digits-1)
     //   and 10^(nerr_digits), then round, and readjust in case rounding
     //   up increased number of digits
 errv=err*exp10(-exponent);
 errint=(long int) floor(errv+0.5);          
 if (errint>=(long int) exp10(ndec)){
    errint/=10;
    exponent++;}

 valv=val*exp10(-exponent);
 valint=(long int) floor(std::abs(valv)+0.5);

 if (exponent>0){
    valint*=exp10(exponent);
    errint*=exp10(exponent);
    }
 if (exponent>=0){
      // make into a string
    ch << valint << "(" << errint << ")";
    res=ch.str();
    }
 else{
    ch << valint << "(" << errint << ")";
    res=ch.str();
    int i,pos=0;
    while (res[pos]!='(') pos++;  // point to "("
    for (i=pos;i<=-exponent;i++) res="0"+res; // pad with 0's on left
    pos=0;
    while (res[pos]!='(') pos++;  // point to "("
    pos+=exponent;
    res=res.insert(pos,1,'.');
    if (ndec+exponent>0){
       pos=0;
       while (res[pos]!=')') pos++;  // point to ")"
       pos+=exponent;
       res=res.insert(pos,1,'.');
       }
    }

    // include negative sign if required
 if ((val<=0.0)&&(valint>0)) res="-"+res;
 return res;
}




int main(){

double val=342.34567;
int nprec=2;

ostringstream out;
out.precision(nprec);
out.setf(std::ios::fixed);
out << val<<endl;

val=0.00000453;
out.setf(std::ios::scientific);
out << val<<endl;


val=3248765436.3235;
out << val<<endl;


cout << out.str()<<endl;

cout << double_to_string(val,12)<<endl;

cout << "Errors now"<<endl;

cout << measure_to_string(1,0.34567,0.003456)<<endl;
cout << measure_to_string(2,0.34567,0.003456)<<endl;
cout << measure_to_string(3,0.34567,0.003456)<<endl;

cout << measure_to_string(1,7845.34567,0.003456)<<endl;
cout << measure_to_string(2,1234.34567,0.003456)<<endl;
cout << measure_to_string(3,7832.34567,0.003456)<<endl;

cout << measure_to_string(1,7845.34567,0.3456)<<endl;
cout << measure_to_string(2,7845.34567,0.3456)<<endl;
cout << measure_to_string(3,7845.34567,0.3456)<<endl;
cout << measure_to_string(1,1234.34567,2.3456)<<endl;
cout << measure_to_string(2,1234.34567,2.3456)<<endl;
cout << measure_to_string(3,1234.34567,2.3456)<<endl;
cout << measure_to_string(1,7832.34567,12.3456)<<endl;
cout << measure_to_string(2,7832.34567,12.3456)<<endl;
cout << measure_to_string(3,7832.34567,12.3456)<<endl;

cout << measure_to_string(1,783287934354.34567,12435.3456)<<endl;
cout << measure_to_string(2,783287934354.34567,12435.3456)<<endl;
cout << measure_to_string(3,783287934354.34567,12435.3456)<<endl;

return 0;
}
