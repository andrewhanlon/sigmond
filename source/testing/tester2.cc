#include <vector>
#include <iostream>
#include <list>
#include <complex>

using namespace std;


    //  returns  x^k  where k = integer

double ipow(double x, int k)
{
 if (k==0) return 1.0;
 else if (k==1) return x;
 else if (k==-1) return 1.0/x;
 double xx=(k>0)?x:1.0/x;
 unsigned int kk=(k>0)?k:-k;
 double val=1.0;
 double w=xx;
 while (kk){
    if (kk&1u) val*=w;
    w*=w; kk>>=1;}
 return val;
}


   //   Base class for evaluating a function and its derivative of
   //   a single variable.  The crucial operator is
   //       F(var,funcval,derivval);

class FuncAndDerivSingleVar
{
 public:
    virtual ~FuncAndDerivSingleVar(){};
    virtual void operator()(double var, double& funcval, double& derivval) = 0;
};


  //   Using a combination of Newton-Raphson and bisection, finds the
  //   root of a function f(x) bracketed between x1 and x2.  The root, returned
  //   as the function value, will be refined until its accuracy is known
  //   within +/- xacc.  "funcd" is a user-supplied routine that
  //   returns both the function value and its derivative.  The object of type T
  //   must have a member function  .eval_func_and_deriv(double,double&,double&)

double rtsafe(FuncAndDerivSingleVar& funcd, double x1, double x2, double xacc, unsigned int maxit)
{
 int j;
 double df,dx,dxold,f,fh,fl;
 double temp,xh,xl,rts;

 funcd(x1,fl,df);
 funcd(x2,fh,df);
 if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    throw(std::string("Root must be bracketed in rtsafe"));
 if (fl == 0.0) return x1;
 if (fh == 0.0) return x2;
 if (fl < 0.0) {
    xl=x1; xh=x2;} 
 else {
    xh=x1; xl=x2;}
 rts=0.5*(x1+x2);
 dxold=std::abs(x2-x1);
 dx=dxold;
 funcd(rts,f,df);
 for (j=1;j<=maxit;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
       || (std::abs(2.0*f) > std::abs(dxold*df))) {
       dxold=dx;
       dx=0.5*(xh-xl);
       rts=xl+dx;
       if (xl == rts) return rts;
    } else {
       dxold=dx;
       dx=f/df;
       temp=rts;
       rts -= dx;
       if (temp == rts) return rts;
    }
    if (std::abs(dx) < xacc) return rts;
    funcd(rts,f,df);
    if (f < 0.0)
       xl=rts;
    else
       xh=rts;
 }
 throw(std::string("Maximum number of iterations exceeded in rtsafe"));
 return 0.0;
}






class PeriodicExpFuncDeriv : public FuncAndDerivSingleVar
{
   double m_r;    //    C(t+step)/C(t)
   int m_step;    //  effective energy step
   int m_k;       //  Textent-2*t

 public:

   PeriodicExpFuncDeriv(double in_r, int in_step, int in_k)
        :  m_r(in_r), m_step(in_step), m_k(in_k) {}
   PeriodicExpFuncDeriv(const PeriodicExpFuncDeriv& in)
        :  m_r(in.m_r), m_step(in.m_step), m_k(in.m_k) {}
   PeriodicExpFuncDeriv& operator=(const PeriodicExpFuncDeriv& in)
    {m_r=in.m_r; m_step=in.m_step; m_k=in.m_k;
     return *this;}

   virtual void operator()(double s, double& f, double& df)
   {
    double s1=ipow(s,m_k-1), s2=ipow(s,m_step-1), s3=ipow(s,m_k-m_step-1);
    df=s1*m_k*m_r-s2*m_step-s3*(m_k-m_step);
    f=(1.0+s1*s)*m_r-s*(s2+s3);
   }

};


/*
class EffectiveEnergyCalculator
{
   unsigned int step;      // effective energy from C(t) and C(t+step)
   unsigned int Textent;   // temporal extent of lattice
   unsigned int type;      // 0 => A*exp(-m*t),  1 => A*(exp(-m*t)+exp(-m*(T-t)))

   mutable int m_current_k;            //  Textent-2*t
   mutable double m_current_r;         //  C(t+step)/C(t)

   EffectiveEnergyCalculator() = delete;

 public:

   EffectiveEnergyCalculator(unsigned int in_step, unsigned int in_Textent,
                             unsigned int in_type);

   EffectiveEnergyCalculator(const EffectiveEnergyCalculator& in)
      :  step(in.step), Textent(in.Textent), type(in.type) {}

   EffectiveEnergyCalculator& operator=(const EffectiveEnergyCalculator& in)
    {step=in.step; Textent=in.Textent; type=in.type;
     return *this;}

   bool calculate(int tvalue, double corr, double corrstep, double& value);

   void eval_func_and_deriv(double s, double& f, double& df);
 

};


EffectiveEnergyCalculator::EffectiveEnergyCalculator(
               unsigned int in_step, unsigned int in_Textent,
               unsigned int in_type)
      :  step(in_step), Textent(in_Textent), type(in_type)
{
 if (Textent<8) 
    throw(std::string("Invalid Textent in EffectiveEnergyCalculator"));
 if ((step<1)||(step>Textent/4))
    throw(std::string("Invalid time step in EffectiveEnergyCalculator"));
 if (type>1)
    throw(std::string("Invalid type in EffectiveEnergyCalculator"));
}


bool EffectiveEnergyCalculator::calculate(int tvalue, double corr, double corrstep, 
                                          double& value)
{
 value=-1.0;
 if (type==0){             // C(t) = A*exp(-m*t)
    double r=corrstep/corr;
    if (r<=0.0) return false;
    value=-log(r)/double(step);
    return true;}
 if (type==1){             // C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
    if ((tvalue<0)||(2*(tvalue+step)>Textent)) return false;
    m_current_r=corrstep/corr;
    if ((m_current_r<0.0)||(m_current_r>=1.0)) return false;
    m_current_k=int(Textent)-2*tvalue;

    double sb=std::pow(m_current_r,1.0/double(step)); // initial guess
    double f,df;
    eval_func_and_deriv(sb,f,df);
    cout << "f = "<<f<<"  df = "<<df<<endl;

    cout << "Start periodic eff 2: r = "<<m_current_r<<" sb = "<<sb<<endl;
    cout << "f = "<<f<<"  df = "<<df<<endl;
    double bstep=-f/df;
    double sa=sb-2*bstep;
    double fnext,dfnext;
    eval_func_and_deriv(sa,fnext,dfnext);
    cout << "f = "<<fnext<<"  df = "<<dfnext<<endl;
    

    double acc=1e-10;
    unsigned int maxit=200;
 
    void (EffectiveEnergyCalculator::*funcd)(double, double&, double&) = 
         &EffectiveEnergyCalculator::eval_func_and_deriv;

    (this->*funcd)(sa,fnext,dfnext);

  void (*fff)(double, double&, double&)=this->*funcd;

  //  &(this->periodic_eff_energy_func);

//static_cast<(void*)(double,double&,double&)>(
//          &(this->periodic_eff_energy_func);
    try{
       double s=rtsafe(*this,sa,sb,acc,maxit);
       value=-1.0/double(step)*log(s);
       return true;}
    catch(...){
       return false;}}

}




void EffectiveEnergyCalculator::eval_func_and_deriv(double s, double& f, double& df)
{
 double s1=ipow2(s,m_current_k-1), s2=ipow2(s,step-1), s3=ipow2(s,m_current_k-step-1);
 df=s1*m_current_k*m_current_r-s2*step-s3*(m_current_k-step);
 f=(1.0+s1*s)*m_current_r-s*(s2+s3);
}


double periodic_eff_energy(double corrt, double corrtplusd, unsigned int d,
                           unsigned int t, unsigned int Textent)
{
 if (2.0*double(t+d)<=double(Textent)){
    double r=corrtplusd/corrt;
    if ((r<0.0)||(r>=1.0)){
       cout << "cannot solve"<<endl; return -1.0;}
    cout << "r = "<<r<<endl;
//    double k = double(Textent-2*t)/double(d);
    double k = double(Textent-2*t)/double(d)-1.0;
    double y, s=0.0;
    double snext=r;
    
    double tol=1e-12;
    unsigned int count=0, countmax=40;
    while ((std::abs(snext-s)>tol)&&(count<countmax)){
       s=snext; count++; 
       y=std::pow(s,k);
//       snext=0.5*r*(1.0+y)*(1.0+sqrt(1.0-4.0*y/(std::pow(r*(1.0+y),2))));
       snext=(r-y)/(1.0-y*r);
    cout << "s = "<<s<<" count = "<<count<<"  snext = "<<snext<<endl;}
    double m=-1.0/double(d)*log(s);
    return m;}
 else{
    cout << "too close to midpoint"<<endl;
    return -1.0;}

}


 

double periodic_eff_energy2(double corrt, double corrtplusd, unsigned int d,
                            unsigned int t, unsigned int Textent)
{
 if (2*(t+d)<=Textent){
    double r=corrtplusd/corrt;
    if ((r<=0.0)||(r>=1.0)){
       cout << "cannot solve"<<endl; return -1.0;}
    cout << "r = "<<r<<endl;
//    double k = double(Textent-2*t)/double(d);
    unsigned int k = Textent-2*t;

    double s=std::pow(r,1.0/double(d));
    double f,df;
    periodic_eff_energy_func(s,f,df,r,k,d);
    cout << "f = "<<f<<"  df = "<<df<<endl;

    cout << "Start periodic eff 2: r = "<<r<<" s = "<<s<<endl;
    cout << "f = "<<f<<"  df = "<<df<<endl;
    double step=-f/df;
    s+=step;
    double fnext,dfnext;
    periodic_eff_energy_func(s,fnext,dfnext,r,k,d);
    cout << "f = "<<fnext<<"  df = "<<dfnext<<endl;
    
    double m=-1.0/double(d)*log(s);
    return m;}
 else{
    cout << "too close to midpoint"<<endl;
    return -1.0;}

}

Define  r = C(t+d)/C(t)

then solve        exp(-m*(t+d)) + exp(-m*(T-t-d))
            r = ----------------------------------
                      exp(-m*t) + exp(-m*(T-t))
*/

class EffectiveEnergyCalculator
{
   unsigned int step;      // effective energy from C(t) and C(t+step)
   unsigned int Textent;   // temporal extent of lattice
   unsigned int type;      // 0 => A*exp(-m*t),  1 => A*(exp(-m*t)+exp(-m*(T-t)))

   EffectiveEnergyCalculator() = delete;

 public:

   EffectiveEnergyCalculator(unsigned int in_step, unsigned int in_Textent,
                             unsigned int in_type);

   EffectiveEnergyCalculator(const EffectiveEnergyCalculator& in)
      :  step(in.step), Textent(in.Textent), type(in.type) {}

   EffectiveEnergyCalculator& operator=(const EffectiveEnergyCalculator& in)
    {step=in.step; Textent=in.Textent; type=in.type;
     return *this;}

   bool calculate(int tvalue, double corr, double corrstep, double& value);

};


EffectiveEnergyCalculator::EffectiveEnergyCalculator(
               unsigned int in_step, unsigned int in_Textent,
               unsigned int in_type)
      :  step(in_step), Textent(in_Textent), type(in_type)
{
 if (Textent<8) 
    throw(std::string("Invalid Textent in EffectiveEnergyCalculator"));
 if ((step<1)||(step>Textent/4))
    throw(std::string("Invalid time step in EffectiveEnergyCalculator"));
 if (type>1)
    throw(std::string("Invalid type in EffectiveEnergyCalculator"));
}


bool EffectiveEnergyCalculator::calculate(int tvalue, double corr, double corrstep, 
                                          double& value)
{
 value=-1.0;
 if (type==0){             // C(t) = A*exp(-m*t)
    double r=corrstep/corr;
    if (r<=0.0) return false;
    value=-log(r)/double(step);
    return true;}
 if (type==1){             // C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
    if ((tvalue<0)||(2*(tvalue+step)>Textent)) return false;  // only try for solution if t+step <= T/2
    double r=corrstep/corr;
    if ((r<0.0)||(r>=1.0)) return false;
    int k=int(Textent)-2*tvalue;
    PeriodicExpFuncDeriv funcd(r,step,k);

    double sb=std::pow(r,1.0/double(step)); // initial guess
    double f,df;
    funcd(sb,f,df);
    double bstep=-f/df;
    double sa=sb+2.0*bstep;
    double fnext,dfnext;
    funcd(sa,fnext,dfnext);
    int bcount=0;
      //  try to bracket the solution
    while (f*fnext>0){
       sa+=bstep; 
       if ((sa<=0.0)||(sa>=1.0)||(bcount>=30)) return false;  // could not bracket
       funcd(sa,fnext,dfnext);}
    double acc=1e-10;
    unsigned int maxit=200;
    try{
       double s=rtsafe(funcd,sa,sb,acc,maxit);
       value=-log(s);
       return true;}
    catch(...){
       return false;}}
}





int main(){


double A=1.23; 
double m=1.3428;
unsigned int Textent=128;
unsigned int d=2;

cout.precision(12);
EffectiveEnergyCalculator meff(d,Textent,1);
for (unsigned int t=0;t<Textent;t++){
   double corrt=A*(exp(-m*t)+exp(-m*(Textent-t)));
   double corrtplusd=A*(exp(-m*(t+d))+exp(-m*(Textent-t-d)));
   double msolve;
   if (meff.calculate(t,corrt,corrtplusd,msolve))
      cout << "t="<<t<<" msolve = "<<msolve<<endl;
   else
      cout << "t="<<" could not solve"<<endl;}
   //msolve=periodic_eff_energy2(corrt,corrtplusd,d,t,Textent);
   //cout << "t="<<t<<" msolve = "<<msolve<<endl;}

/* for (int k=0;k<32;k++){
   cout << "ipow(0.35,"<<k<<") = " <<ipow(0.35,k) <<"  "<<ipow2(0.35,k) <<endl;
   cout << "ipow(0.35,"<<-k<<") = " <<ipow(0.35,-k) <<"  "<<ipow2(0.35,-k) <<endl;} */

return 0;
}
