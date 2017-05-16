#include <vector>
#include <iostream>
#include <list>
#include <complex>
#include <stdexcept>
using namespace std;

typedef  unsigned int   uint;


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
 for (j=1;j<=int(maxit);j++) {
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






template <typename T>
class PeriodicExpFuncDeriv : public FuncAndDerivSingleVar
{
   double m_r;    //    C(t+step)/C(t)
   int m_step;    //  effective energy step
   T m_k;       //  Textent-2*t

 public:

   PeriodicExpFuncDeriv(double in_r, int in_step, T in_k)
        :  m_r(in_r), m_step(in_step), m_k(in_k) {}
   PeriodicExpFuncDeriv(const PeriodicExpFuncDeriv& in)
        :  m_r(in.m_r), m_step(in.m_step), m_k(in.m_k) {}
   PeriodicExpFuncDeriv& operator=(const PeriodicExpFuncDeriv& in)
    {m_r=in.m_r; m_step=in.m_step; m_k=in.m_k;
     return *this;}

   virtual void operator()(double b, double& f, double& df)
   {
    double b1=std::pow(b,m_k-1), b2=std::pow(b,m_step-1), b3=std::pow(b,m_k-m_step-1);
    df=b1*m_k*m_r-b2*m_step-b3*(m_k-m_step);
    f=(1.0+b1*b)*m_r-b*(b2+b3);
   }

};


template <typename T>
class PeriodicExp2FuncDeriv : public FuncAndDerivSingleVar
{
   double m_r;    //    (C(t+step)-C(t))/(C(t)-C(t-step))
   int m_step;    //  effective energy step
   T m_k;       //  Textent-2*t

 public:

   PeriodicExp2FuncDeriv(double in_r, int in_step, T in_k)
        :  m_r(in_r), m_step(in_step), m_k(in_k) {}
   PeriodicExp2FuncDeriv(const PeriodicExp2FuncDeriv& in)
        :  m_r(in.m_r), m_step(in.m_step), m_k(in.m_k) {}
   PeriodicExp2FuncDeriv& operator=(const PeriodicExp2FuncDeriv& in)
    {m_r=in.m_r; m_step=in.m_step; m_k=in.m_k;
     return *this;}

   virtual void operator()(double b, double& f, double& df)
   {
    double b1=std::pow(b,m_k-1), b2=std::pow(b,m_step-1);
    df=-(m_k+m_step)*b1*b2*b*m_r-m_step*b2+m_k*b1;
    f=(1.0-b1*b2*b*b)*m_r-b*(b2-b1);
   }

};


// ****************************************************************************
// *                                                                          *
// *   "EffectiveEnergyCalculator" is used for computing effective energies   *
// *   from diagonal temporal correlators.  The constructor takes three       *
// *   arguments:                                                             *
// *                                                                          *
// *        step => solves for energy using C(t+step), C(t),                  *
// *                       and C(t-step) when added constant present          *
// *        Textent => temporal extent of the lattice                         *
// *        type =>  0 means use C(t) = A*exp(-m*t),                          *
// *                 1 means use C(t) = A*(exp(-m*t)+exp(-m*(T-t)))           *
// *                 2 means use C(t) = A*exp(-m*t) + B0                      *
// *                 3 means use C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0      *
// *                                                                          *
// *   The key member is                                                      *
// *       bool calculate(double& value, int tvalue, double corr,             *
// *                      double corrforwardstep, double corrbackstep=0);     *
// *                                                                          *
// *   Effective energy returned in "value"; returns true if successful,      *
// *   false is returned if the effective energy could not be computed.       *
// *                                                                          *
// ****************************************************************************


class EffectiveEnergyCalculator
{
   unsigned int step;      // effective energy from C(t) and C(t+step)  (and C(t-step) for added const)
   unsigned int Textent;   // temporal extent of lattice
   unsigned int type;      // 0 => A*exp(-m*t),  1 => A*(exp(-m*t)+exp(-m*(T-t)))
                           // 2 => A*exp(-m*t) + B0,  3 = A*(exp(-m*t)+exp(-m*(T-t))) + B0

#ifndef NO_CXX11
   EffectiveEnergyCalculator() = delete;
#else
   EffectiveEnergyCalculator();
#endif

 public:

   EffectiveEnergyCalculator(unsigned int in_step, unsigned int in_Textent,
                             unsigned int in_type);

   EffectiveEnergyCalculator(const EffectiveEnergyCalculator& in)
      :  step(in.step), Textent(in.Textent), type(in.type) {}

   EffectiveEnergyCalculator& operator=(const EffectiveEnergyCalculator& in)
    {step=in.step; Textent=in.Textent; type=in.type;
     return *this;}

   bool needsBackStep() const {return (type>1);}

   bool calculate(double& value, int tvalue, double corr, double corrforwardstep, 
                  double corrbackstep=0);

   bool calculate(double& value, uint tvalue, double corr, double corrforwardstep, 
                  double corrbackstep=0);

   bool calculate(double& value, double tvalue, double corr, double corrforwardstep, 
                  double corrbackstep=0);

 private:


   bool forward_effcalc(double corr, double corrstep, uint step, double& effenergy);

   bool forward_effcalc_with_const(double corr, double corrforwardstep, double corrbackstep, 
                                   uint step, double& effenergy);

   template <typename T>
   bool timesym_effcalc(double corr, double corrstep, uint step, 
                        T tval, uint Textent,  double& effenergy);

   template <typename T>
   bool timesym_effcalc_with_const(double corr, double corrforwardstep, double corrbackstep, 
                                   uint step, T tval, uint Textent,  double& effenergy);

};


// ************************************************************



EffectiveEnergyCalculator::EffectiveEnergyCalculator(
               unsigned int in_step, unsigned int in_Textent,
               unsigned int in_type)
      :  step(in_step), Textent(in_Textent), type(in_type)
{
 if (Textent<8) 
    throw(std::invalid_argument("Invalid Textent in EffectiveEnergyCalculator"));
 if ((step<1)||(step>Textent/4))
    throw(std::invalid_argument("Invalid time step in EffectiveEnergyCalculator"));
 if (type>3)
    throw(std::invalid_argument("Invalid type in EffectiveEnergyCalculator"));
}


bool EffectiveEnergyCalculator::calculate(double& value, int tvalue, double corr, 
                                          double corrstep, double corrbackstep)
{
 value=-1.0;
 if (type==0){            
                   // C(t) = A*exp(-m*t)
    return forward_effcalc(corr,corrstep,step,value);}
 else if (type==1){             
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
    return timesym_effcalc(corr,corrstep,step,tvalue,Textent,value);}
 else if (type==2){
                   // C(t) = A*exp(-m*t) + B0
    return forward_effcalc_with_const(corr,corrstep,corrbackstep,step,value);}
 else if (type==3){
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0
    return timesym_effcalc_with_const(corr,corrstep,corrbackstep,step,tvalue,Textent,value);}
 return false;
}

bool EffectiveEnergyCalculator::calculate(double& value, uint tvalue, double corr, 
                                          double corrstep, double corrbackstep)
{
 value=-1.0;
 if (type==0){            
                   // C(t) = A*exp(-m*t)
    return forward_effcalc(corr,corrstep,step,value);}
 else if (type==1){             
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
    return timesym_effcalc(corr,corrstep,step,int(tvalue),Textent,value);}
 else if (type==2){
                   // C(t) = A*exp(-m*t) + B0
    return forward_effcalc_with_const(corr,corrstep,corrbackstep,step,value);}
 else if (type==3){
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0
    return timesym_effcalc_with_const(corr,corrstep,corrbackstep,step,int(tvalue),Textent,value);}
 return false;
}

bool EffectiveEnergyCalculator::calculate(double& value, double tvalue, double corr, 
                                          double corrstep, double corrbackstep)
{
 value=-1.0;
 if (type==0){            
                   // C(t) = A*exp(-m*t)
    return forward_effcalc(corr,corrstep,step,value);}
 else if (type==1){             
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
    return timesym_effcalc(corr,corrstep,step,tvalue,Textent,value);}
 else if (type==2){
                   // C(t) = A*exp(-m*t) + B0
    return forward_effcalc_with_const(corr,corrstep,corrbackstep,step,value);}
 else if (type==3){
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0
    return timesym_effcalc_with_const(corr,corrstep,corrbackstep,step,tvalue,Textent,value);}
 return false;
}

    //  Routines below compute the "effective energy" for different
    //  assumed correlator forms.  Result is returned in "effenergy"
    //  and routines return "true" if successful, "false" otherwise

               // C(t) = A*exp(-m*t):   corr = C(t),  corrstep = C(t+step)

bool EffectiveEnergyCalculator::forward_effcalc(
                                double corr, double corrstep, uint step, double& effenergy)
{
 double r=corrstep/corr;
 if (r<=0.0) return false;
 effenergy=-log(r)/double(step);
 return true;
}

               // C(t) = A*exp(-m*t) + B0:   corr = C(t),  corrforwardstep = C(t+step)
               //                            corrbackstep = C(t-step)

bool EffectiveEnergyCalculator::forward_effcalc_with_const(
                                double corr, double corrforwardstep, double corrbackstep, 
                                uint step, double& effenergy)
{
 double r=(corrforwardstep-corr)/(corr-corrbackstep);
 if (r<=0.0) return false;
 effenergy=-log(r)/double(step);
 return true;
}

      // C(t) = A*(exp(-m*t)+exp(-m*(T-t))):   corr = C(tval),  corrstep = C(tval+step)
      //
      // Method:  must solve       C(t+s)      A*(exp(-m*(t+s))+exp(-m*(T-t-s)))
      //                      r =  ------  =   ---------------------------------
      //                            C(t)          A*(exp(-m*t)+exp(-m*(T-t)))
      //
      //  Define b = exp(-m)  and K = T-2*t, then solve for b below:
      //
      //                   (1+b^K)*r-b^s-b^(K-s) = 0

template <typename T>
bool EffectiveEnergyCalculator::timesym_effcalc(
                                double corr, double corrstep, uint step, 
                                T tvalue, uint Textent,  double& effenergy)
{
 if ((tvalue<0)||(tvalue>=(int(Textent)-int(step)))) return false;
 double r; T k;
 if (tvalue<(int(Textent)/2)){
    k=T(Textent)-2*tvalue;
    r=corrstep/corr;}
 else{
    T tt=Textent-tvalue-step;
    k=T(Textent)-2*tt;
    r=corr/corrstep;}
 if ((r<0.0)||(r>=1.0)) return false;
 PeriodicExpFuncDeriv<T> funcd(r,step,k);

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
    bcount++;
    funcd(sa,fnext,dfnext);}
 double acc=1e-10;
 unsigned int maxit=200;
 try{
    double s=rtsafe(funcd,sa,sb,acc,maxit);
    effenergy=-log(s);
    return true;}
 catch(const std::exception& xp){
    return false;}
}



      // C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0:   corr = C(t),  corrforwardstep = C(t+step)
      //                                            corrbackstep = C(t-step)  
      //
      // Method:  must solve       C(t+s)-C(t)  
      //                      r =  ------------ 
      //                            C(t)-C(t-s) 
      //
      //  Define b = exp(-m)  and K = T-2*t,  then solve for b below:
      //
      //                   (1-b^(K+s))*r-b^s+b^K = 0

template <typename T>
bool EffectiveEnergyCalculator::timesym_effcalc_with_const(
                                double corr, double corrforwardstep, double corrbackstep, 
                                uint step, T tvalue, uint Textent,  double& effenergy)
{
 if ((tvalue<0)||(tvalue>=(int(Textent)-2*int(step)))) return false;
 double r; T k;
 if (tvalue<(int(Textent)/2)){
    k=T(Textent)-2*tvalue;
    r=(corrforwardstep-corr)/(corr-corrbackstep);}
 else{
    T tt=Textent-tvalue;
    k=T(Textent)-2*tt;
    r=(corrbackstep-corr)/(corr-corrforwardstep);}
 if ((r<0.0)||(r>=1.0)) return false;
 PeriodicExp2FuncDeriv<T> funcd(r,step,k);

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
    bcount++;
    funcd(sa,fnext,dfnext);}
 double acc=1e-10;
 unsigned int maxit=200;
 try{
    double s=rtsafe(funcd,sa,sb,acc,maxit);
    effenergy=-log(s);
    return true;}
 catch(const std::exception& xp){
    return false;}
}


// ***************************************************************************************















// ****************************************************************************
// *                                                                          *
// *   "EffectiveEnergyCalculator" is used for computing effective energies   *
// *   from diagonal temporal correlators.  The constructor takes three       *
// *   arguments:                                                             *
// *                                                                          *
// *        step => solves for energy using C(t+step), C(t),                  *
// *                       and C(t-step) when added constant present          *
// *        Textent => temporal extent of the lattice                         *
// *        type =>  0 means use C(t) = A*exp(-m*t),                          *
// *                 1 means use C(t) = A*(exp(-m*t)+exp(-m*(T-t)))           *
// *                 2 means use C(t) = A*exp(-m*t) + B0                      *
// *                 3 means use C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0      *
// *                                                                          *
// *   The key member is                                                      *
// *       bool calculate(double& value, int tvalue, double corr,             *
// *                      double corrforwardstep, double corrbackstep=0);     *
// *                                                                          *
// *   Effective energy returned in "value"; returns true if successful,      *
// *   false is returned if the effective energy could not be computed.       *
// *                                                                          *
// ****************************************************************************


class EffEnCalc
{
   unsigned int step;      // effective energy from C(t) and C(t+step)  (and C(t-step) for added const)
   unsigned int Textent;   // temporal extent of lattice
   unsigned int type;      // 0 => A*exp(-m*t),  1 => A*(exp(-m*t)+exp(-m*(T-t)))
                           // 2 => A*exp(-m*t) + B0,  3 = A*(exp(-m*t)+exp(-m*(T-t))) + B0

#ifndef NO_CXX11
   EffEnCalc() = delete;
#else
   EffEnCalc();
#endif

 public:

   EffEnCalc(unsigned int in_step, unsigned int in_Textent,
                             unsigned int in_type);

   EffEnCalc(const EffEnCalc& in)
      :  step(in.step), Textent(in.Textent), type(in.type) {}

   EffEnCalc& operator=(const EffEnCalc& in)
    {step=in.step; Textent=in.Textent; type=in.type;
     return *this;}

   bool calculate(double& value, int tvalue, double corr, double corrforwardstep, 
                  double corrbackstep=0);

 private:


   bool forward_effcalc(double corr, double corrstep, uint step, double& effenergy);

   bool forward_effcalc_with_const(double corr, double corrforwardstep, double corrbackstep, 
                                   uint step, double& effenergy);

   bool timesym_effcalc(double corr, double corrstep, uint step, 
                        uint tval, uint Textent,  double& effenergy);

   bool timesym_effcalc_with_const(double corr, double corrforwardstep, double corrbackstep, 
                                   uint step, uint tval, uint Textent,  double& effenergy);

};


EffEnCalc::EffEnCalc(
               unsigned int in_step, unsigned int in_Textent,
               unsigned int in_type)
      :  step(in_step), Textent(in_Textent), type(in_type)
{
 if (Textent<8) 
    throw(std::invalid_argument("Invalid Textent in EffEnCalc"));
 if ((step<1)||(step>Textent/4))
    throw(std::invalid_argument("Invalid time step in EffEnCalc"));
 if (type>3)
    throw(std::invalid_argument("Invalid type in EffEnCalc"));
}


bool EffEnCalc::calculate(double& value, int tvalue, double corr, 
                                          double corrstep, double corrbackstep)
{
 value=-1.0;
 if (type==0){            
                   // C(t) = A*exp(-m*t)
    return forward_effcalc(corr,corrstep,step,value);}
 else if (type==1){             
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
    return timesym_effcalc(corr,corrstep,step,tvalue,Textent,value);}
 else if (type==2){
                   // C(t) = A*exp(-m*t) + B0
    return forward_effcalc_with_const(corr,corrstep,corrbackstep,step,value);}
 else if (type==3){
                   // C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0
    return timesym_effcalc_with_const(corr,corrstep,corrbackstep,step,tvalue,Textent,value);}
 return false;
}

    //  Routines below compute the "effective energy" for different
    //  assumed correlator forms.  Result is returned in "effenergy"
    //  and routines return "true" if successful, "false" otherwise

               // C(t) = A*exp(-m*t):   corr = C(t),  corrstep = C(t+step)

bool EffEnCalc::forward_effcalc(
                                double corr, double corrstep, uint step, double& effenergy)
{
 double r=corrstep/corr;
 if (r<=0.0) return false;
 effenergy=-log(r)/double(step);
 return true;
}

               // C(t) = A*exp(-m*t) + B0:   corr = C(t),  corrforwardstep = C(t+step)
               //                            corrbackstep = C(t-step)

bool EffEnCalc::forward_effcalc_with_const(
                                double corr, double corrforwardstep, double corrbackstep, 
                                uint step, double& effenergy)
{
 double r=(corrforwardstep-corr)/(corr-corrbackstep);
 if (r<=0.0) return false;
 effenergy=-log(r)/double(step);
 return true;
}


      // C(t) = A*(exp(-m*t)+exp(-m*(T-t))):   corr = C(tval),  corrstep = C(tval+step)
      //
      // Method:  must solve       C(t+s)      A*(exp(-m*(t+s))+exp(-m*(T-t-s)))
      //                      r =  ------  =   ---------------------------------
      //                            C(t)          A*(exp(-m*t)+exp(-m*(T-t)))
      //
      //  Define b = exp(-m)  and K = T-2*t-s, then solve for b below:
      //
      //                f(b) =   (1+b^(K+s))*r-b^s-b^K = 0
      //
      //  Solve by Newton-Raphson: delta = - f(b)/ f'(b)   b+=delta 

bool EffEnCalc::timesym_effcalc(
                                double corr, double corrstep, uint step, 
                                uint tvalue, uint Textent,  double& effenergy)
{
 if ((tvalue<0)||(tvalue>=(Textent-step))) return false;
 double r; int k;
 if (tvalue<(Textent/2)){
    k=int(Textent)-2*tvalue-step;
    r=corrstep/corr;}
 else{
    int tt=Textent-tvalue-step;
    k=int(Textent)-2*tt-step;
    r=corr/corrstep;}
 if ((r<0.0)||(r>=1.0)) return false;
 double s=double(step);
 double K=double(k);
 double delta=1.0;
 double acc=1e-12;
 double b=std::pow(r,1.0/s); // initial guess
 int niter=100; int iter=0;
 while ((std::abs(delta)>acc)&&(iter<niter)){
    double bs=std::pow(b,step); 
    double bK=std::pow(b,k);
    delta=b*(bs+bK-r*(1.0+bK*bs))/(r*(K+s)*bK*bs-K*bK-s*bs);
    b+=delta; iter++;}
 if (iter>=niter) return false;
 effenergy=-log(b);
 return true;
}


      // C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0:   corr = C(t),  corrforwardstep = C(t+step)
      //                                            corrbackstep = C(t-step)  
      //
      // Method:  must solve       C(t+s)-C(t)  
      //                      r =  ------------ 
      //                            C(t)-C(t-s) 
      //
      //  Define b = exp(-m)  and K = T-2*t,  then solve for b below:
      //
      //                   (1-b^(K+s))*r-b^s+b^K = 0

bool EffEnCalc::timesym_effcalc_with_const(
                                double corr, double corrforwardstep, double corrbackstep, 
                                uint step, uint tvalue, uint Textent,  double& effenergy)
{
 if ((tvalue<0)||(tvalue>=(Textent-2*step))) return false;
 double r; int k;
 if (tvalue<(Textent/2)){
    k=int(Textent)-2*tvalue;
    r=(corrforwardstep-corr)/(corr-corrbackstep);}
 else{
    int tt=Textent-tvalue-step;
    k=int(Textent)-2*tt;
    r=(corrbackstep-corr)/(corr-corrforwardstep);}
 if ((r<0.0)||(r>=1.0)) return false;
 double s=double(step);
 double K=double(k);
 double delta=1.0;
 double acc=1e-12;
 double b=std::pow(r,1.0/s); // initial guess
 int niter=200; int iter=0;
 while ((std::abs(delta)>acc)&&(iter<niter)){
    double bs=std::pow(b,step); 
    double bK=std::pow(b,k);
    delta=b*(bs-bK-r*(1.0-bK*bs))/(r*(K+s)*bK*bs-K*bK+s*bs);
    b+=delta; iter++;}
 if (iter>=niter) return false;
 effenergy=-log(b);
 return true;
}





int main(){


double A=201.23; 
double m=1.3428;
unsigned int Textent=128;
unsigned int d=2;

cout.precision(12);
EffEnCalc meff(d,Textent,1);
EffectiveEnergyCalculator meff2(d,Textent,1);
for (unsigned int t=0;t<Textent;t++){
   double corrt=A*(exp(-m*t)+exp(-m*(Textent-t)));
   double corrtplusd=A*(exp(-m*(t+d))+exp(-m*(Textent-t-d)));
   double msolve,msolveB;
   if (meff.calculate(msolve,t,corrt,corrtplusd)){
      cout << "t="<<t<<" msolve = "<<msolve<<endl;}
   else
      cout << "t="<<t<<" could not solve"<<endl;
   if (meff2.calculate(msolve,t,corrt,corrtplusd)){
      meff2.calculate(msolveB,double(t),corrt,corrtplusd);
      cout << "t="<<t<<" msolve2 = "<<msolve<<"   "<<msolveB<<endl;}
   else
      cout << "t="<<t<<" could not solve"<<endl;
   //msolve=periodic_eff_energy2(corrt,corrtplusd,d,t,Textent);
   //cout << "t="<<t<<" msolve = "<<msolve<<endl;
   double tf=t;  tf+=0.32;
   corrt=A*(exp(-m*tf)+exp(-m*(double(Textent)-tf)));
   corrtplusd=A*(exp(-m*(tf+d))+exp(-m*(double(Textent)-tf-d)));
   if (meff2.calculate(msolve,int(tf),corrt,corrtplusd)){
      meff2.calculate(msolveB,tf,corrt,corrtplusd);
      cout << "tf="<<tf<<" msolve2 = "<<msolve<<"   "<<msolveB<<endl;}
   else
      cout << "tf="<<tf<<" could not solve"<<endl;

   }
/*
 for (int k=0;k<64;k++){
   cout << "ipow(0.35,"<<k<<") = " <<ipow(0.35,k) <<"  "<<std::pow(0.35,k) <<endl;
   cout << "ipow(0.35,"<<-k<<") = " <<ipow(0.35,-k) <<"  "<<std::pow(0.35,-k) <<endl;} 
*/



cout <<endl<<endl<< " *********************** "<<endl<<endl;


A=201.23; 
m=0.07034;
Textent=128;
double c0=-45.67;
d=2;

cout.precision(12);
EffEnCalc meff3(d,Textent,3);
EffectiveEnergyCalculator meff4(d,Textent,3);
for (unsigned int t=0;t<Textent;t++){
   double corrt=c0+A*(exp(-m*t)+exp(-m*(Textent-t)));
   double corrtplusd=c0+A*(exp(-m*(t+d))+exp(-m*(Textent-t-d)));
   double corrtbackd=c0+A*(exp(-m*(t-d))+exp(-m*(Textent-t+d)));
   double msolve,msolveB;
   cout << "-----------------------"<<endl;
   if (meff3.calculate(msolve,t,corrt,corrtplusd,corrtbackd)){
      cout << "t="<<t<<" msolve = "<<msolve<<"   "<<msolveB<<endl;}
   else
      cout << "t="<<t<<" could not solve"<<endl;
   if (meff4.calculate(msolve,t,corrt,corrtplusd,corrtbackd)){
      meff3.calculate(msolveB,double(t),corrt,corrtplusd,corrtbackd);
      cout << "t="<<t<<" msolve2 = "<<msolve<<endl;}
   else
      cout << "t="<<t<<" could not solve"<<endl;
   }


cout <<endl<<endl<< " *********************** "<<endl<<endl;

m=0.0693197469181;
A=153.199351061;
double DD=0.492212686542;
double B=0.361669109146;
d=3;

cout.precision(12);
EffectiveEnergyCalculator approach(d,Textent,1);
for (unsigned int t=0;t<Textent;t++){
   double tf=double(t);
   double tb=double(Textent)-tf;
   double corrt=A*(exp(-m*tf)*(1.0+B*exp(-DD*DD*tf))+exp(-m*tb)*(1.0+B*exp(-DD*DD*tb)));
   tf=double(t+d);
   tb=double(Textent)-tf;
   double corrtplusd=A*(exp(-m*tf)*(1.0+B*exp(-DD*DD*tf))+exp(-m*tb)*(1.0+B*exp(-DD*DD*tb)));
   double msolve,msolveB;
   if (approach.calculate(msolve,t,corrt,corrtplusd)){
      approach.calculate(msolveB,double(t),corrt,corrtplusd);
      cout << "t="<<t<<" msolve2 = "<<msolve<<"     "<<msolveB<<endl;}
   else
      cout << "t="<<" could not solve"<<endl;}



return 0;
}
