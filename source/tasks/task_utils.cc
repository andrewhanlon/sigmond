#include "task_utils.h"
#include "stopwatch.h"
using namespace std;


// *************************************************************************


bool read_arg_type(XMLHandler& xmlin, ComplexArg& arg)
{
 arg=RealPart;
 XMLHandler xmlt(xmlin);
 int na=xml_child_tag_count(xmlt,"Arg");
 if (na>1) throw(std::invalid_argument("Invalid Arg tag"));
 else if (na==0) return false;
 string reply;
 if (xmlreadifchild(xmlt,"Arg",reply)){
    if ((reply=="ImaginaryPart")||(reply=="Im")) arg=ImaginaryPart;
    else if ((reply=="RealPart")||(reply=="Re")) arg=RealPart;
    else throw(std::invalid_argument("Invalid Arg tag"));
    return true;}
 return false;
}


  //   Using a combination of Newton-Raphson and bisection, finds the
  //   root of a function f(x) bracketed between x1 and x2.  The root, returned
  //   as the function value, will be refined until its accuracy is known
  //   within +/- xacc.  "funcd" is a user-supplied routine that
  //   returns both the function value and its derivative.  The object of type T
  //   must have a member function  .eval_func_and_deriv(double,double&,double&)

double rtsafe(FuncAndDerivSingleVar& funcd, double x1, double x2, double xacc,
              unsigned int maxit)
{
 int j;
 double df,dx,dxold,f,fh,fl;
 double temp,xh,xl,rts;

 funcd(x1,fl,df);
 funcd(x2,fh,df);
 if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    throw(std::invalid_argument("Root must be bracketed in rtsafe"));
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
 throw(std::invalid_argument("Maximum number of iterations exceeded in rtsafe"));
 return 0.0;
}


// ****************************************************




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

   //  Reads temporal correlator data and returns vector time separations
   //  that have all Monte Carlo measurements available

void getCorrelatorAvailableTimes(MCObsHandler *moh,
                                 set<uint>& timesavailable,
                                 const CorrelatorInfo& corr, bool hermitian,
                                 ComplexArg arg)
{
 timesavailable.clear();
 CorrelatorAtTimeInfo corrt(corr,0,hermitian);
 for (uint tval=0;tval<moh->getLatticeTimeExtent();tval++){
    corrt.resetTimeSeparation(tval);
    if (moh->queryFullAndSamplings(MCObsInfo(corrt,arg))) timesavailable.insert(tval);}
}




void getCorrelatorEstimates(MCObsHandler *moh, const CorrelatorInfo& corr,
                  bool hermitian, bool subtract_vev, ComplexArg arg,
                  SamplingMode mode, map<double,MCEstimate>& results)
{
 results.clear();
 CorrelatorAtTimeInfo corrtv(corr,0,hermitian,subtract_vev);
 for (uint tval=0;tval<moh->getLatticeTimeExtent();tval++){
    corrtv.resetTimeSeparation(tval);
    MCObsInfo obskey(corrtv,arg);
    try{
       if (moh->queryFullAndSamplings(obskey,mode)){
          MCEstimate est=moh->getEstimate(obskey,mode);
          results.insert(make_pair(double(tval),est));}}
    catch(const std::exception& xp){}} 
}


 // ******************************************************************


#ifdef COMPLEXNUMBERS

void getHermCorrelatorMatrixAtTime_CurrentSampling(MCObsHandler *moh,
                  const CorrelatorMatrixInfo* cormat, uint timeval,
                  ComplexHermitianMatrix& cormat_estimates,
                  const CorrelatorMatrixInfo* orig_cormat,
                  const TransMatrix* orig_trans)
{
 try{
 const set<OperatorInfo>& corrops=orig_cormat->getOperators();
 uint nops=orig_cormat->getNumberOfOperators();
 bool herm=orig_cormat->isHermitian();
 if (!herm){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for this case"));}
 bool subtract_vevs=orig_cormat->subtractVEV();
 cormat_estimates.resize(nops);
 int row=0;
 for (set<OperatorInfo>::const_iterator snk=corrops.begin();snk!=corrops.end();snk++,row++){
    int col=row;
    for (set<OperatorInfo>::const_iterator src=snk;src!=corrops.end();src++,col++){
       CorrelatorAtTimeInfo corrt(*snk,*src,timeval,herm,subtract_vevs);
       MCObsInfo obskey(corrt,RealPart);
       if (src==snk){
          double cor_re=moh->getCurrentSamplingValue(obskey);  // reads data, subtracts vevs
          cormat_estimates.put(row,col,std::complex<double>(cor_re,0.0));}
       else{
          double cor_re=moh->getCurrentSamplingValue(obskey);
          obskey.setToImaginaryPart();
          double cor_im=moh->getCurrentSamplingValue(obskey);
          cormat_estimates.put(row,col,std::complex<double>(cor_re,cor_im));}}}
 if (orig_cormat!=cormat){
    doMatrixRotation(cormat_estimates,*orig_trans);
      // put bins of correlator matrix of improved operators into memory
    const set<OperatorInfo>& corriops=cormat->getOperators();
    row=0;
    for (set<OperatorInfo>::const_iterator snk=corriops.begin();snk!=corriops.end();snk++,row++){
       int col=row;
       for (set<OperatorInfo>::const_iterator src=snk;src!=corriops.end();src++,col++){
          CorrelatorAtTimeInfo corrt(*snk,*src,timeval,herm,subtract_vevs);
          MCObsInfo obskey(corrt,RealPart);
          if (src==snk){
             moh->putCurrentSamplingValue(obskey,cormat_estimates(row,col).real());
             obskey.setToImaginaryPart();
             moh->putCurrentSamplingValue(obskey,0.0);}
          else{
             moh->putCurrentSamplingValue(obskey,cormat_estimates(row,col).real());
             obskey.setToImaginaryPart();
             moh->putCurrentSamplingValue(obskey,cormat_estimates(row,col).imag());}}}}}
 catch(const std::exception& errmsg){
    cormat_estimates.clear();
    throw(std::invalid_argument(string("Error in getHermCorrelatorMatrixAtTime_CurrentSampling: ")
            +string(errmsg.what())));}
}


void getHermCorrelatorMatrixVEVs_CurrentSampling(MCObsHandler *moh,
                  const CorrelatorMatrixInfo* cormat, CVector& vevs,
                  const CorrelatorMatrixInfo* orig_cormat,
                  const TransMatrix* orig_trans)
{
 try{
 const set<OperatorInfo>& corrops=orig_cormat->getOperators();
 uint nops=orig_cormat->getNumberOfOperators();
 bool subtract_vevs=orig_cormat->subtractVEV();
 if (!subtract_vevs){
    throw(std::invalid_argument("CorrelatorMatrix must have VEV subtractions for this case"));}
 vevs.resize(nops);
 uint count=0;
 for (set<OperatorInfo>::const_iterator it=corrops.begin();it!=corrops.end();it++,count++){
    MCObsInfo obskey(*it,RealPart);
    double vev_re=moh->getCurrentSamplingValue(obskey);
    obskey.setToImaginaryPart();
    double vev_im=moh->getCurrentSamplingValue(obskey);
    vevs[count]=complex<double>(vev_re,vev_im);}
 if (orig_cormat!=cormat){
    doVectorRotation(vevs,*orig_trans);
      // put bins of vevs of improved operators into memory
    const set<OperatorInfo>& corriops=cormat->getOperators();
    uint count=0;
    for (set<OperatorInfo>::const_iterator it=corriops.begin();it!=corriops.end();it++,count++){
       MCObsInfo obskey(*it,RealPart);
       moh->putCurrentSamplingValue(obskey,vevs[count].real());
       obskey.setToImaginaryPart();
       moh->putCurrentSamplingValue(obskey,vevs[count].imag());}}}
 catch(const std::exception& errmsg){
    vevs.clear();
    throw(std::invalid_argument(string("Error in getHermCorrelatorMatrixVEVs_CurrentSampling: ")
            +string(errmsg.what())));}
}


void setImaginaryPartsToZero(ComplexHermitianMatrix& cormat_estimates)
{
 cormat_estimates.setAllImagToZero();
}


#else


void getHermCorrelatorMatrixAtTime_CurrentSampling(MCObsHandler *moh,
                  const CorrelatorMatrixInfo* cormat, uint timeval,
                  RealSymmetricMatrix& cormat_estimates,
                  const CorrelatorMatrixInfo* orig_cormat,
                  const TransMatrix* orig_trans)
{
 try{
 const set<OperatorInfo>& corrops=orig_cormat->getOperators();
 uint nops=orig_cormat->getNumberOfOperators();
 bool herm=orig_cormat->isHermitian();
 if (!herm){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for this case"));}
 bool subtract_vevs=orig_cormat->subtractVEV();
 cormat_estimates.resize(nops);
 int row=0;
 for (set<OperatorInfo>::const_iterator snk=corrops.begin();snk!=corrops.end();snk++,row++){
    int col=row;
    for (set<OperatorInfo>::const_iterator src=snk;src!=corrops.end();src++,col++){
       CorrelatorAtTimeInfo corrt(*snk,*src,timeval,herm,subtract_vevs);
       MCObsInfo obskey(corrt,RealPart);
       cormat_estimates(row,col)=moh->getCurrentSamplingValue(obskey);}}
 if (orig_cormat!=cormat){
    doMatrixRotation(cormat_estimates,*orig_trans);
      // put bins of correlator matrix of improved operators into memory
    const set<OperatorInfo>& corriops=cormat->getOperators();
    row=0;
    for (set<OperatorInfo>::const_iterator snk=corriops.begin();snk!=corriops.end();snk++,row++){
       int col=row;
       for (set<OperatorInfo>::const_iterator src=snk;src!=corriops.end();src++,col++){
          CorrelatorAtTimeInfo corrt(*snk,*src,timeval,herm,subtract_vevs);
          MCObsInfo obskey(corrt,RealPart);
          moh->putCurrentSamplingValue(obskey,cormat_estimates(row,col));}}}}
 catch(const std::exception& errmsg){
    cormat_estimates.clear();
    throw(std::invalid_argument(string("Error in getRealSymCorrelatorMatrixAtTime_CurrentSampling: ")
            +string(errmsg.what())));}
}


void getHermCorrelatorMatrixVEVs_CurrentSampling(MCObsHandler *moh,
                  const CorrelatorMatrixInfo* cormat, RVector& vevs,
                  const CorrelatorMatrixInfo* orig_cormat,
                  const TransMatrix* orig_trans)
{
 try{
 const set<OperatorInfo>& corrops=orig_cormat->getOperators();
 uint nops=orig_cormat->getNumberOfOperators();
 bool subtract_vevs=orig_cormat->subtractVEV();
 if (!subtract_vevs){
    throw(std::invalid_argument("CorrelatorMatrix must have VEV subtractions for this case"));}
 vevs.resize(nops);
 uint count=0;
 for (set<OperatorInfo>::const_iterator it=corrops.begin();it!=corrops.end();it++,count++){
    MCObsInfo obskey(*it,RealPart);
    vevs[count]=moh->getCurrentSamplingValue(obskey);}
 if (orig_cormat!=cormat){
    doVectorRotation(vevs,*orig_trans);
      // put bins of vevs of improved operators into memory
    const set<OperatorInfo>& corriops=cormat->getOperators();
    uint count=0;
    for (set<OperatorInfo>::const_iterator it=corriops.begin();it!=corriops.end();it++,count++){
       MCObsInfo obskey(*it,RealPart);
       moh->putCurrentSamplingValue(obskey,vevs[count]);}}}
 catch(const std::exception& errmsg){
    vevs.clear();
    throw(std::invalid_argument(string("Error in getRealSymCorrelatorMatrixVEVs_CurrentSampling: ")
            +string(errmsg.what())));}
}


void setImaginaryPartsToZero(RealSymmetricMatrix& cormat_estimates)
{
}

#endif


void eraseHermCorrelatorMatrixAtTime(MCObsHandler *moh,
                  const CorrelatorMatrixInfo& cormat, uint timeval)
{
 try{
 const set<OperatorInfo>& corrops=cormat.getOperators();
 bool herm=cormat.isHermitian();
 if (!herm){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for this case"));}
 for (set<OperatorInfo>::const_iterator snk=corrops.begin();snk!=corrops.end();snk++){
    for (set<OperatorInfo>::const_iterator src=snk; src!=corrops.end();src++){
       CorrelatorAtTimeInfo corrt(*snk,*src,timeval,herm,false);  // does not erase VEVs
       MCObsInfo obskey(corrt,RealPart);
       moh->eraseData(obskey);
#ifdef COMPLEXNUMBERS
       obskey.setToImaginaryPart();
       moh->eraseData(obskey);
#endif
       }}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Error in eraseHermCorrelatorMatrixAtTime: ")
            +string(errmsg.what())));}
}


void eraseHermCorrelatorMatrixVEVs(MCObsHandler *moh,
                  const CorrelatorMatrixInfo& cormat)
{
 try{
 const set<OperatorInfo>& corrops=cormat.getOperators();
 bool subtract_vevs=cormat.subtractVEV();
 if (!subtract_vevs){
    throw(std::invalid_argument("CorrelatorMatrix must have VEV subtractions for this case"));}
 for (set<OperatorInfo>::const_iterator it=corrops.begin();it!=corrops.end();it++){
    MCObsInfo obskey(*it,RealPart);
    moh->eraseData(obskey);
#ifdef COMPLEXNUMBERS
    obskey.setToImaginaryPart();
    moh->eraseData(obskey);
#endif
    }}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Error in eraseHermCorrelatorMatrixVEVs: ")
            +string(errmsg.what())));}
}



  // ***************


void getDiagonalCorrelatorsAtTimeEstimates(MCObsHandler *moh,
                  const CorrelatorMatrixInfo& cormat, uint timeval,
                  vector<MCEstimate>& corrdiag_estimates)
{
 try{
 const set<OperatorInfo>& corrops=cormat.getOperators();
 uint nops=cormat.getNumberOfOperators();
 bool herm=cormat.isHermitian();
 if (!herm){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for this case"));}
 bool subtract_vevs=cormat.subtractVEV();
 corrdiag_estimates.resize(nops);
 int row=0;
 for (set<OperatorInfo>::const_iterator snk=corrops.begin();snk!=corrops.end();snk++,row++){
    CorrelatorAtTimeInfo corrt(*snk,*snk,timeval,herm,subtract_vevs);
    MCObsInfo obskey(corrt,RealPart);
    corrdiag_estimates[row]=moh->getEstimate(obskey);}}
 catch(const std::exception& errmsg){
    corrdiag_estimates.clear();
    throw(std::invalid_argument(string("Error in getDiagonalCorrelatorsAtTimeEstimates: ")
            +string(errmsg.what())));}
}



  // *****************************************************************************



   //  Evaluates estimates for the effective energy for all available
   //  times.  Subtract VEVs if "subtract_vev" is input true.
   //  If subtract VEV is requested, an exception is thrown if the
   //  VEV date is not available.  Results are returned in "results"
   //  which is a map, with key given by time separation.  The
   //  effective energy parameters are
   //      step => solves for energy using C(t+step), C(t), and possibly C(t-step)
   //      efftype =>  0 means use C(t) = A*exp(-m*t),
   //                  1 means use C(t) = A*(exp(-m*t)+exp(-m*(T-t)))
   //                  2 means use C(t) = A*exp(-m*t) + B0,
   //                  3 means use C(t) = A*(exp(-m*t)+exp(-m*(T-t))) + B0
   //  You can also provide a constant to subtract from the correlator
   //  before the effective energy is calculated (which is most useful
   //  with efftype 0 and 1, and somewhat redundant with efftypes 2,3).


void getEffectiveEnergy(MCObsHandler *moh, const CorrelatorInfo& corr,
                  bool hermitian, bool subtract_vev, ComplexArg arg,
                  SamplingMode mode, uint step,
                  uint efftype, map<double,MCEstimate>& results,
                  double subtract_const)
{
 results.clear();
 EffectiveEnergyCalculator effcalc(step,moh->getLatticeTimeExtent(),efftype);
 string effname("EffEn_");
 effname+=currDateTimeString();
 MCObsInfo effkey(effname);
 moh->setSamplingMode(mode);

           //  get correlators into memory

 if (!subtract_vev){

    CorrelatorAtTimeInfo corrt(corr,0,hermitian,false);
    for (uint tval=0;tval<moh->getLatticeTimeExtent();tval++){
       corrt.resetTimeSeparation(tval);
       MCObsInfo obskey(corrt,arg);
       try{
          if (moh->queryBins(obskey)){
             moh->getBins(obskey);}}   // read into memory
       catch(const std::exception& xp){}} }

 else{

    set<int> tavail;
    CorrelatorAtTimeInfo corrt(corr,0,hermitian,false);
    CorrelatorAtTimeInfo corrtv(corr,0,hermitian,true);
    MCObsInfo src_re_info(corr.getSource(),RealPart);
    MCObsInfo snk_re_info(corr.getSink(),RealPart);
#ifdef COMPLEXNUMBERS
    MCObsInfo src_im_info(corr.getSource(),ImaginaryPart);
    MCObsInfo snk_im_info(corr.getSink(),ImaginaryPart);
    if (((!moh->queryBins(src_re_info))||(!moh->queryBins(src_im_info))
        ||(!moh->queryBins(snk_re_info))||(!moh->queryBins(snk_im_info)))
	&& ((!moh->queryFullAndSamplings(src_re_info))||(!moh->queryFullAndSamplings(snk_re_info))
	    ||(!moh->queryFullAndSamplings(src_im_info))||(!moh->queryFullAndSamplings(snk_im_info))))
       return;
#else
    if (((!moh->queryBins(src_re_info))||(!moh->queryBins(snk_re_info)))
	&& ((!moh->queryFullAndSamplings(src_re_info))||(!moh->queryFullAndSamplings(snk_re_info))))
      return;
#endif
    for (uint tval=0;tval<moh->getLatticeTimeExtent();tval++){
       corrt.resetTimeSeparation(tval);
       if (moh->queryBins(MCObsInfo(corrt,arg))) tavail.insert(tval);}
    for (moh->begin();!moh->end();++(*moh)){
       double vev=0.0;
       double src_re=moh->getCurrentSamplingValue(src_re_info);
       double snk_re=moh->getCurrentSamplingValue(snk_re_info);
#ifdef COMPLEXNUMBERS
       double src_im=moh->getCurrentSamplingValue(src_im_info);
       double snk_im=moh->getCurrentSamplingValue(snk_im_info);
       if (arg==RealPart)
          vev=snk_re*src_re+snk_im*src_im;
       else
          vev=snk_im*src_re-snk_re*src_im;
#else
       vev=snk_re*src_re;
#endif
       for (set<int>::const_iterator it=tavail.begin();it!=tavail.end();it++){
          corrt.resetTimeSeparation(*it);
          corrtv.resetTimeSeparation(*it);
          MCObsInfo obskey(corrt,arg);
          double corrval=moh->getCurrentSamplingValue(obskey)-vev;
          moh->putCurrentSamplingValue(MCObsInfo(corrtv,arg),corrval,true);}}}

   //  now compute the effective energy

 CorrelatorAtTimeInfo corrt(corr,0,hermitian,subtract_vev);
 CorrelatorAtTimeInfo corrtstep(corr,0,hermitian,subtract_vev);
 double effenergy;
 if (efftype<2){
    for (uint tval=0;tval<moh->getLatticeTimeExtent();tval++){
       corrt.resetTimeSeparation(tval);
       corrtstep.resetTimeSeparation(tval+step);
       MCObsInfo obskey1(corrt,arg);
       MCObsInfo obskey2(corrtstep,arg);
       effkey.resetObsIndex(tval);
       if (moh->queryFullAndSamplings(obskey1)&&moh->queryFullAndSamplings(obskey2)){
        try{
          for (moh->begin();!moh->end();++(*moh)){
             double cval1=moh->getCurrentSamplingValue(obskey1)-subtract_const;
             double cval2=moh->getCurrentSamplingValue(obskey2)-subtract_const;
             if (effcalc.calculate(effenergy,tval,cval1,cval2))
                moh->putCurrentSamplingValue(effkey,effenergy,true);
             else throw(std::runtime_error("Could not compute effective energy"));}
          MCEstimate est=moh->getEstimate(effkey);
          results.insert(make_pair(double(tval)+0.5*double(step),est));}
        catch(const std::exception& xp){}}
       moh->eraseData(effkey);}}
 else{
    CorrelatorAtTimeInfo corrtbackstep(corr,0,hermitian,subtract_vev);
    for (uint tval=step;tval<moh->getLatticeTimeExtent();tval++){
       corrt.resetTimeSeparation(tval);
       corrtstep.resetTimeSeparation(tval+step);
       corrtbackstep.resetTimeSeparation(tval-step);
       MCObsInfo obskey1(corrt,arg);
       MCObsInfo obskey2(corrtstep,arg);
       MCObsInfo obskey3(corrtbackstep,arg);
       effkey.resetObsIndex(tval);
       if (moh->queryFullAndSamplings(obskey1)&&moh->queryFullAndSamplings(obskey2)){
        try{
          for (moh->begin();!moh->end();++(*moh)){
             double cval1=moh->getCurrentSamplingValue(obskey1)-subtract_const;
             double cval2=moh->getCurrentSamplingValue(obskey2)-subtract_const;
             double cval3=moh->getCurrentSamplingValue(obskey3)-subtract_const;
             if (effcalc.calculate(effenergy,tval,cval1,cval2,cval3))
                moh->putCurrentSamplingValue(effkey,effenergy,true);
             else throw(std::runtime_error("Could not compute effective energy"));}
       MCEstimate est=moh->getEstimate(effkey);
       results.insert(make_pair(double(tval),est));}
        catch(const std::exception& xp){}}
       moh->eraseData(effkey);}}

}




// ***************************************************************************************


  // Prototypes of routines in LAPACK library--to call Fortran
  // routines from a C++ program, use extern "C" to tell the
  // compiler that the external routine is a C routine; then
  // add an underscore to the end of the routine name since
  // the routine is in Fortran.  All parameters must be passed
  // as pointers and two-dimensional arrays must follow
  // Fortran conventions.

extern "C"{
   void dpotrf_(char*,int*,double*,int*,int*);
   void dpotri_(char*,int*,double*,int*,int*);
   void dsygv_(int *itype, char *jobz, char *uplo, int *n, double *a,
               int *lda, double *b, int *ldb,
               double *w, double *work, int *lwork, int *info);
   void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
   void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
               double *w, double *work, int *lwork, int *info);
   void zheev_(char *jobz, char *uplo, int *n, double *a, int *lda,
               double *w, double *work, int *lwork, double *rwork,
               int *info);
}

// ***************************************************************


   //  Takes a Hermitian matrix "H" and returns the eigenvalues in
   //  ascending order in "eigvals" and the associated eigenvectors
   //  in the columns of "eigvecs".  Throws an exception if fails.

void Diagonalizer::diagonalize(const RealSymmetricMatrix& H, RVector& eigvals,
                               RMatrix& eigvecs, bool calceigvecs)
{
 int n=H.size();
 if (n==0){
   eigvals.clear();
   eigvecs.clear();return;}
 int lwork=5*n;
 RVector work(lwork);
 eigvals.resize(n);
 int info;
 char jobz=(calceigvecs)?'V':'N';
 char uplo='U';

    // load H (upper triangle) into matf fortran format
    //    (column major; row index changes fastest)
 vector<double> matf(n*n);
 for (int col=0;col<n;++col)
 for (int row=0;row<=col;++row)
    matf[row+n*col]=H(row,col);

 dsyev_(&jobz,&uplo,&n,&matf[0],&n,&eigvals[0],&work[0],&lwork,&info);
 if (info<0){
    throw(std::invalid_argument(" bad arguments in diagonalize"));}
 else if (info>0){
    throw(std::invalid_argument(" no convergence in diagonalize"));}

 if (calceigvecs){
    eigvecs.resize(n,n);
    for (int col=0;col<n;++col)
    for (int row=0;row<n;++row)
       eigvecs(row,col)=matf[row+n*col];}
}


void Diagonalizer::getEigenvectors(const RealSymmetricMatrix& H,
                                   RVector& eigvals, RMatrix& eigvecs)
{
 diagonalize(H,eigvals,eigvecs,true);
}


void Diagonalizer::getEigenvalues(const RealSymmetricMatrix& H,
                                  RVector& eigvals)
{
 RMatrix eigvecs;
 diagonalize(H,eigvals,eigvecs,false);
}




void Diagonalizer::diagonalize(const ComplexHermitianMatrix& H,
                               RVector& eigvals, CMatrix& eigvecs, bool calceigvecs)
{
 int n=H.size();
 if (n==0){
   eigvals.clear();
   eigvecs.clear();return;}
 int lwork=4*n;
 RVector work(2*lwork);
 RVector rwork(3*n);
 eigvals.resize(n);
 int info;
 char jobz=(calceigvecs)?'V':'N';
 char uplo='U';

    // load H (upper triangle) into matf fortran format
    //    (column major; row index changes fastest)
    //    complex stored as real,imag contiguous in fortran
 vector<double> matf(2*n*n);
 for (int col=0;col<n;col++)
 for (int row=0;row<=col;row++){
    int index=2*(row+n*col);
    const complex<double>& z=H(row,col);
    matf[index]=z.real();
    matf[index+1]=z.imag();}

 zheev_(&jobz,&uplo,&n,&matf[0],&n,&eigvals[0],&work[0],&lwork,&rwork[0],&info);
 if (info<0){
    throw(std::invalid_argument(" bad arguments in diagonalize"));}
 else if (info>0){
    throw(std::invalid_argument(" no convergence in diagonalize"));}

 if (calceigvecs){
    eigvecs.resize(n,n);
    for (int col=0;col<n;col++)
    for (int row=0;row<n;row++){
       int index=2*(row+n*col);
       eigvecs(row,col)=complex<double>(matf[index],matf[index+1]);}}
}


void Diagonalizer::getEigenvectors(const ComplexHermitianMatrix& H,
                                   RVector& eigvals, CMatrix& eigvecs)
{
 diagonalize(H,eigvals,eigvecs,true);
}


void Diagonalizer::getEigenvalues(const ComplexHermitianMatrix& H, RVector& eigvals)
{
 CMatrix eigvecs;
 diagonalize(H,eigvals,eigvecs,false);
}


// ****************************************************************

   //  Computes all the eigenvalues and the eigenvectors of a generalized
   //  eigenproblem, of the form   A*y=(lambda)*B*y.  Here, A and B are
   //  assumed to be NxN real symmetric or complex Hermitian, and both A and
   //  B must be positive semidefinite with the null space of B being
   //  entirely contained in the null space of A.  With these properties,
   //  a matrix Z exists such that A = Z Lambda Z^dagger, and if the
   //  null spaces of A and B are the same, then B = Z Z^dagger as well.
   //
   //  Let N0 be the rank of B, and NP be the rank of A, where we must
   //  have NP <= N0 <= N.  Objects of this class compute the NP eigenvalues
   //  in the diagonal matrix Lambda and the NxNP matrices X, Y, and Z
   //  which satisfy
   //
   //      Y^dag B Y = [I]_(NPxNP)    Y^dag A Y = Lambda
   //        X^dag X = [I]_(NPxNP)      X = B^(1/2) Y
   //               A = Z Lambda Z^dag
   //               B = Z Z^dag (if null(A)=null(B))
   //
   //  The matrices A and B are input only and are not destroyed.
   //  Lambda is returned as "eigvals", Y is returned as "eigvecs"
   //  which is useful for rotating the correlation matrix,
   //  X is returned as "orthovecs" whose columns are useful for level
   //  pinning, and Z is returned as "Zmat" which is useful for
   //  estimating operator overlap factors.

   //  If B is NOT positive definite, then the routine solves the
   //  eigensystem in the subspace of B which IS positive definite.
   //  Let lambda_max = the largest magnitude of the eigenvalues, then
   //  eigenvectors whose eigenvalues have magnitude smaller than
   //  lambda_max * min_inv_cond_num are removed. "min_inv_cond_num" is
   //  the minimum inverse condition number.  Recall that the condition
   //  number is the magnitude of the ratio of the largest eigenvalue
   //  over the smallest eigenvalue. If "A" restricted to the positive
   //  definite subspace of "B" is also NOT positive definite, then the
   //  eigenvectors associated with the negative (or small, based on
   //  min_inv_cond_num) eigenvalues are also discarded.

   //  The class checks to see if the null space of B is entirely
   //  contained in the null space of A.  If this is not true,
   //  the X and Y matrices are still correct, but the Z matrix
   //  will be incorrect.  If exceptions are turned on using
   //  setExceptionsOn(), an object of this class throws an exception
   //  if any warning or error is encountered.   With exceptions off,
   //  no exceptions are thrown, but empty results are returned.

   //  Usage:
   //      HermDiagonalizerWithMetric DG(min_inv_condnum);
   //      ComplexHermitianMatrix A, B;
   //      DG.setMetric(B);
   //      DG.getMetricEigenvalues(RVector& Beigvals);
   //      DG.getMetrixRank();
   //      DG.setMatrix(A);
   //      DG.getMatrixRank();
   //      DG.isNullBInNullA();
   //      DG.getEigenvalues(RVector& eigvals);
   //      DG.getEigenvectors(CMatrix& eigvecs);
   //      DG.getOrthovectors(CMatrix& orthovecs);
   //      DG.getZMatrix(CMatrix& Zmat);

   //  Two crucial routines are "setMetric" and "setMatrix".
   //  These are the routines which do the actual diagonalization.
   //  Each returns an integer code which summarizes the success
   //  of the diagonalization.

   //  setMetric(B):
   //    Sets the metric B.  Returns 0 if successful,  -1 if B is not
   //    positive semidefinite, -2 if diagonalizing B failed for some reason,
   //    -3 if B is trivial or the null space is the dimension of B.

   //  setMatrix(A):
   //    Sets the matrix A.  Diagonalizes G = B^(-1/2) A B^(-1/2), checks for small
   //    and negative eigenvalues.  Returns 0 if successful,  -1 if B is not
   //    set, -2 if size of A not same as B, -3 if the null space
   //    is the dimension of A, -4 if diagonalization failed for some reason,
   //    -5 if the null space of A does not contain the entire null space of B,
   //    -6 if A is not positive semidefinite



  //  Method:
  //
  //     Let A and B be n x n Hermitian matrices.
  //     Begin by writing the metric  B = U0 LB U0^dag, where U0 = unitary matrix
  //     and LB is diagonal.  The columns of U0 are the eigenvectors of B.
  //     LB is diagonal with ascending values, so the first few might be
  //     be negative or small.  We want to work in the subspace spanned by
  //     the "good" eigenvectors of B with positive and sufficiently large
  //     eigenvalues.  Call the lower right n0 x n0 square of LB as  Btilde.
  //     Put the good eigenvectors (columns of U0) into the columns of P0.
  //     P0 has n rows and n0 columns.
  //
  //     Now consider  Btilde = P0^dag B P0  and   Atilde = P0^dag A P0, which are
  //     now n0 x n0 Hermitian matrices, and Btilde is guaranteed to be
  //     positive definite.   In fact, Btilde is diagonal.  We then compute
  //              Gtilde = Btilde^(-1/2) Atilde Btilde^(-1/2)
  //     then solve   Gtilde Xtilde = Xtilde D  where Xtilde is unitary and
  //     D are the eigenvalues.  The columns of Xtilde are the orthonormal
  //     eigenvectors returned by the LAPACK solver.  We then compute
  //
  //         orthovecs:   X = P0 * Xtilde
  //         eigvecs:     Y = P0 * Btilde^(-1/2) * Xtilde
  //         Zmat         Z = P0 * Btilde^(1/2) * Xtilde




HermDiagonalizerWithMetric::HermDiagonalizerWithMetric()
   : mininvcondnum(0.0), n(0), n0(0), np(0), xon(true), Bset(false),
     Aset(false), nullB_in_nullA(false), negeigalarm(0.0)
{}


HermDiagonalizerWithMetric::HermDiagonalizerWithMetric(double min_inv_cond_num,
                                                       double negative_eigval_alarm)
   : mininvcondnum(min_inv_cond_num), n(0), n0(0), np(0),
     xon(true), Bset(false), Aset(false), nullB_in_nullA(false),
     negeigalarm(negative_eigval_alarm)
{
 setMinInvCondNum(min_inv_cond_num);
 setNegativeEigenvalueAlarm(negative_eigval_alarm);
}


HermDiagonalizerWithMetric::HermDiagonalizerWithMetric(double min_inv_cond_num)
   : mininvcondnum(min_inv_cond_num), n(0), n0(0), np(0),
     xon(true), Bset(false), Aset(false), nullB_in_nullA(false),
     negeigalarm(0.0)
{
 setMinInvCondNum(min_inv_cond_num);
 setNegativeEigenvalueAlarm(-5.0*min_inv_cond_num);
}


HermDiagonalizerWithMetric::~HermDiagonalizerWithMetric()
{
 clear();
}


void HermDiagonalizerWithMetric::clear()
{
 matb.clear();
 matg.clear();
 Beigvals.clear();
 Geigvals.clear();
 n=0; n0=0; np=0;
 Bset=false;
 Aset=false;
 nullB_in_nullA=false;
}


void HermDiagonalizerWithMetric::clearMatrix()
{
 matg.clear();
 Geigvals.clear();
 np=0;
 Aset=false;
 nullB_in_nullA=false;
}



void HermDiagonalizerWithMetric::setMinInvCondNum(double min_inv_cond_num)
{
// clear();
 if (min_inv_cond_num>=0.0)
    mininvcondnum=min_inv_cond_num;
 else
    if (xon) throw(std::invalid_argument("Min inv cond number must not be negative in HermDiagonalizerWithMetric"));
}


void HermDiagonalizerWithMetric::setNegativeEigenvalueAlarm(
            double negative_eigval_alarm)
{
// clear();
 if (negative_eigval_alarm<=0.0)
    negeigalarm=negative_eigval_alarm;
 else
    if (xon) throw(std::invalid_argument("Negative eigenvalue alarm must not be positive in HermDiagonalizerWithMetric"));
}


void HermDiagonalizerWithMetric::setExceptionsOn()
{
 xon=true;
}


void HermDiagonalizerWithMetric::setExceptionsOff()
{
 xon=false;
}


     //  Sets the metric B.  Diagonalizes B, checks for small and negative
     //  eigenvalues.  Returns 0 if successful,  -1 if B is not
     //  positive semidefinite, -2 if diagonalizing B failed for some reason,
     //  -3 if B is trivial or the null space is the dimension of B.
     //  This routine computes "matb" which contains the eigenvectors of B,
     //  "Beigvals" containing the eigenvalues of B, "n" the size of B,
     //  "n0" the rank of B, and "Bset".

int HermDiagonalizerWithMetric::setMetric(const ComplexHermitianMatrix& B,
                                          LogHelper& xmlout)
{
 clear();
 n=B.size();
 if (n==0) return -3;

 int info;
 char jobz='V';
 char uplo='U';
 int lwork=4*n;
 Beigvals.resize(n);
 RVector work(2*lwork);
 RVector rwork(3*n);

    // load B (upper triangle) into matb fortran format
    //    (column major; row index changes fastest)
 matb.resize(2*n*n);
 for (int col=0;col<n;col++)
 for (int row=0;row<=col;row++){
    int index=2*(row+n*col);
    const complex<double>& z=B(row,col);
    matb[index]=z.real();
    matb[index+1]=z.imag();}

    // solve for eigenvectors and eigenvalues of Hermitian B
    // eigenvectors returned in matb, eigenvalues in Beigvals
 zheev_(&jobz,&uplo,&n,&matb[0],&n,&Beigvals[0],&work[0],&lwork,&rwork[0],&info);

 if (info<0){
    clear();
    if (xon) throw(std::invalid_argument(" bad arguments in HermDiagonalizerWithMetric::setMetric"));
    else return -2;}
 else if (info>0){
    clear();
    if (xon) throw(std::invalid_argument(" no convergence in HermDiagonalizerWithMetric::setMetric"));
    else return -2;}

 xmlout.reset("AnalyzeMetric");
 xmlout.putInt("NumberOfOperators",n);
 xmlout.putString("ReturnCode",getRotateMetricCode(info));
 if (info==0){
    LogHelper xmleig("MetricAllEigenvalues");
    for (uint k=0;k<Beigvals.size();k++)
       xmleig.putReal("Value",Beigvals[k]);
    xmlout.putItem(xmleig);}

   // Beigvals returned in ascending order. Beigvals[n-1] is largest eigenvalue.
   // Eigenvectors whose eigenvalues have magnitude smaller than
   //     Beigvals[n-1] * min_inv_cond_num are removed. "min_inv_cond_num" is
   //  the minimum inverse condition number.  

 double cutoff=std::abs(mininvcondnum*Beigvals[n-1]);
// if (Beigvals[n-1]<mininvcondnum) cutoff=mininvcondnum;
 if (Beigvals[0]<negeigalarm){
    clear();
    xmlout.putString("Error","Metric not positive semidefinite in HermDiagonalizerWithMetric::setMetric");
    if (xon) throw(std::invalid_argument(" Metric not positive semidefinite in HermDiagonalizerWithMetric::setMetric"));
    else return -1;}

   // discard eigenvectors associated with small eigenvalues,
   // error condition if negative eigenvalues
 int Bnum_removed=0;
 while ((Bnum_removed<n)&&(Beigvals[Bnum_removed]<cutoff)) Bnum_removed++;
 if (Bnum_removed==n){
    clear();
    xmlout.putString("Warning","Null space is dim of Metric in HermDiagonalizerWithMetric::setMetric");
    if (xon) throw(std::invalid_argument("Null space is dim of Metric in HermDiagonalizerWithMetric::setMetric"));
    else return -3;}
 n0=n-Bnum_removed;
 Bset=true;
 xmlout.putInt("MetricRank",n0);
 LogHelper xmlmet("MetricRetainedEigenvalues");
 for (int k=0;k<n0;k++)
    xmlmet.putReal("Value",Beigvals[k+Bnum_removed]);
 xmlout.put(xmlmet);
 return 0;
}


int HermDiagonalizerWithMetric::setMetric(const ComplexHermitianMatrix& B)
{
 LogHelper xmlout;
 return setMetric(B,xmlout);
}


void HermDiagonalizerWithMetric::getMetricEigenvalues(RVector& metric_eigvals)
{
 if (!Bset){
    metric_eigvals.clear();
    if (xon) throw(std::invalid_argument("Metric not set in HermDiagonalizerWithMetric::getMetricEigenvalues"));
    return;}
 metric_eigvals=Beigvals;
}


     //  Sets the matrix A.  Diagonalizes G = B^(-1/2) A B^(-1/2), checks for small
     //  and negative eigenvalues.  Returns 0 if successful,  -1 if B is not
     //  set, -2 if size of A not same as B, -3 if the null space
     //  is the dimension of A, -4 if diagonalization failed for some reason,
     //  -5 if the null space of A does not contain the entire null space of B,
     //  -6 if A is not positive semidefinite

int HermDiagonalizerWithMetric::setMatrix(const ComplexHermitianMatrix& A,
                                          LogHelper& xmlout, bool checkNullSpace)
{
 clearMatrix();
 if (!Bset){
    if (xon) throw(std::invalid_argument("cannot set Matrix since Metric NOT set in HermDiagonalizerWithMetric"));
    return -1;}
 if (int(A.size())!=n){
    if (xon) throw(std::invalid_argument("Matrices A and B must be same size in HermDiagonalizerWithMetric::setMatrix"));
    return -2;}
 xmlout.reset("AnalyzeMatrix");

   // check that null space of B is entirely contained with the
   // null space of A:
   //  For each discarded (null space) eigenvector |d> of B, we want to
   //  check that |d> can be written as a linear superposition of
   //  the null space eigenvectors |a> of A.  In other words, that
   //                sum_a <d|a><a|d> ~ 1.0

 int Bnulldim=n-n0;
 if ((Bnulldim>0)&&(checkNullSpace)){
    LogHelper xmlnull("CheckNullSpaceCommonality");
    Diagonalizer DG;
    RVector Aeigvals;
    CMatrix Aeigvecs;
    DG.getEigenvectors(A,Aeigvals,Aeigvecs);
    LogHelper xmla("AMatrixEigenvalues");
    for (int k=0;k<n;k++){
       xmla.putReal("Value",Aeigvals[k]);}
    xmlnull.put(xmla);
    vector<CVector> Anullvecs;
    double Acutoff=std::abs(mininvcondnum*Aeigvals[n-1]);
    int Anulldim=0;
    while ((Anulldim<n)&&(Aeigvals[Anulldim]<Acutoff)){
       Anulldim++;
       CVector Atemp(n);
       for (int j=0;j<n;j++) Atemp[j]=Aeigvecs(j,Anulldim);
       Anullvecs.push_back(Atemp);}
    xmlnull.putUInt("DimensionNullSpaceAMatrix",Anulldim);
    if (Anulldim<Bnulldim)
       xmlnull.putString("WARNING","Anull smaller than B null");
    Aeigvecs.clear();
    if (Anulldim==n){
       clear();
       xmlout.putString("Warning","Null space is dim of Matrix in HermDiagonalizerWithMetric::setMetric");
       if (xon) throw(std::invalid_argument("Null space is dim of Matrix in HermDiagonalizerWithMetric::setMetric"));
       else return -3;}
    double ovmin=1.0;
    for (int bvec=0;bvec<Bnulldim;bvec++){
       CVector Bnullvec(n);
       int ind=2*n*bvec;
       for (int j=0;j<n;j++){
          double q0r=matb[ind++];
          double q0i=matb[ind++];
          Bnullvec(j)=complex<double>(q0r,q0i);}
       LogHelper xmlbvec("MetricNullEigenvectorTest");
       xmlbvec.putUInt("Index",bvec);
       CVector Anullvec;
       multiply(Anullvec,A,Bnullvec); // Anullvec = A * nullvec
       xmlbvec.putReal("MagnitudeOfMatrixElementOfAMatrix",
                   dotProductMagnitude(Bnullvec,Anullvec));
       double ov=0.0;
       for (int k=0;k<Anulldim;k++){
          ov+=dotProductMagnitudeSquared(Anullvecs[k],Bnullvec);}
       xmlbvec.putReal("NormProjectionIntoAMatrixNullSpace",ov);
       xmlnull.put(xmlbvec);
       if (ov<ovmin) ovmin=ov;}
    nullB_in_nullA=(ovmin>0.9);
    if (nullB_in_nullA)
       xmlnull.putString("MetricNullSpace"," SubsetOfMatrixNullSpace ");
    else
       xmlnull.putString("MetricNullSpace"," NOTSubsetOfMatrixNullSpace ");
    xmlout.put(xmlnull);}


 int info;
 char jobz='V';
 char uplo='U';
 int lwork=4*n;
 RVector work(2*lwork);
 RVector rwork(3*n);
 int Bnull=n-n0;
 RVector Btildeinvsqrt(n0);
 RVector ev(n0);
 for (int i=0;i<n0;i++)
    Btildeinvsqrt[i]=1.0/sqrt(Beigvals[i+Bnull]);

     // make the matrix  matg = Btilde^(-1/2) Atilde Btilde^(-1/2)
 matg.resize(2*n0*n0);
 for (int col=0;col<n0;++col){
    int fcol=col+Bnull;
    for (int row=0;row<=col;++row){
       int frow=row+Bnull;
       double tmpr=0.0,tmpi=0.0;
       for (int k=0;k<n;k++){
          int ind1=2*(k+n*fcol);
          double br=matb[ind1], bi=matb[ind1+1];
          for (int l=0;l<n;l++){
             int ind2=2*(l+n*frow);
             double ar=matb[ind2], ai=matb[ind2+1];
             double abr=ar*br+ai*bi;
             double abi=ar*bi-ai*br;
             double Are=real(A(l,k)),Aim=imag(A(l,k));
             tmpr+=abr*Are-abi*Aim;
             tmpi+=abr*Aim+abi*Are;}}
       int index=2*(row+n0*col);
       matg[index]=tmpr*Btildeinvsqrt[row]*Btildeinvsqrt[col];
       matg[index+1]=tmpi*Btildeinvsqrt[row]*Btildeinvsqrt[col];}}

     // diagonalize, orthonormal eigenvectors in columns of matg,
     // eigenvalues in "ev"
 zheev_(&jobz,&uplo,&n0,&matg[0],&n0,&ev[0],&work[0],&lwork,&rwork[0],&info);

 if (info<0){
    clearMatrix();
    if (xon) throw(std::invalid_argument(" bad arguments in HermDiagonalizerWithMetric::setMatrix"));
    return -4;}
 else if (info>0){
    clearMatrix();
    if (xon) throw(std::invalid_argument(" no convergence in HermDiagonalizerWithMetric::setMatrix"));
    return -4;}

 xmlout.putString("ReturnCode",getRotateMatrixCode(info));
 if ((info==0)||(info==-5)){
    LogHelper xmleig("GMatrixAllEigenvalues");
    for (uint k=0;k<ev.size();k++)
       xmleig.putReal("Value",ev[k]);
    xmlout.putItem(xmleig);}

 int Aremove=0;
 double cutoff=std::abs(mininvcondnum*ev[n0-1]);
//  if (ev[n0-1]<mininvcondnum) cutoff=mininvcondnum;
 if (ev[0]<negeigalarm){
   clearMatrix();
    if (xon) throw(std::invalid_argument(" A Matrix not positive semidefinite in HermDiagonalizerWithMetric::setMatrix"));
    else return -6;}
 if(strip_eigenvectors) while ((Aremove<n0)&&(ev[Aremove]<cutoff)) Aremove++;
 invcondnum = ev[Aremove]/ev[n0-1];
 if (Aremove==n0){
    clearMatrix();
    if (xon) throw(std::invalid_argument("Null space is dim of Matrix in HermDiagonalizerWithMetric::setMatrix"));
    else return -3;}

   //  put final retained eigenvalues in "eigvals"
 np=n0-Aremove;
 Geigvals.resize(np);
 for (int l=0;l<np;l++){
    Geigvals[np-l-1]=ev[l+Aremove];}
 xmlout.putInt("GMatrixRank",np);
 LogHelper xmlmat("GMatrixRetainedEigenvalues");
 for (int k=0;k<np;k++)
    xmlmat.putReal("Value",Geigvals[k]);
 xmlout.put(xmlmat);
 Aset=true;

 if (!checkNullSpace) return 0;
 return (nullB_in_nullA)?0:-5;
}




int HermDiagonalizerWithMetric::setMatrix(const ComplexHermitianMatrix& A,
                                          bool checkNullSpace)
{
 LogHelper xmlout;
 return setMatrix(A,xmlout,checkNullSpace);
}





void HermDiagonalizerWithMetric::getEigenvalues(RVector& eigvals)
{
 if (!Aset){
    eigvals.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for HermDiagonalizerWithMetric::getEigenvalues"));
    return;}
 eigvals=Geigvals;
}


void HermDiagonalizerWithMetric::getEigenvectors(CMatrix& eigvecs)
{
 if (!Aset){
    eigvecs.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for HermDiagonalizerWithMetric::getEigenvectors"));
    return;}

 int Bnull=n-n0;
 int Aremove=n0-np;
 RVector Btildeinvsqrt(n0);
 for (int i=0;i<n0;i++)
    Btildeinvsqrt[i]=1.0/sqrt(Beigvals[i+Bnull]);

   //  compute the generalized eigenvectors P0 Btilde^(-1/2) matg

 eigvecs.resize(n,np);
 for (int l=0;l<np;l++)
 for (int i=0;i<n;i++){
    double tmpr=0.0, tmpi=0.0; int jj=Bnull;
    for (int j=0;j<n0;j++){
       int ind1=2*(i+n*jj);
       double ar=matb[ind1],ai=matb[ind1+1];
       int ind2=2*(j+n0*(l+Aremove));
       double br=matg[ind2],bi=matg[ind2+1];
       tmpr+=(ar*br-ai*bi)*Btildeinvsqrt[j];
       tmpi+=(ar*bi+ai*br)*Btildeinvsqrt[j];
       jj++;}
    eigvecs.put(i,np-l-1,complex<double>(tmpr,tmpi));}
}



void HermDiagonalizerWithMetric::getOrthovectors(CMatrix& orthovecs)
{
 if (!Aset){
    orthovecs.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for HermDiagonalizerWithMetric::getOrthovectors"));
    return;}

 int Bnull=n-n0;
 int Aremove=n0-np;

   //  compute the orthonormal eigenvectors P0 matg

 orthovecs.resize(n,np);
 for (int l=0;l<np;l++)
 for (int i=0;i<n;i++){
    double tmpr=0.0, tmpi=0.0; int jj=Bnull;
    for (int j=0;j<n0;j++){
       int ind1=2*(i+n*jj);
       double ar=matb[ind1],ai=matb[ind1+1];
       int ind2=2*(j+n0*(l+Aremove));
       double br=matg[ind2],bi=matg[ind2+1];
       tmpr+=(ar*br-ai*bi);
       tmpi+=(ar*bi+ai*br);
       jj++;}
    orthovecs.put(i,np-l-1,complex<double>(tmpr,tmpi));}
}


void HermDiagonalizerWithMetric::getZMatrix(CMatrix& Zmat)
{
 if (!Aset){
    Zmat.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for HermDiagonalizerWithMetric::getZMatrix"));
    return;}

 int Bnull=n-n0;
 int Aremove=n0-np;
 RVector Btildesqrt(n0);
 for (int i=0;i<n0;i++)
    Btildesqrt[i]=sqrt(Beigvals[i+Bnull]);

   //  compute Zmat = P0 Btilde^(1/2) matg

 Zmat.resize(n,np);
 for (int l=0;l<np;l++)
 for (int i=0;i<n;i++){
    double tmpr=0.0, tmpi=0.0; int jj=Bnull;
    for (int j=0;j<n0;j++){
       int ind1=2*(i+n*jj);
       double ar=matb[ind1],ai=matb[ind1+1];
       int ind2=2*(j+n0*(l+Aremove));
       double br=matg[ind2],bi=matg[ind2+1];
       tmpr+=(ar*br-ai*bi)*Btildesqrt[j];
       tmpi+=(ar*bi+ai*br)*Btildesqrt[j];
       jj++;}
    Zmat.put(i,np-l-1,complex<double>(tmpr,tmpi));}
}


string HermDiagonalizerWithMetric::getRotateMetricCode(int info)
{
 if (info==0) return string(" Success ");
 else if (info==-1) return string(" Not positive semidefinite ");
 else if (info==-2) return string(" Diagonalization failed ");
 else if (info==-3) return string(" Trivial or total null space ");
 else return string(" Unknown ");
}

string HermDiagonalizerWithMetric::getRotateMatrixCode(int info)
{
 if (info==0) return string(" Success ");
 else if (info==-1) return string(" Metric not set ");
 else if (info==-2) return string(" Metric and matrix size mismatch ");
 else if (info==-3) return string(" Trivial or total null space ");
 else if (info==-4) return string(" Diagonalization failed ");
 else if (info==-5) return string(" Null space of metric not subset of null space of matrix ");
 else if (info==-6) return string(" Not positive semidefinite ");
 else return string(" Unknown ");
}




// ******************************************************



RealSymDiagonalizerWithMetric::RealSymDiagonalizerWithMetric()
   : mininvcondnum(0.0), n(0), n0(0), np(0), xon(true), Bset(false),
     Aset(false), nullB_in_nullA(false), negeigalarm(0.0)
{}


RealSymDiagonalizerWithMetric::RealSymDiagonalizerWithMetric(double min_inv_cond_num,
                                                             double negative_eigval_alarm)
   : mininvcondnum(min_inv_cond_num), n(0), n0(0), np(0),
     xon(true), Bset(false), Aset(false), nullB_in_nullA(false),
     negeigalarm(negative_eigval_alarm)
{
 setMinInvCondNum(min_inv_cond_num);
 setNegativeEigenvalueAlarm(negative_eigval_alarm);
}


RealSymDiagonalizerWithMetric::RealSymDiagonalizerWithMetric(double min_inv_cond_num)
   : mininvcondnum(min_inv_cond_num), n(0), n0(0), np(0),
     xon(true), Bset(false), Aset(false), nullB_in_nullA(false),
     negeigalarm(0.0)
{
 setMinInvCondNum(min_inv_cond_num);
 setNegativeEigenvalueAlarm(-5.0*min_inv_cond_num);
}


RealSymDiagonalizerWithMetric::~RealSymDiagonalizerWithMetric()
{
 clear();
}


void RealSymDiagonalizerWithMetric::clear()
{
 matb.clear();
 matg.clear();
 Beigvals.clear();
 Geigvals.clear();
 n=0; n0=0; np=0;
 Bset=false;
 Aset=false;
 nullB_in_nullA=false;
}


void RealSymDiagonalizerWithMetric::clearMatrix()
{
 matg.clear();
 Geigvals.clear();
 np=0;
 Aset=false;
 nullB_in_nullA=false;
}



void RealSymDiagonalizerWithMetric::setMinInvCondNum(double min_inv_cond_num)
{
// clear();
 if (min_inv_cond_num>=0.0)
    mininvcondnum=min_inv_cond_num;
 else
    if (xon) throw(std::invalid_argument("Min inv cond number must not be negative in RealSymDiagonalizerWithMetric"));
}


void RealSymDiagonalizerWithMetric::setNegativeEigenvalueAlarm(
            double negative_eigval_alarm)
{
// clear();
 if (negative_eigval_alarm<=0.0)
    negeigalarm=negative_eigval_alarm;
 else
    if (xon) throw(std::invalid_argument("Negative eigenvalue alarm must not be positive in RealSymDiagonalizerWithMetric"));
}


void RealSymDiagonalizerWithMetric::setExceptionsOn()
{
 xon=true;
}


void RealSymDiagonalizerWithMetric::setExceptionsOff()
{
 xon=false;
}


     //  Sets the metric B.  Diagonalizes B, checks for small and negative
     //  eigenvalues.  Returns 0 if successful,  -1 if B is not
     //  positive semidefinite, -2 if diagonalizing B failed for some reason,
     //  -3 if B is trivial or the null space is the dimension of B.
     //  This routine computes "matb" which contains the eigenvectors of B,
     //  "Beigvals" containing the eigenvalues of B, "n" the size of B,
     //  "n0" the rank of B, and "Bset".

int RealSymDiagonalizerWithMetric::setMetric(const RealSymmetricMatrix& B,
                                             LogHelper& xmlout)
{
 clear();
 n=B.size();
 if (n==0) return -3;

 int info;
 char jobz='V';
 char uplo='U';
 int lwork=5*n;
 RVector work(lwork);
 Beigvals.resize(n);

    // load B (upper triangle) into matb fortran format
    //    (column major; row index changes fastest)
 matb.resize(n*n);
 for (int col=0;col<n;++col)
 for (int row=0;row<=col;++row)
    matb[row+n*col]=B(row,col);

    // solve for eigenvectors and eigenvalues of Hermitian B
    // eigenvectors returned in matb, eigenvalues in Beigvals
 dsyev_(&jobz,&uplo,&n,&matb[0],&n,&Beigvals[0],&work[0],&lwork,&info);

 if (info<0){
    clear();
    if (xon) throw(std::invalid_argument(" bad arguments in RealSymDiagonalizerWithMetric::setMetric"));
    else return -2;}
 else if (info>0){
    clear();
    if (xon) throw(std::invalid_argument(" no convergence in RealSymDiagonalizerWithMetric::setMetric"));
    else return -2;}

 xmlout.reset("AnalyzeMetric");
 xmlout.putInt("NumberOfOperators",n);
 xmlout.putString("ReturnCode",getRotateMetricCode(info));
 if (info==0){
    LogHelper xmleig("MetricAllEigenvalues");
    for (uint k=0;k<Beigvals.size();k++)
       xmleig.putReal("Value",Beigvals[k]);
    xmlout.putItem(xmleig);}

   // Beigvals returned in ascending order
 double cutoff=std::abs(mininvcondnum*Beigvals[n-1]);
 if (Beigvals[n-1]<mininvcondnum) cutoff=mininvcondnum;
 if (Beigvals[0]<negeigalarm){
    clear();
    if (xon) throw(std::invalid_argument(" Metric not positive semidefinite in RealSymDiagonalizerWithMetric::setMetric"));
    else return -1;}

   // discard eigenvectors associated with small eigenvalues,
   // error condition if negative eigenvalues
 int Bnum_removed=0;
 while ((Bnum_removed<n)&&(Beigvals[Bnum_removed]<cutoff)) Bnum_removed++;
 if (Bnum_removed==n){
    clear();
    if (xon) throw(std::invalid_argument("Null space is dim of Metric in RealSymDiagonalizerWithMetric::setMetric"));
    else return -3;}
 n0=n-Bnum_removed;
 Bset=true;
 xmlout.putInt("MetricRank",n0);
 LogHelper xmlmet("MetricRetainedEigenvalues");
 for (int k=0;k<n0;k++)
    xmlmet.putReal("Value",Beigvals[k+Bnum_removed]);
 xmlout.put(xmlmet);
 return 0;
}


int RealSymDiagonalizerWithMetric::setMetric(const RealSymmetricMatrix& B)
{
 LogHelper xmlout;
 return setMetric(B,xmlout);
}


void RealSymDiagonalizerWithMetric::getMetricEigenvalues(RVector& metric_eigvals)
{
 if (!Bset){
    metric_eigvals.clear();
    if (xon) throw(std::invalid_argument("Metric not set in RealSymDiagonalizerWithMetric::getMetricEigenvalues"));
    return;}
 metric_eigvals=Beigvals;
}


     //  Sets the matrix A.  Diagonalizes G = B^(-1/2) A B^(-1/2), checks for small
     //  and negative eigenvalues.  Returns 0 if successful,  -1 if B is not
     //  set, -2 if size of A not same as B, -3 if the null space
     //  is the dimension of A, -4 if diagonalization failed for some reason,
     //  -5 if the null space of A does not contain the entire null space of B,
     //  -6 if A is not positive semidefinite

int RealSymDiagonalizerWithMetric::setMatrix(const RealSymmetricMatrix& A,
                                             LogHelper& xmlout, bool checkNullSpace)
{
 clearMatrix();
 if (!Bset){
    if (xon) throw(std::invalid_argument("cannot set Matrix since Metric NOT set in RealSymDiagonalizerWithMetric"));
    return -1;}
 if (int(A.size())!=n){
    if (xon) throw(std::invalid_argument("Matrices A and B must be same size in RealSymDiagonalizerWithMetric::setMatrix"));
    return -2;}
 xmlout.reset("AnalyzeMatrix");

   // check that null space of B is entirely contained with the
   // null space of A:
   //  For each discarded (null space) eigenvector |d> of B, we want to
   //  check that |d> can be written as a linear superposition of
   //  the null space eigenvectors |a> of A.  In other words, that
   //                sum_a <d|a><a|d> ~ 1.0

 int Bnulldim=n-n0;
 if ((Bnulldim>0)&&(checkNullSpace)){
    LogHelper xmlnull("CheckNullSpaceCommonality");
    Diagonalizer DG;
    RVector Aeigvals;
    RMatrix Aeigvecs;
    DG.getEigenvectors(A,Aeigvals,Aeigvecs);
    LogHelper xmla("AMatrixEigenvalues");
    for (int k=0;k<n;k++){
       xmla.putReal("Value",Aeigvals[k]);}
    xmlnull.put(xmla);
    vector<RVector> Anullvecs;
    double Acutoff=std::abs(mininvcondnum*Aeigvals[n-1]);
    int Anulldim=0;
    while ((Anulldim<n)&&(Aeigvals[Anulldim]<Acutoff)){
       Anulldim++;
       RVector Atemp(n);
       for (int j=0;j<n;j++) Atemp[j]=Aeigvecs(j,Anulldim);
       Anullvecs.push_back(Atemp);}
    xmlnull.putUInt("DimensionNullSpaceAMatrix",Anulldim);
    if (Anulldim<Bnulldim)
       xmlnull.putString("WARNING","Anull smaller than B null");
    Aeigvecs.clear();
    if (Anulldim==n){
       clear();
       xmlout.putString("Warning","Null space is dim of Matrix in HermDiagonalizerWithMetric::setMetric");
       if (xon) throw(std::invalid_argument("Null space is dim of Matrix in HermDiagonalizerWithMetric::setMetric"));
       else return -3;}
    double ovmax=0.0;
    for (int bvec=0;bvec<Bnulldim;bvec++){
       RVector Bnullvec(n);
       int ind=n*bvec;
       for (int j=0;j<n;j++){
          Bnullvec[j]=matb[ind++];}
       LogHelper xmlbvec("MetricNullEigenvectorTest");
       xmlbvec.putUInt("Index",bvec);
       RVector Anullvec;
       multiply(Anullvec,A,Bnullvec); // Anullvec = A * nullvec
       xmlbvec.putReal("MagnitudeOfMatrixElementOfAMatrix",
                   dotProductMagnitude(Bnullvec,Anullvec));
       double ov=0.0;
       for (int k=0;k<Anulldim;k++){
          ov+=dotProductMagnitudeSquared(Anullvecs[k],Bnullvec);}
       xmlbvec.putReal("NormProjectionIntoAMatrixNullSpace",ov);
       xmlnull.put(xmlbvec);
       if (ov>ovmax) ovmax=ov;}
    nullB_in_nullA=(ovmax>0.9);
    if (nullB_in_nullA)
       xmlnull.putString("MetricNullSpace"," SubsetOfMatrixNullSpace ");
    else
       xmlnull.putString("MetricNullSpace"," NOTSubsetOfMatrixNullSpace ");
    xmlout.put(xmlnull);}


 int info;
 char jobz='V';
 char uplo='U';
 int lwork=5*n;
 RVector work(lwork);
 int Bnull=n-n0;
 RVector Btildeinvsqrt(n0);
 RVector ev(n0);
 for (int i=0;i<n0;i++)
    Btildeinvsqrt[i]=1.0/sqrt(Beigvals[i+Bnull]);

     // make the matrix  matg = Btilde^(-1/2) Atilde Btilde^(-1/2)
 matg.resize(n0*n0);
 for (int col=0;col<n0;++col){
    int fcol=col+Bnull;
    for (int row=0;row<=col;++row){
       int frow=row+Bnull;
       double tmp=0.0;
       for (int k=0;k<n;k++){
          double bk=matb[k+n*fcol];
          for (int l=0;l<n;l++)
             tmp+=matb[l+n*frow]*A(l,k)*bk;}
       matg[row+n0*col]=tmp*Btildeinvsqrt[row]*Btildeinvsqrt[col];}}

     // diagonalize, orthonormal eigenvectors in columns of matg,
     // eigenvalues in "ev"
 dsyev_(&jobz,&uplo,&n0,&matg[0],&n0,&ev[0],&work[0],&lwork,&info);

 if (info<0){
    clearMatrix();
    if (xon) throw(std::invalid_argument(" bad arguments in RealSymDiagonalizerWithMetric::setMatrix"));
    return -4;}
 else if (info>0){
    clearMatrix();
    if (xon) throw(std::invalid_argument(" no convergence in RealSymDiagonalizerWithMetric::setMatrix"));
    return -4;}

 xmlout.putString("ReturnCode",getRotateMatrixCode(info));
 if ((info==0)||(info==-5)){
    LogHelper xmleig("GMatrixAllEigenvalues");
    for (uint k=0;k<ev.size();k++)
       xmleig.putReal("Value",ev[k]);
    xmlout.putItem(xmleig);}

 int Aremove=0;
 double cutoff=std::abs(mininvcondnum*ev[n0-1]);
 if (ev[n0-1]<mininvcondnum) cutoff=mininvcondnum;
 if (ev[0]<negeigalarm){
    clearMatrix();
    if (xon) throw(std::invalid_argument(" Matrix not positive semidefinite in RealSymDiagonalizerWithMetric::setMatrix"));
    else return -6;}
 while ((Aremove<n0)&&(ev[Aremove]<cutoff)) Aremove++;
 if (Aremove==n0){
    clearMatrix();
    if (xon) throw(std::invalid_argument("Null space is dim of Matrix in RealSymDiagonalizerWithMetric::setMatrix"));
    else return -3;}

   //  put final retained eigenvalues in "eigvals"
 np=n0-Aremove;
 Geigvals.resize(np);
 for (int l=0;l<np;l++){
    Geigvals[np-l-1]=ev[l+Aremove];}
 xmlout.putInt("GMatrixRank",np);
 LogHelper xmlmat("GMatrixRetainedEigenvalues");
 for (int k=0;k<np;k++)
    xmlmat.putReal("Value",Geigvals[k]);
 xmlout.put(xmlmat);
 Aset=true;

 if (!checkNullSpace) return 0;
 return (nullB_in_nullA)?0:-5;
}


int RealSymDiagonalizerWithMetric::setMatrix(const RealSymmetricMatrix& A,
                                             bool checkNullSpace)
{
 LogHelper xmlout;
 return setMatrix(A,xmlout,checkNullSpace);
}


void RealSymDiagonalizerWithMetric::getEigenvalues(RVector& eigvals)
{
 if (!Aset){
    eigvals.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for RealSymDiagonalizerWithMetric::getEigenvalues"));
    return;}
 eigvals=Geigvals;
}


void RealSymDiagonalizerWithMetric::getEigenvectors(RMatrix& eigvecs)
{
 if (!Aset){
    eigvecs.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for RealSymDiagonalizerWithMetric::getEigenvectors"));
    return;}

 int Bnull=n-n0;
 int Aremove=n0-np;
 RVector Btildeinvsqrt(n0);
 for (int i=0;i<n0;i++)
    Btildeinvsqrt[i]=1.0/sqrt(Beigvals[i+Bnull]);

   //  compute the generalized eigenvectors P0 Btilde^(-1/2) matg

 eigvecs.resize(n,np);
 for (int l=0;l<np;l++)
 for (int i=0;i<n;i++){
    double tmp=0.0; int jj=Bnull;
    for (int j=0;j<n0;j++){
       tmp+=matb[i+n*jj]*Btildeinvsqrt[j]*matg[j+n0*(l+Aremove)];
       jj++;}
    eigvecs(i,np-l-1)=tmp;}
}



void RealSymDiagonalizerWithMetric::getOrthovectors(RMatrix& orthovecs)
{
 if (!Aset){
    orthovecs.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for RealSymDiagonalizerWithMetric::getOrthovectors"));
    return;}

 int Bnull=n-n0;
 int Aremove=n0-np;

   //  compute the orthonormal eigenvectors P0 matg

 orthovecs.resize(n,np);
 for (int l=0;l<np;l++)
 for (int i=0;i<n;i++){
    double tmp=0.0; int jj=Bnull;
    for (int j=0;j<n0;j++){
       tmp+=matb[i+n*jj]*matg[j+n0*(l+Aremove)];
       jj++;}
    orthovecs(i,np-l-1)=tmp;}
}


void RealSymDiagonalizerWithMetric::getZMatrix(RMatrix& Zmat)
{
 if (!Aset){
    Zmat.clear();
    if (xon) throw(std::invalid_argument("Matrix not yet set for RealSymDiagonalizerWithMetric::getZMatrix"));
    return;}

 int Bnull=n-n0;
 int Aremove=n0-np;
 RVector Btildesqrt(n0);
 for (int i=0;i<n0;i++)
    Btildesqrt[i]=sqrt(Beigvals[i+Bnull]);

   //  compute Zmat = P0 Btilde^(1/2) matg

 Zmat.resize(n,np);
 for (int l=0;l<np;l++)
 for (int i=0;i<n;i++){
    double tmp=0.0; int jj=Bnull;
    for (int j=0;j<n0;j++){
       tmp+=matb[i+n*jj]*Btildesqrt[j]*matg[j+n0*(l+Aremove)];
       jj++;}
    Zmat(i,np-l-1)=tmp;}
}


string RealSymDiagonalizerWithMetric::getRotateMetricCode(int info)
{
 return HermDiagonalizerWithMetric::getRotateMetricCode(info);
}

string RealSymDiagonalizerWithMetric::getRotateMatrixCode(int info)
{
 return HermDiagonalizerWithMetric::getRotateMatrixCode(info);
}


// **************************************************************


void analyzeHermCorrelatorMatrix(MCObsHandler *moh,
                  const CorrelatorMatrixInfo& cormat, uint timeval,
                  vector<MCEstimate>& corr_diag_estimates,
                  vector<MCEstimate>& eigenvalues)
{
 try{
 const set<OperatorInfo>& corrops=cormat.getOperators();
 uint nops=cormat.getNumberOfOperators();
 bool herm=cormat.isHermitian();
 if (!herm){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for this case"));}
 bool subvevs=cormat.subtractVEV();
 corr_diag_estimates.resize(nops);
 eigenvalues.resize(nops);

 vector<Vector<double> > jackvals(nops*nops);  //  jackknife samplings
 uint count=0, icol=0;
 set<OperatorInfo>::const_iterator itrow,itcol;
 for (itcol=corrops.begin();itcol!=corrops.end();itcol++,icol++){
    for (itrow=corrops.begin();itrow!=itcol;itrow++){
       CorrelatorAtTimeInfo corrtv(*itrow,*itcol,timeval,true,subvevs);
       MCObsInfo obskey(corrtv,RealPart);
       jackvals[count++]=moh->getJackknifeSamplingValues(obskey);
       obskey.setToImaginaryPart();
       jackvals[count++]=moh->getJackknifeSamplingValues(obskey);
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,timeval,true,false);
       moh->eraseData(MCObsInfo(corrt,RealPart));
       moh->eraseData(MCObsInfo(corrt,ImaginaryPart));}
    CorrelatorAtTimeInfo corrtv(*itcol,*itcol,timeval,true,subvevs);
    MCObsInfo diagkey(corrtv,RealPart);
    jackvals[count++]=moh->getJackknifeSamplingValues(diagkey);
    corr_diag_estimates[icol]=moh->getEstimate(diagkey);
    CorrelatorAtTimeInfo corrt(*itcol,*itcol,timeval,true,false);
    moh->eraseData(MCObsInfo(corrt,RealPart));}
 if (subvevs){
    for (itcol=corrops.begin();itcol!=corrops.end();itcol++){
       moh->eraseData(MCObsInfo(*itcol,RealPart));
       moh->eraseData(MCObsInfo(*itcol,ImaginaryPart));}}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument("Could not calculate jackknife sampling"));}
}



// **************************************************************


// Given a positive definite real symmetric matrix "A",
// this routine constructs its Cholesky decomposition
//                 A = L * transpose(L)
// where L is lower triangular.  Throws an exception if not successful.


void CholeskyDecomposer::getCholesky(const RealSymmetricMatrix& A,
                                     LowerTriangularMatrix<double>& L)
{
 int n=A.size();
 if (n==0){
    L.clear(); return;}
 int info;
 char uplo='L';

    // load A into lower triangle of mata in fortran format
    //    (column major; row index changes fastest)
 vector<double> mata(n*n);
 for (int row=0;row<n;++row)
 for (int col=0;col<=row;++col)
    mata[row+n*col]=A(row,col);

 dpotrf_(&uplo,&n,&mata[0],&n,&info);
 if (info<0){
    L.clear();
    throw(std::invalid_argument(" bad arguments in cholesky"));}
 else if (info>0){
    L.clear();
    throw(std::invalid_argument(" matrix not positive definite in cholesky"));}

 L.resize(n);
 for (int row=0;row<n;++row)
 for (int col=0;col<=row;++col)
    L(row,col)=mata[row+n*col];
}

// *************************************************************

   // Given a positive definite real symmetric matrix "A",
   // this routine constructs the Cholesky decomposition of the
   // inverse of A:
   //                 A^(-1) = transpose(L) * L
   // where L is lower triangular. Throws an exception if not successful.


void CholeskyDecomposer::getCholeskyOfInverse(const RealSymmetricMatrix& A,
                                      LowerTriangularMatrix<double>& L)
{
 try{
    getCholesky(A,L);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Failure in cholesky_of_inverse: ")
              +string(errmsg.what())));}
 int n=A.size();
 for (int i=0;i<n;i++){
    L(i,i)=1.0/L(i,i);
    for (int j=i+1;j<n;j++){
       double sum=0.0;
       for (int k=i;k<j;k++) sum-=L(j,k)*L(k,i);
       L(j,i)=sum/L(j,j);}}
}


// ************************************************************



template <>
double VectorPinner<double>::dot_prod(const double *avec,
                                      const double *bvec) const
{
 double res=0.0;
 for (uint i=0;i<m_veclength;++i)
    res+=avec[i]*bvec[i];
 return res;
}

template <>
std::complex<double> VectorPinner<std::complex<double> >::dot_prod(
                                      const std::complex<double> *avec,
                                      const std::complex<double> *bvec) const
{
 std::complex<double> res(0.0,0.0);
 for (uint i=0;i<m_veclength;++i)
    res+=conjugate(avec[i])*bvec[i];
 return res;
}


// *****************************************************************************

    //   Rescales the matrix "cormat" using the diagonal elements
    //   of the matrix "mat_scales" according to
    //
    //     cormat(i,j) / sqrt( |mat_scales(i,i)|*|mat_scales(j,j)| )


void doRescaleByDiagonals(ComplexHermitianMatrix& cormat,
                          const ComplexHermitianMatrix& mat_scales)
{
 int n=cormat.size();
 if (int(mat_scales.size())!=n)
    throw(std::invalid_argument("Size mismatch in doRescaleByDiagonals"));
 RVector scales(n);
 for (int k=0;k<n;k++)
    scales[k]=1.0/sqrt(std::abs(mat_scales(k,k).real()));
 for (int row=0;row<n;row++)
 for (int col=row;col<n;col++)
    cormat.put(row,col,cormat(row,col)*scales[row]*scales[col]);
}


void doRescaleByDiagonals(RealSymmetricMatrix& cormat,
                          const RealSymmetricMatrix& mat_scales)
{
 int n=cormat.size();
 if (int(mat_scales.size())!=n)
    throw(std::invalid_argument("Size mismatch in doRescaleByDiagonals"));
 RVector scales(n);
 for (int k=0;k<n;k++)
    scales[k]=1.0/sqrt(std::abs(mat_scales(k,k)));
 for (int row=0;row<n;row++)
 for (int col=row;col<n;col++)
    cormat(row,col)*=scales[row]*scales[col];
}


// *****************************************************************************

    //   Rescales the transformation matrix "R" using the
    //   diagonal elements of the matrix "mat_scales" according to
    //
    //     R(i,j) / sqrt( |mat_scales(i,i)| )


void doRescaleTransformation(CMatrix& R,
              const ComplexHermitianMatrix& mat_scales)
{
 int n=R.size(0);
 if (int(mat_scales.size())!=n)
    throw(std::invalid_argument("Size mismatch in doRescaleTransformation"));
 int np=R.size(1);
 RVector scales(n);
 for (int k=0;k<n;k++)
    scales[k]=1.0/sqrt(std::abs(mat_scales(k,k).real()));
 for (int row=0;row<n;row++)
 for (int col=0;col<np;col++)
    R.put(row,col,R(row,col)*scales[row]);
}

void doRescaleTransformation(RMatrix& R,
              const RealSymmetricMatrix& mat_scales)
{
 int n=R.size(0);
 if (int(mat_scales.size())!=n)
    throw(std::invalid_argument("Size mismatch in doRescaleTransformation"));
 int np=R.size(1);
 RVector scales(n);
 for (int k=0;k<n;k++)
    scales[k]=1.0/sqrt(std::abs(mat_scales(k,k)));
 for (int row=0;row<n;row++)
 for (int col=0;col<np;col++)
    R(row,col)*=scales[row];
}


// *************************************************************

     //  Takes Hermitian matrix "A" and replaces it with the "rotated"
     //  matrix   R^dagger A  R.  There is also a version that just
     //  evaluates the diagonal elements of the rotated matrix.

void doMatrixRotation(ComplexHermitianMatrix& A, const CMatrix& R)
{
 int n=R.size(0);
 int np=R.size(1);
 if (int(A.size())!=n)
    throw(std::invalid_argument("Matrix multiply size mismatch"));

 CMatrix AR(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<n;k++)
       tmp+=A(i,k)*R(k,j);
    AR(i,j)=tmp;}

 A.resize(np);
 for (int i=0;i<np;i++)
 for (int j=i;j<np;j++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<n;k++)
       tmp+=conjugate(R(k,i))*AR(k,j);
    if (i!=j) A.put(i,j,tmp);
    else A.put(i,j,complex<double>(tmp.real(),0.0));}
}


void doMatrixRotation(const ComplexHermitianMatrix& A, const CMatrix& R,
                      RVector& Ardiag)
{
 int n=R.size(0);
 int np=R.size(1);
 if (int(A.size())!=n)
    throw(std::invalid_argument("Matrix multiply size mismatch"));

 CMatrix AR(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<n;k++)
       tmp+=A(i,k)*R(k,j);
    AR(i,j)=tmp;}

 Ardiag.resize(np);
 for (int i=0;i<np;i++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=R(k,i).real()*AR(k,i).real()
           +R(k,i).imag()*AR(k,i).imag();
    Ardiag[i]=tmp;}
}


void doMatrixRotation(RealSymmetricMatrix& A, const RMatrix& R)
{
 int n=R.size(0);
 int np=R.size(1);
 if (int(A.size())!=n)
    throw(std::invalid_argument("Matrix multiply size mismatch"));

 RMatrix AR(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=A(i,k)*R(k,j);
    AR(i,j)=tmp;}

 A.resize(np);
 for (int i=0;i<np;i++)
 for (int j=i;j<np;j++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=R(k,i)*AR(k,j);
    A(i,j)=tmp;}
}


void doMatrixRotation(const RealSymmetricMatrix& A, const RMatrix& R,
                      RVector& Ardiag)
{
 int n=R.size(0);
 int np=R.size(1);
 if (int(A.size())!=n)
    throw(std::invalid_argument("Matrix multiply size mismatch"));

 RMatrix AR(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=A(i,k)*R(k,j);
    AR(i,j)=tmp;}

 Ardiag.resize(np);
 for (int i=0;i<np;i++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=R(k,i)*AR(k,i);
    Ardiag[i]=tmp;}
}

  // version with non Hermitian but square matrix A

void doMatrixRotation(CMatrix& A, const CMatrix& R)
{
 int n=R.size(0);
 int np=R.size(1);
 if ((int(A.size(0))!=n)||(int(A.size(1))!=n))
    throw(std::invalid_argument("Matrix multiply size mismatch"));

 CMatrix AR(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<n;k++)
       tmp+=A(i,k)*R(k,j);
    AR(i,j)=tmp;}

 A.resize(np,np);
 for (int i=0;i<np;i++)
 for (int j=0;j<np;j++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<n;k++)
       tmp+=conjugate(R(k,i))*AR(k,j);
    if (i!=j) A.put(i,j,tmp);
    else A.put(i,j,complex<double>(tmp.real(),0.0));}
}


void doMatrixRotation(RMatrix& A, const RMatrix& R)
{
 int n=R.size(0);
 int np=R.size(1);
 if ((int(A.size(0))!=n)||(int(A.size(1))!=n))
    throw(std::invalid_argument("Matrix multiply size mismatch"));

 RMatrix AR(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=A(i,k)*R(k,j);
    AR(i,j)=tmp;}

 A.resize(np,np);
 for (int i=0;i<np;i++)
 for (int j=0;j<np;j++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=R(k,i)*AR(k,j);
    A(i,j)=tmp;}
}

// ********************************************************************

     //  Takes vector "V" and replaces it with the rotated
     //  matrix   R^dagger V.


void doVectorRotation(CVector& V, const CMatrix& R)
{
 int n=R.size(0);
 int np=R.size(1);
 if (int(V.size())!=n)
    throw(std::invalid_argument("Matrix-vector multiply size mismatch"));
 CVector res(np);
 for (int i=0;i<np;i++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<n;k++)
       tmp+=conjugate(R(k,i))*V[k];
    res[i]=tmp;}
 V=res;
}

void doVectorRotation(RVector& V, const RMatrix& R)
{
 int n=R.size(0);
 int np=R.size(1);
 if (int(V.size())!=n)
    throw(std::invalid_argument("Matrix-vector multiply size mismatch"));
 RVector res(np);
 for (int i=0;i<np;i++){
    double tmp=0.0;
    for (int k=0;k<n;k++)
       tmp+=R(k,i)*V[k];
    res[i]=tmp;}
 V=res;
}


// ********************************************************************

     //  Replaces V by R*V


void doMatrixMultiply(const CMatrix& R, CMatrix& V)
{
 int n=R.size(0);
 int np=V.size(1);
 int nk=R.size(1);
 if (int(V.size(0))!=nk)
    throw(std::invalid_argument("Matrix multiply size mismatch"));
 CMatrix RV(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    complex<double> tmp(0.0,0.0);
    for (int k=0;k<nk;k++)
       tmp+=R(i,k)*V(k,j);
    RV(i,j)=tmp;}
 V=RV;
}


void doMatrixMultiply(const RMatrix& R, RMatrix& V)
{
 int n=R.size(0);
 int np=V.size(1);
 int nk=R.size(1);
 if (int(V.size(0))!=nk)
    throw(std::invalid_argument("Matrix multiply size mismatch"));
 RMatrix RV(n,np);
 for (int i=0;i<n;i++)
 for (int j=0;j<np;j++){
    double tmp=0.0;
    for (int k=0;k<nk;k++)
       tmp+=R(i,k)*V(k,j);
    RV(i,j)=tmp;}
 V=RV;
}


// ********************************************************************

     //   outvec = inmat * invec

void multiply(CVector& outvec, const ComplexHermitianMatrix& inmat, const CVector& invec)
{
 uint nn=inmat.size();
 if (nn!=invec.size()) throw(std::invalid_argument("bad matrix-vector multiplying sizes"));
 outvec.resize(nn);
 for (uint row=0;row<nn;row++){
    complex<double> z(0.0,0.0);
    for (uint k=0;k<nn;k++)
       z+=inmat(row,k)*invec[k];
    outvec[row]=z;}
}



void multiply(RVector& outvec, const RealSymmetricMatrix& inmat, const RVector& invec)
{
 uint nn=inmat.size();
 if (nn!=invec.size()) throw(std::invalid_argument("bad matrix-vector multiplying sizes"));
 outvec.resize(nn);
 for (uint row=0;row<nn;row++){
    double z=0.0;
    for (uint k=0;k<nn;k++)
       z+=inmat(row,k)*invec[k];
    outvec[row]=z;}
}

     //   inner product of vectors

complex<double> dotProduct(const CVector& lvec, const CVector& rvec)
{
 uint nv=rvec.size();
 if (nv!=lvec.size()) throw(std::invalid_argument("Mismatch in vector sizes in dotProduct"));
 complex<double> res(0.0,0.0);
 for (uint k=0;k<nv;k++)
    res+=conjugate(lvec[k])*rvec[k];
 return res;
}


double dotProduct(const RVector& lvec, const RVector& rvec)
{
 uint nv=rvec.size();
 if (nv!=lvec.size()) throw(std::invalid_argument("Mismatch in vector sizes in dotProduct"));
 double res=0.0;
 for (uint k=0;k<nv;k++)
    res+=lvec[k]*rvec[k];
 return res;
}

     //   magnitude of inner product of vectors

double dotProductMagnitude(const CVector& lvec, const CVector& rvec)
{
 return std::abs(dotProduct(lvec,rvec));
}

double dotProductMagnitude(const RVector& lvec, const RVector& rvec)
{
 return std::abs(dotProduct(lvec,rvec));
}


     //   magnitude squared of inner product of vectors

double dotProductMagnitudeSquared(const CVector& lvec, const CVector& rvec)
{
 return std::norm(dotProduct(lvec,rvec));
}

double dotProductMagnitudeSquared(const RVector& lvec, const RVector& rvec)
{
 return sqr(dotProduct(lvec,rvec));
}


// ********************************************************************


void array_to_matrix(const Array<complex<double> >& in, CMatrix& out)
{
 if (in.numDimensions()!=2)
    throw(std::invalid_argument("Invalid array to matrix conversion"));
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void array_to_matrix(const Array<complex<float> >& in, CMatrix& out)
{
 if (in.numDimensions()!=2)
    throw(std::invalid_argument("Invalid array to matrix conversion"));
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void array_to_matrix(const Array<double>& in, RMatrix& out)
{
 if (in.numDimensions()!=2)
    throw(std::invalid_argument("Invalid array to matrix conversion"));
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void array_to_matrix(const Array<float>& in, RMatrix& out)
{
 if (in.numDimensions()!=2)
    throw(std::invalid_argument("Invalid array to matrix conversion"));
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void array_to_matrix(const Array<double>& in, CMatrix& out)
{
 if ((in.numDimensions()!=2)||((in.size(1)%2)!=0))
    throw(std::invalid_argument("Invalid array to matrix conversion"));
 uint nrow=in.size(0);
 uint ncol=in.size(1)/2;
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=complex<double>(in(row,col),in(row,col+ncol));
}

void array_to_matrix(const Array<float>& in, CMatrix& out)
{
 if ((in.numDimensions()!=2)||((in.size(1)%2)!=0))
    throw(std::invalid_argument("Invalid array to matrix conversion"));
 uint nrow=in.size(0);
 uint ncol=in.size(1)/2;
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=complex<double>(in(row,col),in(row,col+ncol));
}

void matrix_to_array(const CMatrix& in, Array<complex<double> >& out)
{
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void matrix_to_array(const CMatrix& in, Array<complex<float> >& out)
{
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void matrix_to_array(const RMatrix& in, Array<double>& out)
{
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}

void matrix_to_array(const RMatrix& in, Array<float>& out)
{
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++)
    out(row,col)=in(row,col);
}


void matrix_to_array(const CMatrix& in, Array<double>& out)
{
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,2*ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++){
    out(row,col)=in(row,col).real();
    out(row,col+ncol)=in(row,col).imag();}
}

void matrix_to_array(const CMatrix& in, Array<float>& out)
{
 uint nrow=in.size(0);
 uint ncol=in.size(1);
 out.resize(nrow,2*ncol);
 for (uint row=0;row<nrow;row++)
 for (uint col=0;col<ncol;col++){
    out(row,col)=in(row,col).real();
    out(row,col+ncol)=in(row,col).imag();}
}

void array_to_vector(const Array<double>& in, std::vector<double>& out)
{
 uint n=in.size();
 out.resize(n);
 for (uint k=0;k<n;k++)
    out[k]=in[k];
}

void array_to_RVector(const Array<double>& in, RVector& out)
{
 uint n=in.size();
 out.resize(n);
 for (uint k=0;k<n;k++)
    out[k]=in[k];
}

void vector_to_array(const std::vector<double>& in, Array<double>& out)
{
 uint n=in.size();
 out.resize(n);
 for (uint k=0;k<n;k++)
    out[k]=in[k];
}

void RVector_to_array(const RVector& in, Array<double>& out)
{
 uint n=in.size();
 out.resize(n);
 for (uint k=0;k<n;k++)
    out[k]=in[k];
}


   //  returns a vector of integer times in ascending order from
   //  tmin to tmax, excluding any values contained in "texclude"

vector<uint> form_tvalues(uint tmin, uint tmax,
                          const vector<int>& texclude)
{
 set<uint> tvals;  // values will automatically be sorted in set
 for (uint tt=tmin;tt<=tmax;tt++){
    tvals.insert(tt);}
 for (uint k=0;k<texclude.size();k++){
    tvals.erase(texclude[k]);}
 vector<uint> result(tvals.begin(),tvals.end());
 //if (result.size()<4) throw(std::invalid_argument("Time range too small"));
    // should be sorted already, but double check that it is sorted
 for (uint k=1;k<result.size();k++){
    if (result[k-1]>=result[k])
       throw(std::invalid_argument("Not sorted but should be!"));}
 return result;
}


// ********************************************************************

void doCopyByBins(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 const Vector<double>& inbins=moh.getBins(obs_in);
 moh.putBins(obs_out,inbins);
}


void doCopyBySamplings(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double val=moh.getCurrentSamplingValue(obs_in);
    moh.putCurrentSamplingValue(obs_out,val);}
}


void doSquareByBins(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 const Vector<double>& inbins=moh.getBins(obs_in);
 int nbins=inbins.size();
 Vector<double> sqvalues(nbins);
 for (int bin=0;bin<nbins;bin++)
    sqvalues[bin]=inbins[bin]*inbins[bin];
 moh.putBins(obs_out,sqvalues);
}


void doSquareBySamplings(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double val=moh.getCurrentSamplingValue(obs_in);
    moh.putCurrentSamplingValue(obs_out,val*val);}
}


void doSquareRootByBins(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 const Vector<double>& inbins=moh.getBins(obs_in);
 int nbins=inbins.size();
 Vector<double> sqvalues(nbins);
 for (int bin=0;bin<nbins;bin++)
   sqvalues[bin]=sqrt(inbins[bin]);
 moh.putBins(obs_out,sqvalues);
}


void doSquareRootBySamplings(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double val=moh.getCurrentSamplingValue(obs_in);
    moh.putCurrentSamplingValue(obs_out,sqrt(val));}
}


void doLogByBins(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 const Vector<double>& inbins=moh.getBins(obs_in);
 int nbins=inbins.size();
 Vector<double> logvalues(nbins);
 for (int bin=0;bin<nbins;bin++)
   logvalues[bin]=log(inbins[bin]);
 moh.putBins(obs_out,logvalues);
}


void doLogBySamplings(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double val=moh.getCurrentSamplingValue(obs_in);
    moh.putCurrentSamplingValue(obs_out,log(val));}
}

void doExpByBins(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 const Vector<double>& inbins=moh.getBins(obs_in);
 int nbins=inbins.size();
 Vector<double> expvalues(nbins);
 for (int bin=0;bin<nbins;bin++)
   expvalues[bin]=exp(inbins[bin]);
 moh.putBins(obs_out,expvalues);
}


void doExpBySamplings(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double val=moh.getCurrentSamplingValue(obs_in);
    moh.putCurrentSamplingValue(obs_out,exp(val));}
}


void doRatioByBins(MCObsHandler& moh, const MCObsInfo& obs_numer, const MCObsInfo& obs_denom,
                   const MCObsInfo& obs_ratio)
{
 const Vector<double>& numerbins=moh.getBins(obs_numer);
 const Vector<double>& denombins=moh.getBins(obs_denom);
 int nbins=numerbins.size();
 Vector<double> ratiovalues(nbins);
 for (int bin=0;bin<nbins;bin++)
    ratiovalues[bin]=numerbins[bin]/denombins[bin];
 moh.putBins(obs_ratio,ratiovalues);
}


void doRatioBySamplings(MCObsHandler& moh, const MCObsInfo& obs_numer, const MCObsInfo& obs_denom,
                        const MCObsInfo& obs_ratio)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double ratiovalue=moh.getCurrentSamplingValue(obs_numer)
                     /moh.getCurrentSamplingValue(obs_denom);
    moh.putCurrentSamplingValue(obs_ratio,ratiovalue);}
}

void doGMOByBins(MCObsHandler& moh, const CorrelatorInfo& L, const CorrelatorInfo& S,
                 const CorrelatorInfo& N, const CorrelatorInfo& X,
                   const CorrelatorInfo& GMO)
{
    CorrelatorAtTimeInfo corrL(L,0);
    CorrelatorAtTimeInfo corrS(S,0);
    CorrelatorAtTimeInfo corrN(N,0);
    CorrelatorAtTimeInfo corrX(X,0);
    CorrelatorAtTimeInfo corrGMO(GMO,0);
    double one_third = 1.0/3.0;
    double two_thirds = 2.0/3.0;
     for (uint tval=0;tval<moh.getLatticeTimeExtent();tval++){
        corrL.resetTimeSeparation(tval);
        corrS.resetTimeSeparation(tval);
        corrN.resetTimeSeparation(tval);
        corrX.resetTimeSeparation(tval);
        corrGMO.resetTimeSeparation(tval);
        MCObsInfo obsL(corrL);
        MCObsInfo obsS(corrS);
        MCObsInfo obsN(corrN);
        MCObsInfo obsX(corrX);
        MCObsInfo obsGMO(corrGMO);
        if( moh.queryBins(obsL) ){
            const Vector<double>& Lbins=moh.getBins(obsL);
            const Vector<double>& Sbins=moh.getBins(obsS);
            const Vector<double>& Nbins=moh.getBins(obsN);
            const Vector<double>& Xbins=moh.getBins(obsX);
            int nbins=Lbins.size();
            Vector<double> GMObins(nbins);
            for (int bin=0;bin<nbins;bin++) 
                GMObins[bin]=Lbins[bin]*pow(Sbins[bin],one_third)/( pow(Nbins[bin],two_thirds) * pow(Xbins[bin],two_thirds) );
            moh.putBins(obsGMO,GMObins);
        }
     }
}


void doGMOBySamplings(MCObsHandler& moh, const CorrelatorInfo& L, const CorrelatorInfo& S,
                 const CorrelatorInfo& N, const CorrelatorInfo& X,
                        const CorrelatorInfo& GMO)
{
    CorrelatorAtTimeInfo corrL(L,0);
    CorrelatorAtTimeInfo corrS(S,0);
    CorrelatorAtTimeInfo corrN(N,0);
    CorrelatorAtTimeInfo corrX(X,0);
    CorrelatorAtTimeInfo corrGMO(GMO,0);
    double one_third = 1.0/3.0;
    double two_thirds = 2.0/3.0;
     for (uint tval=0;tval<moh.getLatticeTimeExtent();tval++){
        corrL.resetTimeSeparation(tval);
        corrS.resetTimeSeparation(tval);
        corrN.resetTimeSeparation(tval);
        corrX.resetTimeSeparation(tval);
        corrGMO.resetTimeSeparation(tval);
        MCObsInfo obsL(corrL);
        MCObsInfo obsS(corrS);
        MCObsInfo obsN(corrN);
        MCObsInfo obsX(corrX);
        MCObsInfo obsGMO(corrGMO);
        if( moh.queryBins(obsL) ){
         for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
            double gmovalue=moh.getCurrentSamplingValue(obsL)*pow(moh.getCurrentSamplingValue(obsS),one_third)
                /( pow(moh.getCurrentSamplingValue(obsN),two_thirds) * pow(moh.getCurrentSamplingValue(obsX),two_thirds) );
            moh.putCurrentSamplingValue(obsGMO,gmovalue);
         }
        }
     }
}

void doGMOBySamplings(MCObsHandler& moh, const MCObsInfo& L, const MCObsInfo& S,
                 const MCObsInfo& N, const MCObsInfo& X,
                        const MCObsInfo& GMO)
{
    double one_third = 1.0/3.0;
    double two_thirds = 2.0/3.0;
     for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
        double gmovalue=moh.getCurrentSamplingValue(L)+moh.getCurrentSamplingValue(S)*one_third
            -moh.getCurrentSamplingValue(N)*two_thirds - moh.getCurrentSamplingValue(X)*two_thirds;
        moh.putCurrentSamplingValue(GMO,gmovalue);
     }
}

void doLinearSuperpositionByBins(MCObsHandler& moh, std::vector<MCObsInfo>& suminfos,
                   std::vector<double>& sumcoefs, const MCObsInfo& obs_superposition)
{
 int nsummands=suminfos.size();
 vector<const Vector<double>* > bins(nsummands);
 for (int k=0;k<nsummands;k++)
    bins[k]=&(moh.getBins(suminfos[k]));
 int nbins=bins[0]->size();
 Vector<double> result(nbins);
 for (int bin=0;bin<nbins;bin++){
    double temp=0.0;
    for (int k=0;k<nsummands;k++)
       temp+=sumcoefs[k]*(*(bins[k]))[bin];
    result[bin]=temp;}
 moh.putBins(obs_superposition,result);
}


void doLinearSuperpositionBySamplings(MCObsHandler& moh, std::vector<MCObsInfo>& suminfos,
                   std::vector<double>& sumcoefs, const MCObsInfo& obs_superposition)
{
 int nsummands=suminfos.size();
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double result=0.0;
    for (int k=0;k<nsummands;k++)
       result+=sumcoefs[k]*moh.getCurrentSamplingValue(suminfos[k]);
    moh.putCurrentSamplingValue(obs_superposition,result);}
}


void doAnisoDispersionBySamplings(MCObsHandler& moh, const MCObsInfo& anisotropy_key,
                             const MCObsInfo& restmasssquared_key, double psqfactor,
                             const MCObsInfo& Esqinfo)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double m0sq=moh.getCurrentSamplingValue(restmasssquared_key);
    double xi=moh.getCurrentSamplingValue(anisotropy_key);
    double Esq=m0sq+psqfactor/(xi*xi);
    moh.putCurrentSamplingValue(Esqinfo,Esq);}
}

void doCoeffDispersionBySamplings(MCObsHandler& moh, const MCObsInfo& coeff_key,
                             const MCObsInfo& restmasssquared_key, double psqfactor,
                             const MCObsInfo& Esqinfo)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double m0sq=moh.getCurrentSamplingValue(restmasssquared_key);
    double c=moh.getCurrentSamplingValue(coeff_key);
    double Esq=m0sq+c*psqfactor;
    moh.putCurrentSamplingValue(Esqinfo,Esq);}
}


void doBoostBySamplings(MCObsHandler& moh, const MCObsInfo& restmass_key,
			const MCObsInfo& anisotropy_key, double psqfactor,
			const MCObsInfo& Eboosted)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double m0=moh.getCurrentSamplingValue(restmass_key);
    double xi=moh.getCurrentSamplingValue(anisotropy_key);
    double Esq=m0*m0+psqfactor/(xi*xi);
    moh.putCurrentSamplingValue(Eboosted,sqrt(Esq));}
}


void doBoostBySamplings(MCObsHandler& moh, const MCObsInfo& restmass_key,
			double psqfactor, const MCObsInfo& Eboosted)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double m0=moh.getCurrentSamplingValue(restmass_key);
    double Esq=m0*m0+psqfactor;
    moh.putCurrentSamplingValue(Eboosted,sqrt(Esq));}
}


void doCorrelatorMatrixTimeDifferencesByBins(MCObsHandler& moh,
                const list<OperatorInfo>& origops,
                const list<OperatorInfo>& newops, bool herm, uint tmin, uint tmax,
                set<MCObsInfo>& obskeys, bool erase_orig)
{
 list<OperatorInfo>::const_iterator oldrow,oldcol,newrow,newcol;
 newcol=newops.begin();
 obskeys.clear();
 for (oldcol=origops.begin();oldcol!=origops.end();oldcol++,newcol++){
    newrow=(herm?newcol:newops.begin());
    for (oldrow=(herm?oldcol:origops.begin());oldrow!=origops.end();oldrow++,newrow++){
       CorrelatorAtTimeInfo origcorr(*oldrow,*oldcol,0,herm,false);
       CorrelatorAtTimeInfo newcorr(*newrow,*newcol,0,herm,false);
       for (uint t=tmin;t<tmax;t++){
          origcorr.resetTimeSeparation(t);
          const RVector& bins1=moh.getBins(MCObsInfo(origcorr));
          origcorr.resetTimeSeparation(t+1);
          const RVector& bins2=moh.getBins(MCObsInfo(origcorr));
          RVector newbins(bins1);
          newbins-=bins2;
          newcorr.resetTimeSeparation(t);
          MCObsInfo newkey(newcorr);
          moh.putBins(newkey,newbins);
          obskeys.insert(newkey);
          origcorr.resetTimeSeparation(t);
          if (erase_orig) moh.eraseData(MCObsInfo(origcorr));
#ifdef COMPLEXNUMBERS
          origcorr.resetTimeSeparation(t);
          const RVector& ibins1=moh.getBins(MCObsInfo(origcorr,ImaginaryPart));
          origcorr.resetTimeSeparation(t+1);
          const RVector& ibins2=moh.getBins(MCObsInfo(origcorr,ImaginaryPart));
          newbins=ibins1;
          newbins-=ibins2;
          newkey.setToImaginaryPart();
          moh.putBins(newkey,newbins);
          obskeys.insert(newkey);
          origcorr.resetTimeSeparation(t);
          if (erase_orig) moh.eraseData(MCObsInfo(origcorr,ImaginaryPart));
#endif
          }}}
}


void doCorrelatorMatrixTimeDifferencesBySamplings(MCObsHandler& moh,
                const list<OperatorInfo>& origops,
                const list<OperatorInfo>& newops, bool herm, uint tmin, uint tmax,
                set<MCObsInfo>& obskeys, bool erase_orig)
{
 list<OperatorInfo>::const_iterator oldrow,oldcol,newrow,newcol;
 newcol=newops.begin();
 obskeys.clear();
 for (oldcol=origops.begin();oldcol!=origops.end();oldcol++,newcol++){
    newrow=(herm?newcol:newops.begin());
    for (oldrow=(herm?oldcol:origops.begin());oldrow!=origops.end();oldrow++,newrow++){
       CorrelatorAtTimeInfo origcorr(*oldrow,*oldcol,0,herm,false);
       CorrelatorAtTimeInfo origcorrp(*oldrow,*oldcol,0,herm,false);
       CorrelatorAtTimeInfo newcorr(*newrow,*newcol,0,herm,false);
       for (uint t=tmin;t<tmax;t++){
          origcorr.resetTimeSeparation(t);
          origcorrp.resetTimeSeparation(t+1);
          newcorr.resetTimeSeparation(t);
          MCObsInfo okey(origcorr);
          MCObsInfo okeyp(origcorrp);
          MCObsInfo newkey(newcorr);
          for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
             double res1=moh.getCurrentSamplingValue(okey);
             double res2=moh.getCurrentSamplingValue(okeyp);
             moh.putCurrentSamplingValue(newkey,res1-res2);}
          obskeys.insert(newkey);
          if (erase_orig){
             moh.eraseSamplings(okey);
             if (t==(tmax-1)) moh.eraseSamplings(okeyp);}
#ifdef COMPLEXNUMBERS
          okey.setToImaginaryPart();
          okeyp.setToImaginaryPart();
          newkey.setToImaginaryPart();
          for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
             double res1=moh.getCurrentSamplingValue(okey);
             double res2=moh.getCurrentSamplingValue(okeyp);
             moh.putCurrentSamplingValue(newkey,res1-res2);}
          obskeys.insert(newkey);
          if (erase_orig){
             moh.eraseSamplings(okey);
             if (t==(tmax-1)) moh.eraseSamplings(okeyp);}
#endif
          }}}
}


void doCorrelatorMatrixSuperpositionByBins(MCObsHandler& moh,
             const list<vector<pair<OperatorInfo,double> > >& superposition,
             const vector<OperatorInfo>& resultops, bool herm,
             uint tmin, uint tmax, set<MCObsInfo>& obskeys, bool erase_orig,
             bool ignore_missing)
{
 obskeys.clear();
 uint nops=resultops.size();
 uint nterms=superposition.size();
 for (list<vector<pair<OperatorInfo,double> > >::const_iterator
      st=superposition.begin();st!=superposition.end();st++){
   if (st->size()!=nops) throw(std::invalid_argument("Invalid superposition in doCorrelatorMatrixSuperpositionByBins"));}
 for (uint row=0;row<nops;row++){
    for (uint col=(herm?row:0);col<nops;col++){
       CorrelatorAtTimeInfo resultcorr(resultops[row],resultops[col],0,herm,false);
       vector<CorrelatorAtTimeInfo> corrterms; 
       vector<double> wts;
       for (list<vector<pair<OperatorInfo,double> > >::const_iterator
            st=superposition.begin();st!=superposition.end();st++){
          corrterms.push_back(CorrelatorAtTimeInfo((*st)[row].first,(*st)[col].first,0,herm,false));
          wts.push_back((*st)[row].second*(*st)[col].second);}
       for (uint t=tmin;t<=tmax;t++){
          //cout << "row = "<<row<<" col = "<<col<<" time sep = "<<t<<endl;
          resultcorr.resetTimeSeparation(t);
          MCObsInfo newkey(resultcorr);
          try{
          for (uint k=0;k<nterms;k++){
             corrterms[k].resetTimeSeparation(t);}
          const RVector& bins1=moh.getBins(MCObsInfo(corrterms[0]));
          RVector newbins(bins1);
          newbins*=wts[0];
          for (uint k=1;k<nterms;k++){
             const RVector& binsk=moh.getBins(MCObsInfo(corrterms[k]));
             RVector addbins(binsk);
             addbins*=wts[k];
             newbins+=addbins;}
          moh.putBins(newkey,newbins);
          obskeys.insert(newkey);
          if (erase_orig)
             for (uint k=0;k<nterms;k++){
                moh.eraseData(MCObsInfo(corrterms[k]));}}  // erase original data
          catch(const std::exception& xp){
             if (!ignore_missing)
                throw(std::runtime_error(string("Failure in doCorrelatorMatrixSuperpositionByBins: ")+xp.what()));}
#ifdef COMPLEXNUMBERS
          try{
          const RVector& ibins1=moh.getBins(MCObsInfo(corrterms[0],ImaginaryPart));
          RVector newbins(ibins1);
          newbins*=wts[0];
          for (uint k=1;k<nterms;k++){
             const RVector& ibinsk=moh.getBins(MCObsInfo(corrterms[k],ImaginaryPart));
             RVector addbins(ibinsk);
             addbins*=wts[k];
             newbins+=addbins;}
          newkey.setToImaginaryPart();
          moh.putBins(newkey,newbins);
          obskeys.insert(newkey);
          if (erase_orig)
             for (uint k=0;k<nterms;k++){
                moh.eraseData(MCObsInfo(corrterms[k],ImaginaryPart));}}  // erase original data
          catch(const std::exception& xp){
             if (!ignore_missing)
                throw(std::runtime_error(string("Failure in doCorrelatorMatrixSuperpositionByBins: ")+xp.what()));}
#endif
          }}}
}


void doCorrelatorMatrixSuperpositionBySamplings(MCObsHandler& moh,
             const list<vector<pair<OperatorInfo,double> > >& superposition,
             const vector<OperatorInfo>& resultops, bool herm,
             uint tmin, uint tmax, set<MCObsInfo>& obskeys, bool erase_orig,
             bool ignore_missing)
{
 obskeys.clear();
 uint nops=resultops.size();
 uint nterms=superposition.size();
 for (list<vector<pair<OperatorInfo,double> > >::const_iterator
      st=superposition.begin();st!=superposition.end();st++){
   if (st->size()!=nops) throw(std::invalid_argument("Invalid superposition in doCorrelatorMatrixSuperpositionByBins"));}
 for (uint row=0;row<nops;row++){
    for (uint col=(herm?row:0);col<nops;col++){
       CorrelatorAtTimeInfo resultcorr(resultops[row],resultops[col],0,herm,false);
       vector<CorrelatorAtTimeInfo> corrterms; 
       vector<double> wts;
       for (list<vector<pair<OperatorInfo,double> > >::const_iterator
            st=superposition.begin();st!=superposition.end();st++){
          corrterms.push_back(CorrelatorAtTimeInfo((*st)[row].first,(*st)[col].first,0,herm,false));
          wts.push_back((*st)[row].second*(*st)[col].second);}
       for (uint t=tmin;t<=tmax;t++){
          //cout << "row = "<<row<<" col = "<<col<<" time sep = "<<t<<endl;
          resultcorr.resetTimeSeparation(t);
          MCObsInfo newkey(resultcorr);
          try{
          for (uint k=0;k<nterms;k++){
             corrterms[k].resetTimeSeparation(t);}
          for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
             double res1=moh.getCurrentSamplingValue(MCObsInfo(corrterms[0]));
             double newres=res1*wts[0];
             for (uint k=1;k<nterms;k++){
                double resk=moh.getCurrentSamplingValue(MCObsInfo(corrterms[k]));
                newres+=resk*wts[k];}
             moh.putCurrentSamplingValue(newkey,newres);}
          obskeys.insert(newkey);
          if (erase_orig)
             for (uint k=0;k<nterms;k++){
                moh.eraseSamplings(MCObsInfo(corrterms[k]));}}  // erase original data
          catch(const std::exception& xp){
             if (!ignore_missing)
                throw(std::runtime_error(string("Failure in doCorrelatorMatrixSuperpositionByBins: ")+xp.what()));}
#ifdef COMPLEXNUMBERS
          try{
          for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
             double ires1=moh.getCurrentSamplingValue(MCObsInfo(corrterms[0],ImaginaryPart));
             double newres=ires1*wts[0];
             for (uint k=1;k<nterms;k++){
                double iresk=moh.getCurrentSamplingValue(MCObsInfo(corrterms[k],ImaginaryPart));
                newres+=iresk*wts[k];}
             newkey.setToImaginaryPart();
             moh.putCurrentSamplingValue(newkey,newres);}
          obskeys.insert(newkey);
          if (erase_orig)
             for (uint k=0;k<nterms;k++){
                moh.eraseSamplings(MCObsInfo(corrterms[k],ImaginaryPart));}}  // erase original data
          catch(const std::exception& xp){
             if (!ignore_missing)
                throw(std::runtime_error(string("Failure in doCorrelatorMatrixSuperpositionByBins: ")+xp.what()));}
#endif
          }}}
}


    //  In the pairs below, the bool specify is a vev is to be subtracted.
    //  A ratio is formed: numerator = interactionOpInfo, denominator is
    //  product of individual freeSingleOpInfos.  Each of the individual
    //  correlators is assumed to be real.

void doCorrelatorInteractionRatioBySamplings(MCObsHandler& moh,
                const pair<OperatorInfo,bool>& interactingOpInfo,
                const vector<pair<OperatorInfo,bool> >& freeSingleOpInfos,
                uint tmin, uint tmax, const OperatorInfo& ratio_op, 
                set<MCObsInfo>& obskeys, bool erase_orig)
{
 obskeys.clear();
 uint nterms=freeSingleOpInfos.size();
 if (nterms<2) throw(std::invalid_argument("Two or more FreeSingleOperators required"));
 bool herm=true;
 bool overwrite=true;
 CorrelatorAtTimeInfo ratiocorr(ratio_op,ratio_op,0,herm,false);
 CorrelatorAtTimeInfo intcorr(interactingOpInfo.first,interactingOpInfo.first,0,herm,interactingOpInfo.second);
 vector<CorrelatorAtTimeInfo> freecorrs;
 for (vector<pair<OperatorInfo,bool> >::const_iterator
        st=freeSingleOpInfos.begin();st!=freeSingleOpInfos.end();st++){
    freecorrs.push_back(CorrelatorAtTimeInfo(st->first,st->first,0,herm,st->second));}
 for (uint t=tmin;t<=tmax;t++){
    try{
    ratiocorr.resetTimeSeparation(t);
    intcorr.resetTimeSeparation(t);
    MCObsInfo intkey(intcorr);
    MCObsInfo ratiokey(ratiocorr);
    vector<MCObsInfo> freekeys;
    for (uint k=0;k<nterms;k++){
       freecorrs[k].resetTimeSeparation(t);
       freekeys.push_back(MCObsInfo(freecorrs[k]));}
    for (moh.begin();!moh.end();++moh){
       double ratioval=moh.getCurrentSamplingValue(intkey);
       for (uint k=0;k<nterms;k++){
          double denomval=moh.getCurrentSamplingValue(freekeys[k]);
          ratioval/=denomval;}
       moh.putCurrentSamplingValue(ratiokey,ratioval,overwrite);}
    obskeys.insert(ratiokey);
    if (erase_orig){
       moh.eraseSamplings(intkey);
       for (uint k=0;k<nterms;k++){
          moh.eraseSamplings(freekeys[k]);}}}
  catch(const std::exception& xp){}}
}


void doReconstructEnergyBySamplings(MCObsHandler& moh, const MCObsInfo& energy_diff_key,
			            const MCObsInfo& anisotropy_key, 
                                    const list<pair<MCObsInfo,double> >& scattering_particles,
			            const MCObsInfo& energy_res)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double xi=moh.getCurrentSamplingValue(anisotropy_key);
    double e_diff=moh.getCurrentSamplingValue(energy_diff_key);
    double energy = e_diff;
    for (list<pair<MCObsInfo,double> >::const_iterator it=scattering_particles.begin();
         it!=scattering_particles.end();it++){
      double scat_mass=moh.getCurrentSamplingValue(it->first);
      energy += sqrt(scat_mass*scat_mass + it->second/(xi*xi));}
    moh.putCurrentSamplingValue(energy_res,energy);}
}


void doReconstructEnergyBySamplings(MCObsHandler& moh, const MCObsInfo& energy_diff_key,
                                    const list<pair<MCObsInfo,double> >& scattering_particles, 
                                    const MCObsInfo& energy_res)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double e_diff=moh.getCurrentSamplingValue(energy_diff_key);
    double energy = e_diff;
    for (list<pair<MCObsInfo,double> >::const_iterator it=scattering_particles.begin();
         it!=scattering_particles.end();it++){
      double scat_mass=moh.getCurrentSamplingValue(it->first);
      energy += sqrt(scat_mass*scat_mass + it->second);}
    moh.putCurrentSamplingValue(energy_res,energy);}
}

void doReconstructAmplitudeBySamplings(MCObsHandler& moh, const MCObsInfo& energy_diff_amp_key,
                                       const std::list<MCObsInfo>& scattering_particles_amps, 
                                       const MCObsInfo& amp_res)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double e_diff_amp=moh.getCurrentSamplingValue(energy_diff_amp_key);
    double amplitude = e_diff_amp;
    for (list<MCObsInfo>::const_iterator it=scattering_particles_amps.begin();it!=scattering_particles_amps.end();it++){
      double scat_amp=moh.getCurrentSamplingValue(*it);
      amplitude *= scat_amp;}
    moh.putCurrentSamplingValue(amp_res,amplitude);}
}

void doEnergyDifferenceBySamplings(MCObsHandler& moh, const MCObsInfo& energy_key,
			            const MCObsInfo& anisotropy_key, 
                                    const list<pair<MCObsInfo,double> >& scattering_particles,
			            const MCObsInfo& energy_diff_res)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double xi=moh.getCurrentSamplingValue(anisotropy_key);
    double energy=moh.getCurrentSamplingValue(energy_key);
    double e_diff = energy;
    for (list<pair<MCObsInfo,double> >::const_iterator it=scattering_particles.begin();
         it!=scattering_particles.end();it++){
      double scat_mass=moh.getCurrentSamplingValue(it->first);
      e_diff -= sqrt(scat_mass*scat_mass + it->second/(xi*xi));}
    moh.putCurrentSamplingValue(energy_diff_res,e_diff);}
}


void doEnergyDifferenceBySamplings(MCObsHandler& moh, const MCObsInfo& energy_key,
                                    const list<pair<MCObsInfo,double> >& scattering_particles, 
                                    const MCObsInfo& energy_diff_res)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double energy=moh.getCurrentSamplingValue(energy_key);
    double e_diff = energy;
    for (list<pair<MCObsInfo,double> >::const_iterator it=scattering_particles.begin();
         it!=scattering_particles.end();it++){
      double scat_mass=moh.getCurrentSamplingValue(it->first);
      e_diff -= sqrt(scat_mass*scat_mass + it->second);}
    moh.putCurrentSamplingValue(energy_diff_res,e_diff);}
}

void doCorrelatedDifferenceBySamplings(MCObsHandler& moh, const MCObsInfo& obs1,
                                       const MCObsInfo& obs2, const MCObsInfo& diff)
{
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    double obs1_samp_val=moh.getCurrentSamplingValue(obs1);
    double obs2_samp_val=moh.getCurrentSamplingValue(obs2);
    moh.putCurrentSamplingValue(diff,obs1_samp_val-obs2_samp_val);}
}

   //  Given a list of CorrelatorInfo objects, this returns a list of index pairs
   //  specifying where the CorrelatorInfo objects occur in a CorrelatorMatrixInfo.

std::list<std::pair<uint,uint> > getCorrelatorMatrixIndices(const CorrelatorMatrixInfo& cormat, 
                                         const std::list<CorrelatorInfo>& corinfos, bool herm)
{
 std::list<std::pair<uint,uint> > elem_indices;
 const std::set<OperatorInfo>& opset=cormat.getOperators();
 std::vector<OperatorInfo> opvec(opset.begin(),opset.end());
 std::vector<OperatorInfo>::iterator vt1,vt2;
 for (list<CorrelatorInfo>::const_iterator ct=corinfos.begin();ct!=corinfos.end();++ct){
    vt1=std::find(opvec.begin(),opvec.end(), ct->getSink());
    vt2=std::find(opvec.begin(),opvec.end(), ct->getSource());
    if ((vt1!=opvec.end())&&(vt2!=opvec.end())){
       uint row=std::distance(opvec.begin(),vt1);
       uint col=std::distance(opvec.begin(),vt2);
       elem_indices.push_back(make_pair(row,col));
       if ((herm)&&(row!=col))
          elem_indices.push_back(make_pair(col,row));}}
 return elem_indices;
}

// ********************************************************************

     // Reads an original correlator matrix into memory, does the transformation
     // and puts results into memory.  If "eraseorig" true, the original is erased from
     // memory.  All keys associated with transformed correlator matrix and its vevs
     // that have been put in memory are returned in "obskeys".  This could then be
     // used to write to file.  The columns of T contain the coefficients of the
     // transformed **creation** operators in terms of the original creation operators.

void doTransformedCorrMatrixByBins(MCObsHandler& moh, const vector<OperatorInfo>& origopsvec,
                                   const vector<OperatorInfo>& transopsvec, bool herm, bool subvev,
                                   uint tmin, uint tmax, const TransMatrix& T, bool eraseorig,
                                   set<MCObsInfo>& obskeys)
{
 uint norig=origopsvec.size();
 uint ntrans=transopsvec.size();
 if ((norig<ntrans)||(T.size(0)!=norig)||(T.size(1)!=ntrans))
    throw(std::invalid_argument("Bad sizes in getTransformedCorrMatrixByBins"));
 uint nbins=moh.getNumberOfBins();
 vector<const RVector*> origbins_re(norig*norig);
 vector<RVector> transbins_re(ntrans*ntrans,RVector(nbins));
#ifdef COMPLEXNUMBERS
 vector<const RVector*> origbins_im(norig*norig);
 vector<RVector> transbins_im(ntrans*ntrans,RVector(nbins));
 CMatrix cormat;
#else
 RMatrix cormat;
#endif
 for (uint t=tmin;t<=tmax;t++){
    uint k=0;
    for (uint snk=0;snk<norig;++snk)
    for (uint src=0;src<norig;++src,++k){
       CorrelatorAtTimeInfo origcorr(origopsvec[snk],origopsvec[src],t,herm,false);
#ifdef COMPLEXNUMBERS
       origbins_im[k]=&(moh.getBins(MCObsInfo(origcorr,ImaginaryPart))); 
#endif
       origbins_re[k]=&(moh.getBins(MCObsInfo(origcorr)));}  // load into memory
    for (uint bin=0;bin<nbins;++bin){
       cormat.resize(norig,norig);
       k=0;
       for (uint snk=0;snk<norig;++snk)
       for (uint src=0;src<norig;++src,++k){
#ifdef COMPLEXNUMBERS
          if ((!herm)||(snk!=src))
             cormat.put(snk,src,complex<double>((*origbins_re[k])[bin],(*origbins_im[k])[bin]));
          else
             cormat.put(snk,src,complex<double>((*origbins_re[k])[bin],0.0));
#else
          cormat(snk,src)=(*origbins_re[k])[bin];
#endif
          }
       doMatrixRotation(cormat,T);
       k=0;
       for (uint snk=0;snk<ntrans;++snk)
       for (uint src=0;src<ntrans;++src,++k){
#ifdef COMPLEXNUMBERS
          transbins_re[k][bin]=cormat(snk,src).real();
          if ((!herm)||(src!=snk))
             transbins_im[k][bin]=cormat(snk,src).imag();
          else
             transbins_im[k][bin]=0.0;
#else                
          transbins_re[k][bin]=cormat(snk,src);
#endif
          }}
    k=0;
    for (uint snk=0;snk<ntrans;++snk)
    for (uint src=0;src<ntrans;++src,++k){
       if ((!herm)||(src<=snk)){
          CorrelatorAtTimeInfo transcorr(transopsvec[snk],transopsvec[src],t,herm,false);
          MCObsInfo key(transcorr);
          moh.putBins(key,transbins_re[k]);
          obskeys.insert(key);
#ifdef COMPLEXNUMBERS
          MCObsInfo imkey(transcorr,ImaginaryPart);
          moh.putBins(imkey,transbins_im[k]);
          obskeys.insert(imkey);
#endif
       }}
    if (eraseorig){
       for (uint snk=0;snk<norig;++snk)
       for (uint src=0;src<norig;++src){
          CorrelatorAtTimeInfo origcorr(origopsvec[snk],origopsvec[src],t,herm,false);
#ifdef COMPLEXNUMBERS
          moh.eraseData(MCObsInfo(origcorr,ImaginaryPart));
#endif
          moh.eraseData(MCObsInfo(origcorr));}}}
 if (!subvev) return;

    // rotate vevs if available
 origbins_re.resize(norig);
 transbins_re.resize(ntrans,RVector(nbins));
#ifdef COMPLEXNUMBERS
 origbins_im.resize(norig);
 transbins_im.resize(ntrans,RVector(nbins));
 CVector vevs;
#else
 RVector vevs;
#endif
 for (uint k=0;k<norig;++k){
#ifdef COMPLEXNUMBERS
    origbins_im[k]=&(moh.getBins(MCObsInfo(origopsvec[k],ImaginaryPart)));
#endif
    origbins_re[k]=&(moh.getBins(MCObsInfo(origopsvec[k])));}  // load into memory
 for (uint bin=0;bin<nbins;++bin){
    vevs.resize(norig);
    for (uint k=0;k<norig;++k){
#ifdef COMPLEXNUMBERS
       vevs.put(k,complex<double>((*origbins_re[k])[bin],(*origbins_im[k])[bin]));
#else
       vevs[k]=(*origbins_re[k])[bin];
#endif
       }
    doVectorRotation(vevs,T);
    for (uint k=0;k<ntrans;++k){
#ifdef COMPLEXNUMBERS
       transbins_re[k][bin]=vevs[k].real();
       transbins_im[k][bin]=vevs[k].imag();
#else                
       transbins_re[k][bin]=vevs[k];
#endif
       }}
 for (uint k=0;k<ntrans;++k){
    MCObsInfo key(transopsvec[k]);
    moh.putBins(key,transbins_re[k]);
    obskeys.insert(key);
#ifdef COMPLEXNUMBERS
    MCObsInfo imkey(transopsvec[k],ImaginaryPart);
    moh.putBins(imkey,transbins_im[k]);
    obskeys.insert(imkey);
#endif
    }
 if (eraseorig){
    for (uint k=0;k<norig;++k){
#ifdef COMPLEXNUMBERS
       moh.eraseData(MCObsInfo(origopsvec[k],ImaginaryPart));
#endif
       moh.eraseData(MCObsInfo(origopsvec[k]));}}
}



void doTransformedCorrMatrixBySamplings(MCObsHandler& moh, const vector<OperatorInfo>& origopsvec,
                                        const vector<OperatorInfo>& transopsvec, bool herm, bool subvev,
                                        uint tmin, uint tmax, const TransMatrix& T, bool eraseorig,
                                        set<MCObsInfo>& obskeys, bool separate_vevs)
{
 uint norig=origopsvec.size();
 uint ntrans=transopsvec.size();
 if ((norig<ntrans)||(T.size(0)!=norig)||(T.size(1)!=ntrans))
    throw(std::invalid_argument("Bad sizes in getTransformedCorrMatrixBySamplings"));
 vector<bool> svevs;
 if (!subvev) svevs.push_back(false);
 else if (!separate_vevs) svevs.push_back(true);
 else{ svevs.push_back(false); svevs.push_back(true);}
#ifdef COMPLEXNUMBERS
 CMatrix cormat;
#else
 RMatrix cormat;
#endif
 for (uint t=tmin;t<=tmax;t++){
  for (uint vevind=0;vevind<svevs.size();++vevind){  
    bool svev=svevs[vevind];
    for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
       cormat.resize(norig,norig);
       for (uint snk=0;snk<norig;++snk)
       for (uint src=0;src<norig;++src){
          CorrelatorAtTimeInfo origcorr(origopsvec[snk],origopsvec[src],t,herm,svev);
#ifdef COMPLEXNUMBERS
          if ((!herm)||(snk!=src))
             cormat.put(snk,src,complex<double>(moh.getCurrentSamplingValue(MCObsInfo(origcorr)),
                                                moh.getCurrentSamplingValue(MCObsInfo(origcorr,ImaginaryPart))));
          else
             cormat.put(snk,src,complex<double>(moh.getCurrentSamplingValue(MCObsInfo(origcorr)),0.0));
#else
          cormat(snk,src)=moh.getCurrentSamplingValue(MCObsInfo(origcorr));
#endif
          }
       doMatrixRotation(cormat,T);
       for (uint snk=0;snk<ntrans;++snk)
       for (uint src=0;src<ntrans;++src)
          if ((!herm)||(src<=snk)){
             CorrelatorAtTimeInfo transcorr(transopsvec[snk],transopsvec[src],t,herm,svev);
             MCObsInfo key(transcorr);
             obskeys.insert(key);
#ifdef COMPLEXNUMBERS
             moh.putCurrentSamplingValue(key,cormat(snk,src).real());
             MCObsInfo imkey(transcorr,ImaginaryPart);
             if ((!herm)||(src!=snk))
                moh.putCurrentSamplingValue(imkey,cormat(snk,src).imag());
             else
                moh.putCurrentSamplingValue(imkey,0.0);
             obskeys.insert(imkey);
#else
             moh.putCurrentSamplingValue(key,cormat(snk,src));
#endif
          }}
    if (eraseorig){
       for (uint snk=0;snk<norig;++snk)
       for (uint src=0;src<norig;++src){
          CorrelatorAtTimeInfo origcorr(origopsvec[snk],origopsvec[src],t,herm,svev);
#ifdef COMPLEXNUMBERS
          moh.eraseSamplings(MCObsInfo(origcorr,ImaginaryPart));
#endif
          moh.eraseSamplings(MCObsInfo(origcorr));}}}}
 if ((!subvev)||(!separate_vevs)) return;

#ifdef COMPLEXNUMBERS
 CVector vevs;
#else
 RVector vevs;
#endif
    // rotate vevs if available
 for (moh.setSamplingBegin();!moh.isSamplingEnd();moh.setSamplingNext()){
    vevs.resize(norig);
    for (uint k=0;k<norig;++k){
#ifdef COMPLEXNUMBERS
       vevs.put(k,complex<double>(moh.getCurrentSamplingValue(MCObsInfo(origopsvec[k])),
                                  moh.getCurrentSamplingValue(MCObsInfo(origopsvec[k],ImaginaryPart))));
#else
       vevs[k]=moh.getCurrentSamplingValue(MCObsInfo(origopsvec[k]));
#endif
       }
    doVectorRotation(vevs,T);
    for (uint k=0;k<ntrans;++k){
       MCObsInfo key(transopsvec[k]);
       obskeys.insert(key);
#ifdef COMPLEXNUMBERS
       moh.putCurrentSamplingValue(key,vevs[k].real());
       MCObsInfo imkey(transopsvec[k],ImaginaryPart);
       moh.putCurrentSamplingValue(imkey,vevs[k].imag());
       obskeys.insert(imkey);
#else
       moh.putCurrentSamplingValue(key,vevs[k]);
#endif
       }}
 if (eraseorig){
    for (uint k=0;k<norig;++k){
#ifdef COMPLEXNUMBERS
       moh.eraseSamplings(MCObsInfo(origopsvec[k],ImaginaryPart));
#endif
       moh.eraseSamplings(MCObsInfo(origopsvec[k]));}}
}

// ******************************************************************************************

void setUpRatioFit(MCObsHandler& m_obs, XMLHandler& xmlf, XMLHandler& xmltf, bool write_output, XMLHandler& xmlout, bool erase_orig){
     XMLHandler xmlres(xmlf,"Ratio");
     OperatorInfo ratio_op(xmlres);
     XMLHandler xmlint(xmlf,"InteractingOperator");
     bool numvev=(xmlint.count("SubtractVEV")>0) ? true: false;
     pair<OperatorInfo,bool> numerator=make_pair(OperatorInfo(xmlint),numvev);
     vector<pair<OperatorInfo,bool> > denominator;
     list<XMLHandler> denomxml=xmlf.find_among_children("NonInteractingOperator");
     for (list<XMLHandler>::iterator it=denomxml.begin();it!=denomxml.end();++it){
       OperatorInfo opinfo(*it);
       bool subvev=(it->count("SubtractVEV")>0) ? true: false;
       denominator.push_back(make_pair(opinfo,subvev));}
     uint nterms=denominator.size();
     if (nterms<2) throw(std::invalid_argument("Two or more NonInteractingOperators required"));
     uint tmin,tmax;
     xmlreadchild(xmlf,"MinimumTimeSeparation",tmin);
     xmlreadchild(xmlf,"MaximumTimeSeparation",tmax);
     XMLHandler xmlo, xmldp;
     xmlo.set_root("Ratio");
     ratio_op.output(xmldp);
     xmlo.put_child(xmldp);
     if(write_output) xmlout.put_child(xmlo);
     xmlo.set_root("InteractingOperator");
     numerator.first.output(xmldp);
     xmlo.put_child(xmldp);
     if(write_output) xmlout.put_child(xmlo);
     xmlo.set_root("NonInteractingOperators");
     for (vector<pair<OperatorInfo,bool> >::const_iterator
            it=denominator.begin();it!=denominator.end();it++){
        XMLHandler xmloo; it->first.output(xmloo); 
         xmlo.put_child(xmloo);
     }
     if(write_output) xmlout.put_child(xmlo);
     set<MCObsInfo> obskeys;
     doCorrelatorInteractionRatioBySamplings(m_obs,numerator,denominator,
                                             0,(tmax<64)?64:tmax,ratio_op,obskeys,erase_orig);
     xmltf.rename_tag("TemporalCorrelatorFit");
     XMLHandler xmlro; ratio_op.output(xmlro);
     xmltf.put_child(xmlro); 
}
