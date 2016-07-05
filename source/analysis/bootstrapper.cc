#include "bootstrapper.h"
using namespace std;


// *************************************************************


Bootstrapper::Bootstrapper(unsigned int num_objects, unsigned int num_samples, 
                           unsigned long seed, unsigned int skip_value, 
                           bool precompute)
       :  nobjects(num_objects), nsamples(num_samples),
          rngseed(seeder(seed)), nskip(skip_value), U(rngseed), 
          bitmask(calc_bitmask(num_objects)), counter(0), currbits(0), 
          neednew(true), m_precompute(precompute), m_current(0),
          m_samplings(((precompute)&&(num_samples<65536))?num_samples:1,
                       Vector<unsigned int>((num_objects<65536)?num_objects:1))
{
 try{
    if (nobjects>=65536) throw(std::invalid_argument("Too many objects in Bootstrapper"));
    if (nobjects<=4) throw(std::invalid_argument("Too few objects in Bootstrapper"));
    if (nsamples<=4) throw(std::invalid_argument("Too few resamplings in Bootstrapper"));
    if (nsamples>=65536) throw(std::invalid_argument("Too many resamplings in Bootstrapper"));}
 catch(const std::exception& errmsg){
    cout << "Bootstrapper error: "<<errmsg.what()<<endl;
    throw(std::invalid_argument((string("Invalid bootstrapper: ")+string(errmsg.what())).c_str()));}
 if (precompute){
    m_funcptr=&Bootstrapper::get_resampling_precalc;
    for (unsigned int k=0;k<nsamples;k++){
       m_current=&(m_samplings[k]);
       generate_resampling();}}
 else{
    m_funcptr=&Bootstrapper::get_resampling_onfly;
    m_current=&(m_samplings[0]);
    generate_resampling();}
}


Bootstrapper& Bootstrapper::reset(unsigned int num_objects, unsigned int num_samples, 
                                  unsigned long seed, unsigned int skip_value,
                                  bool precompute)
{
 m_samplings.clear();
 nobjects=num_objects;
 nsamples=num_samples;
 try{
    if (nobjects>=65536) throw(std::invalid_argument("Too many objects in Bootstrapper"));
    if (nobjects<=4) throw(std::invalid_argument("Too few objects in Bootstrapper"));
    if (nsamples<=4) throw(std::invalid_argument("Too few samples in Bootstrapper"));
    if (nsamples>=65536) throw(std::invalid_argument("Too many resamplings in Bootstrapper"));}
 catch(const std::exception& errmsg){
    cout << "Bootstrapper error: "<<errmsg.what()<<endl;
    throw(std::invalid_argument((string("Invalid bootstrapper reset")+string(errmsg.what())).c_str()));}
 rngseed=seeder(seed);
 nskip=skip_value;
 U.reseed(rngseed);
 bitmask=calc_bitmask(num_objects);
 counter=0;
 currbits=0;
 neednew=true;
 m_precompute=precompute;
 if (precompute){
    m_funcptr=&Bootstrapper::get_resampling_precalc;
    m_samplings.resize(nsamples,Vector<unsigned int>(nobjects));
    for (unsigned int k=0;k<nsamples;k++){
       m_current=&(m_samplings[k]);
       generate_resampling();}}
 else{
    m_funcptr=&Bootstrapper::get_resampling_onfly;
    m_samplings.resize(1,Vector<unsigned int>(nobjects));
    m_current=&(m_samplings[0]);
    generate_resampling();}
 return *this;
}


Bootstrapper::~Bootstrapper() {}


unsigned int Bootstrapper::getNumberOfObjects() const
{
 return nobjects;
}

unsigned int Bootstrapper::getNumberOfResamplings()  const
{
 return nsamples;
}

unsigned long Bootstrapper::getRNGSeed() const
{
 return rngseed;
}

unsigned int Bootstrapper::getSkipValue() const
{
 return nskip;
}

unsigned int Bootstrapper::getCurrentResamplingCount() const
{
 return counter;
}

bool Bootstrapper::isPrecomputeMode() const
{
 return m_precompute;
}

unsigned long Bootstrapper::seeder(unsigned long in_seed)
{
 unsigned long s(in_seed);
 if (s==0){
    s=time(NULL);
    s^=s<<8; s^=s<<8; s^=s<<8;}
 return s;
}

unsigned int Bootstrapper::calc_bitmask(unsigned int num_objects)
{
 unsigned int mask=1;
 while (mask<num_objects) mask*=2;
 return mask-1;
}

const Vector<unsigned int>& Bootstrapper::getCurrentResampling() const
{
 return *m_current;
}
    
const Vector<unsigned int>& Bootstrapper::getResampling(int ind)
{
 return (this->*m_funcptr)(ind);
}

const Vector<unsigned int>& Bootstrapper::get_resampling_precalc(int ind)
{
 if ((ind<0)||(ind>=int(nsamples)))
    throw(std::invalid_argument("Invalid sample index in Bootstrapper"));
 counter=ind;
 return m_samplings[ind];
}
 
const Vector<unsigned int>& Bootstrapper::get_resampling_onfly(int ind)
{
 if ((ind<0)||(ind>=int(nsamples)))
    throw(std::invalid_argument("Invalid sample index in Bootstrapper"));
 if (ind==int(counter)) return *m_current;
 int nadvance=ind-counter-1;
 if (int(counter)>ind){
    U.reseed(rngseed);
    neednew=true;
    nadvance=ind;}
 for (int a=0;a<nadvance;++a) advance_resampling();
 generate_resampling();
 counter=ind;
 return *m_current;
}


const Vector<unsigned int>& Bootstrapper::getNextResampling()
{
 unsigned int next=counter+1;
 if (next==nsamples) next=0;
 return getResampling(next);
}

unsigned int Bootstrapper::generate_with_mask()
{
 if (neednew){
    currbits=U.generate();
    neednew=false;}
 else{
    currbits>>=16;
    neednew=true;}
 return bitmask & currbits;
}

unsigned int Bootstrapper::generate()
{
 unsigned int result=generate_with_mask();
 while (result>=nobjects)
    result=generate_with_mask();
 return result;
}

void Bootstrapper::advance_resampling()
{
 for (unsigned int k=0;k<(nobjects+nskip);++k)
    generate();
}

void Bootstrapper::generate_resampling()
{
 for (unsigned int k=0;k<nobjects;++k)
    (*m_current)[k]=generate();
 for (unsigned int k=0;k<nskip;++k) generate();
}


// *************************************************************


UniformDeviate32::UniformDeviate32(uint32 seed)
{
 reseed(seed);
}

 // Set initial "state" using "seed".  If seed==0, then 
 // initialization is done using the system clock.
 // NOTES on hexidecimal numbers:  a constant beginning
 // with "0x" indicates a hexidecimal number (base 16).
 // The "digits" of a hexidecimal number are 
 // 0-9,a,b,c,d,e,f (or upppercase).  A "u" or "U" then
 // indicates it is unsigned, and an "L" indicates it is
 // a long integer.

void UniformDeviate32::reseed(uint32 seed)
{ 
 uint32 s=seed;
 int j;
 state[0]= s & 0xffffffffUL;  // 2^32
 for (j=1;j<N;j++) {
    state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
    state[j] &= 0xffffffffUL;  // for >32 bit machines
    }
 left = 0;
}

// *************************************************************
