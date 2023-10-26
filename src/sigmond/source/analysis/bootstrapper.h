#ifndef BOOTSTRAPPER_H
#define BOOTSTRAPPER_H
#include "matrix.h"
#include <vector>


typedef unsigned long   uint32;
class UniformDeviate32;


  // ******************************************************************
  // *                                                                *
  // *   The class "Bootstrapper" is defined in this file.  The       *
  // *   constructor requires the following information:              *
  // *                                                                *
  // *      - number of objects (such as Monte Carlo measurements)    *
  // *      - number of bootstrap resamplings to be done              *
  // *      - a random number generator seed (32 bit unsigned int)    *
  // *      - a skip unsigned integer (0 to 1024 usually)             *
  // *      - a boolean "precompute"                                  *
  // *                                                                *
  // *   A maximum number of 65536 objects are allowed by this class. *
  // *   For a given bootstrap resampling, the Mersenne twister is    *
  // *   called in sequence to generate the sample. Before generating *
  // *   the next resampling, the Mersenne twister is called "skip"   *
  // *   number of times.  The "skip" parameter allows more variety   *
  // *   in how the bootstrap resamplings are generated.              *
  // *                                                                *
  // *   Once defined, the main member to be used is "getResampling". *
  // *   This member takes an unsigned integer index as argument,     *
  // *   which must range from 0 to one less than the number of       *
  // *   bootstrap resamplings.                                       *
  // *                                                                *
  // *   If "precompute" is set to true, then the bootstrap indices   *
  // *   are computed all at once and stored in memory.  If set       *
  // *   to false, the bootstrap indices are computed "on the fly".   *
  // *   The class has defined a unique way of generating the samples *
  // *   that is repeatable, using the Mersenne twister. So a given   *
  // *   resampling can be regenerated as needed. Setting "precompute"*
  // *   to true uses more memory.  Try both methods to see which is  *
  // *   faster, since speed depends on memory access speed.          *
  // *                                                                *
  // *   Example usage:                                               *
  // *                                                                *
  // *    unsigned int num_objects=1584;   // MC estimates            *
  // *    unsigned int num_samples=1024;   // bootstrap samples       *
  // *    unsigned long rngseed=65231;                                *
  // *    unsigned int skip=723;                                      *
  // *    bool precompute=false;                                      *
  // *                                                                *
  // *    Bootstrapper B(num_objects,num_samples,                     *
  // *                   rngseed,skip,precompute);                    *
  // *                                                                *
  // *    for (int k=0;k<B.getNumberOfResamplings();++k){             *
  // *                                                                *
  // *       const Vector<uint>& sample=B.getResampling(k);           *
  // *                                                                *
  // *        ....use sample... }                                     *
  // *                                                                *
  // *   This class maintains the memory for a given sample as a      *
  // *   vector of unsigned integers.  Other notable members:         *
  // *                                                                *
  // *      B.getNumberOfObjects()                                    *
  // *      B.getNumberOfReSamplings()                                *
  // *      B.getRNGSeed()                                            *
  // *      B.getSkipValue()                                          *
  // *      B.getCurrentSampleIndex()                                 *
  // *                                                                *
  // *                                                                *
  // ******************************************************************



  // ****************************************************************
  // *                                                              *
  // *   Random number generator:  the 32-bit Mersenne Twister      *
  // *   Using objects of class "UniformDeviate32":                 *
  // *                                                              *
  // *   UniformDeviate32 Q;        => defines Q (time-based seed)  *
  // *   UniformDeviate32 Q(seed);  => or use explicit seed         *
  // *   Q.reseed(seed);  => explicit reseeding                     *
  // *   Q.generate();    => returns random 32-bit unsigned int     *
  // *                                                              *
  // ****************************************************************

// The celebrated Mersenne Twister.
// Mersenne Twister random number generator.
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v0.8  24 March 2002  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.
// It was designed with consideration of the flaws in various other
// generators. The period, 2^19937-1, and the order of equidistribution,
// 623 dimensions, are far greater.  The generator is also fast; it
// avoids multiplication and division, and it benefits from caches and
// pipelines. 

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.



class UniformDeviate32
{

 private:

    typedef unsigned long uint32;  // unsigned integer type, at least 32 bits
    static const int N = 624;      // length of state vector
    static const int M = 397;      // period parameter
    static const uint32 MATRIX_A = 0x9908b0dfUL;  // constant vector a
    static const uint32 UMASK = 0x80000000UL;     // most significant w-r bits
    static const uint32 LMASK = 0x7fffffffUL;     // least significant r bits
                      // hexadecimal constants start with 0x, UL=unsigned long

    uint32 state[N];  // internal state
    uint32 *next;     // next value to get from state
    int left;         // number of values left before next_state needed

 public:
          
    UniformDeviate32(uint32 seed);  // Constructor with explicit seed (0 <= seed < 2^32).
    ~UniformDeviate32(){}           // Destructor.
    void reseed(uint32 seed);       // Explicitly re-seeds the generator (0 <= seed < 2^32   or  seed < 0 to seed by time).  
   
    uint32 generate();              // generates random 32 bit unsigned int

 private:

    uint32 twist( const uint32& m, const uint32& s0, const uint32& s1 ) const
      { return m ^ (((s0 & UMASK) | (s1 & LMASK)) >> 1) ^ (-(s1 & 1UL) & MATRIX_A); }

            // prevent copying
   UniformDeviate32(const UniformDeviate32& indata);
   UniformDeviate32& operator=(const UniformDeviate32& indata);

};



  // ********************************************
  // *                                          *
  // *   The class "Bootstrapper"               *
  // *                                          *
  // ********************************************


class Bootstrapper
{

    unsigned int nobjects;
    unsigned int nsamples;
    unsigned long rngseed;
    unsigned int nskip;

    UniformDeviate32 U;
    unsigned int bitmask;
    unsigned int counter;    // current resampling index
    unsigned long currbits;  // current randomly generated integer value
    bool neednew;            // Mersenne twister generates 32 bits; use 16 at a time
    bool m_precompute;
    Vector<unsigned int> *m_current;
    std::vector<Vector<unsigned int> > m_samplings;

    const Vector<unsigned int>& (Bootstrapper::*m_funcptr)(int);  // resampling function pointer


 public:

    Bootstrapper(unsigned int num_objects, unsigned int num_resamplings, 
                 unsigned long seed=0, unsigned int skip_value=0, bool precompute=false);

    Bootstrapper& reset(unsigned int num_objects, unsigned int num_resamplings, 
                        unsigned long seed=0, unsigned int skip_value=0, 
                        bool precompute=false);
    ~Bootstrapper();

    unsigned int getNumberOfObjects() const;
    unsigned int getNumberOfResamplings()  const;
    unsigned long getRNGSeed() const;
    unsigned int getSkipValue() const;
    unsigned int getCurrentResamplingCount() const;
    bool isPrecomputeMode() const;

    const Vector<unsigned int>& getResampling(int ind);
    const Vector<unsigned int>& getNextResampling();
    const Vector<unsigned int>& getCurrentResampling() const;

 private:

    unsigned long seeder(unsigned long in_seed);
    unsigned int calc_bitmask(unsigned int num_objects);
    unsigned int generate_with_mask();
    unsigned int generate();
    void advance_resampling();
    void generate_resampling();
     
    const Vector<unsigned int>& get_resampling_precalc(int ind);
    const Vector<unsigned int>& get_resampling_onfly(int ind);

            // prevent copying
   Bootstrapper(const Bootstrapper& indata);
   Bootstrapper& operator=(const Bootstrapper& indata);
};


  // ****************************************************************


inline UniformDeviate32::uint32 UniformDeviate32::generate()
{
 if (left==0){
    int i;
    uint32 *p=state;
    left=N;
    next=state;
    for (i=N-M+1;--i;p++)
       *p=twist(p[M],p[0],p[1]);
    for (i=M;--i;p++)
       *p=twist(p[M-N],p[0],p[1]);
    *p=twist(p[M-N],p[0],state[0]);
    }
 left--;
 uint32 s1 = *next++;
 s1 ^= (s1 >> 11);
 s1 ^= (s1 <<  7) & 0x9d2c5680UL;
 s1 ^= (s1 << 15) & 0xefc60000UL;
 s1 ^= (s1 >> 18);
 return s1;
}


// *************************************************************
#endif
