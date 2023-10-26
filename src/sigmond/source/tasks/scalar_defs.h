#ifndef SCALAR_DEFS_H
#define SCALAR_DEFS_H
#include <complex>
#include <string>

        //  Macros below refer to number types in the data files

   //  Either one of these macros must be defined (preferably in makefile)
//#define SINGLEPRECISION
//#define DOUBLEPRECISION

   //  Either one of these macros must be defined (preferably in makefile)
//#define COMPLEXNUMBERS
//#define REALNUMBERS

   //   The typedef "InScalar" defines the type used for basic
   //   values in the input files.  If results stored in the
   //   data files are single precision and complex, uncomment
   //   the line with complex<float>, and so on.  

   //   The typedef "Scalar" defines the type used for basic
   //   values in manipulating the data during analysis.
   //   Double precision is strongly recommended, so choose
   //   between real and complex values.

#if defined DOUBLEPRECISION

#if defined COMPLEXNUMBERS
  typedef std::complex<double>    InScalar;
  typedef std::complex<double>    Scalar;
#elif defined REALNUMBERS
  typedef double             InScalar;
  typedef double             Scalar;
#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif

#elif defined SINGLEPRECISION

#if defined COMPLEXNUMBERS
  typedef std::complex<float>     InScalar;
  typedef std::complex<double>    Scalar;
#elif defined REALNUMBERS
  typedef float              InScalar;
  typedef double             Scalar;
#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif

#else
  #error "Either DOUBLEPRECISION or SINGLEPRECISION must be defined"
#endif



   //  Compares the number type read from file in a
   //  <NumberType> tag to the expected type depending on
   //  version of "sigmond" compiled.  Throws exception on failure.

void check_number_type(const std::string& readtype);


double conjugate(const double& x);
float conjugate(const float& x);
std::complex<double> conjugate(const std::complex<double>& z);
std::complex<float> conjugate(const std::complex<float>& z);

double realpart(const double& x);
double realpart(const float& x);
double imaginarypart(const double& x);
double imaginarypart(const float& x);

double realpart(const std::complex<double>& z);
double realpart(const std::complex<float>& z);
double imaginarypart(const std::complex<double>& z);
double imaginarypart(const std::complex<float>& z);

double sqr(const double& x);
double sqr(const std::complex<double>& z);

enum ComplexArg { RealPart, ImaginaryPart };

typedef unsigned int   uint;

#endif


