#include "scalar_defs.h"
#include <iostream>
#include <stdexcept>

using namespace std;


void check_number_type(const string& readtype)
{
#ifdef SINGLEPRECISION
#ifdef COMPLEXNUMBERS
 if (readtype!=string("ComplexSinglePrecision")){
#endif
#ifdef REALNUMBERS
 if (readtype!=string("RealSinglePrecision")){
#endif
#endif

#ifdef DOUBLEPRECISION
#ifdef COMPLEXNUMBERS
 if (readtype!=string("ComplexDoublePrecision")){
#endif
#ifdef REALNUMBERS
 if (readtype!=string("RealDoublePrecision")){
#endif
#endif

   cout << endl<<endl;
   cout << "Number type in files is "<<readtype
        <<" does not match compiled version of sigmond"<<endl;
   cout << "Change macros in scalar_defs.h and recompile"<<endl;
   cout << endl<<endl;
   throw(std::runtime_error("Mismatch in number type"));}
}   


double conjugate(const double& x)
{
 return x;
}

float conjugate(const float& x)
{
 return x;
}

std::complex<double> conjugate(const std::complex<double>& z)
{
 return conj(z);
}

std::complex<float> conjugate(const std::complex<float>& z)
{
 return conj(z);
}

double realpart(const double& x)
{
 return x;
}

double realpart(const float& x)
{
 return double(x);
}

double imaginarypart(const double& x)
{
 return 0.0;
}

double imaginarypart(const float& x)
{
 return 0.0;
}

double realpart(const complex<double>& z)
{
 return real(z);
}

double realpart(const complex<float>& z)
{
 return double(real(z));
}

double imaginarypart(const complex<double>& z)
{
 return imag(z);
}

double imaginarypart(const complex<float>& z)
{
 return double(imag(z));
}

double sqr(const double& x)
{
 return x*x;
}

double sqr(const std::complex<double>& z)
{
 return std::norm(z);
}

