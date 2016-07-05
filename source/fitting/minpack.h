#ifndef MINPACK_H
#define MINPACK_H
#include <cmath>

namespace MinPack {

void qrfac(int m, int n, double *a, int lda, int pivot, int *ipvt,
           int lipvt, double *rdiag, double *acnorm, double *wa);

void qrsolv(int n, double *r, int ldr, const int *ipvt, const double *diag,
            const double *qtb, double *x, double *sdiag, double *wa);

double dpmpar(int i);

double enorm(int n, double *x);

inline double abs(double x)
{return (x>=0.0) ? x : -x;}

inline double sqr(double x)
{return x*x;}

inline int mod(int a, int b)
{return a%b;}

inline double max(double a, double b)
{return (a>b) ? a : b;}

inline double min(double a, double b)
{return (a<b) ? a : b;}

inline int min(int i, int j)
{return (i<j) ? i : j;}

inline int max(int a, int b)
{return (a>b) ? a : b;}

inline double pow_dd(double *x, double *y)
{return pow(*x,*y);}

const int True=1;
const int False=0;

}
#endif
