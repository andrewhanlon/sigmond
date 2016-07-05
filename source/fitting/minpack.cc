#include "minpack.h"

namespace MinPack {

// ******************************************************************
//
//     function dpmpar
//
//     This function provides double machine parameters
//     when the appropriate set of data statements is activated (by
//     removing the c from column 1) and all other data statements are
//     rendered inactive. Most of the parameter values were obtained
//     from the corresponding Bell Laboratories Port Library function.
//
//     The function statement is
//
//       double function dpmpar[i]
//
//     where
//
//       i is an int input variable set to 1, 2, or 3 which
//         selects the desired machine parameter. If the machine has
//         t base b digits and its smallest and largest exponents are
//         emin and emax, respectively, then these parameters are
//
//         dpmpar(1) = b**(1 - t), the machine precision (mcheps),
//
//         dpmpar(2) = b**(emin - 1), the smallest magnitude (minmag),
//
//         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude (maxmag).
//
//     Argonne National Laboratory. MINPACK Project. June 1983.
//     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
//
//     **********

double dpmpar(int i)
{
 if (i==1) return 2.22e-16;
 else if (i==2) return 2.23e-308;
 else return 1.79e308;              // DEC alpha series.
}

// ********************************************************************
//
//     function enorm
//
//     given an n-vector x, this function calculates the
//     euclidean norm of x.
//
//     the euclidean norm is computed by accumulating the sum of
//     squares in three different sums. the sums of squares for the
//     small and large components are scaled so that no overflows
//     occur. non-destructive underflows are permitted. underflows
//     and overflows do not occur in the computation of the unscaled
//     sum of squares for the intermediate components.
//     the definitions of small, intermediate and large components
//     depend on two constants, rdwarf and rgiant. the main
//     restrictions on these constants are that rdwarf**2 not
//     underflow and rgiant**2 not overflow. the constants
//     given here are suitable for every known computer.
//
//     the function statement is
//
//       double function enorm(n,x)
//
//     where
//
//       n is a positive int input variable.
//
//       x is an input array of length n.
//
//     subprograms called
//
//       fortran-supplied ... fabs,sqrt
//
//     argonne national laboratory. minpack project. march 1980.
//     burton s. garbow, kenneth e. hillstrom, jorge j. more
//
//     **********

double enorm(int n, double *x)
{
 int i;
 double agiant, s1, s2, s3, d1, xabs, x1max, x3max;
 double rdwarf=3.834e-20, rgiant=1.304e19;
 double ans=0.0;
 
 s1 = 0.0;
 s2 = 0.0;
 s3 = 0.0;
 x1max = 0.0;
 x3max = 0.0;
 agiant = rgiant / double(n);
 for (i = 0; i < n; ++i) {
     xabs = fabs(x[i]);
     if (xabs >= agiant) {
          // sum for large components.
         if (xabs > x1max) {
             // Computing 2nd power
             d1 = x1max / xabs;
             s1 = 1. + s1 * (d1 * d1);
             x1max = xabs;
         } else {
             // Computing 2nd power
             d1 = xabs / x1max;
             s1 += d1 * d1;
         }
     } else if (xabs <= rdwarf) {
           // sum for small components. 
         if (xabs > x3max) {
             // Computing 2nd power 
             d1 = x3max / xabs;
             s3 = 1. + s3 * (d1 * d1);
             x3max = xabs;
         } else if (xabs != 0.) {
             // Computing 2nd power 
             d1 = xabs / x3max;
             s3 += d1 * d1;
         }
     } else {
         // sum for intermediate components.
         // Computing 2nd power
         s2 += xabs * xabs;
     }
 }

   //  calculation of norm.

 if (s1 != 0.) {
     ans = x1max * sqrt(s1 + (s2 / x1max) / x1max);
 } else if (s2 != 0.) {
     if (s2 >= x3max) {
         ans = sqrt(s2 * (1. + (x3max / s2) * (x3max * s3)));
     } else {
         ans = sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
     }
 } else {
     ans = x3max * sqrt(s3);
 }
 return ans;
}


// ******************************************************************
//
//     subroutine qrfac
//
//     this subroutine uses householder transformations with column
//     pivoting (optional) to compute a qr factorization of the
//     m by n matrix a. that is, qrfac determines an orthogonal
//     matrix q, a permutation matrix p, and an upper trapezoidal
//     matrix r with diagonal elements of nonincreasing magnitude,
//     such that a*p = q*r. the householder transformation for
//     column k, k = 1,2,...,min(m,n), is of the form
//
//                           t
//           i - (1/u(k))*u*u
//
//     where u has zeros in the first k-1 positions. the form of
//     this transformation and the method of pivoting first
//     appeared in the corresponding linpack subroutine.
//
//     the subroutine statement is
//
//       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
//
//     where
//
//       m is a positive int input variable set to the number
//         of rows of a.
//
//       n is a positive int input variable set to the number
//         of columns of a.
//
//       a is an m by n array. on input a contains the matrix for
//         which the qr factorization is to be computed. on output
//         the strict upper trapezoidal part of a contains the strict
//         upper trapezoidal part of r, and the lower trapezoidal
//         part of a contains a factored form of q (the non-trivial
//         elements of the u vectors described above).
//
//       lda is a positive int input variable not less than m
//         which specifies the leading dimension of the array a.
//
//       pivot is a logical input variable. if pivot is set true,
//         then column pivoting is enforced. if pivot is set false,
//         then no column pivoting is done.
//
//       ipvt is an int output array of length lipvt. ipvt
//         defines the permutation matrix p such that a*p = q*r.
//         column j of p is column ipvt[j] of the identity matrix.
//         if pivot is false, ipvt is not referenced.
//
//       lipvt is a positive int input variable. if pivot is false,
//         then lipvt may be as small as 1. if pivot is true, then
//         lipvt must be at least n.
//
//       rdiag is an output array of length n which contains the
//         diagonal elements of r.
//
//       acnorm is an output array of length n which contains the
//         norms of the corresponding columns of the input matrix a.
//         if this information is not needed, then acnorm can coincide
//         with rdiag.
//
//       wa is a work array of length n. if pivot is false, then wa
//         can coincide with rdiag.
//
//     subprograms called
//
//       minpack-supplied ... dpmpar,enorm
//
//       fortran-supplied ... max,sqrt,min
//
//     argonne national laboratory. minpack project. march 1980.
//     burton s. garbow, kenneth e. hillstrom, jorge j. more
//
//     **********

void qrfac(int m, int n, double *a, int lda, int pivot, int *ipvt,
           int lipvt, double *rdiag, double *acnorm, double *wa)
{
 int i,j,jp1,k,kmax,minmn;
 double d1,ajnorm,epsmch,sum,temp;
 const double p05=0.05;

         //  epsmch is the machine precision.

 epsmch = dpmpar(1);

     // compute the initial column norms and initialize several arrays.

 for (j = 0; j < n; ++j) {
    acnorm[j] = enorm(m, &a[j*lda]);
    rdiag[j] = acnorm[j];
    wa[j] = rdiag[j];
    if (pivot) ipvt[j] = j+1;
    }

//     reduce a to r with householder transformations.

 minmn = min(m,n);
 for (j = 0; j < minmn; ++j) {
    if (pivot) {

//        bring the column of largest norm into the pivot position.

       kmax = j;
       for (k = j; k < n; ++k) {
          if (rdiag[k] > rdiag[kmax]) kmax = k;}
       if (kmax != j) {
          for (i = 0; i < m; ++i) {
             temp = a[i + j * lda];
             a[i + j * lda] = a[i + kmax * lda];
             a[i + kmax * lda] = temp;}
          rdiag[kmax] = rdiag[j];
          wa[kmax] = wa[j];
          k = ipvt[j];
          ipvt[j] = ipvt[kmax];
          ipvt[kmax] = k;}}

//        compute the householder transformation to reduce the
//        j-th column of a to a multiple of the j-th unit vector.

    ajnorm = enorm(m - (j+1) + 1, &a[j + j * lda]);
    if (ajnorm != 0.0) {
       if (a[j + j * lda] < 0.) {
          ajnorm = -ajnorm;}
       for (i = j; i < m; ++i) {
          a[i + j * lda] /= ajnorm;}
       a[j + j * lda] += 1.0;

//        apply the transformation to the remaining columns 
//        and update the norms.

       jp1 = j + 1;
       if (n > jp1) {
          for (k = jp1; k < n; ++k) {
             sum = 0.;
             for (i = j; i < m; ++i) {
                sum += a[i + j * lda] * a[i + k * lda];}
             temp = sum / a[j + j * lda];
             for (i = j; i < m; ++i) {
                a[i + k * lda] -= temp * a[i + j * lda];}
             if (pivot && rdiag[k] != 0.) {
                temp = a[j + k * lda] / rdiag[k];
                        // Computing MAX 
                d1 = 1.0 - temp * temp;
                rdiag[k] *= sqrt((max(0.0,d1)));
                        // Computing 2nd power
                d1 = rdiag[k] / wa[k];
                if (p05 * (d1 * d1) <= epsmch) {
                   rdiag[k] = enorm(m - (j+1), &a[jp1 + k * lda]);
                   wa[k] = rdiag[k];}
             }}
          }
       }
    rdiag[j] = -ajnorm;
    }
}


// *****************************************************************
//
//     subroutine qrsolv
//
//     given an m by n matrix a, an n by n diagonal matrix d,
//     and an m-vector b, the problem is to determine an x which
//     solves the system
//
//           a*x = b ,     d*x = 0 ,
//
//     in the least squares sense.
//
//     this subroutine completes the solution of the problem
//     if it is provided with the necessary information from the
//     qr factorization, with column pivoting, of a. that is, if
//     a*p = q*r, where p is a permutation matrix, q has orthogonal
//     columns, and r is an upper triangular matrix with diagonal
//     elements of nonincreasing magnitude, then qrsolv expects
//     the full upper triangle of r, the permutation matrix p,
//     and the first n components of (q transpose)*b. the system
//     a*x = b, d*x = 0, is then equivalent to
//
//                  t       t
//           r*z = q *b ,  p *d*p*z = 0 ,
//
//     where x = p*z. if this system does not have full rank,
//     then a least squares solution is obtained. on output qrsolv
//     also provides an upper triangular matrix s such that
//
//            t   t               t
//           p *(a *a + d*d)*p = s *s .
//
//     s is computed within qrsolv and may be of separate interest.
//
//     the subroutine statement is
//
//       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
//
//     where
//
//       n is a positive int input variable set to the order of r.
//
//       r is an n by n array. on input the full upper triangle
//         must contain the full upper triangle of the matrix r.
//         on output the full upper triangle is unaltered, and the
//         strict lower triangle contains the strict upper triangle
//         (transposed) of the upper triangular matrix s.
//
//       ldr is a positive int input variable not less than n
//         which specifies the leading dimension of the array r.
//
//       ipvt is an int input array of length n which defines the
//         permutation matrix p such that a*p = q*r. column j of p
//         is column ipvt[j] of the identity matrix.
//
//       diag is an input array of length n which must contain the
//         diagonal elements of the matrix d.
//
//       qtb is an input array of length n which must contain the first
//         n elements of the vector (q transpose)*b.
//
//       x is an output array of length n which contains the least
//         squares solution of the system a*x = b, d*x = 0.
//
//       sdiag is an output array of length n which contains the
//         diagonal elements of the upper triangular matrix s.
//
//       wa is a work array of length n.
//
//     subprograms called
//
//       fortran-supplied ... fabs,sqrt
//
//     argonne national laboratory. minpack project. march 1980.
//     burton s. garbow, kenneth e. hillstrom, jorge j. more
//
//     **********

//      int ipvt(n)
//      double r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)


void qrsolv(int n, double *r, int ldr, const int *ipvt, const double *diag,
            const double *qtb, double *x, double *sdiag, double *wa)
{
 int i,j,k,l,nsing;
 double cos,qtbpj,sin,sum,temp;
 double p5=0.5, p25=0.25;

 for (j = 0; j < n; ++j) {
   for (i = j; i < n; ++i) {
      r[i + j * ldr] = r[j + i * ldr];}
   x[j] = r[j + j * ldr];
   wa[j] = qtb[j];}

//     eliminate the diagonal matrix d using a givens rotation. 

 for (j = 0; j < n; ++j) {

//        prepare the row of d to be eliminated, locating the 
//        diagonal element using p from the qr factorization. 

    l = ipvt[j]-1;
    if (diag[l] != 0.) {
       for (k = j; k < n; ++k) {
          sdiag[k] = 0.;}
       sdiag[j] = diag[l];

//        the transformations to eliminate the row of d 
//        modify only a single element of (q transpose)*b 
//        beyond the first n, which is initially zero. 

       qtbpj = 0.;
       for (k = j; k < n; ++k) {

//           determine a givens rotation which eliminates the 
//           appropriate element in the current row of d. 

          if (sdiag[k] != 0.) {
             if (fabs(r[k + k * ldr]) < fabs(sdiag[k])) {
                double cotan;
                cotan = r[k + k * ldr] / sdiag[k];
                sin = p5 / sqrt(p25 + p25 * (cotan * cotan));
                cos = sin * cotan;} 
             else {
                double tan;
                tan = sdiag[k] / r[k + k * ldr];
                cos = p5 / sqrt(p25 + p25 * (tan * tan));
                sin = cos * tan;}

//           compute the modified diagonal element of r and 
//           the modified element of ((q transpose)*b,0). 

             temp = cos * wa[k] + sin * qtbpj;
             qtbpj = -sin * wa[k] + cos * qtbpj;
             wa[k] = temp;

//           accumulate the tranformation in the row of s. 

             r[k + k * ldr] = cos * r[k + k * ldr] + sin * sdiag[k];
             if (n > k+1) {
                for (i = k+1; i < n; ++i) {
                   temp = cos * r[i + k * ldr] + sin * sdiag[i];
                   sdiag[i] = -sin * r[i + k * ldr] + cos * sdiag[i];
                   r[i + k * ldr] = temp;}}
             }
          }
      }

//        store the diagonal element of s and restore 
//        the corresponding diagonal element of r. 

   sdiag[j] = r[j + j * ldr];
   r[j + j * ldr] = x[j];}

//     solve the triangular system for z. if the system is 
//     singular, then obtain a least squares solution. 

   nsing = n;
   for (j = 0; j < n; ++j) {
      if (sdiag[j] == 0. && nsing == n) {
         nsing = j;}
      if (nsing < n) {
         wa[j] = 0.;}}
   if (nsing >= 1) {
      for (k = 1; k <= nsing; ++k) {
         j = nsing - k;
         sum = 0.;
         if (nsing > j+1) {
             for (i = j+1; i < nsing; ++i) {
                 sum += r[i + j * ldr] * wa[i];}}
         wa[j] = (wa[j] - sum) / sdiag[j];}}

//    permute the components of z back to components of x. 

 for (j = 0; j < n; ++j) {
     l = ipvt[j]-1;
     x[l] = wa[j];}
 return;

} 


// ****************************************************************** 
}

