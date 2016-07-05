#include "minpack.h"
#include <iostream>
#include "chisq_base.h"

using namespace std;


    // Nonlinear least-squares minimization code from 
    // MINPACK-1, translated into C++ and adapted for
    // lattice field theory analysis uses.  
    // Thankfully, hidden inside the "MinPack" namespace.


//     ***********************************************************
//
//     the purpose of lmder is to minimize the sum of the squares of
//     m nonlinear functions in n variables by a modification of the
//     levenberg-marquardt algorithm. this is done by using the more
//     general least-squares solver lmder. the user must provide a
//     subroutine which calculates the functions and the jacobian.
//
//     the subroutine statement is
//
//       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
//                        ipvt,wa,lwa)
//
//     where
//
//       fcn is the name of the user-supplied subroutine which
//         calculates the functions and the jacobian. fcn must
//         be declared in an external statement in the user
//         calling program, and should be written as follows.
//
//         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
//         int m,n,ldfjac,iflag
//         double x(n),fvec(m),fjac(ldfjac,n)
//         ----------
//         if iflag = 1 calculate the functions at x and
//         return this vector in fvec. do not alter fjac.
//         if iflag = 2 calculate the jacobian at x and
//         return this matrix in fjac. do not alter fvec.
//         ----------
//         return
//         end
//
//         the value of iflag should not be changed by fcn unless
//         the user wants to terminate execution of lmder1.
//         in this case set iflag to a negative int.
//
//       m is a positive int input variable set to the number
//         of functions.
//
//       n is a positive int input variable set to the number
//         of variables. n must not exceed m.
//
//       x is an array of length n. on input x must contain
//         an initial estimate of the solution vector. on output x
//         contains the final estimate of the solution vector.
//
//       fvec is an output array of length m which contains
//         the functions evaluated at the output x.
//
//       fjac is an output m by n array. the upper n by n submatrix
//         of fjac contains an upper triangular matrix r with
//         diagonal elements of nonincreasing magnitude such that
//
//                t     t           t
//               p *(jac *jac)*p = r *r,
//
//         where p is a permutation matrix and jac is the final
//         calculated jacobian. column j of p is column ipvt[j]
//         (see below) of the identity matrix. the lower trapezoidal
//         part of fjac contains information generated during
//         the computation of r.
//
//       ldfjac is a positive int input variable not less than m
//         which specifies the leading dimension of the array fjac.
//
//       tol is a nonnegative input variable. termination occurs
//         when the algorithm estimates either that the relative
//         error in the sum of squares is at most tol or that
//         the relative error between x and the solution is at
//         most tol.
//
//       info is an int output variable. if the user has
//         terminated execution, info is set to the (negative)
//         value of iflag. see description of fcn. otherwise,
//         info is set as follows.
//
//         info = 0  improper input parameters.
//
//         info = 1  algorithm estimates that the relative error
//                   in the sum of squares is at most tol.
//
//         info = 2  algorithm estimates that the relative error
//                   between x and the solution is at most tol.
//
//         info = 3  conditions for info = 1 and info = 2 both hold.
//
//         info = 4  fvec is orthogonal to the columns of the
//                   jacobian to machine precision.
//
//         info = 5  number of calls to fcn with iflag = 1 has
//                   reached 100*(n+1).
//
//         info = 6  tol is too small. no further reduction in
//                   the sum of squares is possible.
//
//         info = 7  tol is too small. no further improvement in
//                   the approximate solution x is possible.
//
//       ipvt is an int output array of length n. ipvt
//         defines a permutation matrix p such that jac*p = q*r,
//         where jac is the final calculated jacobian, q is
//         orthogonal (not stored), and r is upper triangular
//         with diagonal elements of nonincreasing magnitude.
//         column j of p is column ipvt[j] of the identity matrix.
//
//       wa is a work array of length lwa.
//
//       lwa is a positive int input variable not less than 5*n+m.
//
//     subprograms called
//
//       user-supplied ...... fcn
//
//       minpack-supplied ... lmder
//
//     argonne national laboratory. minpack project. march 1980.
//     burton s. garbow, kenneth e. hillstrom, jorge j. more
//
//     **********

namespace MinPack {

// ******************************************************************

void lmpar(int n, double *r, int ldr, int *ipvt, double *diag,
           double *qtb, double delta, double *par, double *x,
           double *sdiag, double *wa1, double* wa2);


// *****************************************************************
//
//     subroutine lmder
//
//     the purpose of lmder is to minimize the sum of the squares of
//     m nonlinear functions in n variables by a modification of
//     the levenberg-marquardt algorithm. the user must provide a
//     subroutine which calculates the functions and the jacobian.
//
//     the subroutine statement is
//
//       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
//                        maxfev,diag,mode,factor,nprint,info,nfev,
//                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
//
//     where
//
//       fcn is the name of the user-supplied subroutine which
//         calculates the functions and the jacobian. fcn must
//         be declared in an external statement in the user
//         calling program, and should be written as follows.
//
//         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
//         int m,n,ldfjac,iflag
//         double x(n),fvec(m),fjac(ldfjac,n)
//         ----------
//         if iflag = 1 calculate the functions at x and
//         return this vector in fvec. do not alter fjac.
//         if iflag = 2 calculate the jacobian at x and
//         return this matrix in fjac. do not alter fvec.
//         ----------
//         return
//         end
//
//         the value of iflag should not be changed by fcn unless
//         the user wants to terminate execution of lmder.
//         in this case set iflag to a negative int.
//
//       m is a positive int input variable set to the number
//         of functions.
//
//       n is a positive int input variable set to the number
//         of variables. n must not exceed m.
//
//       x is an array of length n. on input x must contain
//         an initial estimate of the solution vector. on output x
//         contains the final estimate of the solution vector.
//
//       fvec is an output array of length m which contains
//         the functions evaluated at the output x.
//
//       fjac is an output m by n array. the upper n by n submatrix
//         of fjac contains an upper triangular matrix r with
//         diagonal elements of nonincreasing magnitude such that
//
//                t     t           t
//               p *(jac *jac)*p = r *r,
//
//         where p is a permutation matrix and jac is the final
//         calculated jacobian. column j of p is column ipvt[j]
//         (see below) of the identity matrix. the lower trapezoidal
//         part of fjac contains information generated during
//         the computation of r.
//
//       ldfjac is a positive int input variable not less than m
//         which specifies the leading dimension of the array fjac.
//
//       ftol is a nonnegative input variable. termination
//         occurs when both the actual and predicted relative
//         reductions in the sum of squares are at most ftol.
//         therefore, ftol measures the relative error desired
//         in the sum of squares.
//
//       xtol is a nonnegative input variable. termination
//         occurs when the relative error between two consecutive
//         iterates is at most xtol. therefore, xtol measures the
//         relative error desired in the approximate solution.
//
//       gtol is a nonnegative input variable. termination
//         occurs when the cosine of the angle between fvec and
//         any column of the jacobian is at most gtol in absolute
//         value. therefore, gtol measures the orthogonality
//         desired between the function vector and the columns
//         of the jacobian.
//
//       maxfev is a positive int input variable. termination
//         occurs when the number of calls to fcn with iflag = 1
//         has reached maxfev.
//
//       diag is an array of length n. if mode = 1 (see
//         below), diag is internally set. if mode = 2, diag
//         must contain positive entries that serve as
//         multiplicative scale factors for the variables.
//
//       mode is an int input variable. if mode = 1, the
//         variables will be scaled internally. if mode = 2,
//         the scaling is specified by the input diag. other
//         values of mode are equivalent to mode = 1.
//
//       factor is a positive input variable used in determining the
//         initial step bound. this bound is set to the product of
//         factor and the euclidean norm of diag*x if nonzero, or else
//         to factor itself. in most cases factor should lie in the
//         interval (.1,100.).100. is a generally recommended value.
//
//       nprint is an int input variable that enables controlled
//         printing of iterates if it is positive. in this case,
//         fcn is called with iflag = 0 at the beginning of the first
//         iteration and every nprint iterations thereafter and
//         immediately prior to return, with x, fvec, and fjac
//         available for printing. fvec and fjac should not be
//         altered. if nprint is not positive, no special calls
//         of fcn with iflag = 0 are made.
//
//       info is an int output variable. if the user has
//         terminated execution, info is set to the (negative)
//         value of iflag. see description of fcn. otherwise,
//         info is set as follows.
//
//         info = 0  improper input parameters.
//
//         info = 1  both actual and predicted relative reductions
//                   in the sum of squares are at most ftol.
//
//         info = 2  relative error between two consecutive iterates
//                   is at most xtol.
//
//         info = 3  conditions for info = 1 and info = 2 both hold.
//
//         info = 4  the cosine of the angle between fvec and any
//                   column of the jacobian is at most gtol in
//                   absolute value.
//
//         info = 5  number of calls to fcn with iflag = 1 has
//                   reached maxfev.
//
//         info = 6  ftol is too small. no further reduction in
//                   the sum of squares is possible.
//
//         info = 7  xtol is too small. no further improvement in
//                   the approximate solution x is possible.
//
//         info = 8  gtol is too small. fvec is orthogonal to the
//                   columns of the jacobian to machine precision.
//
//       nfev is an int output variable set to the number of
//         calls to fcn with iflag = 1.
//
//       njev is an int output variable set to the number of
//         calls to fcn with iflag = 2.
//
//       ipvt is an int output array of length n. ipvt
//         defines a permutation matrix p such that jac*p = q*r,
//         where jac is the final calculated jacobian, q is
//         orthogonal (not stored), and r is upper triangular
//         with diagonal elements of nonincreasing magnitude.
//         column j of p is column ipvt[j] of the identity matrix.
//
//       qtf is an output array of length n which contains
//         the first n elements of the vector (q transpose)*fvec.
//
//       wa1, wa2, and wa3 are work arrays of length n.
//
//       wa4 is a work array of length m.
//
//     subprograms called
//
//       user-supplied ...... fcn
//
//       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
//
//       fortran-supplied ... fabs,max,min,sqrt,mod
//
//     argonne national laboratory. minpack project. march 1980.
//     burton s. garbow, kenneth e. hillstrom, jorge j. more
//
//     **********

int lmder(ChiSquare *p, int m, int n, vector<double>& params, 
          vector<double>& residuals, RMatrix& gradients, 
          int ldfjac, double ftol, double xtol, double gtol, int maxfev, 
          int mode, double factor, int nprint, int &nfev, 
          int &njev, vector<double>& vdiag, vector<double>& vwa1, 
          vector<double>& vwa2, vector<double>& vwa3, vector<double>& vqtf,
          vector<double> vwa4, vector<int> vipvt, ostringstream& outlog)
{     
 int i, j, l, iter=0, iflag, info;
 double d1, d2, par, sum, temp, temp1, temp2, ratio;
 double delta = 0.0;
 double fnorm=0.0, gnorm, pnorm, xnorm = 0., fnorm1, actred, dirder, 
        epsmch, prered;
 double p1=0.1, p5=0.5, p25=0.25, p75=0.75, p0001=1.0e-4;

 double *x=&params[0];
 double *fvec=&residuals[0];
 double *fjac=&gradients(0,0);
 double *diag=&vdiag[0];
 double *wa1=&vwa1[0];
 double *wa2=&vwa2[0]; 
 double *wa3=&vwa3[0];
 double *qtf=&vqtf[0];
 double *wa4=&vwa4[0];
 int *ipvt=&vipvt[0];

//     epsmch is the machine precision.

 epsmch = dpmpar(1);

 info = 0;
 iflag = 0;
 nfev= 0;
 njev = 0;

    //     check the input parameters for errors. 

 if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. || 
         gtol < 0. || maxfev <= 0 || factor <= 0.) {
     goto TERMINATE;
 }
 if (mode == 2) {
     for (j = 0; j < n; ++j) {
         if (diag[j] <= 0.) {
             goto TERMINATE;
         }
     }
 }

    //     evaluate the function at the starting point 
    //     and calculate its norm. 

 p->evalResiduals(params,residuals);
 nfev= 1;
 fnorm = enorm(m, fvec);

    //     initialize levenberg-marquardt parameter and iteration counter. 

 par = 0.;
 iter = 1;

    //     beginning of the outer loop. 

 for (;;) {

    //        calculate the jacobian matrix. 

     p->evalResGradients(params,gradients);
     ++njev;

    //        if requested, call fcn to enable printing of iterates. 

     if (nprint > 0) {
         iflag = 0;
         if ((iter - 1) % nprint == 0) {
            outlog << " Iteration "<<iter-1<<":   chi_sq = "
                 <<sqr(fnorm)<<endl;
            for (i=0;i<n;i++)
                outlog << "     param["<<i<<"] = "<<x[i]<<endl;
         }
         if (iflag < 0) {
             goto TERMINATE;
         }
     }

    //        compute the qr factorization of the jacobian. 

     qrfac(m, n, fjac, ldfjac, True, ipvt, n, wa1, wa2, wa3);

    //        on the first iteration and if mode is 1, scale according 
    //        to the norms of the columns of the initial jacobian. 

     if (iter == 1) {
         if (mode != 2) {
             for (j = 0; j < n; ++j) {
                 diag[j] = wa2[j];
                 if (wa2[j] == 0.) {
                     diag[j] = 1.;
                 }
             }
         }

    //        on the first iteration, calculate the norm of the scaled x 
    //        and initialize the step bound delta. 

         for (j = 0; j < n; ++j) {
             wa3[j] = diag[j] * x[j];
         }
         xnorm = enorm(n, wa3);
         delta = factor * xnorm;
         if (delta == 0.) {
             delta = factor;
         }
     }

    //        form (q transpose)*fvec and store the first n components in 
    //        qtf. 

     for (i = 0; i < m; ++i) {
         wa4[i] = fvec[i];
     }
     for (j = 0; j < n; ++j) {
         if (fjac[j + j * ldfjac] != 0.) {
             sum = 0.;
             for (i = j; i < m; ++i) {
                 sum += fjac[i + j * ldfjac] * wa4[i];
             }
             temp = -sum / fjac[j + j * ldfjac];
             for (i = j; i < m; ++i) {
                 wa4[i] += fjac[i + j * ldfjac] * temp;
             }
         }
         fjac[j + j * ldfjac] = wa1[j];
         qtf[j] = wa4[j];
     }

    //        compute the norm of the scaled gradient. 

     gnorm = 0.;
     if (fnorm != 0.) {
         for (j = 0; j < n; ++j) {
             l = ipvt[j]-1;
             if (wa2[l] != 0.) {
                 sum = 0.;
                 for (i = 0; i <= j; ++i) {
                     sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
                 }
                        // Computing MAX 
                 d1 = fabs(sum / wa2[l]);
                 gnorm = max(gnorm,d1);
             }
         }
     }

    //        test for convergence of the gradient norm. 

     if (gnorm <= gtol) {
         info = 4;
     }
     if (info != 0) {
         goto TERMINATE;
     }

    //        rescale if necessary. 

     if (mode != 2) {
         for (j = 0; j < n; ++j) {
                    // Computing MAX 
             d1 = diag[j], d2 = wa2[j];
             diag[j] = max(d1,d2);
         }
     }

    //        beginning of the inner loop. 

     do {

    //           determine the levenberg-marquardt parameter. 

         lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta,
               &par, wa1, wa2, wa3, wa4);

    //           store the direction p and x + p. calculate the norm of p. 

         for (j = 0; j < n; ++j) {
             wa1[j] = -wa1[j];
             wa2[j] = x[j] + wa1[j];
             wa3[j] = diag[j] * wa1[j];
         }
         pnorm = enorm(n, wa3);

    //           on the first iteration, adjust the initial step bound. 

         if (iter == 1) {
             delta = min(delta,pnorm);
         }

    //           evaluate the function at x + p and calculate its norm. 

         p->evalResiduals(vwa2,vwa4);
         ++nfev;
         if (iflag < 0) {
             goto TERMINATE;
         }
         fnorm1 = enorm(m, wa4);

    //           compute the scaled actual reduction. 

         actred = -1.;
         if (p1 * fnorm1 < fnorm) {
                    // Computing 2nd power 
             d1 = fnorm1 / fnorm;
             actred = 1. - d1 * d1;
         }

    //           compute the scaled predicted reduction and 
    //           the scaled directional derivative. 

         for (j = 0; j < n; ++j) {
             wa3[j] = 0.;
             l = ipvt[j]-1;
             temp = wa1[l];
             for (i = 0; i <= j; ++i) {
                 wa3[i] += fjac[i + j * ldfjac] * temp;
             }
         }
         temp1 = enorm(n, wa3) / fnorm;
         temp2 = (sqrt(par) * pnorm) / fnorm;
         prered = temp1 * temp1 + temp2 * temp2 / p5;
         dirder = -(temp1 * temp1 + temp2 * temp2);

    //           compute the ratio of the actual to the predicted 
    //           reduction. 

         ratio = 0.;
         if (prered != 0.) {
             ratio = actred / prered;
         }

    //           update the step bound. 

         if (ratio <= p25) {
             if (actred >= 0.) {
                 temp = p5;
             } else {
                 temp = p5 * dirder / (dirder + p5 * actred);
             }
             if (p1 * fnorm1 >= fnorm || temp < p1) {
                 temp = p1;
             }
                    // Computing MIN 
             d1 = pnorm / p1;
             delta = temp * min(delta,d1);
             par /= temp;
         } else {
             if (par == 0. || ratio >= p75) {
                 delta = pnorm / p5;
                 par = p5 * par;
             }
         }

    //           test for successful iteration. 

         if (ratio >= p0001) {

    //           successful iteration. update x, fvec, and their norms. 

             for (j = 0; j < n; ++j) {
                 x[j] = wa2[j];
                 wa2[j] = diag[j] * x[j];
             }
             for (i = 0; i < m; ++i) {
                 fvec[i] = wa4[i];
             }
             xnorm = enorm(n, wa2);
             fnorm = fnorm1;
             ++iter;
         }

    //           tests for convergence. 

      //   if ((fabs(actred) <= ftol) && (prered <= ftol) && (p5 * ratio <= 1.)) {
      //       info = 1;
      //   }
         if (delta <= xtol * xnorm) {
             info = 2;
         }
         if ((fabs(actred) <= ftol) && (prered <= ftol) 
             && (p5 * ratio <= 1.) && (info == 2)) {
             info = 3;
         }
         if (info != 0) {
             goto TERMINATE;
         }

    //           tests for termination and stringent tolerances. 

         if (nfev>= maxfev) {
             info = 5;
         }
         if (fabs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1.) {
             info = 6;
         }
         if (delta <= epsmch * xnorm) {
             info = 7;
         }
         if (gnorm <= epsmch) {
             info = 8;
         }
         if (info != 0) {
             goto TERMINATE;
         }

    //           end of the inner loop. repeat if iteration unsuccessful. 

     } while (ratio < p0001);

    //        end of the outer loop. 

 }
TERMINATE:

       //     termination, either normal or user imposed. 

 if (iflag < 0) {
     info = iflag;
 }
 if (nprint > 0) {
    outlog << " Iteration "<<iter<<" (FINAL):   chi_sq = "
         << sqr(fnorm)<<endl;
    for (i=0;i<n;i++)
       outlog << "     param["<<i<<"] = "<<x[i]<<endl;
 }
 return info;

}


// ******************************************************************
//
//     subroutine lmpar
//
//     given an m by n matrix a, an n by n nonsingular diagonal
//     matrix d, an m-vector b, and a positive number delta,
//     the problem is to determine a value for the parameter
//     par such that if x solves the system
//
//           a*x = b ,     sqrt(par)*d*x = 0 ,
//
//     in the least squares sense, and dxnorm is the euclidean
//     norm of d*x, then either par is zero and
//
//           (dxnorm-delta)<=0.1*delta ,
//
//     or par is positive and
//
//           abs(dxnorm-delta)<=0.1*delta .
//
//     this subroutine completes the solution of the problem
//     if it is provided with the necessary information from the
//     qr factorization, with column pivoting, of a. that is, if
//     a*p = q*r, where p is a permutation matrix, q has orthogonal
//     columns, and r is an upper triangular matrix with diagonal
//     elements of nonincreasing magnitude, then lmpar expects
//     the full upper triangle of r, the permutation matrix p,
//     and the first n components of (q transpose)*b. on output
//     lmpar also provides an upper triangular matrix s such that
//
//            t   t                   t
//           p *(a *a + par*d*d)*p = s *s .
//
//     s is employed within lmpar and may be of separate interest.
//
//     only a few iterations are generally needed for convergence
//     of the algorithm. if, however, the limit of 10 iterations
//     is reached, then the output par will contain the best
//     value obtained so far.
//
//     the subroutine statement is
//
//       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
//                        wa1,wa2)
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
//       delta is a positive input variable which specifies an upper
//         bound on the euclidean norm of d*x.
//
//       par is a nonnegative variable. on input par contains an
//         initial estimate of the levenberg-marquardt parameter.
//         on output par contains the final estimate.
//
//       x is an output array of length n which contains the least
//         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
//         for the output par.
//
//       sdiag is an output array of length n which contains the
//         diagonal elements of the upper triangular matrix s.
//
//       wa1 and wa2 are work arrays of length n.
//
//     subprograms called
//
//       minpack-supplied ... dpmpar,enorm,qrsolv
//
//       fortran-supplied ... fabs,max,min,sqrt
//
//     argonne national laboratory. minpack project. march 1980.
//     burton s. garbow, kenneth e. hillstrom, jorge j. more
//
//     **********


void lmpar(int n, double *r, int ldr, int *ipvt, double *diag,
           double *qtb, double delta, double *par, double *x,
           double *sdiag, double *wa1, double* wa2)
{
 int iter,j,l,nsing;
 double dxnorm,dwarf,fp,gnorm,parc,parl,paru,temp,d1,d2;
 double p1=0.1, p001=0.001;

//     dwarf is the smallest positive magnitude.

 dwarf = dpmpar(2);

 nsing = n;
 for (j = 0; j < n; ++j) {
     wa1[j] = qtb[j];
     if (r[j + j * ldr] == 0. && nsing == n) {
         nsing = j;
     }
     if (nsing < n) {
         wa1[j] = 0.;
     }
 }
 if (nsing >= 1) {
     int k;
     for (k = 1; k <= nsing; ++k) {
         j = nsing - k;
         wa1[j] /= r[j + j * ldr];
         temp = wa1[j];
         if (j >= 1) {
             int i;
             for (i = 0; i < j; ++i) {
                 wa1[i] -= r[i + j * ldr] * temp;
             }
         }
     }
 }
 for (j = 0; j < n; ++j) {
     l = ipvt[j]-1;
     x[l] = wa1[j];
 }

//     initialize the iteration counter. 
//     evaluate the function at the origin, and test 
//     for acceptance of the gauss-newton direction. 

 iter = 0;
 for (j = 0; j < n; ++j) {
     wa2[j] = diag[j] * x[j];
 }
 dxnorm = enorm(n, wa2);
 fp = dxnorm - delta;
 if (fp <= p1 * delta) {
     goto TERMINATE;
 }

//     if the jacobian is not rank deficient, the newton 
//     step provides a lower bound, parl, for the zero of 
//     the function. otherwise set this bound to zero. 

 parl = 0.;
 if (nsing >= n) {
     for (j = 0; j < n; ++j) {
         l = ipvt[j]-1;
         wa1[j] = diag[l] * (wa2[l] / dxnorm);
     }
     for (j = 0; j < n; ++j) {
         double sum = 0.;
         if (j >= 1) {
             int i;
             for (i = 0; i < j; ++i) {
                 sum += r[i + j * ldr] * wa1[i];
             }
         }
         wa1[j] = (wa1[j] - sum) / r[j + j * ldr];
     }
     temp = enorm(n, wa1);
     parl = fp / delta / temp / temp;
 }

//     calculate an upper bound, paru, for the zero of the function. 

 for (j = 0; j < n; ++j) {
     double sum;
     int i;
     sum = 0.;
     for (i = 0; i <= j; ++i) {
         sum += r[i + j * ldr] * qtb[i];
     }
     l = ipvt[j]-1;
     wa1[j] = sum / diag[l];
 }
 gnorm = enorm(n, wa1);
 paru = gnorm / delta;
 if (paru == 0.) {
     paru = dwarf / min(delta,double(p1));
 }

//     if the input par lies outside of the interval (parl,paru), 
//     set par to the closer endpoint. 

 *par = max(*par,parl);
 *par = min(*par,paru);
 if (*par == 0.) {
     *par = gnorm / dxnorm;
 }

//     beginning of an iteration. 

 for (;;) {
     ++iter;

//        evaluate the function at the current value of par. 

     if (*par == 0.) {
         // Computing MAX 
         d1 = dwarf, d2 = p001 * paru;
         *par = max(d1,d2);
     }
     temp = sqrt(*par);
     for (j = 0; j < n; ++j) {
         wa1[j] = temp * diag[j];
     }
     qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
     for (j = 0; j < n; ++j) {
         wa2[j] = diag[j] * x[j];
     }
     dxnorm = enorm(n, wa2);
     temp = fp;
     fp = dxnorm - delta;

//        if the function is small enough, accept the current value 
//        of par. also test for the exceptional cases where parl 
//        is zero or the number of iterations has reached 10. 

     if (fabs(fp) <= p1 * delta || (parl == 0. && fp <= temp && temp < 0.) || iter == 10) {
         goto TERMINATE;
     }

//        compute the newton correction. 

     for (j = 0; j < n; ++j) {
         l = ipvt[j]-1;
         wa1[j] = diag[l] * (wa2[l] / dxnorm);
     }
     for (j = 0; j < n; ++j) {
         wa1[j] /= sdiag[j];
         temp = wa1[j];
         if (n > j+1) {
             int i;
             for (i = j+1; i < n; ++i) {
                 wa1[i] -= r[i + j * ldr] * temp;
             }
         }
     }
     temp = enorm(n, wa1);
     parc = fp / delta / temp / temp;

//        depending on the sign of the function, update parl or paru. 

     if (fp > 0.) {
         parl = max(parl,*par);
     }
     if (fp < 0.) {
         paru = min(paru,*par);
     }

//        compute an improved estimate for par. 

        // Computing MAX 
     d1 = parl, d2 = *par + parc;
     *par = max(d1,d2);

//        end of an iteration. 

 }
TERMINATE:

//     termination. 

 if (iter == 0) {
     *par = 0.;
 }


} 


// ******************************************************************
}

