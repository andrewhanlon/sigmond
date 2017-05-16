#ifndef TASK_UTILS_H
#define TASK_UTILS_H

#include <iostream>
#include <string>
#include "xml_handler.h"
#include "scalar_defs.h"
#include "matrix.h"
#include "mcobs_handler.h"
#include "chisq_base.h"
#include "minimizer.h"
#include "correlator_matrix_info.h"
#include "array.h"
#include "log_helper.h"

// *******************************************************************
// *                                                                 *
// *   Various important utility routines and classes are defined    *
// *   in this file, mostly related to linear algebra.  Matrix       *
// *   diagonalizations, Cholesky decompositions, and vector         *
// *   pinnings are included.                                        *
// *                                                                 *
// *******************************************************************

 
bool read_arg_type(XMLHandler& xmlin, ComplexArg& arg);


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
   //   within +/- xacc.  "funcd" is a user-supplied functor that
   //   returns both the function value and its derivative.  

double rtsafe(FuncAndDerivSingleVar& funcd, double x1, double x2, double xacc, 
              unsigned int maxit);


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

   //  Reads temporal correlator data and returns vector time separations 
   //  that have all Monte Carlo measurements available

void getCorrelatorAvailableTimes(MCObsHandler *moh,
                                 std::set<uint>& timesavailable,
                                 const CorrelatorInfo& corr, bool hermitian,
                                 ComplexArg arg);

 // ******************************************************************

   //  Reads and computes correlator estimates, returning estimates
   //  in the "results" map.  The integer key in this map is the
   //  time separation.  Note that results are also "put" into 
   //  the memory of the MObsHandler pointed to by "moh".

void getCorrelatorEstimates(MCObsHandler *moh, const CorrelatorInfo& corr, 
                            bool hermitian, bool subtract_vev, ComplexArg arg, 
                            SamplingMode mode, std::map<int,MCEstimate>& results);

 // ******************************************************************

   //  Reads and computes complete correlator matrix estimates for
   //  the current sampling for a particular time separation "timeval".
   //  Results are put into "cormat_estimates".  If the bins of each
   //  element of the correlation matrix are not already in memory,
   //  they are read from file and put into memory, including any
   //  VEVs, if necessary.  Then the current sampling can be computed.
   //  The bins and the current sampling estimate are **left** in 
   //  memory, as well as the current sampling being returned in 
   //  "cormat_estimates".  Hence, a later "eraseHermCorrelatorMatrixAtTime"
   //  is needed to free up memory.  The bins are left in memory in
   //  case other resamplings need to be computed.  The "erase"
   //  routines delete both the bins and all resamplings.


#ifdef COMPLEXNUMBERS

void getHermCorrelatorMatrixAtTime_CurrentSampling(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat, uint timeval,
                  ComplexHermitianMatrix& cormat_estimates);

void getHermCorrelatorMatrixVEVs_CurrentSampling(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat,
                  CVector& vev_estimates);

#else 

void getHermCorrelatorMatrixAtTime_CurrentSampling(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat, uint timeval,
                  RealSymmetricMatrix& cormat_estimates);

void getHermCorrelatorMatrixVEVs_CurrentSampling(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat,
                  RVector& vev_estimates);

#endif

         // routine below does not erase the VEVs
void eraseHermCorrelatorMatrixAtTime(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat, uint timeval);

void eraseHermCorrelatorMatrixVEVs(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat);


    //  Routine below gets the diagonal elements as MCEstimates

void getDiagonalCorrelatorsAtTimeEstimates(MCObsHandler *moh, 
                  const CorrelatorMatrixInfo& cormat, uint timeval,
                  std::vector<MCEstimate>& corrdiag_estimates);


 // ******************************************************************


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
   //  before the effective energy is calculated.


void getEffectiveEnergy(MCObsHandler *moh, const CorrelatorInfo& corr, 
                        bool hermitian, bool subtract_vev, ComplexArg arg, 
                        SamplingMode mode, uint step, uint efftype, 
                        std::map<int,MCEstimate>& results,
                        double subtract_const=0.0);


// ***************************************************************************************

   //  Takes a Hermitian matrix "H" and returns the eigenvalues in
   //  ascending order in "eigvals" and the associated eigenvectors
   //  in the columns of "eigvecs".  Throws an exception if fails.
   //  Versions for only the eigenvalues are also available.

class Diagonalizer
{
 public:

    Diagonalizer(){}
    ~Diagonalizer(){}

    void getEigenvalues(const RealSymmetricMatrix& H, RVector& eigvals);
    void getEigenvectors(const RealSymmetricMatrix& H, RVector& eigvals,
                         RMatrix& eigvecs);
    void getEigenvalues(const ComplexHermitianMatrix& H, RVector& eigvals);
    void getEigenvectors(const ComplexHermitianMatrix& H, 
                         RVector& eigvals, CMatrix& eigvecs);

 private:

    void diagonalize(const RealSymmetricMatrix& H, RVector& eigvals,
                     RMatrix& eigvecs, bool calceigvecs);
    void diagonalize(const ComplexHermitianMatrix& H, 
                     RVector& eigvals, CMatrix& eigvecs, bool calceigvecs);
};


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
   //      DG.isNullMetricInNullMatrix();
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

   //  For each discarded (null space) eigenvector |d> of B, we want to 
   //  check that |d> can be written as a linear superposition of 
   //  the null space eigenvectors |a> of A.  In other words, that
   //                sum_a <d|a><a|d> ~ 1.0


class HermDiagonalizerWithMetric
{
    double mininvcondnum;
    std::vector<double> matb, matg;
    RVector Beigvals, Geigvals;
    int n, n0, np;
    bool xon, Bset, Aset, nullB_in_nullA;
    double negeigalarm;

 public:

    HermDiagonalizerWithMetric();

    HermDiagonalizerWithMetric(double min_inv_cond_num,
                               double negative_eigval_alarm);

    HermDiagonalizerWithMetric(double min_inv_cond_num);

    ~HermDiagonalizerWithMetric();

    void clear();

    void clearMatrix();

    void setMinInvCondNum(double min_inv_cond_num);

    void setNegativeEigenvalueAlarm(double negative_eigval_alarm);

    double getMinInvCondNum() const {return mininvcondnum;}

    double getNegativeEigenvalueAlarm() const {return negeigalarm;}

    void setExceptionsOn();

    void setExceptionsOff();



    int setMetric(const ComplexHermitianMatrix& B, LogHelper& xmlout);

    int setMetric(const ComplexHermitianMatrix& B);

    bool isMetricSet() const {return Bset;}

    void getMetricEigenvalues(RVector& metric_eigvals);

    int getMetricRank() const {return n0;}



    int setMatrix(const ComplexHermitianMatrix& A, LogHelper& xmlout, 
                  bool checkNullSpace=false);
 
    int setMatrix(const ComplexHermitianMatrix& A, bool checkNullSpace=false);

    bool isMatrixSet() const {return Aset;}

    int getMatrixRank() const {return np;}

    bool isNullMetricInNullMatrix() const {return nullB_in_nullA;}

    void getEigenvalues(RVector& eigvals);

    void getEigenvectors(CMatrix& eigvecs);

    void getOrthovectors(CMatrix& orthovecs);

    void getZMatrix(CMatrix& Zmat);

 private:

    static std::string getRotateMetricCode(int info);
    static std::string getRotateMatrixCode(int info);
    friend class RealSymDiagonalizerWithMetric;
    friend class SinglePivotOfCorrMat;
    friend class RollingPivotOfCorrMat;

};



class RealSymDiagonalizerWithMetric
{
    double mininvcondnum;
    std::vector<double> matb, matg;
    RVector Beigvals, Geigvals;
    int n, n0, np;
    bool xon, Bset, Aset, nullB_in_nullA;
    double negeigalarm;

 public:

    RealSymDiagonalizerWithMetric();

    RealSymDiagonalizerWithMetric(double min_inv_cond_num,
                                  double negative_eigval_alarm);

    RealSymDiagonalizerWithMetric(double min_inv_cond_num);

    ~RealSymDiagonalizerWithMetric();

    void clear();

    void clearMatrix();

    void setMinInvCondNum(double min_inv_cond_num);

    void setNegativeEigenvalueAlarm(double negative_eigval_alarm);

    double getMinInvCondNum() const {return mininvcondnum;}

    double getNegativeEigenvalueAlarm() const {return negeigalarm;}

    void setExceptionsOn();

    void setExceptionsOff();



    int setMetric(const RealSymmetricMatrix& B, LogHelper& xmlout);

    int setMetric(const RealSymmetricMatrix& B);

    bool isMetricSet() const {return Bset;}

    void getMetricEigenvalues(RVector& metric_eigvals);

    int getMetricRank() const {return n0;}



    int setMatrix(const RealSymmetricMatrix& A, LogHelper& xmlout,
                  bool checkNullSpace=false);
 
    int setMatrix(const RealSymmetricMatrix& A, bool checkNullSpace=false);

    bool isMatrixSet() const {return Aset;}

    int getMatrixRank() const {return np;}

    bool isNullMetricInNullMatrix() const {return nullB_in_nullA;}

    void getEigenvalues(RVector& eigvals);

    void getEigenvectors(RMatrix& eigvecs);

    void getOrthovectors(RMatrix& orthovecs);

    void getZMatrix(RMatrix& Zmat);

 private:

    static std::string getRotateMetricCode(int info);
    static std::string getRotateMatrixCode(int info);
    friend class SinglePivotOfCorrMat;
    friend class RollingPivotOfCorrMat;

};




// **************************************************************

   //  Given a positive definite real symmetric matrix "A",
   //  an object of this class can construct its Cholesky decomposition
   //                 A = L * transpose(L)
   //  where L is lower triangular, or the Cholesky decomposition of the
   //  inverse of A:
   //                 A^(-1) = transpose(L) * L
   //  where L is lower triangular.  Throws an exception if not 
   //  successful.

class CholeskyDecomposer
{
 public:

    CholeskyDecomposer(){}
    ~CholeskyDecomposer(){}

    void getCholesky(const RealSymmetricMatrix& A, 
                     LowerTriangularMatrix<double>& L);
    void getCholeskyOfInverse(const RealSymmetricMatrix& A, 
                              LowerTriangularMatrix<double>& L);
};




// **************************************************************


   //   An object of class "VectorPinner" stores an indexed list of
   //   reference vectors, which may or may not be mutually orthogonal.
   //   When given a bunch of vectors, it "pins" them to the reference 
   //   vectors.  In other words, it establishes a mapping of the vectors
   //   onto the reference vectors based on maximal inner products (overlaps).
   //   The routine takes each vector and finds which reference vector has
   //   the maximum overlap with it.  If any maximum overlap is smaller than
   //   "warn_fraction", a warning is returned.  If repeated references
   //   are not allowed, then each vector, in the order presented to
   //   the routine, is compared only to the reference vectors not already
   //   matched.  If repeated references are allowed, the routine returns
   //   whether or not a repeated reference has occurred.


template <typename T>
class VectorPinner
{
    std::vector<Vector<T> > m_ref_vecs;
    double m_warn_fraction;
    uint m_numrefs;
    uint m_veclength;
    bool m_repeats;

#ifndef NO_CXX11
    VectorPinner(const VectorPinner& copy) = delete;
    VectorPinner& operator=(const VectorPinner& copy) = delete;
#else
    VectorPinner(const VectorPinner& copy);
    VectorPinner& operator=(const VectorPinner& copy);
#endif

 public:

    VectorPinner(double warn_frac=0.7);
    VectorPinner(const std::vector<Vector<T> >& ref_vecs, double warn_frac=0.7);
    VectorPinner(const Matrix<T>& ref_columns, double warn_frac=0.7);
    ~VectorPinner(){}

    void addReferenceVector(const Vector<T>& ref_vec);
    void resetReferenceVectors(const Matrix<T>& new_ref_columns);

    void setWarningFraction(double warn_frac);
    void setOffRepeatedPinnings() { m_repeats=false;}
    void setOnRepeatedPinnings() { m_repeats=true;}

    double getWarningFraction() const {return m_warn_fraction;}
    unsigned int getNumberRefVectors() const {return m_numrefs;}
    unsigned int getVectorLengths() const {return m_veclength;}
    bool areReferencesOrthogonal() const;
    bool areRepeatsAllowed() const {return m_repeats;}

    void getPinnings(const std::vector<Vector<T> >& vecs,
                     std::vector<uint>& pinnings, bool& repeat_occurred,
                     unsigned int& warnings);
    void getPinnings(const Matrix<T>& col_vecs,
                     std::vector<uint>& pinnings, bool& repeat_occurred,
                     unsigned int& warnings);

 private:

    T dot_prod(const T *avec, const T *bvec) const;
    void do_pin(const T *avec, std::list<uint>& ref_inds,
                uint& ref_pin, double& overlap, bool update) const;
    bool check_for_repeats(const std::vector<uint>& pinnings) const;

};


template <typename T>
VectorPinner<T>::VectorPinner(double warn_frac) 
       : m_warn_fraction(warn_frac), m_numrefs(0), m_veclength(0),
         m_repeats(true)
{
 if ((warn_frac<=0.0)||(warn_frac>=1.0))
    throw(std::invalid_argument("Invalid warning fraction in VectorPinner"));
}


template <typename T>
VectorPinner<T>::VectorPinner(const std::vector<Vector<T> >& ref_vecs, 
                              double warn_frac)
    : m_ref_vecs(ref_vecs), m_warn_fraction(warn_frac), 
      m_numrefs(ref_vecs.size()), m_repeats(true)
{
 if ((warn_frac<=0.0)||(warn_frac>=1.0))
    throw(std::invalid_argument("Invalid warning fraction in VectorPinner"));
 m_veclength=(m_numrefs>0)?ref_vecs[0].size():0;
 for (uint v=1;v<m_numrefs;++v)
    if (ref_vecs[v].size()!=m_veclength)
       throw(std::invalid_argument("Reference vectors must all be same length in VectorPinner"));
 for (uint v=0;v<m_numrefs;++v){
    double rescale=1.0/sqrt(std::abs(dot_prod(&m_ref_vecs[v][0],
                                              &m_ref_vecs[v][0])));
    for (uint k=0;k<m_veclength;++k)
       m_ref_vecs[v][k]*=rescale;}
}


template <typename T>
VectorPinner<T>::VectorPinner(const Matrix<T>& ref_columns, 
                              double warn_frac)
    : m_ref_vecs(ref_columns.size(1)), m_warn_fraction(warn_frac), 
      m_numrefs(ref_columns.size(1)), m_veclength(ref_columns.size(0)), 
      m_repeats(true)
{
 if ((warn_frac<=0.0)||(warn_frac>=1.0))
    throw(std::invalid_argument("Invalid warning fraction in VectorPinner"));
 for (uint v=0;v<m_numrefs;++v){
    m_ref_vecs[v].resize(m_veclength);
    for (uint k=0;k<m_veclength;++k)
       m_ref_vecs[v][k]=ref_columns(k,v);
    double rescale=1.0/sqrt(std::abs(dot_prod(&m_ref_vecs[v][0],
                                              &m_ref_vecs[v][0])));
    for (uint k=0;k<m_veclength;++k)
       m_ref_vecs[v][k]*=rescale;}
}


template <typename T>
void VectorPinner<T>::addReferenceVector(const Vector<T>& ref_vec)
{
 if (ref_vec.size()==0) return;
 if (m_numrefs==0){
    m_veclength=ref_vec.size();}
 else{
    if (ref_vec.size()!=m_veclength)
       throw(std::invalid_argument("Reference vectors must all be same length in VectorPinner"));}
 double rescale=1.0/sqrt(std::abs(dot_prod(&ref_vec[0],&ref_vec[0])));
 m_ref_vecs.push_back(ref_vec);
 for (uint k=0;k<m_veclength;++k)
    m_ref_vecs[m_numrefs][k]*=rescale;
 m_numrefs++;
}


template <typename T>
void VectorPinner<T>::resetReferenceVectors(const Matrix<T>& new_ref_columns)
{
 if ((new_ref_columns.size(1)!=m_numrefs)||(new_ref_columns.size(0)!=m_veclength)){
    throw(std::invalid_argument("Mismatch in updating references in VectorPinner"));}
 for (uint v=0;v<m_numrefs;++v){
    m_ref_vecs[v].resize(m_veclength);
    for (uint k=0;k<m_veclength;++k)
       m_ref_vecs[v][k]=new_ref_columns(k,v);
    double rescale=1.0/sqrt(std::abs(dot_prod(&m_ref_vecs[v][0],
                                              &m_ref_vecs[v][0])));
    for (uint k=0;k<m_veclength;++k)
       m_ref_vecs[v][k]*=rescale;}
}


template <typename T>
void VectorPinner<T>::setWarningFraction(double warn_frac)
{
 if ((warn_frac<=0.0)||(warn_frac>=1.0))
    throw(std::invalid_argument("Invalid warning fraction in VectorPinner"));
 m_warn_fraction=warn_frac;
}


template <typename T>
bool VectorPinner<T>::areReferencesOrthogonal() const
{
 for (uint j=0;j<m_numrefs;++j)
 for (uint k=0;k<j;++k){
    if (std::abs(dot_prod(&m_ref_vecs[j][0],&m_ref_vecs[k][0]))>1e-10){
       return false;}}
 return true;
}



template <typename T>
void VectorPinner<T>::getPinnings(const std::vector<Vector<T> >& vecs,
                                  std::vector<uint>& pinnings, bool& repeat_occurred,
                                  unsigned int& warnings)
{
 if (vecs.size()>m_numrefs)
    throw(std::invalid_argument("Number of vectors exceeds references in VectorPinner"));
 std::list<uint> ref_inds;
 for (uint v=0;v<m_numrefs;++v) ref_inds.push_back(v);
 warnings=0;
 repeat_occurred=false;
 pinnings.resize(vecs.size());
 double overlap;
 for (uint k=0;k<vecs.size();++k){
    if (vecs[k].size()!=m_veclength)
       throw(std::invalid_argument("vector has incorrect length in VectorPinning::getPinnings"));
    do_pin(&vecs[k][0],ref_inds,pinnings[k],overlap,!m_repeats);
    if (overlap<m_warn_fraction) warnings++;}
 if (m_repeats)
    repeat_occurred=check_for_repeats(pinnings);
}


template <typename T>
void VectorPinner<T>::getPinnings(const Matrix<T>& mat, std::vector<uint>& pinnings, 
                                  bool& repeat_occurred, unsigned int& warnings)
{
 if (mat.size(1)>m_numrefs)
    throw(std::invalid_argument("Number of vectors exceeds references in VectorPinner"));
 if (mat.size(0)!=m_veclength)
    throw(std::invalid_argument("vectors have incorrect length in VectorPinning::getPinnings"));
 std::list<uint> ref_inds;
 for (uint v=0;v<m_numrefs;++v) ref_inds.push_back(v);
 warnings=0;
 repeat_occurred=false;
 pinnings.resize(mat.size(1));
 double overlap;
 for (uint k=0;k<mat.size(1);++k){
    do_pin(&mat.get(0,k),ref_inds,pinnings[k],overlap,!m_repeats);
    if (overlap<m_warn_fraction) warnings++;}
 if (m_repeats)
    repeat_occurred=check_for_repeats(pinnings);
}


template <typename T>
void VectorPinner<T>::do_pin(const T *avec, std::list<uint>& ref_inds,
                             uint& ref_pin, double& overlap, bool update) const
{
 std::list<uint>::iterator pin=ref_inds.end();
 overlap=0.0;
 double a_rescale=1.0/sqrt(std::abs(dot_prod(avec,avec)));
 for (std::list<uint>::iterator it=ref_inds.begin();it!=ref_inds.end();++it){
    double cv=std::abs(dot_prod(&m_ref_vecs[*it][0],avec))*a_rescale;
    if (cv>overlap){
       pin=it; overlap=cv;}}
 if (pin==ref_inds.end())
    throw(std::invalid_argument("No suitable overlap found in VectorPinner"));
 ref_pin=*pin;
 if (update)
    ref_inds.erase(pin);
}


template <typename T>
bool VectorPinner<T>::check_for_repeats(const std::vector<uint>& pinnings) const
{
 for (uint j=0;j<pinnings.size();++j)
 for (uint k=0;k<j;++k)
    if (pinnings[j]==pinnings[k]) return true;
 return false;
}

// *************************************************************

    //   Rescales the matrix "cormat" using the diagonal elements
    //   of the matrix "mat_scales" according to
    //
    //     cormat(i,j) / sqrt( |mat_scales(i,i)|*mat_scales(j,j)| )


void doRescaleByDiagonals(ComplexHermitianMatrix& cormat,
                          const ComplexHermitianMatrix& mat_scales);

void doRescaleByDiagonals(RealSymmetricMatrix& cormat,
                          const RealSymmetricMatrix& mat_scales);


// *************************************************************

    //   Rescales the transformation matrix "R" using the 
    //   diagonal elements of the matrix "mat_scales" according to
    //
    //     R(i,j) / sqrt( |mat_scales(i,i)| )


void doRescaleTransformation(CMatrix& R,
              const ComplexHermitianMatrix& mat_scales);

void doRescaleTransformation(RMatrix& R,
              const RealSymmetricMatrix& mat_scales);


// *************************************************************

     //  Takes Hermitian matrix "A" and replaces it with the rotated
     //  matrix   R^dagger A  R.  This is also a version that just
     //  evaluates the diagonal elements of the rotated matrix.

void doMatrixRotation(ComplexHermitianMatrix& A, const CMatrix& R);

void doMatrixRotation(const ComplexHermitianMatrix& A, const CMatrix& R,
                      RVector& Ardiag);


void doMatrixRotation(RealSymmetricMatrix& A, const RMatrix& R);

void doMatrixRotation(const RealSymmetricMatrix& A, const RMatrix& R,
                      RVector& Ardiag);


// ********************************************************************

     //  Takes vector "V" and replaces it with the rotated
     //  matrix   R^dagger V.  


void doVectorRotation(CVector& V, const CMatrix& R);

void doVectorRotation(RVector& V, const RMatrix& R);


// ********************************************************************

     //   outvec = inmat * invec

void multiply(CVector& outvec, const ComplexHermitianMatrix& inmat, const CVector& invec);
void multiply(RVector& outvec, const RealSymmetricMatrix& inmat, const RVector& invec);

     //   inner product of vectors

std::complex<double> dotProduct(const CVector& lvec, const CVector& rvec);
double dotProduct(const RVector& lvec, const RVector& rvec);

     //   magnitude of inner product of vectors

double dotProductMagnitude(const CVector& lvec, const CVector& rvec);
double dotProductMagnitude(const RVector& lvec, const RVector& rvec);

     //   magnitude squared of inner product of vectors

double dotProductMagnitudeSquared(const CVector& lvec, const CVector& rvec);
double dotProductMagnitudeSquared(const RVector& lvec, const RVector& rvec);


// ********************************************************************

void array_to_matrix(const Array<std::complex<double> >& in, CMatrix& out);
void array_to_matrix(const Array<std::complex<float> >& in, CMatrix& out);
void array_to_matrix(const Array<double>& in, RMatrix& out);
void array_to_matrix(const Array<float>& in, RMatrix& out);
void array_to_matrix(const Array<double>& in, CMatrix& out);
void array_to_matrix(const Array<float>& in, CMatrix& out);

void matrix_to_array(const CMatrix& in, Array<std::complex<double> >& out);
void matrix_to_array(const CMatrix& in, Array<std::complex<float> >& out);
void matrix_to_array(const RMatrix& in, Array<double>& out);
void matrix_to_array(const RMatrix& in, Array<float>& out);
void matrix_to_array(const CMatrix& in, Array<double>& out);
void matrix_to_array(const CMatrix& in, Array<float>& out);

void array_to_vector(const Array<double>& in, std::vector<double>& out);
void array_to_RVector(const Array<double>& in, RVector& out);

void vector_to_array(const std::vector<double>& in, Array<double>& out);
void RVector_to_array(const RVector& in, Array<double>& out);

std::vector<uint> form_tvalues(uint tmin, uint tmax, 
                               const std::vector<int>& texclude);

// ********************************************************************

void doSquareByBins(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out);

void doSquareBySamplings(MCObsHandler& moh, const MCObsInfo& obs_in, const MCObsInfo& obs_out);

void doRatioByBins(MCObsHandler& moh, const MCObsInfo& obs_numer, const MCObsInfo& obs_denom,
                   const MCObsInfo& obs_ratio);

void doRatioBySamplings(MCObsHandler& moh, const MCObsInfo& obs_numer, 
                        const MCObsInfo& obs_denom, const MCObsInfo& obs_ratio);

void doLinearSuperpositionByBins(MCObsHandler& moh, std::vector<MCObsInfo>& suminfos,
                   std::vector<double>& sumcoefs, const MCObsInfo& obs_superposition);

void doLinearSuperpositionBySamplings(MCObsHandler& moh, std::vector<MCObsInfo>& suminfos,
                   std::vector<double>& sumcoefs, const MCObsInfo& obs_superposition);

            // evaluates energy_squared = rest_mass_squared + psqfactor / (xi*xi)
            //    where psqfactor = (2*Pi/L)^2*nsq

void doDispersionBySamplings(MCObsHandler& moh, const MCObsInfo& anisotropy_key, 
                             const MCObsInfo& restmasssquared_key,  double psqfactor,
                             const MCObsInfo& Esqinfo);

// ********************************************************************

inline bool level_compare(const std::pair<double,uint>& a, const std::pair<double,uint>& b)
{
 return (a.first<b.first);
}

// ********************************************************************
#endif
