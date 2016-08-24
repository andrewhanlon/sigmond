#ifndef SINGLE_PIVOT_H
#define SINGLE_PIVOT_H

#include "mcobs_handler.h"
#include "correlator_matrix_info.h"
#include "task_utils.h"
#include "task_handler.h"
#include "diag_corr_set.h"

#if defined COMPLEXNUMBERS
  typedef CMatrix                                TransMat;
  typedef HermDiagonalizerWithMetric             DiagonalizerWithMetric; 
  typedef Array<std::complex<double> >           ArrayBuf;
  typedef CVector                                VVector;
#elif defined REALNUMBERS
  typedef RMatrix                           TransMat;
  typedef RealSymDiagonalizerWithMetric     DiagonalizerWithMetric;
  typedef Array<double>                     ArrayBuf;
  typedef RVector                           VVector;
#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif


// ***********************************************************************************
// *                                                                                 *
// *   This implements the Single Pivot method.  For given Hermitian correlation     *
// *   matrix, three time slices are chosen:  tauN <= tau0 < tauD.  The              *
// *   matrix at rescaling time "tauN" is used to rescale the correlation matrix:    *
// *                                                                                 *
// *      C[i,j](t) = rawC[i,j](t) / sqrt(  rawC[i,i](tauN) rawC[j,j](tauN) )        *
// *                                                                                 *
// *   For metric time "tau0" and diagonalization time "tauD", the transformation    *
// *   matrix "V" whose columns are the eigenvectors that satisfy                    *
// *                                                                                 *
// *          C(tauD) V =   C(tau0) V  Lambda                                        *
// *                                                                                 *
// *   where "Lambda" is a diagonal matrix of eigenvalues, are determined.  The      *
// *   rescaled correlation matrix C(t) is then "rotated"                            *
// *                                                                                 *
// *         Crot(t) = V^dagger C(t) V                                               *
// *                                                                                 *
// *   and we have                                                                   *
// *                                                                                 *
// *         Crot(tauD) = V^dagger C(tauD) V = Lambda                                *
// *         Crot(tau0) = V^dagger C(tau0) V = I (identity)                          *
// *                                                                                 *
// *   If there are VEVs, these are also rotated and subtracted.  An                 *
// *   additional rephasing is done so the rotated VEVs are real and positive.       *
// *   In the single pivot method, this rotation is done bin-by-bin.                 *
// *   Only the diagonal elements of the "rotated" correlation matrix are            *
// *   determined since these are the only ones needed (currently).                  *
// *                                                                                 *
// *   Use the static member "initiateSinglePivot" to initiate a single pivot.       *
// *   Based on the input XML, it either calculates and creates a new object         *
// *   (new memory), reads a previously calculated object from file (new memory),    *
// *   or gets a pointer to a previously calculated object saved in the task         *
// *   handler data map (no new memory).  If the input XML requests to store the     *
// *   object in the task handler data map, then "keep_in_task_map" is returned      *
// *   true, and the end user need not worry about deallocating the memory since     *
// *   this will be automatically done by the destructor of the task handler.        *
// *   If "keep_in_task_map" is returned false, then the user is responsible for     *
// *   deleting the object.                                                          *
// *                                                                                 *
// *   Input XML to create a new pivot:                                              *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <RotatedCorrelator>                                                     *
// *           <GIOperator>...</GIOperator>                                          *
// *         </RotatedCorrelator>                                                    *
// *         <AssignName>PivTester</AssignName>  (optional)                          *
// *         <NormTime>3</NormTime>                                                  *
// *         <MetricTime>6</MetricTime>                                              *
// *         <DiagonalizeTime>12</DiagonalizeTime>                                   *
// *         <MinimumInverseConditionNumber>0.01</MinimumInverseConditionNumber>     *
// *         <NegativeEigenvalueAlarm>-0.01</NegativeEigenvalueAlarm>  (optional)    *
// *         <CheckMetricErrors/>    (optional)                                      *
// *         <CheckCommonMetricMatrixNullSpace/>    (optional)                       *
// *         <WritePivotToFile>    (optional)                                        *
// *            <PivotFileName>pivot_test</PivotFileName>                            *
// *            <Overwrite/>                                                         *
// *         </WritePivotToFile>                                                     *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in a file:               *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <ReadPivotFromFile>                                                     *
// *            <PivotFileName>pivot_file</PivotFileName>                            *
// *         </ReadPivotFromFile>                                                    *
// *         ... other tasks ...                                                     *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in memory:               *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <GetFromMemory>                                                         *
// *            <IDName>PivTester</IDName>                                           *
// *         </GetFromMemory>                                                        *
// *         ... other tasks ...                                                     *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *                                                                                 *
// *   Input XML for tasks:                                                          *
// *                                                                                 *
// *      <WriteRotatedCorrToFile>    (optional)                                     *
// *         <RotatedCorrFileName>rotated_corr_bins</RotatedCorrFileName>            *
// *         <Overwrite/>                                                            *
// *      </WriteRotatedCorrToFile>                                                  *
// *                                                                                 *
// ***********************************************************************************



class SinglePivotOfCorrMat : public TaskHandlerData
{

   MCObsHandler *m_moh;
   const CorrelatorMatrixInfo *m_cormat_info;
   GenIrrepOperatorInfo *m_rotated_info;
   const TransMatrix *m_Zmat, *m_transmat;
   uint m_tauN, m_tau0, m_tauD;
   double m_min_inv_condnum;
   double m_neg_eig_alarm;
   DiagonalCorrelatorSet *m_rotcorset;

#ifndef NO_CXX11
    SinglePivotOfCorrMat() = delete;
    SinglePivotOfCorrMat(const SinglePivotOfCorrMat& copy) = delete;
    SinglePivotOfCorrMat& operator=(const SinglePivotOfCorrMat& copy) = delete;
#else
    SinglePivotOfCorrMat();
    SinglePivotOfCorrMat(const SinglePivotOfCorrMat& copy);
    SinglePivotOfCorrMat& operator=(const SinglePivotOfCorrMat& copy);
#endif

 public:

   SinglePivotOfCorrMat(TaskHandler& taskhandler, ArgsHandler& xmlin,
                        LogHelper& xmlout);
   ~SinglePivotOfCorrMat();

   static SinglePivotOfCorrMat* initiateSinglePivot(
                   TaskHandler& taskhandler, ArgsHandler& xmlin,
                   LogHelper& xmlout, bool& keep_in_task_map);


   uint getNumberOfOperators() const;

   uint getNumberOfLevels() const;

   GenIrrepOperatorInfo getRotatedOperator() const;

   bool isVEVsubtracted() const;


   void doRotation(uint tmin, uint tmax, LogHelper& xmllog);
 
   void writeRotated(uint tmin, uint tmax, const std::string& corrfile,
                     bool overwrite, LogHelper& xmlout);

   void computeZMagnitudesSquared(Matrix<MCEstimate>& ZMagSq);

 private:

   static SinglePivotOfCorrMat* initiateFromMemory(TaskHandler& taskhandler, 
                  ArgsHandler& xml_in, LogHelper& xmlout);

   static bool putInMemory(TaskHandler& taskhandler, ArgsHandler& xmlin,
                           LogHelper& xmlout, SinglePivotOfCorrMat* pivot);

   void initiate_new(ArgsHandler& xml_in, LogHelper& xmlout);
   void initiate_from_file(ArgsHandler& xml_in, LogHelper& xmlout);
   void clear();

   void create_pivot(LogHelper& xmllog, bool checkMetricErrors, 
                     bool checkCommonNullSpace);
   void do_vev_rotation();
   void do_corr_rotation(uint timeval, bool diagonly);
   void write_to_file(const std::string& fname, bool overwrite);

};


#endif
