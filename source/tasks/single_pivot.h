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
// *   This implements the Single Pivot method.  To apply this method, one first     *
// *   creates a single pivot using the "<SinglePivotInitiate>" tag.  The result     *
// *   of this initiation can be saved to file with a "<WritePivotToFile>" tag for   *
// *   subsequent use in other runs, or the pivot can be directly used.  Later       *
// *   runs can initiate a pivot by reading from file with a "<ReadPivotFromFile>"   *
// *   tag, or if it is already stored in persistent memory from a previous task, a  *
// *   "<GetFromMemory>" tag can be used.  Regardless of how the pivot is            *
// *   initialized, the second step is to use the pivot to rotate a correlation      *
// *   matrix using the "<WriteRotatedCorrToFile>" tag.  The 3rd step is to perform  *
// *   fits to the diagonal elements of the rotated correlation matrix.  These fit   *
// *   results can be saved to file, especially the amplitude factors.  The last     *
// *   step is to use the pivot information and the fit amplitudes to determine the  *
// *   operator overlap "Z" factors using a  "<DoCorrMatrixZMagSquares>" tag.        *
// *   The initial ordering of the levels is based on the diagonalization at time    *
// *   separation "tauD".  However, this ordering may not agree with that from       *
// *   the final fit energies. By inserting fit energy information for all levels,   *
// *   the level ordering can be changed to agree with increasing fit energy.        *
// *                                                                                 *
// *   For given Hermitian correlation matrix, three time slices are chosen:         *
// *   tauN <= tau0 < tauD.  The matrix at rescaling time "tauN" is used to rescale  *
// *   the correlation matrix:                                                       *
// *                                                                                 *
// *      C[i,j](t) = rawC[i,j](t) / sqrt(  rawC[i,i](tauN) rawC[j,j](tauN) )        *
// *                                                                                 *
// *   For metric time "tau0" and diagonalization time "tauD", the following         *
// *   procedure is followed:                                                        *
// *                                                                                 *
// *   (1) The eigenvalues and eigenvectors of the NxN matrix C(tau0) are            *
// *       determined using the full ensemble of configurations.  Let "L0max"        *
// *       denote the eigenvalue of largest magnitude.  Put the N0 <= N eigenvectors *
// *       associated with eigenvalues greater than                                  *
// *       "L0max" * "MinimumInverseConditionNumber" into the columns of a matrix    *
// *       named "P0".  Define the N0 x N0 matrices                                  *
// *          Ctilde(tau0) = P0^dag C(tau0) P0                                       *
// *          Ctilde(t) = P0^dag C(t) P0                                             *
// *          Gtilde(t) = Ctilde(tau0)^(-1/2) Ctilde(t) Ctilde(tau0)^(-1/2)          *
// *                                                                                 *
// *   (2) Solve for the eigenvalues and eigenvectors of Gtilde(tauD) using the      *
// *       full ensemble only.  Let "Ltmax" denote the eigenvalue of largest         *
// *       magnitude.  Put the NP <= N0 eigenvectors associated with eigenvalues     *
// *       greater than "Ltmax" * "MinimumInverseConditionNumber" into the columns   *
// *       of a matrix called VtildeD.                                               *
// *                                                                                 *
// *   (3) Evaluate the diagonal elements of the "rotated" correlation matrix        *
// *          Dtilde(t) = VtildeD^dag Gtilde(t) VtildeD,                             *
// *       on the individual bins.  On the full ensemble, we will have               *
// *           Dtilde(tau0) = 1,  Dtilde(tauD) = diagonal.                           *
// *       For different resamplings, the above relations will not be true.          *
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
// *         <CorrelatorMatrixInfo> ... </CorrelatorMatrixInfo>                      *
// *         <ImprovedOperators> ... </ImprovedOperators>  (optional)                *
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
// *         <PrintTransformationMatrix/> (optional)                                 *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *   The <RotatedCorrelator> tag specifies the name to give the rotated            *
// *   operators.  Any integer index specified is ignored.  If the matrix to         *
// *   be rotated is "N" x "N", then the N rotated operators will all have the       *
// *   same isospin and irrep labels and ID name, but the ID index will vary         *
// *   from 0 to N-1.                                                                *
// *                                                                                 *
// *   The <CorrelatorMatrixInfo> tag specifies the correlator matrix of operators   *
// *   to be rotated. The tag <HermitianMatrix> must be present.  If the tag         *
// *   <ImprovedOperators> is present, this means that the operators are linear      *
// *   combinations of another set of operators.  The linear combinations must then  *
// *   be given within this tag in the form:                                         *
// *                                                                                 *
// *    <ImprovedOperators>                                                          *
// *      <ImprovedOperator>                                                         *
// *         <OpName>                                                                *
// *          <GIOperatorString>isotriplet P=(0,0,0) A1gp_1 RotTester 0</GIOperatorString>
// *         </OpName>                                                               *
// *         <OpTerm>                                                                *
// *           <BLOperatorString>pion P=(0,0,0) A1gp_1 SD_0</BLOperatorString>       *
// *           <Coefficient>(-0.0593735752248,0.0421528577847)</Coefficient>         *
// *         </OpTerm>                                                               *
// *          ...                                                                    *
// *      </ImprovedOperator>                                                        *
// *       ....                                                                      *
// *    </ImprovedOperators>                                                         *
// *                                                                                 *
// *   If <AssignName> is present, the pivot is inserted into the task handler       *
// *   data map, and can be accessed using this ID tag.  If not present, the         *
// *   pivot is not stored in persistent memory.  If the pivot is written to         *
// *   file, this ID name is NOT put into the file, allowing subsequent programs     *
// *   to assign whatever name they wish.                                            * 
// *                                                                                 *
// *   The tag "MinimumInverseConditionNumber" has already been explained above,     *
// *   but its purpose is to remove noisy states.  States that are not sufficiently  *
// *   independent of the other states can become dominated by noise.                *
// *   The fractional errors in the diagonal elements of rawC(tau0), rawC(tauD) are  *
// *   printed out for informational purposes.                                       *
// *                                                                                 *
// *   If the tag "CheckMetricErrors" is present, the eigenvalues of the largest     *
// *   and smallest magnitudes of C(tau0) are determined by jackknife and the        *
// *   condition number with errors is output.  This can help in determining         *
// *   the "MinimumInverseConditionNumber" to use.                                   *
// *                                                                                 *
// *   If the tag "NegativeEigenvalueAlarm" is set to a negative value, then         *
// *   if any eigenvalues are less than this value, this fact is reported.           *
// *                                                                                 *
// *   If the tag "CheckCommonMetricMatrixNullSpace" is present, then in             *
// *   Gtilde(tauD), the null space of Ctilde(tauD) is checked to see that           *
// *   it contains the entire null space of Ctilde(tau0).  This is a desirable       *
// *   property when removing noisy eigenvectors.  Finding this false indicates      *
// *   caution in interpreting the overlap factors.                                  *
// *                                                                                 *
// *   If the "WritePivotToFile" tag is assigned, the information in the pivot       *
// *   is written out to a file so that it can be input by later sigmond runs.       *
// *   The pivot file is an IOMap (a single integer key) with an XML header string   *
// *   and binary data. The XML header string contains                               *
// *      - the correlator matrix info                                               *
// *      - the rotated correlator ID info                                           *
// *      - the norm time, metric time, diagonalize time                             *
// *      - the minimum inverse condition number                                     *
// *   The binary data contains                                                      *
// *      - the rotation (transformation) matrix (key = 0)                           *
// *      - the matrix needed to compute the overlap Zmag squares (key = 1)          *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in a file:               *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <ReadPivotFromFile>                                                     *
// *            <PivotFileName>pivot_file</PivotFileName>                            *
// *         </ReadPivotFromFile>                                                    *
// *         <AssignName>PivTester</AssignName>  (optional)                          *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in memory:               *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <GetFromMemory>                                                         *
// *            <IDName>PivTester</IDName>                                           *
// *         </GetFromMemory>                                                        *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *                                                                                 *
// *   Input XML for writing rotated correlators to file (as bins or samplings):     *
// *                                                                                 *
// *      <WriteRotatedCorrToFile>    (optional)                                     *
// *         <RotatedCorrFileName>rotated_corr_bins</RotatedCorrFileName>            *
// *         <Type>bins</Type>  (or samplings)                                       *
// *         <Overwrite/>                                                            *
// *      </WriteRotatedCorrToFile>                                                  *
// *                                                                                 *
// *   The member "computeZMagnitudesSquared" computes the overlap factors,          *
// *   returning them in a matrix of MCEstimates  |Z(opindex,level)|^2 for all       *
// *   operators for all levels.                                                     *
// *                                                                                 *
// *                                                                                 *
// ***********************************************************************************


// ***********************************************************************************
// *                                                                                 *
// *  Implementation notes:                                                          *
// *                                                                                 *
// *  - The original correlator information is stored in the object pointed to       *
// *    by "m_cormat_info".                                                          *
// *  - Only the diagonal elements of the rotated correlator are dealt with.         *
// *    These are specified by the information stored in the object pointed to       *
// *    by "m_rotated_info".  All of the diagonal elements have the same             *
// *    name, by differ in their ID index.  The diagonal elements are ordered        *
// *    according to their effective energy at the diagonalization time "tauD".      *
// *    This ordering may not be exactly the same as the eventually order of         *
// *    energies as determined from large time fits.                                 *
// *  - "m_tauN" is the time separation used for rescaling by norms.                 *
// *  - "m_tau0" is the time separation used for the "metric" C(tau0)                *
// *  - "m_tauD" is the time separation used for the diagonalization                 *
// *  - "m_transmat" points to the transformation matrix "R" that gives the          *
// *    rotated correlator matrix in terms of the original raw matrix using          *
// *       R^dagger C_orig R  where                                                  *
// *           R[i,j]= 1/sqrt(C_orig[i,i]) P0[i,k] Ctilde(tau0)^(-1/2)[k,j]          *
// *           Ctilde(t) = P0^dag C(t) P0                                            *
// *           C[i,j](t) = C_orig[i,j](t)/sqrt(C_orig[i,i](tauN)*C_orig[j,j](tauN))  *
// *           columns of P0 are retained eigenvectors of C(tau0)                    *
// *  - "m_Zmat" points to the matrix needed to compute the |Z|^2 overlaps.          *
// *    Once fits to each diagonal rotated correlator are done using                 *
// *            |Ztilde[level]|^2 exp(-E[level]*t)                                   *
// *    the overlaps of the original operators are obtained from                     *
// *       |Z[opindex,level]|^2 = |(*m_Zmat)(opindex,level)|^2 * |Ztilde[level]|^2   *
// *  - "m_ampkeys" is a map that stores the MCObsInfo keys for getting the          *
// *    fit values for |Ztilde[level]|^2.                                            *
// *  - "m_min_inv_condnum" is the minimum inverse condition number to use           *
// *    is removing noise.                                                           *
// *  - "m_neg_eig_alarm" is the threshold (negative value) for reporting a          *
// *    negative value for an eigenvalue of a correlation matrix.                    *
// *                                                                                 *
// *  - The constructor initializes all quantities, except m_ampkeys. The            *
// *    constructor does the diagonalizations needed and evaluates and stores        *
// *    the transformation and Z calculation matrices.  It also writes these         *
// *    quantities into a pivot file, if requested.  The constructor does NOT        *
// *    carry out the correlator rotations.                                          *
// *  - To keep the pivot in persistent memory, do not use the constructor           *
// *    directly, but instead call the static routine "initiateSinglePivot"          *
// *    which inserts the pivot in the TaskHandlerData map.                          *
// *  - The computation of the diagonal elements of the rotated correlation          *
// *    matrix is accomplished by the "doRotation" member, which puts the            *
// *    diagonal elements of the rotated correlation matrix in memory.               *
// *  - The member "writeRotated" can be used to store the diagonal elements         *
// *    of the rotated correlation matrix in a file.                                 *
// *  - The member "insertEnergyFitInfo" can be used for each level to indicate      *
// *    the MCObsInfo key for the fit values of the level energies.                  *
// *  - Once all of the fit energies are inserted, the member                        *
// *    "reorderLevelsByFitEnergy" can be used to reorder the level indices to       *
// *    match increasing fit energy values.                                          *
// *  - The member "insertAmplitudeFitInfo" can be used for each level to            *
// *    indicate the MCObsInfo key for the fit values of |Ztilde[level]|^2.          *
// *  - Once all of these infos are inserted, the overlap factors of the             *
// *    original operators can be computed using the member routine                  *
// *    "computeZMagnitudesSquared" which returns the overlaps in                    *
// *    a  Matrix<MCEstimate>& ZMagSq.                                               *
// *                                                                                 *
// *                                                                                 *
// ***********************************************************************************



class SinglePivotOfCorrMat : public TaskHandlerData
{

   MCObsHandler *m_moh;
   const CorrelatorMatrixInfo *m_cormat_info, *m_orig_cormat_info;
   GenIrrepOperatorInfo *m_rotated_info;
   const TransMatrix *m_Zmat, *m_transmat, *m_imp_trans;
   uint m_tauN, m_tau0, m_tauD;
   double m_min_inv_condnum;
   double m_neg_eig_alarm;
   std::map<uint,MCObsInfo> m_ampkeys;
   std::map<uint,MCObsInfo> m_energykeys;
   std::vector<uint> m_reorder;

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

   const std::set<OperatorInfo>& getOperators() const;

   GenIrrepOperatorInfo getRotatedOperator() const;

   bool subtractVEV() const;


   void doRotation(uint tmin, uint tmax, LogHelper& xmllog);
 
   void writeRotated(uint tmin, uint tmax, const std::string& corrfile,
                     bool overwrite, LogHelper& xmlout, bool bins);


   void insertAmplitudeFitInfo(uint level, const MCObsInfo& ampinfo);

   MCObsInfo getAmplitudeKey(uint level) const;

   bool allAmplitudeFitInfoAvailable() const
    {return (m_ampkeys.size()>0)&&(m_ampkeys.size()==getNumberOfLevels());}

   void insertEnergyFitInfo(uint level, const MCObsInfo& ampinfo);

   MCObsInfo getEnergyKey(uint level) const;

   bool allEnergyFitInfoAvailable() const
    {return (m_energykeys.size()>0)&&(m_energykeys.size()==getNumberOfLevels());}


   void reorderLevelsByFitEnergy();

   void clearReordering();

   bool areLevelsReordered() const
    {return !(m_reorder.empty());}

         //  get |Z(opindex,level)|^2 for all operators for all levels

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
   void print_trans(LogHelper& xmllog);
   void read_trans(ArgsHandler& xmlin, 
          std::map<OperatorInfo,std::map<OperatorInfo,Scalar> >& trans);
   void setup_improved_operators(const std::map<OperatorInfo,
                std::map<OperatorInfo,Scalar> >& trans);

};


#endif
