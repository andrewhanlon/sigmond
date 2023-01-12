#ifndef ROLLING_PIVOT_H
#define ROLLING_PIVOT_H

#include "mcobs_handler.h"
#include "correlator_matrix_info.h"
#include "task_utils.h"
#include "task_handler.h"
#include "diag_corr_set.h"

// using namespace std;

#if defined COMPLEXNUMBERS
  typedef CMatrix                                TransMat;
  typedef HermDiagonalizerWithMetric             DiagonalizerWithMetric; 
  typedef Array<double>                          RArrayBuf;
  typedef CVector                                VVector;
  typedef VectorPinner<std::complex<double> >    LevelPinner;
#elif defined REALNUMBERS
  typedef RMatrix                           TransMat;
  typedef RealSymDiagonalizerWithMetric     DiagonalizerWithMetric;
  typedef Array<double>                     RArrayBuf;
  typedef RVector                           VVector;
  typedef VectorPinner<double>              LevelPinner;
#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif


// ***********************************************************************************
// *                                                                                 *
// *   This implements the Rolling Pivot method.  To apply this method, one first    *
// *   creates a rolling pivot using the "<RollingPivotInitiate>" tag.  The result   *
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
// *   separation "tauZ".  However, this ordering may not agree with that from       *
// *   the final fit energies. By inserting fit energy information for all levels,   *
// *   the level ordering can be changed to agree with increasing fit energy.        *
// *                                                                                 *
// *   For given Hermitian correlation matrix, three time slices are chosen:         *
// *   tauN <= tau0 < tauZ.  The matrix at rescaling time "tauN" is used to rescale  *
// *   the correlation matrix:                                                       *
// *                                                                                 *
// *      C[i,j](t) = rawC[i,j](t) / sqrt(  rawC[i,i](tauN) rawC[j,j](tauN) )        *
// *                                                                                 *
// *   For metric time "tau0" and Zmag time "tauZ", the following                    *
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
// *   (2) Solve for the eigenvalues and eigenvectors of Gtilde(tauZ) using the      *
// *       full ensemble only.  Let "Ltmax" denote the eigenvalue of largest         *
// *       magnitude.  Put the NP <= N0 eigenvectors associated with eigenvalues     *
// *       greater than "Ltmax" * "MinimumInverseConditionNumber" into the columns   *
// *       of a matrix called VtildeD.                                               *
// *                                                                                 *
// *   (3) Evaluate the diagonal elements of the "rotated" correlation matrix        *
// *          Dtilde(tauZ) = VtildeD^dag Gtilde(tauZ) VtildeD,                       *
// *       on the individual bins.                                                   *
// *                                                                                 *
// *   (4) Rotation matrix for tauZ will be stored for ZMag analysis.                *
// *
// *   (5) Steps (2) and (3) will be repeated for every time slice beginning with    *
// *       t = tauZ - 1 (if exists) and decremented until the code reaches the       *
// *       first time slice. The VectorPinner will be used to match and reorder the  *
// *       eigenvectors between time slices in accordance with the pivot at tauZ.    *
// *       When it reaches the first time slice, the vector pinner is reset to the   *
// *       tauZ pivot, then continues diagonalizing the correlation matrix for       *
// *       times t = tauZ+1 to the final time slice in the same manner. On the full  *
// *       ensemble, we will have  Dtilde(tau0) = 1,  Dtilde(t) = diagonal.          *
// *       For different resamplings, the above relations will not be true.          *
// *                                                                                 *
// *   If there are VEVs, these are not implemented correctly and need additional    *
// *   development and testing.                                                      *
// *                                                                                 *
// *   Use the static member "initiateRollingPivot" to initiate a pivot.             *
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
// *      <RollingPivotInitiate>                                                     *
// *         <RotatedCorrelator>                                                     *
// *           <GIOperator>...</GIOperator>                                          *
// *         </RotatedCorrelator>                                                    *
// *         <AssignName>PivTester</AssignName>  (optional)                          *
// *         <CorrelatorMatrixInfo> ... </CorrelatorMatrixInfo>                      *
// *         <NormTime>3</NormTime>                                                  *
// *         <MetricTime>6</MetricTime>                                              *
// *         <ZMatrixTime>12</ZMatrixTime>                                           *
// *         <MinimumInverseConditionNumber>0.01</MinimumInverseConditionNumber>     *
// *         <NegativeEigenvalueAlarm>-0.01</NegativeEigenvalueAlarm>  (optional)    *
// *         <WarningFraction>0.7</WarningFraction>  (optional)                      *
// *         <CheckMetricErrors/>    (optional)                                      *
// *         <CheckCommonMetricMatrixNullSpace/>    (optional)                       *
// *         <WritePivotToFile>    (optional)                                        *
// *            <PivotFileName>pivot_test</PivotFileName>                            *
// *            <Overwrite/>                                                         *
// *         </WritePivotToFile>                                                     *
// *      </RollingPivotInitiate>                                                    *
// *                                                                                 *
// *   The <RotatedCorrelator> tag specifies the name to give the rotated            *
// *   operators.  Any integer index specified is ignored.  If the matrix to         *
// *   be rotated is "N" x "N", then the N rotated operators will all have the       *
// *   same isospin and irrep labels and ID name, but the ID index will vary         *
// *   from 0 to N-1.                                                                *
// *                                                                                 *
// *   The <CorrelatorMatrixInfo> tag specifies the original correlator matrix       *
// *   of operators to be rotated. The tag <HermitianMatrix> must be present.        *
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
// *                                                                                 *
// *   If the tag "CheckMetricErrors" is present, the eigenvalues of the largest     *
// *   and smallest magnitudes of C(tau0) are determined by jackknife and the        *
// *   condition number with errors is output.  This can help in determining         *
// *   the "MinimumInverseConditionNumber" to use.                                   *
// *                                                                                 *
// *   If the tag "NegativeEigenvalueAlarm" is set to a negative value, then         *
// *   if any eigenvalues are less than this value, this fact is reported.           *
// *                                                                                 *
// *   If the tag "WarningFraction" is set then the vector pinner will use that      *
// *   value as the minimum accepted overlap for vector pinner rather than the       *
// *   default 0.7.                                                                  *
// *                                                                                 *
// *   If the tag "CheckCommonMetricMatrixNullSpace" is present, then in             *
// *   Gtilde(tauZ), the null space of Ctilde(tauZ) is checked to see that           *
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
// *      <RollingPivotInitiate>                                                     *
// *         <ReadPivotFromFile>                                                     *
// *            <PivotFileName>pivot_file</PivotFileName>                            *
// *         </ReadPivotFromFile>                                                    *
// *         <AssignName>PivTester</AssignName>  (optional)                          *
// *      </RollingPivotInitiate>                                                    *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in memory:               *
// *                                                                                 *
// *      <RollingPivotInitiate>                                                     *
// *         <GetFromMemory>                                                         *
// *            <IDName>PivTester</IDName>                                           *
// *         </GetFromMemory>                                                        *
// *      </RollingPivotInitiate>                                                    *
// *                                                                                 *
// *                                                                                 *
// *   Input XML for tasks:                                                          *
// *                                                                                 *
// *      <WriteRotatedCorrToFile>    (optional)                                     *
// *         <RotatedCorrFileName>rotated_corr_bins</RotatedCorrFileName>            *
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
// *    according to their effective energy at the diagonalization time "tauZ".      *
// *    This ordering may not be exactly the same as the eventually order of         *
// *    energies as determined from large time fits.                                 *
// *  - "m_tauN" is the time separation used for rescaling by norms.                 *
// *  - "m_tau0" is the time separation used for the "metric" C(tau0)                *
// *  - "m_tauZ" is the time separation used for the zmag calculation and            *
// *       eigenvector ordering                                                      *
// *  - "m_refstart" points to the transformation matrix "R" that gives the          *
// *    rotated correlator matrix in terms of the original matrix using              *
// *       R^dagger C_orig R                                                         *
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
// *  - "m_vecpin" tracks the eigenvalue overlaps and reorders them according to     *
// *       the tauZ pivot                                                            *
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



class RollingPivotOfCorrMat : public TaskHandlerData
{

   MCObsHandler *m_moh;
   const CorrelatorMatrixInfo *m_cormat_info; //, *m_orig_cormat_info;
   GenIrrepOperatorInfo *m_rotated_info;
   DiagonalizerWithMetric *m_diag;
   const TransMatrix *m_refstart, *m_Zmat; //, *m_transmat, *m_imp_trans;
   uint m_taurecent;
   TransMatrix m_refrecent,m_phase_matrix,m_vev_rotator;
   LevelPinner m_vecpin; 
   uint m_tauN, m_tau0, m_tauZ;
   double m_min_inv_condnum;
   double m_neg_eig_alarm;
   std::map<uint,MCObsInfo> m_ampkeys;
   std::map<uint,MCObsInfo> m_energykeys;
   std::vector<uint> m_reorder;
   bool m_vevs_avail;

#ifndef NO_CXX11
    RollingPivotOfCorrMat() = delete;
    RollingPivotOfCorrMat(const RollingPivotOfCorrMat& copy) = delete;
    RollingPivotOfCorrMat& operator=(const RollingPivotOfCorrMat& copy) = delete;
#else
    RollingPivotOfCorrMat();
    RollingPivotOfCorrMat(const RollingPivotOfCorrMat& copy);
    RollingPivotOfCorrMat& operator=(const RollingPivotOfCorrMat& copy);
#endif

 public:

   RollingPivotOfCorrMat(TaskHandler& taskhandler, ArgsHandler& xmlin,
                        LogHelper& xmlout);
   ~RollingPivotOfCorrMat();

   static RollingPivotOfCorrMat* initiateRollingPivot(
                   TaskHandler& taskhandler, ArgsHandler& xmlin,
                   LogHelper& xmlout, bool& keep_in_task_map);
   uint getTauN() const
    {return m_tauN;}

   uint getTau0() const
    {return m_tau0;}

   uint getTauZ() const
    {return m_tauZ;}

   uint getNumberOfOperators() const;

   uint getNumberOfLevels() const;

   const std::set<OperatorInfo>& getOperators() const;

   GenIrrepOperatorInfo getRotatedOperator() const;

   bool subtractVEV() const;


   void doRotation(uint tmin, uint tmax, LogHelper& xmllog);
 
   void writeRotated(uint tmin, uint tmax, bool remove_off_diag, const std::string& corrfile, 
                       WriteMode wmode, LogHelper& xmlout, char mode, char file_format);


   void insertAmplitudeFitInfo(uint level, const MCObsInfo& ampinfo);

   MCObsInfo getAmplitudeKey(uint level) const;

   bool allAmplitudeFitInfoAvailable() const
    {return (m_ampkeys.size()>0)&&(m_ampkeys.size()==getNumberOfLevels());}


   void insertEnergyFitInfo(uint level, const MCObsInfo& ampinfo);

   MCObsInfo getEnergyKey(uint level) const;

   bool allEnergyFitInfoAvailable() const
    {return (m_energykeys.size()>0)&&(m_energykeys.size()==getNumberOfLevels());}

   void reorderLevelsByFitEnergy(LogHelper& xmllog);

         //  get |Z(opindex,level)|^2 for all operators for all levels

   void computeZMagnitudesSquared(Matrix<MCEstimate>& ZMagSq);
   
   std::string type(){return "RollingPivotOfCorrMat";}

 private:

   static RollingPivotOfCorrMat* initiateFromMemory(TaskHandler& taskhandler, 
                  ArgsHandler& xml_in, LogHelper& xmlout);

   static bool putInMemory(TaskHandler& taskhandler, ArgsHandler& xmlin,
                           LogHelper& xmlout, RollingPivotOfCorrMat* pivot);

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
