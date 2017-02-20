#ifndef DIAG_CORR_SET_H
#define DIAG_CORR_SET_H

#include <vector>
#include <map>
#include "operator_info.h"
#include "model_tcorr.h"

// *****************************************************************************
// *                                                                           *
// *   This is a simple class useful for storing information about a set of    *
// *   diagonal temporal correlators.  The CorrelatorInfo's are available,     *
// *   as well as information related to a fit to each correlator.             *
// *   The first XML construction below can be used regardless of the order    *
// *   of the operators in the XML.  The second construction with the          *
// *   <Sequential> tag takes the operators in the order given in the XML.     *
// *                                                                           *
// *   Input XML:                                                              *
// *                                                                           *
// *    <DiagonalCorrelatorSet>                                                *
// *       <DiagonalCorrelator>                                                *
// *          <OperatorString>pion P=(0,0,0) A1um_1 SS_0</OperatorString>      *
// *          <OperatorIndex>0</OperatorIndex>                                 *
// *       </DiagonalCorrelator>                                               *
// *       <DiagonalCorrelator>                                                *
// *          <OperatorString>kaon P=(0,0,0) A1u_1 SS_0</OperatorString>       *
// *          <OperatorIndex>1</OperatorIndex>                                 *
// *       </DiagonalCorrelator>                                               *
// *       <SubtractVEV/>   (optional)                                         *
// *    </DiagonalCorrelatorSet>                                               *
// *                                                                           *
// *    or                                                                     *
// *                                                                           *
// *    <DiagonalCorrelatorSet>                                                *
// *      <Sequential>                                                         *
// *          <OperatorString>pion P=(0,0,0) A1um_1 SS_0</OperatorString>      *
// *          <OperatorString>kaon P=(0,0,0) A1u_1 SS_0</OperatorString>       *
// *      </Sequential>                                                        *
// *       <SubtractVEV/>   (optional)                                         *
// *    </DiagonalCorrelatorSet>                                               *
// *                                                                           *
// *                                                                           *
// *    or                                                                     *
// *                                                                           *
// *    <DiagonalCorrelatorSet>                                                *
// *      <RotatedSequential>                                                  *
// *         <GIOperatorString>....</GIOperatorString> (can omit ID index)     *
// *         <NumberOfLevels>34</NumberOfLevels>                               *
// *      </RotatedSequential>                                                 *
// *       <SubtractVEV/>   (optional)                                         *
// *    </DiagonalCorrelatorSet>                                               *
// *                                                                           *
// *                                                                           *
// *    or                                                                     *
// *                                                                           *
// *    <DiagonalCorrelatorSet>                                                *
// *      <CorrelatorMatrixInfo>                                               *
// *      ....                                                                 *
// *      </CorrelatorMatrixInfo>                                              *
// *    </DiagonalCorrelatorSet>                                               *
// *                                                                           *
// *                                                                           *
// *****************************************************************************



class DiagonalCorrelatorSet
{

   std::vector<OperatorInfo> m_opset;
   std::map<OperatorInfo,TCorrFitInfo> m_fitinfos;
   bool m_subvev;

 public:

   DiagonalCorrelatorSet() : m_subvev(false) {}

   DiagonalCorrelatorSet(XMLHandler& xmlin);

   DiagonalCorrelatorSet(const DiagonalCorrelatorSet& cor) 
        : m_opset(cor.m_opset),  m_fitinfos(cor.m_fitinfos),
          m_subvev(cor.m_subvev) {}

   DiagonalCorrelatorSet& operator=(const DiagonalCorrelatorSet& cor)
    {m_opset=cor.m_opset;
     m_fitinfos=cor.m_fitinfos;
     m_subvev=cor.m_subvev;
     return *this;}

   ~DiagonalCorrelatorSet(){}

   void outputOperatorInfos(XMLHandler& xmlout);

   void output(XMLHandler& xmlout);

   void addCorrelator(const OperatorInfo& anop);

   void insertFitResult(const OperatorInfo& anop, 
                        const TCorrFitInfo& fitresult);

   void insertFitResult(uint opnum, 
                        const TCorrFitInfo& fitresult);

   void setVEVsubtractedOn() {m_subvev=true;}

   void setVEVsubtractedOff() {m_subvev=false;}

   void setVEVsubtracted(bool svev) {m_subvev=svev;}

   uint getNumberOfCorrelators() const
    {return m_opset.size();}

   CorrelatorInfo getCorrelatorInfo(uint opnum) const;

   bool isVEVsubtracted() const
    {return m_subvev;}

   MCObsInfo getEnergyKey(uint opnum) const;

   MCObsInfo getAmplitudeKey(uint opnum) const;

   const TCorrFitInfo& getFitInfo(uint opnum) const;

   bool allFitInfoAvailable() const
    {return (m_fitinfos.size()==m_opset.size());}

};


#endif
