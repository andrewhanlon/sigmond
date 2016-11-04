#ifndef CORRELATOR_MATRIX_INFO_H
#define CORRELATOR_MATRIX_INFO_H

#include <map>
#include <vector>
#include "operator_info.h"
#include "correlator_info.h"

class CorrMatrixIterator;

// *******************************************************************
// *                                                                 *
// *   The class "CorrelatorMatrixInfo" is defined in this file.     *
// *   This class stores information about the operators that        *
// *   define a temporal correlation matrix, including whether       *
// *   VEV subtraction is needed and whether or not the matrix is    *
// *   Hermitian. "CorrelatorMatrixInfo" construction by XML content *
// *   requires XML in the following format:                         *
// *                                                                 *
// *       <CorrelatorMatrixInfo>                                    *
// *           <Operator>...</Operator>                              *
// *           <Operator>...</Operator>                              *
// *        or   <OperatorString>...</OperatorString>                *
// *                ...                                              *
// *           <HermitianMatrix/>    (optional)                      *
// *           <SubtractVEV/>    (optional)                          *
// *           <AssignName>descriptive_name</AssignName>             *
// *       </CorrelatorMatrixInfo>                                   *
// *                                                                 *
// *   or by                                                         *
// *                                                                 *
// *       <CorrelatorMatrixInfo>                                    *
// *           <Name>descriptive_name</Name>                         *
// *       </CorrelatorMatrixInfo>                                   *
// *                                                                 *
// *   See the file "operator_info.h" for more information about     *
// *   the XML format needed inside the <Operator> or                *
// *   <OperatorString> tags.                                        *
// *                                                                 *
// *   The <AssignName> tag can be used to associate a string        *
// *   (with no blanks) with a CorrelatorMatrixInfo object.          *
// *   Then this correlation matrix can be referred to in            *
// *   subsequent tasks using just a <Name> tag.  A static map is    *
// *   used to keep track of all named correlation matrix objects.   *
// *                                                                 *
// *   This class provides a means of iterating over all correlators *
// *   in the matrix:                                                *
// *                                                                 *
// *     CorrelatorMatrixInfo cormat(....);                          *
// *     for (cormat.begin();!cormat.end();cormat++){                *
// *          cormat.getCurrentCorrelatorInfo();    }                *
// *                                                                 *
// *                                                                 * 
// *******************************************************************


class CorrelatorMatrixInfo
{

   std::set<OperatorInfo> m_opinfos;
   bool m_hermitian;
   bool m_vevsubt;

   CorrelatorInfo *m_current;
   std::set<OperatorInfo>::const_iterator src_it, snk_it;

 public:

   CorrelatorMatrixInfo(XMLHandler& xml_in);

   CorrelatorMatrixInfo(const std::set<OperatorInfo>& inops, bool herm, bool subvev)
       : m_opinfos(inops), m_hermitian(herm), m_vevsubt(subvev), m_current(0) {}

   CorrelatorMatrixInfo(const CorrelatorMatrixInfo& cor) : m_opinfos(cor.m_opinfos),
         m_hermitian(cor.m_hermitian), m_vevsubt(cor.m_vevsubt), m_current(0) {}

   CorrelatorMatrixInfo& operator=(const CorrelatorMatrixInfo& cor)
    {m_opinfos=cor.m_opinfos;
     m_hermitian=cor.m_hermitian;
     m_vevsubt=cor.m_vevsubt;
     m_current=0;
     return *this;}

   ~CorrelatorMatrixInfo()
    {delete m_current;}

   void setName(const std::string& corrmatname);

   CorrelatorMatrixInfo& begin()
    {src_it=m_opinfos.begin(); snk_it=m_opinfos.begin();
     if (m_current==0) m_current=new CorrelatorInfo(*snk_it,*src_it);
     else *m_current=CorrelatorInfo(*snk_it,*src_it);
     return *this;}

   bool end() const
    {return ((src_it==m_opinfos.end())&&(snk_it==m_opinfos.end()));}

   CorrelatorMatrixInfo& operator++()
    {
     if (m_current==0) 
        throw(std::invalid_argument("Cannot increment uninitialized iterator in CorrelatorMatrixInfo"));
     src_it++;
     if (src_it==m_opinfos.end()){
        snk_it++;
        if (m_hermitian) src_it=snk_it;
        else src_it=m_opinfos.begin();}
     if (snk_it==m_opinfos.end()){
        src_it=m_opinfos.end();
        return *this;}
     *m_current=CorrelatorInfo(*snk_it,*src_it);
     return *this;}

   CorrelatorMatrixInfo operator++(int k)
    {return ++(*this);}

   const CorrelatorInfo& getCurrentCorrelatorInfo() const
    {if (m_current==0) throw(std::invalid_argument("No current CorrelatorInfo"));
     return *m_current;}


   bool isHermitian() const
     { return m_hermitian; }

   void setHermitian()
     { m_hermitian = true; }

   void setNoHermitian()
     { m_hermitian = false; }

   bool isVEVSubtracted() const
     { return m_vevsubt; }

   void setNoVEVSubtracted()
     { m_vevsubt=false;}

   void setVEVSubtracted()
     { m_vevsubt=true;}

   const std::set<OperatorInfo>& getOperators() const
     {return m_opinfos;}

   uint getNumberOfOperators() const
     {return m_opinfos.size();}

   std::string output(bool longform=false, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=false) const;  // XML output


   bool operator==(const CorrelatorMatrixInfo& rhs) const;

   bool operator!=(const CorrelatorMatrixInfo& rhs) const;

   bool operator<(const CorrelatorMatrixInfo& rhs) const;

 private:

   static std::map<std::string,CorrelatorMatrixInfo> m_cormat_namemap;

};



// **************************************************************
#endif  
