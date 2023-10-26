#ifndef LAPH_OPERATOR_INFO_H
#define LAPH_OPERATOR_INFO_H

#include <vector>
#include "momenta.h"
#include "args_handler.h"
#include "basic_laph_operator_info.h"
#include "gen_irrep_operator_info.h"


// *******************************************************************
// *                                                                 *
// *   Objects of class "OperatorInfo" store identifying information *
// *   about one particular QCD stationary-state operator.  An       *
// *   "OperatorInfo" object is either a "BasicLapHOperatorInfo"     *
// *   or a "GenIrrepOperatorInfo".  The XML used to create a        *
// *   "BasicLapHOperatorInfo" or a "GenIrrepOperatorInfo" is used.  *
// *   Implementation note: the right most 3 bits of icode[0]        *
// *   specify whether the object is a "BasicLapHOperatorInfo"       *
// *   or a "GenIrrepOperatorInfo": if these 3 bits are all 1, then  *
// *   the object is of type "GenIrrepOperatorInfo", otherwise it    *
// *   is of type "BasicLapHOperatorInfo".                           *
// *                                                                 *
// *******************************************************************


class OperatorInfo
{

   std::vector<unsigned int> icode;

 public:

   enum OpKind { BasicLapH, GenIrrep };

   OperatorInfo();

   OperatorInfo(XMLHandler& xml_in);

   OperatorInfo(const std::string& opstring, OperatorInfo::OpKind opkind = BasicLapH);

   OperatorInfo(const OperatorInfo& B) 
         : icode(B.icode) {}

   OperatorInfo& operator=(const OperatorInfo& B)
    {icode=B.icode; return *this;}

   OperatorInfo(const BasicLapHOperatorInfo& B) 
         : icode(B.icode) {}

   OperatorInfo(const GenIrrepOperatorInfo& B) 
         : icode(B.icode) {}

   ~OperatorInfo(){}


    // output functions
    
   bool isBasicLapH() const;
   
   bool isGenIrrep() const;

   bool isBackwards() const;

   void setBackwards();

   void setForwards();

   BasicLapHOperatorInfo getBasicLapH() const;  // throws if not basic lapH

   GenIrrepOperatorInfo getGenIrrep() const;    // throws if not gen irrep

   OperatorInfo& resetGenIrrepIDIndex(uint index);


   std::string output(bool longform=false, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=false) const;  // XML output

   std::string short_output() const;

   bool operator==(const OperatorInfo& rhs) const;

   bool operator!=(const OperatorInfo& rhs) const;

   bool operator<(const OperatorInfo& rhs) const;


   friend class MCObsInfo;
   friend class CorrelatorInfo;
   friend class CorrelatorAtTimeInfo;

 private:

   OperatorInfo(std::vector<unsigned int>::const_iterator inbegin,
                std::vector<unsigned int>::const_iterator inend) 
          : icode(inbegin,inend) {}


};


// **************************************************
#endif  
