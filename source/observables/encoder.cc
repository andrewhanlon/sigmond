#include "encoder.h"
#include <stdexcept>

using namespace std;

typedef unsigned int  uint;


void encode_string_to_uints(const string& astr, unsigned int maxchar,
                            vector<unsigned int>& icode)
{
 uint nchar=astr.length();
 if (nchar>maxchar){
    throw(std::invalid_argument("String for encoding is too long"));}
 if (nchar==0){
    throw(std::invalid_argument("String for encoding cannot be empty"));}
 for (uint k=0;k<nchar;k++){
    if (isspace(astr[k]))
       throw(std::invalid_argument("String for encoding cannot contain white space"));}
 uint nchar_per_int=sizeof(unsigned int); // should be four
 string buffer(astr);
 uint npad=nchar%nchar_per_int;
 if (npad>0){
    npad=nchar_per_int-npad;
    buffer+=string(npad,' ');}
 nchar=buffer.length();
 uint nint=nchar/nchar_per_int;
 icode.resize(nint);
 uint ic=0;
 for (uint ind=0;ind<nint;ind++){
   uint tt=(unsigned int)(buffer[ic++]);
   tt<<=8; tt|=(unsigned int)(buffer[ic++]);
   tt<<=8; tt|=(unsigned int)(buffer[ic++]); 
   tt<<=8; tt|=(unsigned int)(buffer[ic++]);
   icode[ind]=tt;}
}



string decode_uints_to_string(const vector<unsigned int>& icode)
{
 uint nint=icode.size();
 uint nchar_per_int=sizeof(unsigned int); // should be four
 vector<char> buf(nint*nchar_per_int);
 for (uint ind=0;ind<nint;ind++){
    uint tcode=icode[ind];
    uint ic=(ind+1)*nchar_per_int;
    buf[--ic]=tcode&255u;
    tcode>>=8; buf[--ic]=tcode&255u; 
    tcode>>=8; buf[--ic]=tcode&255u; 
    tcode>>=8; buf[--ic]=tcode&255u;}
  string astr(buf.begin(),buf.end());
  std::size_t found=astr.find_last_not_of(" ");
  if (found!=std::string::npos)
     astr.erase(found+1);
 return astr;
}


// *************************************************************
