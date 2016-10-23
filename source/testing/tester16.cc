#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
using namespace std;

typedef unsigned int   uint;

      //  Extract a value from a string, throwing an exception on failure.

template <typename T>
void primitive_extract_from_string(const std::string& pstring, T& result, 
                                   const char* ptype) 
{
 try{
    std::istringstream is(pstring);
    std::ios::iostate state=is.exceptions();
    is.exceptions(std::ios::failbit);
    is.setf(std::ios_base::boolalpha);
       // try to read the type from the istringstream. 
       //   bool's should be "true" or "false" strings 
    is  >> result;
       // turn off exceptions on failure
    is.exceptions(state);
       // look for non white spaces in remainder of stream
    char c;
    while (is.get(c),!is.eof()){
       if (!isspace(c)) throw(std::invalid_argument("err"));}}
 catch(const std::exception& xp){ 
    std::ostringstream err;
    err << "Error: Failed to convert string: \"" 
        << pstring << "\" to " << ptype;
    throw(std::invalid_argument(err.str().c_str()));}
}

inline void extract_from_string(const std::string& pstring, std::string& result) 
{ result=pstring;}

inline void extract_from_string(const std::string& pstring, int& result) 
{ primitive_extract_from_string<int>(pstring, result, "int");}

inline void extract_from_string(const std::string& pstring, unsigned int& result) 
{ primitive_extract_from_string<unsigned int>(pstring, result, "unsigned int");}
    
inline void extract_from_string(const std::string& pstring, short int& result) 
{ primitive_extract_from_string<short int>(pstring, result,"short int");}

inline void extract_from_string(const std::string& pstring, unsigned short int& result) 
{ primitive_extract_from_string<unsigned short int>(pstring, result, "unsigned short int");}

inline void extract_from_string(const std::string& pstring, long int& result) 
{ primitive_extract_from_string<long int>(pstring, result, "long int");}
    
inline void extract_from_string(const std::string& pstring, unsigned long int& result) 
{ primitive_extract_from_string<unsigned long int>(pstring, result, "unsigned long int");}
    
inline void extract_from_string(const std::string& pstring, float& result) 
{ primitive_extract_from_string<float>(pstring, result, "float");}

inline void extract_from_string(const std::string& pstring, double& result) 
{ primitive_extract_from_string<double>(pstring, result, "double");}

inline void extract_from_string(const std::string& pstring, bool& result) 
{ primitive_extract_from_string<bool>(pstring, result, "bool");}



bool parse(const string& idstr)
{
 vector<size_t> pos;
 for (uint k=0;k<idstr.length();k++){
    if (idstr[k]=='|') pos.push_back(k);}
 if (pos.size()!=7) return false;
 unsigned int t1,t2,t3,t4,t5,t6,t7;
 try{
 extract_from_string(idstr.substr(pos[0]+1,pos[1]-pos[0]-1),t1);
 extract_from_string(idstr.substr(pos[1]+1,pos[2]-pos[1]-1),t2);
 extract_from_string(idstr.substr(pos[2]+1,pos[3]-pos[2]-1),t3);
 extract_from_string(idstr.substr(pos[3]+1,pos[4]-pos[3]-1),t4);
 extract_from_string(idstr.substr(pos[4]+1,pos[5]-pos[4]-1),t5);
 extract_from_string(idstr.substr(pos[5]+1,pos[6]-pos[5]-1),t6);
 extract_from_string(idstr.substr(pos[6]+1,string::npos),t7);
 string res;
 res=idstr.substr(0,pos[0]);}
 catch(const std::exception& xp){
   return false;}
 return true;
}

int main(){

bool flag=parse("abcdef|01|234|5|678|9a|01|123");

cout << "flag = "<<flag<<endl;

return 0;
}
