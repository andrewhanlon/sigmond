#include <string>
#include <iostream>
using namespace std;

/*
enum ComplexArg { RealPart, ImaginaryPart };

void secondary_encode(const std::string& instr, unsigned int index,
                      bool nonsimple, ComplexArg arg,
                      std::vector<unsigned int> icode)
{
}

void secondary_decode(const std::string& instr, unsigned int index,
                      bool nonsimple, ComplexArg arg,
                      std::vector<unsigned int> icode)
{
}
*/



unsigned int four_char_to_uint(const char *cstr)
{
 unsigned int code=unsigned int(cstr[0]);
 code<<=8;
 code|=unsigned int(cstr[1]);
 code<<=8;
 code|=unsigned int(cstr[2]);
 code<<=8;
 code|=unsigned int(cstr[3]);
 return code;
}



int main(){

string istr("ThisIsAString0");

unsigned int n=istr.length();
cout << "length of string = "<<n<<endl;

cout << "size of int = "<<sizeof(unsigned int)<<endl;

unsigned int irem=n%sizeof(unsigned int);
cout << "irem = "<<irem<<endl;
if (irem>0){
   for (unsigned int k=0;k<(sizeof(unsigned int)-irem);k++)
      istr+=" ";}

cout << "istr = <"<<istr<<">"<<endl;

cout << "length of istr now = "<<istr.length()<<endl;
n=istr.length();

{union { char si[32]; unsigned int ui[8]; };

for (unsigned int k=0;k<n;++k)
   si[k]=istr[k];

for (unsigned int k=0;k<n;k++)
   cout << "si["<<k<<"] = "<<si[k]<<endl;

for (unsigned int k=0;k<n/4;k++)
   cout << "ui["<<k<<"] = "<<ui[k]<<endl;}

return 0;
}
