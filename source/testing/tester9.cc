#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;

  // Trims leading and trailing white space, then checks
  // to make sure each character is alphanumeric, underscore, 
  // or period. If name is invalid, an empty string is returned.

string tidyName(const string& str)
{
 string tmp;
 size_t len=str.length();
 if (len==0) return tmp;
 size_t start=0;
 while ((isspace(str[start]))&&(start<len)) start++;
 if (start==len) return tmp;
 size_t stop=len-1;
 while ((isspace(str[stop]))&&(stop>start)) stop--;
 for (size_t i=start;i<=stop;i++){
    char c=str[i];
    if (isalnum(c)||(c=='_')||(c=='.'))
       tmp.push_back(c);
    else return string("");}
 return tmp;
}


int main(){

cout << "|"<<tidyName("aname_rabbit")<<"|"<<endl;
cout << "|"<<tidyName("aname_rab88bit")<<"|"<<endl;
cout << "|"<<tidyName("  aname5_rabbit")<<"|"<<endl;
cout << "|"<<tidyName("aname_rabbit  ")<<"|"<<endl;
cout << "|"<<tidyName("   aname.rabbit  ")<<"|"<<endl;
cout << "|"<<tidyName("   aname./rabbit  ")<<"|"<<endl;
cout << "|"<<tidyName("    ")<<"|"<<endl;
cout << "|"<<tidyName("")<<"|"<<endl;
cout << "|"<<tidyName("  \t  \n aname_\rweb\n\t")<<"|"<<endl;
cout << "|"<<tidyName("  \t  \n aname_rabbit\n\t\r  ")<<"|"<<endl;

return 0;
}
