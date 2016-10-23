#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <list>
using namespace std;

  // Trims leading and trailing white space, then checks
  // to make sure each character is alphanumeric, underscore, 
  // or period. If name is invalid, an empty string is returned.

string tidyName(const string& str)
{
 string tmp;
 size_t len=str.length();
 if (len==0)
    throw(string("Invalid tag name"));
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
 if (tmp.empty())
    throw(string("Invalid tag name"));
 if (tmp=="__sub__args___")
    throw(string("Cannot use __sub__args___ as tag name in ArgsHandler"));
 return tmp;
}


class ArgBase
{
 public:
   ArgBase(){}
   virtual ~ArgBase() {}
};

template <typename T>
class ArgPointer : public ArgBase
{
   T* data;

 public:

  ArgPointer(T* inptr) {data=inptr;}
  ~ArgPointer() {}
  virtual T* getPtr() {return data;}
  virtual const T* getConstPtr() { return data;}

};


//  An object of the class "XMLHandler" handles only textual information.
//  An object of the class "ArgsHandler" bundles together, gets, and
//  echos data (argments) of all different types.  It associated XML
//  tag names with pointers to the data.  The data is not allocated nor
//  destroyed.  The calling procedure must handle all of the memory
//  management.



class ArgsHandler
{

    struct ArgsItem
    {
      std::string name;
      ArgBase *ptr;
      bool optional;
    };

    std::list<ArgsItem> m_items;
    std::string m_rootname;
    bool m_enforceroot;

 public:

    ArgsHandler(const std::string& roottag, bool enforceroot=true)
     : m_rootname(tidyName(roottag)), m_enforceroot(enforceroot){}

    ~ArgsHandler(){clear();}

//    void addItem(const std::string& tagname, uint& inref, bool optional=false);
//    void addItem(const std::string& tagname, int& inref, bool optional=false);
//    void addItem(const std::string& tagname, std::string& inref, bool optional=false);
//    void addItem(const std::string& tagname, double& inref, bool optional=false);
//    void addItem(const std::string& tagname, bool& inref, bool optional=false);
//    void addItem(const std::string& tagname, ArgsHandler& subargs);

    template <typename T>
    void addItem(const std::string& tagname, T& inref, bool optional=false)
    {
     ArgsItem a;
     a.name=tidyName(tagname);
     a.ptr=new ArgPointer<T>(&inref);
     a.optional=optional;
     m_items.push_back(a);
     }

    void addItem(ArgsHandler& subhandler, bool optional=false)
    {
     addItem("__sub__args___",subhandler,optional);
    }

  //  void get(XMLHandler& xmlin);
  //  void output(XMLHandler& xmlecho); 
    void clear()
    {
     list<ArgsItem>::iterator it;
     for (it=m_items.begin();it!=m_items.end();it++)
        delete it->ptr;
     m_items.clear();
    }

};








int main(){

ArgsHandler gin("Root");
int k;
gin.addItem("first_child",k);
double g;
gin.addItem("second_child",g);
string str;
gin.addItem("third_child",str);

k=3;
g=21.76;
str="a string";

ArgsHandler g2("Filio");
double f=6.7;
g2.addItem("bango",f);
gin.addItem(g2);


return 0;
}
