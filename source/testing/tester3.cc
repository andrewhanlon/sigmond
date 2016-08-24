#include <string>
#include <iostream>
using namespace std;

class Dummy
{
   int m_value;

  public:

   Dummy() : m_value(0) {}
   Dummy(int inval) : m_value(inval) {}
   ~Dummy() { cout << "Destructor called"<<endl;}

   int getValue() const {return m_value;}

};

void testfunc(int inval)
{
 Dummy A;
 Dummy B(inval);
 if (B.getValue()<A.getValue()) throw(string("THROW"));
 cout << "No throw"<<endl;
}


int main(){

cout << "This routine tests whether a destructor is called for local var"<<endl;
cout << "when an exception thrown in a subroutine"<<endl;

cout <<endl<<"Try 1"<<endl;
try{
   testfunc(5);}
catch(...){
   cout << "Caught exception"<<endl;}

cout <<endl<<"Try 2"<<endl;
try{
   testfunc(2);}
catch(...){
   cout << "Caught exception"<<endl;}

cout <<endl<<"Try 3"<<endl;
try{
   testfunc(-3);}
catch(...){
   cout << "Caught exception"<<endl;}

cout <<endl<<"Try 4"<<endl;
try{
   testfunc(15);}
catch(...){
   cout << "Caught exception"<<endl;}

cout <<endl<<"Try 5"<<endl;
try{
   testfunc(-2);}
catch(...){
   cout << "Caught exception"<<endl;}

cout <<endl<<"Done"<<endl;

return 0;
}
