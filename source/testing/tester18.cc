#include <iostream>
using namespace std;

class Dummy
{
  int m_value;

  Dummy();  // no default

 public:

  Dummy(int k) : m_value(k) {}

  int getValue() const {return m_value;}

  static int numints() { return 1;}

};

int main()
{

 Dummy x(5);

 cout << x.getValue()<<endl;

 cout << Dummy::numints()<<endl;

 cout << x.numints()<<endl;

 return 0;
}
