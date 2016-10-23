#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;


class Base
{
 protected:

  int m_base_val;

 public:

  Base() : m_base_val(0) {}
  Base(int v) : m_base_val(v) {}

  virtual int getBValue() const = 0;
  virtual ~Base() {}
};

class Derived : public Base
{
 
  double m_derived_val;

 public:

  Derived() : m_derived_val(0), Base(0) {}
  Derived(double f) : m_derived_val(f), Base(0) {}
  Derived(double f, int v) : m_derived_val(f), Base(v) {}

  virtual ~Derived() {}

  double getDValue() const {return m_derived_val;}
  virtual int getBValue() const {return m_base_val;}
};

class Derived2 : public Base
{
 
  double m_derived_val;

 public:

  Derived2() : m_derived_val(0), Base(0) {}
  Derived2(double f) : m_derived_val(f), Base(0) {}
  Derived2(double f, int v) : m_derived_val(f), Base(v) {}

  virtual ~Derived2() {}

  double getDValue() const {return m_derived_val;}
  virtual int getBValue() const {return m_base_val;}
};

int main(){

Derived* p=new Derived(5.5,2);
cout << p->getDValue()<<endl;
cout << p->getBValue()<<endl;

Base* pp=p;
//cout << pp->getDValue()<<endl;
cout << pp->getBValue()<<endl;

Base* q=new Derived(8.8,3);
cout << q->getBValue()<<endl;

Derived2* qq=dynamic_cast<Derived*>(q);
cout << qq->getDValue()<<endl;


return 0;
}
