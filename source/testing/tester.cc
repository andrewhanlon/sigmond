#include <vector>
#include <iostream>
#include <list>
#include <complex>

using namespace std;

double conjugate(const double& x)
{
 return x;
}

float conjugate(const float& x)
{
 return x;
}

std::complex<double> conjugate(const std::complex<double>& z)
{
 return conj(z);
}

std::complex<float> conjugate(const std::complex<float>& z)
{
 return conj(z);
}

double realpart(const double& x)
{
 return x;
}

double realpart(const float& x)
{
 return double(x);
}

double imaginarypart(const double& x)
{
 return 0.0;
}

double imaginarypart(const float& x)
{
 return 0.0;
}

double realpart(const complex<double>& z)
{
 return real(z);
}

double realpart(const complex<float>& z)
{
 return double(real(z));
}

double imaginarypart(const complex<double>& z)
{
 return imag(z);
}

double imaginarypart(const complex<float>& z)
{
 return double(imag(z));
}





int main(){

vector<int> a(73890000);
for (int i=0;i<a.size();i++) a[i]=i;
//for (int i=0;i<a.size();i++) a[i]*=2;
//for (int i=0;i<a.size();i++) a[i]+=5;
//for (int i=0;i<a.size();i++) a[i]*=3;
//for (int i=0;i<a.size();i++) a[i]-=7;

//for (vector<int>::iterator it=a.begin();it!=a.end();it++) *it*=2;
//for (vector<int>::iterator it=a.begin();it!=a.end();it++) *it+=5;
//for (vector<int>::iterator it=a.begin();it!=a.end();it++) *it*=3;
//for (vector<int>::iterator it=a.begin();it!=a.end();it++) *it-=7;

vector<int> b(a.begin()+88,a.begin()+92);
for (int k=0;k<b.size();k++)
   cout <<" b["<<k<<"] = "<<b[k]<<endl;


list<int> ga;
//ga.push_back(1);
//ga.push_back(2);
list<int> gb;
//gb.push_back(3);
//gb.push_back(4);

ga.splice(ga.end(),gb);

cout << "ga list:";
for (list<int>::const_iterator it=ga.begin();it!=ga.end();it++) cout <<"  "<<*it;
cout << endl;
cout << "gb list:";
for (list<int>::const_iterator it=gb.begin();it!=gb.end();it++) cout <<"  "<<*it;
cout << endl;


double *ptr=0;//new double(5);

//cout << *ptr<<endl;

delete ptr;



complex<double> z(1.3,5.2);
cout << realpart(z)<<endl;
cout << imaginarypart(z) <<endl;


double x=5.6;
cout << realpart(x)<<endl;
cout << imaginarypart(x)<<endl;


return 0;
}
