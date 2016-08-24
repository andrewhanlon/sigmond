#include "testing.h"
#include "multi_compare.h"
using namespace std;


   // check=0 if a==b, check=1 if a<b, check=2 if a>b
void doamulticompare(const std::vector<int>& a, const std::vector<int>& b,
                     int check, unsigned int& successcount,
                     unsigned int& failcount)
{
 bool lessthan=(check==1);
 bool greaterthan=(check==2);
 bool equal=(check==0);
 bool notequal=!equal;
 
 bool flag=multiLessThan(a,b);
 if (flag==lessthan) successcount++; else failcount++;
 flag=multiLessThan(b,a);
 if (flag==greaterthan) successcount++; else failcount++;
 flag=multiEqual(a,b);
 if (flag==equal) successcount++; else failcount++;
 flag=multiNotEqual(a,b);
 if (flag==notequal) successcount++; else failcount++;
}

void doamulticompare2(const std::vector<int>& a, const std::vector<int>& b,
                      unsigned int& successcount,
                      unsigned int& failcount)
{
 if ((a.size()<7)||(b.size()<7)){
    failcount++;
    return;}

 vector<int> u,v;
 int ndim=2;
 u.resize(ndim); v.resize(ndim);
 for (int k=0;k<ndim;++k){
    u[k]=a[k]; v[k]=b[k];}

 if (multiLessThan(u[0],v[0],u[1],v[1])==multiLessThan(u,v)) successcount++;
 else failcount++;
 if (   multiEqual(u[0],v[0],u[1],v[1])==multiEqual(u,v)) successcount++;
 else failcount++;
 if (multiNotEqual(u[0],v[0],u[1],v[1])==multiNotEqual(u,v)) successcount++;
 else failcount++;

 ndim=3;
 u.resize(ndim); v.resize(ndim);
 for (int k=0;k<ndim;++k){
    u[k]=a[k]; v[k]=b[k];}
 if (multiLessThan(u[0],v[0],u[1],v[1],u[2],v[2])==multiLessThan(u,v)) successcount++;
 else failcount++;
 if (   multiEqual(u[0],v[0],u[1],v[1],u[2],v[2])==multiEqual(u,v)) successcount++;
 else failcount++;
 if (multiNotEqual(u[0],v[0],u[1],v[1],u[2],v[2])==multiNotEqual(u,v)) successcount++;
 else failcount++;

 ndim=4;
 u.resize(ndim); v.resize(ndim);
 for (int k=0;k<ndim;++k){
    u[k]=a[k]; v[k]=b[k];}
 if (multiLessThan(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3])==multiLessThan(u,v)) successcount++;
 else failcount++;
 if (   multiEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3])==multiEqual(u,v)) successcount++;
 else failcount++;
 if (multiNotEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3])==multiNotEqual(u,v)) successcount++;
 else failcount++;

 ndim=5;
 u.resize(ndim); v.resize(ndim);
 for (int k=0;k<ndim;++k){
    u[k]=a[k]; v[k]=b[k];}
 if (multiLessThan(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4])==multiLessThan(u,v)) successcount++;
 else failcount++;
 if (   multiEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4])==multiEqual(u,v)) successcount++;
 else failcount++;
 if (multiNotEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4])==multiNotEqual(u,v)) successcount++;
 else failcount++;

 ndim=6;
 u.resize(ndim); v.resize(ndim);
 for (int k=0;k<ndim;++k){
    u[k]=a[k]; v[k]=b[k];}
 if (multiLessThan(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4],u[5],v[5])==multiLessThan(u,v)) successcount++;
 else failcount++;
 if (   multiEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4],u[5],v[5])==multiEqual(u,v)) successcount++;
 else failcount++;
 if (multiNotEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4],u[5],v[5])==multiNotEqual(u,v)) successcount++;
 else failcount++;

 ndim=7;
 u.resize(ndim); v.resize(ndim);
 for (int k=0;k<ndim;++k){
    u[k]=a[k]; v[k]=b[k];}
 if (multiLessThan(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4],u[5],v[5],u[6],v[6])==multiLessThan(u,v)) successcount++;
 else failcount++;
 if (   multiEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4],u[5],v[5],u[6],v[6])==multiEqual(u,v)) successcount++;
 else failcount++;
 if (multiNotEqual(u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4],u[5],v[5],u[6],v[6])==multiNotEqual(u,v)) successcount++;
 else failcount++;
}

   // check=0 if a==b, check=1 if a<b, check=2 if a>b
void doamulticompare3(const std::set<int>& a, const std::set<int>& b,
                      int check, unsigned int& successcount,
                      unsigned int& failcount)
{
 bool lessthan=(check==1);
 bool greaterthan=(check==2);
 bool equal=(check==0);
 bool notequal=!equal;
 
 bool flag=multiLessThan(a,b);
 if (flag==lessthan) successcount++; else failcount++;
 flag=multiLessThan(b,a);
 if (flag==greaterthan) successcount++; else failcount++;
 flag=multiEqual(a,b);
 if (flag==equal) successcount++; else failcount++;
 flag=multiNotEqual(a,b);
 if (flag==notequal) successcount++; else failcount++;
}



void testmulticompare(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMultiCompare")==0)
 return;

 xml_in.seek_unique("TestMultiCompare");
 cout << "Beginning test of multi compare"<<endl<<endl;
 
 unsigned int successcount=0;
 unsigned int failcount=0;
 
 vector<int> a,b;
 doamulticompare(a,b,0,successcount,failcount);
  
 a.resize(1); a[0]=3;
 doamulticompare(a,b,2,successcount,failcount);
 doamulticompare(a,a,0,successcount,failcount);

 a.resize(2); a[0]=1; a[1]=-5;
 doamulticompare(a,b,2,successcount,failcount);
 doamulticompare(a,a,0,successcount,failcount);
 doamulticompare(b,a,1,successcount,failcount);

 b=a;
 doamulticompare(a,b,0,successcount,failcount);
 b[0]=7;
 doamulticompare(a,b,1,successcount,failcount);
 b[0]=-2;
 doamulticompare(a,b,2,successcount,failcount);

 b.resize(4); b[0]=0; b[1]=1; b[2]=3; b[3]=8;
 doamulticompare(a,b,1,successcount,failcount);
 doamulticompare(b,a,2,successcount,failcount);

 a=b;
 doamulticompare(a,b,0,successcount,failcount);
 b[2]=11;
 doamulticompare(a,b,1,successcount,failcount);
 b[0]=-3;
 doamulticompare(a,b,2,successcount,failcount);
 

 a.resize(7), b.resize(7);
 b[0]=0; b[1]=1; b[2]=3; b[3]=8; b[4]=5; b[5]=0; b[6]=12;
 for (a[0]=-3;a[0]<=4;a[0]++)
 for (a[1]=-1;a[1]<=1;a[1]++)
 for (a[2]=0;a[2]<=4;a[2]++)
 for (a[3]=-3;a[3]<=15;a[3]+=4)
 for (a[4]=4;a[4]<=6;a[4]++)
 for (a[5]=-1;a[5]<=1;a[5]++)
 for (a[6]=9;a[6]<=13;a[6]+=2){
    doamulticompare2(a,b,successcount,failcount);
    doamulticompare2(b,a,successcount,failcount);}

 if (allEqual(1,1,1)==true) successcount++; else failcount++;
 if (allEqual(2,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,1,2)==false) successcount++; else failcount++;
 if (allEqual(2,2,1)==false) successcount++; else failcount++;
    
 if (allEqual(1,1,1,1)==true) successcount++; else failcount++;
 if (allEqual(2,1,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,1,2,1)==false) successcount++; else failcount++;
 if (allEqual(2,2,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,2,2,1)==false) successcount++; else failcount++;
    
 if (allEqual(1,1,1,1,1)==true) successcount++; else failcount++;
 if (allEqual(2,1,1,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,1,2,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,2,1,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,2,2,1,1)==false) successcount++; else failcount++;
 if (allEqual(2,2,2,2,1)==false) successcount++; else failcount++;
 if (allEqual(2,2,2,1,2)==false) successcount++; else failcount++;

 if (allNotEqual(1,1,1)==false) successcount++; else failcount++;
 if (allNotEqual(2,1,1)==false) successcount++; else failcount++;
 if (allNotEqual(2,1,2)==false) successcount++; else failcount++;
 if (allNotEqual(2,2,1)==false) successcount++; else failcount++;
 if (allNotEqual(2,3,1)==true) successcount++; else failcount++;

 if (allNotEqual(1,1,1,1)==false) successcount++; else failcount++;
 if (allNotEqual(2,1,1,3)==false) successcount++; else failcount++;
 if (allNotEqual(2,1,2,4)==false) successcount++; else failcount++;
 if (allNotEqual(2,2,1,4)==false) successcount++; else failcount++;
 if (allNotEqual(2,3,1,4)==true) successcount++; else failcount++;

 if (allNotEqual(1,1,1,1,5)==false) successcount++; else failcount++;
 if (allNotEqual(2,1,1,3,5)==false) successcount++; else failcount++;
 if (allNotEqual(2,1,2,4,5)==false) successcount++; else failcount++;
 if (allNotEqual(2,2,1,4,5)==false) successcount++; else failcount++;
 if (allNotEqual(2,3,1,4,5)==true) successcount++; else failcount++;


 set<int> sa,sb;
 doamulticompare3(sa,sb,0,successcount,failcount);
  
 sa.insert(3);
 doamulticompare3(sa,sb,2,successcount,failcount);
 doamulticompare3(sa,sa,0,successcount,failcount);
 doamulticompare3(sb,sa,1,successcount,failcount);

 sa.insert(1); sa.insert(-5);
 doamulticompare3(sa,sb,2,successcount,failcount);
 doamulticompare3(sa,sa,0,successcount,failcount);
 doamulticompare3(sb,sa,1,successcount,failcount);

 sb=sa;
 doamulticompare3(sa,sb,0,successcount,failcount);
 sb.insert(7);
 doamulticompare3(sa,sb,1,successcount,failcount);
 doamulticompare3(sb,sa,2,successcount,failcount);
 sa.insert(12);
 doamulticompare3(sa,sb,2,successcount,failcount);

 sb.insert(1); sb.insert(8);
 doamulticompare3(sa,sb,1,successcount,failcount);
 doamulticompare3(sb,sa,2,successcount,failcount);

 sa=sb;
 doamulticompare3(sa,sb,0,successcount,failcount);
 sb.insert(0); sa.insert(9);
 doamulticompare3(sa,sb,2,successcount,failcount);
 sa.insert(-13); sb.insert(13);
 doamulticompare3(sa,sb,1,successcount,failcount);



 cout << endl;
 cout << " Number of successful tests = "<<successcount<<endl;
 cout << "     Number of FAILED tests = "<<failcount<<endl<<endl;    
}


// ***********************************************
