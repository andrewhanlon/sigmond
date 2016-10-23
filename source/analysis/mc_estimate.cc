#include "mc_estimate.h"
using namespace std;


MCEstimate::MCEstimate() : m_store(6), m_mode(Bootstrap)
{}

MCEstimate::MCEstimate(SamplingMode inmode) 
           : m_store((inmode==Bootstrap)?6:3), m_mode(inmode)
{}

MCEstimate::MCEstimate(const MCEstimate& in) 
           : m_store(in.m_store), m_mode(in.m_mode)
{}

MCEstimate& MCEstimate::operator=(const MCEstimate& in)
{
 m_store=in.m_store;
 m_mode=in.m_mode;
 return *this;
}

double MCEstimate::getFullEstimate() const
{
 return m_store[0];
}

double MCEstimate::getAverageEstimate() const
{
 return m_store[1];
}

double MCEstimate::getSymmetricError() const
{
 return m_store[2];
}

double MCEstimate::getLowerConfLimit() const
{
 if (m_mode!=Bootstrap) 
    throw(std::invalid_argument("Bootstrap mode required for getDownError in MCEstimate")); 
 return m_store[3];
}

double MCEstimate::getMedian() const
{
 if (m_mode!=Bootstrap) 
    throw(std::invalid_argument("Bootstrap mode required for getMedian in MCEstimate")); 
 return m_store[4];
}

double MCEstimate::getUpperConfLimit() const
{
 if (m_mode!=Bootstrap) 
    throw(std::invalid_argument("Bootstrap mode required for getUpError in MCEstimate")); 
 return m_store[5];
}


bool MCEstimate::isJackknifeMode() const
{
 return (m_mode==Jackknife);
}


bool MCEstimate::isBootstrapMode() const
{
 return (m_mode==Bootstrap);
}



string MCEstimate::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}


string MCEstimate::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void MCEstimate::output(XMLHandler& xmlout) const
{
 xmlout.set_root("MCEstimate");
 if (m_mode==Jackknife)
    xmlout.put_child("ResamplingMode","Jackknife");
 else
    xmlout.put_child("ResamplingMode","Bootstrap");
 xmlout.put_child("FullEstimate",make_string(getFullEstimate()));
 xmlout.put_child("SampleAverage",make_string(getAverageEstimate()));
 xmlout.put_child("SymmetricError",make_string(getSymmetricError()));
 if (m_mode==Bootstrap){
   xmlout.put_child("LowerConfLimit",make_string(getLowerConfLimit()));
   xmlout.put_child("Median",make_string(getMedian()));
   xmlout.put_child("UpperConfLimit",make_string(getUpperConfLimit()));}
}


void MCEstimate::rescale(double r)
{
 for (uint k=0;k<m_store.size();k++)
    m_store[k]*=r;
}


void MCEstimate::jackassign(double full, double avg, double error)
{
 m_store.resize(3);
 m_store[0]=full;
 m_store[1]=avg;
 m_store[2]=error;
 m_mode=Jackknife;
}


void MCEstimate::bootassign(double full, double avg, double error,
                            double low, double med, double upp)
{
 m_store.resize(6);
 m_store[0]=full;
 m_store[1]=avg;
 m_store[2]=error;
 m_store[3]=low;
 m_store[4]=med;
 m_store[5]=upp;
 m_mode=Bootstrap;
}


// ********************************************************************


string SimpleMCEstimate::str(unsigned int nerr_digits) const
{
 ostringstream ch;
 long int exponent,ndec,errint,valint;
 double errv,valv;
 string res;

 if (nerr_digits<1) ndec=1;
 else ndec=nerr_digits;

 exponent=(int) floor(log10(m_stddev));
 exponent-=ndec-1;
     //   rescale so that error is now a number between 10^(nerr_digits-1)
     //   and 10^(nerr_digits), then round, and readjust in case rounding
     //   up increased number of digits
 errv=m_stddev*pow10(-exponent);
 errint=(long int) floor(errv+0.5);          
 if (errint>=(long int) pow10(ndec)){
    errint/=10;
    exponent++;}

 valv=m_mean*pow10(-exponent);
 valint=(long int) floor(std::abs(valv)+0.5);

 if (exponent>0){
    valint*=pow10(exponent);
    errint*=pow10(exponent);
    }
 if (exponent>=0){
      // make into a string
    ch << valint << "(" << errint << ")";
    res=ch.str();
    }
 else{
    ch << valint << "(" << errint << ")";
    res=ch.str();
    int i,pos=0;
    while (res[pos]!='(') pos++;  // point to "("
    for (i=pos;i<=-exponent;i++) res="0"+res; // pad with 0's on left
    pos=0;
    while (res[pos]!='(') pos++;  // point to "("
    pos+=exponent;
    res=res.insert(pos,1,'.');
    if (ndec+exponent>0){
       pos=0;
       while (res[pos]!=')') pos++;  // point to ")"
       pos+=exponent;
       res=res.insert(pos,1,'.');
       }
    }

    // include negative sign if required
 if ((m_mean<=0.0)&&(valint>0)) res="-"+res;
 return res;
}
