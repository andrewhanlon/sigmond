#include "histogram.h"
#include <algorithm>
#include <list>
using namespace std;

// *****************************************************************
// *                                                               *
// *   Given a vector of Monte Carlo data in "mcdata" (either the  *
// *   real parts or the imaginary parts to produce double         *
// *   precision data), this routine returns in "outliers" a       *
// *   vector of indices corresponding to data points that are     *
// *   significantly outside of the range of the data.  This helps *
// *   to identify Monte Carlo measurements that may have been     *
// *   corrupted.  Outliers are identified as follows: the data is *
// *   sorted, then the data range "a..b" for the middle half of   *
// *   the data points is found.  Let the mid-range be "m=(a+b)/2" *
// *   then outliers are points outside the range  "m-v .. m+v"    *
// *   where "v = outlier_scale * (b-a)/2 ".  By default,          *
// *   "outlier_scale" is input as 5.0, but this can be changed    *
// *   by the calling code.                                        *
// *                                                               *
// *****************************************************************
  

void getOutliers(const Vector<double>& mcdata, Vector<uint>& outliers,
                 double outlier_scale)
{
 outliers.clear();
 uint ndata=mcdata.size();
 if (ndata<32) return;
 vector<double> temp(mcdata.c_vector());
 std::sort(temp.begin(),temp.end());
 double lo=temp[ndata/4]; 
 double hi=temp[3*ndata/4];
 double m=(lo+hi)/2.0;
 double v=outlier_scale*(hi-lo)/2.0;
 lo=m-v;
 hi=m+v;
 std::list<uint> ol;
 for (uint k=0;k<ndata;k++)
   if ((mcdata[k]<lo)||(mcdata[k]>hi)) ol.push_back(k);
 if (ol.empty()) return;
 outliers.resize(uint(ol.size()));
 uint k=0;
 for (list<uint>::const_iterator it=ol.begin();it!=ol.end();it++,k++)
    outliers[k]=*it;
}


// *****************************************************************
// *                                                               *
// *   This class takes a vector of Monte Carlo data in "mcdata"   *
// *   (either the real parts or the imaginary parts to produce    *
// *   double precision data) and computes a histogram of the      *
// *   results.  The number of histogram bars can be chosen,       *
// *   as well as the lower and upper limits of a range of values  *
// *   to focus on.  The bar index is zero offset.                 *
// *                                                               *
// *****************************************************************
  
Histogram::Histogram(const Vector<double>& mcdata, double lowerlimit, 
                     double upperlimit, uint nbars)
{
 initialize(mcdata,lowerlimit,upperlimit,nbars,false);
}

Histogram::Histogram(const Vector<double>& mcdata, uint nbars)
{
 initialize(mcdata,0.0,0.0,nbars,true);
}

Histogram::Histogram(const Histogram& inhisto) : 
                    m_lower(inhisto.m_lower),  m_upper(inhisto.m_upper),
                    m_nbars(inhisto.m_nbars),  m_histo(inhisto.m_histo)
{}

Histogram& Histogram::operator=(const Histogram& inhisto)
{
 m_lower=inhisto.m_lower, 
 m_upper=inhisto.m_upper;
 m_nbars=inhisto.m_nbars;
 m_histo=inhisto.m_histo;
 return *this;
}

double Histogram::getRangeLowerLimit() const
{
 return m_lower;
}

double Histogram::getRangeUpperLimit() const
{
 return m_upper;
}

double Histogram::getBarWidth() const
{
 return (m_upper-m_lower)/m_nbars;
}

uint Histogram::getNumberOfBars() const
{
 return m_nbars;
}

double Histogram::getBarRangeLowerLimit(uint jbar) const
{
 if (jbar>=m_nbars) throw(std::invalid_argument("Histogram:getBarRangeLowerLimit index out of range"));
 return m_lower+((m_upper-m_lower)/m_nbars)*jbar;
}

double Histogram::getBarRangeUpperLimit(uint jbar) const
{
 if (jbar>=m_nbars) throw(std::invalid_argument("Histogram:getBarRangeUpperLimit index out of range"));
 if (jbar==(m_nbars-1)) return m_upper;
 return m_lower+((m_upper-m_lower)/m_nbars)*(jbar+1);
}

double Histogram::getBarMiddleLocation(uint jbar) const
{
 if (jbar>=m_nbars) throw(std::invalid_argument("Histogram:getBarMiddleLocation index out of range"));
 return m_lower+((m_upper-m_lower)/m_nbars)*(double(jbar)+0.5);
}

uint Histogram::getBarHeight(uint jbar) const
{
 if (jbar>=m_nbars) throw(std::invalid_argument("Histogram:getBarHeight index out of range"));
 return m_histo[jbar];
}


uint Histogram::getMaxHeight() const
{
 uint hmax=0;
 for (uint jbar=0;jbar<m_nbars;jbar++)
    if (hmax<m_histo[jbar]) hmax=m_histo[jbar];
 return hmax;
}


void Histogram::initialize(const Vector<double>& mcdata, double lowerlimit, 
                           double upperlimit, uint nbars, bool usedatalimits)
{
 if (mcdata.size()<12) throw(std::invalid_argument("Too little data for Histogram"));
 if (nbars<4) throw(std::invalid_argument("At least 4 bars needed for Histogram"));
 if (nbars>65536) throw(std::invalid_argument("Too many bars requested for Histogram"));
 m_nbars=nbars;
 if (usedatalimits){
    m_lower=mcdata[0];
    m_upper=mcdata[0];
    for (uint k=1;k<mcdata.size();k++){
       if (mcdata[k]<m_lower) m_lower=mcdata[k];
       if (mcdata[k]>m_upper) m_upper=mcdata[k];}
    double barwidth=(m_upper-m_lower)/m_nbars;
    if (barwidth==0.0) barwidth=0.1;
    m_lower-=0.25*barwidth;
    m_upper+=0.25*barwidth;}
 else{
    if (lowerlimit>=upperlimit) throw(std::invalid_argument("Invalid range for Histogram"));
    m_lower=lowerlimit;
    m_upper=upperlimit;}

 vector<double> temp(mcdata.c_vector());
 std::sort(temp.begin(),temp.end());

 m_histo.resize(m_nbars);
 for (uint k=0;k<m_nbars;k++) m_histo[k]=0;

 vector<double>::const_iterator it=temp.begin();
 while ((it!=temp.end())&&(*it<m_lower)) ++it;

 double barwidth=(m_upper-m_lower)/m_nbars;
 double barlo=m_lower;
 double barup=barlo+barwidth;
 for (uint k=0;k<m_nbars;k++){
    while ((it!=temp.end())&&(*it>=barlo)&&(*it<barup)){
       m_histo[k]++; it++;}
    barlo=barup;
    barup=barlo+barwidth;}

}

// **************************************************************
