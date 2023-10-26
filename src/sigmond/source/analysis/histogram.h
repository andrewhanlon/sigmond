#ifndef HISTOGRAM_H
#define HISTOGRAM_H
#include "matrix.h"
#include "scalar_defs.h"

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
                 double outlier_scale=5.0);



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
  
class Histogram {

  double m_lower, m_upper;
  uint m_nbars;
  Vector<uint> m_histo;

 public:

  Histogram(const Vector<double>& mcdata, double lowerlimit, 
            double upperlimit, uint nbars=40);
  Histogram(const Vector<double>& mcdata, uint nbars=40);
  Histogram(const Histogram& inhisto);
  Histogram& operator=(const Histogram& inhisto);

  double getRangeLowerLimit() const;
  double getRangeUpperLimit() const;
  double getBarWidth() const;
  uint getNumberOfBars() const;
  double getBarRangeLowerLimit(uint jbar) const;
  double getBarRangeUpperLimit(uint jbar) const;
  double getBarMiddleLocation(uint jbar) const;
  uint getBarHeight(uint jbar) const;
  uint getMaxHeight() const;

 private:

  void initialize(const Vector<double>& mcdata, double lowerlimit, 
                  double upperlimit, uint nbars, bool usedatalimits);

};



// **************************************************************
#endif
