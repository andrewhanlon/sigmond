#ifndef MOMENTA_H
#define MOMENTA_H

#include "multi_compare.h"
#include <string>


// *******************************************************************
// *                                                                 *
// *   Defines a class "Momentum".   The class "Momentum" is just a  *
// *   convenient struct for storing a three-momentum defined as     *
// *   three integers (since momentum is quantized on a toroid).     *
// *                                                                 *
// *                                                                 *
// *******************************************************************


struct Momentum
{

 int x,y,z;

 Momentum(){}
 Momentum(int px, int py, int pz) : x(px), y(py), z(pz) {}
 Momentum(const Momentum& rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}
 Momentum& operator=(const Momentum& rhs)
  {x=rhs.x; y=rhs.y; z=rhs.z; return *this;}

 bool operator<(const Momentum& rhs) const
  {return multiLessThan(x,rhs.x,  y,rhs.y,  z,rhs.z);}
 bool operator==(const Momentum& rhs) const
  {return multiEqual(x,rhs.x,  y,rhs.y,  z,rhs.z);}
 bool operator!=(const Momentum& rhs) const
  {return multiNotEqual(x,rhs.x,  y,rhs.y,  z,rhs.z);}

   //   Returns true if (px,py,pz) is an allowed momentum ray:
   //     000  +00  0+0  00+  -00  0-0  00-
   //     ++0  +-0  +0+  +0-  0++  0+-
   //     --0  -+0  -0-  -0+  0--  0-+
   //     +++  ++-  +-+  +--  ---  --+  -+-  -++
   
 bool isAllowedRay(std::string& ray) const;

   //   Returns the momentum direction. This will either be
   //   the momentum itself, or the momentum rescaled.  For example,
   //   (0,0,2) will be rescaled to (0,0,1).

 Momentum getMomentumDirection() const;

   //   Returns the 3-character string associated with this
   //   momentum direction.  Recall that "+" refers to +1, "-" refers
   //   to -1, "#" refers to +2, and "=" refers to -2.
   //   The momentum is suitably rescaled.  For example,
   //   (0,0,2) will be rescaled to (0,0,1) -> "00+".

 std::string getMomentumType() const;
 
   //   Returns the 3-character string associated with this
   //   momentum direction.  Recall that "+" refers to +1, "-" refers
   //   to -1, "#" refers to +2, and "=" refers to -2.
   //   The momentum is NOT rescaled.  

 std::string getMomentumString() const;

 private:

 unsigned int get_gcd(unsigned int ia, unsigned int ib) const;
 
 const static char charcodes[5];
 
};


// **************************************************
#endif  
