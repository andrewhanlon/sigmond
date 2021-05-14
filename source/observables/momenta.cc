#include "momenta.h"
#include <cstdlib>
#include <cmath>

using namespace std;


   //   Allowed momentum rays:
   //     000  +00  0+0  00+  -00  0-0  00-
   //     ++0  +-0  +0+  +0-  0++  0+-
   //     --0  -+0  -0-  -0+  0--  0-+
   //     +++  ++-  +-+  +--  ---  --+  -+-  -++

bool Momentum::isAllowedRay(std::string& ray) const
{
 ray.clear();
 int X=(x<0)?-x:x;  // absolute values
 int Y=(y<0)?-y:y;
 int Z=(z<0)?-z:z;
 int n=0,nz[3];
 if (X>0) nz[n++]=X;
 if (Y>0) nz[n++]=Y;
 if (Z>0) nz[n++]=Z;
 if ((n==2)&&(nz[0]!=nz[1])) return false;
 if ((n==3)&&((nz[0]!=nz[1])||(nz[1]!=nz[2]))) return false;

 char dir[3]={'-','0','+'};
 X=(x>0)?1:((x<0)?-1:0);  // get zero or sign
 Y=(y>0)?1:((y<0)?-1:0);
 Z=(z>0)?1:((z<0)?-1:0);
 ++X; ++Y; ++Z;
 ray="mom_ray_"; ray+=dir[X]; ray+=dir[Y]; ray+=dir[Z];
 return true;
}


unsigned int Momentum::get_gcd(unsigned int ia, unsigned int ib) const
{
 unsigned int b=ib, a=ia, t;
 while (b!=0){
    t=b;
    b=a % t;
    a=t;}
 return a;
}

   //   Returns the momentum direction. This will either be
   //   the momentum itself, or the momentum rescaled.  For example,
   //   (0,0,2) will be rescaled to (0,0,1).

Momentum Momentum::getMomentumDirection() const
{
 int gcd=get_gcd(get_gcd(abs(x),abs(y)),abs(z));
 if (gcd<=1) return Momentum(x,y,z);
 else return Momentum(x/gcd,y/gcd,z/gcd);
}


   //   Returns the 3-character string associated with this
   //   momentum direction.  Recall that "+" refers to +1, "-" refers
   //   to -1, "#" refers to +2, and "=" refers to -2.
   //   The momentum is suitably rescaled.  For example,
   //   (0,0,2) will be rescaled to (0,0,1) -> "00+".
   //   A component value larger than 2 will be returned as "X".
   
std::string Momentum::getMomentumType() const
{
 Momentum dir=getMomentumDirection();
 char codestr[3];
 codestr[0]=((dir.x>=-2)&&(dir.x<=2))?charcodes[dir.x+2]:'X';
 codestr[1]=((dir.y>=-2)&&(dir.y<=2))?charcodes[dir.y+2]:'X';
 codestr[2]=((dir.z>=-2)&&(dir.z<=2))?charcodes[dir.z+2]:'X';
 return string(codestr,3);
}


   //   Returns the 3-character string associated with this
   //   momentum direction.  Recall that "+" refers to +1, "-" refers
   //   to -1, "#" refers to +2, and "=" refers to -2.
   //   The momentum is NOT rescaled.  
   //   A component value larger than 2 will be returned as "X".
   
std::string Momentum::getMomentumString() const
{
 char codestr[3];
 codestr[0]=((x>=-2)&&(x<=2))?charcodes[x+2]:'X';
 codestr[1]=((y>=-2)&&(y<=2))?charcodes[y+2]:'X';
 codestr[2]=((z>=-2)&&(z<=2))?charcodes[z+2]:'X';
 return string(codestr,3);
}

uint Momentum::getPsq() const
{
  return pow(x, 2.) + pow(y, 2.) + pow(z, 2.);
}

double Momentum::getPsq(int spatial_extent) const
{
  return pow(2.*M_PI, 2.)*(pow(double(x)/spatial_extent, 2.) + pow(double(y)/spatial_extent, 2.) + pow(double(z)/spatial_extent, 2.));
}

double Momentum::getPsq(int x_extent, int y_extent, int z_extent) const
{
  return pow(2.*M_PI, 2.)*(pow(double(x)/x_extent, 2.) + pow(double(y)/y_extent, 2.) + pow(double(z)/z_extent, 2.));
}

const char Momentum::charcodes[5]={'=','-','0','+','#'};

// **************************************************
