#ifndef MULTI_BOOLEAN_H
#define MULTI_BOOLEAN_H
#include <vector>
#include <set>
#include <cstring>

// **********************************************************************
// *                                                                    *
// *   Objects used as keys in C++ maps require a less than             *
// *   operator.  The templated code here is useful for performing      *
// *   the less than operation on objects constructed from              *
// *   several other objects.  The equality and inequality operations   *
// *   are also given.                                                  *
// *                                                                    *
// *   Basically, two vector quantities u and v are compared,           *
// *   component by component.  u[0] and v[0] are compared first,       *
// *   then u[1] and v[1], and so on.   Many of the routines below      *
// *   require the vectors be input as   u[0], v[0], u[1], v[1], ...    *
// *                                                                    *
// **********************************************************************


template <typename T1, typename T2>
bool multiLessThan(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2)
{
 return   ((u1<v1) || ((u1==v1)
        && (u2<v2) ));
}

template <typename T1, typename T2>
bool multiEqual(const T1& u1, const T1& v1, 
                const T2& u2, const T2& v2)
{
 return (u1==v1)&&(u2==v2);
}
 
template <typename T1, typename T2>
bool multiNotEqual(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2)
{
 return (u1!=v1)||(u2!=v2);
}


// ******************************************************************


template <typename T1, typename T2, typename T3>
bool multiLessThan(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3)
{
 return   ((u1<v1) || ((u1==v1) 
       && ((u2<v2) || ((u2==v2)
       &&  (u3<v3) ))));
}

template <typename T1, typename T2, typename T3>
bool multiEqual(const T1& u1, const T1& v1, 
                const T2& u2, const T2& v2,
                const T3& u3, const T3& v3)
{
 return (u1==v1)&&(u2==v2)&&(u3==v3);
}
 
template <typename T1, typename T2, typename T3>
bool multiNotEqual(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3)
{
 return (u1!=v1)||(u2!=v2)||(u3!=v3);
}

// ******************************************************************


template <typename T1, typename T2, typename T3, typename T4>
bool multiLessThan(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4)
{
 return   ((u1<v1) || ((u1==v1) 
       && ((u2<v2) || ((u2==v2)
       && ((u3<v3) || ((u3==v3)
       &&  (u4<v4) ))))));
}


template <typename T1, typename T2, typename T3, typename T4>
bool multiEqual(const T1& u1, const T1& v1, 
                const T2& u2, const T2& v2,
                const T3& u3, const T3& v3,
                const T4& u4, const T4& v4)
{
 return (u1==v1)&&(u2==v2)&&(u3==v3)&&(u4==v4);
}

 
template <typename T1, typename T2, typename T3, typename T4>
bool multiNotEqual(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4)
{
 return (u1!=v1)||(u2!=v2)||(u3!=v3)||(u4!=v4);
}


// ******************************************************************


template <typename T1, typename T2, typename T3, typename T4, typename T5>
bool multiLessThan(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4,
                   const T5& u5, const T5& v5)
{
 return   ((u1<v1) || ((u1==v1) 
       && ((u2<v2) || ((u2==v2)
       && ((u3<v3) || ((u3==v3)
       && ((u4<v4) || ((u4==v4)
       &&  (u5<v5) ))))))));
}


template <typename T1, typename T2, typename T3, typename T4, typename T5>
bool multiEqual(const T1& u1, const T1& v1, 
                const T2& u2, const T2& v2,
                const T3& u3, const T3& v3,
                const T4& u4, const T4& v4,
                const T5& u5, const T5& v5)
{
 return (u1==v1)&&(u2==v2)&&(u3==v3)&&(u4==v4)&&(u5==v5);
}

 
template <typename T1, typename T2, typename T3, typename T4, typename T5>
bool multiNotEqual(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4,
                   const T5& u5, const T5& v5)
{
 return (u1!=v1)||(u2!=v2)||(u3!=v3)||(u4!=v4)||(u5!=v5);
}


// ******************************************************************


template <typename T1, typename T2, typename T3, typename T4, 
          typename T5, typename T6>
bool multiLessThan(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4,
                   const T5& u5, const T5& v5,
                   const T6& u6, const T6& v6)
{
 return   ((u1<v1) || ((u1==v1) 
       && ((u2<v2) || ((u2==v2)
       && ((u3<v3) || ((u3==v3)
       && ((u4<v4) || ((u4==v4)
       && ((u5<v5) || ((u5==v5)
       &&  (u6<v6) ))))))))));
}


template <typename T1, typename T2, typename T3, typename T4, 
          typename T5, typename T6>
bool multiEqual(const T1& u1, const T1& v1, 
                const T2& u2, const T2& v2,
                const T3& u3, const T3& v3,
                const T4& u4, const T4& v4,
                const T5& u5, const T5& v5,
                const T6& u6, const T6& v6)
{
 return (u1==v1)&&(u2==v2)&&(u3==v3)&&(u4==v4)&&(u5==v5)&&(u6==v6);
}

 
template <typename T1, typename T2, typename T3, typename T4, 
          typename T5, typename T6>
bool multiNotEqual(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4,
                   const T5& u5, const T5& v5,
                   const T6& u6, const T6& v6)
{
 return (u1!=v1)||(u2!=v2)||(u3!=v3)||(u4!=v4)||(u5!=v5)||(u6!=v6);
}


// ******************************************************************


template <typename T1, typename T2, typename T3, typename T4, 
          typename T5, typename T6, typename T7>
bool multiLessThan(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4,
                   const T5& u5, const T5& v5,
                   const T6& u6, const T6& v6,
                   const T7& u7, const T7& v7)
{
 return   ((u1<v1) || ((u1==v1) 
       && ((u2<v2) || ((u2==v2)
       && ((u3<v3) || ((u3==v3)
       && ((u4<v4) || ((u4==v4)
       && ((u5<v5) || ((u5==v5)
       && ((u6<v6) || ((u6==v6)
       &&  (u7<v7) ))))))))))));
}


template <typename T1, typename T2, typename T3, typename T4, 
          typename T5, typename T6, typename T7>
bool multiEqual(const T1& u1, const T1& v1, 
                const T2& u2, const T2& v2,
                const T3& u3, const T3& v3,
                const T4& u4, const T4& v4,
                const T5& u5, const T5& v5,
                const T6& u6, const T6& v6,
                const T7& u7, const T7& v7)
{
 return (u1==v1)&&(u2==v2)&&(u3==v3)&&(u4==v4)&&(u5==v5)&&(u6==v6)&&(u7==v7);
}

 
template <typename T1, typename T2, typename T3, typename T4, 
          typename T5, typename T6, typename T7>
bool multiNotEqual(const T1& u1, const T1& v1, 
                   const T2& u2, const T2& v2,
                   const T3& u3, const T3& v3,
                   const T4& u4, const T4& v4,
                   const T5& u5, const T5& v5,
                   const T6& u6, const T6& v6,
                   const T7& u7, const T7& v7)
{
 return (u1!=v1)||(u2!=v2)||(u3!=v3)||(u4!=v4)||(u5!=v5)||(u6!=v6)||(u7!=v7);
}


// ******************************************************************


template <typename T>
bool multiLessThan(const std::vector<T>& u, const std::vector<T>& v)
{
 if (u.size()<v.size()) return true;
 if (u.size()>v.size()) return false;
 for (size_t k=0;k<u.size();++k){
    if (u[k]<v[k]) return true;
    if (u[k]!=v[k]) return false;}
 return false;
}

template <typename T>
bool multiEqual(const std::vector<T>& u, const std::vector<T>& v)
{
 if (u.size()!=v.size()) return false;
 for (size_t k=0;k<u.size();++k){
    if (u[k]!=v[k]) return false;}
 return true;
}
 
template <typename T>
bool multiNotEqual(const std::vector<T>& u, const std::vector<T>& v)
{
 if (u.size()!=v.size()) return true;
 for (size_t k=0;k<u.size();++k){
    if (u[k]!=v[k]) return true;}
 return false;
}


template <typename T>
bool multiLessThan(const std::set<T>& u, const std::set<T>& v)
{
 if (u.size()<v.size()) return true;
 if (u.size()>v.size()) return false;
 typename std::set<T>::const_iterator ut,vt;
 for (ut=u.begin(),vt=v.begin();(ut!=u.end())&&(vt!=v.end());ut++,vt++){
    if (*ut<*vt) return true;
    if (*ut!=*vt) return false;}
 return false;
}

template <typename T>
bool multiEqual(const std::set<T>& u, const std::set<T>& v)
{
 if (u.size()!=v.size()) return false;
 typename std::set<T>::const_iterator ut,vt;
 for (ut=u.begin(),vt=v.begin();(ut!=u.end())&&(vt!=v.end());ut++,vt++){
    if (*ut!=*vt) return false;}
 return true;
}
 
template <typename T>
bool multiNotEqual(const std::set<T>& u, const std::set<T>& v)
{
 if (u.size()!=v.size()) return true;
 typename std::set<T>::const_iterator ut,vt;
 for (ut=u.begin(),vt=v.begin();(ut!=u.end())&&(vt!=v.end());ut++,vt++){
    if (*ut!=*vt) return true;}
 return false;
}


// ******************************************************************

template <typename T>
bool allEqual(const T& u1, const T& u2, const T& u3)
{
 return (u1==u2)&&(u2==u3);
}

template <typename T>
bool allEqual(const T& u1, const T& u2, const T& u3, const T& u4)
{
 return (u1==u2)&&(u2==u3)&&(u3==u4);
}

template <typename T>
bool allEqual(const T& u1, const T& u2, const T& u3, const T& u4, const T& u5)
{
 return (u1==u2)&&(u2==u3)&&(u3==u4)&&(u4==u5);
}

template <typename T>
bool allNotEqual(const T& u1, const T& u2, const T& u3)
{
 return (u1!=u2)&&(u1!=u3)&&(u2!=u3);
}

template <typename T>
bool allNotEqual(const T& u1, const T& u2, const T& u3, const T& u4)
{
 return (u1!=u2)&&(u1!=u3)&&(u1!=u4)&&(u2!=u3)&&(u2!=u4)&&(u3!=u4);
}

template <typename T>
bool allNotEqual(const T& u1, const T& u2, const T& u3, const T& u4, const T& u5)
{
 return (u1!=u2)&&(u1!=u3)&&(u1!=u4)&&(u1!=u5)&&(u2!=u3)&&(u2!=u4)&&(u2!=u5)
        &&(u3!=u4)&&(u3!=u5)&&(u4!=u5);
}

// ******************************************************************

#endif
