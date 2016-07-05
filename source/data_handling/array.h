#ifndef ARRAY_H
#define ARRAY_H
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>

// *************************************************************
// *                                                           *
// *   The class "Array" is defined in this file.  The data    *
// *   is stored on disk in terms of "Array" structures, so    *
// *   the IO handlers/maps need access to this class.         *
// *                                                           *
// *************************************************************


namespace LaphEnv {

// **************************************************************************
// *                                                                        *
// *  Class "Array" is a simple container for an N-dimensional array.  The  *
// *  dimension is not templated so that different dimensions are handled   *
// *  by a single type. This is useful as a data type in an operator end    *
// *  line since the dimension needs to be different depending on how       *
// *  many quark lines are always contracted.  Storage is column-major      *
// *  (indices on left are fastest varying).  Zero-offset is used.          *
// *                                                                        *
// **************************************************************************

  // Use the #define below to include range checking code.
  // Do not define and range checking code is omitted (faster).
//#define SAFETY_FLAG



template <typename T>
class Array {

   std::vector<unsigned int> m_sizes;
   std::vector<T> m_store;

 public:

   Array();
   Array(int size1);
   Array(int size1, int size2);
   Array(int size1, int size2, int size3);
   Array(int size1, int size2, int size3, int size4);
   Array(int size1, int size2, int size3, int size4, 
         int size5);
   Array(int size1, int size2, int size3, int size4, 
         int size5, int size6);
   Array(int size1, int size2, int size3, int size4, 
         int size5, int size6, int size7);
   Array(int size1, int size2, int size3, int size4, 
         int size5, int size6, int size7, int size8);
   Array(const std::vector<unsigned int>& sizevals);
   Array(const std::vector<int>& sizevals);
   Array(const std::vector<unsigned int>& sizevals, const T& initial_value);
   Array(const std::vector<int>& sizevals, const T& initial_value);
   Array(const Array<T>& incoming);


   ~Array();
   Array<T>& clear(); 

   Array<T>& operator=(const T& val);
   Array<T>& operator+=(const T& val);
   Array<T>& operator-=(const T& val);
   Array<T>& operator*=(const T& val);
   Array<T>& operator/=(const T& val);
   Array<T>& operator=(const Array<T>& incoming); 
   Array<T>& operator+=(const Array<T>& incoming);
   Array<T>& operator-=(const Array<T>& incoming);
   Array<T>& operator*=(const Array<T>& incoming);
   Array<T>& operator/=(const Array<T>& incoming);

   T& operator()(int i);             
   const T& operator()(int i) const; 
   T& operator[](int i);             
   const T& operator[](int i) const; 
   T& operator()(int i0, int i1);    
   const T& operator()(int i0, int i1) const; 
   T& operator()(int i0, int i1, int i2);     
   const T& operator()(int i0, int i1, int i2) const;
   T& operator()(int i0, int i1, int i2, int i3);  
   const T& operator()(int i0, int i1, int i2, int i3) const; 
   T& operator()(int i0, int i1, int i2, int i3,
                 int i4);  
   const T& operator()(int i0, int i1, int i2, int i3,
                       int i4) const; 
   T& operator()(int i0, int i1, int i2, int i3,
                 int i4, int i5);  
   const T& operator()(int i0, int i1, int i2, int i3,
                       int i4, int i5) const; 
   T& operator()(int i0, int i1, int i2, int i3,
                 int i4, int i5, int i6);  
   const T& operator()(int i0, int i1, int i2, int i3,
                       int i4, int i5, int i6) const; 
   T& operator()(int i0, int i1, int i2, int i3,
                 int i4, int i5, int i6, int i7);  
   const T& operator()(int i0, int i1, int i2, int i3,
                       int i4, int i5, int i6, int i7) const; 
   T& operator()(const std::vector<int>& ind);             
   const T& operator()(const std::vector<int>& ind) const; 
   T& operator()(const std::vector<unsigned int>& ind);
   const T& operator()(const std::vector<unsigned int>& ind) const;

   uint size() const {return m_store.size();}
   uint size(int i) const;     
   uint numDimensions() const {return m_sizes.size();}
   const std::vector<uint>& sizes() const {return m_sizes;}

   Array<T>& resize();
   Array<T>& resize(int size1);
   Array<T>& resize(int size1, int size2);
   Array<T>& resize(int size1, int size2, int size3);
   Array<T>& resize(int size1, int size2, int size3, int size4);
   Array<T>& resize(int size1, int size2, int size3, int size4,
                    int size5);
   Array<T>& resize(int size1, int size2, int size3, int size4,
                    int size5, int size6);
   Array<T>& resize(int size1, int size2, int size3, int size4,
                    int size5, int size6, int size7);
   Array<T>& resize(int size1, int size2, int size3, int size4,
                    int size5, int size6, int size7, int size8);
   Array<T>& resize(const std::vector<int>& sizevals);
   Array<T>& resize(const std::vector<unsigned int>& sizevals);



 private:

   unsigned int prod(const std::vector<uint>& ivec) const;
   int prod(const std::vector<int>& ivec) const;

   friend class IOHandler;

};



template <typename T>
Array<T>::Array()
{}

template <typename T>
Array<T>::Array(int size1) 
         : m_sizes(1), m_store(abs(size1))
{
 if (size1<1) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
}

template <typename T>
Array<T>::Array(int size1, int size2) 
         : m_sizes(2), m_store(abs(size1*size2))
{
 if ((size1<1)||(size2<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
}

template <typename T>
Array<T>::Array(int size1, int size2, int size3) 
         : m_sizes(3), m_store(abs(size1*size2*size3))
{
 if ((size1<1)||(size2<1)||(size3<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
}

template <typename T>
Array<T>::Array(int size1, int size2, int size3, int size4) 
         : m_sizes(4), m_store(abs(size1*size2*size3*size4))
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
}

template <typename T>
Array<T>::Array(int size1, int size2, int size3, int size4,
                int size5) 
         : m_sizes(5), m_store(abs(size1*size2*size3*size4*size5))
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
}

template <typename T>
Array<T>::Array(int size1, int size2, int size3, int size4,
                int size5, int size6) 
         : m_sizes(6), m_store(abs(size1*size2*size3*size4*size5*size6))
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)||(size6<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_sizes[5]=size6;
}

template <typename T>
Array<T>::Array(int size1, int size2, int size3, int size4,
                int size5, int size6, int size7) 
         : m_sizes(7), m_store(abs(size1*size2*size3*size4*size5*size6*size7))
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)||(size6<1)||(size7<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_sizes[5]=size6;
 m_sizes[6]=size7;
}

template <typename T>
Array<T>::Array(int size1, int size2, int size3, int size4,
                int size5, int size6, int size7, int size8) 
         : m_sizes(8), m_store(abs(size1*size2*size3*size4*size5*size6*size7*size8))
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)||(size6<1)||(size7<1)||(size8<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_sizes[5]=size6;
 m_sizes[6]=size7;
 m_sizes[7]=size8;
}

template <typename T>
Array<T>::Array(const std::vector<int>& sizevals) 
         : m_sizes(sizevals.begin(),sizevals.end()), m_store(prod(sizevals))
{
 for (std::vector<int>::const_iterator
    it=sizevals.begin();it!=sizevals.end();++it){
    if (*it<1) throw(std::invalid_argument("Invalid Array size"));}
}

template <typename T>
Array<T>::Array(const std::vector<unsigned int>& sizevals) 
         : m_sizes(sizevals), m_store(prod(sizevals))
{
 for (std::vector<unsigned int>::const_iterator
    it=sizevals.begin();it!=sizevals.end();++it){
    if (*it<1) throw(std::invalid_argument("Invalid Array size"));}
}

template <typename T>
Array<T>::Array(const std::vector<int>& sizevals, const T& initial_value) 
         : m_sizes(sizevals.begin(),sizevals.end()), m_store(prod(sizevals),initial_value)
{
 for (std::vector<int>::const_iterator
    it=sizevals.begin();it!=sizevals.end();++it){
    if (*it<1) throw(std::invalid_argument("Invalid Array size"));}
}

template <typename T>
Array<T>::Array(const std::vector<unsigned int>& sizevals, const T& initial_value) 
         : m_sizes(sizevals), m_store(prod(sizevals),initial_value)
{
 for (std::vector<unsigned int>::const_iterator
    it=sizevals.begin();it!=sizevals.end();++it){
    if (*it<1) throw(std::invalid_argument("Invalid Array size"));}
}


template <typename T>
Array<T>::Array(const Array<T>& incoming) 
         : m_sizes(incoming.m_sizes), m_store(incoming.m_store)
{}



template <typename T>
Array<T>::~Array()
{}


template <typename T>
Array<T>& Array<T>::clear()
{
 m_store.clear();
 m_sizes.clear();
 return *this;
}



template <typename T>
inline Array<T>& Array<T>::operator=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]=val;
 return *this;
}


template <typename T>
inline Array<T>& Array<T>::operator+=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]+=val;
 return *this;
}

template <typename T>
inline Array<T>& Array<T>::operator-=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]-=val;
 return *this;
}

template <typename T>
inline Array<T>& Array<T>::operator*=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]*=val;
 return *this;
}

template <typename T>
inline Array<T>& Array<T>::operator/=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]/=val;
 return *this;
}


template <typename T>
inline Array<T>& Array<T>::operator=(const Array<T>& incoming)
{
 if (this==&incoming) return *this;
 m_store=incoming.m_store;
 m_sizes=incoming.m_sizes;
 return *this;
}



template <typename T>
inline Array<T>& Array<T>::operator+=(const Array<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_sizes!=incoming.m_sizes) 
    throw(std::invalid_argument("Array size/shape mismatch"));
#endif
 for (uint i=0;i<m_store.size();i++) m_store[i]+=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Array<T>& Array<T>::operator-=(const Array<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_sizes!=incoming.m_sizes) 
    throw(std::invalid_argument("Array size/shape mismatch"));
#endif
 for (uint i=0;i<m_store.size();i++) m_store[i]-=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Array<T>& Array<T>::operator*=(const Array<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_sizes!=incoming.m_sizes) 
    throw(std::invalid_argument("Array size/shape mismatch"));
#endif
 for (uint i=0;i<m_store.size();i++) m_store[i]*=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Array<T>& Array<T>::operator/=(const Array<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_sizes!=incoming.m_sizes) 
    throw(std::invalid_argument("Array size/shape mismatch"));
#endif
 for (uint i=0;i<m_store.size();i++) m_store[i]/=incoming.m_store[i];
 return *this;
}


template <typename T>
inline T& Array<T>::operator[](int i)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=1) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i<0)||(i>=int(m_sizes[0]))) throw(std::invalid_argument("Array: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline const T& Array<T>::operator[](int i) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=1) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i<0)||(i>=int(m_sizes[0]))) throw(std::invalid_argument("Array: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline T& Array<T>::operator()(int i)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=1) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i<0)||(i>=int(m_sizes[0]))) throw(std::invalid_argument("Array: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline const T& Array<T>::operator()(int i) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=1) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i<0)||(i>=int(m_sizes[0]))) throw(std::invalid_argument("Array: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline T& Array<T>::operator()(int i0, int i1)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=2) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
#endif
 return m_store[i0+m_sizes[0]*i1];
}

template <typename T>
inline const T& Array<T>::operator()(int i0, int i1) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=2) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
#endif
 return m_store[i0+m_sizes[0]*i1];
}


template <typename T>
inline T& Array<T>::operator()(int i0, int i1, int i2)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=3) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*i2)];
}


template <typename T>
inline const T& Array<T>::operator()(int i0, int i1, int i2) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=3) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*i2)];
}

template <typename T>
inline T& Array<T>::operator()(int i0, int i1, int i2, int i3)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=4) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*i3))];
}


template <typename T>
inline const T& Array<T>::operator()(int i0, int i1, int i2, int i3) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=4) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*i3))];
}


template <typename T>
inline T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                               int i4)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=5) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]*i4)))];
}


template <typename T>
inline const T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                                     int i4) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=5) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]*i4)))];
}

template <typename T>
inline T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                               int i4, int i5)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=6) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
 if ((i5<0)||(i5>=int(m_sizes[5]))) throw(std::invalid_argument("Array: 6th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]
                *(i4+m_sizes[4]*i5))))];
}


template <typename T>
inline const T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                                     int i4, int i5) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=6) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
 if ((i5<0)||(i5>=int(m_sizes[5]))) throw(std::invalid_argument("Array: 6th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]
                *(i4+m_sizes[4]*i5))))];
}

template <typename T>
inline T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                               int i4, int i5, int i6)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=7) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
 if ((i5<0)||(i5>=int(m_sizes[5]))) throw(std::invalid_argument("Array: 6th index out of range"));
 if ((i6<0)||(i6>=int(m_sizes[6]))) throw(std::invalid_argument("Array: 7th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]*(i4+m_sizes[4]
                *(i5+m_sizes[5]*i6)))))];
}


template <typename T>
inline const T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                                     int i4, int i5, int i6) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=7) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
 if ((i5<0)||(i5>=int(m_sizes[5]))) throw(std::invalid_argument("Array: 6th index out of range"));
 if ((i6<0)||(i6>=int(m_sizes[6]))) throw(std::invalid_argument("Array: 7th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]*(i4+m_sizes[4]
                *(i5+m_sizes[5]*i6)))))];
}

template <typename T>
inline T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                               int i4, int i5, int i6, int i7)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=8) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
 if ((i5<0)||(i5>=int(m_sizes[5]))) throw(std::invalid_argument("Array: 6th index out of range"));
 if ((i6<0)||(i6>=int(m_sizes[6]))) throw(std::invalid_argument("Array: 7th index out of range"));
 if ((i7<0)||(i7>=int(m_sizes[7]))) throw(std::invalid_argument("Array: 8th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]
                *(i4+m_sizes[4]*(i5+m_sizes[5]*(i6+m_sizes[6]*i7))))))];
}


template <typename T>
inline const T& Array<T>::operator()(int i0, int i1, int i2, int i3,
                                     int i4, int i5, int i6, int i7) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=8) throw(std::invalid_argument("Array: improper rank of array indices"));
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Array: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Array: 2nd index out of range"));
 if ((i2<0)||(i2>=int(m_sizes[2]))) throw(std::invalid_argument("Array: 3rd index out of range"));
 if ((i3<0)||(i3>=int(m_sizes[3]))) throw(std::invalid_argument("Array: 4th index out of range"));
 if ((i4<0)||(i4>=int(m_sizes[4]))) throw(std::invalid_argument("Array: 5th index out of range"));
 if ((i5<0)||(i5>=int(m_sizes[5]))) throw(std::invalid_argument("Array: 6th index out of range"));
 if ((i6<0)||(i6>=int(m_sizes[6]))) throw(std::invalid_argument("Array: 7th index out of range"));
 if ((i7<0)||(i7>=int(m_sizes[7]))) throw(std::invalid_argument("Array: 8th index out of range"));
#endif
 return m_store[i0+m_sizes[0]*(i1+m_sizes[1]*(i2+m_sizes[2]*(i3+m_sizes[3]
                *(i4+m_sizes[4]*(i5+m_sizes[5]*(i6+m_sizes[6]*i7))))))];
}

template <typename T>
inline T& Array<T>::operator()(const std::vector<int>& ind)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=ind.size()) throw(std::invalid_argument("Array: improper rank of array indices"));
 for (uint i=0; i<m_sizes.size(); ++i)
    if ((ind[i]<0)||(ind[i]>=int(m_sizes[i])))
       throw(std::invalid_argument("Array: an index is out of range"));
#endif
 if (m_sizes.empty()) throw(std::invalid_argument("Empty Array...element access not possible"));
 uint off=ind[m_sizes.size()-1];
 for (int i=int(m_sizes.size())-2; i>=0; --i)
    off=off*m_sizes[i]+ind[i];
 return m_store[off];
}

template <typename T>
inline const T& Array<T>::operator()(const std::vector<int>& ind) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=ind.size()) throw(std::invalid_argument("Array: improper rank of array indices"));
 for (uint i=0; i<m_sizes.size(); ++i)
    if ((ind[i]<0)||(ind[i]>=int(m_sizes[i])))
       throw(std::invalid_argument("Array: an index is out of range"));
#endif
 if (m_sizes.empty()) throw(std::invalid_argument("Empty Array...element access not possible"));
 uint off=ind[m_sizes.size()-1];
 for (int i=int(m_sizes.size())-2; i>=0; --i)
    off=off*m_sizes[i]+ind[i];
 return m_store[off];
}


template <typename T>
inline T& Array<T>::operator()(const std::vector<uint>& ind)
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=ind.size()) throw(std::invalid_argument("Array: improper rank of array indices"));
 for (uint i=0; i<m_sizes.size(); ++i)
    if ((ind[i]<0)||(ind[i]>=m_sizes[i]))
       throw(std::invalid_argument("Array: an index is out of range"));
#endif
 if (m_sizes.empty()) throw(std::invalid_argument("Empty Array...element access not possible"));
 uint off=ind[m_sizes.size()-1];
 for (int i=int(m_sizes.size())-2; i>=0; --i)
    off=off*m_sizes[i]+ind[i];
 return m_store[off];
}


template <typename T>
inline const T& Array<T>::operator()(const std::vector<uint>& ind) const
{
#ifdef SAFETY_FLAG
 if (m_sizes.size()!=ind.size()) throw(std::invalid_argument("Array: improper rank of array indices"));
 for (uint i=0; i<m_sizes.size(); ++i)
    if ((ind[i]<0)||(ind[i]>=m_sizes[i]))
       throw(std::invalid_argument("Array: an index is out of range"));
#endif
 if (m_sizes.empty()) throw(std::invalid_argument("Empty Array...element access not possible"));
 uint off=ind[m_sizes.size()-1];
 for (int i=int(m_sizes.size())-2; i>=0; --i)
    off=off*m_sizes[i]+ind[i];
 return m_store[off];
}


template <typename T>
inline uint Array<T>::size(int i) const
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_sizes.size()))) throw(std::invalid_argument("Array size index out of range"));
#endif
 return m_sizes[i];
}


template <typename T>
Array<T>& Array<T>::resize()
{
 return clear();
}

template <typename T>
Array<T>& Array<T>::resize(int size1)
{
 if (size1<1) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(1);
 m_sizes[0]=size1;
 m_store.resize(size1);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(int size1, int size2)
{
 if ((size1<1)||(size2<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(2);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_store.resize(size1*size2);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(int size1, int size2, int size3)
{
 if ((size1<1)||(size2<1)||(size3<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(3);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_store.resize(size1*size2*size3);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(int size1, int size2, int size3, int size4)
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(4);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_store.resize(size1*size2*size3*size4);
 return *this;
}


template <typename T>
Array<T>& Array<T>::resize(int size1, int size2, int size3, int size4,
                           int size5) 
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(5);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_store.resize(size1*size2*size3*size4*size5);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(int size1, int size2, int size3, int size4,
                           int size5, int size6) 
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)||(size6<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(6);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_sizes[5]=size6;
 m_store.resize(size1*size2*size3*size4*size5*size6);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(int size1, int size2, int size3, int size4,
                           int size5, int size6, int size7) 
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)||(size6<1)||(size7<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(7);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_sizes[5]=size6;
 m_sizes[6]=size7;
 m_store.resize(size1*size2*size3*size4*size5*size6*size7);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(int size1, int size2, int size3, int size4,
                           int size5, int size6, int size7, int size8) 
{
 if ((size1<1)||(size2<1)||(size3<1)||(size4<1)
    ||(size5<1)||(size6<1)||(size7<1)||(size8<1)) 
    throw(std::invalid_argument("Invalid Array size"));
 m_sizes.resize(8);
 m_sizes[0]=size1;
 m_sizes[1]=size2;
 m_sizes[2]=size3;
 m_sizes[3]=size4;
 m_sizes[4]=size5;
 m_sizes[5]=size6;
 m_sizes[6]=size7;
 m_sizes[7]=size8;
 m_store.resize(size1*size2*size3*size4*size5*size6*size7*size8);
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(const std::vector<int>& sizevals)
{
 for (std::vector<int>::const_iterator
    it=sizevals.begin();it!=sizevals.end();++it){
    if (*it<1) throw(std::invalid_argument("Invalid Array size"));}
 m_sizes.assign(sizevals.begin(),sizevals.end());
 m_store.resize(prod(sizevals));
 return *this;
}

template <typename T>
Array<T>& Array<T>::resize(const std::vector<unsigned int>& sizevals)
{
 for (std::vector<unsigned int>::const_iterator
    it=sizevals.begin();it!=sizevals.end();++it){
    if (*it<1) throw(std::invalid_argument("Invalid Array size"));}
 m_sizes=sizevals;
 m_store.resize(prod(sizevals));
 return *this;
}


template <typename T>
unsigned int Array<T>::prod(const std::vector<uint>& ivec) const
{
 unsigned int res=1;
 for (std::vector<unsigned int>::const_iterator
    it=ivec.begin();it!=ivec.end();++it) res*=(*it);
 if (res<0) res=-res;
 return res;
}

template <typename T>
int Array<T>::prod(const std::vector<int>& ivec) const
{
 int res=1;
 for (std::vector<int>::const_iterator
    it=ivec.begin();it!=ivec.end();++it) res*=(*it);
 if (res<0) res=-res;
 return res;
}


// **************************************************************
}
#endif
