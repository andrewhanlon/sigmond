#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <complex>
#include "scalar_defs.h"
#include <stdexcept>

// ***************************************************************
// *                                                             *
// *   The classes "Vector" and "Matrix" are defined in          *
// *   this file.  These are used by the data analysis           *
// *   routines.  "Vector" and "Matrix" are 1-dimensional and    *
// *   2-dimensional, respectively, versions of "Array".         *
// *   The class "LowerTriangularMatrix" is also defined,        *
// *   which is useful for Cholesky decompositions and only      *
// *   stores the nonzero elements in the lower triangle.        *
// *   These classes are templated.  Individual elements         *
// *   can be assigned and accessed with parentheses and         *
// *   square brackets, in the case of "Vector".  However,       *
// *   "put" and "get" methods are also available.               *
// *   "LowerTriangularMatrix" assumes a square matrix,          *
// *   and throws an exception when attempting to access an      *
// *   element in the upper right triangle.                      *
// *                                                             *
// *   All indices are ZERO offset.                              *
// *                                                             *
// *   Two other important classes defined here are              *
// *   "RealSymmetricMatrix" and "ComplexHermitianMatrix".       *
// *   Only the upper right triangle is stored to save memory.   *
// *   Assignment of elements for "ComplexHermitianMatrix"       *
// *   MUST be done using the "put" member, so that assignment   *
// *   to diagonal elements can be checked for a zero imaginary  *
// *   part.                                                     *
// *                                                             *
// *   The typedefs "RVector", "CVector", "RMatrix", "CMatrix",  *
// *   and "LTMatrix" are also defined in this file.             *
// *                                                             *
// *                                                             *
// ***************************************************************


  // Use the #define below to include range checking code.
  // Do not define and range checking code is omitted (faster).
#define SAFETY_FLAG

// **************************************************************



template <typename T>
class Vector {

   std::vector<T> m_store;

 public:

   Vector();
   Vector(int size);
   Vector(uint size);
   Vector(int size, const T& initial_value);
   Vector(uint size, const T& initial_value);
   Vector(const Vector<T>& incoming);
   Vector(const std::vector<T>& incoming);

   ~Vector();
   Vector<T>& clear(); 

   Vector<T>& operator=(const T& val);
   Vector<T>& operator+=(const T& val);
   Vector<T>& operator-=(const T& val);
   Vector<T>& operator*=(const T& val);
   Vector<T>& operator/=(const T& val);
   Vector<T>& operator=(const Vector<T>& incoming); 
   Vector<T>& operator+=(const Vector<T>& incoming);
   Vector<T>& operator-=(const Vector<T>& incoming);
   Vector<T>& operator*=(const Vector<T>& incoming);
   Vector<T>& operator/=(const Vector<T>& incoming);

   T& operator()(int i);             
   T& operator[](int i);             
   const T& operator[](int i) const; 
   const T& operator()(int i) const; 
   void put(int i, const T& inputvalue);             
   const T& get(int i) const; 

   T& operator()(uint i);             
   T& operator[](uint i);             
   const T& operator[](uint i) const; 
   const T& operator()(uint i) const; 
   void put(uint i, const T& inputvalue);             
   const T& get(uint i) const; 

   uint size() const {return m_store.size();}

   Vector<T>& resize();
   Vector<T>& resize(int size);
   Vector<T>& resize(uint size);
   std::vector<T> c_vector() const;

   std::vector<T>& c_vector_ref();  // use with caution

};



template <typename T>
Vector<T>::Vector()
{}


template <typename T>
Vector<T>::Vector(int size)  : m_store((size>0)?size:0)
{
 if (size<1) throw(std::invalid_argument("Invalid Vector size"));
}


template <typename T>
Vector<T>::Vector(uint size)  : m_store(size)
{}


template <typename T>
Vector<T>::Vector(int size, const T& initial_value) 
         : m_store((size>0)?size:0,initial_value)
{
 if (size<1) throw(std::invalid_argument("Invalid Vector size"));
}


template <typename T>
Vector<T>::Vector(uint size, const T& initial_value) 
         : m_store(size,initial_value)
{}


template <typename T>
Vector<T>::Vector(const Vector<T>& incoming) 
         : m_store(incoming.m_store)
{}


template <typename T>
Vector<T>::Vector(const std::vector<T>& incoming) 
         : m_store(incoming)
{}


template <typename T>
Vector<T>::~Vector()
{}


template <typename T>
Vector<T>& Vector<T>::clear()
{
 m_store.clear();
 return *this;
}



template <typename T>
inline Vector<T>& Vector<T>::operator=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]=val;
 return *this;
}


template <typename T>
inline Vector<T>& Vector<T>::operator+=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]+=val;
 return *this;
}

template <typename T>
inline Vector<T>& Vector<T>::operator-=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]-=val;
 return *this;
}

template <typename T>
inline Vector<T>& Vector<T>::operator*=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]*=val;
 return *this;
}

template <typename T>
inline Vector<T>& Vector<T>::operator/=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]/=val;
 return *this;
}


template <typename T>
inline Vector<T>& Vector<T>::operator=(const Vector<T>& incoming)
{
 if (this==&incoming) return *this;
 m_store=incoming.m_store;
 return *this;
}



template <typename T>
inline Vector<T>& Vector<T>::operator+=(const Vector<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_store.size()!=incoming.m_store.size())
    throw(std::invalid_argument("Vector size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]+=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Vector<T>& Vector<T>::operator-=(const Vector<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_store.size()!=incoming.m_store.size())
    throw(std::invalid_argument("Vector size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]-=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Vector<T>& Vector<T>::operator*=(const Vector<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_store.size()!=incoming.m_store.size())
    throw(std::invalid_argument("Vector size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]*=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Vector<T>& Vector<T>::operator/=(const Vector<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_store.size()!=incoming.m_store.size())
    throw(std::invalid_argument("Vector size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]/=incoming.m_store[i];
 return *this;
}


template <typename T>
inline T& Vector<T>::operator[](int i)
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_store.size()))) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline const T& Vector<T>::operator[](int i) const
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_store.size()))) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline T& Vector<T>::operator()(int i)
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_store.size()))) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline const T& Vector<T>::operator()(int i) const
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_store.size()))) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline void Vector<T>::put(int i, const T& inputvalue)
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_store.size()))) throw(std::invalid_argument("Vector: index out of range"));
#endif
 m_store[i]=inputvalue;
}

template <typename T>
inline const T& Vector<T>::get(int i) const
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>=int(m_store.size()))) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline T& Vector<T>::operator[](uint i)
{
#ifdef SAFETY_FLAG
 if (i>=m_store.size()) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline const T& Vector<T>::operator[](uint i) const
{
#ifdef SAFETY_FLAG
 if (i>=m_store.size()) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline T& Vector<T>::operator()(uint i)
{
#ifdef SAFETY_FLAG
 if (i>=m_store.size()) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline const T& Vector<T>::operator()(uint i) const
{
#ifdef SAFETY_FLAG
 if (i>=m_store.size()) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
inline void Vector<T>::put(uint i, const T& inputvalue)
{
#ifdef SAFETY_FLAG
 if (i>=m_store.size()) throw(std::invalid_argument("Vector: index out of range"));
#endif
 m_store[i]=inputvalue;
}

template <typename T>
inline const T& Vector<T>::get(uint i) const
{
#ifdef SAFETY_FLAG
 if (i>=m_store.size()) throw(std::invalid_argument("Vector: index out of range"));
#endif
 return m_store[i];
}


template <typename T>
Vector<T>& Vector<T>::resize()
{
 return clear();
}


template <typename T>
Vector<T>& Vector<T>::resize(int size)
{
 if (size<1) clear();
 m_store.resize(size);
 return *this;
}


template <typename T>
Vector<T>& Vector<T>::resize(uint size)
{
 m_store.resize(size);
 return *this;
}


template <typename T>
std::vector<T> Vector<T>::c_vector() const
{
 return m_store;
}



template <typename T>
std::vector<T>& Vector<T>::c_vector_ref()
{
 return m_store;
}


// **************************************************************



template <typename T>
class Matrix {

   std::vector<unsigned int> m_sizes;
   std::vector<T> m_store;

 public:

   Matrix();
   Matrix(int size0, int size1);
   Matrix(uint size0, uint size1);
   Matrix(int size0, int size1, const T& initial_value);
   Matrix(uint size0, uint size1, const T& initial_value);
   Matrix(const std::vector<unsigned int>& sizevals);
   Matrix(const std::vector<int>& sizevals);
   Matrix(const std::vector<unsigned int>& sizevals, const T& initial_value);
   Matrix(const std::vector<int>& sizevals, const T& initial_value);
   Matrix(const Matrix<T>& incoming);

   ~Matrix();
   Matrix<T>& clear(); 

   Matrix<T>& operator=(const T& val);
   Matrix<T>& operator+=(const T& val);
   Matrix<T>& operator-=(const T& val);
   Matrix<T>& operator*=(const T& val);
   Matrix<T>& operator/=(const T& val);
   Matrix<T>& operator=(const Matrix<T>& incoming); 
   Matrix<T>& operator+=(const Matrix<T>& incoming);
   Matrix<T>& operator-=(const Matrix<T>& incoming);

   T& operator()(int i0, int i1);    
   T& operator()(uint i0, uint i1);    
   T& operator()(int i0, uint i1) {return operator()(i0,int(i1));}
   T& operator()(uint i0, int i1) {return operator()(int(i0),i1);}
   T& operator()(const std::vector<int>& ind);             
   T& operator()(const std::vector<unsigned int>& ind);
   const T& operator()(uint i0, uint i1) const; 
   const T& operator()(int i0, int i1) const; 
   const T& operator()(const std::vector<int>& ind) const; 
   const T& operator()(const std::vector<unsigned int>& ind) const;

   void put(int i0, int i1, const T& inputvalue);    
   void put(uint i0, uint i1, const T& inputvalue);    
   void put(int i0, uint i1, const T& inputvalue) {put(i0,int(i1),inputvalue);}   
   void put(uint i0, int i1, const T& inputvalue) {put(int(i0),i1,inputvalue);}    
   void put(const std::vector<int>& ind, const T& inputvalue);             
   void put(const std::vector<unsigned int>& ind, const T& inputvalue);
   const T& get(int i0, int i1) const; 
   const T& get(uint i0, uint i1) const; 
   const T& get(int i0, uint i1) const {return get(i0,int(i1));}
   const T& get(uint i0, int i1) const {return get(int(i0),i1);} 
   const T& get(const std::vector<int>& ind) const; 
   const T& get(const std::vector<unsigned int>& ind) const;

   uint size() const {return m_store.size();}
   uint size(int i) const;     
   uint size(uint i) const;     
   const std::vector<uint>& sizes() const {return m_sizes;}

   Matrix<T>& resize();
   Matrix<T>& resize(int size1, int size2);
   Matrix<T>& resize(uint size1, uint size2);
   Matrix<T>& resize(const std::vector<int>& sizevals);
   Matrix<T>& resize(const std::vector<unsigned int>& sizevals);


 private:

   unsigned int get_index(int ind0, int ind1) const;
   unsigned int get_index(uint ind0, uint ind1) const;
   unsigned int get_index(const std::vector<int>& ind) const;
   unsigned int get_index(const std::vector<uint>& ind) const;

   unsigned int set_sizes(int size0, int size1);
   unsigned int set_sizes(uint size0, uint size1);
   unsigned int set_sizes(const std::vector<uint>& ivec);
   unsigned int set_sizes(const std::vector<int>& ivec);

};



template <typename T>
unsigned int Matrix<T>::set_sizes(int size0, int size1)
{
 if ((size0<0)||(size1<0)){
    m_sizes[0]=m_sizes[1]=0;
    throw(std::invalid_argument("Invalid Matrix size"));}
 if ((size0==0)||(size1==0)){
    m_sizes[0]=m_sizes[1]=0;
    return 0;}
 m_sizes[0]=size0;
 m_sizes[1]=size1;
 return size0*size1;
}

template <typename T>
unsigned int Matrix<T>::set_sizes(uint size0, uint size1)
{
 if ((size0==0)||(size1==0)){
    m_sizes[0]=m_sizes[1]=0;
    return 0;}
 m_sizes[0]=size0;
 m_sizes[1]=size1;
 return size0*size1;
}

template <typename T>
unsigned int Matrix<T>::set_sizes(const std::vector<uint>& sizevals)
{
 if (sizevals.size()!=2) throw(std::invalid_argument("Invalid Matrix sizes"));
 return set_sizes(sizevals[0],sizevals[1]);
}

template <typename T>
unsigned int Matrix<T>::set_sizes(const std::vector<int>& sizevals)
{
 if (sizevals.size()!=2) throw(std::invalid_argument("Invalid Matrix sizes"));
 return set_sizes(sizevals[0],sizevals[1]);
}



template <typename T>
Matrix<T>::Matrix()  : m_sizes(2)
{
 m_sizes[0]=m_sizes[1]=0;
}

template <typename T>
Matrix<T>::Matrix(int size0, int size1) 
         : m_sizes(2), m_store(set_sizes(size0,size1)) 
{}

template <typename T>
Matrix<T>::Matrix(uint size0, uint size1) 
         : m_sizes(2), m_store(set_sizes(size0,size1)) 
{}

template <typename T>
Matrix<T>::Matrix(int size0, int size1, const T& initial_value) 
         : m_sizes(2), m_store(set_sizes(size0,size1),initial_value)
{}

template <typename T>
Matrix<T>::Matrix(uint size0, uint size1, const T& initial_value) 
         : m_sizes(2), m_store(set_sizes(size0,size1),initial_value)
{}

template <typename T>
Matrix<T>::Matrix(const std::vector<int>& sizevals) 
         : m_sizes(2), m_store(set_sizes(sizevals))
{}

template <typename T>
Matrix<T>::Matrix(const std::vector<unsigned int>& sizevals) 
         : m_sizes(2), m_store(set_sizes(sizevals))
{}

template <typename T>
Matrix<T>::Matrix(const std::vector<int>& sizevals, const T& initial_value) 
         : m_sizes(2), m_store(set_sizes(sizevals),initial_value)
{}

template <typename T>
Matrix<T>::Matrix(const std::vector<unsigned int>& sizevals, const T& initial_value) 
         : m_sizes(2), m_store(set_sizes(sizevals),initial_value)
{}


template <typename T>
Matrix<T>::Matrix(const Matrix<T>& incoming) 
         : m_sizes(incoming.m_sizes), m_store(incoming.m_store)
{}



template <typename T>
Matrix<T>::~Matrix()
{}


template <typename T>
Matrix<T>& Matrix<T>::clear()
{
 m_store.clear();
 m_sizes[0]=m_sizes[1]=0;
 return *this;
}



template <typename T>
inline Matrix<T>& Matrix<T>::operator=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]=val;
 return *this;
}


template <typename T>
inline Matrix<T>& Matrix<T>::operator+=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]+=val;
 return *this;
}

template <typename T>
inline Matrix<T>& Matrix<T>::operator-=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]-=val;
 return *this;
}

template <typename T>
inline Matrix<T>& Matrix<T>::operator*=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]*=val;
 return *this;
}

template <typename T>
inline Matrix<T>& Matrix<T>::operator/=(const T& val)
{
 for (uint i=0;i<m_store.size();i++) m_store[i]/=val;
 return *this;
}


template <typename T>
inline Matrix<T>& Matrix<T>::operator=(const Matrix<T>& incoming)
{
 if (this==&incoming) return *this;
 m_store=incoming.m_store;
 m_sizes=incoming.m_sizes;
 return *this;
}



template <typename T>
inline Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_sizes!=incoming.m_sizes) 
    throw(std::invalid_argument("Matrix size/shape mismatch"));
#endif
 for (uint i=0;i<m_store.size();i++) m_store[i]+=incoming.m_store[i];
 return *this;
}

template <typename T>
inline Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_sizes!=incoming.m_sizes) 
    throw(std::invalid_argument("Matrix size/shape mismatch"));
#endif
 for (uint i=0;i<m_store.size();i++) m_store[i]-=incoming.m_store[i];
 return *this;
}



template <typename T>
inline unsigned int Matrix<T>::get_index(int i0, int i1) const
{
#ifdef SAFETY_FLAG
 if ((i0<0)||(i0>=int(m_sizes[0]))) throw(std::invalid_argument("Matrix: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_sizes[1]))) throw(std::invalid_argument("Matrix: 2nd index out of range"));
#endif
 return i0+m_sizes[0]*i1;
}


template <typename T>
inline unsigned int Matrix<T>::get_index(uint i0, uint i1) const
{
#ifdef SAFETY_FLAG
 if (i0>=m_sizes[0]) throw(std::invalid_argument("Matrix: 1st index out of range"));
 if (i1>=m_sizes[1]) throw(std::invalid_argument("Matrix: 2nd index out of range"));
#endif
 return i0+m_sizes[0]*i1;
}


template <typename T>
inline unsigned int Matrix<T>::get_index(const std::vector<int>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("Matrix: improper rank of Matrix indices"));
#endif
 return get_index(ind[0],ind[1]);
}


template <typename T>
inline unsigned int Matrix<T>::get_index(const std::vector<uint>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("Matrix: improper rank of Matrix indices"));
#endif
 return get_index(ind[0],ind[1]);
}


template <typename T>
inline T& Matrix<T>::operator()(int i0, int i1)
{
 return m_store[get_index(i0,i1)];
}


template <typename T>
inline T& Matrix<T>::operator()(uint i0, uint i1)
{
 return m_store[get_index(i0,i1)];
}


template <typename T>
inline T& Matrix<T>::operator()(const std::vector<int>& ind)
{
 return m_store[get_index(ind)];
}


template <typename T>
inline T& Matrix<T>::operator()(const std::vector<uint>& ind)
{
 return m_store[get_index(ind)];
}


template <typename T>
inline const T& Matrix<T>::operator()(int i0, int i1) const
{
 return m_store[get_index(i0,i1)];
}


template <typename T>
inline const T& Matrix<T>::operator()(uint i0, uint i1) const
{
 return m_store[get_index(i0,i1)];
}


template <typename T>
inline const T& Matrix<T>::operator()(const std::vector<int>& ind) const
{
 return m_store[get_index(ind)];
}


template <typename T>
inline const T& Matrix<T>::operator()(const std::vector<uint>& ind) const
{
 return m_store[get_index(ind)];
}




template <typename T>
inline void Matrix<T>::put(int i0, int i1, const T& inputvalue)
{
 m_store[get_index(i0,i1)]=inputvalue;
}


template <typename T>
inline void Matrix<T>::put(uint i0, uint i1, const T& inputvalue)
{
 m_store[get_index(i0,i1)]=inputvalue;
}


template <typename T>
inline void Matrix<T>::put(const std::vector<int>& ind, const T& inputvalue)
{
 m_store[get_index(ind)]=inputvalue;
}


template <typename T>
inline void Matrix<T>::put(const std::vector<uint>& ind, const T& inputvalue)
{
 m_store[get_index(ind)]=inputvalue;
}


template <typename T>
inline const T& Matrix<T>::get(int i0, int i1) const
{
 return m_store[get_index(i0,i1)];
}


template <typename T>
inline const T& Matrix<T>::get(uint i0, uint i1) const
{
 return m_store[get_index(i0,i1)];
}


template <typename T>
inline const T& Matrix<T>::get(const std::vector<int>& ind) const
{
 return m_store[get_index(ind)];
}


template <typename T>
inline const T& Matrix<T>::get(const std::vector<uint>& ind) const
{
 return m_store[get_index(ind)];
}




template <typename T>
inline uint Matrix<T>::size(int i) const
{
#ifdef SAFETY_FLAG
 if ((i<0)||(i>1)) throw(std::invalid_argument("Matrix size index out of range"));
#endif
 return m_sizes[i];
}


template <typename T>
inline uint Matrix<T>::size(uint i) const
{
#ifdef SAFETY_FLAG
 if (i>1) throw(std::invalid_argument("Matrix size index out of range"));
#endif
 return m_sizes[i];
}


template <typename T>
Matrix<T>& Matrix<T>::resize()
{
 return clear();
}

template <typename T>
Matrix<T>& Matrix<T>::resize(int size0, int size1)
{
 m_store.resize(set_sizes(size0,size1));
 return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::resize(uint size0, uint size1)
{
 m_store.resize(set_sizes(size0,size1));
 return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::resize(const std::vector<int>& sizevals)
{
 m_store.resize(set_sizes(sizevals));
 return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::resize(const std::vector<unsigned int>& sizevals)
{
 m_store.resize(set_sizes(sizevals));
 return *this;
}





// **************************************************************




template <typename T>
class LowerTriangularMatrix {

   std::vector<T> m_store;
   unsigned int m_size;

 public:

   LowerTriangularMatrix();
   LowerTriangularMatrix(int insize);
   LowerTriangularMatrix(uint insize);
   LowerTriangularMatrix(int insize, const T& initial_value);
   LowerTriangularMatrix(uint insize, const T& initial_value);
   LowerTriangularMatrix(const LowerTriangularMatrix<T>& incoming);

   ~LowerTriangularMatrix();
   LowerTriangularMatrix<T>& clear(); 

   LowerTriangularMatrix<T>& operator=(const T& val);
   LowerTriangularMatrix<T>& operator+=(const T& val);
   LowerTriangularMatrix<T>& operator-=(const T& val);
   LowerTriangularMatrix<T>& operator*=(const T& val);
   LowerTriangularMatrix<T>& operator/=(const T& val);
   LowerTriangularMatrix<T>& operator=(const LowerTriangularMatrix<T>& incoming); 
   LowerTriangularMatrix<T>& operator+=(const LowerTriangularMatrix<T>& incoming);
   LowerTriangularMatrix<T>& operator-=(const LowerTriangularMatrix<T>& incoming);

   T& operator()(int i0, int i1);    
   T& operator()(uint i0, uint i1);    
   T& operator()(int i0, uint i1) {return operator()(i0,int(i1));}   
   T& operator()(uint i0, int i1) {return operator()(int(i0),i1);}   
   T& operator()(const std::vector<int>& ind);             
   T& operator()(const std::vector<unsigned int>& ind);
   const T& operator()(int i0, int i1) const; 
   const T& operator()(uint i0, uint i1) const; 
   const T& operator()(int i0, uint i1) const {return operator()(i0,int(i1));}
   const T& operator()(uint i0, int i1) const {return operator()(int(i0),i1);}
   const T& operator()(const std::vector<int>& ind) const; 
   const T& operator()(const std::vector<unsigned int>& ind) const;

   void put(int i0, int i1, const T& inputvalue);    
   void put(uint i0, uint i1, const T& inputvalue);    
   void put(int i0, uint i1, const T& inputvalue) {put(i0,int(i1),inputvalue);}    
   void put(uint i0, int i1, const T& inputvalue) {put(int(i0),i1,inputvalue);} 
   void put(const std::vector<int>& ind, const T& inputvalue);             
   void put(const std::vector<unsigned int>& ind, const T& inputvalue);
   const T& get(int i0, int i1) const; 
   const T& get(uint i0, uint i1) const; 
   const T& get(int i0, uint i1) const {return get(i0,int(i1));}
   const T& get(uint i0, int i1) const {return get(int(i0),i1);}
   const T& get(const std::vector<int>& ind) const; 
   const T& get(const std::vector<unsigned int>& ind) const;

   uint size() const {return m_size;}

   LowerTriangularMatrix<T>& resize();
   LowerTriangularMatrix<T>& resize(int size);
   LowerTriangularMatrix<T>& resize(uint size);


 private:

   unsigned int get_index(int ind0, int ind1) const;
   unsigned int get_index(uint ind0, uint ind1) const;
   unsigned int get_index(const std::vector<int>& ind) const;
   unsigned int get_index(const std::vector<uint>& ind) const;

   unsigned int set_size(int insize);
   unsigned int set_size(uint insize);

};


template <typename T>
unsigned int LowerTriangularMatrix<T>::set_size(int insize)
{
 if (insize<0) throw(std::invalid_argument("Invalid LowerTriangularMatrix size"));
 m_size=insize;
 return (m_size*(m_size+1))/2;
}


template <typename T>
unsigned int LowerTriangularMatrix<T>::set_size(uint insize)
{
 m_size=insize;
 return (m_size*(m_size+1))/2;
}


template <typename T>
LowerTriangularMatrix<T>::LowerTriangularMatrix()
{m_size=0;}


template <typename T>
LowerTriangularMatrix<T>::LowerTriangularMatrix(int insize)  
                   : m_store(set_size(insize))
{}


template <typename T>
LowerTriangularMatrix<T>::LowerTriangularMatrix(uint insize)  
                   : m_store(set_size(insize))
{}


template <typename T>
LowerTriangularMatrix<T>::LowerTriangularMatrix(int insize, const T& initial_value) 
                   : m_store(set_size(insize),initial_value)
{}


template <typename T>
LowerTriangularMatrix<T>::LowerTriangularMatrix(uint insize, const T& initial_value) 
                   : m_store(set_size(insize),initial_value)
{}


template <typename T>
LowerTriangularMatrix<T>::LowerTriangularMatrix(const LowerTriangularMatrix<T>& incoming) 
         : m_store(incoming.m_store), m_size(incoming.m_size)
{}


template <typename T>
LowerTriangularMatrix<T>::~LowerTriangularMatrix()
{}


template <typename T>
LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::clear()
{
 m_store.clear();
 m_size=0;
 return *this;
}



template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]=val;
 return *this;
}


template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator+=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]+=val;
 return *this;
}

template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator-=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]-=val;
 return *this;
}

template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator*=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]*=val;
 return *this;
}

template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator/=(const T& val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]/=val;
 return *this;
}


template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator=(const LowerTriangularMatrix<T>& incoming)
{
 if (this==&incoming) return *this;
 m_store=incoming.m_store;
 m_size=incoming.m_size;
 return *this;
}



template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator+=(const LowerTriangularMatrix<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_size!=incoming.m_size)
    throw(std::invalid_argument("LowerTriangularMatrix size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]+=incoming.m_store[i];
 return *this;
}

template <typename T>
inline LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::operator-=(const LowerTriangularMatrix<T>& incoming)
{
#ifdef SAFETY_FLAG
 if (m_size!=incoming.m_size)
    throw(std::invalid_argument("LowerTriangularMatrix size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]-=incoming.m_store[i];
 return *this;
}




template <typename T>
inline unsigned int LowerTriangularMatrix<T>::get_index(int i0, int i1) const
{
#ifdef SAFETY_FLAG
 if ((i0<0)||(i0>=int(m_size))) throw(std::invalid_argument("LowerTriangularMatrix: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_size))) throw(std::invalid_argument("LowerTriangularMatrix: 2nd index out of range"));
 if (i0<i1) throw(std::invalid_argument("LowerTriangularMatrix: 2nd index cannot exceed 1st index"));
#endif
 return (i0*(i0+1))/2+i1;
}


template <typename T>
inline unsigned int LowerTriangularMatrix<T>::get_index(uint i0, uint i1) const
{
#ifdef SAFETY_FLAG
 if (i0>=m_size) throw(std::invalid_argument("LowerTriangularMatrix: 1st index out of range"));
 if (i1>=m_size) throw(std::invalid_argument("LowerTriangularMatrix: 2nd index out of range"));
 if (i0<i1) throw(std::invalid_argument("LowerTriangularMatrix: 2nd index cannot exceed 1st index"));
#endif
 return (i0*(i0+1))/2+i1;
}


template <typename T>
inline unsigned int LowerTriangularMatrix<T>::get_index(const std::vector<int>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("LowerTriangularMatrix: improper rank of LowerTriangularMatrix indices"));
#endif
 return get_index(ind[0],ind[1]);
}


template <typename T>
inline unsigned int LowerTriangularMatrix<T>::get_index(const std::vector<uint>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("LowerTriangularMatrix: improper rank of LowerTriangularMatrix indices"));
#endif
 return get_index(ind[0],ind[1]);
}



template <typename T>
inline T& LowerTriangularMatrix<T>::operator()(int i0, int i1)
{
 return m_store[get_index(i0,i1)];
}

template <typename T>
inline T& LowerTriangularMatrix<T>::operator()(uint i0, uint i1)
{
 return m_store[get_index(i0,i1)];
}

template <typename T>
inline T& LowerTriangularMatrix<T>::operator()(const std::vector<int>& ind)
{
 return m_store[get_index(ind)];
}

template <typename T>
inline T& LowerTriangularMatrix<T>::operator()(const std::vector<uint>& ind)
{
 return m_store[get_index(ind)];
}


template <typename T>
inline const T& LowerTriangularMatrix<T>::operator()(int i0, int i1) const
{
 return m_store[get_index(i0,i1)];
}

template <typename T>
inline const T& LowerTriangularMatrix<T>::operator()(uint i0, uint i1) const
{
 return m_store[get_index(i0,i1)];
}

template <typename T>
inline const T& LowerTriangularMatrix<T>::operator()(const std::vector<int>& ind) const
{
 return m_store[get_index(ind)];
}

template <typename T>
inline const T& LowerTriangularMatrix<T>::operator()(const std::vector<uint>& ind) const
{
 return m_store[get_index(ind)];
}



template <typename T>
inline void LowerTriangularMatrix<T>::put(int i0, int i1, const T& inputvalue)
{
 m_store[get_index(i0,i1)]=inputvalue;
}

template <typename T>
inline void LowerTriangularMatrix<T>::put(uint i0, uint i1, const T& inputvalue)
{
 m_store[get_index(i0,i1)]=inputvalue;
}

template <typename T>
inline void LowerTriangularMatrix<T>::put(const std::vector<int>& ind, const T& inputvalue)
{
 m_store[get_index(ind)]=inputvalue;
}

template <typename T>
inline void LowerTriangularMatrix<T>::put(const std::vector<uint>& ind, const T& inputvalue)
{
 m_store[get_index(ind)]=inputvalue;
}


template <typename T>
inline const T& LowerTriangularMatrix<T>::get(int i0, int i1) const
{
 return m_store[get_index(i0,i1)];
}

template <typename T>
inline const T& LowerTriangularMatrix<T>::get(uint i0, uint i1) const
{
 return m_store[get_index(i0,i1)];
}

template <typename T>
inline const T& LowerTriangularMatrix<T>::get(const std::vector<int>& ind) const
{
 return m_store[get_index(ind)];
}

template <typename T>
inline const T& LowerTriangularMatrix<T>::get(const std::vector<uint>& ind) const
{
 return m_store[get_index(ind)];
}



template <typename T>
LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::resize()
{
 return clear();
}

template <typename T>
LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::resize(int insize)
{
 m_store.resize(set_size(insize));
 return *this;
}


template <typename T>
LowerTriangularMatrix<T>& LowerTriangularMatrix<T>::resize(uint insize)
{
 m_store.resize(set_size(insize));
 return *this;
}


// **************************************************************


class RealSymmetricMatrix {

   std::vector<double> m_store;
   unsigned int m_size;

 public:

   RealSymmetricMatrix();
   RealSymmetricMatrix(int insize);
   RealSymmetricMatrix(uint insize);
   RealSymmetricMatrix(int insize, double initial_value);
   RealSymmetricMatrix(uint insize, double initial_value);
   RealSymmetricMatrix(const RealSymmetricMatrix& incoming);

   ~RealSymmetricMatrix();
   RealSymmetricMatrix& clear(); 

   RealSymmetricMatrix& operator=(double val);
   RealSymmetricMatrix& operator+=(double val);
   RealSymmetricMatrix& operator-=(double val);
   RealSymmetricMatrix& operator*=(double val);
   RealSymmetricMatrix& operator/=(double val);
   RealSymmetricMatrix& operator=(const RealSymmetricMatrix& incoming); 
   RealSymmetricMatrix& operator+=(const RealSymmetricMatrix& incoming);
   RealSymmetricMatrix& operator-=(const RealSymmetricMatrix& incoming);

   double& operator()(int i0, int i1);    
   double& operator()(uint i0, uint i1);    
   double& operator()(int i0, uint i1) {return operator()(i0,int(i1));}    
   double& operator()(uint i0, int i1) {return operator()(int(i0),i1);}
   double& operator()(const std::vector<int>& ind);             
   double& operator()(const std::vector<unsigned int>& ind);
   const double operator()(int i0, int i1) const; 
   const double operator()(uint i0, uint i1) const; 
   const double operator()(int i0, uint i1) const  {return operator()(i0,int(i1));} 
   const double operator()(uint i0, int i1) const  {return operator()(int(i0),i1);}
   const double operator()(const std::vector<int>& ind) const; 
   const double operator()(const std::vector<unsigned int>& ind) const;

   void put(int i0, int i1, double inputvalue);    
   void put(uint i0, uint i1, double inputvalue);    
   void put(int i0, uint i1, double inputvalue) {put(i0,int(i1),inputvalue);}    
   void put(uint i0, int i1, double inputvalue) {put(int(i0),i1,inputvalue);}
   void put(const std::vector<int>& ind, double inputvalue);             
   void put(const std::vector<unsigned int>& ind, double inputvalue);
   const double get(int i0, int i1) const; 
   const double get(uint i0, uint i1) const; 
   const double get(int i0, uint i1) const {return get(i0,int(i1));}
   const double get(uint i0, int i1) const {return get(int(i0),i1);}
   const double get(const std::vector<int>& ind) const; 
   const double get(const std::vector<unsigned int>& ind) const;

   uint size() const {return m_size;}

   RealSymmetricMatrix& resize();
   RealSymmetricMatrix& resize(int size);
   RealSymmetricMatrix& resize(uint size);


 private:

   unsigned int get_index(int ind0, int ind1) const;
   unsigned int get_index(uint ind0, uint ind1) const;
   unsigned int get_index(const std::vector<int>& ind) const;
   unsigned int get_index(const std::vector<uint>& ind) const;

   unsigned int set_size(int insize);
   unsigned int set_size(uint insize);

};


inline unsigned int RealSymmetricMatrix::set_size(int insize)
{
 if (insize<0) throw(std::invalid_argument("Invalid RealSymmetricMatrix size"));
 m_size=insize;
 return (m_size*(m_size+1))/2;
}


inline unsigned int RealSymmetricMatrix::set_size(uint insize)
{
 m_size=insize;
 return (m_size*(m_size+1))/2;
}


inline unsigned int RealSymmetricMatrix::get_index(int i0, int i1) const
{
#ifdef SAFETY_FLAG
 if ((i0<0)||(i0>=int(m_size))) throw(std::invalid_argument("RealSymmetricMatrix: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_size))) throw(std::invalid_argument("RealSymmetricMatrix: 2nd index out of range"));
#endif
 return (i0<=i1) ? (i1*(i1+1))/2+i0 : (i0*(i0+1))/2+i1;
}


inline unsigned int RealSymmetricMatrix::get_index(uint i0, uint i1) const
{
#ifdef SAFETY_FLAG
 if (i0>=m_size) throw(std::invalid_argument("RealSymmetricMatrix: 1st index out of range"));
 if (i1>=m_size) throw(std::invalid_argument("RealSymmetricMatrix: 2nd index out of range"));
#endif
 return (i0<=i1) ? (i1*(i1+1))/2+i0 : (i0*(i0+1))/2+i1;
}


inline unsigned int RealSymmetricMatrix::get_index(const std::vector<int>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("RealSymmetricMatrix: improper rank of RealSymmetricMatrix indices"));
#endif
 return get_index(ind[0],ind[1]);
}


inline unsigned int RealSymmetricMatrix::get_index(const std::vector<uint>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("RealSymmetricMatrix: improper rank of RealSymmetricMatrix indices"));
#endif
 return get_index(ind[0],ind[1]);
}



inline double& RealSymmetricMatrix::operator()(int i0, int i1)
{
 return m_store[get_index(i0,i1)];
}

inline double& RealSymmetricMatrix::operator()(uint i0, uint i1)
{
 return m_store[get_index(i0,i1)];
}

inline double& RealSymmetricMatrix::operator()(const std::vector<int>& ind)
{
 return m_store[get_index(ind)];
}

inline double& RealSymmetricMatrix::operator()(const std::vector<uint>& ind)
{
 return m_store[get_index(ind)];
}


inline const double RealSymmetricMatrix::operator()(int i0, int i1) const
{
 return m_store[get_index(i0,i1)];
}

inline const double RealSymmetricMatrix::operator()(uint i0, uint i1) const
{
 return m_store[get_index(i0,i1)];
}

inline const double RealSymmetricMatrix::operator()(const std::vector<int>& ind) const
{
 return m_store[get_index(ind)];
}

inline const double RealSymmetricMatrix::operator()(const std::vector<uint>& ind) const
{
 return m_store[get_index(ind)];
}




inline void RealSymmetricMatrix::put(int i0, int i1, double inputvalue)
{
 m_store[get_index(i0,i1)]=inputvalue;
}

inline void RealSymmetricMatrix::put(uint i0, uint i1, double inputvalue)
{
 m_store[get_index(i0,i1)]=inputvalue;
}

inline void RealSymmetricMatrix::put(const std::vector<int>& ind, double inputvalue)
{
 m_store[get_index(ind)]=inputvalue;
}

inline void RealSymmetricMatrix::put(const std::vector<uint>& ind, double inputvalue)
{
 m_store[get_index(ind)]=inputvalue;
}


inline const double RealSymmetricMatrix::get(int i0, int i1) const
{
 return m_store[get_index(i0,i1)];
}

inline const double RealSymmetricMatrix::get(uint i0, uint i1) const
{
 return m_store[get_index(i0,i1)];
}

inline const double RealSymmetricMatrix::get(const std::vector<int>& ind) const
{
 return m_store[get_index(ind)];
}

inline const double RealSymmetricMatrix::get(const std::vector<uint>& ind) const
{
 return m_store[get_index(ind)];
}



// **************************************************************




class ComplexHermitianMatrix {

   std::vector<std::complex<double> > m_store;
   unsigned int m_size;

 public:

   ComplexHermitianMatrix();
   ComplexHermitianMatrix(int size);
   ComplexHermitianMatrix(uint size);
   ComplexHermitianMatrix(int size, double initial_value);
   ComplexHermitianMatrix(uint size, double initial_value);
   ComplexHermitianMatrix(const ComplexHermitianMatrix& incoming);

   ~ComplexHermitianMatrix();
   ComplexHermitianMatrix& clear(); 

   ComplexHermitianMatrix& operator=(double val);
   ComplexHermitianMatrix& operator+=(double val);
   ComplexHermitianMatrix& operator-=(double val);
   ComplexHermitianMatrix& operator*=(double val);
   ComplexHermitianMatrix& operator/=(double val);
   ComplexHermitianMatrix& operator=(const ComplexHermitianMatrix& incoming); 
   ComplexHermitianMatrix& operator+=(const ComplexHermitianMatrix& incoming);
   ComplexHermitianMatrix& operator-=(const ComplexHermitianMatrix& incoming);

   const std::complex<double> operator()(int i0, int i1) const; 
   const std::complex<double> operator()(uint i0, uint i1) const; 
   const std::complex<double> operator()(int i0, uint i1) const {return operator()(i0,int(i1));}
   const std::complex<double> operator()(uint i0, int i1) const {return operator()(int(i0),i1);}
   const std::complex<double> operator()(const std::vector<int>& ind) const; 
   const std::complex<double> operator()(const std::vector<unsigned int>& ind) const;

   void put(int i0, int i1, const std::complex<double>& inputvalue);    
   void put(uint i0, uint i1, const std::complex<double>& inputvalue);    
   void put(int i0, uint i1, const std::complex<double>& inputvalue) {put(i0,int(i1),inputvalue);}   
   void put(uint i0, int i1, const std::complex<double>& inputvalue) {put(int(i0),i1,inputvalue);}   
   void put(const std::vector<int>& ind, const std::complex<double>& inputvalue);             
   void put(const std::vector<unsigned int>& ind, const std::complex<double>& inputvalue);

   const std::complex<double> get(int i0, int i1) const; 
   const std::complex<double> get(uint i0, uint i1) const; 
   const std::complex<double> get(int i0, uint i1) const {return get(i0,int(i1));}
   const std::complex<double> get(uint i0, int i1) const {return get(int(i0),i1);}
   const std::complex<double> get(const std::vector<int>& ind) const; 
   const std::complex<double> get(const std::vector<unsigned int>& ind) const;

   uint size() const {return m_size;}

   ComplexHermitianMatrix& resize();
   ComplexHermitianMatrix& resize(int size);
   ComplexHermitianMatrix& resize(uint size);


 private:

   unsigned int get_index(int ind0, int ind1) const;
   unsigned int get_index(uint ind0, uint ind1) const;
   unsigned int set_size(int insize);
   unsigned int set_size(uint insize);


};


inline unsigned int ComplexHermitianMatrix::set_size(int insize)
{
 if (insize<0) throw(std::invalid_argument("Invalid ComplexHermitianMatrix size"));
 m_size=insize;
 return (m_size*(m_size+1))/2;
}

inline unsigned int ComplexHermitianMatrix::set_size(uint insize)
{
 m_size=insize;
 return (m_size*(m_size+1))/2;
}

                 //  (i0<=i1) assumed below
inline unsigned int ComplexHermitianMatrix::get_index(int i0, int i1) const
{
#ifdef SAFETY_FLAG
 if ((i0<0)||(i0>=int(m_size))) throw(std::invalid_argument("ComplexHermitianMatrix: 1st index out of range"));
 if ((i1<0)||(i1>=int(m_size))) throw(std::invalid_argument("ComplexHermitianMatrix: 2nd index out of range"));
#endif
 return (i1*(i1+1))/2+i0;
}

                 //  (i0<=i1) assumed below
inline unsigned int ComplexHermitianMatrix::get_index(uint i0, uint i1) const
{
#ifdef SAFETY_FLAG
 if (i0>=m_size) throw(std::invalid_argument("ComplexHermitianMatrix: 1st index out of range"));
 if (i1>=m_size) throw(std::invalid_argument("ComplexHermitianMatrix: 2nd index out of range"));
#endif
 return (i1*(i1+1))/2+i0;
}



inline const std::complex<double> ComplexHermitianMatrix::operator()(int i0, int i1) const
{
 return get(i0,i1);
}

inline const std::complex<double> ComplexHermitianMatrix::operator()(uint i0, uint i1) const
{
 return get(i0,i1);
}

inline const std::complex<double> ComplexHermitianMatrix::operator()(const std::vector<int>& ind) const
{
 return get(ind);
}

inline const std::complex<double> ComplexHermitianMatrix::operator()(const std::vector<uint>& ind) const
{
 return get(ind);
}




inline void ComplexHermitianMatrix::put(int i0, int i1, const std::complex<double>& inputvalue)
{
 if ((i0==i1)&&(inputvalue.imag()!=0.0))
    throw(std::invalid_argument("ComplexHermitianMatrix: nonzero imaginary part in diagonal element assignment"));
 if (i0<=i1)
    m_store[get_index(i0,i1)]=inputvalue;
 else
    m_store[get_index(i1,i0)]=conjugate(inputvalue);
}


inline void ComplexHermitianMatrix::put(uint i0, uint i1, const std::complex<double>& inputvalue)
{
 if ((i0==i1)&&(inputvalue.imag()!=0.0))
    throw(std::invalid_argument("ComplexHermitianMatrix: nonzero imaginary part in diagonal element assignment"));
 if (i0<=i1)
    m_store[get_index(i0,i1)]=inputvalue;
 else
    m_store[get_index(i1,i0)]=conjugate(inputvalue);
}


inline void ComplexHermitianMatrix::put(const std::vector<uint>& ind, const std::complex<double>& inputvalue)
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("ComplexHermitianMatrix: improper rank of ComplexHermitianMatrix indices"));
#endif
 put(ind[0],ind[1],inputvalue);
}

inline void ComplexHermitianMatrix::put(const std::vector<int>& ind, const std::complex<double>& inputvalue)
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("ComplexHermitianMatrix: improper rank of ComplexHermitianMatrix indices"));
#endif
 put(ind[0],ind[1],inputvalue);
}




inline const std::complex<double> ComplexHermitianMatrix::get(int i0, int i1) const
{
 return (i0<=i1) ? m_store[get_index(i0,i1)] : conjugate(m_store[get_index(i1,i0)]);
}

inline const std::complex<double> ComplexHermitianMatrix::get(uint i0, uint i1) const
{
 return (i0<=i1) ? m_store[get_index(i0,i1)] : conjugate(m_store[get_index(i1,i0)]);
}

inline const std::complex<double> ComplexHermitianMatrix::get(const std::vector<int>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("ComplexHermitianMatrix: improper rank of ComplexHermitianMatrix indices"));
#endif
 return (ind[0]<=ind[1]) ? m_store[get_index(ind[0],ind[1])] : conjugate(m_store[get_index(ind[1],ind[0])]);
}

inline const std::complex<double> ComplexHermitianMatrix::get(const std::vector<uint>& ind) const
{
#ifdef SAFETY_FLAG
 if (ind.size()!=2) throw(std::invalid_argument("ComplexHermitianMatrix: improper rank of ComplexHermitianMatrix indices"));
#endif
 return (ind[0]<=ind[1]) ? m_store[get_index(ind[0],ind[1])] : conjugate(m_store[get_index(ind[1],ind[0])]);
}



// **************************************************************

typedef Vector<double>                  RVector;
typedef Vector<std::complex<double> >   CVector;

typedef Matrix<double>                  RMatrix;
typedef Matrix<std::complex<double> >   CMatrix;

typedef LowerTriangularMatrix<double>   LTMatrix;

#if defined COMPLEXNUMBERS
typedef ComplexHermitianMatrix          HermMatrix;
typedef Matrix<std::complex<double> >   TransMatrix;
typedef Vector<std::complex<double> >   SVector;

#elif defined REALNUMBERS
typedef RealSymmetricMatrix             HermMatrix;
typedef Matrix<double>                  TransMatrix;
typedef Vector<double>                  SVector;

#else
  #error "Either COMPLEXNUMBERS or REALNUMBERS must be defined"
#endif

// **************************************************************


RVector copyRVectorWithoutZerothElement(const RVector& in);


// **************************************************************
#endif
