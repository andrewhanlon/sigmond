#include "matrix.h"


// **************************************************************


RealSymmetricMatrix::RealSymmetricMatrix()
{m_size=0;}


RealSymmetricMatrix::RealSymmetricMatrix(int insize)  
                   : m_store(set_size(insize))
{}

RealSymmetricMatrix::RealSymmetricMatrix(uint insize)  
                   : m_store(set_size(insize))
{}

RealSymmetricMatrix::RealSymmetricMatrix(int insize, double initial_value) 
                   : m_store(set_size(insize),initial_value)
{}

RealSymmetricMatrix::RealSymmetricMatrix(uint insize, double initial_value) 
                   : m_store(set_size(insize),initial_value)
{}

RealSymmetricMatrix::RealSymmetricMatrix(const RealSymmetricMatrix& incoming) 
         : m_store(incoming.m_store), m_size(incoming.m_size)
{}


RealSymmetricMatrix::~RealSymmetricMatrix()
{}


RealSymmetricMatrix& RealSymmetricMatrix::clear()
{
 m_store.clear();
 m_size=0;
 return *this;
}



RealSymmetricMatrix& RealSymmetricMatrix::operator=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]=val;
 return *this;
}


RealSymmetricMatrix& RealSymmetricMatrix::operator+=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]+=val;
 return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator-=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]-=val;
 return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator*=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]*=val;
 return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator/=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]/=val;
 return *this;
}


RealSymmetricMatrix& RealSymmetricMatrix::operator=(const RealSymmetricMatrix& incoming)
{
 if (this==&incoming) return *this;
 m_store=incoming.m_store;
 m_size=incoming.m_size;
 return *this;
}



RealSymmetricMatrix& RealSymmetricMatrix::operator+=(const RealSymmetricMatrix& incoming)
{
#ifdef SAFETY_FLAG
 if (m_size!=incoming.m_size)
    throw(std::invalid_argument("RealSymmetricMatrix size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]+=incoming.m_store[i];
 return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator-=(const RealSymmetricMatrix& incoming)
{
#ifdef SAFETY_FLAG
 if (m_size!=incoming.m_size)
    throw(std::invalid_argument("RealSymmetricMatrix size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]-=incoming.m_store[i];
 return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::resize()
{
 return clear();
}

RealSymmetricMatrix& RealSymmetricMatrix::resize(int insize)
{
 m_store.resize(set_size(insize));
 return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::resize(uint insize)
{
 m_store.resize(set_size(insize));
 return *this;
}


// **************************************************************


ComplexHermitianMatrix::ComplexHermitianMatrix()
{m_size=0;}


ComplexHermitianMatrix::ComplexHermitianMatrix(int insize)  
                   : m_store(set_size(insize))
{}

ComplexHermitianMatrix::ComplexHermitianMatrix(uint insize)  
                   : m_store(set_size(insize))
{}

ComplexHermitianMatrix::ComplexHermitianMatrix(int insize, double initial_value) 
                   : m_store(set_size(insize),std::complex<double>(initial_value,0.0))
{}

ComplexHermitianMatrix::ComplexHermitianMatrix(uint insize, double initial_value) 
                   : m_store(set_size(insize),std::complex<double>(initial_value,0.0))
{}

ComplexHermitianMatrix::ComplexHermitianMatrix(const ComplexHermitianMatrix& incoming) 
         : m_store(incoming.m_store), m_size(incoming.m_size)
{}


ComplexHermitianMatrix::~ComplexHermitianMatrix()
{}


ComplexHermitianMatrix& ComplexHermitianMatrix::clear()
{
 m_store.clear();
 m_size=0;
 return *this;
}



ComplexHermitianMatrix& ComplexHermitianMatrix::operator=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]=val;
 return *this;
}


ComplexHermitianMatrix& ComplexHermitianMatrix::operator+=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]+=val;
 return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator-=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]-=val;
 return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator*=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]*=val;
 return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator/=(double val)
{
 for (uint i=0;i<m_store.size();++i) m_store[i]/=val;
 return *this;
}


ComplexHermitianMatrix& ComplexHermitianMatrix::operator=(const ComplexHermitianMatrix& incoming)
{
 if (this==&incoming) return *this;
 m_store=incoming.m_store;
 m_size=incoming.m_size;
 return *this;
}



ComplexHermitianMatrix& ComplexHermitianMatrix::operator+=(const ComplexHermitianMatrix& incoming)
{
#ifdef SAFETY_FLAG
 if (m_size!=incoming.m_size)
    throw(std::invalid_argument("ComplexHermitianMatrix size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]+=incoming.m_store[i];
 return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator-=(const ComplexHermitianMatrix& incoming)
{
#ifdef SAFETY_FLAG
 if (m_size!=incoming.m_size)
    throw(std::invalid_argument("ComplexHermitianMatrix size mismatch"));
#endif
 for (uint i=0;i<m_store.size();++i) m_store[i]-=incoming.m_store[i];
 return *this;
}



ComplexHermitianMatrix& ComplexHermitianMatrix::resize()
{
 return clear();
}

ComplexHermitianMatrix& ComplexHermitianMatrix::resize(int insize)
{
 m_store.resize(set_size(insize));
 return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::resize(uint insize)
{
 m_store.resize(set_size(insize));
 return *this;
}


// ***************************************************************************


RVector copyRVectorWithoutZerothElement(const RVector& in)
{
 uint n=in.size()-1;
 RVector buffer(n);
 for (uint k=0;k<n;k++)
    buffer[k]=in[k+1];
 return buffer;
}


// **************************************************************
