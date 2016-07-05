#ifndef VEC_UTILS_H
#define VEC_UTILS_H
#include <list>
#include <vector>

// *******************************************************************************
// *                                                                             *
// *    The following utility routines are also defined here:                    *
// *          "make_list",  "make_vector"                                        *
// *                                                                             *
// *******************************************************************************

  //   Simple routine to help make a list of objects

template <typename T>
std::list<T> make_list()
{
 std::list<T> result;
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1)
{
 std::list<T> result;
 result.push_back(v1);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 result.push_back(v6);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6, const T& v7)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 result.push_back(v6);
 result.push_back(v7);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6, const T& v7, const T& v8)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 result.push_back(v6);
 result.push_back(v7);
 result.push_back(v8);
 return result;
}

  //   Simple routine to help make a vector of objects

template <typename T>
std::vector<T> make_vector()
{
 std::vector<T> result;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1)
{
 std::vector<T> result(1);
 result[0]=v1;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2)
{
 std::vector<T> result(2);
 result[0]=v1;
 result[1]=v2;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3)
{
 std::vector<T> result(3);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4)
{
 std::vector<T> result(4);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5)
{
 std::vector<T> result(5);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6)
{
 std::vector<T> result(6);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 result[5]=v6;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6, const T& v7)
{
 std::vector<T> result(7);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 result[5]=v6;
 result[6]=v7;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6, const T& v7, const T& v8)
{
 std::vector<T> result(8);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 result[5]=v6;
 result[6]=v7;
 result[7]=v8;
 return result;
}

 
// ****************************************************************************

#endif
