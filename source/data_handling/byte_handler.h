#ifndef BYTEORDER_H
#define BYTEORDER_H

#include<vector>
#include <cstdio>

// **************************************************
// *                                                *
// *   This class handles byte swapping to          *
// *   convert endian, and evaluating check sums.   *
// *                                                *
// *    1  = little-endian    (alpha, x86_64, ...)  *
// *    2  = big-endian       (sun, ibm, hp, ...)   *
// *                                                *
// **************************************************

class ByteHandler
{
       // disable copying
    ByteHandler(const ByteHandler&);
    ByteHandler& operator=(const ByteHandler&);
    
 public:

    typedef unsigned int        n_uint32_t;
    typedef unsigned short int  n_uint16_t;    


    ByteHandler();

         // Is the native byte order big endian?
    bool big_endian();

         // Byte-swap an array of data each of size nmemb
    void byte_swap(void *ptr, size_t size, size_t nmemb);

         // Return a check sum given an input checksum and buffer data
    n_uint32_t get_checksum(n_uint32_t crc, const unsigned char *buf, size_t len);
 
    n_uint32_t get_checksum(n_uint32_t crc, const char *buf, size_t len);


 private:
 
    std::vector<n_uint32_t> crc_table;

};

// *************************************************************
#endif
