#ifndef ENDIAN_H
#define ENDIAN_H

#include <cstdint>

enum Endian { LittleEndian, BigEndian, UnknownEndian };

Endian endianness()
{
 union {
    uint32_t i;
    char c[4];
 } bint = {0x01020304};

 if (bint.c[0] == 1u)
    return BigEndian;
 else if (bint.c[0] == 4u)
    return LittleEndian;
 else
    return UnknownEndian;
}

#endif
