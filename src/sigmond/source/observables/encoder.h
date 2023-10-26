#ifndef ENCODER_H
#define ENCODER_H

#include <string>
#include <vector>


// *********************************************************************
// *                                                                   *
// *   These utility routines convert a string to a vector of          *
// *   unsigned integers and vice versa.  A "maxchar" is required.     *
// *   The string cannot contain any white space and cannot be empty.  *
// *                                                                   *
// *********************************************************************


void encode_string_to_uints(const std::string& astr, unsigned int maxchar,
                            std::vector<unsigned int>& icode);

std::string decode_uints_to_string(const std::vector<unsigned int>& icode);


#endif
