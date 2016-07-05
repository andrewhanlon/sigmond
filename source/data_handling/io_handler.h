#ifndef IO_HANDLER_H
#define IO_HANDLER_H

#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <complex>
#include "qdp_byteorder.h"
#include "array.h"

#ifndef NO_CXX11
#include <type_traits>
#endif

namespace LaphEnv {


typedef std::complex<double> dcmplx;
typedef std::complex<float>  fcmplx;


 // *********************************************************************************
 // *                                                                               *
 // *       class IOHandler:     random access input/output                         *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *   This class provides binary input/output. It allows the user to choose       *
 // *   between big-endian and little-endian format, and it allows the user to turn *
 // *   check sums on or off.  Files written by this class will start with a single *
 // *   character 'B' or 'L' to indicate endian-ness.  Then an ID string is output, *
 // *   which always has the length "ID_string_length"=32 (spaces are padded if     *
 // *   needed).  Files read by this class expect this starting behavior.  An error *
 // *   occurs if the ID string given in the open command on an existing file does  *
 // *   not match that in the file.                                                 *
 // *                                                                               *
 // *   All errors are considered fatal in this class, so the end user need not     *
 // *   perform any error checking.  The ID string is useful for ensuring that      *
 // *   the file contains the kind of information that the user is expecting.       *
 // *                                                                               *
 // *   Files can be opened in one of four modes:                                   *
 // *     (1) ReadOnly  -- fails if the file does not exist or read not allowed     *
 // *     (2) ReadWriteFailIfExists --  fails if the file exists                    *
 // *     (3) ReadWriteEraseIfExists  -- erases the file if it exists               *
 // *     (4) ReadWriteUpdateIfExists  -- updates existing file or creates new      *
 // *   Modes 2,3,4 will create a new file if it does not exist.  When creating     *
 // *   a new file, the endian format choices are 'B', 'L', or 'N' (native).        *
 // *   Random access to the file is employed.                                      *
 // *                                                                               *
 // *   There are separate open routines for the different modes, or you can        *
 // *   use the OpenMode enum in the general open routine.  Alternatively, you      *
 // *   can use iostream openmodes as follows:                                      *
 // *      ReadOnly                  =   ios::in                                    *
 // *      ReadWriteFailIfExists     =   ios::in | ios::out                         *
 // *      ReadWriteEraseIfExists    =   ios::in | ios::out | ios::trunc            *
 // *                                                                               *
 // *   There are seek and tell members for moving around in the file.  It is an    *
 // *   error to seek before the start of the data, but you can seek past the end   *
 // *   of the file.  You can seek relative to the start of the data (33 characters *
 // *   past the start of the file), the current location, or the end of file.      *
 // *   Use negative offsets relative to the end of file to go backwards.           *
 // *                                                                               *
 // *   Input is done with read(..) commands, and output is done with write(..)     *
 // *   commands.  The data types currently supported for read/write are            *
 // *                                                                               *
 // *       - basic data types (int, bool, float, string, etc.)                     *
 // *       - complex numbers using the STL complex class                           *
 // *       - vector, Array of basic/complex data types                             *
 // *                                                                               *
 // *   Quantities above are contiguous in memory and regular (each element is      *
 // *   the same size) and these facts have been used to speed up I/O.  The file    *
 // *   pointer is advanced by the size of the object after the I/O operation.      *
 // *                                                                               *
 // *   Examples of writing and reading:                                            *
 // *                                                                               *
 // *      IOHandler io;                                                            *
 // *      io.open(....);                                                           *
 // *      int k=5;     io.seekFromStart(...); write(io,k);                         *
 // *      float x=5.4; io.seekFromStart(...); write(io,x);                         *
 // *      int j; float y;  io.seekFromStart(...); read(io,j,y); // in sequence     *
 // *      string str("ABC"); io.seekFromStart(...); write(io,str);                 *
 // *                                                                               *
 // *   If check sums are turned on, objects of this class maintain a check sum.    *
 // *   Doing an explicit seek resets the checksum; changing from a read to a       *
 // *   write or vice versa resets the checksum; or an explicit reset can also      *
 // *   be done.  The checksum is updated as bytes are read or written.             *
 // *   Successive reads or successive writes update the checksum.                  *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************

 
    //  The C++ standard requires that elements of a "vector" are
    //  stored in contiguous memory, and this class uses this fact.  A general "struct"
    //  may or may not be stored in contiguous memory; padding characters
    //  could be inserted between the struct elements (the GNU compiler
    //  allows an "__attribute__ ((packed))" to prevent such padding).
    //  This is why we must treat complex numbers differently than basic floats, etc.


class IOHandler
{

   std::fstream fh;

   bool read_only;
   bool openflag;
   bool read_mode;
   
   char endian_format;            // 'B' for big-endian, 'L' for little-endian
   bool endian_convert;

   bool checksum_on;
   std::string m_filename;
   bool is_new_file;
 
   QDPUtil::n_uint32_t checksum;

      // disallow copying
   IOHandler(const IOHandler&);
   IOHandler(IOHandler&);
   IOHandler& operator=(const IOHandler&);
   IOHandler& operator=(IOHandler&);


 public:
 
   enum OpenMode { ReadOnly, ReadWriteFailIfExists, ReadWriteEraseIfExists, 
                   ReadWriteUpdateIfExists };

   explicit IOHandler();

   IOHandler(const std::string& filename, OpenMode mode=ReadOnly,
             const std::string& filetype_id="", char endianness='N',
             bool turn_on_checksum=false);

   IOHandler(const std::string& filename, std::ios_base::openmode mode,
             const std::string& filetype_id="", char endianness='N',
             bool turn_on_checksum=false);

   void open(const std::string& filename, OpenMode mode=ReadOnly,
             const std::string& filetype_id="", char endianness='N',
             bool turn_on_checksum=false);

   void open(const std::string& filename, std::ios_base::openmode mode,
             const std::string& filetype_id="", char endianness='N',
             bool turn_on_checksum=false);

   void openReadOnly(const std::string& filename, const std::string& filetype_id="",
                     bool turn_on_checksum=false);

   void openNew(const std::string& filename, bool fail_if_exists=true,
                const std::string& filetype_id="", char endianness='N',
                bool turn_on_checksum=false);

   void openUpdate(const std::string& filename, const std::string& filetype_id="", 
                   char endianness='N', bool turn_on_checksum=false);

   ~IOHandler();

          // closes current file if open, otherwise no action taken
   void close();


          // informational routines
          
   bool isOpen() const { return openflag; }
   
   bool isNewFile() const { return is_new_file; }

   std::string getFileName() const { return m_filename; }
   
   bool isChecksumOn() const { return checksum_on; }
   
   bool isEndianConversionOn() const { return endian_convert; }
   
   bool isFileLittleEndian() const { return (endian_format=='L'); }
   
   bool isFileBigEndian() const { return (endian_format=='B'); }

   bool isReadOnly() const { return read_only; }
   

   typedef std::iostream::pos_type   pos_type;  // position in buffer
   typedef std::iostream::off_type   off_type;  // offset in buffer
   typedef std::iostream::seekdir    whence_type;
   typedef std::iostream::openmode   openmode_type;

         //  Set the file pointer relative to start of data in file,
         //  end of file (use negative to go backward), or current 
         //  location.  Checksum is reset if in use.

         //  Caution: make sure to convert sizeof(...) quantities
         //  to an IOHandler::off_type(sizeof(...)) if you wish to use 
         //  negative offsets.  sizeof(...) returns an unsigned 
         //  integer type.

   void seekFromStart(off_type offset);  // start means start of data (not file)
   void seekFromCurr(off_type offset);
   void seekFromEnd(off_type offset);
    
   void seek(pos_type offset);          // from start of data
   void seekBegin(off_type offset);     // from start of data
   void seekRelative(off_type offset);
   void seekEnd(off_type offset);
   void rewind();   // puts pointer at start of data in file

         //  Getting the file pointer location in bytes from start of data

   pos_type tell();
   pos_type currentPosition();


         //  Check summing
   
   void turnOnChecksum();
   void turnOffChecksum();
   void resetChecksum();   
   QDPUtil::n_uint32_t getChecksum();
   void printFileID();

     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened and its file type matches 
     // "filetype_id", returns false otherwise.  The string is returned in 
     // "stringvalue".  This routine is used by objects in data_io_handler.h 
     // when the multi-file handlers build up their maps of file keys.
           
   bool peekString(std::string& stringvalue, unsigned int position, 
                   const std::string& filename,
                   const std::string& filetype_id="");

     // Open file, read ID string at beginning of file, then close. Returns
     // true if file can be opened and read; false other.
     // Returns ID string in "stringvalue".

   bool peekID(std::string& stringvalue, const std::string& filename);

         // write routines  

   void write(const std::string& output);
   void write(const char& output);
   void write(const int& output);
   void write(const unsigned int& output);
   void write(const long int& output);
   void write(const unsigned long int& output);
   void write(const long long int& output);
   void write(const unsigned long long int& output);
   void write(const float& output);
   void write(const double& output);
   void write(const bool& output);
   void write(const fcmplx& output);
   void write(const dcmplx& output);
   void write(const std::vector<char>& output);
   void write(const std::vector<int>& output);
   void write(const std::vector<unsigned int>& output);
   void write(const std::vector<long int>& output);
   void write(const std::vector<unsigned long int>& output);
   void write(const std::vector<float>& output);
   void write(const std::vector<double>& output);
   void write(const std::vector<fcmplx>& output);
   void write(const std::vector<dcmplx>& output);
   void write(const Array<char>& output);
   void write(const Array<int>& output);
   void write(const Array<unsigned int>& output);
   void write(const Array<long int>& output);
   void write(const Array<unsigned long int>& output);
   void write(const Array<float>& output);
   void write(const Array<double>& output);
   void write(const Array<fcmplx>& output);
   void write(const Array<dcmplx>& output);
   void multi_write(const char* output, int n);
   void multi_write(const int* output, int n);
   void multi_write(const unsigned int* output, int n);
   void multi_write(const long int* output, int n);
   void multi_write(const unsigned long int* output, int n);
   void multi_write(const long long int* output, int n);
   void multi_write(const unsigned long long int* output, int n);
   void multi_write(const float* output, int n);
   void multi_write(const double* output, int n);
   void multi_write(const bool* output, int n);
   void multi_write(const fcmplx* output, int n);
   void multi_write(const dcmplx* output, int n);


          // read routines

   void read(std::string& input);
   void read(char& input);
   void read(int& input);
   void read(unsigned int& input);
   void read(long int& input);
   void read(unsigned long int& input);
   void read(long long int& input);
   void read(unsigned long long int& input);
   void read(float& input);
   void read(double& input);
   void read(bool& input);
   void read(fcmplx& input);
   void read(dcmplx& input);
   void read(std::vector<char>& input);
   void read(std::vector<int>& input);
   void read(std::vector<unsigned int>& input);
   void read(std::vector<long int>& input);
   void read(std::vector<unsigned long int>& input);
   void read(std::vector<float>& input);
   void read(std::vector<double>& input);
   void read(std::vector<fcmplx>& input);
   void read(std::vector<dcmplx>& input);
   void read(Array<char>& input);
   void read(Array<int>& input);
   void read(Array<unsigned int>& input);
   void read(Array<long int>& input);
   void read(Array<unsigned long int>& input);
   void read(Array<float>& input);
   void read(Array<double>& input);
   void read(Array<fcmplx>& input);
   void read(Array<dcmplx>& input);
             // routines below assume memory has already been allocated
   void multi_read(char* input, int n);
   void multi_read(int* input, int n);
   void multi_read(unsigned int* input, int n);
   void multi_read(long int* input, int n);
   void multi_read(unsigned long int* input, int n);
   void multi_read(long long int* input, int n);
   void multi_read(unsigned long long int* input, int n);
   void multi_read(float* input, int n);
   void multi_read(double* input, int n);
   void multi_read(bool* input, int n);
   void multi_read(fcmplx* input, int n);
   void multi_read(dcmplx* input, int n);
 

         //  number of bytes routines (will be useful by IOMap for
         //  verifying read/writes and determining whether overwrites
         //  should occur)

   size_t numbytes(const std::string& data) const;
   size_t numbytes(const char& data) const;
   size_t numbytes(const int& data) const;
   size_t numbytes(const unsigned int& data) const;
   size_t numbytes(const long int& data) const;
   size_t numbytes(const unsigned long int& data) const;
   size_t numbytes(const long long int& data) const;
   size_t numbytes(const unsigned long long int& data) const;
   size_t numbytes(const float& data) const;
   size_t numbytes(const double& data) const;
   size_t numbytes(const bool& data) const;
   size_t numbytes(const fcmplx& data) const;
   size_t numbytes(const dcmplx& data) const;
   size_t numbytes(const std::vector<char>& data) const;
   size_t numbytes(const std::vector<int>& data) const;
   size_t numbytes(const std::vector<unsigned int>& data) const;
   size_t numbytes(const std::vector<long int>& data) const;
   size_t numbytes(const std::vector<unsigned long int>& data) const;
   size_t numbytes(const std::vector<float>& data) const;
   size_t numbytes(const std::vector<double>& data) const;
   size_t numbytes(const std::vector<fcmplx>& data) const;
   size_t numbytes(const std::vector<dcmplx>& data) const;
   size_t numbytes(const Array<char>& data) const;
   size_t numbytes(const Array<int>& data) const;
   size_t numbytes(const Array<unsigned int>& data) const;
   size_t numbytes(const Array<long int>& data) const;
   size_t numbytes(const Array<unsigned long int>& data) const;
   size_t numbytes(const Array<float>& data) const;
   size_t numbytes(const Array<double>& data) const;
   size_t numbytes(const Array<fcmplx>& data) const;
   size_t numbytes(const Array<dcmplx>& data) const;


 private:

      // private utility routines
       
   void open_existing_file(const std::string& filetype_id,
                           openmode_type access_mode);

   void open_new_file(const std::string& filetype_id, char endianness);

   void clear();

   bool fileExists();

   void writeCommon(const char *data, size_t nbytes);

   void readCommon(char *data, size_t nbytes);

   void check_for_failure(int errcode, const std::string& mesg="");

   std::string int_to_string(int intval);

   std::string tidyString(const std::string& str);  

   void readIDstring(char &endianness, std::string& ID_string);

   void writeIDstring(const std::string& ID_string);

   void write_common(const char* output, size_t element_size, size_t nelements);

   void read_common(char* output, size_t element_size, size_t nelements);
  

   void delete_file();
   void file_open(openmode_type access_mode, const std::string& errmsg);
   void file_close();
   void file_seek(off_type offset, whence_type whence);
   pos_type file_tell();

   bool peeker(std::string& stringvalue, unsigned int position, 
               const std::string& filename,
               const std::string& filetype_id="");

       // constants
       
   static const int ID_string_length;
   static const pos_type data_start_pos;

   static const int IO_ERR_NO_SUCH_FILE;
   static const int IO_ERR_ACCESS;
   static const int IO_ERR_OTHER;
   static const int IO_SUCCESS;

   static const openmode_type IO_MODE_RDONLY;
   static const openmode_type IO_MODE_RDWR;
   static const openmode_type IO_MODE_CREATE;
   static const openmode_type IO_MODE_EXCL;

   static const whence_type IO_SEEK_BEG;
   static const whence_type IO_SEEK_CUR;
   static const whence_type IO_SEEK_END;

#ifdef NO_CXX11
          // for static (compile time) assertion
   template <bool b>
   void static__assert()
   { typedef char asserter[b?1:-1]; }
#endif

          // private write members

   template <typename T>
   void write_basic(const T& output);

   template <typename T>
   void write_complex(const std::complex<T>& output);

   template <typename T>
   void multi_write_basic(const T* output, int n);

   template <typename T>
   void multi_write_complex(const std::complex<T>* output, int n);

   template <typename T>
   void write_basic_vector(const std::vector<T>& output);

   template <typename T>
   void write_complex_vector(const std::vector<std::complex<T> >& output);

   template <typename T>
   void write_array(const Array<T>& output);


          // private read members

   template <typename T>
   void read_basic(T& output);

   template <typename T>
   void read_complex(std::complex<T>& output);

   template <typename T>
   void multi_read_basic(T* output, int n);

   template <typename T>
   void multi_read_complex(std::complex<T>* output, int n);

   template <typename T>
   void read_basic_vector(std::vector<T>& output);

   template <typename T>
   void read_complex_vector(std::vector<std::complex<T> >& output);

   template <typename T>
   void read_array(Array<T>& output);



          // private numbytes members

   template <typename T>
   size_t numbytes_basic(const T& data) const;

   template <typename T>
   size_t numbytes_complex(const std::complex<T>& data) const;

   template <typename T>
   size_t numbytes_basic_vector(const std::vector<T>& data) const;

   template <typename T>
   size_t numbytes_complex_vector(const std::vector<std::complex<T> >& data) const;

   template <typename T>
   size_t numbytes_array(const Array<T>& data) const;


};


 // ***************************************************************

            // private write members

template <typename T>
void IOHandler::write_basic(const T& output)
{ 
 write_common((const char*)&output, sizeof(T), 1);
}

template <typename T>
void IOHandler::write_complex(const std::complex<T>& output)
{
 T data[2]; data[0]=real(output); data[1]=imag(output);
 write_common((const char*)&data, sizeof(T), 2);
}

template <typename T>
void IOHandler::multi_write_basic(const T* output, int n)
{ 
 write_common((const char*)output, sizeof(T), n);
}

template <typename T>
void IOHandler::multi_write_complex(const std::complex<T>* output, int n)
{
 std::vector<T> data(2*n);
 int count=0;
 const std::complex<T> *ptr=output; 
 for (int k=0;k<n;++k){
    data[count++]=real(*ptr);
    data[count++]=imag(*ptr);
    ++ptr;}
 write_common((const char*)&data[0],sizeof(T),size_t(2*n));
}

template <typename T>
void IOHandler::write_basic_vector(const std::vector<T>& output)
{
 int n=output.size();        // vector guarantees contiguous memory allocation
 write_basic<int>(n);
 if (n==0) return;
 write_common((const char*)&(output[0]),sizeof(T),size_t(n)); 
}

template <typename T>
void IOHandler::write_complex_vector(const std::vector<std::complex<T> >& output)
{
 int n=output.size();
 write_basic<int>(n);
 if (n==0) return;
 std::vector<T> data(2*n);
 int count=0; for (int k=0;k<n;++k){
    data[count++]=real(output[k]);
    data[count++]=imag(output[k]);}
 write_common((const char*)&data[0],sizeof(T),size_t(2*n));
}

template <typename T>
void IOHandler::write_array(const Array<T>& output)
{
 unsigned int n=output.numDimensions();
 write_basic<unsigned int>(n);
 if (n==0) return;
 write(output.m_sizes);
 write(output.m_store);
}


// ************************************************************  
 
            // private read members

template <typename T>
void IOHandler::read_basic(T& input)
{ 
 read_common((char*)&input, sizeof(T), 1);
}

template <typename T>
void IOHandler::read_complex(std::complex<T>& input)
{
 T data[2]; 
 read_common((char*)&data, sizeof(T), 2);
 input=std::complex<T>(data[0],data[1]);
}

template <typename T>
void IOHandler::multi_read_basic(T* input, int n)
{ 
 read_common((char*)input, sizeof(T), n);
}

template <typename T>
void IOHandler::multi_read_complex(std::complex<T>* input, int n)
{
 std::vector<T> data(2*n);
 read_common((char*)&data[0],sizeof(T),size_t(2*n));
 int count=0;
 std::complex<T> *ptr=input; 
 for (int k=0;k<n;++k){
    *ptr=std::complex<T>(data[count],data[count+1]);
    count+=2; ++ptr;}
}

template <typename T>
void IOHandler::read_basic_vector(std::vector<T>& input)
{
 int n;
 read_basic<int>(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "vector too large...bad file location?");
 input.resize(n);
 read_common((char*)&input[0],sizeof(T),size_t(n));
}

template <typename T>
void IOHandler::read_complex_vector(std::vector<std::complex<T> >& input)
{
 int n;
 read_basic<int>(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "vector too large...bad file location?");
 input.resize(n);
 std::vector<T> data(2*n);
 read_common((char*)&data[0],sizeof(T),size_t(2*n));
 int count=0; for (int k=0;k<n;++k){
    input[k]=std::complex<T>(data[count],data[count+1]);
    count+=2;}
}

template <typename T>
void IOHandler::read_array(Array<T>& input)
{
 unsigned int n=input.numDimensions();
 read_basic<unsigned int>(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure(n>128, "Array number of dimensions too large...bad file location?");
 read(input.m_sizes);
 read(input.m_store);
}


// ************************************************************  


template <typename T>
size_t IOHandler::numbytes_basic(const T& data) const
{ 
 return sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes_complex(const std::complex<T>& data) const
{ 
 return 2*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes_basic_vector(const std::vector<T>& data) const
{
 return sizeof(int)+data.size()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes_complex_vector(const std::vector<std::complex<T> >& data) const
{
 return sizeof(int)+2*data.size()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes_array(const Array<T>& data) const
{
 return sizeof(unsigned int)+numbytes(data.m_sizes)+numbytes(data.m_store);
}


// **************************************************************

inline void write(IOHandler& ioh, const std::string& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const char& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const int& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const unsigned int& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const long int& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const unsigned long int& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const long long int& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const unsigned long long int& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const float& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const double& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const bool& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const fcmplx& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const dcmplx& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<char>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<unsigned int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<long int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<unsigned long int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<float>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<double>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<fcmplx>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const std::vector<dcmplx>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<char>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<unsigned int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<long int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<unsigned long int>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<float>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<double>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<fcmplx>& output)
 { ioh.write(output); }

inline void write(IOHandler& ioh, const Array<dcmplx>& output)
 { ioh.write(output); }




inline void read(IOHandler& ioh, std::string& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, char& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, int& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, unsigned int& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, long int& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, unsigned long int& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, long long int& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, unsigned long long int& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, float& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, double& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, bool& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, fcmplx& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, dcmplx& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<char>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<unsigned int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<long int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<unsigned long int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<float>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<double>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<fcmplx>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, std::vector<dcmplx>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<char>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<unsigned int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<long int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<unsigned long int>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<float>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<double>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<fcmplx>& input)
 { ioh.read(input); }

inline void read(IOHandler& ioh, Array<dcmplx>& input)
 { ioh.read(input); }


 


inline size_t numbytes(IOHandler& ioh, const std::string& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const char& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const int& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const unsigned int& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const long int& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const unsigned long int& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const long long int& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const unsigned long long int& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const float& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const double& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const bool& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const fcmplx& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const dcmplx& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<char>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<unsigned int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<long int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<unsigned long int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<float>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<double>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<fcmplx>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const std::vector<dcmplx>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<char>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<unsigned int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<long int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<unsigned long int>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<float>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<double>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<fcmplx>& data)
 { return ioh.numbytes(data); }

inline size_t numbytes(IOHandler& ioh, const Array<dcmplx>& data)
 { return ioh.numbytes(data); }


// **************************************************************
}
#endif
