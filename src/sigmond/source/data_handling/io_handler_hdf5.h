#ifndef IO_HANDLER_HDF5_H
#define IO_HANDLER_HDF5_H

#include <vector>
#include <set>
#ifdef XML
#include <unistd.h> 
#endif
#ifdef HDF5
#include <hdf5.h>
#endif
#include <stdexcept>
#include <complex>
#include <iostream>
#include "byte_handler.h"
#include "array.h"

typedef std::complex<double> dcmplx;
typedef std::complex<float>  fcmplx;


 // *********************************************************************************
 // *                                                                               *
 // *       class IOHDF5Handler:     random access input/output                     *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *   This class provides binary input/output using HDF5. It allows the user to   *
 // *   choose between big-endian and little-endian format, and it allows the user  *
 // *   to turn check sums on or off.  Files written by this class will contain     *
 // *   endian-ness information and an ID string. Files read by this class expect   *
 // *   this information.  An error occurs if the ID string given in the open       *
 // *   command on an existing file does not match that in the file.                *
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
 // *   Unlike fstreams, there are no seek functions. Instead, the group or         *
 // *   directory structure of an HDF5 file is used.  Navigation is done very       *
 // *   similarly to linux files: cd to change directory, mkdir to make a new       *
 // *   directory.  The behavior of cd and mkdir differ from the linux commands.    *
 // *   "mkdir" will create all directories (groups) in the request path, and       *
 // *   "cd" has an option to create the directory if it does not already exist.    *
 // *                                                                               *
 // *   Input is done with read(..) commands, and output is done with write(..)     *
 // *   commands.  The data types currently supported for read/write are            *
 // *                                                                               *
 // *       - basic data types (int, bool, float, string, etc.)                     *
 // *       - complex numbers using the STL complex class                           *
 // *       - vector, Array of basic/complex data types                             *
 // *                                                                               *
 // *   If check sums are turned on, objects of this class maintain a check sum.    *
 // *   Changing from a read to a write or vice versa resets the checksum; or an    *
 // *   explicit reset can also be done.  Navigating using "cd" resets the checksum.*
 // *   The checksum is updated as bytes are read or written. Successive reads or   *
 // *   successive writes update the checksum.                                      *
 // *                                                                               *
 // *   WARNING: sigmond multidimensional Array is column major (indices on left    *
 // *   are fastest varying).  HDF5 stores multidimensional arrays as row major,    *
 // *   which cannot be changed!  Objects of this class deal with this easily,      *
 // *   but if you use separate code to access files produced by this class,        *
 // *   the separate code must be made aware of this!!!!                            *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************


class IOHDF5Handler
{

   hid_t fid;          // the HDF5 file id
   std::string m_filename;
   bool is_new_file;

   bool read_only;
   bool openflag;
   bool read_mode;
   
   char endian_format;            // 'B' for big-endian, 'L' for little-endian
   bool endian_convert;

   bool checksum_on;
   hid_t currid;      // HDF5 identifier for current working group (directory)
   std::vector<std::string> currwd;
 
   ByteHandler::n_uint32_t checksum;

      // disallow copying
   IOHDF5Handler(const IOHDF5Handler&);
   IOHDF5Handler(IOHDF5Handler&);
   IOHDF5Handler& operator=(const IOHDF5Handler&);
   IOHDF5Handler& operator=(IOHDF5Handler&);


 public:
 
   enum OpenMode { ReadOnly, ReadWriteFailIfExists, ReadWriteEraseIfExists, 
                   ReadWriteUpdateIfExists };

   explicit IOHDF5Handler();

   IOHDF5Handler(const std::string& filename, OpenMode mode=ReadOnly,
                 const std::string& filetype_id="", char endianness='N',
                 bool turn_on_checksum=false);

   IOHDF5Handler(const std::string& filename, std::ios_base::openmode mode,
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

   ~IOHDF5Handler();

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

   typedef std::iostream::openmode   openmode_type;

         //  Check summing
   
   void turnOnChecksum();
   void turnOffChecksum();
   void resetChecksum();   
   ByteHandler::n_uint32_t getChecksum();
   void printFileID();

     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened, its file type matches 
     // "filetype_id", and a string value is found; returns false otherwise.  
     // The string is returned in "stringvalue".  This routine is used by 
     // objects in data_io_handler.h when the multi-file handlers build 
     // up their maps of file keys.
           
   bool peekString(std::string& stringvalue, const std::string& path,
                   const std::string& filename,
                   const std::string& filetype_id="") const;

     // Open file, read ID string at beginning of file, then close. Returns
     // true if file can be opened and read; false other.
     // Returns ID string in "stringvalue".

   bool peekID(std::string& stringvalue, const std::string& filename) const;

         // navigational routines

   void mkdir(const std::string& path);
   void cd(const std::string& path, bool makedir=false);
   std::string pwd() const { return form_absolute_path(currwd); }
   bool queryDir(const std::string& path) const;
   bool queryData(const std::string& objname) const;
   std::set<std::string> getAllDataNames() const;
   std::set<std::string> getDataNamesInCurrentDir() const;
   std::set<std::string> getDirNamesInCurrentDir() const;
   std::set<std::string> getAllDirNames() const;

         // write routines  

   void write(const std::string& objname, const std::string& output);
   void write(const std::string& objname, const char& output);
   void write(const std::string& objname, const int& output);
   void write(const std::string& objname, const unsigned int& output);
   void write(const std::string& objname, const long int& output);
   void write(const std::string& objname, const unsigned long int& output);
   void write(const std::string& objname, const long long int& output);
   void write(const std::string& objname, const unsigned long long int& output);
   void write(const std::string& objname, const float& output);
   void write(const std::string& objname, const double& output);
   void write(const std::string& objname, const fcmplx& output);
   void write(const std::string& objname, const dcmplx& output);
   void write(const std::string& objname, const std::vector<char>& output);
   void write(const std::string& objname, const std::vector<int>& output);
   void write(const std::string& objname, const std::vector<unsigned int>& output);
   void write(const std::string& objname, const std::vector<long int>& output);
   void write(const std::string& objname, const std::vector<unsigned long int>& output);
   void write(const std::string& objname, const std::vector<float>& output);
   void write(const std::string& objname, const std::vector<double>& output);
   void write(const std::string& objname, const std::vector<fcmplx>& output);
   void write(const std::string& objname, const std::vector<dcmplx>& output);
   void write(const std::string& objname, const Array<int>& output);
   void write(const std::string& objname, const Array<unsigned int>& output);
   void write(const std::string& objname, const Array<long int>& output);
   void write(const std::string& objname, const Array<unsigned long int>& output);
   void write(const std::string& objname, const Array<float>& output);
   void write(const std::string& objname, const Array<double>& output);
   void write(const std::string& objname, const Array<fcmplx>& output);
   void write(const std::string& objname, const Array<dcmplx>& output);

          // read routines

   void read(const std::string& objname, std::string& input);
   void read(const std::string& objname, char& input);
   void read(const std::string& objname, int& input);
   void read(const std::string& objname, unsigned int& input);
   void read(const std::string& objname, long int& input);
   void read(const std::string& objname, unsigned long int& input);
   void read(const std::string& objname, long long int& input);
   void read(const std::string& objname, unsigned long long int& input);
   void read(const std::string& objname, float& input);
   void read(const std::string& objname, double& input);
   void read(const std::string& objname, fcmplx& input);
   void read(const std::string& objname, dcmplx& input);
   void read(const std::string& objname, std::vector<char>& input);
   void read(const std::string& objname, std::vector<int>& input);
   void read(const std::string& objname, std::vector<unsigned int>& input);
   void read(const std::string& objname, std::vector<long int>& input);
   void read(const std::string& objname, std::vector<unsigned long int>& input);
   void read(const std::string& objname, std::vector<float>& input);
   void read(const std::string& objname, std::vector<double>& input);
   void read(const std::string& objname, std::vector<fcmplx>& input);
   void read(const std::string& objname, std::vector<dcmplx>& input);
   void read(const std::string& objname, Array<int>& input);
   void read(const std::string& objname, Array<unsigned int>& input);
   void read(const std::string& objname, Array<long int>& input);
   void read(const std::string& objname, Array<unsigned long int>& input);
   void read(const std::string& objname, Array<float>& input); 
   void read(const std::string& objname, Array<double>& input);
   void read(const std::string& objname, Array<fcmplx>& input);
   void read(const std::string& objname, Array<dcmplx>& input);

 private:

      // private utility routines

   void open_existing_file(const std::string& filetype_id);
   void open_new_file(const std::string& filetype_id, char endianness);
   void delete_file();
   void file_close();
   void clear();

   bool fileExists();
   void check_for_failure(int errcode, const std::string& mesg="") const;
   void check_for_hid_failure(hid_t id, const std::string& mesg="") const;
   void check_for_hid_failure2(hid_t id, const std::string& mesg1="",const std::string& mesg2="") const;
   void check_for_herr_failure(herr_t status, const std::string& mesg="") const;
   std::string tidyString(const std::string& str, const std::string& leadtrailchars=" ") const;  
   std::vector<std::string> tokenize(const std::string& str, const std::string& delimiter = "/", 
                                     bool allowdot=true, bool allowdotdot=false) const;
   void check_token(const std::string& token, bool allowdot=true, bool allowdotdot=false) const;
   void check_path(const std::string& str, bool allowdot=true, bool allowdotdot=false) const;
   std::string form_absolute_path(const std::vector<std::string>& path_links) const;
   void readIDstring(char &endianness, std::string& ID_string);
   void writeIDstring(const std::string& ID_string);
   char query_obj(const std::string& objname) const;
   void collect_data_names(hid_t loc_id, const std::string& path, const char* name,
                           std::set<std::string>& collected_names) const;
   void collect_dir_names(hid_t loc_id, const std::string& path, const char* name,
                          std::set<std::string>& collected_dir_names) const;
   void collect_data_names_in_currdir(hid_t loc_id,
                           std::set<std::string>& collected_names, H5O_type_t query_type) const;
   void peek_strings(const std::string& filename, const std::vector<std::string>& paths,
                     std::vector<std::string>& stringvalues) const; 

       // constants
       
   static const int IO_ERR_NO_SUCH_FILE;
   static const int IO_ERR_ACCESS;
   static const int IO_ERR_OTHER;
   static const int IO_SUCCESS;

   static const openmode_type IO_MODE_RDONLY;
   static const openmode_type IO_MODE_RDWR;
   static const openmode_type IO_MODE_CREATE;
   static const openmode_type IO_MODE_EXCL;

   static ByteHandler m_bytehandler;

   hid_t dtype_int;
   hid_t dtype_long;
   hid_t dtype_llong;
   hid_t dtype_uint;
   hid_t dtype_ulong;
   hid_t dtype_ullong;
   hid_t dtype_float;
   hid_t dtype_double;

   void assign_dtypes();

          // private write members

   template <typename T>
   void write_atomics(const std::string& objname, const T* output, size_t nT,
                      hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void write_complex_values(const std::string& objname, 
                             const std::complex<T>* output, size_t n,
                             hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void write_array(const std::string& objname, const Array<T>& output,
                    hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void write_complex_array(const std::string& objname, const Array<T>& output,
                            hid_t dtype_id, hid_t mem_type_id);

          // private read members

   template <typename T>
   void read_atomic(const std::string& objname, T& input,
                    hid_t dtype_id, hid_t mem_type_id);
   template <typename T>
   void read_atomics(const std::string& objname, std::vector<T>& input,
                     hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void read_complex(const std::string& objname, std::complex<T>& output,
                     hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void read_complex_vector(const std::string& objname, 
                            std::vector<std::complex<T> >& output,
                            hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void read_array(const std::string& objname, Array<T>& output,
                   hid_t dtype_id, hid_t mem_type_id);

   template <typename T>
   void read_complex_array(const std::string& objname, Array<T>& output,
                           hid_t dtype_id, hid_t mem_type_id);
};


 // ***************************************************************

            // private write members

template <typename T>
void IOHDF5Handler::write_atomics(const std::string& objname, const T* output, size_t nT,
                                  hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write when no open file");}
 if (read_only){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write to read-only file");}
 check_path(objname);
 hsize_t dsize=nT;
 hid_t dataspace_id = H5Screate_simple(1, &dsize, NULL);
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dcreate2(*cwdptr, objname.c_str(), dtype_id, dataspace_id, 
                               H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 check_for_hid_failure(dataset_id,"Could not create object name for writing");
 herr_t status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, output);
 check_for_herr_failure(status,"Could not write data");
 status = H5Sclose(dataspace_id);
 status = H5Dclose(dataset_id);
 if (read_mode){
    read_mode=false; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
             reinterpret_cast<const char *>(output), nT*sizeof(T));
}


template <typename T>
void IOHDF5Handler::write_complex_values(const std::string& objname, 
                                         const std::complex<T>* output, size_t n,
                                         hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write when no open file");}
 if (read_only){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write to read-only file");}
 std::vector<T> buffer(2*n);
 T* dptr=buffer.data();
 const std::complex<T> *ptr=output; 
 for (size_t k=0;k<n;++k){
    *dptr=real(*ptr); ++dptr;
    *dptr=imag(*ptr); ++dptr;
    ++ptr;}
 write_atomics(objname,buffer.data(),2*n,dtype_id,mem_type_id);
 if (read_mode){
    read_mode=false; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
              reinterpret_cast<const char *>(output), 2*n*sizeof(T));
}


template <typename T>
void IOHDF5Handler::write_array(const std::string& objname, const Array<T>& output,
                                hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write when no open file");}
 if (read_only){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write to read-only file");}
 check_path(objname);
 uint rank = output.numDimensions();
 hsize_t* dims=new hsize_t[rank];
 for (uint k=0;k<rank;++k) dims[k]=output.size(k);
 hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dcreate2(*cwdptr, objname.c_str(), dtype_id, dataspace_id, 
                               H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 check_for_hid_failure(dataset_id,"Could not create object name for writing");
 herr_t status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, output.m_store.data());
 check_for_herr_failure(status,"Could not write data");
 status = H5Sclose(dataspace_id);
 status = H5Dclose(dataset_id);
 delete [] dims;
 if (read_mode){
    read_mode=false; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
           reinterpret_cast<const char *>(output.m_store.data()), output.size()*sizeof(T));
}


template <typename T>
void IOHDF5Handler::write_complex_array(const std::string& objname, const Array<T>& output,
                                        hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write when no open file");}
 if (read_only){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write to read-only file");}
 check_path(objname);
 uint rank = output.numDimensions()+1;
 hsize_t* dims=new hsize_t[rank];
 dims[0]=2; // real and imaginary parts
 for (uint k=1;k<rank;++k) dims[k]=output.size(k-1);
 hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dcreate2(*cwdptr, objname.c_str(), dtype_id, dataspace_id, 
                               H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 check_for_hid_failure(dataset_id,"Could not create object name for writing");
 herr_t status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, output.m_store.data());
 check_for_herr_failure(status,"Could not write data");
 status = H5Sclose(dataspace_id);
 status = H5Dclose(dataset_id);
 delete [] dims;
 if (read_mode){
    read_mode=false; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
               reinterpret_cast<const char *>(output.m_store.data()), 2*output.size()*sizeof(T));
 }

// ************************************************************  
 
            // private read members

template <typename T>
void IOHDF5Handler::read_atomic(const std::string& objname, T& input,
                                hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dopen2(*cwdptr,objname.c_str(),H5P_DEFAULT);
 check_for_hid_failure2(dataset_id,"Could not find object name",objname);
 hid_t dtype = H5Dget_type(dataset_id);
 if (!H5Tequal(dtype,dtype_id)) check_for_failure(IO_ERR_OTHER,"Datatype mismatch during read");
 hid_t dataspace_id = H5Dget_space(dataset_id);
 herr_t status = H5Dread(dataset_id,mem_type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,&input);
 herr_t status1 = H5Sclose(dataspace_id);
 status1 = H5Dclose(dataset_id);                 
 check_for_herr_failure(status,"Could not read data"); 
 check_for_herr_failure(status1,"Could not read data"); 
 H5Tclose(dtype);
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
            reinterpret_cast<const char *>(&input), sizeof(T));
}


template <typename T>
void IOHDF5Handler::read_atomics(const std::string& objname, std::vector<T>& input,
                                 hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dopen2(*cwdptr,objname.c_str(),H5P_DEFAULT);
 check_for_hid_failure2(dataset_id,"Could not find object name",objname);
 hid_t dtype = H5Dget_type(dataset_id);
 if (!H5Tequal(dtype,dtype_id)) check_for_failure(IO_ERR_OTHER,"Datatype mismatch during read");
 hsize_t dsize;
 hid_t dataspace_id = H5Dget_space(dataset_id);
 uint rank=H5Sget_simple_extent_ndims(dataspace_id);
 if (rank!=1) check_for_failure(IO_ERR_OTHER,"Datatype rank mismatch during read");
 herr_t status = H5Sget_simple_extent_dims(dataspace_id,&dsize,NULL);
 input.resize(dsize);
 status = H5Dread(dataset_id,mem_type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,input.data());
 herr_t status1 = H5Sclose(dataspace_id);
 status1 = H5Dclose(dataset_id);                 
 check_for_herr_failure(status,"Could not read data"); 
 check_for_herr_failure(status1,"Could not read data"); 
 H5Tclose(dtype);
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
              reinterpret_cast<const char *>(input.data()), input.size()*sizeof(T));
}


template <typename T>
void IOHDF5Handler::read_complex(const std::string& objname, std::complex<T>& input,
                                 hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::vector<T> buffer;
 read_atomics(objname,buffer,dtype_id,mem_type_id);
 if (buffer.size()!=2){
    check_for_failure(IO_ERR_OTHER,"Failure reading a complex value");}
 input=std::complex<T>(buffer[0],buffer[1]);
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
             reinterpret_cast<const char *>(&input), 2*sizeof(T));
}


template <typename T>
void IOHDF5Handler::read_complex_vector(const std::string& objname, 
                                        std::vector<std::complex<T> >& input,
                                        hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::vector<T> buffer;
 read_atomics(objname,buffer,dtype_id,mem_type_id);
 size_t n=buffer.size();
 if ((n%2)!=0){
    check_for_failure(IO_ERR_OTHER,"Failure reading a complex vector");}
 n/=2;
 input.resize(n);
 std::complex<T> *ptr=input.data();
 const T *bptr=buffer.data();
 for (size_t k=0;k<n;++k){
    *ptr=std::complex<T>(*bptr,*(bptr+1));
    bptr+=2; ++ptr;}
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
              reinterpret_cast<const char *>(input.data()), input.size()*sizeof(T));
}


template <typename T>
void IOHDF5Handler::read_array(const std::string& objname, Array<T>& input,
                               hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dopen2(*cwdptr,objname.c_str(),H5P_DEFAULT);
 check_for_hid_failure2(dataset_id,"Could not find object name",objname);
 hid_t dtype = H5Dget_type(dataset_id);
 if (!H5Tequal(dtype,dtype_id)) check_for_failure(IO_ERR_OTHER,"Datatype mismatch during read");
 hid_t dataspace_id = H5Dget_space(dataset_id);
 int rank=H5Sget_simple_extent_ndims(dataspace_id);
 if (rank<=0) check_for_failure(IO_ERR_OTHER,"Datatype rank mismatch during read");
 hsize_t* dims=new hsize_t[rank];
 herr_t status = H5Sget_simple_extent_dims(dataspace_id,dims,NULL);
 std::vector<uint> sizevals(rank);
 for (int k=0;k<rank;++k) sizevals[k]=dims[k];
 delete [] dims;
 input.resize(sizevals);
 status = H5Dread(dataset_id,mem_type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,input.m_store.data());
 herr_t status1 = H5Sclose(dataspace_id);
 status1 = H5Dclose(dataset_id);                 
 check_for_herr_failure(status,"Could not read data"); 
 check_for_herr_failure(status1,"Could not read data"); 
 H5Tclose(dtype);
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
              reinterpret_cast<const char *>(input.m_store.data()), input.size()*sizeof(T));
}


template <typename T>
void IOHDF5Handler::read_complex_array(const std::string& objname, Array<T>& input,
                                       hid_t dtype_id, hid_t mem_type_id)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dopen2(*cwdptr,objname.c_str(),H5P_DEFAULT);
 check_for_hid_failure2(dataset_id,"Could not find object name",objname);
 hid_t dtype = H5Dget_type(dataset_id);
 if (!H5Tequal(dtype,dtype_id)) check_for_failure(IO_ERR_OTHER,"Datatype mismatch during read");
 hid_t dataspace_id = H5Dget_space(dataset_id);
 int rank=H5Sget_simple_extent_ndims(dataspace_id);
 if (rank<=0) check_for_failure(IO_ERR_OTHER,"Datatype rank mismatch during read");
 hsize_t* dims=new hsize_t[rank];
 herr_t status = H5Sget_simple_extent_dims(dataspace_id,dims,NULL);
 std::vector<uint> sizevals(rank-1);
 for (int k=1;k<rank;++k) sizevals[k-1]=dims[k];
 delete [] dims;
 input.resize(sizevals);
 status = H5Dread(dataset_id,mem_type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,input.m_store.data());
 herr_t status1 = H5Sclose(dataspace_id);
 status1 = H5Dclose(dataset_id);                 
 check_for_herr_failure(status,"Could not read data"); 
 check_for_herr_failure(status1,"Could not read data"); 
 H5Tclose(dtype);
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, 
               reinterpret_cast<const char *>(input.m_store.data()), 2*input.size()*sizeof(T));
}

// **************************************************************
#endif
