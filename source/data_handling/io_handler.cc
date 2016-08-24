#include "io_handler.h"
#include <unistd.h>
using namespace std;


const int IOHandler::ID_string_length = 32;
const IOHandler::pos_type IOHandler::data_start_pos = 1+IOHandler::ID_string_length;


// *************************************************************************

            //  serial code


void IOHandler::delete_file()
{
 check_for_failure(remove(m_filename.c_str()),"Failure deleting file");
}

void IOHandler::file_open(IOHandler::openmode_type access_mode, 
                          const std::string& errmsg)
{
 fh.open(m_filename.c_str(), std::ios::binary | access_mode);
 check_for_failure(!fh, errmsg);
}
   
void IOHandler::file_close()
{
 fh.close();
 check_for_failure(!fh, "Failure during close");
}

const int IOHandler::IO_ERR_NO_SUCH_FILE= 1;
const int IOHandler::IO_ERR_ACCESS=       1;
const int IOHandler::IO_ERR_OTHER=        1;
const int IOHandler::IO_SUCCESS=          0;

const IOHandler::openmode_type IOHandler::IO_MODE_RDONLY= ios::in;
const IOHandler::openmode_type IOHandler::IO_MODE_RDWR=   ios::in | ios::out;
const IOHandler::openmode_type IOHandler::IO_MODE_CREATE= ios::trunc;
const IOHandler::openmode_type IOHandler::IO_MODE_EXCL=   ios::trunc;

const IOHandler::whence_type IOHandler::IO_SEEK_BEG  = ios::beg;
const IOHandler::whence_type IOHandler::IO_SEEK_CUR  = ios::cur;
const IOHandler::whence_type IOHandler::IO_SEEK_END  = ios::end;


void IOHandler::check_for_failure(int errcode, const std::string& mesg)
{
 if (errcode==IOHandler::IO_SUCCESS) return;
 std::cout << "IOHandler error with file "<<m_filename<<":"
             <<endl<<"  "<<mesg<<endl;
 exit(errcode);
}


void IOHandler::file_seek(off_type offset, whence_type whence)
{
 fh.seekg(offset, whence);
 check_for_failure(fh.fail(),"Failure during seek");
}

IOHandler::pos_type IOHandler::file_tell()
{
 pos_type current=fh.tellg();
 if (current<0) check_for_failure(1,"Failure during tell");
 return current;
}


void IOHandler::writeCommon(const char *data, size_t nbytes)
{
 pos_type currdisp=fh.tellg();
 fh.seekp(currdisp);  // just to be sure
 fh.write(data,nbytes);
 check_for_failure(fh.fail(),"Failure during common write");
 fh.seekg(currdisp+pos_type(nbytes));
}


void IOHandler::readCommon(char *data, size_t nbytes)
{
 pos_type currdisp=fh.tellg();
 fh.read(data,nbytes);
 check_for_failure(fh.fail(),"Failure during common read");
 fh.seekg(currdisp+pos_type(nbytes)); 
}


// *************************************************************************





IOHandler::IOHandler() : read_only(true), openflag(false), read_mode(true),
                         endian_format('U'), endian_convert(false), 
                         checksum_on(false), is_new_file(false), checksum(0)
{
#ifndef NO_CXX11
 static_assert(sizeof(int)==4,"Invalid int size");
#else
 static__assert<sizeof(int)==4>();
#endif
}


IOHandler::IOHandler(const std::string& filename, OpenMode mode,
                     const std::string& filetype_id, char endianness,
                     bool turn_on_checksum)
{
#ifndef NO_CXX11
 static_assert(sizeof(int)==4,"Invalid int size");
#else
 static__assert<sizeof(int)==4>();
#endif
 openflag=false;
 open(filename,mode,filetype_id,endianness,turn_on_checksum);
}

IOHandler::IOHandler(const std::string& filename, ios_base::openmode mode,
                     const std::string& filetype_id, char endianness,
                     bool turn_on_checksum)
{
#ifndef NO_CXX11
 static_assert(sizeof(int)==4,"Invalid int size");
#else
 static__assert<sizeof(int)==4>();
#endif
 openflag=false;
 open(filename,mode,filetype_id,endianness,turn_on_checksum);
}


void IOHandler::open(const std::string& filename, IOHandler::OpenMode mode,
                     const std::string& filetype_id, char endianness,
                     bool turn_on_checksum)
{
 close();
 read_only=(mode==ReadOnly)?true:false;

 m_filename=tidyString(filename);
 if (m_filename.empty())
     check_for_failure(IO_ERR_NO_SUCH_FILE,"Empty file name");

 bool exists=fileExists();
 if ((mode==ReadOnly)&&(!exists))
    check_for_failure(IO_ERR_NO_SUCH_FILE,
           "Failure during ReadOnly open: file does not exist");
 else if ((mode==ReadWriteFailIfExists)&&(exists))
    check_for_failure(IO_ERR_ACCESS,
           "Failure during ReadWrite open: file exists and FailIfExists mode");
 else if ((mode==ReadWriteEraseIfExists)&&(exists)){
    delete_file();
    exists=false;}
 if (exists){
    IOHandler::openmode_type access=(mode==ReadOnly) ? 
                IO_MODE_RDONLY : IO_MODE_RDWR;
    is_new_file=false;
    open_existing_file(filetype_id,access);}
 else{
    is_new_file=true;
    open_new_file(filetype_id,endianness);}

 checksum_on=turn_on_checksum;
 checksum=0;
 read_mode=true;
}


void IOHandler::open(const std::string& filename, ios_base::openmode mode,
                     const std::string& filetype_id, char endianness,
                     bool turn_on_checksum)
{
 IOHandler::OpenMode iomode;
 if (!(mode & ios_base::out)) iomode=ReadOnly;
 else{
    if (mode & ios_base::trunc) iomode=ReadWriteEraseIfExists;
    else iomode=ReadWriteFailIfExists;}
 open(filename,iomode,filetype_id,endianness,turn_on_checksum);
}


void IOHandler::openReadOnly(const std::string& filename, 
                             const std::string& filetype_id,
                             bool turn_on_checksum)
{
 open(filename,ReadOnly,filetype_id,'N',turn_on_checksum);
}


void IOHandler::openNew(const std::string& filename, bool fail_if_exists,
                        const std::string& filetype_id, char endianness,
                        bool turn_on_checksum)
{
 OpenMode iomode=(fail_if_exists) ? ReadWriteFailIfExists : ReadWriteEraseIfExists;
 open(filename,iomode,filetype_id,endianness,turn_on_checksum);
}


void IOHandler::openUpdate(const std::string& filename,
                           const std::string& filetype_id, char endianness,
                           bool turn_on_checksum)
{
 open(filename,ReadWriteUpdateIfExists,filetype_id,endianness,turn_on_checksum);
}



void IOHandler::open_existing_file(const std::string& filetype_id,
                                   IOHandler::openmode_type access_mode)
{
 file_open(access_mode,"Could not open existing file: "+m_filename);
     // get endian info, and check file id
 openflag=true;
 string ID_string;
 readIDstring(endian_format,ID_string); 
 bool flag;
 flag=true;
 if ((endian_format!='B')&&(endian_format!='L')){
    flag=false; std::cerr << "Invalid endian format"<<std::endl;}
 if (tidyString(filetype_id)!=ID_string){
    flag=false; 
    std::cerr << "File = "<<m_filename<<std::endl;
    std::cerr << "File ID mismatch:"<<std::endl;
    std::cerr << "File contains ID: <"<<ID_string<<">"<<std::endl;
    std::cerr << "ID requested was: <"<<tidyString(filetype_id)<<">"<<std::endl;}
 if (!flag)
    check_for_failure(IO_ERR_OTHER,"Error during open");
 if (QDPUtil::big_endian())
    endian_convert=(endian_format=='L')?true:false;
 else 
    endian_convert=(endian_format=='B')?true:false;
}


void IOHandler::open_new_file(const std::string& filetype_id,
                              char endianness)
{
 if (endianness=='N'){  // native 
    endian_format=QDPUtil::big_endian() ? 'B':'L';
    endian_convert=false;}
 else if ((endianness!='B')&&(endianness!='L'))
    check_for_failure(IO_ERR_OTHER,"Invalid endian format");
 else{
    endian_format=endianness;
    if (QDPUtil::big_endian()) 
       endian_convert=(endian_format=='L')?true:false;
    else 
       endian_convert=(endian_format=='B')?true:false;}

 IOHandler::openmode_type amode = IO_MODE_RDWR | IO_MODE_CREATE;
 
 file_open(amode,"Could not open file "+m_filename);

 openflag=true;
 writeIDstring(filetype_id);
}


       // Destructor

IOHandler::~IOHandler() 
{
 if (openflag) clear(); 
}


     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened and its file type matches 
     // "filetype_id", returns false otherwise.  The string is returned in 
     // "stringvalue".  This routine is used by objects in data_io_handler.h 
     // when the multi-file handlers build up their maps of file keys.
           
bool IOHandler::peekString(std::string& stringvalue, unsigned int byte_offset,
                           const std::string& filename, const std::string& filetype_id)
{
 stringvalue.clear();
 string fname=tidyString(filename);
 if (fname.empty()) return false;
 bool flag;
 flag=peeker(stringvalue,byte_offset,fname,filetype_id);
 return flag;
}



bool IOHandler::peeker(std::string& stringvalue, unsigned int byte_offset,
                       const std::string& fname, const std::string& filetype_id)
{
 ifstream in(fname.c_str(), std::ios::binary | ios::in);
 if (!in) return false;
 string ID_string(ID_string_length+1,' ');
 in.read((char*)&ID_string[0],ID_string_length+1);
 if (in.fail()) return false;
 char endian=ID_string[0];
 if ((endian!='B')&&(endian!='L')) return false;
 ID_string.erase(0,1);
 ID_string=tidyString(ID_string);
 if (tidyString(filetype_id)!=ID_string) return false;
 in.seekg(pos_type(byte_offset),ios::cur);
 if (in.fail()) return false;
 unsigned int n;
 in.read((char*)&n,sizeof(int));
 if (in.fail()) return false;
 if (QDPUtil::big_endian())
    endian_convert=(endian=='L')?true:false;
 else 
    endian_convert=(endian=='B')?true:false;
 if (endian_convert) QDPUtil::byte_swap(&n,sizeof(int),1);
 if (n>16777216) return false;  // too large for string...must be corrupt data
 stringvalue.resize(n);
 in.read((char*)&stringvalue[0],sizeof(char)*n);
 if (in.fail()){ stringvalue.clear(); return false;}
 return true;       // in destructor will close the file
}

     // return true if file can be opened and read; false other.
     // Return ID string in "stringvalue".

bool IOHandler::peekID(std::string& stringvalue, const std::string& filename)
{
 stringvalue.clear();
 string fname=tidyString(filename);
 if (fname.empty()) return false;
 ifstream in(fname.c_str(), std::ios::binary | ios::in);
 if (!in) return false;
 string ID_string(ID_string_length+1,' ');
 in.read((char*)&ID_string[0],ID_string_length+1);
 if (in.fail()) return false;
 char endian=ID_string[0];
 if ((endian!='B')&&(endian!='L')) return false;
 ID_string.erase(0,1);
 stringvalue=tidyString(ID_string);
 return true;       // in destructor will close the file
}



      // Clear, then reset finfo

void IOHandler::close()
{ 
 if (!openflag) return;
 clear();
}

     // Close the file, reset all data members except finfo

void IOHandler::clear()
{ 
 file_close();
 m_filename.clear();
 openflag=false;
 endian_convert=false;
 endian_format='U';  // undefined
 checksum=0;
 read_only=true;
 read_mode=true;
 checksum_on=false;
 is_new_file=false;
}

         //  Set the file pointer relative to start of file,
         //  end of file (positive is backward), or current location

void IOHandler::seekFromStart(off_type offset)
{
 file_seek(offset+data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOHandler::seekFromCurr(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_CUR);
 checksum=0;
}
 
void IOHandler::seekFromEnd(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_END);
 checksum=0;
}

void IOHandler::seek(pos_type offset)
{ 
 file_seek(offset+data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOHandler::seekBegin(off_type offset)
{ 
 file_seek(offset+data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOHandler::seekRelative(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_CUR);
 checksum=0;
}

void IOHandler::seekEnd(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_END);
 checksum=0;
}

void IOHandler::rewind()
{ 
 file_seek(data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}


      //  Get the file pointer location in bytes from start of data

IOHandler::pos_type IOHandler::currentPosition()
{
 return file_tell()-data_start_pos;
}

IOHandler::pos_type IOHandler::tell()
{
 return file_tell()-data_start_pos;
}


void IOHandler::turnOnChecksum()
{
 if (checksum_on) return;
 checksum_on=true;
 checksum=0;
}

void IOHandler::turnOffChecksum()
{
 checksum_on=false;
}

void IOHandler::resetChecksum()
{
 checksum=0;
}

QDPUtil::n_uint32_t IOHandler::getChecksum()
{
 check_for_failure(!checksum_on,"Invalid call to getChecksum since checksums not turned on");
 return checksum;
}



        // Print out file ID. Does not 
        // change file pointers.

void IOHandler::printFileID()
{
 char endian;
 string ID_string;
 pos_type curr=file_tell();
 readIDstring(endian,ID_string);
 file_seek(curr,IOHandler::IO_SEEK_BEG);
 std::cout << "File = "<<m_filename<<std::endl<<"ID string = <"
           <<tidyString(ID_string)<<">"<<std::endl;
}

   // Write ID string at start of file. Length of string
   // is always "ID_string_length".  First character will
   // be 'B' or 'L' to indicate endian format to use.

void IOHandler::writeIDstring(const std::string& ID_string)
{
 check_for_failure(int(ID_string.length())>ID_string_length,
                  "IOHandler file ID string too long: cannot exceed "
                  +int_to_string(ID_string_length)+" characters");
 string buf(1,endian_format);
 buf+=ID_string;
 int nblanks=ID_string_length-ID_string.length();
 if (nblanks>0)
    buf+=string(nblanks,' ');
 file_seek(0,IOHandler::IO_SEEK_BEG);
 writeCommon(&buf[0],ID_string_length+1);
}


   // Reads endian format and ID string from current file.
   // Results are returned in "endianness" and "ID_string".
   // All file pointers updated.

void IOHandler::readIDstring(char &endianness, std::string& ID_string)
{
 ID_string.resize(ID_string_length+1);
 file_seek(0,IOHandler::IO_SEEK_BEG);
 readCommon(&ID_string[0],ID_string_length+1);
 endianness=ID_string[0];
 ID_string.erase(0,1);
 ID_string=tidyString(ID_string);
}




      // Checks if a file exists.
      
bool IOHandler::fileExists()
{
 bool result;
 result = (access(m_filename.c_str(),F_OK) == 0) ? true : false;
 return result;
}


      // Converts an integer to a string
      
string IOHandler::int_to_string(int intval)
{ 
 std::ostringstream oss;
 oss << intval;
 return oss.str();
}

      // Removes leading and trailing blanks in a string

string IOHandler::tidyString(const string& str)   
{
 string tmp;
 for (int i=0;i<int(str.length());i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 int start=tmp.find_first_not_of(" ");
 if (size_t(start)==string::npos) return "";
 int len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}


// *******************************************************************

      //   Main input/output routines that do the byte-swapping (if needed),
      //   update the check sum, and do the read/write.

void IOHandler::write_common(const char* output, size_t elementbytes, size_t nelements)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (read_mode){
    read_mode=false; checksum=0;}
 if (endian_convert){
    QDPUtil::byte_swap(const_cast<char *>(output), elementbytes, nelements);
    if (checksum_on) checksum = QDPUtil::crc32(checksum, output, elementbytes*nelements);
    writeCommon(output, elementbytes*nelements);
    QDPUtil::byte_swap(const_cast<char *>(output), elementbytes, nelements);}
 else{
    if (checksum_on) checksum = QDPUtil::crc32(checksum, output, elementbytes*nelements);
    writeCommon(output, elementbytes*nelements);}
}


void IOHandler::read_common(char* input, size_t elementbytes, size_t nelements)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!read_mode){
    read_mode=true; checksum=0;}
 readCommon(input, elementbytes*nelements);
 if (checksum_on) checksum = QDPUtil::crc32(checksum, input, elementbytes*nelements);
 if (endian_convert) QDPUtil::byte_swap(input, elementbytes, nelements);
}


// *******************************************************************


void IOHandler::write(const std::string& output)
{
 int n=output.length();
 write_common((const char*)&n, sizeof(int), 1);
 write_common(output.data(), sizeof(char), n);
}

void IOHandler::write(const char& output)
{
 write_basic<char>(output);
}

void IOHandler::write(const int& output)
{
 write_basic<int>(output);
}

void IOHandler::write(const unsigned int& output)
{
 write_basic<unsigned int>(output);
}

void IOHandler::write(const long int& output)
{
 write_basic<long int>(output);
}

void IOHandler::write(const unsigned long int& output)
{
 write_basic<unsigned long int>(output);
}

void IOHandler::write(const long long int& output)
{
 write_basic<long long int>(output);
}

void IOHandler::write(const unsigned long long int& output)
{
 write_basic<unsigned long long int>(output);
}

void IOHandler::write(const float& output)
{
 write_basic<float>(output);
}

void IOHandler::write(const double& output)
{
 write_basic<double>(output);
}

void IOHandler::write(const bool& output)
{
 write_basic<bool>(output);
}

void IOHandler::write(const fcmplx& output)
{
 write_complex<float>(output);
}
 
void IOHandler::write(const dcmplx& output)
{
 write_complex<double>(output);
}


void IOHandler::write(const std::vector<char>& output)
{
 write_basic_vector<char>(output);
}

void IOHandler::write(const std::vector<int>& output)
{
 write_basic_vector<int>(output);
}

void IOHandler::write(const std::vector<unsigned int>& output)
{
 write_basic_vector<unsigned int>(output);
}

void IOHandler::write(const std::vector<long int>& output)
{
 write_basic_vector<long int>(output);
}

void IOHandler::write(const std::vector<unsigned long int>& output)
{
 write_basic_vector<unsigned long int>(output);
}

void IOHandler::write(const std::vector<float>& output)
{
 write_basic_vector<float>(output);
}

void IOHandler::write(const std::vector<double>& output)
{
 write_basic_vector<double>(output);
}

void IOHandler::write(const std::vector<fcmplx>& output)
{
 write_complex_vector<float>(output);
}

void IOHandler::write(const std::vector<dcmplx>& output)
{
 write_complex_vector<double>(output);
}

void IOHandler::write(const Array<char>& output)
{
 write_array<char>(output);
}

void IOHandler::write(const Array<int>& output)
{
 write_array<int>(output);
}

void IOHandler::write(const Array<unsigned int>& output)
{
 write_array<unsigned int>(output);
}

void IOHandler::write(const Array<long int>& output)
{
 write_array<long int>(output);
}

void IOHandler::write(const Array<unsigned long int>& output)
{
 write_array<unsigned long int>(output);
}

void IOHandler::write(const Array<float>& output)
{
 write_array<float>(output);
}

void IOHandler::write(const Array<double>& output)
{
 write_array<double>(output);
}

void IOHandler::write(const Array<fcmplx>& output)
{
 write_array<fcmplx>(output);
}

void IOHandler::write(const Array<dcmplx>& output)
{
 write_array<dcmplx>(output);
}


void IOHandler::multi_write(const char* output, int n)
{ 
 multi_write_basic<char>(output,n);
}

void IOHandler::multi_write(const int* output, int n)
{ 
 multi_write_basic<int>(output,n);
}

void IOHandler::multi_write(const unsigned int* output, int n)
{ 
 multi_write_basic<unsigned int>(output,n);
}

void IOHandler::multi_write(const long int* output, int n)
{ 
 multi_write_basic<long int>(output,n);
}

void IOHandler::multi_write(const unsigned long int* output, int n)
{ 
 multi_write_basic<unsigned long int>(output,n);
}

void IOHandler::multi_write(const long long int* output, int n)
{ 
 multi_write_basic<long long int>(output,n);
}

void IOHandler::multi_write(const unsigned long long int* output, int n)
{ 
 multi_write_basic<unsigned long long int>(output,n);
}

void IOHandler::multi_write(const float* output, int n)
{ 
 multi_write_basic<float>(output,n);
}

void IOHandler::multi_write(const double* output, int n)
{ 
 multi_write_basic<double>(output,n);
}

void IOHandler::multi_write(const bool* output, int n)
{ 
 multi_write_basic<bool>(output,n);
}

void IOHandler::multi_write(const fcmplx* output, int n)
{
 multi_write_complex<float>(output,n);
}

void IOHandler::multi_write(const dcmplx* output, int n)
{ 
 multi_write_complex<double>(output,n);
}


// ********************************************************************

void IOHandler::read(std::string& input)
{
 int n;
 read_common((char*)&n, sizeof(int), 1);
    // if reading in wrong location, could get nonsense here,
    // so limit to a 16MB string
 check_for_failure(n>16777216, "string for read too large...bad file location?");
 char* str = new(nothrow) char[n];
 check_for_failure((str==0),"IOHandler::read---unable to allocate memory");
 readCommon(str, sizeof(char)*n);
 if (checksum_on) checksum = QDPUtil::crc32(checksum, str, sizeof(char)*n);
 input.assign(str, n);
 delete[] str;
}

void IOHandler::read(char& input)
{
 read_basic<char>(input);
}

void IOHandler::read(int& input)
{
 read_basic<int>(input);
}

void IOHandler::read(unsigned int& input)
{
 read_basic<unsigned int>(input);
}

void IOHandler::read(long int& input)
{
 read_basic<long int>(input);
}

void IOHandler::read(unsigned long int& input)
{
 read_basic<unsigned long int>(input);
}

void IOHandler::read(long long int& input)
{
 read_basic<long long int>(input);
}

void IOHandler::read(unsigned long long int& input)
{
 read_basic<unsigned long long int>(input);
}

void IOHandler::read(float& input)
{
 read_basic<float>(input);
}

void IOHandler::read(double& input)
{
 read_basic<double>(input);
}

void IOHandler::read(bool& input)
{
 read_basic<bool>(input);
}

void IOHandler::read(fcmplx& input)
{
 read_complex<float>(input);
}
 
void IOHandler::read(dcmplx& input)
{
 read_complex<double>(input);
}


void IOHandler::read(std::vector<char>& input)
{
 read_basic_vector<char>(input);
}

void IOHandler::read(std::vector<int>& input)
{
 read_basic_vector<int>(input);
}

void IOHandler::read(std::vector<unsigned int>& input)
{
 read_basic_vector<unsigned int>(input);
}

void IOHandler::read(std::vector<long int>& input)
{
 read_basic_vector<long int>(input);
}

void IOHandler::read(std::vector<unsigned long int>& input)
{
 read_basic_vector<unsigned long int>(input);
}

void IOHandler::read(std::vector<float>& input)
{
 read_basic_vector<float>(input);
}

void IOHandler::read(std::vector<double>& input)
{
 read_basic_vector<double>(input);
}

void IOHandler::read(std::vector<fcmplx>& input)
{
 read_complex_vector<float>(input);
}

void IOHandler::read(std::vector<dcmplx>& input)
{
 read_complex_vector<double>(input);
}

void IOHandler::read(Array<char>& input)
{
 read_array<char>(input);
}

void IOHandler::read(Array<int>& input)
{
 read_array<int>(input);
}

void IOHandler::read(Array<unsigned int>& input)
{
 read_array<unsigned int>(input);
}

void IOHandler::read(Array<long int>& input)
{
 read_array<long int>(input);
}

void IOHandler::read(Array<unsigned long int>& input)
{
 read_array<unsigned long int>(input);
}

void IOHandler::read(Array<float>& input)
{
 read_array<float>(input);
}

void IOHandler::read(Array<double>& input)
{
 read_array<double>(input);
}

void IOHandler::read(Array<fcmplx>& input)
{
 read_array<fcmplx>(input);
}

void IOHandler::read(Array<dcmplx>& input)
{
 read_array<dcmplx>(input);
}


void IOHandler::multi_read(char* input, int n)
{ 
 multi_read_basic<char>(input,n);
}

void IOHandler::multi_read(int* input, int n)
{ 
 multi_read_basic<int>(input,n);
}

void IOHandler::multi_read(unsigned int* input, int n)
{ 
 multi_read_basic<unsigned int>(input,n);
}

void IOHandler::multi_read(long int* input, int n)
{ 
 multi_read_basic<long int>(input,n);
}

void IOHandler::multi_read(unsigned long int* input, int n)
{ 
 multi_read_basic<unsigned long int>(input,n);
}

void IOHandler::multi_read(long long int* input, int n)
{ 
 multi_read_basic<long long int>(input,n);
}

void IOHandler::multi_read(unsigned long long int* input, int n)
{ 
 multi_read_basic<unsigned long long int>(input,n);
}

void IOHandler::multi_read(float* input, int n)
{ 
 multi_read_basic<float>(input,n);
}

void IOHandler::multi_read(double* input, int n)
{ 
 multi_read_basic<double>(input,n);
}

void IOHandler::multi_read(bool* input, int n)
{ 
 multi_read_basic<bool>(input,n);
}

void IOHandler::multi_read(fcmplx* input, int n)
{
 multi_read_complex<float>(input,n);
}

void IOHandler::multi_read(dcmplx* input, int n)
{ 
 multi_read_complex<double>(input,n);
}

// ********************************************************************


    // the number of bytes that these quantities occupy in an
    // IOHandler file

size_t IOHandler::numbytes(const std::string& data) const
{
 return sizeof(int)+data.length();
}

size_t IOHandler::numbytes(const char& data) const
{
 return numbytes_basic<char>(data);
}

size_t IOHandler::numbytes(const int& data) const
{
 return numbytes_basic<int>(data);
}

size_t IOHandler::numbytes(const unsigned int& data) const
{
 return numbytes_basic<unsigned int>(data);
}

size_t IOHandler::numbytes(const long int& data) const
{
 return numbytes_basic<long int>(data);
}

size_t IOHandler::numbytes(const unsigned long int& data) const
{
 return numbytes_basic<unsigned long int>(data);
}

size_t IOHandler::numbytes(const long long int& data) const
{
 return numbytes_basic<long long int>(data);
}

size_t IOHandler::numbytes(const unsigned long long int& data) const
{
 return numbytes_basic<unsigned long long int>(data);
}

size_t IOHandler::numbytes(const float& data) const
{
 return numbytes_basic<float>(data);
}

size_t IOHandler::numbytes(const double& data) const
{
 return numbytes_basic<double>(data);
}

size_t IOHandler::numbytes(const bool& data) const
{
 return numbytes_basic<bool>(data);
}

size_t IOHandler::numbytes(const fcmplx& data) const
{
 return numbytes_complex<float>(data);
}
 
size_t IOHandler::numbytes(const dcmplx& data) const
{
 return numbytes_complex<double>(data);
}


size_t IOHandler::numbytes(const std::vector<char>& data) const
{
 return numbytes_basic_vector<char>(data);
}

size_t IOHandler::numbytes(const std::vector<int>& data) const
{
 return numbytes_basic_vector<int>(data);
}

size_t IOHandler::numbytes(const std::vector<unsigned int>& data) const
{
 return numbytes_basic_vector<unsigned int>(data);
}

size_t IOHandler::numbytes(const std::vector<long int>& data) const
{
 return numbytes_basic_vector<long int>(data);
}

size_t IOHandler::numbytes(const std::vector<unsigned long int>& data) const
{
 return numbytes_basic_vector<unsigned long int>(data);
}

size_t IOHandler::numbytes(const std::vector<float>& data) const
{
 return numbytes_basic_vector<float>(data);
}

size_t IOHandler::numbytes(const std::vector<double>& data) const
{
 return numbytes_basic_vector<double>(data);
}

size_t IOHandler::numbytes(const std::vector<fcmplx>& data) const
{
 return numbytes_complex_vector<float>(data);
}

size_t IOHandler::numbytes(const std::vector<dcmplx>& data) const
{
 return numbytes_complex_vector<double>(data);
}

size_t IOHandler::numbytes(const Array<char>& data) const
{
 return numbytes_array<char>(data);
}

size_t IOHandler::numbytes(const Array<int>& data) const
{
 return numbytes_array<int>(data);
}

size_t IOHandler::numbytes(const Array<unsigned int>& data) const
{
 return numbytes_array<unsigned int>(data);
}

size_t IOHandler::numbytes(const Array<long int>& data) const
{
 return numbytes_array<long int>(data);
}

size_t IOHandler::numbytes(const Array<unsigned long int>& data) const
{
 return numbytes_array<unsigned long int>(data);
}

size_t IOHandler::numbytes(const Array<float>& data) const
{
 return numbytes_array<float>(data);
}

size_t IOHandler::numbytes(const Array<double>& data) const
{
 return numbytes_array<double>(data);
}

size_t IOHandler::numbytes(const Array<fcmplx>& data) const
{
 return numbytes_array<fcmplx>(data);
}

size_t IOHandler::numbytes(const Array<dcmplx>& data) const
{
 return numbytes_array<dcmplx>(data);
}



// ***************************************************************

