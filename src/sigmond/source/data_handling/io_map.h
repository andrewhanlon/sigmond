#ifndef IO_MAP_H
#define IO_MAP_H

#include "io_map_fstream.h"
#include "io_map_hdf5.h"
#include "xml_handler.h"

 // *********************************************************************************
 // *                                                                               *
 // *       class IOMap:       random access input/output                           *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *  Objects of class IOMap handle the input and output of one type of data       *
 // *  of class "V" using one type of key of class "K". Essentially, an IOMap       *
 // *  object is like a C++ map, but instead of dealing with memory, it deals with  *
 // *  data on a disk.  The steps in using the class are as follows:                *
 // *                                                                               *
 // *   (1) define a class for the keys (see below)                                 *
 // *   (2) define a class for the values (see below)                               *
 // *   (3) create an IOMap object                                                  *
 // *   (4) open a file by calling  openReadOnly, openNew, or openUpdate            *
 // *   (5) use "put" to insert data into the file                                  *
 // *   (6) use "get" to retrieve data from the file                                *
 // *   (7) close the file (or destroy the IOMap object)                            *
 // *                                                                               *
 // *    class RecordType{....};                                                    *
 // *    class DataType{....};                                                      *
 // *    IOMap<RecordType,DataType> iomap;                                          *
 // *    K key(...); V val;                                                         *
 // *    iomap.get(key,val);                                                        *
 // *    iomap.put(key,val);                                                        *
 // *                                                                               *
 // *  More details on the usage are below.                                         *
 // *                                                                               *
 // *  This class accesses the data using a base class IOMapBase.  Several          *
 // *  IO map handlers can then be derived from this base class in order to         *
 // *  allow multiple file formats.  Currently, two different file formats are      *
 // *  available:                                                                   *
 // *      (1) original SigMonD fstreams format (IOFSTRMap with IOFSTRHandler)      *
 // *      (2) HDF5 format (IOHDF5Map with IOHDF5Handler)                           *
 // *  The record key class and values have different requirements for these        *
 // *  different formats.  These differences will be emphasized below.              *
 // *                                                                               *
 // *  The record key class should have the following features:                     *
 // *                                                                               *
 // *    (1) since used in a C++ map, a less than operator must be define           *
 // *              const K& K::operator<(const K& rhs);                             *
 // *    (2) a numbytes(ioh,K) function, where ioh is an IOFSTRHandler, must be     *
 // *         defined to give the number of bytes each key occupies in an           *
 // *         IOFSTRHandler file                                                    *
 // *    (3) a copy constructor K(const K& in) must be defined                      *
 // *          (a default constructor is not needed)                                *
 // *    (4) a multi_read(ioh, vector<K>&,n) must be defined to read n keys         *
 // *    (5) a multi_write(ioh, const vector<K>&) must be defined                   *
 // *    (6) a K.serialize() function member which returns a string to represent    *
 // *         the data, and a constructor K(string) which can assign the            *
 // *         data from the string created by "serialize"                           *
 // *                                                                               *
 // *  The keys are used in a C++ map so keeping the keys small makes for a more    *
 // *  efficient search.  For the original fstreams format, the number of bytes     *
 // *  must be the same for all key values.  For the HDF5 format, the keys          *
 // *  are expressed as serialized strings.                                         *
 // *                                                                               *
 // *  The value type must have the following features:                             *
 // *                                                                               *
 // *   (1) a write(ioh, const V&) must be defined (ioh is an IOFSTRHandler object) *
 // *   (2) a read(ioh, V&) must be defined                                         *
 // *   (3) a numbytes(ioh,V) must be defined giving number of bytes occupied       *
 // *           by V in an IOFSTRHandler file                                       *
 // *                                                                               *
 // *  The size of the value type does not need to be the same for all values. All  *
 // *  of the basic data types, multi1d<T>, multi2d<T>, multi3d<T>, vector<T>,      *
 // *  Array<T> objects already have the above functions defined.                   *
 // *                                                                               *
 // *  IOMap files can be opened using one of three open routines:                  *
 // *                                                                               *
 // *     (1) openReadOnly  -- fails if the file does not exist or read not allowed *
 // *     (2) openNew -- deletes any existing file then creates a new file          *
 // *     (3) openUpdate -- updates existing file or creates new                    *
 // *                                                                               *
 // *  An ID string is used to identify an IOMap file. You should                   *
 // *  choose a string that is based on the key and value types.  During an open    *
 // *  of an existing file, an exact match of the ID string is needed or the open   *
 // *  fails.  For the original fstreams format, the ID string must be              *
 // *  32-characters.  This does not apply for the HDF5 format.                     *
 // *                                                                               *
 // *  During an open of a new file, a header string is written.  This can be of    *
 // *  any length.  During an open of an existing file, the header string is        *
 // *  read and returned (there is one read-only open routine that does not read    *
 // *  the header string).                                                          *
 // *                                                                               *
 // *  During an open of a new file, you can request little endian ('L') format,    *
 // *  big endian ('B') format, or native ('N') format.  You can use check sums     *
 // *  or not.  If check sums are included in a file, they will continue to be      *
 // *  included in any future insertions, even if not used.  During reading, you    *
 // *  can ignore checksums even if they are included in the file.                  *
 // *                                                                               *
 // *  Inserting new data adds the data to the file.  If you attempt to add data    *
 // *  whose key already exists in the file, the data in the file will be           *
 // *  overwritten if overwrites are allowed and if, for the original fstreams      *
 // *  format, the size of the new data cannot larger than the size of the data in  *
 // *  the file (this restriction does not apply to the HDF5 format). To simplify   *
 // *  matters, no erase member is available.  If you really need to erase records, *
 // *  read a file and copy the records you wish to keep to a new file.             *
 // *                                                                               *
 // *  All errors that occur during a "get" or a "put" throw a string (so you can   *
 // *  output more meaning information). All other errors are fatal and cause       *
 // *  an abort.                                                                    *
 // *                                                                               *
 // *  The member "keepKeys" is used to limit attention to a subset of keys.        *
 // *  It is useful if only a few keys are needed, so only those keys are           *
 // *  kept in memory.  The "keepKeys" members returns true if all requested        *
 // *  keys are available, false if some are missing (not available).               *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************



template<typename K, typename V>
class IOMap
{

    char m_file_format;    // 'F' fstreams (default), 'H' hdf5
    IOMapBase<K,V> *m_iomap_ptr;

       // Determine format of file if it exists; if does not exist, set it
       // to "format_if_not_exist".  If "format_if_not_exist" is 'U' for unknown,
       // abort if file does not exist

    void get_file_format(const std::string& filename, char use_format_if_not_exist)
    {
     delete m_iomap_ptr;
     std::string ID;
     std::string fname(filename);
     size_t pos=filename.find("[");
     if (pos!=std::string::npos){
        fname=filename.substr(0,pos);}
     IOFSTRHandler iohA;
     IOHDF5Handler iohB;
     char format_if_not_exist=use_format_if_not_exist;
#ifdef DEFAULT_FSTREAM
     if (use_format_if_not_exist=='D') format_if_not_exist='F';
#else
     if (use_format_if_not_exist=='D') format_if_not_exist='H';
#endif
     if (iohA.peekID(ID,fname)){
        m_file_format='F';  // fstreams format
        m_iomap_ptr=new IOFSTRMap<K,V>;}
     else if (iohB.peekID(ID,fname)){
        m_file_format='H';  // HDF5 format
        m_iomap_ptr=new IOHDF5Map<K,V>;}
     else if (format_if_not_exist=='F'){
        m_file_format='F';  // fstreams format
        m_iomap_ptr=new IOFSTRMap<K,V>;}
     else if (format_if_not_exist=='H'){
        m_file_format='H';  // HDF5 format
        m_iomap_ptr=new IOHDF5Map<K,V>;}
     else{
        std::cout << "IOMap cannot determine file format of "<<filename<<std::endl;
        exit(1);}
    }


    void reset_file_format(const std::string& filename, char use_file_format)
    {
     delete m_iomap_ptr;
     char file_format=use_file_format;
#ifdef DEFAULT_FSTREAM
     if (use_file_format=='D') file_format='F';
#else
     if (use_file_format=='D') file_format='H';
#endif
     if (file_format=='F'){
        m_file_format='F';  // fstreams format
        m_iomap_ptr=new IOFSTRMap<K,V>;}
     else if (file_format=='H'){
        m_file_format='H';  // HDF5 format
        m_iomap_ptr=new IOHDF5Map<K,V>;}
     else{
        std::cout << "Invalid file format requested in IOMap for file "<<filename<<std::endl;
        exit(1);}
    }

  public:
 
#ifdef DEFAULT_FSTREAM
    IOMap() :  m_file_format('F'), m_iomap_ptr(new IOFSTRMap<K,V>) {}
#else
    IOMap() :  m_file_format('H'), m_iomap_ptr(new IOHDF5Map<K,V>) {}
#endif

            // read only open, returns header string
            
    void openReadOnly(const std::string& filename, 
                      const std::string& filetype_id,
                      std::string& header, bool turn_on_checksum=false)
     {get_file_format(filename,'U');
      m_iomap_ptr->openReadOnly(filename,filetype_id,header,turn_on_checksum);}

            // read only open, ignores header string

    void openReadOnly(const std::string& filename, 
                      const std::string& filetype_id,
                      bool turn_on_checksum=false)
     {get_file_format(filename,'U');
      m_iomap_ptr->openReadOnly(filename,filetype_id,turn_on_checksum);}

            // open a new file in read/write mode, writes the header string (fails 
            // if the file exists and "fail_if_exists" is true; if "fail_if_exists"
            // is false, deletes the existing file to start a new file)

    void openNew(const std::string& filename, const std::string& filetype_id, 
                 const std::string& header, bool fail_if_exists=true, char endianness='N',
                 bool turn_on_checksum=false, bool overwrites_allowed=false, 
#ifdef DEFAULT_FSTREAM
                 char file_format='F')
#else
                 char file_format='H')
#endif
     {reset_file_format(filename,file_format);
      m_iomap_ptr->openNew(filename,filetype_id,header,fail_if_exists,endianness,
                           turn_on_checksum,overwrites_allowed);}
                 
            // open a file in read/write mode; if file exists, the header
            // string is read and returned in "header" and writes will update
            // the existing file; otherwise, a new file is created (in which 
            // case, the header string is needed as input so it can be written
            // into the new file)

    void openUpdate(const std::string& filename, const std::string& filetype_id, 
                    std::string& header, char endianness='N', 
                    bool turn_on_checksum=false, bool overwrites_allowed=false, 
#ifdef DEFAULT_FSTREAM
                    char format_if_new='F')
#else
                    char format_if_new='H')
#endif
     {get_file_format(filename,format_if_new);
      m_iomap_ptr->openUpdate(filename,filetype_id,header,endianness,turn_on_checksum, 
                              overwrites_allowed);}

    ~IOMap() {delete m_iomap_ptr;}

    void close() {m_iomap_ptr->close();}



    std::string getHeader() {return m_iomap_ptr->getHeader();}
    
    std::string getFileName() const {return m_iomap_ptr->getFileName();}

    bool isOpen() const {return m_iomap_ptr->isOpen();}

    bool isNewFile() const {return m_iomap_ptr->isNewFile();}

    bool isOverwriteOn() const {return m_iomap_ptr->isOverwriteOn();}


    bool areChecksumsInFile() const {return m_iomap_ptr->areChecksumsInFile();}
      
    bool isFileLittleEndian() const {return m_iomap_ptr->isFileLittleEndian();}
   
    bool isFileBigEndian() const {return m_iomap_ptr->isFileBigEndian();}



    void setHighVerbosity() {m_iomap_ptr->setHighVerbosity();}
    
    void setMediumVerbosity() {m_iomap_ptr->setMediumVerbosity();}
   
    void setNoVerbosity() {m_iomap_ptr->setNoVerbosity();}

    void setDisallowOverwrites() {m_iomap_ptr->setDisallowOverwrites();}

    void setAllowOverwrites() {m_iomap_ptr->setAllowOverwrites();}



    void put(const K& key, const V& val) {m_iomap_ptr->put(key,val);}
    
    void get(const K& key, V& val) {m_iomap_ptr->get(key,val);}

    bool get_maybe(const K& key, V& val) {return m_iomap_ptr->get_maybe(key,val);}

    bool exist(const K& key) const {return m_iomap_ptr->exist(key);} 
    
    void flush() {m_iomap_ptr->flush();}

    std::string get_format() const
    {if (m_file_format=='F') return std::string(" in fstreams format");
     else return std::string(" in HDF5 format");}


    unsigned int size() const {return m_iomap_ptr->size();}

    void getKeys(std::vector<K>& keys) const {m_iomap_ptr->getKeys(keys);}

    void getKeys(std::set<K>& keys) const {m_iomap_ptr->getKeys(keys);}

    bool keepKeys(const std::set<K>& keys_to_keep) {return m_iomap_ptr->keepKeys(keys_to_keep);}

          // "peek" assumes file is not opened yet; it opens the file, peeks, and
          // then closes the file

    bool peekHeader(std::string& header, const std::string& filename, 
                    const std::string& filetype_id)
     {IOFSTRMap<K,V> iomA;
      bool flag=iomA.peekHeader(header,filename,filetype_id);
      if (!flag){
         IOHDF5Map<K,V> iomB;
         flag=iomB.peekHeader(header,filename,filetype_id);}
      return flag;}

};


// ************************************************************************


inline bool IOMapPeekID(std::string& ID, const std::string& filename)
{
 IOFSTRHandler iohA;
 if (iohA.peekID(ID,filename)) return true;
 IOHDF5Handler iohB;
 std::string fname(filename);
 size_t pos=filename.find("[");
 if (pos!=std::string::npos){
    fname=filename.substr(0,pos);}
 if (iohB.peekID(ID,fname)) return true;
 return false;
}


// ************************************************************************

    
     //  A simple integer key class that you can use right out of the box.

class UIntKey
{
   unsigned int value;

 public:

   UIntKey() : value(0) {}
   UIntKey(int inval) : value(inval) {}
   UIntKey(const UIntKey& in) : value(in.value) {}
   UIntKey& operator=(const UIntKey& in) {value=in.value; return *this;}
   UIntKey(XMLHandler &xmlin) {xmlread(xmlin,"KeyValue",value,"UIntKey");}

   UIntKey(const std::string& instr)
   {if (instr.substr(0,8)!="UIntKey_"){
       throw(std::invalid_argument("Invalid string input to UIntKey"));}
    extract_from_string(instr.substr(8,instr.length()-8),value);}

   ~UIntKey() {}

   bool operator<(const UIntKey& rhs) const {return (value<rhs.value);}
   bool operator==(const UIntKey& rhs) const {return (value==rhs.value);}
   bool operator!=(const UIntKey& rhs) const {return (value!=rhs.value);}

   unsigned int getValue() const {return value;}
   void assign(unsigned int ival) {value=ival;}

   void output(XMLHandler& xmlw) const 
   {xmlw.set_root("KeyValue",make_string(getValue()));}

   std::string output() const 
   {XMLHandler xmlout; output(xmlout); return xmlout.output();}
   
   std::string serialize() const
   {return std::string("UIntKey_")+make_string(value);}

   explicit UIntKey(const unsigned int* buf) {value=*buf;}
   static int numints() {return 1;} 
   size_t numbytes() const {return sizeof(unsigned int);}
   void copyTo(unsigned int* buf) const { *buf=value;}

};



inline size_t numbytes(IOFSTRHandler& ioh, const UIntKey& uikey)
{
 return uikey.numbytes();
}


#endif
