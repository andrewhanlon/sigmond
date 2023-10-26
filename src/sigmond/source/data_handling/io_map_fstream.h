#ifndef IO_MAP_LAPH_H
#define IO_MAP_LAPH_H

#include "io_handler_fstream.h"
#include "io_map_base.h"
#include <map>
#include <set>

#ifndef NO_CXX11
#include <type_traits>
#endif
 
#define UINTSIZE     4
//#define SIZETSIZE    4
#define SIZETSIZE    8
#define ULONGLONG    8


 // *********************************************************************************
 // *                                                                               *
 // *       class IOFSTRMap:       random access input/output                       *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *  Objects of class IOFSTRMap handle the input and output of one type of data   *
 // *  of class "V" using one type of key of class "K". Essentially, an IOFSTRMap   *
 // *  object is like a C++ map, but instead of dealing with memory, it deals with  *
 // *  data on a disk.  The steps in using the class are as follows:                *
 // *                                                                               *
 // *   (1) define a class for the keys (see below)                                 *
 // *   (2) define a class for the values (see below)                               *
 // *   (3) create an IOFSTRMap object                                              *
 // *   (4) open a file by calling  openReadOnly, openNew, or openUpdate            *
 // *   (5) use "put" to insert data into the file                                  *
 // *   (6) use "get" to retrieve data from the file                                *
 // *   (7) close the file (or destroy the IOFSTRMap object)                        *
 // *                                                                               *
 // *    class RecordType{....};                                                    *
 // *    class DataType{....};                                                      *
 // *    IOFSTRMap<RecordType,DataType> IOFSTRMap;                                  *
 // *    K key(...); V val;                                                         *
 // *    IOFSTRMap.get(key,val);                                                    *
 // *    IOFSTRMap.put(key,val);                                                    *
 // *                                                                               *
 // *  More details on the usage are below.                                         *
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
 // *                                                                               *
 // *  The number of bytes must be the same for all key values.  The keys are used  *
 // *  in a C++ map so keeping the keys small makes for a more efficient search.    *
 // *  Also, all of the keys stored in the file are read and written in a single    *
 // *  read/write.                                                                  *
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
 // *  IOFSTRMap files can be opened using one of three open routines:              *
 // *                                                                               *
 // *     (1) openReadOnly  -- fails if the file does not exist or read not allowed *
 // *     (2) openNew -- deletes any existing file then creates a new file          *
 // *     (3) openUpdate -- updates existing file or creates new                    *
 // *                                                                               *
 // *  A 32-character ID string is used to identify an IOFSTRMap file. You should   *
 // *  choose a string that is based on the key and value types.  During an open    *
 // *  of an existing file, an exact match of the ID string is needed or the open   *
 // *  fails.                                                                       *
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
 // *  overwritten if overwrites are allowed and if the size of the new data is     *
 // *  not larger than the size of the data in the file. To simplify matters, no    *
 // *  erase member is available.  If you really need to erase records, read a      *
 // *  file and copy the records you wish to keep to a new file.                    *
 // *                                                                               *
 // *  Upon closing a file, the map locations are written to the end of the file.   *
 // *  If the program aborts, the map locations do not get written and the file     *
 // *  is corrupted.  To prevent this, use the "flush" command which immediately    *
 // *  writes the current record locations at the end of the file.                  *
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
 // *  A completed IOFSTRMap file contains (in the following order):                *
 // *   (1) endian character 'L' or 'B'                                             *
 // *   (2) a 32-character ID string                                                *
 // *   (3) location where the map is stored in the file (8 bytes)                  *
 // *   (4) character 'Y' or 'N' describing if check sums stored in file            *
 // *   (5) header string of any length                                             *
 // *   (6) the data records one after another                                      *
 // *   (7) the file map describing locations of the records and keys               *
 // *   (8) an ending character 'E'                                                 *
 // *                                                                               *
 // *  Each record contains (in the following order):                               *
 // *   (1) the size in bytes of the record (excluding checksum)                    *
 // *   (2) the data itself in binary format                                        *
 // *   (3) the checksum of the data (optionally)                                   *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************

/*
   //  Below is a sample class for an IOFSTRMap key.  Use this as a starting
   //  point for defining your own key class.


class TwoSpin
{

    unsigned int s1,s2;    // each value between 0 and 3 (say)
    
  public:
  
    TwoSpin(int in1, int in2);   // no default constructor by design !!
    
    TwoSpin(const TwoSpin& rhs) : s1(rhs.s1), s2(rhs.s2) {}
    
    TwoSpin& operator=(const TwoSpin& rhs)
     {s1=rhs.s1; s2=rhs.s2; return *this;}
    
    std::string output() const;
    
    bool operator<(const TwoSpin& rhs) const
    { return ((s1<rhs.s1) || ( (s1==rhs.s1)&&(s2<rhs.s2) ) ); }
    
    friend void multi_write(IOFSTRHandler& ioh, const vector<TwoSpin>& output);

};

   // no default constructor is needed
   
inline TwoSpin::TwoSpin(int in1, int in2)
{
 if ((in1<0)||(in1>3)||(in2<0)||(in2>3)){
    QDPIO::cerr << "Invalid TwoSpin initialization"<<endl;
    QDP_abort(1);}
 s1=static_cast<unsigned int>(in1);
 s2=static_cast<unsigned int>(in2);
}

inline std::string  TwoSpin::output() const
{
 stringstream oss;
 oss << "("<< s1 <<", "<< s2 <<")";
 return oss.str();
}

    // the size the key occupies in the IOFSTRHandler file (in bytes)
 
size_t numbytes(IOFSTRHandler& ioh, const TwoSpin& ss)
{
 return 2*sizeof(unsigned int);
}

   // copies TwoSpin members into an integer array, then does an
   // integer IOFSTRHandler multi_write (handles byte swapping)

void multi_write(IOFSTRHandler& ioh, const vector<TwoSpin>& output)
{
 int n=output.size();
 if (n<1) return;
 unsigned int *buf=new(nothrow) unsigned int[2*n];
 if (!buf){ 
    QDPIO::cerr << "could not allocate buffer for TwoSpin multiwrite"<<endl; 
    QDP_abort(1);}
 int k=0;
 for (vector<TwoSpin>::const_iterator it=output.begin();it!=output.end();it++){
    buf[k++]=it->s1; buf[k++]=it->s2;}
 ioh.multi_write(buf,2*n);   // int write handles byte swapping if needed
 delete [] buf;
}

   // does an integer IOFSTRHandler multi_read (which handle byte swapping,
   // then copies data into the TwoSpin members

void multi_read(IOFSTRHandler& ioh, vector<TwoSpin>& input, int n)
{
 input.clear();
 if (n<1) return;
 input.reserve(n);  // no default constructor needed here
 unsigned int *buf=new(nothrow) unsigned int[2*n];
 if (!buf){ 
    QDPIO::cerr << "could not allocate buffer for TwoSpin multiread"<<endl; 
    QDP_abort(1);}
 ioh.multi_read(buf,2*n);   // read into ints handles byte swapping if needed
 for (int k=0;k<n;k++)
    input.push_back(TwoSpin(buf[2*k],buf[2*k+1]));
 delete [] buf;
}
*/

// ******************************************************************

         //   Helper routines for simple key classes that are 
         //   a small number of unsigned ints.  Just ensure that
         //   the key class "T" has member functions
         //      static int numints();  <- numints in each key  
         //      void copyTo(unsigned int *buf) const;
         //      T(const unsigned int *buf);   // constructor
         //      size_t numbytes() const <- number of bytes of each key


template <typename T>
void multi_write(IOFSTRHandler& ioh, const std::vector<T>& output)
{
 int n=output.size();
 if (n<1) return;
 int nint=T::numints();
 std::vector<unsigned int> buf(n*nint);
 int k=0;
 for (typename std::vector<T>::const_iterator it=output.begin();it!=output.end();it++){
    it->copyTo(&buf[k]); k+=nint;}
 ioh.multi_write(&buf[0],n*nint);   // int write handles byte swapping if needed
}

template <typename T>
void multi_read(IOFSTRHandler& ioh, std::vector<T>& input, int n)
{
 input.clear();
 if (n<1) return;
 input.reserve(n); 
 int nint=T::numints();
 std::vector<unsigned int> buf(nint*n);
 ioh.multi_read(&buf[0],nint*n);   // read into ints handles byte swapping if needed
 for (int k=0;k<n*nint;k+=nint)
    input.push_back(T(&buf[k]));
}



 // ********************************************************
 // *                                                      *
 // *              The main event: the IOFSTRMap           *
 // *                                                      *
 // ********************************************************


template<typename K, typename V>
class IOFSTRMap  : public IOMapBase<K,V>
{

     typedef IOFSTRHandler::pos_type pos_type; 
     typedef IOFSTRHandler::off_type off_type; 
     typedef unsigned long long  w_pos;
    
     mutable IOFSTRHandler ioh;
     mutable std::map<K, pos_type> file_map;
     bool checksums_in_file;
     bool use_checksums;
     bool allow_overwrites;
     bool verbose1,verbose2;

     w_pos new_write_pos;
     bool write_map_on_close;

      // disallow copying
     IOFSTRMap(const IOFSTRMap&);
     IOFSTRMap(IOFSTRMap&);
     IOFSTRMap& operator=(const IOFSTRMap&);
     IOFSTRMap& operator=(IOFSTRMap&);

  public:
 
    typedef IOFSTRHandler::OpenMode OpenMode;

    IOFSTRMap();

            // read only open, returns header string
            
    virtual void openReadOnly(const std::string& filename, 
                              const std::string& filetype_id,
                              std::string& header, bool turn_on_checksum=false);

            // read only open, ignores header string

    virtual void openReadOnly(const std::string& filename, 
                              const std::string& filetype_id,
                              bool turn_on_checksum=false);

            // open a new file in read/write mode, writes the header string (fails 
            // if the file exists and "fail_if_exists" is true; if "fail_if_exists"
            // is false, deletes the existing file to start a new file)

    virtual void openNew(const std::string& filename, 
                         const std::string& filetype_id, 
                         const std::string& header,  
                         bool fail_if_exists=true, char endianness='N',
                         bool turn_on_checksum=false, bool overwrites_allowed=false);

            // open a file in read/write mode; if file exists, the header
            // string is read and returned in "header" and writes will update
            // the existing file; otherwise, a new file is created (in which 
            // case, the header string is needed as input so it can be written
            // into the new file)

    virtual void openUpdate(const std::string& filename, 
                            const std::string& filetype_id, 
                            std::string& header, 
                            char endianness='N', bool turn_on_checksum=false, 
                            bool overwrites_allowed=false);

    virtual ~IOFSTRMap() {close();}

    virtual void close();



    virtual std::string getHeader();  // file must be open
    
        // Version that assumes file is not open; file closed afterwards.
        // Returns false if file cannot be opened.
        
    bool peekHeader(std::string& header, const std::string& filename, 
                    const std::string& filetype_id);
    
    virtual std::string getFileName() const { return ioh.getFileName(); }

    virtual bool isOpen() const { return ioh.isOpen(); }

    virtual bool isNewFile() const { return ioh.isNewFile(); }

    virtual bool isOverwriteOn() const { return allow_overwrites; }


    virtual bool areChecksumsInFile() const { return checksums_in_file; }
      
    virtual bool isFileLittleEndian() const { return ioh.isFileLittleEndian(); }
   
    virtual bool isFileBigEndian() const { return ioh.isFileBigEndian(); }



    virtual void setHighVerbosity();

    virtual void setMediumVerbosity();
   
    virtual void setNoVerbosity() { verbose1=verbose2=false; }

    virtual void setDisallowOverwrites() { allow_overwrites=false; }

    virtual void setAllowOverwrites() { allow_overwrites=true; }



    virtual void put(const K& key, const V& val);
    
    virtual void get(const K& key, V& val);  // throws exception or aborts if fails

    virtual bool get_maybe(const K& key, V& val);  // returns false is fails, true otherwise

    virtual bool exist(const K& key) const;
    
    virtual void flush();  // puts file in finalized state so no data loss if abort occurs
     


    virtual unsigned int size() const { return file_map.size(); }

    virtual void getKeys(std::vector<K>& keys) const;

    virtual void getKeys(std::set<K>& keys) const;


//    void outputContents(const std::string& logfile, bool printdata);


    virtual bool keepKeys(const std::set<K>& keys_to_keep);   // keep only those keys in "keys_to_keep"
                          // return true if all keys in "keys_to_keep" are available,
                          // false otherwise


  private:

    void initialize(const std::string& filename, OpenMode mode,
                    const std::string& filetype_id, char endianness,
                    bool turn_on_checksum, std::string& header, 
                    bool read_header, bool overwrites);

    void readMap(w_pos mapstart);
    
    void writeMap(w_pos mapstart);
    
    void check_for_failure(bool errcond, const std::string& msg, bool abort=true);

#ifdef NO_CXX11
         // for static (compile time) assertion (produces compiler error if not satisfied)
   template <bool b>
   void staticassert()
   { typedef char asserter[b?1:-1]; }
#endif

};

    
    

template <typename K, typename V>
IOFSTRMap<K,V>::IOFSTRMap() : checksums_in_file(false), use_checksums(false), 
                              allow_overwrites(false), verbose1(false), verbose2(false), 
                              new_write_pos(0), write_map_on_close(false)
{
          // rely on specific sizes in the file so check these sizes
#ifndef NO_CXX11
 static_assert((sizeof(unsigned)==UINTSIZE)&&(sizeof(w_pos)==ULONGLONG)
               &&(sizeof(size_t)==SIZETSIZE),"Invalid data sizes");
#else
 staticassert<(sizeof(unsigned)==UINTSIZE)&&(sizeof(w_pos)==ULONGLONG)
               &&(sizeof(size_t)==SIZETSIZE)>();
#endif
}


template <typename K, typename V>
void IOFSTRMap<K,V>::openReadOnly(const std::string& filename, 
                                  const std::string& filetype_id,
                                  std::string& header, bool turn_on_checksum)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOFSTRMap!");
 OpenMode mode=IOFSTRHandler::ReadOnly;
 initialize(filename,mode,filetype_id,'N',turn_on_checksum, 
            header,true,false);
}



template <typename K, typename V>
void IOFSTRMap<K,V>::openReadOnly(const std::string& filename, 
                                  const std::string& filetype_id,
                                  bool turn_on_checksum)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOFSTRMap!");
 OpenMode mode=IOFSTRHandler::ReadOnly;
 std::string header;
 initialize(filename,mode,filetype_id,'N',turn_on_checksum, 
            header,false,false);
}


template <typename K, typename V>
void IOFSTRMap<K,V>::openNew(const std::string& filename, 
                             const std::string& filetype_id, 
                             const std::string& header,  
                             bool fail_if_exists, char endianness,
                             bool turn_on_checksum, bool overwrites_allowed)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOFSTRMap!");
 OpenMode mode=(fail_if_exists)?IOFSTRHandler::ReadWriteFailIfExists : IOFSTRHandler::ReadWriteEraseIfExists;
 std::string tmp_header(header);
 initialize(filename,mode,filetype_id,endianness,
            turn_on_checksum,tmp_header,false,overwrites_allowed);
}


template <typename K, typename V>
void IOFSTRMap<K,V>::openUpdate(const std::string& filename, 
                                const std::string& filetype_id, 
                                std::string& header, 
                                char endianness,  bool turn_on_checksum, 
                                bool overwrites_allowed)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOFSTRMap!");
 OpenMode mode=IOFSTRHandler::ReadWriteUpdateIfExists;
 initialize(filename,mode,filetype_id,endianness,
            turn_on_checksum,header,true,overwrites_allowed);
}



template <typename K, typename V>
void IOFSTRMap<K,V>::initialize(const std::string& filename, OpenMode mode,
                                const std::string& filetype_id, char endianness,
                                bool turn_on_checksum, std::string& header, 
                                bool read_header, bool overwrites)
{
 ioh.open(filename,mode,filetype_id,endianness,turn_on_checksum);
 allow_overwrites=overwrites;

 if (ioh.isNewFile()){                // this is a new file
    if (verbose1) std::cout << "IOFSTRMap: opening new file "<<filename<<std::endl;
    char cksum=(turn_on_checksum)?'Y':'N';
    checksums_in_file=use_checksums=turn_on_checksum;
    new_write_pos=0;
    ioh.write(new_write_pos);  // placeholder for the map start
    ioh.write(cksum);
    write(ioh,header);
    new_write_pos=ioh.tell();
    write_map_on_close=true;}
 else{                         // this is an existing file
    if (verbose1) std::cout << "IOFSTRMap: opening already existing file "<<filename<<std::endl;
    ioh.rewind();
    ioh.read(new_write_pos);
    char cksum;
    ioh.read(cksum);
    checksums_in_file=(cksum=='Y')?true:false;
    use_checksums=turn_on_checksum && checksums_in_file;
    if (!use_checksums) ioh.turnOffChecksum();
    if (read_header) ioh.read(header);
    write_map_on_close=false;
    readMap(new_write_pos);
    }

 if (verbose2){
    if (ioh.isChecksumOn()) std::cout << "   Checksums are ON"<<std::endl;
    else std::cout << "   Checksums are OFF"<<std::endl;
    if (ioh.isFileLittleEndian()) std::cout << "   File is little endian"<<std::endl;
    else std::cout << "   File is big endian"<<std::endl;
    if (ioh.isEndianConversionOn()) std::cout << "   Endian conversion is ON"<<std::endl;
    else std::cout << "   Endian conversion is OFF"<<std::endl;
    if (ioh.isReadOnly()) std::cout << "   File is opened in read-only mode"<<std::endl;
    else std::cout << "   File is opened in read-write mode"<<std::endl;
    std::cout << "   Number of key-value pairs is "<<file_map.size()<<std::endl;}
}


template <typename K, typename V>
void IOFSTRMap<K,V>::close()
{
 if (!ioh.isOpen()) return;
 if (write_map_on_close) writeMap(new_write_pos);
 if (verbose1) std::cout << "IOFSTRMap: closed file "<<ioh.getFileName()<<std::endl;
 ioh.close();
 file_map.clear();
}



template <typename K, typename V>
void IOFSTRMap<K,V>::flush()
{
 if (!ioh.isOpen()) return;
 if (write_map_on_close) writeMap(new_write_pos);
 if (verbose1) std::cout << "IOFSTRMap: flushed file "<<ioh.getFileName()<<std::endl;
 write_map_on_close=false;
}


template<typename K, typename V>
void IOFSTRMap<K,V>::writeMap(w_pos mapstart)
{
 if (checksums_in_file) ioh.turnOnChecksum();
 ioh.seekFromStart(static_cast<off_type>(mapstart));
 size_t nrecords=file_map.size();
 ioh.write(nrecords);

 size_t keysize;
 if (nrecords==0) keysize=0;
 else keysize=numbytes(ioh,file_map.begin()->first); 
 ioh.write(keysize);
 
 std::vector<K> keybuf; keybuf.reserve(nrecords);
 w_pos* posbuf=new(std::nothrow) w_pos[nrecords];
 check_for_failure((!posbuf),"allocation error in IOFSTRMap");
 w_pos* posptr=posbuf;
 for (typename std::map<K,pos_type>::const_iterator it=file_map.begin();
          it!=file_map.end();it++){
    keybuf.push_back(it->first);
   *posptr=it->second; posptr++;}
 multi_write(ioh,keybuf);
 ioh.multi_write(posbuf,nrecords);
 delete [] posbuf;
 check_for_failure((size_t(ioh.tell())-size_t(mapstart))!=(nrecords*(keysize+sizeof(mapstart))
                    +2*sizeof(size_t)),"bad write of file map");

 if (checksums_in_file){
    ByteHandler::n_uint32_t checksum=ioh.getChecksum();
    ioh.write(checksum);}
 if (!use_checksums) ioh.turnOffChecksum();
 char eof='E';
 ioh.write(eof);
 ioh.rewind(); 
 ioh.write(mapstart);
}




template<typename K, typename V>
void IOFSTRMap<K,V>::readMap(w_pos mapstart)
{
 check_for_failure(mapstart==0,"IOFSTRMap file "+ioh.getFileName()+" corrupted from prior aborted execution");
 ioh.seekFromStart(static_cast<off_type>(mapstart));
 size_t nrecords;
 ioh.read(nrecords);
 check_for_failure(nrecords>16777216,"Too many records during readMap: bad read?");

 size_t keysize;
 ioh.read(keysize);
 check_for_failure(keysize>1024,"Key size too large during readMap: bad read?");

 if (nrecords>0){
    std::vector<K> keybuf; keybuf.reserve(nrecords);
    w_pos* posbuf=new(std::nothrow) w_pos[nrecords]; 
    check_for_failure((!posbuf),"allocation error in IOFSTRMap"); 
    multi_read(ioh,keybuf,nrecords);
    ioh.multi_read(posbuf,nrecords);
    for (unsigned int k=0;k<nrecords;k++){
       file_map.insert(std::make_pair(keybuf[k],static_cast<pos_type>(posbuf[k])));}
    delete [] posbuf;
    check_for_failure(keysize!=numbytes(ioh,file_map.begin()->first),
                     "Bad keysize in reading file map");}
 check_for_failure((size_t(ioh.tell())-size_t(mapstart))!=(nrecords*(keysize+sizeof(mapstart))
                    +2*sizeof(size_t)),"bad read of file map");
 
 if (use_checksums){
    ByteHandler::n_uint32_t checksumA=ioh.getChecksum();
    ByteHandler::n_uint32_t checksumB;
    ioh.read(checksumB);
    check_for_failure(checksumA!=checksumB,"checksum mismatch: bad read of file map");}
 else if (checksums_in_file)
    ioh.seekFromCurr(sizeof(ByteHandler::n_uint32_t));

 char eof;
 ioh.read(eof);
 check_for_failure(eof!='E',"bad read of file map");
}


template<typename K, typename V>
void IOFSTRMap<K,V>::setHighVerbosity() 
{ 
 verbose1=verbose2=true;
}


template<typename K, typename V>
void IOFSTRMap<K,V>::setMediumVerbosity() 
{
 verbose1=true;
 verbose2=false;
}


template<typename K, typename V>
bool IOFSTRMap<K,V>::exist(const K& key) const 
{
 return (file_map.find(key) == file_map.end()) ? false : true;
}
  

template<typename K, typename V>
std::string IOFSTRMap<K,V>::getHeader()
{
 check_for_failure(!ioh.isOpen(),"IOFSTRMap not open: cannot get header");
 std::string tmp;
 ioh.seekFromStart(sizeof(w_pos)+1);
 ioh.read(tmp);
 return tmp;
}


template<typename K, typename V>
bool IOFSTRMap<K,V>::peekHeader(std::string& header,
                                const std::string& filename, 
                                const std::string& filetype_id)
{
 return ioh.peekString(header,sizeof(w_pos)+1,filename,filetype_id);
}
 
 
template<typename K, typename V>
void IOFSTRMap<K,V>::put(const K& key, const V& val) 
{
 check_for_failure((!ioh.isOpen())||(ioh.isReadOnly()),"Write operation disallowed",false);
 w_pos start_write_pos=new_write_pos;
 size_t sz=numbytes(ioh,val);
 bool flag=exist(key);
 if (flag){
    check_for_failure(!allow_overwrites,"key already in map and overwrites disallowed",false);
    start_write_pos=static_cast<w_pos>(file_map[key]);
    ioh.seekFromStart(static_cast<off_type>(start_write_pos));
    size_t cursize; 
    ioh.read(cursize);
    check_for_failure((cursize<sz),"can only overwrite record if new size is not larger",false);
    if (verbose2) std::cout << "IOFSTRMap::put overwriting existing record at file location "
                            <<start_write_pos<<std::endl;
    }
 else{
    if (verbose2) std::cout << "IOFSTRMap::put of new record at file location "
                            <<start_write_pos<<std::endl;
    write_map_on_close=true;
    file_map.insert(std::make_pair(key,start_write_pos));}
 if (checksums_in_file) ioh.turnOnChecksum();
 ioh.seekFromStart(static_cast<off_type>(start_write_pos));
 ioh.write(sz);
 write(ioh,val);
 size_t recordsize=sizeof(size_t)+sz;
 if (checksums_in_file){
    ByteHandler::n_uint32_t checksum=ioh.getChecksum();
    ioh.write(checksum);
    recordsize+=sizeof(checksum);}
 if (!use_checksums) ioh.turnOffChecksum();
 pos_type curpos=ioh.tell();
 check_for_failure((size_t(curpos)-size_t(start_write_pos))!=recordsize,
                   "bad value write in IOFSTRMap",false);
 if (!flag){
    new_write_pos=static_cast<w_pos>(curpos);
    ioh.rewind(); w_pos zero=0;
    ioh.write(zero);}
} 


template<typename K, typename V>
void IOFSTRMap<K,V>::get(const K& key, V& val) 
{
 check_for_failure((!ioh.isOpen())||(!exist(key)),
                   "Read failed: file not open or key not in file",false);
 off_type start = static_cast<off_type>(file_map.find(key)->second); 
 if (verbose2) std::cout << "IOFSTRMap::get at file location "<<start<<std::endl;
 ioh.seekFromStart(start);
 size_t sz;
 ioh.read(sz);
 read(ioh,val);
 check_for_failure(sz!=numbytes(ioh,val),"value size mismatch",false);
 check_for_failure((size_t(ioh.tell())-size_t(start))!=(sizeof(size_t)+sz),
                   "bad value read in IOFSTRMap",false);
 if (use_checksums){
    ByteHandler::n_uint32_t checksumA=ioh.getChecksum();
    ByteHandler::n_uint32_t checksumB;
    ioh.read(checksumB);
    check_for_failure(checksumA!=checksumB,"checksum mismatch",false);}
}


    // returns false if "get" cannot be achieved, true otherwise

template<typename K, typename V>
bool IOFSTRMap<K,V>::get_maybe(const K& key, V& val) 
{
 if (!ioh.isOpen()) return false;
 typename std::map<K, pos_type>::iterator ft=file_map.find(key);
 if (ft==file_map.end()) return false;
 off_type start = static_cast<off_type>(ft->second); 
 if (verbose2) std::cout << "IOFSTRMap::get at file location "<<start<<std::endl;
 ioh.seekFromStart(start);
 size_t sz;
 ioh.read(sz);
 read(ioh,val);
 if (sz!=numbytes(ioh,val)) return false;
 if ((size_t(ioh.tell())-size_t(start))!=(sizeof(size_t)+sz)) return false;
 if (use_checksums){
    ByteHandler::n_uint32_t checksumA=ioh.getChecksum();
    ByteHandler::n_uint32_t checksumB;
    ioh.read(checksumB);
    if (checksumA!=checksumB) return false;}
 return true;
}


template<typename K, typename V>
void IOFSTRMap<K,V>::getKeys(std::vector<K>& keys) const 
{
 keys.clear();
 keys.reserve(file_map.size());
 for (typename std::map<K,pos_type>::const_iterator it  = file_map.begin();
      it != file_map.end(); ++it){
      keys.push_back(it->first);}
}


template<typename K, typename V>
void IOFSTRMap<K,V>::getKeys(std::set<K>& keys) const 
{
 keys.clear();
 for (typename std::map<K,pos_type>::const_iterator it  = file_map.begin();
      it != file_map.end(); ++it){
      keys.insert(it->first);}
}

            // keep only those keys in "keys_to_keep"
            // return true of all keys in "keys_to_keep" are available,
            // false otherwise

template<typename K, typename V>
bool IOFSTRMap<K,V>::keepKeys(const std::set<K>& keys_to_keep)  
{
 std::map<K, pos_type> newmap;
 for (typename std::set<K>::const_iterator it=keys_to_keep.begin();it!=keys_to_keep.end();it++){
    typename std::map<K, pos_type>::iterator mt=file_map.find(*it);
    if (mt!=file_map.end()) newmap.insert(*mt);} 
 file_map=newmap;
 return (file_map.size()==keys_to_keep.size());
}



template<typename K, typename V>
void IOFSTRMap<K,V>::check_for_failure(bool errcond, const std::string& mesg,
                                   bool abort)
{
 if (!errcond) return;
 //std::cerr << "IOFSTRMap error: "<<mesg<<std::endl;
 if (abort) exit(1);
 else throw(std::invalid_argument(mesg));
}


// **************************************************************
#endif
