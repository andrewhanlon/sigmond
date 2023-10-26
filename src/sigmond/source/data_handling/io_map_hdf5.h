#ifndef IO_MAP_HDF5_H
#define IO_MAP_HDF5_H

#include "io_handler_hdf5.h"
#include "io_map_base.h"
#include "xml_handler.h"
#include <map>
#include <set>
#include <algorithm>
 

 // *********************************************************************************
 // *                                                                               *
 // *       class IOHDF5Map:       random access input/output                       *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *  Objects of class IOHDF5Map handle the input and output of one type of data   *
 // *  of class "V" using one type of key of class "K". Essentially, an IOHDF5Map   *
 // *  object is like a C++ map, but instead of dealing with memory, it deals with  *
 // *  data on a disk.  HDF5 files can be very large.  An IOHDF5Map is used to only *
 // *  access a small part of the file.  You MUST specify a root path (group), then *
 // *  only data in that group will be accessed.  In addition to the keys and data, *
 // *  a header string can be stored under each root path.                          *
 // *                                                                               *
 // *  The root path must be specified in the filename by appending [rootpath].     *
 // *  Example:   /raid/channel/filename[chan1/mom2]                                *
 // *                                                                               *
 // *  The steps in using the class are as follows:                                 *
 // *                                                                               *
 // *   (1) define a class for the keys (see below)                                 *
 // *   (2) define a class for the values (see below)                               *
 // *   (3) create an IOHDF5Map object                                              *
 // *   (4) open a file by calling  openReadOnly, openNew, or openUpdate            *
 // *   (5) use "put" to insert data into the file                                  *
 // *   (6) use "get" to retrieve data from the file                                *
 // *   (7) close the file (or destroy the IOHDF5Map object)                        *
 // *                                                                               *
 // *    class RecordType{....};                                                    *
 // *    class DataType{....};                                                      *
 // *    IOHDF5Map<RecordType,DataType> IOHDF5Map;                                  *
 // *    K key(...); V val;                                                         *
 // *    IOHDF5Map.get(key,val);                                                    *
 // *    IOHDF5Map.put(key,val);                                                    *
 // *                                                                               *
 // *  More details on the usage are below.                                         *
 // *                                                                               *
 // *  The record key class should have the following features:                     *
 // *                                                                               *
 // *    (1) since used in a C++ map, a less than operator must be define           *
 // *              const K& K::operator<(const K& rhs);                             *
 // *    (2) a K.serialize() function member which returns a string to represent    *
 // *         the data, and a constructor K(string) which can assign the            *
 // *         data from the string created by "serialize"                           *
 // *    (3) a copy constructor K(const K& in) must be defined                      *
 // *          (a default constructor is not needed)                                *
 // *                                                                               *
 // *  The value type must have the following features:                             *
 // *                                                                               *
 // *   (1) a write(ioh, const V&) must be defined (ioh is an IOHDF5Handler object) *
 // *   (2) a read(ioh, V&) must be defined                                         *
 // *                                                                               *
 // *  The size of the value type does not need to be the same for all values. All  *
 // *  of the basic data types, multi1d<T>, multi2d<T>, multi3d<T>, vector<T>,      *
 // *  Array<T> objects already have the above functions defined.                   *
 // *                                                                               *
 // *  IOHDF5Map files can be opened using one of three open routines:              *
 // *                                                                               *
 // *     (1) openReadOnly  -- fails if the file does not exist or read not allowed *
 // *     (2) openNew -- deletes any existing file then creates a new file          *
 // *     (3) openUpdate -- updates existing file or creates new                    *
 // *                                                                               *
 // *  An ID string is used to identify an IOHDF5MapHDF5 file. You should choose a  *
 // *  string that is based on the key and value types.  During an open of an       *
 // *  existing file, an exact match of the ID string is needed or the open         *
 // *  fails.                                                                       *
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


 // ********************************************************
 // *                                                      *
 // *              The main event: the IOHDF5Map           *
 // *                                                      *
 // ********************************************************


template<typename K, typename V>
class IOHDF5Map  :  public IOMapBase<K,V>
{

     mutable IOHDF5Handler ioh;
     bool checksums_in_file;
     bool use_checksums;
     bool allow_overwrites;
     bool verbose1,verbose2;

      // disallow copying
     IOHDF5Map(const IOHDF5Map&);
     IOHDF5Map(IOHDF5Map&);
     IOHDF5Map& operator=(const IOHDF5Map&);
     IOHDF5Map& operator=(IOHDF5Map&);

  public:
 
    typedef IOHDF5Handler::OpenMode OpenMode;

    IOHDF5Map();

            // read only open, returns header string
            
    virtual void openReadOnly(const std::string& filename_and_rootpath, 
                              const std::string& filetype_id,
                              std::string& header, bool turn_on_checksum=false);

            // read only open, ignores header string

    virtual void openReadOnly(const std::string& filename_and_rootpath, 
                              const std::string& filetype_id,
                              bool turn_on_checksum=false);

            // open a new file in read/write mode, writes the header string (fails 
            // if the file exists and "fail_if_exists" is true; if "fail_if_exists"
            // is false, deletes the existing file to start a new file)

    virtual void openNew(const std::string& filename_and_rootpath, 
                         const std::string& filetype_id,
                         const std::string& header,  
                         bool fail_if_exists=true, char endianness='N',
                         bool turn_on_checksum=false, bool overwrites_allowed=false);

            // open a file in read/write mode; if file exists, the header
            // string is read and returned in "header" and writes will update
            // the existing file; otherwise, a new file is created (in which 
            // case, the header string is needed as input so it can be written
            // into the new file)

    virtual void openUpdate(const std::string& filename_and_rootpath, 
                            const std::string& filetype_id,
                            std::string& header, 
                            char endianness='N', bool turn_on_checksum=false, 
                            bool overwrites_allowed=false);

    virtual ~IOHDF5Map() {close();}

    virtual void close();



    virtual std::string getHeader();  // file must be open
    
        // Version that assumes file is not open; file closed afterwards.
        // Returns false if file cannot be opened.
        
    bool peekHeader(std::string& header, const std::string& filename_and_rootpath, 
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
     


    virtual unsigned int size() const;

    virtual void getKeys(std::vector<K>& keys) const;

    virtual void getKeys(std::set<K>& keys) const;

    virtual bool keepKeys(const std::set<K>& keys_to_keep);   // keep only those keys in "keys_to_keep"
                          // return true if all keys in "keys_to_keep" are available,
                          // false otherwise

    std::set<std::string> getAvailableRootPaths() const;


  private:

    void initialize(const std::string& filename_and_rootpath, OpenMode mode,
                    const std::string& filetype_id, 
                    char endianness,
                    bool turn_on_checksum, std::string& header, 
                    bool read_header, bool overwrites);

    void check_for_failure(bool errcond, const std::string& msg, bool abort=true) const;
    
    uint get_size() const;
 
    char get_filename_and_path(const std::string& filename_and_rootpath, 
                               std::string& filename, std::string& rootpath); 

};


// *************************************************************************************


   //  Take string containing filename and root path of form  filename[rootpath]
   //  and separates them, returning them in "filename" and "rootpath".
   //  Returns 'B' for both if successful, 'F' if only filename could be extracted,
   //  and 'N' for none if neither could be determined.

template <typename K, typename V>
char IOHDF5Map<K,V>::get_filename_and_path(const std::string& filename_and_rootpath, 
                                           std::string& filename, std::string& rootpath)
{
 filename.clear(); rootpath.clear();
 size_t nlbrack = std::count(filename_and_rootpath.begin(), filename_and_rootpath.end(), '[');
 size_t nrbrack = std::count(filename_and_rootpath.begin(), filename_and_rootpath.end(), ']');
 if ((nlbrack==0)&&(nrbrack==0)){
    filename=tidyString(filename_and_rootpath); 
    if (filename.empty()) return 'N';
    else return 'F';}
 if ((nlbrack!=1)||(nrbrack!=1)) return 'N';
 size_t rstart=filename_and_rootpath.find("[");
 size_t rfinish=filename_and_rootpath.find("]");
 if (rfinish<rstart) return 'N';
 filename=tidyString(filename_and_rootpath.substr(0,rstart));
 rootpath=tidyString(filename_and_rootpath.substr(rstart+1,rfinish-rstart-1));
 if (filename.empty()) return 'N';
 if (rootpath.empty()) return 'F';
 return 'B';
}


template <typename K, typename V>
IOHDF5Map<K,V>::IOHDF5Map() : checksums_in_file(false), use_checksums(false), 
                              allow_overwrites(false), verbose1(false), verbose2(false)
{}


template <typename K, typename V>
void IOHDF5Map<K,V>::openReadOnly(const std::string& filename_and_rootpath, 
                                  const std::string& filetype_id,
                                  std::string& header, bool turn_on_checksum)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOHDF5Map!");
 OpenMode mode=IOHDF5Handler::ReadOnly;
 initialize(filename_and_rootpath,mode,filetype_id,'N',turn_on_checksum, 
            header,true,false);
}



template <typename K, typename V>
void IOHDF5Map<K,V>::openReadOnly(const std::string& filename_and_rootpath, 
                                  const std::string& filetype_id,
                                  bool turn_on_checksum)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOHDF5Map!");
 OpenMode mode=IOHDF5Handler::ReadOnly;
 std::string header;
 initialize(filename_and_rootpath,mode,filetype_id,'N',turn_on_checksum, 
            header,false,false);
}


template <typename K, typename V>
void IOHDF5Map<K,V>::openNew(const std::string& filename_and_rootpath, 
                             const std::string& filetype_id, 
                             const std::string& header,  
                             bool fail_if_exists, char endianness,
                             bool turn_on_checksum, bool overwrites_allowed)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOHDF5Map!");
 OpenMode mode=(fail_if_exists)?IOHDF5Handler::ReadWriteFailIfExists : IOHDF5Handler::ReadWriteEraseIfExists;
 std::string tmp_header(header);
 initialize(filename_and_rootpath,mode,filetype_id,endianness,
            turn_on_checksum,tmp_header,false,overwrites_allowed);
}


template <typename K, typename V>
void IOHDF5Map<K,V>::openUpdate(const std::string& filename_and_rootpath, 
                                const std::string& filetype_id, 
                                std::string& header, 
                                char endianness,  bool turn_on_checksum, 
                                bool overwrites_allowed)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOHDF5Map!");
 OpenMode mode=IOHDF5Handler::ReadWriteUpdateIfExists;
 initialize(filename_and_rootpath,mode,filetype_id,endianness,
            turn_on_checksum,header,true,overwrites_allowed);
}



template <typename K, typename V>
void IOHDF5Map<K,V>::initialize(const std::string& filename_and_rootpath, OpenMode mode,
                                const std::string& filetype_id,
                                char endianness, bool turn_on_checksum, 
                                std::string& header, bool read_header, 
                                bool overwrites)
{
 std::string filename,rootpath;
 char flag=get_filename_and_path(filename_and_rootpath,filename,rootpath);
 check_for_failure(flag=='N',"Problem initializing IOHDF5Map: bad file name or root path");
 ioh.open(filename,mode,filetype_id,endianness,turn_on_checksum);
 if (flag=='F'){
    std::set<std::string> rootpaths(getAvailableRootPaths());
    std::cout << "Problem initializing IOHDF5Map: bad root path"<<std::endl;
    std::cout << "Available root paths are as follows:"<<std::endl;
    for (std::set<std::string>::iterator it=rootpaths.begin();it!=rootpaths.end();++it){
       std::cout << *it<<std::endl;}
    ioh.close();
    exit(1);}
 allow_overwrites=overwrites;

 if (ioh.isNewFile()){                // this is a new file
    if (verbose1) std::cout << "IOHDF5Map: opening new file "<<filename<<std::endl;
    char cksum=(turn_on_checksum)?'Y':'N';
    checksums_in_file=use_checksums=turn_on_checksum;
    ioh.mkdir(rootpath);
    ioh.cd(rootpath);
    ioh.mkdir("Values");
    ioh.write("Header",tidyString(header));
    ioh.write("IncludeCKS",cksum);
    if (checksums_in_file){ 
       ioh.mkdir("Checksums");}}
 else{                         // this is an existing file
    if (verbose1) std::cout << "IOHDF5Map: opening already existing file "<<filename<<std::endl;
    if ((ioh.queryDir(rootpath))||(ioh.isReadOnly())){
       ioh.cd(rootpath);
       char cksum;
       ioh.read("IncludeCKS",cksum);
       checksums_in_file=(cksum=='Y')?true:false;
       use_checksums=turn_on_checksum && checksums_in_file;
       if (!use_checksums) ioh.turnOffChecksum();
       if (read_header) ioh.read("Header",header);}
    else if (!(ioh.isReadOnly())){
       ioh.mkdir(rootpath);
       ioh.cd(rootpath);
       ioh.mkdir("Values");
       char cksum=(turn_on_checksum)?'Y':'N';
       checksums_in_file=use_checksums=turn_on_checksum;
       ioh.write("Header",tidyString(header));
       ioh.write("IncludeCKS",cksum);
       if (checksums_in_file){ 
          ioh.mkdir("Checksums");}}}
 if (verbose2){
    if (ioh.isChecksumOn()) std::cout << "   Checksums are ON"<<std::endl;
    else std::cout << "   Checksums are OFF"<<std::endl;
    if (ioh.isFileLittleEndian()) std::cout << "   File is little endian"<<std::endl;
    else std::cout << "   File is big endian"<<std::endl;
    if (ioh.isEndianConversionOn()) std::cout << "   Endian conversion is ON"<<std::endl;
    else std::cout << "   Endian conversion is OFF"<<std::endl;
    if (ioh.isReadOnly()) std::cout << "   File is opened in read-only mode"<<std::endl;
    else std::cout << "   File is opened in read-write mode"<<std::endl;
    std::cout << "   Number of key-value pairs is "<<get_size()<<std::endl;}
}


template <typename K, typename V>
void IOHDF5Map<K,V>::close()
{
 if (!ioh.isOpen()) return;
 if (verbose1) std::cout << "IOHDF5Map: closed file "<<ioh.getFileName()<<std::endl;
 ioh.close();
}



template <typename K, typename V>
void IOHDF5Map<K,V>::flush()
{
 if (!ioh.isOpen()) return;
 if (verbose1) std::cout << "IOHDF5Map: flushed file "<<ioh.getFileName()<<std::endl;
}



template<typename K, typename V>
void IOHDF5Map<K,V>::setHighVerbosity() 
{ 
 verbose1=verbose2=true;
}


template<typename K, typename V>
void IOHDF5Map<K,V>::setMediumVerbosity() 
{
 verbose1=true;
 verbose2=false;
}


template<typename K, typename V>
bool IOHDF5Map<K,V>::exist(const K& key) const 
{
 check_for_failure(!ioh.isOpen(),"IOHDF5Map not open: cannot perform query"); 
 return ioh.queryData(std::string("Values/")+key.serialize());
}
  

template<typename K, typename V>
std::string IOHDF5Map<K,V>::getHeader()
{
 check_for_failure(!ioh.isOpen(),"IOHDF5Map not open: cannot get header");
 std::string tmp;
 ioh.read("Header",tmp);
 return tmp;
}


template<typename K, typename V>
uint IOHDF5Map<K,V>::get_size() const
{
 check_for_failure(!ioh.isOpen(),"IOHDF5Map not open: cannot get size");
 ioh.cd("Values");
 std::set<std::string> datanames(ioh.getDataNamesInCurrentDir());
 ioh.cd("..");
 return datanames.size();
}


template<typename K, typename V>
uint IOHDF5Map<K,V>::size() const
{
 return get_size();
}


template<typename K, typename V>
bool IOHDF5Map<K,V>::peekHeader(std::string& header,
                                const std::string& filename_and_rootpath, 
                                const std::string& filetype_id)
{
 std::string filename,rootpath;
 char flag=get_filename_and_path(filename_and_rootpath,filename,rootpath);
 if (flag!='B') return false;
 std::string ppath(rootpath); ppath+="/Header";
 return ioh.peekString(header,ppath,filename,filetype_id);
}


template<typename K, typename V>
void IOHDF5Map<K,V>::put(const K& key, const V& val)
{
 check_for_failure((!ioh.isOpen())||(ioh.isReadOnly()),"Write operation disallowed",false);
 std::string keyname(key.serialize());
 bool flag=ioh.queryData(std::string("Values/")+keyname);
 if (flag){
    check_for_failure(!allow_overwrites,"key already in map and overwrites disallowed",false);}
 else{
    if (verbose2) std::cout << "IOHDF5Map::put of new record"<<std::endl;}
 if (checksums_in_file) ioh.turnOnChecksum();
 ioh.write(std::string("Values/")+keyname,val);
 if (checksums_in_file){
    ByteHandler::n_uint32_t checksum=ioh.getChecksum();
    ioh.write(std::string("Checksums/")+keyname,checksum);}
 if (!use_checksums) ioh.turnOffChecksum();
}


template<typename K, typename V>
void IOHDF5Map<K,V>::get(const K& key, V& val) 
{
 check_for_failure(!ioh.isOpen(),"Read failed: file not open",false);
 if (verbose2) std::cout << "IOHDF5Map::get "<<std::endl;
 std::string keyname(key.serialize());
 ioh.read(std::string("Values/")+keyname,val);
 if (use_checksums){
    ByteHandler::n_uint32_t checksumA=ioh.getChecksum();
    ByteHandler::n_uint32_t checksumB;
    ioh.read(std::string("Checksums/")+keyname,checksumB);
    check_for_failure(checksumA!=checksumB,"checksum mismatch",false);}
}


    // returns false if "get" cannot be achieved, true otherwise

template<typename K, typename V>
bool IOHDF5Map<K,V>::get_maybe(const K& key, V& val) 
{
 if (!exist(key)) return false;
 if (verbose2) std::cout << "IOHDF5Map::get_maybe "<<std::endl;
 std::string keyname(key.serialize());
 ioh.read(std::string("Values/")+keyname,val);
 if (use_checksums){
    ByteHandler::n_uint32_t checksumA=ioh.getChecksum();
    ByteHandler::n_uint32_t checksumB;
    ioh.read(std::string("Checksums/")+keyname,checksumB);
    if (checksumA!=checksumB) return false;}
 return true;
}


template<typename K, typename V>
void IOHDF5Map<K,V>::getKeys(std::set<K>& keys) const 
{
 keys.clear();
 check_for_failure(!ioh.isOpen(),"IOHDF5Map not open: cannot get keys");
 ioh.cd("Values");
 std::set<std::string> keynames(ioh.getDataNamesInCurrentDir());
 ioh.cd("..");
 for (std::set<std::string>::iterator it=keynames.begin();it!=keynames.end();++it)
    keys.insert(K(*it));
}


template<typename K, typename V>
void IOHDF5Map<K,V>::getKeys(std::vector<K>& keys) const
{
 keys.clear();
 check_for_failure(!ioh.isOpen(),"IOHDF5Map not open: cannot get keys");
 std::set<K> keyset;
 getKeys(keyset);
 std::copy(keyset.begin(),keyset.end(),std::back_inserter(keys));
}

            // keep only those keys in "keys_to_keep"
            // return true of all keys in "keys_to_keep" are available,
            // false otherwise

template<typename K, typename V>
bool IOHDF5Map<K,V>::keepKeys(const std::set<K>& keys_to_keep)  
{
 std::set<K> keys_avail;
 getKeys(keys_avail);
 for (typename std::set<K>::const_iterator it=keys_to_keep.begin();it!=keys_to_keep.end();it++){
    typename std::set<K>::iterator kt=keys_avail.find(*it);
    if (kt==keys_avail.end()) return false;} 
 return true;
}



template <typename K, typename V>
std::set<std::string> IOHDF5Map<K,V>::getAvailableRootPaths() const
{
 check_for_failure(!ioh.isOpen(),"IOHDF5Map not open: cannot get available root paths");
 std::set<std::string> alldirs(ioh.getAllDirNames());
 std::set<std::string> rootpaths;
 for (std::set<std::string>::const_iterator it=alldirs.begin();it!=alldirs.end();++it){
    size_t pos=it->find("Values");
    if (pos!=std::string::npos){
       rootpaths.insert(it->substr(0,pos));}}
 return rootpaths;
}



template<typename K, typename V>
void IOHDF5Map<K,V>::check_for_failure(bool errcond, const std::string& mesg,
                                       bool abort) const
{
 if (!errcond) return;
 std::cerr << "IOHDF5Map error: "<<mesg<<std::endl;
 if (abort) exit(1);
 else throw(std::invalid_argument(mesg));
}


// **************************************************************
#endif
