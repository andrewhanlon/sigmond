#ifndef IO_MAP_BASE_H
#define IO_MAP_BASE_H

#include <set>
#include <string>
#include <vector>


 // ********************************************************
 // *                                                      *
 // *      Purely virtual base class: the IOMapBase        *
 // *                                                      *
 // ********************************************************


template<typename K, typename V>
class IOMapBase
{

  public:
 
 
    IOMapBase() {}

            // read only open, returns header string
            
    virtual void openReadOnly(const std::string& filename, 
                              const std::string& filetype_id,
                              std::string& header, bool turn_on_checksum=false) = 0;

            // read only open, ignores header string

    virtual void openReadOnly(const std::string& filename, 
                              const std::string& filetype_id,
                              bool turn_on_checksum=false) = 0;

            // open a new file in read/write mode, writes the header string (fails 
            // if the file exists and "fail_if_exists" is true; if "fail_if_exists"
            // is false, deletes the existing file to start a new file)

    virtual void openNew(const std::string& filename, 
                         const std::string& filetype_id, 
                         const std::string& header,  
                         bool fail_if_exists=true, char endianness='N',
                         bool turn_on_checksum=false, bool overwrites_allowed=false) = 0;

            // open a file in read/write mode; if file exists, the header
            // string is read and returned in "header" and writes will update
            // the existing file; otherwise, a new file is created (in which 
            // case, the header string is needed as input so it can be written
            // into the new file)

    virtual void openUpdate(const std::string& filename, 
                            const std::string& filetype_id, 
                            std::string& header, 
                            char endianness='N', bool turn_on_checksum=false, 
                            bool overwrites_allowed=false) = 0;

    virtual ~IOMapBase() {}

    virtual void close() = 0;



    virtual std::string getHeader() = 0;  // file must be open
    
        // Version that assumes file is not open; file closed afterwards.
        // Returns false if file cannot be opened.
        
    virtual bool peekHeader(std::string& header, const std::string& filename, 
                            const std::string& filetype_id) = 0;
    
    virtual std::string getFileName() const = 0;

    virtual bool isOpen() const = 0;

    virtual bool isNewFile() const = 0;

    virtual bool isOverwriteOn() const = 0;


    virtual bool areChecksumsInFile() const = 0;
      
    virtual bool isFileLittleEndian() const = 0;
   
    virtual bool isFileBigEndian() const = 0;



    virtual void setHighVerbosity() = 0;

    virtual void setMediumVerbosity() = 0;
   
    virtual void setNoVerbosity() = 0;

    virtual void setDisallowOverwrites() = 0;

    virtual void setAllowOverwrites() = 0;



    virtual void put(const K& key, const V& val) = 0;
    
    virtual void get(const K& key, V& val) = 0;  // throws exception or aborts if fails

    virtual bool get_maybe(const K& key, V& val) = 0;  // returns false is fails, true otherwise

    virtual bool exist(const K& key) const = 0;
    
    virtual void flush() = 0;  // puts file in finalized state so no data loss if abort occurs
     


    virtual unsigned int size() const = 0;

    virtual void getKeys(std::vector<K>& keys) const = 0;

    virtual void getKeys(std::set<K>& keys) const = 0;

    virtual bool keepKeys(const std::set<K>& keys_to_keep) = 0;   // keep only those keys in "keys_to_keep"
                          // return true if all keys in "keys_to_keep" are available,
                          // false otherwise


};

#endif
