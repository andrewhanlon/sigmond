#ifndef DATA_IO_HANDLER_H
#define DATA_IO_HANDLER_H

#include "io_map.h"
#include "xml_handler.h"
#include <algorithm>
#include <set>

 // *********************************************************************************
 // *                                                                               *
 // *   The classes defined in this file that are important are                     *
 // *                                                                               *
 // *         DataGetHandlerSF   (single file)                                      *
 // *         DataPutHandlerSF   (single file)                                      *
 // *         DataGetHandlerMF   (multi-file)                                       *
 // *                                                                               *
 // *                                                                               *
 // *   Theses classes handle inserting and retrieving data to/from files via       *
 // *   IOMap objects.  Use "put" to build up the data in files, and "get"          *
 // *   to access the data in the files.  These classes do NOT store the data       *
 // *   in memory; they just handle reading and writing to file.                    *
 // *                                                                               *
 // *   Objects of these "put" handlers always assume an "updating" mode. Existing  *
 // *   files are never erased, and new files are created as needed.  New records   *
 // *   are added to the files.  If the key of a record to be put already exists    *
 // *   in a file, the put will only occur if "overwrite" is specified AND the      *
 // *   size of the data to be put does not exceed the size of the data already in  *
 // *   the file for that key.                                                      *
 // *                                                                               *
 // *   To use the single-file classes, one needs the following ingredients:        *
 // *                                                                               *
 // *    (a) a file name (string) and an IOHandler file ID string                   *
 // *                                                                               *
 // *    (b) a handler class "H" that has members                                   *
 // *                                                                               *
 // *            bool checkHeader(XMLHandler& xmlin)                                *
 // *                                                                               *
 // *          that checks that the header information in the file is what it       *
 // *          should be for class "H", returning a boolean value, and              *
 // *                                                                               *
 // *            void writeHeader(XMLHandler& xmlout)                               *
 // *                                                                               *
 // *        that writes out the header string to xmlout                            *
 // *                                                                               *
 // *    (c) the record key class must have all of the features of an IOMap key:    *
 // *                                                                               *
 // *       -- since used in a C++ map, a less than operator must be define         *
 // *              const K& K::operator<(const K& rhs);                             *
 // *       -- a numbytes(ioh,K) function, where ioh is an IOHandler, must be       *
 // *           defined to give the number of bytes each key occupies in an         *
 // *           IOHandler file                                                      *
 // *       -- a copy constructor K(const K& in) must be defined                    *
 // *           (a default constructor is not needed)                               *
 // *       -- a multi_read(ioh, vector<K>&,n) must be defined to read n keys       *
 // *       -- a multi_write(ioh, const vector<K>&) must be defined                 *
 // *       -- an "output(XMLHander&)" output member                                *
 // *                                                                               *
 // *    (d) the value type must have the following features:                       *
 // *                                                                               *
 // *       -- a write(ioh, const V&) must be defined (ioh is an IOHandler object)  *
 // *       -- a read(ioh, V&) must be defined                                      *
 // *       -- a numbytes(ioh,V) must be defined giving number of bytes occupied    *
 // *            by V in an IOHandler file                                          *
 // *                                                                               *
 // *   The class "DataGetHandlerMF" extends the get to multiple files.  A linear   *
 // *   search for a key through all files is done, so this is not efficient for    *
 // *   a large number of files.                                                    *
 // *                                                                               *
 // *   The member "keepKeys" is used to limit attention to a subset of keys.       *
 // *   It is useful if only a few keys are needed, so only those keys are          *
 // *   kept in memory.  The "keepKeys" members returns true if all requested       *
 // *   keys are available, false if some are missing (not available).              *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************
 
 
   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerSF                      *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataGetHandlerSF
{

    IOMap<R,D> *iomptr;
    H& handler;
    
 public:

    DataGetHandlerSF(H& in_handler, const std::string& file_name, 
                     const std::string& filetype_id,
                     bool use_checksums=false);

    ~DataGetHandlerSF() {delete iomptr;}

    void close() { iomptr->close();}

    bool keepKeys(const std::set<R>& keys_to_keep);

    std::string getFileName() const {return iomptr->getFileName();}


    bool queryData(const R& rkey);

    void getData(const R& rkey, D& result);

    bool getDataMaybe(const R& rkey, D& result);


    std::set<R> getKeys();

    void outputKeys(XMLHandler& xmlout);
    
    unsigned int size() const {return iomptr->size();}



 private:

    void fail(const R& rkey);

    void fail(const std::string& msg);

          // disallow copies
    DataGetHandlerSF(const DataGetHandlerSF& in);
    DataGetHandlerSF& operator=(const DataGetHandlerSF& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataGetHandlerSF<H,R,D>::DataGetHandlerSF(H& in_handler, const std::string& filename,
                                          const std::string& filetype_id, 
                                          bool use_checksums)
                       :  handler(in_handler)
{
 std::string headerxml; 
 {IOMap<R,D> iom;
 bool exists=iom.peekHeader(headerxml,filename,filetype_id);
 if (!exists) throw(std::invalid_argument(std::string("could not open file ")+filename
                  +std::string(" for reading")));}

 XMLHandler xmlr; xmlr.set_from_string(headerxml);
 if (!handler.checkHeader(xmlr)){
    throw(std::invalid_argument(std::string("Header string in file is \n")+headerxml
         +std::string("\n header info in file ")+filename
         +std::string(" does not match info in current Handler\n ...execution aborted...\n")));}
 try{
    iomptr=new IOMap<R,D>; 
    iomptr->openReadOnly(filename,filetype_id,use_checksums);}
 catch(const std::exception& xp) {
    fail("could not open file "+filename+" for reading");}
}


template <typename H, typename R, typename D>
bool DataGetHandlerSF<H,R,D>::keepKeys(const std::set<R>& keys_to_keep)
{
 return iomptr->keepKeys(keys_to_keep);
} 


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::getData(const R& rkey, D& result)
{
 try {iomptr->get(rkey,result);}
 catch(const std::exception& xp){fail(rkey);}
}


template <typename H, typename R, typename D>
bool DataGetHandlerSF<H,R,D>::getDataMaybe(const R& rkey, D& result)
{
 return iomptr->get_maybe(rkey,result);
}


template <typename H, typename R, typename D>
bool DataGetHandlerSF<H,R,D>::queryData(const R& rkey)
{
 return iomptr->exist(rkey);
}


template <typename H, typename R, typename D>
std::set<R> DataGetHandlerSF<H,R,D>::getKeys()
{
 std::set<R> keys;
 iomptr->getKeys(keys);
 return keys;
} 


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::outputKeys(XMLHandler& xmlout)
{
 xmlout.set_root("AvailableKeys");
 std::set<R> keys;
 iomptr->getKeys(keys);
 for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
    XMLHandler xmlk;
    kt->output(xmlk);
    XMLHandler xmlt("Key");
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);}
}



template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::fail(const R& rkey)
{
 std::cout << "DataGetHandlerSF could not find requested record:"<<std::endl;
 std::cout << " File name: "<< iomptr->getFileName() <<std::endl;
 XMLHandler xmlout;
 xmlout.set_root("RecordKey");
 XMLHandler xmlh;
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 delete iomptr;
 exit(1);
}



template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::fail(const std::string& msg)
{
 delete iomptr; iomptr=0;
 throw(std::invalid_argument(std::string("DataGetHandlerSF error: ")+msg));
}



   // **************************************************************
   // *                                                            *
   // *                      DataPutHandlerSF                      *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataPutHandlerSF
{

    H& handler;
    IOMap<R,D> *iomptr;


 public:

    DataPutHandlerSF(H& inptr, const std::string& file_name, 
                     const std::string& filetype_id,
                     bool overwrite=false, bool use_checksums=false);

    ~DataPutHandlerSF() {delete iomptr;}

    std::string getFileName() const {return iomptr->getFileName();}


    void putData(const R& rkey, const D& data);

    void flush();                   // writes out file map at end of file

    bool queryData(const R& rkey);  // already exists?

    void close() { iomptr->close();}

 private:
 
    void fail(const std::string& msg);

    void fail(const R& rkey);


          // disallow copies and default
    DataPutHandlerSF();
    DataPutHandlerSF(const DataPutHandlerSF& in);
    DataPutHandlerSF& operator=(const DataPutHandlerSF& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataPutHandlerSF<H,R,D>::DataPutHandlerSF(H& in_handler,
                                          const std::string& file_name,
                                          const std::string& filetype_id,
                                          bool overwrite, bool use_checksums)
                     :  handler(in_handler)
{
 XMLHandler headerxml;
 handler.writeHeader(headerxml);  // write header info 
 try{
   iomptr=new IOMap<R,D>;
   std::string header(headerxml.str());
   iomptr->openUpdate(file_name,filetype_id,header,'L',
                      use_checksums,overwrite);
   if (!(iomptr->isNewFile())){
      XMLHandler xmlr; xmlr.set_from_string(header);
      if (!handler.checkHeader(xmlr)){
         fail("Header string in file is \n"+header+"\n header info in file "
              +file_name+" does not match info in current Handler\n ");}}}
 catch(const std::exception& xp) {
    fail("could not open file "+file_name+" for writing"); }
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::putData(const R& rkey, const D& data)
{
 try{
    iomptr->put(rkey,data);}
 catch(const std::exception& xp){
    fail(rkey);}
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::flush()
{
 iomptr->flush();
}


template <typename H, typename R, typename D>
bool DataPutHandlerSF<H,R,D>::queryData(const R& rkey)
{
 return iomptr->exist(rkey);
}



template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::fail(const std::string& msg)
{
// std::cout << "DataPutHandlerSF error: "<<msg<<std::endl;
// delete iomptr;
// exit(1);
 throw(std::runtime_error(std::string("DataPutHandlerSF error: ")+msg));
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::fail(const R& rkey)
{
/* std::cout << "DataPutHandlerSF could not insert requested record:"<<std::endl;
 std::cout << " File name: "<< iomptr->getFileName() <<std::endl;
 XMLHandler xmlout;
 xmlout.set_root("RecordKey");
 XMLHandler xmlh;
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 delete iomptr;
 exit(1); */
 XMLHandler xmlkey; rkey.output(xmlkey);
 throw(std::runtime_error(std::string("DataPutHandlerSF could not insert requested record: ")
         +xmlkey.str()));
}



   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerMF                      *
   // *                                                            *
   // **************************************************************




template <typename H, typename R, typename D>
class DataGetHandlerMF
{

    std::list<DataGetHandlerSF<H,R,D>* > getptrs;
    H& m_handler;
    std::string m_filetype_id;
    bool m_use_checksums;

 public:

    DataGetHandlerMF(H& in_handler, const std::set<std::string>& file_names, 
                     const std::string& filetype_id,
                     bool use_checksums=false);

    ~DataGetHandlerMF();

    void addFile(const std::string& file_name);

    bool addFile(const std::string& file_name, const std::set<R>& keys_to_keep);

    void removeFile(const std::string& file_name);

    bool keepKeys(const std::set<R>& keys_to_keep);

    void clear();
 
    void close();

    std::set<std::string> getFileNames() const;


    bool queryData(const R& rkey);

    void getData(const R& rkey, D& result);

    bool getDataMaybe(const R& rkey, D& result);



    std::set<R> getKeys();

    void outputKeys(XMLHandler& xmlout);
    
    void getFileMap(XMLHandler& xmlout);

    unsigned int size() const;


 private:

          // disallow copies
    DataGetHandlerMF(const DataGetHandlerMF& in);
    DataGetHandlerMF& operator=(const DataGetHandlerMF& in);

};



template <typename H, typename R, typename D>
DataGetHandlerMF<H,R,D>::DataGetHandlerMF(H& in_handler, const std::set<std::string>& file_names, 
                                          const std::string& filetype_id, 
                                          bool use_checksums)
                       :  m_handler(in_handler), m_filetype_id(filetype_id),
                          m_use_checksums(use_checksums)
{
 for (std::set<std::string>::const_iterator it=file_names.begin();it!=file_names.end();it++)
    addFile(*it);
}


template <typename H, typename R, typename D>
DataGetHandlerMF<H,R,D>::~DataGetHandlerMF()
{
 clear();
}


template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::addFile(const std::string& file_name)
{
 std::set<R> keys_to_keep;
 addFile(file_name,keys_to_keep);
}


    // only the keys in "keys_to_keep" will be added (if not empty)

template <typename H, typename R, typename D>
bool DataGetHandlerMF<H,R,D>::addFile(const std::string& file_name, const std::set<R>& keys_to_keep)
{
 bool flag=true;
 try{
    // first check that file_name not already open for getting
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    if (tidyString(file_name)==(*it)->getFileName()){
       throw(std::invalid_argument(std::string("cannot addFile ")+file_name
               +std::string("in DataGetHandlerMF since already open ")));}
    // create the single-file getter
 DataGetHandlerSF<H,R,D>* newget=0;
 try{
    newget=new DataGetHandlerSF<H,R,D>(
                 m_handler,file_name,m_filetype_id,m_use_checksums);
    if (!keys_to_keep.empty())
       flag=newget->keepKeys(keys_to_keep);}
 catch(const std::exception& xp){
    throw(std::invalid_argument(std::string("cannot addFile ")+file_name
          +std::string("in DataGetHandlerMF since single-file open failed ")));}
 std::set<R> newkeys(newget->getKeys());
    // check that all keys in new file are different from all keys already available
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++){
    std::set<R> rkeys((*it)->getKeys());
    std::set<R> intersect;
    std::set_intersection(newkeys.begin(),newkeys.end(),rkeys.begin(),rkeys.end(),
                          std::inserter(intersect,intersect.begin()));
    if (!intersect.empty()){
       delete newget;
       throw(std::invalid_argument(std::string("cannot addFile ")+file_name
              +std::string("in DataGetHandlerMF since duplicate keys ")));}}
 getptrs.push_back(newget);}
 catch(const std::exception& xp){
    std::cout << "Fatal error: addFile failed in DataGetHandlerMF: "<<xp.what()<<std::endl;
    exit(1);}
 return flag;
}


template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::removeFile(const std::string& file_name)
{
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    if ((*it)->getFileName()==tidyString(file_name)){
       delete *it;
       getptrs.erase(it);
       return;}
}



template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::clear()
{
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    delete *it;
 getptrs.clear();
}


template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::close()
{
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    (*it)->close();
}


template <typename H, typename R, typename D>
bool DataGetHandlerMF<H,R,D>::keepKeys(const std::set<R>& keys_to_keep)
{
 uint ksize=0;
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++){
    (*it)->keepKeys(keys_to_keep);
    ksize+=(*it)->size();}
 return (ksize==keys_to_keep.size());
} 


template <typename H, typename R, typename D>
std::set<std::string> DataGetHandlerMF<H,R,D>::getFileNames() const
{
 std::set<std::string> fnames;
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::const_iterator it=getptrs.begin();it!=getptrs.end();it++)
    fnames.insert((*it)->getFileName());
 return fnames;
}


template <typename H, typename R, typename D>
bool DataGetHandlerMF<H,R,D>::queryData(const R& rkey)
{
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    if ((*it)->queryData(rkey)) return true;
 return false;
}


template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::getData(const R& rkey, D& result)
{
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    if ((*it)->getDataMaybe(rkey,result)) return;
 XMLHandler xmlk; rkey.output(xmlk);
 throw(std::runtime_error(std::string("Could not find key ")+xmlk.output()
         +std::string(" in DataGetHandlerMF")));
}



template <typename H, typename R, typename D>
bool DataGetHandlerMF<H,R,D>::getDataMaybe(const R& rkey, D& result)
{
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++)
    if ((*it)->getDataMaybe(rkey,result)) return true;
 return false;
}



template <typename H, typename R, typename D>
std::set<R> DataGetHandlerMF<H,R,D>::getKeys()
{
 std::set<R> keys;
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++){
    std::set<R> nextkeys((*it)->getKeys());
    keys.insert(nextkeys.begin(),nextkeys.end());}
 return keys;
}


template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::outputKeys(XMLHandler& xmlout)
{
 std::set<R> keys(getKeys());
 xmlout.set_root("AvailableKeys");
 for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
    XMLHandler xmlk;
    kt->output(xmlk);
    XMLHandler xmlt("Key");
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);}
}


template <typename H, typename R, typename D>
void DataGetHandlerMF<H,R,D>::getFileMap(XMLHandler& xmlout)
{
 xmlout.set_root("FileMap");
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::iterator it=getptrs.begin();it!=getptrs.end();it++){
    std::set<R> keys((*it)->getKeys());
    std::string fname((*it)->getFileName());
    for (typename std::set<R>::iterator et=keys.begin();et!=keys.end();et++){
       XMLHandler xmle("Entry");
       XMLHandler xmlo; et->output(xmlo);
       xmle.put_child(xmlo);
       xmle.put_child("Name",fname);
       xmlout.put_child(xmle);}}
}


template <typename H, typename R, typename D>
unsigned int DataGetHandlerMF<H,R,D>::size() const
{
 unsigned int tsize=0;
 for (typename std::list<DataGetHandlerSF<H,R,D>* >::const_iterator it=getptrs.begin();it!=getptrs.end();it++)
    tsize+=(*it)->size();
 return tsize;
}


// **************************************************************
#endif
