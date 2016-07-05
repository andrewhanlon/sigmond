#ifndef LAPH_DATA_IO_HANDLER_H
#define LAPH_DATA_IO_HANDLER_H

#include "io_map.h"
#include "filelist_info.h"
#include "xml_handler.h"
#include <set>
#include <list>
#include <sstream>

namespace LaphEnv {


 // *********************************************************************************
 // *                                                                               *
 // *   The class defined in this file is                                           *
 // *                                                                               *
 // *         DataGetHandlerMF   (multi-file)                                       *
 // *                                                                               *
 // *   This is a very important mid-level class for getting access to the LapH     *
 // *   data to analyze.  It utilizes IOMap objects (which utilizes IOHandler       *
 // *   objects) to keep track of a large set of files containing data. Objects of  *
 // *   this class are meant to be used by higher-level "CorrelatorDataHandler"     *
 // *   and "VEVDataHandler" objects.  Objects of this class just handle reading    *
 // *   reading the data from files and do NOT storage the results in memory        *
 // *   in any way.  This class is templated on a file key type, and record key     *
 // *   type, and a data type:                                                      *
 // *                                                                               *
 // *    (a) this class extracts the file key from the header string, so            *
 // *          the file key type must have a constructor that takes only an         *
 // *          XMLHandler and an output(XMLHandler&) member; you need to specify the*
 // *          header tag that the file key will be enclosed in; the file key must  *
 // *          also have a "<" operator and an "==" operator defined                *
 // *                                                                               *
 // *    (b) the record key class must have all of the features of an IOMap key:    *
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
 // *                                                                               *
 // *    (c) the value type must have the following features:                       *
 // *                                                                               *
 // *       -- a write(ioh, const V&) must be defined (ioh is an IOHandler object)  *
 // *       -- a read(ioh, V&) must be defined                                      *
 // *       -- a numbytes(ioh,V) must be defined giving number of bytes occupied    *
 // *            by V in an IOHandler file                                          *
 // *                                                                               *
 // *                                                                               *
 // *   Since this class is meant to be used by "CorrelatorDataHandler" and         *
 // *   "VEVDataHandler" objects, its constructor is somewhat awkward.  These       *
 // *   higher-level classes have the responsibility of constructing the            *
 // *   objects used in the constructor below.  The constructor takes the form:     *
 // *                                                                               *
 // *      DataGetHandlerMF<FileKey,RecordKey,DataType>(                            *
 // *                const std::vector<std::string>& stubs,                         *
 // *                const std::map<F,std::pair<int,int> >& data_files,             *
 // *                const std::string& filetype_id,                                *
 // *                unsigned int in_max_open_files,                                *
 // *                double cleanfrac=0.25,                                         *
 // *                bool use_checksums=false);                                     *
 // *                                                                               *
 // *   "stubs" is a vector containing the FileListInfo stubs for all of the        *
 // *   files to access.  "data_files" contains a map that associates a given       *
 // *   FileKey with a file specified by a pair of integers: the first integer      *
 // *   is the index of the stub as specified in "stubs" and the second integer     *
 // *   is the suffix index.  "filetype_id" is the LapH string expected at the      *
 // *   start of each file, and "use_checksums" specifies whether or not            *
 // *   checksum testing is performed.                                              *
 // *                                                                               *
 // *   The number of files to consider can be quite large (in the thousands),      *
 // *   so to avoid memory exhaustion, this class has a "maxopenfiles" member.      *
 // *   Once this many files are opened, further attempts to open new files         *
 // *   results in a garbage collection: some fraction of open files, specified     *
 // *   by "cleanfrac", are closed.  A list of opened files is maintained,          *
 // *   in order of their opening, and the fraction of open files at the            *
 // *   beginning of this list, presumably the oldest files, are closed.            *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************
 
 

   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which contains a pointer to a stub name, a suffix index, and an IOMap 
   // pointer.  The stub names are stored in "filestubs".  Upon construction, 
   // "fileMap" is assigned from an input map that associates file names
   // with file keys.  This input map is assumed to have been created
   // in such a way that all necessary checks have been done.  No
   // further checks are done by this class, to speed up access times.
   // Initially, no files are open.  As data is accessed, the files are 
   // opened and left open.  "queryData" and "getData" open files as needed,
   // and leave them open.

   // Keep in mind that while a file is open, the IOMap keeps all of the
   // record keys in memory.  Hence, memory exhaustion is an issue for a
   // large number of files.  To prevent memory exhaustion, the constructor
   // requires a maximum number of open files.  A list of open files is
   // maintained in "opened".




template <typename F, typename R, typename D>
class DataGetHandlerMF
{

    typedef std::set<std::string>::const_iterator    StubPtr;

    struct FileMapValue
    {
      StubPtr itstub;
      int suffix;
      IOMap<R,D> *fptr;

      FileMapValue(StubPtr instub, int insuffix) 
            : itstub(instub), suffix(insuffix), fptr(0) {}
      FileMapValue(const FileMapValue& fmv) 
            : itstub(fmv.itstub), suffix(fmv.suffix), fptr(fmv.fptr) {}
      ~FileMapValue() { delete fptr;}
    };

    typedef std::map<F,FileMapValue>   FileMapType;

    std::set<std::string> filestubs;
    FileMapType  fileMap;
    std::string fid;
    bool checksums;
    unsigned int maxopenfiles;
    unsigned int defaultclean;
    std::list<FileMapValue* > opened;


 public:

    DataGetHandlerMF(const std::vector<std::string>& stubs,
                     const std::map<F,std::pair<int,int> >& data_files,
                     const std::string& filetype_id, 
                     unsigned int in_max_open_files,
                     double cleanfrac=0.25,
                     bool use_checksums=false);

    ~DataGetHandlerMF() {}

    bool queryData(const F& fkey, const R& rkey);

    bool queryFile(const F& fkey);

    void getData(const F& fkey, const R& rkey, D& result);

    void close(const F& fkey);

    void close(double fraction);

    void close();


    void getFileMap(XMLHandler& xmlout) const;

    std::set<F> getFileKeys() const;

    std::set<R> getKeys(const F& fkey);

    void outputKeys(XMLHandler& xmlout);


 private:

    IOMap<R,D>* get_file_ptr(const F& fkey);
    
    void open(FileMapValue& fmv);
    
    void close(FileMapValue& fmv);

    void close(unsigned int number_to_close);

    std::string getFileName(const FileMapValue& fmv) const;

    unsigned int getNumberToClose(double cleanfrac, unsigned int total) const;

    void fail(const F& fkey, const R& rkey);
    
    void fail(const F& fkey);

    void fail(const std::string& msg);
    
          // disallow copies
    DataGetHandlerMF(const DataGetHandlerMF& in);
    DataGetHandlerMF& operator=(const DataGetHandlerMF& in);

};



   // Constructor sets the fileMap.

template <typename F, typename R, typename D>
DataGetHandlerMF<F,R,D>::DataGetHandlerMF(
                     const std::vector<std::string>& stubs,
                     const std::map<F,std::pair<int,int> >& in_files,
                     const std::string& filetype_id,
                     unsigned int in_max_open_files,
                     double cleanfrac, bool use_checksums)
           :  fid(tidyString(filetype_id)), checksums(use_checksums), 
              maxopenfiles(in_max_open_files), 
              defaultclean(getNumberToClose(cleanfrac,in_max_open_files))
{
 if ((maxopenfiles<4)||(maxopenfiles>8192))
    throw(std::invalid_argument("Max open files in DataGetHandlerMF is unreasonable"));
 for (typename std::map<F,std::pair<int,int> >::const_iterator 
       it=in_files.begin();it!=in_files.end();it++){
    int stubind=(it->second).first;
    int suffix=(it->second).second;
    if ((stubind<0)||(stubind>=int(stubs.size())))
       throw(std::invalid_argument("Improper stub index in DataGetHandlerMF"));
    StubPtr itstub=filestubs.insert(stubs[stubind]).first;
    fileMap.insert(std::make_pair(it->first, FileMapValue(itstub,suffix)));}
}



template <typename F, typename R, typename D>
bool DataGetHandlerMF<F,R,D>::queryData(const F& fkey, const R& rkey)
{
 IOMap<R,D> *fptr=get_file_ptr(fkey);
 if (fptr==0) return false;
 return fptr->exist(rkey);
}


template <typename F, typename R, typename D>
bool DataGetHandlerMF<F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 return (it!=fileMap.end());
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::getData(const F& fkey, const R& rkey,
                                      D& result)
{
 IOMap<R,D> *fptr=get_file_ptr(fkey);
 if (fptr==0) fail(fkey);
 try {fptr->get(rkey,result);}
 catch(...){fail(fkey,rkey);}
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::close(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it!=fileMap.end()){
    close(it->second);}
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::close(double fraction)
{
 unsigned int num=getNumberToClose(fraction,opened.size());
 close(num);
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::close()
{
 for (typename FileMapType::iterator it=fileMap.begin();it!=fileMap.end();it++){
    delete it->second.fptr; it->second.fptr=0;}
 opened.clear();
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::getFileMap(XMLHandler& xmlout) const
{
 xmlout.set_root("FileMap");
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    XMLHandler xmlt("Entry");
    XMLHandler xmlp;
    it->first.output(xmlp); xmlt.put_child(xmlp);
    xmlt.put_child("Name",getFileName(it->second));
    xmlout.put_child(xmlt);}
}

template <typename F, typename R, typename D>
std::set<F> DataGetHandlerMF<F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 


template <typename F, typename R, typename D>
std::set<R> DataGetHandlerMF<F,R,D>::getKeys(const F& fkey)
{
 std::set<R> keys;
 IOMap<R,D>* fptr=get_file_ptr(fkey);
 if (fptr!=0) fptr->getKeys(keys);
 return keys;
}
 


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::outputKeys(XMLHandler& xmlout)
{
 xmlout.set_root("AvailableKeys");
 for (typename FileMapType::iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    std::set<R> keys;
    if (it->second.fptr!=0) 
       it->second.fptr->getKeys(keys);
    else{
       open(it->second);
       it->second.fptr->getKeys(keys);
       close(it->second);} 
    for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
       XMLHandler xmli,xmlk;
       it->first.output(xmli);
       kt->output(xmlk);
       XMLHandler xmlt("Key");
       xmlt.put_child(xmli);
       xmlt.put_child(xmlk);
       xmlout.put_child(xmlt);}}
}


                //   private members



template <typename F, typename R, typename D>
IOMap<R,D>* DataGetHandlerMF<F,R,D>::get_file_ptr(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it==fileMap.end()) return 0;
 if (it->second.fptr==0) open(it->second);
 return it->second.fptr;
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::open(FileMapValue& fmv)
{
 if (opened.size()>=maxopenfiles){
    close(defaultclean);}
 try {
    fmv.fptr=new IOMap<R,D>;
    fmv.fptr->openReadOnly(getFileName(fmv),fid,checksums);
    opened.push_back(&fmv);}
 catch(...){
    fail("failure opening file "+getFileName(fmv)+" in DataGetHandlerMF");}
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::close(FileMapValue& fmv)
{
 delete fmv.fptr; 
 fmv.fptr=0;
 for (typename std::list<DataGetHandlerMF<F,R,D>::FileMapValue* >::iterator 
      it=opened.begin();it!=opened.end();it++){
    if (*it==&fmv){
       opened.erase(it);
       return;}}
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::close(unsigned int number_to_close)
{
 unsigned int n=number_to_close;
 if (opened.size()<n) n=opened.size();
 for (unsigned int k=0;k<n;k++){
    delete opened.front()->fptr;
    opened.front()->fptr=0;
    opened.pop_front();}
}


template <typename F, typename R, typename D>
std::string DataGetHandlerMF<F,R,D>::getFileName(const FileMapValue& fmv) const
{
 std::stringstream fs;
 fs << *(fmv.itstub) << "." << fmv.suffix;
 return fs.str();
}


template <typename F, typename R, typename D>
unsigned int DataGetHandlerMF<F,R,D>::getNumberToClose(
                   double cleanfrac, unsigned int total) const
{
 if (cleanfrac<=0.0) return 1;
 else if (cleanfrac>=1.0) return total;
 unsigned int toclean=floor(cleanfrac*total);
 if (toclean<1) toclean=1;
 return toclean;
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::fail(const F& fkey, const R& rkey)
{
 XMLHandler xmlout;
 xmlout.set_root("FileRecordKey");
 XMLHandler xmlh;
 fkey.output(xmlh);
 xmlout.put_child(xmlh);
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 throw(std::invalid_argument((std::string("Could not find requested record")+xmlout.str()).c_str()));
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::fail(const F& fkey)
{
 XMLHandler xmlout;
 fkey.output(xmlout);
 throw(std::invalid_argument((std::string("Could not find requested file key")+xmlout.str()).c_str()));
}


template <typename F, typename R, typename D>
void DataGetHandlerMF<F,R,D>::fail(const std::string& msg)
{
 throw(std::invalid_argument((std::string("DataGetHandlerMF error: ")+msg).c_str()));
}


// **************************************************************
}
#endif
