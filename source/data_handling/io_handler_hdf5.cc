#include "io_handler_hdf5.h"
#include <unistd.h>
using namespace std;

// *************************************************************************

const int IOHDF5Handler::IO_ERR_NO_SUCH_FILE= 1;
const int IOHDF5Handler::IO_ERR_ACCESS=       1;
const int IOHDF5Handler::IO_ERR_OTHER=        1;
const int IOHDF5Handler::IO_SUCCESS=          0;

const IOHDF5Handler::openmode_type IOHDF5Handler::IO_MODE_RDONLY= ios::in;
const IOHDF5Handler::openmode_type IOHDF5Handler::IO_MODE_RDWR=   ios::in | ios::out;
const IOHDF5Handler::openmode_type IOHDF5Handler::IO_MODE_CREATE= ios::trunc;
const IOHDF5Handler::openmode_type IOHDF5Handler::IO_MODE_EXCL=   ios::trunc;

ByteHandler IOHDF5Handler::m_bytehandler;


// *************************************************************************


IOHDF5Handler::IOHDF5Handler() : fid(-1), is_new_file(false), read_only(true), openflag(false), 
                                 read_mode(true), endian_format('U'), endian_convert(false), 
                                 checksum_on(false), currid(-1), checksum(0)
{}


IOHDF5Handler::IOHDF5Handler(const std::string& filename, OpenMode mode,
                             const std::string& filetype_id, char endianness,
                             bool turn_on_checksum)
{
 openflag=false;
 open(filename,mode,filetype_id,endianness,turn_on_checksum);
}

IOHDF5Handler::IOHDF5Handler(const std::string& filename, ios_base::openmode mode,
                             const std::string& filetype_id, char endianness,
                             bool turn_on_checksum)
{
 openflag=false;
 open(filename,mode,filetype_id,endianness,turn_on_checksum);
}


void IOHDF5Handler::open(const std::string& filename, IOHDF5Handler::OpenMode mode,
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
    is_new_file=false; 
    open_existing_file(filetype_id);}
 else{
    is_new_file=true;
    open_new_file(filetype_id,endianness);}

 checksum_on=turn_on_checksum;
 checksum=0;
 read_mode=true;
 currid=H5Gopen(fid,"/",H5P_DEFAULT);
 check_for_hid_failure(currid,"Could not open root group");
 assign_dtypes();
}


void IOHDF5Handler::assign_dtypes()
{
 if (endian_format=='L'){ 
    dtype_int=H5T_STD_I32LE;
    if (sizeof(long)==4){
       dtype_long=H5T_STD_I32LE;
       dtype_ulong=H5T_STD_U32LE;}
    else if (sizeof(long)==8){
       dtype_long=H5T_STD_I64LE;
       dtype_ulong=H5T_STD_U64LE;}
    else{
       check_for_failure(IO_ERR_OTHER,"Unsupported long size");}
    dtype_llong=H5T_STD_I64LE;
    dtype_uint=H5T_STD_U32LE;
    dtype_ullong=H5T_STD_U64LE;
    dtype_float=H5T_IEEE_F32LE;
    dtype_double=H5T_IEEE_F64LE;}
 else if (endian_format=='B'){
    dtype_int=H5T_STD_I32BE;
    if (sizeof(long)==4){
       dtype_long=H5T_STD_I32BE;
       dtype_ulong=H5T_STD_U32BE;}
    else if (sizeof(long)==8){
       dtype_long=H5T_STD_I64BE;
       dtype_ulong=H5T_STD_U64BE;}
    else{
       check_for_failure(IO_ERR_OTHER,"Unsupported long size");}
    dtype_llong=H5T_STD_I64BE;
    dtype_uint=H5T_STD_U32BE;
    dtype_ullong=H5T_STD_U64BE;
    dtype_float=H5T_IEEE_F32BE;
    dtype_double=H5T_IEEE_F64BE;}
 else{
    check_for_failure(IO_ERR_OTHER,"Unsupported endian format");}
}



void IOHDF5Handler::open(const std::string& filename, ios_base::openmode mode,
                         const std::string& filetype_id, char endianness,
                         bool turn_on_checksum)
{
 IOHDF5Handler::OpenMode iomode;
 if (!(mode & ios_base::out)) iomode=ReadOnly;
 else{
    if (mode & ios_base::trunc) iomode=ReadWriteEraseIfExists;
    else iomode=ReadWriteFailIfExists;}
 open(filename,iomode,filetype_id,endianness,turn_on_checksum);
}


void IOHDF5Handler::openReadOnly(const std::string& filename, 
                                 const std::string& filetype_id,
                                 bool turn_on_checksum)
{
 open(filename,ReadOnly,filetype_id,'N',turn_on_checksum);
}


void IOHDF5Handler::openNew(const std::string& filename, bool fail_if_exists,
                            const std::string& filetype_id, char endianness,
                            bool turn_on_checksum)
{
 OpenMode iomode=(fail_if_exists) ? ReadWriteFailIfExists : ReadWriteEraseIfExists;
 open(filename,iomode,filetype_id,endianness,turn_on_checksum);
}


void IOHDF5Handler::openUpdate(const std::string& filename,
                               const std::string& filetype_id, char endianness,
                               bool turn_on_checksum)
{
 open(filename,ReadWriteUpdateIfExists,filetype_id,endianness,turn_on_checksum);
}



void IOHDF5Handler::open_existing_file(const std::string& filetype_id)
{
 if (read_only){
    fid=H5Fopen(m_filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    check_for_hid_failure(fid,"Could not open file for reading only");}
 else{
    fid=H5Fopen(m_filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    check_for_hid_failure(fid,"Could not open file for update reading/writing");}

     // get endian info, and check file id
 openflag=true;
 std::string ID_string;
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
 if (m_bytehandler.big_endian())
    endian_convert=(endian_format=='L')?true:false;
 else 
    endian_convert=(endian_format=='B')?true:false;
}


void IOHDF5Handler::open_new_file(const std::string& filetype_id,
                                  char endianness)
{
 if (endianness=='N'){  // native 
    endian_format=m_bytehandler.big_endian() ? 'B':'L';
    endian_convert=false;}
 else if ((endianness!='B')&&(endianness!='L'))
    check_for_failure(IO_ERR_OTHER,"Invalid endian format");
 else{
    endian_format=endianness;
    if (m_bytehandler.big_endian()) 
       endian_convert=(endian_format=='L')?true:false;
    else 
       endian_convert=(endian_format=='B')?true:false;}
 fid=H5Fcreate(m_filename.c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT);
 check_for_hid_failure(fid,"Could not open file for reading/writing if nonexisting");
 openflag=true;
 writeIDstring(filetype_id);
}

       // Destructor

IOHDF5Handler::~IOHDF5Handler() 
{
 clear(); 
}


//  Creates a new "directory" (HDF5 group); does NOT behave the same way as 
//  linux mkdir --- does nothing if already exists, creates all groups above if
//  they do not already exist.   "path" can be absolute or relative
//  Current working directory does NOT change.

void IOHDF5Handler::mkdir(const std::string& path)
{
 if (!openflag) return;
 if (queryDir(path)) return;
 std::string ppath(tidyString(path));
 hid_t* cwdptr=(ppath[0]=='/')? (&fid) : (&currid);
 vector<std::string> tokens(tokenize(path));  // performs some checks
 if (tokens.size()==0) return;
 if (ppath[0]=='/') tokens[0]="/"+tokens[0];
 for (uint k=1;k<tokens.size();++k)
    tokens[k]="/"+tokens[k];
 ppath="";
 for (uint i=0;i<tokens.size();++i){
    ppath+=tokens[i];
    if (!(queryDir(ppath))){
       hid_t gid = H5Gcreate2(*cwdptr,ppath.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
       check_for_hid_failure2(gid,"mkdir failed for path ",ppath);
       H5Gclose(gid);}}
}

//  Changes the current working directory; currid is both input and output
//  a path "." does nothing, ".." goes up one level (unless at root)
//  path can be absolute or relative; fid and filename refer to the file;
//  currid is current group, which will be changed; current path is also
//  specified by "currwd" as a vector of strings which include no "/".
//  If "mkdir" is true, it will create the path if it does not exist,
//  before changing the current working directory.

void IOHDF5Handler::cd(const std::string& path, bool makedir)
{
 if (!openflag) return;
 if (makedir) mkdir(path);
 std::string ppath(tidyString(path));
 if (ppath==".") return;
 vector<std::string> newpath(currwd);
 vector<std::string> upcheck(tokenize(ppath,"/",true,true));
 bool up=(upcheck.size()>0)?true:false;
 for (uint k=0;k<upcheck.size();++k){
    if (upcheck[k]!=".."){
       up=false; k=upcheck.size();}}
 if (up){
    for (uint k=0;k<upcheck.size();++k){
       if (!(newpath.empty())){
          newpath.pop_back();}}}
 else if (ppath[0]=='/'){
    newpath=tokenize(ppath);}
 else{
    vector<std::string> add=tokenize(path);
    newpath.insert(newpath.end(),add.begin(),add.end());}
 hid_t newid=H5Gopen(fid,form_absolute_path(newpath).c_str(),H5P_DEFAULT);
 check_for_hid_failure2(newid,"cd failed for path",path);
 H5Gclose(currid);
 currid=newid;
 currwd=newpath;
 checksum=0;
}

   //  query if the group or directory "path" exists in the file

bool IOHDF5Handler::queryDir(const std::string& path) const
{
 return (query_obj(path)=='G');
}

   //  query if the data object "objname" exists in the file

bool IOHDF5Handler::queryData(const std::string& objname) const
{
 return (query_obj(objname)=='D');
}


   //  Returns 'G' if "objname" exists and is a group
   //  Returns 'D' if "objname" exists and is a dataset
   //  Returns 'N' for all other cases

char IOHDF5Handler::query_obj(const std::string& objname) const
{
 if (!openflag) return 'N';
 std::string ppath(tidyString(objname));
 const hid_t* cwdptr=(ppath[0]=='/')? (&fid) : (&currid);
 vector<std::string> dirlist(tokenize(ppath));
 htri_t exists=H5Lexists(*cwdptr,dirlist[0].c_str(),H5P_DEFAULT);
 if (exists<=0) return 'N';
 exists=H5Oexists_by_name(*cwdptr,dirlist[0].c_str(),H5P_DEFAULT);
 if (exists<=0) return 'N';
 std::string tmpstring=dirlist[0];
 for (uint i=1; i<dirlist.size(); i++){
    tmpstring+="/"+dirlist[i];
    exists=H5Lexists(*cwdptr,tmpstring.c_str(),H5P_DEFAULT);
    if (exists<=0) return 'N';
    exists=H5Oexists_by_name(*cwdptr,tmpstring.c_str(),H5P_DEFAULT);
    if (exists<=0) return 'N';}
 H5O_info_t objinfo;
 herr_t errhandle=H5Oget_info_by_name(*cwdptr,objname.c_str(),&objinfo,H5P_DEFAULT);
 if (errhandle<0){
    check_for_herr_failure(errhandle,"query_obj failed");}
 if (objinfo.type==H5O_TYPE_GROUP) return 'G';
 else if (objinfo.type==H5O_TYPE_DATASET) return 'D';
 return 'N';
}


std::set<std::string> IOHDF5Handler::getAllDataNames() const
{
 set<std::string> result;
 if (openflag){
    std::string path;
    collect_data_names(fid,path,"",result);
    result.erase("/Info/FIdentifier");
    result.erase("/Info/Endianness");}
 return result;
}


std::set<std::string> IOHDF5Handler::getAllDirNames() const
{
 set<std::string> result;
 if (openflag){
    std::string path;
    collect_dir_names(fid,path,"",result);}
 return result;
}


void IOHDF5Handler::collect_data_names(hid_t loc_id, const std::string& path,
                    const char* name, set<std::string>& collected_names) const
{
 H5O_info_t info;
 herr_t status = H5Oget_info(loc_id,&info);
 if (status<0){
    check_for_herr_failure(status,"collect_data_names failed");}
 if (info.type == H5O_TYPE_DATASET){
    std::string newpath(path); std::string add(name);
    if (add.length()>0){ newpath+="/"; newpath+=add;}
    collected_names.insert(newpath); return;}
 if (info.type == H5O_TYPE_GROUP){
    std::string newpath(path); std::string add(name);
    if (add.length()>0){ newpath+="/"; newpath+=add;}
    H5G_info_t ginfo; 
    status = H5Gget_info(loc_id,&ginfo);
    if (status<0){
       check_for_herr_failure(status,"collect_data_names failed");}
    for (hsize_t i=0;i<ginfo.nlinks;i++){
       ssize_t nsize = H5Lget_name_by_idx(loc_id,".",H5_INDEX_NAME,H5_ITER_INC,
                                          i, NULL, 0, H5P_DEFAULT)+1;
       char* nextname=new char[nsize];
       nsize = H5Lget_name_by_idx(loc_id,".",H5_INDEX_NAME,H5_ITER_INC,i,nextname,
                                  size_t(nsize), H5P_DEFAULT);
       hid_t next_id = H5Oopen(loc_id,nextname,H5P_DEFAULT);
       check_for_hid_failure(next_id," failure in collect_data_names");
       collect_data_names(next_id,newpath,nextname,collected_names);
       status = H5Oclose(next_id);
       delete [] nextname;}}
}


void IOHDF5Handler::collect_dir_names(hid_t loc_id, const std::string& path,
                    const char* name, set<std::string>& collected_dir_names) const
{
 H5O_info_t info;
 herr_t status = H5Oget_info(loc_id,&info);
 if (status<0){
    check_for_herr_failure(status,"collect_dir_names failed");}
 if (info.type == H5O_TYPE_GROUP){
    std::string newpath(path); std::string add(name);
    add=tidyString(add);
    if (add.length()>0){ 
       newpath+="/"; newpath+=add;
       collected_dir_names.insert(newpath);}
    H5G_info_t ginfo; 
    status = H5Gget_info(loc_id,&ginfo);
    if (status<0){
       check_for_herr_failure(status,"collect_dir_names failed");}
    for (hsize_t i=0;i<ginfo.nlinks;i++){
       ssize_t nsize = H5Lget_name_by_idx(loc_id,".",H5_INDEX_NAME,H5_ITER_INC,
                                          i, NULL, 0, H5P_DEFAULT)+1;
       char* nextname=new char[nsize];
       nsize = H5Lget_name_by_idx(loc_id,".",H5_INDEX_NAME,H5_ITER_INC,i,nextname,
                                  size_t(nsize), H5P_DEFAULT);
       hid_t next_id = H5Oopen(loc_id,nextname,H5P_DEFAULT);
       check_for_hid_failure(next_id," failure in collect_dir_names");
       collect_dir_names(next_id,newpath,nextname,collected_dir_names);
       status = H5Oclose(next_id);
       delete [] nextname;}}
}


std::set<std::string> IOHDF5Handler::getDataNamesInCurrentDir() const
{
 set<std::string> result;
 if (openflag){
    collect_data_names_in_currdir(currid,result,H5O_TYPE_DATASET);}
 return result;
}


std::set<std::string> IOHDF5Handler::getDirNamesInCurrentDir() const
{
 set<std::string> result;
 if (openflag){
    collect_data_names_in_currdir(currid,result,H5O_TYPE_GROUP);}
 return result;
}


void IOHDF5Handler::collect_data_names_in_currdir(hid_t loc_id,
                         set<std::string>& collected_names, H5O_type_t query_type) const
{
 H5G_info_t ginfo; 
 herr_t status = H5Gget_info(loc_id,&ginfo);
 if (status<0){
    check_for_herr_failure(status,"collect_data_names_in_currdir failed");}
 for (hsize_t i=0;i<ginfo.nlinks;i++){
    ssize_t nsize = H5Lget_name_by_idx(loc_id,".",H5_INDEX_NAME,H5_ITER_INC,
                                       i, NULL, 0, H5P_DEFAULT)+1;
    char* nextname=new char[nsize];
    nsize = H5Lget_name_by_idx(loc_id,".",H5_INDEX_NAME,H5_ITER_INC,i,nextname,
                               size_t(nsize), H5P_DEFAULT);
    hid_t next_id = H5Oopen(loc_id,nextname,H5P_DEFAULT);
    check_for_hid_failure(next_id," failure in collect_data_names");
    H5O_info_t info;
    herr_t status = H5Oget_info(next_id,&info);
    if (status<0){
       check_for_herr_failure(status,"collect_data_names_in_currdir failed");}
    if (info.type == query_type){
       std::string item(nextname);
       if (item.length()>0){
          collected_names.insert(item);}}
    status = H5Oclose(next_id);
    delete [] nextname;}
}


     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened, its file type matches 
     // "filetype_id", and a string value is found; returns false otherwise.  
     // The string is returned in "stringvalue".  This routine is used by 
     // objects in data_io_handler.h when the multi-file handlers build 
     // up their maps of file keys.

bool IOHDF5Handler::peekString(std::string& stringvalue, const std::string& path,
                               const std::string& filename,
                               const std::string& filetype_id) const
{
 stringvalue.clear();
 vector<string> paths(2),strvalues;
 paths[0]="/Info/FIdentifier";
 paths[1]=path;
 peek_strings(filename,paths,strvalues);
 if (filetype_id==strvalues[0]){
    stringvalue=strvalues[1];}
 return !(stringvalue.empty());
}

     // return true if file can be opened and read; false other.
     // Return ID string in "stringvalue".

bool IOHDF5Handler::peekID(std::string& stringvalue, const std::string& filename) const
{
 stringvalue.clear();
 vector<string> paths(1),strvalues;
 paths[0]="/Info/FIdentifier";
 peek_strings(filename,paths,strvalues);
 stringvalue=strvalues[0];
 return !(stringvalue.empty());
}

      // Clear, then reset finfo

void IOHDF5Handler::close()
{ 
 if (!openflag) return;
 clear();
}

     // Close the file, reset all data members except finfo

void IOHDF5Handler::clear()
{ 
 if (openflag) file_close();
 m_filename.clear();
 openflag=false;
 endian_convert=false;
 endian_format='U';  // undefined
 checksum=0;
 read_only=true;
 read_mode=true;
 checksum_on=false;
 is_new_file=false;
 fid=-1;
 currid=-1;
 currwd.clear();
}

void IOHDF5Handler::turnOnChecksum()
{
 if (checksum_on) return;
 checksum_on=true;
 checksum=0;
}

void IOHDF5Handler::turnOffChecksum()
{
 checksum_on=false;
}

void IOHDF5Handler::resetChecksum()
{
 checksum=0;
}

ByteHandler::n_uint32_t IOHDF5Handler::getChecksum()
{
 check_for_failure(!checksum_on,"Invalid call to getChecksum since checksums not turned on");
 return checksum;
}


        // Print out file ID. Does not 
        // change file pointers.

void IOHDF5Handler::printFileID()
{
 char endian;
 std::string ID_string;
 readIDstring(endian,ID_string);
 std::cout << "File = "<<m_filename<<std::endl<<"ID string = <"
           <<tidyString(ID_string)<<">"<<std::endl;
}

   // Write ID string and endianness ('B' or 'L')

void IOHDF5Handler::writeIDstring(const std::string& ID_string)
{
 if (openflag){
    std::string buf(tidyString(ID_string));
    mkdir("/Info");
    write("/Info/FIdentifier",buf);
    write("/Info/Endianness",endian_format);}
 checksum=0;
}


   // Reads endian format and ID string from current file.
   // Results are returned in "endianness" and "ID_string".
   // All file pointers updated.

void IOHDF5Handler::readIDstring(char& endianness, std::string& ID_string)
{
 if (openflag){
    read("/Info/FIdentifier",ID_string);
    read("/Info/Endianness",endianness);}
 else{
    endianness='U';
    ID_string.clear();}
 checksum=0;
}


      // Checks if a file exists.
      
bool IOHDF5Handler::fileExists()
{
 bool result;
 result = (access(m_filename.c_str(),F_OK) == 0) ? true : false;
 return result;
}



// *************************************************************************



void IOHDF5Handler::check_for_failure(int errcode, const std::string& mesg) const
{
 if (errcode==IOHDF5Handler::IO_SUCCESS) return;
 std::cout << "IOHDF5Handler error with file "<<m_filename<<":"
             <<endl<<"  "<<mesg<<endl;
 exit(errcode);
}

void IOHDF5Handler::check_for_hid_failure(hid_t id, const std::string& mesg) const
{
 if (id<0){
    std::cout << "IOHDF5Handler error with file "<<m_filename<<":"
                <<endl<<"  "<<mesg<<endl;
    exit(id);}
}

void IOHDF5Handler::check_for_hid_failure2(hid_t id, const std::string& mesg1, 
                                           const std::string& mesg2) const
{
 if (id<0){
    std::cout << "IOHDF5Handler error with file "<<m_filename<<":"
                <<endl<<"  "<<mesg1<<" "<<mesg2<<endl;
    exit(id);}
}

void IOHDF5Handler::check_for_herr_failure(herr_t status, const std::string& mesg) const
{
 if (status<0){
    std::cout << "IOHDF5Handler error with file "<<m_filename<<":"
                <<endl<<"  "<<mesg<<endl;
    exit(status);}
}

void IOHDF5Handler::delete_file()
{
 check_for_failure(remove(m_filename.c_str()),"Failure deleting file");
}


void IOHDF5Handler::file_close()
{
 herr_t status1=H5Gclose(currid);
 herr_t status2=H5Fclose(fid);
 check_for_herr_failure(status1,"Failure during close");
 check_for_herr_failure(status2,"Failure during close");
}

  // removes tabs, newline, linefeed characters, then trims
  // leading and trailing characters specified in "leadtrailchars"

std::string IOHDF5Handler::tidyString(const std::string& str, 
                                      const std::string& leadtrailchars) const  
{
 std::string tmp;
 for (size_t i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 size_t start=tmp.find_first_not_of(leadtrailchars);
 if (start==std::string::npos) return "";
 size_t len=tmp.find_last_not_of(leadtrailchars)-start+1;
 return tmp.substr(start,len);
}

   // Cuts a string into tokens, which are separated by "delimiter",
   // returning results in a vector of strings which do not contain the
   // delimiting characters.  Tokens which are just "." or empty are discarded.

std::vector<std::string> IOHDF5Handler::tokenize(const std::string& str, 
                            const std::string& delimiter, bool allowdot, 
                            bool allowdotdot) const
{
 std::string s(tidyString(str," /"));
 std::vector<std::string> tokens;
 size_t start = 0;
 size_t end = s.find(delimiter);
 std::string token;
 while (end != std::string::npos) {
    token=s.substr(start, end - start);
    check_token(token,allowdot,allowdotdot); 
    if ((token.size()>0)&&(token!=".")) tokens.push_back(token);
    start = end + delimiter.size();
    end = s.find(delimiter, start);}
 token=s.substr(start, end - start);
 check_token(token,allowdot,allowdotdot); 
 if ((token.size()>0)&&(token!=".")) tokens.push_back(token);
 return tokens;
}


void IOHDF5Handler::check_token(const std::string& token, bool allowdot, 
                                bool allowdotdot) const
{
// if (token.find_first_of(" ")!=std::string::npos)
//    check_for_failure(IO_ERR_OTHER," blank spaces in object names not allowed here!!");
 if ((!allowdotdot)&&(token==".."))
    check_for_failure(IO_ERR_OTHER,".. not allowed for group names");
 if ((!allowdot)&&(token=="."))
    check_for_failure(IO_ERR_OTHER,". not allowed for group names");
 for (uint i=0;i<token.length();++i){
    if ((iswspace(token[i]))&&(token[i]!=' ')){
       check_for_failure(IO_ERR_OTHER," white space other than blanks not allowed in object names!!");}}
}


void IOHDF5Handler::check_path(const std::string& str, bool allowdot, 
                               bool allowdotdot) const
{
 std::string delimiter("/");
 std::string s(tidyString(str," /"));
 size_t start = 0;
 size_t end = s.find(delimiter);
 std::string token;
 while (end != std::string::npos) {
    token=s.substr(start, end - start);
    check_token(token,allowdot,allowdotdot); 
    start = end + delimiter.size();
    end = s.find(delimiter, start);}
 token=s.substr(start, end - start);
 check_token(token,allowdot,allowdotdot);
}

    // take a vector of strings and concatenates them into a single
    // string with a "/" between the vector elements.  Path always starts with
    // a "/" so is an absolute path

std::string IOHDF5Handler::form_absolute_path(const std::vector<std::string>& path_links) const
{
 if (path_links.empty()) return std::string("/");
 std::string result="";
 for (std::vector<std::string>::const_iterator it=path_links.begin();it!=path_links.end();++it)
    result+="/"+(*it);
 return result;
}

// *******************************************************************

void IOHDF5Handler::write(const std::string& objname, const std::string& output)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write when no open file");}
 if (read_only){
    check_for_failure(IO_ERR_ACCESS,"Attempt to write to read-only file");}
 check_path(objname);
 int n=output.length();
 hid_t typid=H5Tcreate(H5T_STRING,n+1); // add 1 for terminating null character
 H5Tset_cset(typid,H5T_CSET_UTF8);
 hid_t spaceid=H5Screate(H5S_SCALAR);
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataid=H5Dcreate2(*cwdptr,objname.c_str(),typid,spaceid, 
                         H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 herr_t status =  H5Dwrite(dataid,typid,H5S_ALL,H5S_ALL,H5P_DEFAULT,output.c_str());
 check_for_herr_failure(status,"Could not write string data");
 H5Tclose(typid);
 H5Dclose(dataid);
 H5Sclose(spaceid);
 if (read_mode){
    read_mode=false; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, output.c_str(), output.length());
}


void IOHDF5Handler::write(const std::string& objname, const char& output)
{
 char str[2]; str[0]=output; str[1]='\0';
 write(objname,std::string(str));
}

void IOHDF5Handler::write(const std::string& objname, const int& output)
{
 write_atomics<int>(objname,&output,1,dtype_int,H5T_NATIVE_INT);
}

void IOHDF5Handler::write(const std::string& objname, const unsigned int& output)
{
 write_atomics<unsigned int>(objname,&output,1,dtype_uint,H5T_NATIVE_UINT);
}

void IOHDF5Handler::write(const std::string& objname, const long int& output)
{
 write_atomics<long int>(objname,&output,1,dtype_long,H5T_NATIVE_LONG);
}

void IOHDF5Handler::write(const std::string& objname, const unsigned long int& output)
{
 write_atomics<unsigned long int>(objname,&output,1,dtype_ulong,H5T_NATIVE_ULONG);
}

void IOHDF5Handler::write(const std::string& objname, const long long int& output)
{
 write_atomics<long long int>(objname,&output,1,dtype_llong,H5T_NATIVE_LLONG);
}

void IOHDF5Handler::write(const std::string& objname, const unsigned long long int& output)
{
 write_atomics<unsigned long long int>(objname,&output,1,dtype_ullong,H5T_NATIVE_ULLONG);
}

void IOHDF5Handler::write(const std::string& objname, const float& output)
{
 write_atomics<float>(objname,&output,1,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::write(const std::string& objname, const double& output)
{
 write_atomics<double>(objname,&output,1,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::write(const std::string& objname, const fcmplx& output)
{
 write_complex_values<float>(objname,&output,1,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::write(const std::string& objname, const dcmplx& output)
{
 write_complex_values<double>(objname,&output,1,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<char>& output)
{
 std::string s(output.begin(), output.end());
 write(objname,s);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<int>& output)
{
 write_atomics<int>(objname,output.data(),output.size(),dtype_int,H5T_NATIVE_INT);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<unsigned int>& output)
{
 write_atomics<unsigned int>(objname,output.data(),output.size(),dtype_uint,H5T_NATIVE_UINT);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<long int>& output)
{
 write_atomics<long int>(objname,output.data(),output.size(),dtype_long,H5T_NATIVE_LONG);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<unsigned long int>& output)
{
 write_atomics<unsigned long int>(objname,output.data(),output.size(),dtype_ullong,H5T_NATIVE_ULLONG);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<float>& output)
{
 write_atomics<float>(objname,output.data(),output.size(),dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<double>& output)
{
 write_atomics<double>(objname,output.data(),output.size(),dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<fcmplx>& output)
{
 write_complex_values<float>(objname,output.data(),output.size(),dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::write(const std::string& objname, const std::vector<dcmplx>& output)
{
 write_complex_values<double>(objname,output.data(),output.size(),dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::write(const std::string& objname, const Array<int>& output)
{
 write_array<int>(objname,output,dtype_int,H5T_NATIVE_INT);
}

void IOHDF5Handler::write(const std::string& objname, const Array<unsigned int>& output)
{
 write_array<unsigned int>(objname,output,dtype_uint,H5T_NATIVE_UINT);
}

void IOHDF5Handler::write(const std::string& objname, const Array<long int>& output)
{
 write_array<long int>(objname,output,dtype_long,H5T_NATIVE_LONG);
}

void IOHDF5Handler::write(const std::string& objname, const Array<unsigned long int>& output)
{
 write_array<unsigned long int>(objname,output,dtype_ulong,H5T_NATIVE_ULONG);
}

void IOHDF5Handler::write(const std::string& objname, const Array<float>& output)
{
 write_array<float>(objname,output,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::write(const std::string& objname, const Array<double>& output)
{
 write_array<double>(objname,output,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::write(const std::string& objname, const Array<fcmplx>& output)
{
 write_complex_array<fcmplx>(objname,output,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::write(const std::string& objname, const Array<dcmplx>& output)
{
 write_complex_array<dcmplx>(objname,output,dtype_double,H5T_NATIVE_DOUBLE);
}


// ********************************************************************

void IOHDF5Handler::read(const std::string& objname, std::string& input)
{
 if (!openflag){
    check_for_failure(IO_ERR_ACCESS,"Attempt to read when no open file");}
 std::string obj(tidyString(objname));
 hid_t* cwdptr=(obj[0]=='/')? (&fid) : (&currid);
 hid_t dataset_id = H5Dopen2(*cwdptr,objname.c_str(),H5P_DEFAULT);
 check_for_hid_failure2(dataset_id,"Could not find object name",objname);
 hid_t dtype = H5Dget_type(dataset_id);
 if (H5Tget_class(dtype)!=H5T_STRING){
    check_for_failure(IO_ERR_OTHER,"Datatype mismatch during string read");}
 hsize_t dsize=H5Tget_size(dtype);
 char* buffer=new char[dsize];
 hid_t nat_type_id=H5Tget_native_type(dtype,H5T_DIR_ASCEND);
 herr_t status = H5Dread(dataset_id,nat_type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,buffer);
 check_for_herr_failure(status,"Could not read string data"); 
 input=std::string(buffer,dsize-1);
 delete [] buffer;
 H5Dclose(dataset_id);
 H5Tclose(nat_type_id);
 H5Tclose(dtype);
 if (!read_mode){
    read_mode=true; checksum=0;}
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, input.c_str(), input.length());
}


void IOHDF5Handler::read(const std::string& objname, char& input)
{
 std::string s;
 read(objname,s);
 if (s.length()!=1){
    check_for_failure(IO_ERR_OTHER,"Failure during a character read");}
 input=s[0];
}

void IOHDF5Handler::read(const std::string& objname, int& input)
{
 read_atomic<int>(objname,input,dtype_int,H5T_NATIVE_INT);
}

void IOHDF5Handler::read(const std::string& objname, unsigned int& input)
{
 read_atomic<unsigned int>(objname,input,dtype_uint,H5T_NATIVE_UINT);
}

void IOHDF5Handler::read(const std::string& objname, long int& input)
{
 read_atomic<long int>(objname,input,dtype_long,H5T_NATIVE_LONG);
}

void IOHDF5Handler::read(const std::string& objname, unsigned long int& input)
{
 read_atomic<unsigned long int>(objname,input,dtype_ulong,H5T_NATIVE_ULONG);
}

void IOHDF5Handler::read(const std::string& objname, long long int& input)
{
 read_atomic<long long int>(objname,input,dtype_llong,H5T_NATIVE_LLONG);
}

void IOHDF5Handler::read(const std::string& objname, unsigned long long int& input)
{
 read_atomic<unsigned long long int>(objname,input,dtype_ullong,H5T_NATIVE_ULLONG);
}

void IOHDF5Handler::read(const std::string& objname, float& input)
{
 read_atomic<float>(objname,input,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::read(const std::string& objname, double& input)
{
 read_atomic<double>(objname,input,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::read(const std::string& objname, fcmplx& input)
{
 read_complex<float>(objname,input,dtype_float,H5T_NATIVE_FLOAT);
}
 
void IOHDF5Handler::read(const std::string& objname, dcmplx& input)
{
 read_complex<double>(objname,input,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<char>& input)
{
 std::string s;
 read(objname,s);
 input=std::vector<char>(s.begin(), s.end());
}

void IOHDF5Handler::read(const std::string& objname, std::vector<int>& input)
{
 read_atomics<int>(objname,input,dtype_int,H5T_NATIVE_INT);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<unsigned int>& input)
{
 read_atomics<unsigned int>(objname,input,dtype_uint,H5T_NATIVE_UINT);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<long int>& input)
{
 read_atomics<long int>(objname,input,dtype_long,H5T_NATIVE_LONG);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<unsigned long int>& input)
{
 read_atomics<unsigned long int>(objname,input,dtype_ulong,H5T_NATIVE_ULONG);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<float>& input)
{
 read_atomics<float>(objname,input,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<double>& input)
{
 read_atomics<double>(objname,input,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<fcmplx>& input)
{
 read_complex_vector<float>(objname,input,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::read(const std::string& objname, std::vector<dcmplx>& input)
{
 read_complex_vector<double>(objname,input,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::read(const std::string& objname, Array<int>& input)
{
 read_array<int>(objname,input,dtype_int,H5T_NATIVE_INT);
}

void IOHDF5Handler::read(const std::string& objname, Array<unsigned int>& input)
{
 read_array<unsigned int>(objname,input,dtype_uint,H5T_NATIVE_UINT);
}

void IOHDF5Handler::read(const std::string& objname, Array<long int>& input)
{
 read_array<long int>(objname,input,dtype_long,H5T_NATIVE_LONG);
}

void IOHDF5Handler::read(const std::string& objname, Array<unsigned long int>& input)
{
 read_array<unsigned long int>(objname,input,dtype_ulong,H5T_NATIVE_ULONG);
}

void IOHDF5Handler::read(const std::string& objname, Array<float>& input)
{
 read_array<float>(objname,input,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::read(const std::string& objname, Array<double>& input)
{
 read_array<double>(objname,input,dtype_double,H5T_NATIVE_DOUBLE);
}

void IOHDF5Handler::read(const std::string& objname, Array<fcmplx>& input)
{
 read_complex_array<fcmplx>(objname,input,dtype_float,H5T_NATIVE_FLOAT);
}

void IOHDF5Handler::read(const std::string& objname, Array<dcmplx>& input)
{
 read_complex_array<dcmplx>(objname,input,dtype_double,H5T_NATIVE_DOUBLE);
}

// ********************************************************************

     // Open file, read strings at particular paths, then close. If
     // any path cannot be found, an empty string is returned for that path.
     // Note that the paths must be absolute path names.

void IOHDF5Handler::peek_strings(const std::string& filename,
                                 const std::vector<std::string>& paths,
                                 std::vector<std::string>& stringvalues) const
{
 stringvalues.clear();
 for (uint k=0;k<paths.size();++k){
    stringvalues.push_back("");}
 std::string fname=tidyString(filename);
 if (fname.empty()) return;

    // Suppress HDF5 error messages: save current error handler, then turn off 
 H5E_auto2_t hdf5error_func;
 void *hdf5error_data;
 H5Eget_auto(H5E_DEFAULT, &hdf5error_func, &hdf5error_data);
 H5Eset_auto(H5E_DEFAULT, NULL, NULL);

 htri_t fflag=H5Fis_hdf5(fname.c_str());
 if (fflag<=0) return;
 hid_t file_id=H5Fopen(fname.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
 if (file_id<0) return;

 for (uint k=0;k<paths.size();k++){
    std::string ppath(tidyString(paths[k]));
    if (ppath.empty()){ continue;}
    hid_t dataset_id = H5Dopen2(file_id,ppath.c_str(),H5P_DEFAULT);
    if (dataset_id<0){ continue;}
    hid_t dtype = H5Dget_type(dataset_id);
    if (H5Tget_class(dtype)!=H5T_STRING){ 
       H5Tclose(dtype); H5Dclose(dataset_id); continue;}
    hsize_t dsize=H5Tget_size(dtype);
    char* buffer=new char[dsize];
    hid_t nat_type_id=H5Tget_native_type(dtype,H5T_DIR_ASCEND);
    herr_t status = H5Dread(dataset_id,nat_type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,buffer);
    if (status>=0){
       stringvalues[k]=std::string(buffer,dsize-1);}
    delete [] buffer;
    H5Dclose(dataset_id);
    H5Tclose(nat_type_id);
    H5Tclose(dtype);}

   // Restore previous error handler
 H5Eset_auto(H5E_DEFAULT, hdf5error_func, hdf5error_data);
 H5Fclose(file_id);
}

// ********************************************************************
