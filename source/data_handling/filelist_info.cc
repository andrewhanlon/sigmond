#include "filelist_info.h"
#include <stdexcept>
using namespace std;

namespace LaphEnv {

 // *************************************************************


FileListInfo::FileListInfo(XMLHandler& xml_in)
{ 
 XMLHandler xmlr(xml_in,"FileListInfo");
 set_info(xmlr);
}


FileListInfo::FileListInfo(XMLHandler& xml_in, const string& outertag)
{
 XMLHandler xmlq(xml_in,outertag);
 XMLHandler xmlr(xmlq,"FileListInfo");
 set_info(xmlr);
}


FileListInfo::FileListInfo(const std::string& stub, int min_suffix, 
                           int max_suffix, bool over_write)
{
 try{set_info(stub,min_suffix,max_suffix,over_write);}
 catch(const std::exception& msg){
    cerr << "invalid FileListInfo construction"<<endl;
    cerr << msg.what()<<endl;
    exit(1);}
}


void FileListInfo::set_info(XMLHandler& xmlr)
{
 string stub;
 xmlread(xmlr,"FileNameStub",stub,"FileListInfo");

 int min_suffix=0;
 if (xml_tag_count(xmlr,"MinFileNumber")==1)
    xmlread(xmlr,"MinFileNumber",min_suffix,"FileListInfo");

 int max_suffix;
 xmlread(xmlr,"MaxFileNumber",max_suffix,"FileListInfo");

 bool overwrite = false;  // protect mode
 if (xml_tag_count(xmlr,"FileMode")==1){
    string fmode;
    xmlread(xmlr,"FileMode",fmode,"FileListInfo");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") overwrite=true;}

 try{set_info(stub,min_suffix,max_suffix,overwrite);}
 catch(const std::exception& msg){
    xml_cerr(xmlr,string(msg.what())); xmlreadfail(xmlr,"FileListInfo");}
}


void FileListInfo::set_info(const std::string& stub, int min_suffix, 
                            int max_suffix, bool over_write)
{
 m_file_stub=tidyString(stub);
 if (m_file_stub.empty()){
    throw(std::invalid_argument("Blank file name in FileListInfo"));}
 m_max_file_number=max_suffix;
 m_overwrite_mode = over_write;  // protect mode
 m_min_file_number=min_suffix;
 if ((m_min_file_number>m_max_file_number)||(m_min_file_number<0)){
    throw(std::invalid_argument("minimum file number > maximum file number in FileListInfo"));}
} 


FileListInfo::FileListInfo(const FileListInfo& fin)
   : m_file_stub(fin.m_file_stub), 
     m_max_file_number(fin.m_max_file_number),
     m_min_file_number(fin.m_min_file_number),
     m_overwrite_mode(fin.m_overwrite_mode) {}


FileListInfo& FileListInfo::operator=(const FileListInfo& fin)
{
 m_file_stub=fin.m_file_stub;
 m_max_file_number=fin.m_max_file_number;
 m_min_file_number=fin.m_min_file_number;
 m_overwrite_mode=fin.m_overwrite_mode;
 return *this;
}


std::string FileListInfo::getFileName(int suffix) const
{
 stringstream fs;
 fs << m_file_stub << "." << suffix;
 return fs.str();
}

int FileListInfo::getFirstAvailableSuffix() const
{
 for (int suffix=m_min_file_number;suffix<=m_max_file_number;suffix++){
    string filename=getFileName(suffix);
    if (!fileExists(filename)) return suffix;}
 cerr << "no suffix numbers are available for writing"<<endl;
 cerr << " ... increase maxFilenumber"<<endl;
 throw(std::invalid_argument("no suffix numbers available for writing"));
}


bool FileListInfo::operator==(const FileListInfo& in) const
{
 return  ((m_file_stub==in.m_file_stub)           
        &&(m_max_file_number==in.m_max_file_number)
        &&(m_min_file_number==in.m_min_file_number));
}

void FileListInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("FileListInfo");
 xmlout.put_child("FileNameStub",m_file_stub);
 xmlout.put_child("MinFileNumber",make_string(m_min_file_number));
 xmlout.put_child("MaxFileNumber",make_string(m_max_file_number));
 if (m_overwrite_mode) xmlout.put_child("FileMode","overwrite");
 else xmlout.put_child("FileMode","protect");
}

string FileListInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}

string FileListInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

// ***************************************************************
}

