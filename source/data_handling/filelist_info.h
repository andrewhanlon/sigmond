#ifndef FILE_LIST_INFO_H
#define FILE_LIST_INFO_H

#include "xml_handler.h"

namespace LaphEnv {


// ****************************************************************
// *                                                              *
// *  "FileListInfo" stores information about a list of files     *
// *  having a common stub and a numerical suffix from a minimum  *
// *  value to a maximum value.  A file mode is also stored       *
// *  which indicates whether files can be overwritten or not.    *
// *                                                              *
// *  Required XML input for setting the handler info:            *
// *                                                              *
// *   <FileListInfo>                                             *
// *      <FileNameStub>  ...  </FileNameStub>                    *
// *      <MaxFileNumber> ...  </MaxFileNumber>                   *
// *      <MinFileNumber> ...  </MinFileNumber> (default=0)       *
// *      <FileMode>      ...  </FileMode>   (optional)           *
// *   </FileListInfo>                                            *
// *                                                              *
// *  FileMode can be omitted (no overwrite) or can be set to     *
// *  "overwrite".                                                *
// *                                                              *
// ****************************************************************


class FileListInfo
{

   std::string m_file_stub;            // stub of files to handle
   int m_max_file_number;
   int m_min_file_number;
   bool m_overwrite_mode;         


 public:


   FileListInfo(XMLHandler& xml_in);

   FileListInfo(XMLHandler& xml_in, const std::string& outertag);

   FileListInfo(const std::string& stub, int min_suffix, int max_suffix,
                bool over_write=false);

   FileListInfo(const FileListInfo& fin);

   FileListInfo& operator=(const FileListInfo& fin);

   void setOverWrite() {m_overwrite_mode=true;}

   void setNoOverWrite() {m_overwrite_mode=false;}

   ~FileListInfo(){}



   std::string getFileStub() const { return m_file_stub; }

   int getMaxFileNumber() const { return m_max_file_number; }
  
   int getMinFileNumber() const { return m_min_file_number; }

   bool isModeOverwrite() const { return m_overwrite_mode; }

   bool operator==(const FileListInfo& rhs) const;



   std::string getFileName(int suffix) const;

   int getFirstAvailableSuffix() const;

   void output(XMLHandler& xmlout) const;

   std::string output(int indent=0) const;

   std::string str() const;


 private:

   void set_info(XMLHandler& xml_in);
   void set_info(const std::string& stub, int min_suffix, 
                 int max_suffix, bool over_write);

};


// ***************************************************************
}
#endif  
