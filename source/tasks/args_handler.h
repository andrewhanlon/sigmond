#ifndef ARGS_HANDLER_H
#define ARGS_HANDLER_H

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include "xml_handler.h"

// *************************************************************************
// *                                                                       *
// *   A common task in "sigmond" is to read XML and extract data values,  *
// *   such as integers, floats, and complicated objects, then often       *
// *   these input data values are echoed back in log files.  The class    *
// *   "ArgsHandler" is meant for this purpose.  With a minimal of         *
// *   code, a large variety of quantities can be extracted from XML,      *
// *   checked for correctness, and an XMLHandler is built up for echoing. *
// *   If an incorrect value is encountered, an exception is thrown.       *
// *   Optional arguments can also be handled.                             *
// *                                                                       *
// *   There are two main constructors: The first just uses an XMLHandler  *
// *   as its single argument.  The current location of this handler       *
// *   becomes the root tag for the input xml.  All data is searched       *
// *   for in the children and descendents of the new root tag.            *
// *   The second constructor takes an XMLHandler and a tag string,        *
// *   plus an optional boolean (false by default).  If the boolean        *
// *   is false or absent, a seek_unique() is done to find the tag         *
// *   name, which becomes the new root tag.  The search is done           *
// *   among all descendents of the current location.  If the              *
// *   boolean is true, the search for the tag name is done only to root   *
// *   and among the children of the current location.  A string exception *
// *   is thrown is an error occurs, such as not being able to find        *
// *   the tag name, or encountering multiple instances of the tag         *
// *   name.  A third constructor is used for getting sub-level arguments  *
// *   (see below).                                                        *
// *                                                                       *
// *   There are two constructors which take a **set** of roottags.        *
// *   If the handler finds just ONE occurrence of ONE of these roottags,  *
// *   it sets up fine (seeking in child nodes only or deeper, depending   *
// *   on "seek_to_child_only"). The user can then query with              *
// *   "getInputRootTag()" to see which root tag was found.  An exception  *
// *   is thrown if none of the roottags is found, or multiple instances   *
// *   are found.                                                          *
// *                                                                       *
// *   All arguments read are placed in an output XMLHandler which         *
// *   can be echoed once all data have been input.  By default,           *
// *   all input is echoed, but echoing can be turned off using            *
// *   "setOffEcho", then turned back on later using "setOnEcho".          *
// *   In this way, the user can selectively echo the input.               *
// *                                                                       *
// *   In searching for the tag name for some requested data, only the     *
// *   child nodes are searched!!  The tag name must occur only once       *
// *   among the children, or an exception is thrown.  An exception is     *
// *   thrown if not found, unless an "optional" get is used.              *
// *                                                                       *
// *   The member "queryTag" searches for the given tag name among the     *
// *   current node and all of its children nodes (not grandchildren).     *
// *   It returns false if not found, true if found once, and throws       *
// *   an exception if multiple instances are found.                       *
// *                                                                       *
// *   To search farther down in descendents, a new ArgsHandler object     *
// *   must be created with one of the children as the new root, then      *
// *   that object is used to get data from the children of that child     *
// *   (the grandchildren of the original root).  A constructor which      *
// *   takes the original ArgsHandler as an argument is available to       *
// *   facilitate this. One then uses an "insert" into the top ArgsHandler *
// *   to facilitate the XML echo.                                         *
// *                                                                       *
// *   Usage:                                                              *
// *                                                                       *
// *     int k;      // define various objects for input                   *
// *     double g;                                                         *
// *     string str;                                                       *
// *     int m=3;    // an optional argument with a default value          *
// *                                                                       *
// *     try{                                                              *
// *        ArgsHandler gin(xml1);   // get input from xml1                *
// *        gin.getInt("Tag1",k);                                          *
// *        gin.getReal("Tag2",g);                                         *
// *        gin.getString("Tag3",str);                                     *
// *        gin.getOptionalInt("Tag4",m);                                  *
// *        int p=gin.getInt("Tag5");                                      *
// *        OperatorInfo op1; gin.getItem("OpTag1",op1);                   *
// *        OperatorInfo op2(gin.getItem<OperatorInfo>("OpTag2"));         *
// *                                                                       *
// *        if (gin.queryTag("Tag6")){                                     *
// *           ArgsHandler gin2(gin,"Tag6");                               *
// *           int q=gin2.getInt("Grandchild",q);                          *
// *           gin.insert(gin2);}   // to include for echo                 *
// *                                                                       *
// *        XMLHandler xmlout;                                             *
// *        gin.echo(xmlout);}                                             *
// *     catch(const std::exception& errmsg){                              *
// *        throw;}                                                        *
// *                                                                       *
// *                                                                       *
// *************************************************************************



class ArgsHandler
{

    mutable XMLHandler m_xmlin;
    XMLHandler m_xmlout;
    bool m_echo;

 public:

           //  current location in xmlin becomes root of m_xmlin

    ArgsHandler(XMLHandler& xmlin) 
       : m_xmlin(xmlin), m_xmlout(xmlin.get_tag_name()), m_echo(true) 
     {m_xmlin.set_exceptions_on();}

           //  if "seek_to_child_only" is false, does seek_unique in xmlin so 
           //  "roottag" is new root of m_xmlin, seeking in children, 
           //  grand children, etc. recursively; if "seek_to_child_only" is 
           //  true, it seeks for "roottag" only in root and among the children
           //  of the current location in xmlin

    ArgsHandler(XMLHandler& xmlin, const std::string& roottag,
                bool seek_to_child_only=true)
        : m_xmlin(xmlin,roottag,seek_to_child_only), m_xmlout(roottag), 
          m_echo(true)
     {m_xmlin.set_exceptions_on();}


    ArgsHandler(ArgsHandler& gtop, const std::string& roottag,
                bool seek_to_child_only=true)
        : m_xmlin(gtop.m_xmlin,roottag,seek_to_child_only), m_xmlout(roottag), 
          m_echo(true) {}

          // The two constructors below take a **set** of roottags.
          // If the handler finds just ONE occurrence of ONE of these
          // roottags, it sets up fine (seeking in child nodes only
          // or deeper, depending on "seek_to_child_only"). The user
          // can then query with "getInputRootTag()" to see which
          // root tag was found.  An exception is thrown if none
          // of the roottags is found, or multiple instances are found.

    ArgsHandler(XMLHandler& xmlin, const std::set<std::string>& roottags,
                bool seek_to_child_only=true)
     {assign_with_choices(xmlin,roottags,seek_to_child_only);
      m_xmlin.set_exceptions_on();}

    ArgsHandler(ArgsHandler& gtop, const std::set<std::string>& roottags,
                bool seek_to_child_only=true)
     {assign_with_choices(gtop.m_xmlin,roottags,seek_to_child_only);}


    void setOnEcho() {m_echo=true;}

    void setOffEcho() {m_echo=false;}


    void echo(XMLHandler& xmlout) const
     {xmlout=m_xmlout;}

    std::string echo() const
     {return m_xmlout.output();}

    void getInput(XMLHandler& xmlin) const
     {xmlin=m_xmlin;}

    std::string getInput() const
     {return m_xmlin.output();}

    std::string getInputRootTag() const
     {m_xmlin.seek_root(); return m_xmlin.get_node_name();}



    void insert(ArgsHandler& subhandler)
     {if (m_echo) m_xmlout.put_child(subhandler.m_xmlout);}


    bool queryTag(const std::string& tagname)
     {return m_xmlin.query_unique_to_among_children(tagname);}


    void getInt(const std::string& tagname, int& val)
     {get_basic_item(tagname,val,true);}

    int getInt(const std::string& tagname)
     {return get_basic_item<int>(tagname);}

    void getOptionalInt(const std::string& tagname, int& val)
     {get_basic_item(tagname,val,false);}



    void getUInt(const std::string& tagname, uint& val)
     {get_basic_item(tagname,val,true);}

    uint getUInt(const std::string& tagname)
     {return get_basic_item<uint>(tagname);}

    void getOptionalUInt(const std::string& tagname, uint& val)
     {get_basic_item(tagname,val,false);}



    void getReal(const std::string& tagname, double& val)
     {get_basic_item(tagname,val,true);}

    double getReal(const std::string& tagname)
     {return get_basic_item<double>(tagname);}

    void getOptionalReal(const std::string& tagname, double& val)
     {get_basic_item(tagname,val,false);}


          // string is trimmed of leading and trailing white space
          // exception thrown if string is entirely white space;
          // for "getNameItem", string must be valid name (each character 
          // is alphanumeric, underscore, or period)

    void getString(const std::string& tagname, std::string& val)
     {get_string_item(tagname,val,true);}

    std::string getString(const std::string& tagname)
     {return get_string_item(tagname,false);}

    std::string getName(const std::string& tagname)
     {return get_string_item(tagname,true);}

    void getOptionalString(const std::string& tagname, std::string& val)
     {get_string_item(tagname,val,false);}


          //  "true" value assumed and echoed if "tagname" is empty tag or
          //  <tagname>true</tagname> encountered;  "false" value
          //  assumed if "tagname" not found as child (no echo) or
          //  <tagname>false</tagname> encountered (echoed)
 
    void getBool(const std::string& tagname, bool& val)
     {get_bool_item(tagname,val,true);}

    void getOptionalBool(const std::string& tagname, bool& val)
     {get_bool_item(tagname,val,false);}

    bool getBool(const std::string& tagname)
     {bool b; get_bool_item(tagname,b,true); return b;}



    void getIntVector(const std::string& tagname, std::vector<int>& val)
     {get_basic_item(tagname,val,true);}

    std::vector<int> getIntVector(const std::string& tagname)
     {return get_basic_item<std::vector<int> >(tagname);}



          //  "tagname" is used only when outputting error message if failure
    template <typename T>
    void getItem(const std::string& tagname, T& val)
    {
     try{
        val=T(m_xmlin);
        m_xmlin.seek_root();
        XMLHandler xmlt; val.output(xmlt);
        if (m_echo) m_xmlout.put_child(xmlt);}
     catch(const std::exception& errmsg){
        throw(std::invalid_argument(error_msg(tagname,errmsg.what())));}
    }


          //  "tagname" is used only when outputting error message if failure
    template <typename T>
    T getItem(const std::string& tagname)
    {
     try{
        T val(m_xmlin);
        m_xmlin.seek_root();
        XMLHandler xmlt; val.output(xmlt);
        if (m_echo) m_xmlout.put_child(xmlt);
        return val;}
     catch(const std::exception& errmsg){
        throw(std::invalid_argument(std::string(tagname).c_str()));}
    }


 private:


    void assign_with_choices(XMLHandler& xmlin, const std::set<std::string>& roottags,
                             bool seek_to_child_only)
    {
     XMLHandler xmlt(xmlin);  // current location in xmlin is now root
     uint count=0, ct;
     std::set<std::string>::const_iterator it,rt;
     for (it=roottags.begin();it!=roottags.end();it++){
        if (seek_to_child_only) ct=xmlt.count_to_among_children(*it);
        else ct=xmlt.count(*it);
        count+=ct;
        if (ct==1) rt=it;}
     if (count!=1) throw(std::invalid_argument("ArgsHandler construction failure: no single root tag amoung choices found"));
     if (!seek_to_child_only) xmlt.seek_unique(*rt);
     else xmlt.seek_unique_to_child(*rt);
     m_xmlin.set(xmlt);
     m_xmlout.set_root(*rt);
     m_echo=true;
    }

    template <typename T>
    void get_basic_item(const std::string& tagname, T& val, bool required)
    {
     try{
        m_xmlin.seek_unique_to_child(tagname);
        std::string content=m_xmlin.get_text_content();
        m_xmlin.seek_root();
        extract_from_string(content,val);
        if (m_echo) m_xmlout.put_child(tagname,make_string(val));}
     catch(const std::exception& errmsg){
        if (required) throw(std::invalid_argument(error_msg(tagname,errmsg.what())));
        m_xmlin.seek_root();
        if (m_echo) m_xmlout.put_child(tagname,make_string(val));}
    }

    template <typename T>
    T get_basic_item(const std::string& tagname)
    {
     try{
        m_xmlin.seek_unique_to_child(tagname);
        std::string content=m_xmlin.get_text_content();
        m_xmlin.seek_root();
        T val;
        extract_from_string(content,val);
        if (m_echo) m_xmlout.put_child(tagname,make_string(val));
        return val;}
     catch(const std::exception& errmsg){
        throw(std::invalid_argument(error_msg(tagname,errmsg.what())));}
    }

    void get_string_item(const std::string& tagname, std::string& val, bool required)
    {
     try{
        m_xmlin.seek_unique_to_child(tagname);
        val=m_xmlin.get_text_content();
        m_xmlin.seek_root();
        if (m_echo) m_xmlout.put_child(tagname,val);}
     catch(const std::exception& errmsg){
        if (required) throw(std::invalid_argument(error_msg(tagname,errmsg.what())));
        m_xmlin.seek_root();
        if (m_echo) m_xmlout.put_child(tagname,val);}
    }

    std::string get_string_item(const std::string& tagname, bool is_name)
    {
     try{
        m_xmlin.seek_unique_to_child(tagname);
        std::string val=m_xmlin.get_text_content();   // white space throws exception
        m_xmlin.seek_root();
        if (is_name){
           val=tidyName(val);
           if (val.empty()) throw(std::invalid_argument("Invalid name string"));}
        if (m_echo) m_xmlout.put_child(tagname,val);
        return val;}
     catch(const std::exception& errmsg){
        throw(std::invalid_argument(error_msg(tagname,errmsg.what())));}
    }

          //  "true" value assumed and echoed if "tagname" is empty tag or
          //  <tagname>true</tagname> encountered;  "false" value
          //  assumed if "tagname" not found as child (no echo) or
          //  <tagname>false</tagname> encountered (echoed)


    void get_bool_item(const std::string& tagname, bool& val, bool required)
    {
     try{
        m_xmlin.seek_unique_to_child(tagname);}
     catch(const std::exception& xp){
        m_xmlin.seek_root();
        if (required) val=false;
        return;}
     if (m_xmlin.is_empty_tag()){
        val=true;
        m_xmlin.seek_root();
        if (m_echo) m_xmlout.put_child(tagname);
        return;}
     std::string content=m_xmlin.get_text_content();
     m_xmlin.seek_root();
     if (content=="true") val=true;
     else if (content=="false") val=false;
     else throw(std::invalid_argument(error_msg(tagname,"Invalid boolean value")));
     if (m_echo) m_xmlout.put_child(tagname,content);
    }


    const char* error_msg(const std::string& tagname, const std::string& errmsg)
    {
     m_xmlin.seek_root(); 
     return (std::string("Could not read <")+tagname
         +std::string("> when input root tag is <")+m_xmlin.get_node_name()
         +std::string("> Message: ")+errmsg).c_str();
    }
    
};


#endif
