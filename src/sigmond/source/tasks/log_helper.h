#ifndef LOG_HELPER_H
#define LOG_HELPER_H

#include <string>
#include <vector>
#include "xml_handler.h"
#include "args_handler.h"


// *************************************************************************
// *                                                                       *
// *   A common task in "sigmond" is to output XML describing a variety    *
// *   data, such as integers, floats, and complicated objects.  The       *
// *   class "LogHelper" is meant for this purpose.  With a minimal of     *
// *   code, a large variety of quantities can be output to XML.           *
// *                                                                       *
// *   Usage:                                                              *
// *                                                                       *
// *     int k;      // define various objects for output                  *
// *     double g;                                                         *
// *     string str;                                                       *
// *     bool m=false;                                                     *
// *                                                                       *
// *     try{                                                              *
// *        LogHelper gout("Roottag");                                     *
// *        gout.putInt("Tag1",k);                                         *
// *        gout.putReal("Tag2",g);                                        *
// *        gout.putString("Tag3",str);                                    *
// *        gout.putBool("Tag4",m);                                        *
// *        int p=gout.putInt("Tag5");                                     *
// *        OperatorInfo op1; gout.putItem(op1);                           *
// *                                                                       *
// *        LogHelper gout2("Tag6");                                       *
// *        gout2.putInt("Vege",4);                                        *
// *        gout.put(gout2);   // to include sub helpers                   *
// *                                                                       *
// *        XMLHandler xmlout;                                             *
// *        gout.output(xmlout);          // output into XMLHandler        *
// *        cout << gout.output()<<endl;  // output as string              *
// *                                                                       *
// *                                                                       *
// *************************************************************************



class LogHelper
{

    XMLHandler m_xmlout;

 public:

    LogHelper(const std::string& roottag)
        : m_xmlout(roottag) {}

    LogHelper() 
        : m_xmlout("LogHelper") {}   // default root tag

    LogHelper(XMLHandler& xmllog)    // will have same root
        : m_xmlout(xmllog,XMLHandler::pointer)
     {if ((xmllog.empty())||(!xmllog.good()))
         throw(std::invalid_argument("Invalid XML to create LogHelper"));}

    void output(XMLHandler& xmlout) const
     {xmlout=m_xmlout;}

    std::string output() const
     {return m_xmlout.output();}


    void reset(const std::string& roottag)   // clears and resets root tag
     {m_xmlout.set_root(roottag);}


    void put(LogHelper& subhandler)
     {m_xmlout.put_child(subhandler.m_xmlout);}

    void put(XMLHandler& xmlg)
     {m_xmlout.put_child(xmlg);}

    void putEcho(ArgsHandler& argh)
     {XMLHandler xmlecho; argh.echo(xmlecho);
      m_xmlout.put_child(xmlecho);}

    void putEcho(ArgsHandler& argh, const std::string& toptag)
     {XMLHandler xmlecho; argh.echo(xmlecho);
      xmlecho.rename_tag(toptag);
      m_xmlout.put_child(xmlecho);}

    void putInt(const std::string& tagname, const int& val)
     {m_xmlout.put_child(tagname,make_string(val));}

    void putIntVector(const std::string& tagname, const std::vector<int>& val)
     {m_xmlout.put_child(tagname,make_string(val));}

    void putUInt(const std::string& tagname, const uint& val)
     {m_xmlout.put_child(tagname,make_string(val));}

    void putReal(const std::string& tagname, const double& val)
     {m_xmlout.put_child(tagname,make_string(val));}

    void putString(const std::string& tagname, const std::string& val)
     {m_xmlout.put_child(tagname,val);}

           // writes out "true" or "false" for bool
    void putBool(const std::string& tagname, const bool& val)
     {m_xmlout.put_child(tagname,make_string(val));}

           // writes empty tag "<tagname/>" if true, nothing otherwise
    void putBoolAsEmpty(const std::string& tagname, const bool& val)
     {if (val) m_xmlout.put_child(tagname);}

           // puts the XML output for "val" as child
    template <typename T>
    void putItem(const T& val)
    {
     XMLHandler xt; val.output(xt);
     m_xmlout.put_child(xt);
    }

           // puts the XML output for "val" as child inside a tag
           // with name "tagname"
    template <typename T>
    void putItem(const std::string& tagname, const T& val)
    {
     XMLHandler xt; val.output(xt);
     XMLHandler xtt(tagname); xtt.put_child(xt);
     m_xmlout.put_child(xtt);
    }

};


#endif
