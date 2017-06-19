#ifndef XML_HANDLER_H
#define XML_HANDLER_H
#include <iostream>
#include <set>
#include <list>
#include <string>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <stack>
#include <unistd.h>
#include <stdexcept>

// *************************************************************************
// *                                                                       *
// *                      Brief description of XML                         *
// *                                                                       *
// *  Extensible Markup Language (XML) is a markup language that defines   *
// *  a set of rules for encoding documents in a format that is both       *
// *  human-readable and machine-readable. The characters which make up    *
// *  an XML document are divided into "markup" and "content". All strings *
// *  which constitute markup begin with the character < and end           *
// *  with a >. Strings of characters which are not markup are content.    *
// *                                                                       *
// *  A "tag" is a markup construct that begins with < and ends with >     *
// *  and contains a valid tag name and zero or one / characters before    *
// *  or after the tag name.  Tags come in three flavors:                  *
// *        - start-tags; for example: <section>                           *
// *        - end-tags; for example: </section>                            *
// *        - empty-element tags; for example: <line-break/>               *
// *  There must not be any white space between the starting < character   *
// *  and the tag name.  However, white space can occur between the end    *
// *  of the tag name and the ending > character.  There are certain       *
// *  rules for determining if a character string is a valid tag name.     *
// *  Here, we use the following rules:                                    *
// *      -- the first character must be a letter or underscore            *
// *      -- all other characters must be a letter, numerical digit,       *
// *               underscore, period, hyphen, exclamation mark            *
// *      -- first 3 characters cannot be "xml" (case insensitive)         *
// *                                                                       *
// *  An XML "comment" starts with <!-- and ends with -->. These           *
// *  characters and all characters between them are ignored as            *
// *  comments. The string "--" (double-hyphen) is not allowed             *
// *  inside comments.  Comments may not appear inside markup text;        *
// *  they may appear anywhere in a document outside other markup.         *
// *                                                                       *
// *  An "element" is a character string that either begins with a         *
// *  start-tag and ends with a matching end-tag or consists only of       *
// *  an empty-element tag.  The element includes the start and            *
// *  end tags.  The "name" of the element is the name of the start        *
// *  tag, and the name of the end tag must be the same.  The              *
// *  characters between the start- and end-tags, if any, are the          *
// *  element's "content", and may contain text and/or markup,             *
// *  including other elements, which are called child elements.           *
// *  An element is called "simple" if its content is textual only.        *
// *                                                                       *
// *  XML elements can contain one or more "attributes".  Attributes       *
// *  are superfluous, in my opinion. Any information in an attribute      *
// *  can be included as an additional child tag.  Hence, the class        *
// *  defined here will not deal with XML containing attributes.           *
// *  The one exception is in the declaration element.                     *
// *                                                                       *
// *  An XML document may begin by declaring some information about        *
// *  itself, as in the following example:                                 *
// *         <?xml version="1.0" encoding="UTF-8" standalone="yes" ?>      *
// *  This is called a "declaration".  It must have the above form,        *
// *  and it must start at the very first character.                       *
// *                                                                       *
// *  The XML specification defines an XML document as a text that         *
// *  is well-formed, i.e. it satisfies a list of syntax rules provided    *
// *  in the specification. The list is fairly lengthy; some key points    *
// *  are:                                                                 *
// *    -- The begin, end, and empty-element tags that delimit the         *
// *       elements are correctly nested, with none missing and none       *
// *       overlapping.                                                    *
// *    -- The element tags are case-sensitive; the beginning and end      *
// *       tags must match exactly. Tag names must be valid as described   *
// *       above.                                                          *
// *    -- There is a single "root" element that contains all the other    *
// *       elements.                                                       *
// *                                                                       *
// *                                                                       *
// *                                                                       *
// *                       XML Document Object Module                      *
// *                                                                       *
// *  The W3C Document Object Model (DOM) is a platform and                *
// *  language-neutral interface that allows programs and scripts to       *
// *  dynamically access and update the content, structure, and style of   *
// *  an XML document. The XML DOM defines a standard way for accessing    *
// *  and manipulating XML documents. The DOM presents an XML document     *
// *  as a tree-structure.                                                 *
// *                                                                       *
// *  The DOM introduces the concept of an XML "node". An XML document is  *
// *  represented by a tree of nodes.  Each node has one parent node, zero *
// *  or more child nodes, and zero or more sibling nodes.  A node can be  *
// *  an element node or a text node.  For example, the following XML      *
// *  document                                                             *
// *     <root>                                                            *
// *        <tag1>fox <tag2>dog</tag2> cat</tag1>                          *
// *        <tag3>goose</tag3>                                             *
// *        <tag4>lion <tag5/></tag4>                                      *
// *     </root>                                                           *
// *  contains one root node "root" whose parent is NULL and who has no    *
// *  siblings, and three children "tag1", "tag3", "tag4". The node "tag1" *
// *  has parent "root", siblings "tag3" and "tag4", and children nodes    *
// *  "fox", "tag2", "cat".  Note that node "fox" is a text node which     *
// *  has no siblings and no children.  Node "tag2" has parent "tag1" and  *
// *  siblings "fox" and "cat" and one child node "dog".  NOTE: Text is    *
// *  always stored in text nodes. A common error in DOM processing is to  *
// *  navigate to an element node and expect it to contain the text.       *
// *  However, even the simplest element node has a text node under it.    *
// *  For example, in <year>2005</year>, there is an element node (year),  *
// *  and a text node under it, which contains the text (2005).            *
// *                                                                       *
// *************************************************************************



// *******************************************************************
// *                                                                 *
// *  This file defines the class XMLHandler which is used to        *
// *  store and manipulate XML documents.  An XMLHandler is          *
// *  designed like an fstream.  A pointer to a current location     *
// *  in the document is maintained, like a file pointer.  Use       *
// *  "seek" to move around the document, like random access.  You   *
// *  can move around by tag name, like random access, or            *
// *  sequentially by next sibling, first child, parent, etc.  The   *
// *  status of the last seek can be determined with the "good" or   *
// *  "fail" members. Use "get" and "put" to extract and insert      *
// *  textual data.  "get" and "put" are done relative to the        *
// *  current location.  "seek" methods change the current           *
// *  location. Like an fstream, you can turn exceptions on or       *
// *  off.  When on, exceptions are thrown when the state of         *
// *  the object is no longer good.  An XMLHandler in some sense     *
// *  can be viewed as a glorified pointer to XML data.              *
// *                                                                 *
// *  NOTE: XMLHandler does NOT allow attributes in the XML          *
// *  document.  Attributes are superfluous, in my opinion. Any      *
// *  information in an attribute can be included as an additional   *
// *  XML tag.  Hence, to simplify this handler, attributes are      *
// *  disallowed (the one exception is in the XML declaration tag).  *
// *  XMLHandler is not meant to be a full-fledged XML document      *
// *  object module.  It is programmed with simplicity in mind,      *
// *  not efficiency.                                                *
// *                                                                 *
// *  An XMLHandler maintains a pointer to the actual data, which    *
// *  is an object of the class XMLDoc.  XMLDoc objects cannot be    *
// *  created by the end user; XMLDoc objects are only created by    *
// *  XMLHandler objects.  A given XMLDoc object can be shared by    *
// *  several XMLHandler objects.                                    *
// *                                                                 *
// *                                                                 *
// *  Constructors:  There is a default constructor that sets all    *
// *  pointers to null.  Another constructor takes one or two        *
// *  strings and creates a single root tag (empty or simple).       *
// *                                                                 *
// *     XMLHandler xmlh0;   // default constructor                  *
// *     XMLHandler xmlh1("roottag");  // create empty root element  *
// *     XMLHandler xmlh2("tag","text"); // simple root element      *
// *                                                                 *
// *  Various "set" members are available to a previously declared   *
// *  handler.  This allows assignment from an XML string, or a      *
// *  string residing in a file.  Previous content in the handler    *
// *  is cleared.                                                    *
// *                                                                 *
// *     xmlh.set("roottag","text");                                 *
// *     xmlh.set_root("roottag","text"); // same routine            *
// *     xmlh.set_from_string(xmlstring);  // parses string          *
// *     xmlh.set_from_file(filename);                               *
// *                                                                 *
// *  Copy constructor:   The XMLHandler copy constructor takes an   *
// *  optional enum "copymode".  The allowed values of mode below    *
// *  are "subtree_pointer" (default), "subtree_copy", "pointer",    *
// *  or "copy".                                                     *  
// *                                                                 *
// *     XMLHandler(const XMLHandler& xmlin, copymode mode);         *
// *                                                                 *
// *  If "mode" has value "subtree_pointer" or "subtree_copy", the   *
// *  current location of "xmlin" becomes the root element of this   *
// *  handler (and so, the new handler cannot access any of the      *
// *  ancestor nodes).  Otherwise, the root of this handler will be  *
// *  the same as "xmlin".  If "mode" has value "subtree_copy" or    *
// *  "copy", then a completely new deep copy is made and the root   *
// *  element is either the current location or the root of "xmlin"; *
// *  otherwise, this is just a new pointer to the same data in      *
// *  memory and the current element is the same as in "xmlin".      *
// *  When the assignment operator is used, no deep copy is done and *
// *  the entire document is accessible (equivalent to "pointer").   *
// *                                                                 *
// *     XMLHandler xmlh4(xmlh2);                                    *
// *     XMLHandler xmlh5(xmlh2,false,true);                         *
// *     XMLHandler xmlh6;  xmlh6=xml4;                              *
// *                                                                 *
// *  There is also a "set" routine that does a shallow or deep      *
// *  copy of an entire document or a subtree:                       *
// *                                                                 *
// *     xmlh.set(xmlhin, copymode mode);                            *
// *                                                                 *
// *  Nodes can be inserted and appended using "put_child" and       *
// *  "put_sibling" members.  "put_sibling" inserts a new node as    *
// *  the next sibling of the current node, then the current node    *
// *  pointer is updated to the top node of the inserted node(s).    *
// *  Fails if the current node is the root node.  "put_child"       *
// *  inserts as the last child of the current node, but the         *
// *  current node pointer is NOT changed.  Fails if the current     *
// *  node is a text node.  For put operations, failure always       *
// *  throws an exception.                                           *
// *                                                                 *
// *     xmlh.put_sibling(xmlin);  // from root of xmlin             *
// *     xmlh.put_sibling("tagname","textcontent");                  *
// *     xmlh.put_sibling_text_node("text");                         *
// *                                                                 *
// *     xmlh.put_child(xmlin);  // from root of xmlin               *
// *     xmlh.put_child("tagname","textcontent");                    *
// *     xmlh.put_child_text_node("text");                           *
// *                                                                 *
// *  Various "seek" members can be used for moving around the       *
// *  document.  You can seek incrementally by nodes, or by          *
// *  searching for a particular tag name:                           *
// *                                                                 *
// *     xmlh.seek_root();                                           *
// *     xmlh.seek_next_sibling();  // fails if no more siblings     *
// *     xmlh.seek_first_child();   // fails if no children          *
// *     xmlh.seek_child("tagname");                                 *
// *     xmlh.seek_parent();                                         *
// *     xmlh.seek_next_node();  // children first, then siblings    *
// *     xmlh.seek_next("tagname"); // children first, then siblings *
// *     xmlh.seek_unique("tagname"); // fails if absent,not unique  *
// *     xmlh.find("tagname"); // return list<XMLHandler>            *
// *     xmlh.count("tagname");  // returns integer                  *
// *                                                                 *
// *  Some useful routines:                                          *
// *                                                                 *
// *     xml_root_assert(xmlh,"Name");  // check root tag is "Name"  *  
// *     xmlreadchild(xmlh,"UniqueName",val);                        *
// *                                                                 *
// *                                                                 *
// *******************************************************************




class XMLHandler;   // forward declaration



   // XMLDoc is a tree of XMLNode objects.
   // Each node has a name string, a boolean indicating if it
   // is a text node or not, and pointers to its parent node,
   // first child node, and its next sibling node.  A text node
   // cannot have a child node, and the root node cannot have
   // any siblings.  A node without a child is either a text node
   // or an empty tag. Comments in the XML string are ignored.  
   // Any declaration tag is kept.


class XMLDoc
{

 public:
   
   struct XMLNode
   {
      std::string name;
      bool text;
      XMLNode *parent, *firstchild, *nextsibling;
      XMLNode(const std::string& inname, XMLNode *inparent, bool is_text=false)
                  : name(inname), text(is_text), parent(inparent), 
                    firstchild(0), nextsibling(0) {}
      ~XMLNode(){}
    private:
      XMLNode(const XMLNode&);  // disable copying
      XMLNode& operator=(const XMLNode&);
      XMLNode& operator=(XMLNode&);
   };


         // constructors and destructors

   explicit XMLDoc(const std::string& xmlstr);  // parse string

   explicit XMLDoc(const std::string& tagname, 
                   const std::string& text_content);  // single simple tag

   XMLDoc(XMLNode* top, const std::string& indecl);

   XMLDoc(const XMLDoc& xmldoc);             // undefined to prevent copying
   XMLDoc& operator=(const XMLDoc& xmldoc);  // undefined to prevent copying

   ~XMLDoc();

   //   <?   =  startdecl     ?>  =  enddecl    <     =  starttag     >  =  endtag
   //   </   =  startclose    />  =  endempty   <!--  =  startcom   -->  =  endcom

   enum XMLEvent {  startdecl, enddecl, starttag, endtag,
                    startclose, endempty, startcom, endcom, none };

       //  creation members

   XMLNode* create_by_parsing(const std::string& xmlstr, size_t start, 
                              size_t stop);
   XMLNode* create_simple_tag(const std::string& tagname, const std::string& text,
                              XMLNode* parent);
   XMLNode* create_text_node(const std::string& text, XMLNode *parent);
   XMLNode* create_copy(XMLNode* top);

       //  destruction members

   void destroy(XMLNode*& top);  // top and descendents destroyed
   void clear();                 // deallocate all memory

       //  connect/disconnect members

   void connect_as_lastchild(XMLNode *into, XMLNode* add);
   void connect_as_nextsibling(XMLNode *into, XMLNode* add);
   void disconnect(XMLNode *top);


       //  utility routines

   std::string get_tag_name(const std::string& instr, size_t charstart, 
                            size_t charstop);
   std::string trim(const std::string& str);   
   bool nonwhitespace(const std::string& instr, size_t charstart, size_t charstop);
   void find_next_xml_event(const std::string& xmlstr, size_t start, 
                            size_t in_stop, XMLEvent& type,  size_t& pos, 
                            bool incomment=false);
   void find_next_xml_tag_event(const std::string& xmlstr, size_t start, size_t stop,
                                XMLEvent& type, size_t& pos, std::string& textpassed);
   void find_root_start_tag(const std::string& xmlstr, size_t start, size_t in_stop,
                            size_t& pos);
   XMLNode* parse_root_tag(const std::string& xmlstr, size_t& start, size_t stop, 
                           XMLEvent& type, size_t& pos, 
                           std::stack<XMLDoc::XMLNode*>& ancestors,
                           std::stack<XMLDoc::XMLNode*>& lastchildren);
   void parse_to_next_event(const std::string& xmlstr, size_t& curr, size_t stop, 
                            XMLEvent& lastevent, std::stack<XMLNode*>& ancestors,
                            std::stack<XMLNode*>& lastchildren);
   size_t get_declaration(const std::string& instr);
   bool get_attribute(const std::string& instr, const std::string& attrkey,
                      std::string& attrvalue, size_t& start, size_t& stop);
   std::ostream& output(std::ostream& os, XMLNode* top, int indent=0) const;
   void delete_node(XMLNode*& node);
   void add_reference(XMLHandler *ref);
   bool remove_reference(XMLHandler *ref);

   friend class XMLHandler;
   friend class XMLContentComparer;

      //  data members

   XMLNode *root;
   std::string declaration;
   std::set< XMLHandler* > refset; // the XMLHandler objects that point to this

};



// *********************************************************


class XMLHandler
{

   XMLDoc *content;              // pointer to XML data
   XMLDoc::XMLNode *root;        // root element for this handler
   XMLDoc::XMLNode *current;     // current location in document
   bool exceptions;

 public:

   enum copymode { subtree_pointer, subtree_copy, pointer, copy };
  
   explicit XMLHandler(int excpts=1) 
        : content(0), root(0), current(0), exceptions(bool(excpts)) {}

   XMLHandler(const std::string& roottagname, 
              const std::string& content=std::string());

   XMLHandler(const char* roottagname, const char* rootcontent);

   explicit XMLHandler(const char* roottagname);

       // Copy constructor. If "mode" has value "subtree_pointer" or
       // "subtree_copy", the current location of "xmlin" becomes the root 
       // element of this handler (cannot access ancestor nodes).  Otherwise,
       // the root of this handler will be the same as "xmlin".  If "mode" has
       // value "subtree_copy" or "copy", then a completely new deep copy is 
       // made and the root element is either the current location or the root
       // of "xmlin"; otherwise, this is just a new pointer to the same data 
       // in memory and the current element is the same as in "xmlin".

   XMLHandler(const XMLHandler& xmlin, copymode mode=subtree_pointer);

       // Copy constructor.  First does a "seek_unique" looking for
       // the tag with name specified in "newroottag".  Fails
       // if not found or multiple instances found.  This creates a new
       // pointer to a subtree of the same XML document.

   XMLHandler(const XMLHandler& xmlin, const std::string& newroottag,
              bool seek_to_child_only=false);

   XMLHandler(const XMLHandler& xmlin, const char* newroottag);


          // no full copy, entire document retained
   XMLHandler& operator=(const XMLHandler& xmlin); 

   XMLHandler& set(const std::string& roottagname, 
                   const std::string& rootcontent=std::string());

   XMLHandler& set(const char* roottagname, const char* rootcontent);

   XMLHandler& set(const char* roottagname);

   XMLHandler& set_root(const std::string& roottagname, 
                        const std::string& rootcontent=std::string());

   XMLHandler& set_root(const char* roottagname, const char* rootcontent);

   XMLHandler& set_root(const char* roottagname);

   XMLHandler& set_from_string(const std::string& xmlin);

   XMLHandler& set_from_file(const std::string& filename);

   XMLHandler& set(const XMLHandler& xmlin, copymode mode=subtree_pointer);

   ~XMLHandler();

   XMLHandler& clear();


   std::string output(int indent=0) const;  // formatted with newlines, indentations

   std::string str() const;   // just a simple string, no newlines nor indentations

   std::string output_current(int indent=0) const;


   bool good( ) const {return (current!=0);}

   bool fail( ) const {return (current==0);}


   void set_exceptions_on() {exceptions=true;}

   void set_exceptions_off() {exceptions=false;}

   void set_exceptions(bool state) {exceptions=state;}

   bool get_exceptions() const {return exceptions;}

   bool empty() const {return (content==0);}

             //  seek members (if seek fails, state of handler
             //  is set to "bad", and if exceptions are turned
             //  on, a string exception is thrown)

   void seek_root();

   void seek_next_sibling();

   void seek_first_child();

   void seek_child(const std::string& tagname);  // must be unique

   void seek_parent();

   void seek_next_node();

   void seek_next(const std::string& tagname);

   void seek_unique(const std::string& tagname);  // starting root and all descendents

   void seek_unique_to_child(const std::string& tagname);  // seek in root and children only

   void seek_unique_child(const std::string& tagname);  // seek among children only (not root)

   std::list<XMLHandler> find(const std::string& tagname) const; // starting from root

   std::list<XMLHandler> find(const std::list<std::string>& tagnames) const;

   int count(const std::string& tagname) const;  // starting from root

   std::list<XMLHandler> find_among_children(const std::string& tagname) const; 

   std::list<XMLHandler> find_among_children(const std::list<std::string>& tagnames) const;

   int count_among_children(const std::string& tagname) const;

   int count_to_among_children(const std::string& tagname) const;  // includes current

   int count_children() const;

       // returns false if "tagname" not current or among children of current,
       // true if one occurrence, throws exception if multiple occurrences

   bool query_unique_to_among_children(const std::string& tagname) const; 

             // various get members

   std::string get_node_name() const;

   std::string get_node_value() const;

   bool is_text_valued() const;

   std::string get_tag_name() const;     // returns empty string if invalid

   std::string get_text_content() const; // returns empty string if invalid

   bool is_simple_element() const;

   bool is_empty_tag() const;


             // Insertion. "put_sibling" inserts as the
             // next sibling of the current node, then the
             // current node pointer is updated to the top
             // node of the inserted node(s).  Fails if the
             // current node is the root node. Exception thrown
             // if failure occurs.  All XMLHanders pointing to
             // the same content are changed.

   void put_sibling(const XMLHandler& xmlin);  // from root of xmlin

   void put_sibling(const std::string& tagname);

   void put_sibling(const std::string& tagname, 
                    const std::string& text_content);

   void put_sibling_text_node(const std::string& text);


             // Insertion. "put_child" inserts as the
             // last child of the current node, but the
             // current node pointer is NOT changed.  Fails
             // if the current node is a text node.  Exception 
             // thrown if failure occurs. All XMLHanders pointing
             // to the same content are changed.

   void put_child(const XMLHandler& xmlin); // from root of xmlin

   void put_child(const std::string& tagname);

   void put_child(const std::string& tagname, 
                  const std::string& text_content);

   void put_child_text_node(const std::string& text);


             // Inserts or changes the textual content
             // of the current node (if a text node).
             // Exception throw if failure occurs.  All XMLHanders
             // pointing to the same content are changed.

   void set_text_content(const std::string& in_text_content);


             // Renames the current tag.
             // Exception throw if failure occurs.  All XMLHanders
             // pointing to the same content are changed.

   void rename_tag(const std::string& new_tag_name);

             // Erases current element (if not a text node).
             // Exception thrown if failure occurs.  All XMLHanders 
             // pointing to the same content are changed.

   void erase_current_element();


             // Erases current node if a text node.
             // Exception thrown if failure occurs. All XMLHanders 
             // pointing to the same content are changed.

   void erase_text_node();

             // Compares the XML content in current handler with that in 
             // "xmlh_cmp".  In doing the comparison, the
             // order of sibling nodes is irrelevant, and textual content
             // is compared token by token.  Integer tokens are compared as
             // integers, and floating-point tokens are compared as 
             // floats.  If the difference between floats is less than
             // "float_rel_tol", they are considered the same.

   bool IsEqualContent(XMLHandler& xmlh_cmp, float float_rel_tol=1e-4);


 private:
  
   void do_set(const std::string& roottagname, 
               const std::string& rootcontent);
   void do_set(const XMLHandler& xmlin, bool subtree=true, 
               bool makecopy=false);
   void do_set_from_string(const std::string& xmlin);
   bool seek_next_sib_or_parent();
   std::string empty_string() const;

   friend class XMLDoc;
   friend class XMLContentComparer;

};


        //  class used by XMLHandler for content comparison


class XMLContentComparer
{

   float m_float_rel_tol;
   XMLContentComparer(float tol=1e-4) : m_float_rel_tol(std::abs(tol)) {}

 public:
   XMLContentComparer(const XMLContentComparer& in) : m_float_rel_tol(in.m_float_rel_tol) {}
   XMLContentComparer& operator=(const XMLContentComparer& in) 
      {m_float_rel_tol=in.m_float_rel_tol; return *this;}
   ~XMLContentComparer() {}

 private:
   typedef XMLDoc::XMLNode*  xmlNodePtr;
   bool IsEqual(XMLHandler& xmlh1, XMLHandler& xmlh2);
   void get_next_token(const char*& start, int& nchar);
   int token_strcmp(const char *a, int na, const char *b, int nb);
   bool is_integer_token(const char *const a, int nchar, int& ivalue);
   bool is_float_token(const char *const a, int nchar, float& fvalue);
   int token_content_cmp(const char *a, int na, const char *b, int nb);
   int content_cmp(const char *a, const char *b);
   int sibling_node_cmp(xmlNodePtr a_node, xmlNodePtr b_node);
   int node_cmp(xmlNodePtr a_node, xmlNodePtr b_node);

   friend class XMLDoc;
   friend class XMLHandler;

 public:
   bool operator()(xmlNodePtr a_node, xmlNodePtr b_node); // less than op

};


 // ******************************************************************

         //  Helper routines

inline std::ostream& operator<<(std::ostream& os, const XMLHandler& xmlh)
{
 os << xmlh.output();
 return os;
}


// *********************************************************

      //  Extract a value from a string, throwing an exception on failure.

template <typename T>
void primitive_extract_from_string(const std::string& pstring, T& result, 
                                   const char* ptype) 
{
 try{
    std::istringstream is(pstring);
    std::ios::iostate state=is.exceptions();
    is.exceptions(std::ios::failbit);
    is.setf(std::ios_base::boolalpha);
       // try to read the type from the istringstream. 
       //   bool's should be "true" or "false" strings 
    is  >> result;
       // turn off exceptions on failure
    is.exceptions(state);
       // look for non white spaces in remainder of stream
    char c;
    while (is.get(c),!is.eof()){
       if (!isspace(c)) throw(std::invalid_argument("err"));}}
 catch(const std::exception& xp){ 
    std::ostringstream err;
    err << "Error: Failed to convert string: \"" 
        << pstring << "\" to " << ptype;
    throw(std::invalid_argument(err.str()));}
}

inline void extract_from_string(const std::string& pstring, std::string& result) 
{ result=pstring;}

inline void extract_from_string(const std::string& pstring, int& result) 
{ primitive_extract_from_string<int>(pstring, result, "int");}

inline void extract_from_string(const std::string& pstring, unsigned int& result) 
{ primitive_extract_from_string<unsigned int>(pstring, result, "unsigned int");}
    
inline void extract_from_string(const std::string& pstring, short int& result) 
{ primitive_extract_from_string<short int>(pstring, result,"short int");}

inline void extract_from_string(const std::string& pstring, unsigned short int& result) 
{ primitive_extract_from_string<unsigned short int>(pstring, result, "unsigned short int");}

inline void extract_from_string(const std::string& pstring, long int& result) 
{ primitive_extract_from_string<long int>(pstring, result, "long int");}
    
inline void extract_from_string(const std::string& pstring, unsigned long int& result) 
{ primitive_extract_from_string<unsigned long int>(pstring, result, "unsigned long int");}
    
inline void extract_from_string(const std::string& pstring, float& result) 
{ primitive_extract_from_string<float>(pstring, result, "float");}

inline void extract_from_string(const std::string& pstring, double& result) 
{ primitive_extract_from_string<double>(pstring, result, "double");}

inline void extract_from_string(const std::string& pstring, bool& result) 
{ primitive_extract_from_string<bool>(pstring, result, "bool");}




template <typename T>
void primitive_vector_extract(const std::string& pstring, std::vector<T>& result, 
                              const char* ptype) 
{
 try{
    result.clear();
    std::istringstream is(pstring);
    is.exceptions(std::ios::badbit);
    is.setf(std::ios_base::boolalpha);
    T tokenvalue;
    while (is >> tokenvalue) result.push_back(tokenvalue);
    if ((!is.eof()) && is.fail()){
       throw(std::invalid_argument(std::string("Error in reading vector ")+pstring));}}
 catch(const std::exception& xp){ 
    std::ostringstream err;
    err << "Error: Failed to convert string: \"" 
        << pstring << "\" to vector of " << ptype;
    throw(std::invalid_argument(err.str()));}
}


inline void extract_from_string(const std::string& pstring, std::vector<int>& result)
{ primitive_vector_extract(pstring, result, "int"); }

inline void extract_from_string(const std::string& pstring, std::vector<unsigned int>& result)
{ primitive_vector_extract(pstring, result, "unsigned int"); }

inline void extract_from_string(const std::string& pstring, std::vector<short int>& result)
{ primitive_vector_extract(pstring, result, "short int"); }

inline void extract_from_string(const std::string& pstring, std::vector<unsigned short int>& result)
{ primitive_vector_extract(pstring, result, "unsigned short int"); }

inline void extract_from_string(const std::string& pstring, std::vector<long int>& result)
{ primitive_vector_extract(pstring, result, "long int"); }

inline void extract_from_string(const std::string& pstring, std::vector<unsigned long int>& result)
{ primitive_vector_extract(pstring, result, "unsigned long int"); }

inline void extract_from_string(const std::string& pstring, std::vector<float>& result)
{ primitive_vector_extract(pstring, result, "float"); }

inline void extract_from_string(const std::string& pstring, std::vector<double>& result)
{ primitive_vector_extract(pstring, result, "double"); }

inline void extract_from_string(const std::string& pstring, std::vector<bool>& result)
{ primitive_vector_extract(pstring, result, "bool"); }



// **********************************************************************


      //  Conversion of numerical data to string.

template <typename T>
std::string numerical_to_string(const T& data) 
{
 try{
    std::ostringstream oss;
    oss.setf(std::ios_base::boolalpha);
    oss.precision(12);
    oss << data;
    return oss.str();}
 catch(const std::exception& xp){ 
    std::ostringstream err;
    err << "Error: Failed to convert numerical to string: data = "<<data;
    throw(std::invalid_argument(err.str()));}
}


inline std::string make_string(const int& data) 
{
 return numerical_to_string<int>(data);
}

inline std::string make_string(const unsigned int& data) 
{
 return numerical_to_string<unsigned int>(data);
}
    
inline std::string make_string(const short int& data) 
{
 return numerical_to_string<short int>(data);
}

inline std::string make_string(const unsigned short int& data) 
{
 return numerical_to_string<unsigned short int>(data);
}

inline std::string make_string(const long int& data) 
{
 return numerical_to_string<long int>(data);
}
    
inline std::string make_string(const unsigned long int& data) 
{
 return numerical_to_string<unsigned long int>(data);
}
    
inline std::string make_string(const float& data) 
{
 return numerical_to_string<float>(data);
}

inline std::string make_string(const double& data) 
{
 return numerical_to_string<double>(data);
}


inline std::string make_string(const bool& data) 
{
 return numerical_to_string<bool>(data);
}

inline std::string make_string(const std::complex<double>& data) 
{
 return numerical_to_string<std::complex<double> >(data);
}

inline std::string make_string(const std::complex<float>& data) 
{
 return numerical_to_string<std::complex<float> >(data);
}


template <typename T>
std::string make_string(const std::vector<T>& datavec)
{
 std::string tmp;
 if (datavec.size()==0) return tmp;
 tmp=make_string(datavec[0]);
 for (typename std::vector<T>::size_type k=1;k<datavec.size();++k)
    tmp+=" "+make_string(datavec[k]);
 return tmp; 
}


inline std::string scientific_to_string(const float& data) 
{
 std::ostringstream oss;
 oss.precision(10);
 oss.setf(std::ios::scientific);
 oss << data;
 return oss.str();
}

inline std::string scientific_to_string(const double& data) 
{
 std::ostringstream oss;
 oss.precision(18);
 oss.setf(std::ios::scientific);
 oss << data;
 return oss.str();
}

// ******************************************************************

  // Removes tabs and newline characters, then trims
  // leading and trailing blanks.

std::string tidyString(const std::string& str);    


  // first tidies the string in "str", then removes any
  // leading path/subdirectory information in the file name

std::string tidyFileName(const std::string& str);


  // Trims leading and trailing white space, then checks
  // to make sure each character is alphanumeric, underscore, 
  // or period. If name is invalid, an empty string is returned.

std::string tidyName(const std::string& str);


  // converts an integer to a string

std::string int_to_string(int ival);



// *********************************************************

    // This routine calls an xml read inside a try block
    // and outputs an informative message if the read fails.
    // The added parameter is a string which should be the
    // name of the class which called this function.
    // Note that a **tag name** should be used.

template <typename T>
void xmlread(XMLHandler& xmlin, const std::string& tagname, T& val,
             const std::string& callingClass)
{
 XMLHandler xmlh(xmlin);
 bool state=xmlh.get_exceptions();
 try{
    xmlh.set_exceptions_on();
    xmlh.seek_unique(tagname);
    std::string content=xmlh.get_text_content();
    extract_from_string(content,val);
    xmlh.set_exceptions(state);}
 catch(const std::exception& err_msg){
    xmlh.set_exceptions(state);
//    std::cerr <<"Invalid read of "<<tagname<<" in "<<callingClass
//              <<":  "<<err_msg<<std::endl;
    throw;}
}

template <typename T>
void xmlreadchild(XMLHandler& xmlin, const std::string& tagname, T& val,
                  const std::string& callingClass)
{
 XMLHandler xmlh(xmlin);
 bool state=xmlh.get_exceptions();
 try{
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content=xmlh.get_text_content();
    extract_from_string(content,val);
    xmlh.set_exceptions(state);}
 catch(const std::exception& err_msg){
    xmlh.set_exceptions(state);
//    std::cerr <<"Invalid read of "<<tagname<<" in "<<callingClass
//              <<":  "<<err_msg<<std::endl;
    throw;}
}

template <typename T>
void xmlreadchild(XMLHandler& xmlin, const std::string& tagname, T& val)
{
 XMLHandler xmlh(xmlin);
 bool state=xmlh.get_exceptions();
 try{
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content=xmlh.get_text_content();
    extract_from_string(content,val);
    xmlh.set_exceptions(state);}
 catch(const std::exception& err_msg){
    xmlh.set_exceptions(state);
    throw;}
}

// *********************************************************

    // This returns the number of times that the tag "tagname"
    // is found in the descendents of the current context.
    // A **tag name** should be input.

int xml_tag_count(XMLHandler& xmlh, const std::string& tagname);

int xml_child_tag_count(XMLHandler& xmlh, const std::string& tagname);

    // Check for existence of a unique tag "tagname" and if not
    // present, throws exception.  The found tag becomes current location.

void xml_tag_assert(XMLHandler& xmlh, const std::string& tagname);

void xml_tag_assert(XMLHandler& xmlh, const std::string& tagname,
                    const std::string& infoname);

    // Check for existence of a unique tag "tagname1" or "tagname2" 
    // and if neither present, throw exception.
  
void xml_either_tag_assert(XMLHandler& xmlh, const std::string& tagname1,
                           const std::string& tagname2);


    // Check for existence of a unique tag "tagname" as child of
    // current node, and if not present, throws exception.
    // The found tag becomes current location.

void xml_child_assert(XMLHandler& xmlh, const std::string& tagname,
                      const std::string& infoname);

void xml_child_assert(XMLHandler& xmlh, const std::string& tagname);

    // Check that root tag is "tagname", throws exception if not.

void xml_root_assert(XMLHandler& xmlh, const std::string& tagname,
                     const std::string& infoname);

void xml_root_assert(XMLHandler& xmlh, const std::string& tagname);



void xml_cerr(XMLHandler& xmlh, const std::string& msg);


// *********************************************************

template <typename T>
bool xmlreadif(XMLHandler& xmlh, const std::string& tagname, T& val,
               const std::string& callingClass)
{
 if (xml_tag_count(xmlh,tagname)==1){
    xmlread(xmlh,tagname,val,callingClass);
    return true;}
 return false;
}

template <typename T>
bool xmlreadifchild(XMLHandler& xmlh, const std::string& tagname, T& val)
{
 try{
    xmlreadchild(xmlh,tagname,val);
    return true;}
 catch(const std::exception& xp){
    return false;}
}

// *********************************************************

template <typename T>
void assertEqual(const T& obj1, const T& obj2, const std::string& callingClass)
{
 try{
    obj1.checkEqual(obj2);}
 catch(const std::exception& err_msg){
    std::cerr << err_msg.what() <<" in "<<callingClass<<std::endl;
    exit(1);}
}
 
// *********************************************************

   // outputs the current context for a failed xml read
   // (typically in an Info constructor)

void xmlreadfail(XMLHandler& xmlh, const std::string& infoname);


// *********************************************************



  // Compares the XML content in "doc1" and "doc2" and returns
  // "true" if they are the same.  In doing the comparison, the
  // order of sibling nodes is irrelevant, and textual content
  // is compared token by token.  Integer tokens are compared as
  // integers, and floating-point tokens are compared as 
  // floats.  If the difference between floats is less than
  // "float_rel_tol", they are considered the same.


bool xmlContentIsEqual(const std::string& doc1, const std::string& doc2, 
                       float float_rel_tol = 1e-6);


  // Same as above, but applied to XMLHandlers "xmlh1" and "xmlh2".

bool xmlContentIsEqual(XMLHandler& xmlh1, XMLHandler& xmlh2, 
                       float float_rel_tol = 1e-6);

  // Same as above, but first tries a straight string comparison.
  // If strings are the same, returns true; otherwise, an
  // XML content comparison is made.

bool headerMatch(const std::string& doc1, const std::string& doc2, 
                 float float_rel_tol = 1e-6);


// *********************************************************

    //  tests if file having name "file_name" exists

bool fileExists(const std::string& file_name);


// *********************************************************
#endif
