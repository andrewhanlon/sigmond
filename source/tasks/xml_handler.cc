#include "xml_handler.h"
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cstring>
using namespace std;


 // *****************************************************

    //  Constructor from a string.  First extract the XML 
    //  declaration, if there is one.  Then parse the XML, 
    //  then check that there is just one top or root node. 
    //  Throw exception if anything in the string is invalid.

XMLDoc::XMLDoc(const string& xmlstr)
{
 root=0;
 size_t start=get_declaration(xmlstr);
 try{
    root=create_by_parsing(xmlstr,start,xmlstr.length());}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("Invalid XML"));} 
}

    //  Constructor as a single simple element.

XMLDoc::XMLDoc(const string& tagname, const string& text_content)
{
 root=0;
 try{
    root=create_simple_tag(tagname,text_content,0);}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("Invalid XML root tag"));}
}


     //  Copy constructor.  Copies only those nodes that
     //  are descendents from "top".

XMLDoc::XMLDoc(XMLNode* top, const string& indecl) : declaration(indecl)
{
 root=0;
 if (top==0) throw(std::invalid_argument("Subtree is not valid XML doc on its own"));
 if (top->text) throw(std::invalid_argument("Subtree is not valid XML doc on its own"));
 root=create_copy(top);
}


      //  Actual routines that do the memory allocation

XMLDoc::XMLNode* XMLDoc::create_by_parsing(const string& xmlstr, 
                                           size_t start, size_t stop)
{
 XMLNode *top=0;
 try{
    stack<XMLDoc::XMLNode*> ancestors;
    stack<XMLDoc::XMLNode*> lastchildren;
    XMLEvent lastevent;
    size_t curr;
    top=parse_root_tag(xmlstr,start,stop,lastevent,curr,ancestors,lastchildren);
    while (lastevent!=none)
        parse_to_next_event(xmlstr,curr,stop,lastevent,ancestors,lastchildren);
    if (!ancestors.empty()) throw(std::invalid_argument("Invalid XML"));
    return top;}
 catch(const std::exception& xp){
    destroy(top); throw(std::invalid_argument("Invalid XML or allocation error"));} 
}

  // exception thrown if "tagname" is invalid

XMLDoc::XMLNode* XMLDoc::create_simple_tag(const string& tagname, 
                                           const std::string& in_text,
                                           XMLDoc::XMLNode* parent)
{
 XMLNode *top=0;
 string tag(get_tag_name(tagname,0,tagname.length()));
 top=new XMLNode(tag,parent);
 string text(trim(in_text));
 if (!text.empty())
    top->firstchild=new XMLNode(text,top,true);
 return top;
}

  // returns pointer to text node; null pointer if text is empty

XMLDoc::XMLNode* XMLDoc::create_text_node(const string& in_text, 
                                          XMLDoc::XMLNode *parent)
{
 string text(trim(in_text));
 XMLNode *top=0;
 if (!text.empty())
    top=new XMLNode(text,parent,true);
 return top;
}


          // copies top and its descendents

XMLDoc::XMLNode* XMLDoc::create_copy(XMLDoc::XMLNode* top)
{
 if (top==0){
    throw(std::invalid_argument("XMLDoc create copy failed: nothing to copy"));}
 XMLNode *curr=top;
 XMLNode *newtop=0;
 newtop=new XMLNode(top->name,top->parent,top->text);
 XMLNode *nodecopy=newtop;
 bool ascending=false;
 bool done=(top->firstchild==0);  // if empty root, done
 while (!done){
    if ((curr->firstchild!=0)&&(!ascending)){
       curr=curr->firstchild;
       nodecopy->firstchild=new XMLNode(curr->name,nodecopy,curr->text);
       nodecopy=nodecopy->firstchild;}
    else if (curr->nextsibling!=0){
       ascending=false;
       curr=curr->nextsibling;
       nodecopy->nextsibling=new XMLNode(curr->name,nodecopy->parent,curr->text);
       nodecopy=nodecopy->nextsibling;}
    else{
       curr=curr->parent;
       nodecopy=nodecopy->parent; 
       ascending=true;
       if (curr==top) done=true;}
    }
 return newtop;
}

     // top and all descendents pruned from tree;
     // if the "current" pointers of any of the XMLHandler
     // objects point to one of the deleted nodes, those
     // current pointers must be nullified!

void XMLDoc::destroy(XMLNode*& top)
{
 if (top==0) return;
 XMLNode *target;
 XMLNode *curr=top->firstchild;
 bool ascending=false;
 while (curr){
    target=curr;
    if ((!ascending)&&(curr->firstchild!=0)){
       curr=curr->firstchild;}
    else if (curr->nextsibling!=0){
       ascending=false;
       curr=curr->nextsibling;
       delete_node(target);}
    else{
       curr=curr->parent;
       ascending=true;
       delete_node(target);
       if (curr==top) curr=0;}
    }
 delete_node(top);
 top=0;
}

  // This deletes a single node, updating the "current" pointers
  // of any XMLHandlers that point to this node.  This is only meant
  // to be used by the "destroy" member.

void XMLDoc::delete_node(XMLNode*& node)
{
    // before removing node, find all XMLHandlers that point to this
    // document, nullify all "current" pointers that point to this node.
    // If the root pointer of an XMLHander points to this node, then
    // that entire XMLHandler must be nullified and removed from "refset".
 set<XMLHandler*> eraseset;  // keeps track of XMLHandler to remove from "refset"
 for (set< XMLHandler* >::iterator it=refset.begin();it!=refset.end();++it){
    if ((*it)->current==node){
       (*it)->current=0;}
    if ((*it)->root==node){
       (*it)->root=0; (*it)->content=0;
       eraseset.insert(*it);}}
   // now remove elements from "refset"
 for (set<XMLHandler*>::iterator it=eraseset.begin();it!=eraseset.end();++it)
    refset.erase(*it);
   // finally it is safe to delete the node
 delete node;
 node=0;
}


void XMLDoc::clear()
{
 destroy(root);
}


XMLDoc::~XMLDoc()
{
 clear();
}

     //  "add" below must have no siblings (but children are okay)

void XMLDoc::connect_as_lastchild(XMLNode *into, XMLNode* add)
{
 if (add==0) return;
 if ((into==0)||(add->nextsibling!=0))
    throw(std::invalid_argument("connection cannot be made in XMLDoc"));
 if (into->firstchild==0){
    into->firstchild=add;}
 else{
    XMLNode *x=into->firstchild; 
    while (x->nextsibling!=0) x=x->nextsibling;
    x->nextsibling=add;}
 add->parent=into;
}

     //  "add" below must have no siblings (but children are okay)

void XMLDoc::connect_as_nextsibling(XMLNode *into, XMLNode* add)
{
 if (add==0) return;
 if ((into==0)||(add->nextsibling!=0))
    throw(std::invalid_argument("connection cannot be made in XMLDoc"));
 XMLNode *nxsb=into->nextsibling;
 into->nextsibling=add;
 add->nextsibling=nxsb;
 add->parent=into->parent;
}


void XMLDoc::disconnect(XMLNode *top)
{
 if ((top==0)||(top->parent==0)) return;
 XMLNode *parent=top->parent;
 if (parent->firstchild==top){
    parent->firstchild=top->nextsibling;}
 else{
    XMLNode *x=parent->firstchild;
    while (x->nextsibling!=top) x=x->nextsibling;
    x->nextsibling=top->nextsibling;}
 top->parent=0;
 top->nextsibling=0;
}


    // This routine takes the characters in "instr", starting at
    // position "charstart" and stopping at "charstop" (does not
    // include this last character) and returns the XML tag name.
    // An exception is thrown if the tag name is invalid.  There
    // can be NO leading blanks, but trailing blanks are allowed,
    // but removed.

string XMLDoc::get_tag_name(const string& instr, size_t charstart, 
                            size_t charstop) 
{
 if (charstart>=charstop) 
    throw(std::invalid_argument("Invalid XML tag name"));
 string tagName(instr,charstart,charstop-charstart);
    // remove any trailing blanks
 size_t pos=tagName.find_last_not_of(" ");
 if (pos<(tagName.length()-1))
    tagName.erase(pos+1,tagName.length()-pos-1);
   // check valid XML tag name
   // no other blanks allowed
 pos=tagName.find_first_of(" ");
 bool success=true;
 if ((pos!=string::npos)||(tagName.empty()))
    success=false;
   // first 3 characters cannot be xml (case insensitive)
 if ((success)&&(tagName.length()>=3)){
    if (((tagName[0]=='x')||(tagName[0]=='X'))
      &&((tagName[1]=='m')||(tagName[1]=='M'))
      &&((tagName[2]=='l')||(tagName[2]=='L')))
        success=false;}
   // first character must be letter or underscore
 if (success){
    char c=tagName[0];
    if (!(((c>='a')&&(c<='z'))||((c>='A')&&(c<='Z'))||(c=='_')))
       success=false;
   // all other characters must be letter, numerical digit,
   // underscore, period, hyphen, exclamation mark
    unsigned int k=1;
    while ((success)&&(k<tagName.length())){
       c=tagName[k];
       if (!(((c>='a')&&(c<='z'))||((c>='A')&&(c<='Z'))
           ||((c>='0')&&(c<='9'))||(c=='_')||(c=='-')||(c=='.')||(c=='!')))
          success=false;
       k++;}}
 if (!success)
    throw(std::invalid_argument(string("Invalid XML tag name: <")+tagName+string(">")));
 return tagName;
}


  // removes tabs, newline, linefeed characters, then trims
  // leading and trailing blanks.

string XMLDoc::trim(const string& str)   
{
 string tmp;
 for (size_t i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 size_t start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 size_t len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}

   // checks the characters in "instr" starting at location "charstart"
   // and stopping at location "charstop" (not inclusive) and returns
   // "true" if any of them are NOT "white space" (blanks, line feeds, 
   // tabs, new line, etc.)

bool XMLDoc::nonwhitespace(const string& instr, size_t charstart, size_t charstop)
{
 bool nonwhite=false;
 size_t stop=(charstop<instr.length())?charstop:instr.length();
 for (size_t k=charstart;k<stop;++k){
    char c=instr[k];
    if ((c!=' ')&&(c!='\n')&&(c!='\t')&&(c!='\r')&&(c!='\v')&&(c!='\f')){
       nonwhite=true; break;}}
 return nonwhite;
}

    // Searches the string in "xmlstr" starting at character location "start"
    // and stopping at "in_stop" (not inclusive) for the first occurrence of
    // an XML event (see enum XMLEvent). It returns the type of event found
    // in "type" and its location in "pos" (the location of the "<" or ">").
    // Note: if "incomment" is true, then only "-->" will be searched for.
    // If no event is found, the type "none" is returned and pos will be
    // set to "in_stop".  Note that if invalid XML is encountered, a string
    // exception will be thrown.

void XMLDoc::find_next_xml_event(const string& xmlstr, size_t start, size_t in_stop,
                                 XMLDoc::XMLEvent& type, size_t& pos, bool incomment)
{
 size_t stop=min(in_stop,xmlstr.length());
 type=none; pos=stop;
 if (start>=stop) return;
 
 if (incomment){
    pos=xmlstr.find("--",start);
    if ((pos==string::npos)||(pos>(in_stop-3))) throw(std::invalid_argument("Invalid XML"));
    pos+=2; 
    if (xmlstr[pos]!='>') throw(std::invalid_argument("Invalid XML"));
    type=endcom;
    return;
    }
 for (pos=start;pos<stop;pos++){
    if (xmlstr[pos]=='<'){
       char next=(pos<(stop-1))?xmlstr[pos+1]:' ';
       if (next=='/'){
          type=startclose; return;}
       else if (next=='?'){
          type=startdecl; return;}
       else if (next=='!'){
          if (((pos+3)<stop)&&(xmlstr[pos+2]=='-')&&(xmlstr[pos+3]=='-')){
             type=startcom; return;}
          else
             throw(std::invalid_argument("Invalid XML"));}
       type=starttag; return;}
    else if (xmlstr[pos]=='>'){
       char prev=(pos>start)?xmlstr[pos-1]:' ';
       if (prev=='/'){
          type=endempty; return;}
       else if (prev=='?'){
          type=enddecl; return;}
       else if (prev=='-'){
          if (((pos-2)>=start)&&(xmlstr[pos-2]=='-')){
             type=endcom; return;}
          else
             throw(std::invalid_argument("Invalid XML"));}
       type=endtag; return;}
    }
 pos=stop;
}

    // Searches the string in "xmlstr" starting at character location "start"
    // and stopping at "in_stop" (not inclusive) for the first occurrence of
    // an XML event, skipping comments and accumulating any non-markup 
    // text that is passed, returning it in "textpassed".  It returns the type 
    // of event found in "type" and its location in "pos" (the location of the 
    // "<" or ">").  If no event is found, the type "none" is returned and pos 
    // will be set to "in_stop".  Throws an exception if a start declaration 
    // or end declaration is encountered, or if invalid comment events are found.

void XMLDoc::find_next_xml_tag_event(const string& xmlstr, size_t start, size_t in_stop,
                                     XMLDoc::XMLEvent& type, size_t& pos, 
                                     string& textpassed)
{
 textpassed.clear();
 size_t curr=start;
 bool incomment=false;
 bool notdone=true;
 while (notdone){
    find_next_xml_event(xmlstr,curr,in_stop,type,pos,incomment);
    if ((type==startdecl)||(type==enddecl)) throw(std::invalid_argument("Invalid XML"));
    if (type==endcom){
       if (!incomment) throw(std::invalid_argument("Invalid XML"));
       incomment=false;
       curr=pos+1;}
    else{
       if ((pos>curr)&&(type!=endtag)&&(type!=endempty)){
          textpassed.append(xmlstr,curr,pos-curr);}
       if (type==startcom){
          incomment=true; curr=pos+4;}
       else
          notdone=false;}}
 textpassed=trim(textpassed);
}

    // Searches the string in "xmlstr" starting at character location "start"
    // and stopping at "in_stop" (not inclusive) for the first occurrence of
    // an XML start tag event, skipping comments and testing for any non-whitespace 
    // non-markup text that is passed.  Throws an exception if non-whitespace text 
    // is encountered, if invalid comment events are found, or any event other
    // than comments or a start tag are found.

void XMLDoc::find_root_start_tag(const string& xmlstr, size_t start, size_t in_stop,
                                 size_t& pos)
{
 size_t curr=start;
 bool incomment=false;
 bool notdone=true;
 XMLEvent type;
 while (notdone){
    find_next_xml_event(xmlstr,curr,in_stop,type,pos,incomment);
    if ((type!=startcom)&&(type!=endcom)&&(type!=starttag))
       throw(std::invalid_argument("Invalid XML"));
    if (type==endcom){
       if (!incomment) throw(std::invalid_argument("Invalid XML"));
       incomment=false;
       curr=pos+1;}
    else{
       if (nonwhitespace(xmlstr,curr,pos)){
          throw(std::invalid_argument("Invalid XML"));}
       if (type==startcom){
          incomment=true; curr=pos+4;}
       else
          notdone=false;}}
}

   //  This routine finds the root tag and starts filling the
   //  "ancestors" stack if the root tag is not an empty tag.

XMLDoc::XMLNode* XMLDoc::parse_root_tag(const string& xmlstr, size_t& start, size_t stop, 
                                        XMLDoc::XMLEvent& type, size_t& pos, 
                                        std::stack<XMLDoc::XMLNode*>& ancestors,
                                        std::stack<XMLDoc::XMLNode*>& lastchildren)
{
 size_t pos1;
 find_root_start_tag(xmlstr,start,stop,pos1);
 find_next_xml_event(xmlstr,pos1+1,stop,type,pos,false);
 string roottag;
 if (type==endtag){
    roottag=get_tag_name(xmlstr,pos1+1,pos);}
 else if (type==endempty){
    roottag=get_tag_name(xmlstr,pos1+1,pos-1);}
 else
    throw(std::invalid_argument("Invalid XML"));
 pos++;
 XMLDoc::XMLNode* top=new XMLNode(roottag,0);
 if (type==endtag){
    ancestors.push(top);
    lastchildren.push(0);}
 return top;
}


  //  This routine parses the string "xmlstr" from the current location 
  //  until the next XML event is encountered, and builds up the XMLDoc nodes.
  //  On entering this routine, the following have already been set up
  //    ancestors  ->  stack of ancestor nodes of the current node
  //                     (must not be empty)
  //    curr       ->  current location in "xmlstr" string
  //    lastevent  ->  type of the last XML event (must be one of
  //                      starttag, endtag, startclose, endempty)
  //  The routine updates "curr", "lastevent", and "ancestors" while building
  //  up the XMLDoc nodes.

void XMLDoc::parse_to_next_event(const string& xmlstr, size_t& curr, size_t stop, 
                                 XMLEvent& lastevent, stack<XMLNode*>& ancestors,
                                 stack<XMLNode*>& lastchildren)
{
 XMLEvent nextevent;
 size_t pos;
 string textpassed;
 find_next_xml_tag_event(xmlstr,curr,stop,nextevent,pos,textpassed);

 if (ancestors.empty()){
    if ((nextevent!=none)||(!textpassed.empty()))
       throw(std::invalid_argument("Invalid XML"));
    lastevent=nextevent;
    return;}

 if (lastevent==starttag){
    if (nextevent==endempty){
        string tagname=get_tag_name(xmlstr,curr,pos-1);
        XMLNode *newtag=new XMLNode(tagname,ancestors.top(),false);
        if (lastchildren.top()==0) ancestors.top()->firstchild=newtag;
        else lastchildren.top()->nextsibling=newtag;
        lastchildren.top()=newtag;
        lastevent=endempty;
        curr=pos+1; return;
        }
    else if (nextevent==endtag){
        string tagname=get_tag_name(xmlstr,curr,pos);
        XMLNode *newtag=new XMLNode(tagname,ancestors.top(),false);
        if (lastchildren.top()==0) ancestors.top()->firstchild=newtag;
        else lastchildren.top()->nextsibling=newtag;
        lastchildren.top()=newtag;
        lastevent=endtag;
        curr=pos+1; 
        ancestors.push(newtag);
        lastchildren.push(0);
        return;}
    throw(std::invalid_argument("Invalid XML"));}
 if (lastevent==startclose){
    if (nextevent==endtag){
       string tagname=get_tag_name(xmlstr,curr,pos);
       if (tagname!=ancestors.top()->name) throw(std::invalid_argument("Invalid XML"));
       ancestors.pop();
       lastchildren.pop();
       lastevent=endtag;
       curr=pos+1;
       return;}
    throw(std::invalid_argument("Invalid XML"));}
 if (!textpassed.empty()){
    XMLNode *textnode=new XMLNode(textpassed,ancestors.top(),true);
    if (lastchildren.top()==0) ancestors.top()->firstchild=textnode;
    else lastchildren.top()->nextsibling=textnode;
    lastchildren.top()=textnode;}
 if (nextevent==starttag){
    curr=pos+1; lastevent=nextevent;
    return;}
 else if (nextevent==startclose){
    curr=pos+2; lastevent=nextevent;
    return;}
 else if (nextevent==none){
    curr=stop; lastevent=nextevent;
    return;}
 throw(std::invalid_argument("Invalid XML"));
}


    //  Routine used in checking the declaration string.  The
    //  attribute with key "attrkey" is searched for, and if found,
    //  its value is put into "attrvalue".  The location of the
    //  attribute element is returned in "start" (first character
    //  of the element) and "stop" (first character after the element).
    //  Return true if found, false otherwise.

bool XMLDoc::get_attribute(const string& instr, const string& attrkey,
                           string& attrvalue, size_t& start, size_t& stop)
{
 start=instr.find(string(" ")+attrkey);
 if (start==string::npos) return false;
 size_t curr=start+attrkey.length()+1;
 if ((instr.length()<(curr+2))||(instr[curr]!='=')
    ||(instr[curr+1]!='\"')) return false;
 curr=curr+2;
 stop=instr.find("\"",curr);
 if (stop==string::npos) return false;
 attrvalue.assign(instr,curr,stop-curr);
 stop++;
 return true;
}

    //  Extract the leading XML declaration tag from "instr",
    //  throwing an exception is invalid.  The tag is stored
    //  in the class member "declaration".  Returns the
    //  character location one past the end of the declaration tag.

size_t XMLDoc::get_declaration(const string& instr)
{
 XMLEvent type;
 size_t pos;
 find_next_xml_event(instr,0,instr.length(),type,pos);
 if (type!=startdecl) return 0;
 if (pos!=0) throw(std::invalid_argument("Invalid XML: declaration must start at first character"));
 if (instr.length()<8) throw(std::invalid_argument("Invalid XML declaration"));
 if ((instr[2]!='x')||(instr[3]!='m')||(instr[4]!='l')) throw(std::invalid_argument("Invalid XML"));
 find_next_xml_event(instr,4,instr.length(),type,pos);
 if (type!=enddecl) throw(std::invalid_argument("Invalid XML declaration"));
 declaration.assign(instr,0,pos+1);
      // check the validity of the XML declaration
 size_t verA,verB,encA,encB,stdA,stdB;
 string ver,enc,std;
 bool flag;
 flag=get_attribute(declaration,"version",ver,verA,verB);
 if ((!flag)||(ver!="1.0")) throw(std::invalid_argument("Invalid XML declaration"));

 flag=get_attribute(declaration,"encoding",enc,encA,encB);
 if (flag){
    if (!( (enc=="UTF-8")||(enc=="UTF-16")||(enc=="ISO-8859-1")) 
        ||(encA<verB)) throw(std::invalid_argument("Invalid XML declaration"));
    verA=encA; verB=encB;}
 flag=get_attribute(declaration,"standalone",std,stdA,stdB);
 if (flag){
    if (!( (std=="yes")||(std=="no"))||(stdA<verB))
       throw(std::invalid_argument("Invalid XML declaration"));
    verA=stdA; verB=stdB;}
 verA=declaration.length();
 if ((verA!=(verB+2))&&(verA!=(verB+3)))
    throw(std::invalid_argument("Invalid XML declaration"));

 return pos+1;
}

    //  "initial_indent" negative means "no indenting, no newlines"
    //  which is useful for file headers

ostream& XMLDoc::output(ostream& os, XMLNode* top, int initial_indent) const
{
 if (top==0) return os;
 string tab("   ");
 int indent=(initial_indent>0)?initial_indent:0;
 bool indflag=(initial_indent>=0)?true:false;
 XMLNode *curr=top;
 bool ascending=false;
 bool newline=indflag;
 bool done=(top->firstchild==0);
 if ((!declaration.empty())&&(curr==root)) os << declaration;
 while (!done){
    if (newline){
       os << endl;
       for (int k=0;k<indent;k++) os <<tab;}
          // will move down to child (curr cannot be empty nor text)
    if ((curr->firstchild!=0)&&(!ascending)){
       os << "<"<<curr->name<<">";
       if ((curr->firstchild->text)&&(curr->firstchild->firstchild==0)
            &&(curr->firstchild->nextsibling==0)){
          newline=false;}   // if curr+child = simple tag
       else{
          newline=indflag; indent++;}
       curr=curr->firstchild;}
           // will move to next sibling (either no child or ascending)
    else if (curr->nextsibling!=0){
       if (curr->text) os <<curr->name;
       else if (ascending) os << "</"<<curr->name<<">";
       else os << "<"<<curr->name<<"/>";
       ascending=false;
       newline=indflag;
       curr=curr->nextsibling;}
              // will move up to parent (no next sibling 
    else{     //        and (ascending or no first child)  )
       if (curr->text) os << curr->name;
       else if (ascending){
          os << "</"<<curr->name<<">"; newline=indflag;}
       else{
          os << "<"<<curr->name<<"/>"; newline=indflag;}
       if (newline) indent--;
       curr=curr->parent;
       ascending=true;
       if (curr==top) done=true;}}
 if (newline){
    os << endl;
    for (int k=0;k<indent;k++) os <<tab;}
 if (ascending) os << "</"<<top->name<<">"<<endl;
 else os << "<"<<top->name<<"/>"<<endl;
 return os;
}


void XMLDoc::add_reference(XMLHandler* ref)
{
 refset.insert(ref);
}

   // Removes the pointer "ref" from refset.  Returns
   // true if refset is empty afterwards; false otherwise.

bool XMLDoc::remove_reference(XMLHandler* ref)
{
 refset.erase(ref);
 return refset.empty();
}



   // *********************************************************



XMLHandler::XMLHandler(const string& roottagname, 
                       const string& rootcontent)
      : content(0), root(0), current(0), exceptions(true)
{
 do_set(roottagname,rootcontent);
}


XMLHandler::XMLHandler(const char* roottagname, const char* rootcontent)
      : content(0), root(0), current(0), exceptions(true)
{
 do_set(string(roottagname),string(rootcontent));
}


XMLHandler::XMLHandler(const char* roottagname)
      : content(0), root(0), current(0), exceptions(true)
{
 do_set(string(roottagname),string());
}


     // Copy constructor. If "mode" has value "subtree_pointer" or
     // "subtree_copy", the current location of "xmlin" becomes the root 
     // element of this handler (cannot access ancestor nodes).  Otherwise,
     // the root of this handler will be the same as "xmlin".  If "mode" has
     // value "subtree_copy" or "copy", then a completely new deep copy is 
     // made and the root element is either the current location or the root
     // of "xmlin"; otherwise, this is just a new pointer to the same data 
     // in memory and the current element is the same as in "xmlin".

XMLHandler::XMLHandler(const XMLHandler& xmlin, 
                       XMLHandler::copymode mode)
           : content(0), root(0), current(0), exceptions(xmlin.exceptions)
{
 bool subtree= (mode==subtree_pointer)||(mode==subtree_copy);
 bool makecopy= (mode==subtree_copy)||(mode==copy);
 do_set(xmlin,subtree,makecopy);
}

       // Copy constructor.  First does a "seek_unique" looking for
       // the tag with name specified in "newroottag".  Throws exception
       // if not found or multiple instances found.  This creates a new
       // pointer to a subtree of the same XML document.

XMLHandler::XMLHandler(const XMLHandler& xmlin, const std::string& newroottag,
                       bool seek_to_child_only)
           : content(0), root(0), current(0), exceptions(true)
{
 do_set(xmlin,false,false); 
 try{
    if (!seek_to_child_only) seek_unique(newroottag);
    else seek_unique_to_child(newroottag);}
 catch(const std::exception& msg){
    clear(); throw;}
 if (current!=0){
    root=current;
    exceptions=xmlin.exceptions;}
}

XMLHandler::XMLHandler(const XMLHandler& xmlin, const char* newroottag)
           : content(0), root(0), current(0), exceptions(true)
{
 do_set(xmlin,false,false); 
 try{
    seek_unique(newroottag);}
 catch(const std::exception& msg){
    clear(); throw;}
 if (current!=0){
    root=current;
    exceptions=xmlin.exceptions;}
}

        // no full copy, entire document retained

XMLHandler& XMLHandler::operator=(const XMLHandler& xmlin)  
{
 if (&xmlin==this) return *this;  // check for self-assignment
 clear();
 content=xmlin.content;
 current=xmlin.current;
 root=xmlin.root;
 if (content!=0){
    content->add_reference(this);}
 exceptions=xmlin.exceptions;
 return *this;
}


XMLHandler& XMLHandler::set(const string& roottagname, 
                            const string& rootcontent)
{
 clear();
 do_set(roottagname,rootcontent);
 return *this;
}


XMLHandler& XMLHandler::set(const char* roottagname, 
                            const char* rootcontent)
{
 clear();
 do_set(string(roottagname),string(rootcontent));
 return *this;
}


XMLHandler& XMLHandler::set(const char* roottagname)
{
 clear();
 do_set(string(roottagname),string());
 return *this;
}


XMLHandler& XMLHandler::set_root(const string& roottagname, 
                                 const string& rootcontent)
{
 clear();
 do_set(roottagname,rootcontent);
 return *this;
}


XMLHandler& XMLHandler::set_root(const char* roottagname, 
                                 const char* rootcontent)
{
 clear();
 do_set(string(roottagname),string(rootcontent));
 return *this;
}


XMLHandler& XMLHandler::set_root(const char* roottagname)
{
 clear();
 do_set(string(roottagname),string());
 return *this;
}


XMLHandler& XMLHandler::set_from_string(const std::string& xmlin)
{
 clear();
 do_set_from_string(xmlin);
 return *this;
}


XMLHandler& XMLHandler::set_from_file(const std::string& filename)
{
 clear();
 string xmlstr;
 try{
    ifstream fin(filename.c_str());
    if (!fin)
       throw(std::invalid_argument(string("Could not open file ")+filename));
    string line;
    while ((!fin.eof())&&(fin.good())){
       getline(fin,line);
       xmlstr+=line;}
    if (!fin.eof())
       throw(std::invalid_argument(string("Problem occurred while reading file ")+filename));
    fin.close();}
 catch(const std::exception& err){
    cout << "XML document file error: "<<err.what()<<endl;
    current=0;
    if (exceptions) throw;}
 do_set_from_string(xmlstr);
 return *this;
}


XMLHandler& XMLHandler::set(const XMLHandler& xmlin, 
                            XMLHandler::copymode mode)
{
 clear();
 bool subtree= (mode==subtree_pointer)||(mode==subtree_copy);
 bool makecopy= (mode==subtree_copy)||(mode==copy);
 do_set(xmlin,subtree,makecopy);
 return *this;
}


void XMLHandler::do_set(const string& roottagname, 
                        const string& rootcontent)
{
 try{
    content=new XMLDoc(roottagname,rootcontent);
    content->add_reference(this);
    root=content->root;
    current=root;}
 catch(const std::exception& xp){
    cout << "XML document error: "<<xp.what()<<endl;
    clear();
    if (exceptions) throw;}
}


void XMLHandler::do_set(const XMLHandler& xmlin, bool subtree, bool makecopy)
{
 XMLDoc::XMLNode *newroot=(subtree ? xmlin.current : xmlin.root);
 if (newroot->text)
    throw(std::runtime_error("XMLHandler copy failed: current node is text valued, cannot be root tag"));
 if ((makecopy)&&(newroot!=0)){
    try{
       content=new XMLDoc(newroot,xmlin.content->declaration);
       root=content->root;
       current=root;
       content->add_reference(this);}
    catch(const std::exception& xp){
       cout << "XML document error: "<<xp.what()<<endl;
       clear();
       if (exceptions) throw;}
    }
 else if (newroot!=0){
    content=xmlin.content;
    root=newroot;
    //  cannot set root->nextsibling=0; since affects original XMLdoc
    current=xmlin.current;
    if (content!=0){
       content->add_reference(this);}}
 else{
    if (exceptions) throw(std::invalid_argument("XML set failed"));}
}


void XMLHandler::do_set_from_string(const std::string& xmlin)
{
 try{
    content=new XMLDoc(xmlin);
    content->add_reference(this);
    root=content->root;
    current=root;}
 catch(const std::exception& xp){
    cout << "XML document error: "<<xp.what()<<endl;
    clear();
    if (exceptions) throw;}
}


XMLHandler::~XMLHandler()
{
 clear();
}


XMLHandler& XMLHandler::clear()
{
 if (!content) return *this;
 if (content->remove_reference(this)){
    delete content;}
 content=0; root=0; current=0;
 return *this;
}


string XMLHandler::output(int indent) const
{
 if (!content) return string();
 ostringstream os;
 content->output(os,root,indent);
 return os.str();
}

string XMLHandler::str() const
{
 if (!content) return string();
 ostringstream os;
 content->output(os,root,-1);   // no indenting
 return os.str();
}

string XMLHandler::output_current(int indent) const
{
 if (!content) return string();
 ostringstream os;
 content->output(os,current,indent);
 return os.str();
}



      //  seek members (if seek fails, state of handler
      //  is set to "bad", and if exceptions are turned
      //  on, a string exception is thrown)

void XMLHandler::seek_root()
{
 current=root;
 if ((current==0)&&(exceptions))
    throw(std::invalid_argument("Seek root failed: no document loaded"));
}

void XMLHandler::seek_next_sibling()
{
 if (current==root) current=0;
 else if (current) current=current->nextsibling;
 if ((current==0)&&(exceptions))
    throw(std::invalid_argument("Seek next sibling failed"));
}

void XMLHandler::seek_first_child()
{
 if (current) current=current->firstchild;
 if ((current==0)&&(exceptions))
    throw(std::invalid_argument("Seek first child failed"));
}

void XMLHandler::seek_child(const string& tagname)
{
 int ncount=0;
 if (current){
    XMLHandler tagfinder(*this,XMLHandler::subtree_pointer);
    tagfinder.seek_first_child();
    tagfinder.set_exceptions_off();
    while ((tagfinder.good())&&(ncount<2)){
       if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
          current=tagfinder.current;
          ncount++;}
       tagfinder.seek_next_sibling();}}
 if (ncount!=1){
    current=0;
    if (exceptions) 
       throw(std::invalid_argument(string("seek child tag ")+tagname+string(" failed: none or plural")));}
}

void XMLHandler::seek_parent()
{
 if (current==root) current=0;
 else if (current) current=current->parent;
 if ((current==0)&&(exceptions))
    throw(std::invalid_argument("Seek parent failed"));
}

void XMLHandler::seek_next_node()
{
 if (current){
    if (current->firstchild!=0){
       current=current->firstchild;}
    else if ((current->nextsibling!=0)&&(current!=root)){
       current=current->nextsibling;}
    else{
       current=current->parent;
       while (seek_next_sib_or_parent());}}
 if ((current==0)&&(exceptions))
    throw(std::invalid_argument("Seek next node failed"));
}

bool XMLHandler::seek_next_sib_or_parent()
{
 if (!current) return false;
 else if (current==root){
    current=0; return false;}
 else if (current->nextsibling!=0){
    current=current->nextsibling;
    return false;}
 else{
    current=current->parent;
    if (current==root->parent) current=0;
    return true;}
}


void XMLHandler::seek_next(const std::string& tagname)
{
 seek_next_node();
 while ((current!=0)&&((current->text)||(current->name!=tagname)))
    seek_next_node();
 if ((current==0)&&(exceptions))
    throw(std::invalid_argument(string("Seek next node by tag name ")
         +tagname+string(" failed")));
}


void XMLHandler::seek_unique(const string& tagname)
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 while ((tagfinder.good())&&(ncount<2)){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       current=tagfinder.current;
       ncount++;}
    tagfinder.seek_next_node();}
 if (ncount!=1){
    current=0;
    if (exceptions) 
       throw(std::invalid_argument(string("seek unique of tag ")+tagname
           +string(" failed: none or plural")));}
}

void XMLHandler::seek_unique_to_child(const string& tagname)
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 bool ischild=false;
 while ((tagfinder.good())&&(ncount<2)){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       current=tagfinder.current;
       ncount++;}
    if (ischild) tagfinder.seek_next_sibling();
    else{ tagfinder.seek_first_child(); ischild=true;}}
 if (ncount!=1){
    current=0;
    if (exceptions) 
       throw(std::invalid_argument(string("seek unique of tag ")+tagname
             +string(" failed: none or plural")));}
}

void XMLHandler::seek_unique_child(const string& tagname)
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 tagfinder.seek_first_child(); 
 while ((tagfinder.good())&&(ncount<2)){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       current=tagfinder.current;
       ncount++;}
    tagfinder.seek_next_sibling();}
 if (ncount!=1){
    current=0;
    if (exceptions) 
       throw(std::invalid_argument(string("seek unique of tag ")+tagname
            +string(" failed: none or plural")));}
}

list<XMLHandler> XMLHandler::find(const string& tagname) const
{
 list<XMLHandler> found;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 while (tagfinder.good()){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       XMLHandler xmlt(tagfinder);
       found.push_back(xmlt);}
    tagfinder.seek_next_node();}
 return found;
}

list<XMLHandler> XMLHandler::find(const list<string>& tagnames) const
{
 list<XMLHandler> found;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 while (tagfinder.good()){
    if (!(tagfinder.current->text)){
       for (list<string>::const_iterator it=tagnames.begin();it!=tagnames.end();it++){
          if (tagfinder.current->name==*it){
             XMLHandler xmlt(tagfinder);
             found.push_back(xmlt); break;}}}
    tagfinder.seek_next_node();}
 return found;
}


int XMLHandler::count(const string& tagname) const
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 while (tagfinder.good()){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       ncount++;}
    tagfinder.seek_next_node();}
 return ncount;
}


list<XMLHandler> XMLHandler::find_among_children(const string& tagname) const
{
 list<XMLHandler> found;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.set_exceptions_off();
 tagfinder.seek_first_child();
 while (tagfinder.good()){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       XMLHandler xmlt(tagfinder);
       found.push_back(xmlt);}
    tagfinder.seek_next_sibling();}
 return found;
}


list<XMLHandler> XMLHandler::find_among_children(const list<string>& tagnames) const
{
 list<XMLHandler> found;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.set_exceptions_off();
 tagfinder.seek_first_child();
 while (tagfinder.good()){
    if (!(tagfinder.current->text)){
       for (list<string>::const_iterator it=tagnames.begin();it!=tagnames.end();it++){
          if (tagfinder.current->name==*it){
             XMLHandler xmlt(tagfinder);
             found.push_back(xmlt); break;}}}
    tagfinder.seek_next_sibling();}
 return found;
}


int XMLHandler::count_among_children(const string& tagname) const
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.set_exceptions_off();
 tagfinder.seek_first_child();
 while (tagfinder.good()){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       ncount++;}
    tagfinder.seek_next_sibling();}
 return ncount;
}

int XMLHandler::count_to_among_children(const string& tagname) const
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
// tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 bool ischild=false;
 while (tagfinder.good()){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       ncount++;}
    if (ischild) tagfinder.seek_next_sibling();
    else{ tagfinder.seek_first_child(); ischild=true;}}
 return ncount;
}

bool XMLHandler::query_unique_to_among_children(const string& tagname) const
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
// tagfinder.seek_root();
 tagfinder.set_exceptions_off();
 bool ischild=false;
 while (tagfinder.good()){
    if ((!(tagfinder.current->text))&&(tagfinder.current->name==tagname)){
       ncount++;}
    if (ncount>1)
       throw(std::invalid_argument(std::string("Multiple occurrences of ")
            +tagname+std::string(" when one or none required")));
    if (ischild) tagfinder.seek_next_sibling();
    else{ tagfinder.seek_first_child(); ischild=true;}}
 return ncount;
}

int XMLHandler::count_children() const
{
 int ncount=0;
 XMLHandler tagfinder(*this,XMLHandler::pointer);
 tagfinder.seek_first_child();
 tagfinder.set_exceptions_off();
 while (tagfinder.good()){
    ncount++;
    tagfinder.seek_next_sibling();}
 return ncount;
}


std::string XMLHandler::get_node_name() const
{
 return (current) ? ((current->text) ? string("#text"): current->name)
          : empty_string();
}


std::string XMLHandler::get_node_value() const
{
 return is_text_valued() ? current->name : empty_string();
}


bool XMLHandler::is_text_valued() const
{
 return (current) ? current->text : false;
}


string XMLHandler::get_tag_name() const     // returns empty string if invalid
{
 return ((current)&&(!(current->text))) ? current->name : empty_string();
}

 
string XMLHandler::get_text_content() const // returns empty string if invalid
{
 return is_simple_element() ? current->firstchild->name : empty_string();
}


bool XMLHandler::is_simple_element() const
{
 return (current)&&(!(current->text))&&(current->firstchild)
        &&(current->firstchild->text)&&(current->firstchild->firstchild==0)
        &&(current->firstchild->nextsibling==0);
}


bool XMLHandler::is_empty_tag() const
{
 return (current)&&(!(current->text))&&(current->firstchild==0);
}

string XMLHandler::empty_string() const
{
 if (exceptions) throw(std::invalid_argument("Error"));
 return string();
}


             // Insertion. "put_sibling" inserts as the
             // next sibling of the current node, then the
             // current node pointer is updated to the top
             // node of the inserted node(s).  Fails if the
             // current node is the root node or current
             // state is not good.

void XMLHandler::put_sibling(const XMLHandler& xmlin)
{
 try{
    if ((fail())||(xmlin.fail())||(current==root))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_copy(xmlin.root);
    content->connect_as_nextsibling(current,branch);
    seek_next_sibling();}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_sibling failed"));}
}


void XMLHandler::put_sibling(const string& tagname)
{
 try{
    if ((fail())||(current==root))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_simple_tag(tagname,string(),0);
    content->connect_as_nextsibling(current,branch);
    seek_next_sibling();}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_sibling failed"));}
}


void XMLHandler::put_sibling(const string& tagname, 
                             const string& text_content)
{
 try{
    if ((fail())||(current==root))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_simple_tag(tagname,text_content,0);
    content->connect_as_nextsibling(current,branch);
    seek_next_sibling();}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_sibling failed"));}
}


void XMLHandler::put_sibling_text_node(const string& text)
{
 try{
    if ((fail())||(current==root))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_text_node(text,0);
    if (branch==0) throw(std::invalid_argument("Error"));
    content->connect_as_nextsibling(current,branch);
    seek_next_sibling();}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_sibling_text_node failed"));}
}


          // Insertion. "put_child" inserts as the
          // last child of the current node, but the
          // current node pointer is NOT changed.  Fails
          // if the current node is a text node.

void XMLHandler::put_child(const XMLHandler& xmlin)
{
 try{
    if ((fail())||(xmlin.fail())||(current->text))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_copy(xmlin.root);
    content->connect_as_lastchild(current,branch);}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_child failed"));}
}


void XMLHandler::put_child(const string& tagname)
{
 try{
    if ((fail())||(current->text))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_simple_tag(tagname,string(),0);
    content->connect_as_lastchild(current,branch);}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_child failed"));}
}


void XMLHandler::put_child(const string& tagname, 
                           const string& text_content)
{
 try{
    if ((fail())||(current->text))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_simple_tag(tagname,text_content,0);
    content->connect_as_lastchild(current,branch);}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_child failed"));}
}


void XMLHandler::put_child_text_node(const string& text)
{
 try{
    if ((fail())||(current->text))
       throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *branch=content->create_text_node(text,0);
    if (branch==0) throw(std::invalid_argument("Error"));
    content->connect_as_lastchild(current,branch);}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("put_child_text_node failed"));}
}


          // Inserts or changes the textual content
          // of the current node (if a text node).

void XMLHandler::set_text_content(const string& in_text_content)
{
 if (good()){
    string text_content(content->trim(in_text_content));
    if ((current->text)&&(!text_content.empty())){
       current->name=text_content; return;}}
 clear(); throw(std::invalid_argument("set_text_content failed"));
}

          // Renames the current tag.
          // Exception throw if failure occurs.  All XMLHanders
          // pointing to the same content are changed.

void XMLHandler::rename_tag(const string& newtagname)
{
 string tag(content->get_tag_name(newtagname,0,newtagname.length()));
 try{
    if ((fail())||(current->text))
       throw(std::invalid_argument("Error"));
    current->name=tag;}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("rename tag failed"));}
}


          // Erases current element (if not a text node).
          // Current pointer will point to parent node afterwards.

void XMLHandler::erase_current_element()
{
 try{
    if ((fail())||(current->text)) throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *p=current->parent;
    content->disconnect(current);
    content->destroy(current);
    current=p;}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("erase_current_element failed"));}
}

          // Erases current node if a text node.
          // Current pointer will point to parent node afterwards.

void XMLHandler::erase_text_node()
{
 try{
    if ((fail())||(!current->text)) throw(std::invalid_argument("Error"));
    XMLDoc::XMLNode *p=current->parent;
    content->disconnect(current);
    content->destroy(current);
    current=p;}
 catch(const std::exception& xp){
    clear(); throw(std::invalid_argument("erase_text_node failed"));}
}


bool XMLHandler::IsEqualContent(XMLHandler& xmlh_cmp, float float_rel_tol)
{
 XMLContentComparer temp(float_rel_tol);
 return temp.IsEqual(*this,xmlh_cmp);
}


// ***************************************************************


  // Beginning with character at "start", moves the "character"
  // pointer to the first non-ignorable character (or NULL if
  // no such characters found).  "nchar" will contain the 
  // number of non-ignorable characters in the next token. 

void XMLContentComparer::get_next_token(const char*& start, int& nchar)
{
 char delimit[] = " \t\n\r";
 start += strspn(start, delimit); 
 if (*start=='\0'){ start=0; nchar=0; return;}
 const char *stop = strpbrk(start, delimit);
 if (stop==0) nchar=strlen(start);
 else{
    nchar=0; const char *tmp=start;
    while (tmp!=stop){ tmp++; nchar++;}}
}

   //  Generalization of strcmp but not using null characters
   //  as string ending; uses numbers of characters in na,nb
   //  Return 0 if na-character string in "a" equals nb-character
   //  string in "b".  Negative return if a-string precedes b-string,
   //  and >0 return if a-string comes after b-string.

int XMLContentComparer::token_strcmp(const char *a, int na, 
                                     const char *b, int nb)
{
 if (na==nb) return strncmp(a,b,na);
 else if (na>nb){
    int k=strncmp(a,b,nb);
    if (k!=0) return k;
    else return 1;}
 else{
    int k=strncmp(a,b,na);
    if (k!=0) return k;
    else return -1;}
}

   // Check if the "nchar"-character token starting at "a"
   // is a valid integer; if so, return integer value in "ivalue".
   // We assume the token has no blank characters.

bool XMLContentComparer::is_integer_token(const char *const a, int nchar, 
                                          int& ivalue)
{
 if (nchar==0) return false;
 const char* ap=a;
 int k=0;
 if ((*ap=='+')||(*ap=='-')){ap++; k++;}
 if (k==nchar) return false;
 while (k<nchar){
    if (!isdigit(*ap)) return false;
    ap++; k++;}
 stringstream oss; oss << "%" << nchar << "d";
 return sscanf(a,oss.str().c_str(),&ivalue);
}

   // Check if the "nchar"-character token starting at "a"
   // is a valid floating-point number; if so, return value in "fvalue".
   // We assume the token has no blank characters.

bool XMLContentComparer::is_float_token(const char *const a, int nchar,
                                        float& fvalue)
{
 if (nchar==0) return false;
 const char* ap=a;
 int k=0;
 if ((*ap=='+')||(*ap=='-')){ap++; k++;}
 if (k==nchar) return false;
 bool dotflag=false, eflag=false;
 while (k<(nchar-1)){
    if (*ap=='.'){
       if (dotflag) return false;
       dotflag=true;}
    else if ((*ap=='e')||(*ap=='E')){
       if (eflag) return false;
       dotflag=true;
       eflag=true;  ap++; k++;
       if ((*ap!='+')&&(*ap!='-')){ ap--; k--;}
       else if (k==nchar) return false;}
    else if (!isdigit(*ap)) return false;
    ap++; k++;
    }
 if (!((isdigit(*ap))||((*ap=='.')&&(!dotflag)))) return false; 
     // last character must be digit or . (if not dots yet)
 stringstream oss; oss << "%" << nchar << "g";
 return sscanf(a,oss.str().c_str(),&fvalue);
}

    // This routine compares the content in the token strings "a" and "b"
    // having "na" and "nb" number of characters.  First the token are
    // compared as integers, then as floats (to within "float_rel_tol"),
    // then finally as strings.  A zero value is returned if they are
    // equal, a number > 0 if a>b and a number <0 is a<b.

int XMLContentComparer::token_content_cmp(const char *a, int na, 
                                          const char *b, int nb)
{
 int ia,ib;
 float fa,fb;
 if (is_integer_token(a,na,ia)&&is_integer_token(b,nb,ib)){
    if (ia>ib) return 1;
    else if (ia==ib){ return 0;}
    else return -1;}
 else if (is_float_token(a,na,fa) && is_float_token(b,nb,fb)){
    float favg=0.5*abs(fa+fb);
    if ((fa-fb)>(favg*m_float_rel_tol)) return 1;
    else if ((fa-fb)<-(favg*m_float_rel_tol)) return -1;
    else return 0;} 
 return token_strcmp(a,na,b,nb);
}


int XMLContentComparer::content_cmp(const char *a, const char *b)
{
 int result;
 const char *ap=a; int na;
 const char *bp=b; int nb;
 get_next_token(ap,na);
 get_next_token(bp,nb);
 while ((ap!=0)&&(bp!=0)){
    result=token_content_cmp(ap,na,bp,nb);
    if (result!=0) return result;
    ap+=na; get_next_token(ap,na);
    bp+=nb; get_next_token(bp,nb);
    }
 if ((ap==0)&&(bp==0)) return 0;
 if (ap!=0) return 1;
 else return -1;
}


int XMLContentComparer::node_cmp(xmlNodePtr a_node,
                                 xmlNodePtr b_node)
{
    // compare node types first
 if (a_node->text < b_node->text){ return -1;}
 else if (a_node->text > b_node->text){ return 1;}

    // if text node, compare by textual content
 if (a_node->text){
    if (!(b_node->text)) exit(1); // something went wrong
    return content_cmp(a_node->name.c_str(), b_node->name.c_str());}

    // compare node names next
 int result=strcmp(a_node->name.c_str(), b_node->name.c_str());
 if (result!=0) return result;

    // now compare by children
 return sibling_node_cmp(a_node->firstchild,b_node->firstchild);
}


int XMLContentComparer::sibling_node_cmp(xmlNodePtr a_node,
                                         xmlNodePtr b_node)
{
 xmlNodePtr cur_node;

 list<xmlNodePtr> alist;
 for (cur_node=a_node; cur_node; cur_node=cur_node->nextsibling)
    alist.push_back(cur_node);
 alist.sort(*this);

 list<xmlNodePtr> blist;
 for (cur_node=b_node; cur_node; cur_node=cur_node->nextsibling)
    blist.push_back(cur_node);
 blist.sort(*this);

 int na=alist.size();
 int nb=blist.size();
 int n = (na>=nb)?nb:na;

 list<xmlNodePtr>::const_iterator at=alist.begin();
 list<xmlNodePtr>::const_iterator bt=blist.begin();

 for (int k=0;k<n;k++){
    int result=node_cmp(*at,*bt);
    if (result!=0) return result;
    at++; bt++;}
 if ((at==alist.end())&&(bt==blist.end())) return 0;
 if (at==alist.end()) return -1;
 else return 1; 
}


bool XMLContentComparer::operator()(xmlNodePtr a_node,
                                    xmlNodePtr b_node)
{
 return (node_cmp(a_node,b_node)<0);
}

bool XMLContentComparer::IsEqual(XMLHandler& xmlh1, XMLHandler& xmlh2)
{
    //Get the root element nodes
 if ((xmlh1.fail())||(xmlh2.fail())) 
    throw(std::invalid_argument("Comparison of invalid XMLHandler"));

 int result=sibling_node_cmp(xmlh1.root,xmlh2.root);
 return (result==0);
}

// *************************************************************

  // removes tabs, newline, linefeed characters, then trims
  // leading and trailing blanks.

string tidyString(const string& str)   
{
 string tmp;
 for (size_t i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 size_t start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 size_t len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}


  // first tidies the string in "str", then removes any
  // leading path/subdirectory information in the file name

string tidyFileName(const string& str)   
{
 string tmp(tidyString(str));
 size_t pos=tmp.rfind('/');
 if (pos==string::npos) return tmp;
 else return tmp.substr(pos+1,tmp.length()-pos);
}

  // Trims leading and trailing white space, then checks
  // to make sure each character is alphanumeric, underscore, 
  // or period or backslash. If name is invalid, an empty string is returned.

string tidyName(const string& str)
{
 string tmp;
 size_t len=str.length();
 if (len==0) return tmp;
 size_t start=0;
 while ((isspace(str[start]))&&(start<len)) start++;
 if (start==len) return tmp;
 size_t stop=len-1;
 while ((isspace(str[stop]))&&(stop>start)) stop--;
 for (size_t i=start;i<=stop;i++){
    char c=str[i];
    if (isalnum(c)||(c=='_')||(c=='.')||(c=='/'))
       tmp.push_back(c);
    else return string("");}
 return tmp;
}


// *********************************************************

    // This returns the number of times that the tag "tagname"
    // is found in the XML document "xmlh".
    // A **tag name** should be input.

int xml_tag_count(XMLHandler& xmlh, const string& tagname)
{
 return xmlh.count(tagname);
}

int xml_child_tag_count(XMLHandler& xmlh, const string& tagname)
{
 return xmlh.count_among_children(tagname);
}


    // Check for existence of a unique tag "tagname" and if not
    // present, throws exception.
  
void xml_tag_assert(XMLHandler& xmlh, const std::string& tagname)
{
 if (xml_tag_count(xmlh,tagname)!=1){
    cout << "Bad XML input:"<<endl;
    cout << "Expected one <"<<tagname<<"> tag"<<endl;
    xmlreadfail(xmlh,tagname);}
 xmlh.seek_unique(tagname);
}


    // Check for existence of a unique tag "tagname1" or "tagname2" 
    // and if neither present, throw exception.
  
void xml_either_tag_assert(XMLHandler& xmlh, const std::string& tagname1,
                           const std::string& tagname2)
{
 if (xml_tag_count(xmlh,tagname1)==1){
    xmlh.seek_unique(tagname1); return;}
 else if (xml_tag_count(xmlh,tagname2)==1){
    xmlh.seek_unique(tagname2); return;}
 cout << "Bad XML input:"<<endl;
 cout << "Expected one <"<<tagname1<<"> or <"<<tagname2<<"> tag"<<endl;
 xmlreadfail(xmlh,tagname1);
}


void xml_tag_assert(XMLHandler& xmlh, const std::string& tagname,
                    const std::string& infoname)
{
 if (xml_tag_count(xmlh,tagname)!=1){
    cout << "Bad XML input to "<<infoname<<endl;
    cout << "Expected one <"<<tagname<<"> tag"<<endl;
    xmlreadfail(xmlh,infoname);}
 xmlh.seek_unique(tagname);
}


void xml_child_assert(XMLHandler& xmlh, const std::string& tagname,
                      const std::string& infoname)
{
 if (xml_child_tag_count(xmlh,tagname)!=1){
    cout << "Bad XML input:"<<endl;
    cout << "Expected one <"<<tagname<<"> tag among children"<<endl;
    xmlreadfail(xmlh,tagname);}
 xmlh.seek_unique_child(tagname);
}


void xml_child_assert(XMLHandler& xmlh, const std::string& tagname)
{
 if (xml_child_tag_count(xmlh,tagname)!=1){
    cout << "Bad XML input:"<<endl;
    cout << "Expected one <"<<tagname<<"> tag among children"<<endl;
    xmlreadfail(xmlh,tagname);}
 xmlh.seek_unique_child(tagname);
}


void xml_root_assert(XMLHandler& xmlh, const std::string& tagname,
                     const std::string& infoname)
{
 try{
    xmlh.seek_root();
    if (xmlh.get_node_name()!=tagname) throw(std::invalid_argument("assert"));}
 catch(const std::exception& xp){
    cout << "Bad XML input to "<<infoname<<endl;
    cout << "Expected root tag to be <"<<tagname<<">"<<endl<<endl;
    xmlreadfail(xmlh,infoname);}
}


void xml_root_assert(XMLHandler& xmlh, const std::string& tagname)
{
 try{
    xmlh.seek_root();
    if (xmlh.get_node_name()!=tagname) throw(std::invalid_argument("assert"));}
 catch(const std::exception& xp){
    cout << "Bad XML input:"<<endl;
    cout << "Expected root tag to be <"<<tagname<<">"<<endl<<endl;
    xmlreadfail(xmlh,tagname);}
}


void xml_cerr(XMLHandler& xmlh, const std::string& msg)
{
 cout << msg << endl;
}

   // outputs the current context for a failed xml read
   // (typically in an Info constructor) then throws exception.

void xmlreadfail(XMLHandler& xmlh, const std::string& infoname)
{
 cout << "*****ERROR**** "<<infoname<<" constructor failed"
      << " on XML content:"<<xmlh.output_current()<<endl;
 throw(std::invalid_argument("XML Read failure"));
}


// *************************************************************


bool xmlContentIsEqual(const string& doc1, const string& doc2, 
                       float float_rel_tol)
{
 XMLHandler xmlh1,xmlh2;
 xmlh1.set_from_string(doc1);
 xmlh2.set_from_string(doc2);
 return xmlh1.IsEqualContent(xmlh2,float_rel_tol);
}


bool xmlContentIsEqual(XMLHandler& xmlh1, XMLHandler& xmlh2, 
                       float float_rel_tol)
{
 return xmlh1.IsEqualContent(xmlh2,float_rel_tol);
}


bool headerMatch(const string& doc1, const string& doc2, 
                 float float_rel_tol)
{
 if (doc1==doc2) return true;
 return xmlContentIsEqual(doc1,doc2,float_rel_tol);
}

// ************************************************************

    //  tests if file having name "file_name" exists

bool fileExists(const std::string& file_name)
{
 return (access(file_name.c_str(),F_OK) == 0) ? true : false;
}


// **********************************************************
