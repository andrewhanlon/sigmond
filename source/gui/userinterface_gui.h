#ifndef USER_INTERFACE_GUI_H
#define USER_INTERFACE_GUI_H

#include <gtkmm.h>
#include <list>
//#include "InfoList.h"
//#include "LogHandler.h"
//#include "TaskList.h"
//using namespace std;


   // *********************************************************************
   //
   //  The classes declared in this file define the graphical user interface
   //  (GUI) for analyzing Monte Carlo data.  Such analysis needs input from
   //  the user: which tasks to perform and which information to use.  This 
   //  GUI is built using the gtkmm 3.0 toolkit.
   //
   //  The important classes defined in this file are as follows:
   //
   //    (1) UserInterface
   //          - a singleton which controls all windows
   //          - creates, destroys, and keeps track of all UINode objects
   //    (2) UINode
   //          - an abstract base class
   //    (3) TaskMenu
   //          - derived from UINode
   //          - offers choices of tasks that can be done
   //          - an object of class "TaskList" specifies the menus to build
   //          - output from the tasks is shown via a "LogDisplay"
   //          - contains a "HelpHandler" to display helpful info to user
   //    (4) InfoMenu
   //          - derived from UINode
   //          - prompts the user for needed information to carry out tasks
   //          - an object of class "InfoList" specifies the needed info
   //          - contains a "HelpHandler" to display helpful info to user
   //    (5) HelpHandler
   //          - controls the display of a help window for a given UINode
   //          - ensures only one window is shown for a given UINode
   //          - kills the help window if the parent UINode is closed
   //    (6) HelpDisplay
   //          - derived from UINode
   //          - shows help for the current menu
   //          - accessed only through a "HelpHandler" object
   //    (7) LogDisplay
   //          - derived from Gtk::Frame
   //          - a log file chapter-by-chapter browser used by "TaskMenu"
   //    (8) ProgressInterrupter
   //          - derived from UINode
   //          - allows the user to interrupt the computation
   //
   //  The analysis software will also create and handle pop-up dialogs as 
   //  needed.  These are implemented as Gtk::Dialog's and not as Gtk::Window's.
   //  The above classes cannot be accessed by the programmer directly.
   //  The code is designed for access only through "UserInterface" (see below).
   //
   //  The base class "InfoMenuItem" and derived classes "StringMenuItem", etc.
   //  are also defined to facilitate the creation of "InfoMenu" objects.
   //
   // *****************************************************************

   // forward declarations

class UINode; 
//class TaskMenu;
//class InfoMenu;
//class HelpHandler;
//class HelpDisplay;
//class ProgressInterrupter;

//class InfoMenuItem;
//class StringMenuItem;
//class OutputFileMenuItem;
//class InputFileMenuItem;
//class TextSelectionSetMenuItem;
//class RealNumberMenuItem;


//typedef sigc::slot<Real> ComputeEngine;

   // *****************************************************************
   //
   //  "UserInterface" is a singleton that handles all UINodes (Monte 
   //  Carlo Data Analysis Windows).  It creates and destroys them, keeps 
   //  track of them, and runs and quits the main loop.  Before creating
   //  and running any windows, it is necessary to first call 
   //
   //       UserInterface::Initiate(argc,argv,helpdir)
   //
   //  with the standard argc,argv arguments to main().  This requirement
   //  comes from the necessity to create a Gtk::Main object as dictated
   //  by gtkmm.  Failure to call this routine first will result in a
   //  fatal run-time error.  "helpdir" should be the directory which
   //  contains all files used to display help messages.  Generally, this
   //  is passed in by the script that runs the program so the end user
   //  never needs to worry about this.
   //
   //  A "TaskMenu" and "InfoMenu" are created and run using the static
   //  members
   //
   //       UserInterface::CreateRun_TaskMenu(....)
   //       UserInterface::CreateRun_InfoMenu(....)
   //
   //  See comments before "TaskMenu" and "InfoMenu" below for more information.
   //  A progress bar, allowing computation interruption, is created and run
   //  using the static member
   //
   //       UserInterface::Run_DisplayProgress(...)
   //
   //  The three routines above do everything that is needed in gtkmm to create 
   //  and show the window, as well as starting the main event loop.  No call
   //  to Gtk::Main::run() is needed.  Windows are destroyed using the
   //  static member function
   //
   //       UserInterface::DestroyWindow(UINode *w)
   //
   //  which destroys the window pointed to by "w". UserInterface maintains
   //  a list of pointers to all windows, and when this list becomes empty,
   //  UserInterface automatically quits the main event loop.  When a given
   //  window is destroyed, all of its inner windows are also destroyed (see
   //  later).
   //
   //  In "gtkmm", windows are generally C++ objects, with their constructors
   //  and destructors creating and destroying the windows.  Although one
   //  usually has some control over constructing objects, destructors are
   //  only called when the object goes out of scope.  The class "UserInterface"
   //  is used to circumvent the problems of relying on scope/constructors
   //  to control windows.  "UserInterface" uses the free store (heap) to
   //  create and destroy windows at the exact moments desired by the user.
   //  Hence, all windows essentially have global scope.  Also, when a window
   //  is destroyed, "delete" is used so all of its resources are freed up.
   //  Merely using gtkmm's "hide" to remove a window as an object eliminates the
   //  display of the window, but its resources are still present, taking up
   //  memory until the destructor is automatically called when the window
   //  object goes out of scope.  
   //
   //  Since an object cannot "delete" itself, the use of "new" and "delete"
   //  to control window management necessitates an external object, such as
   //  "UserInterface", to manage the creation and destruction of window
   //  objects.
   //
   //  Other helpful routines defined in "UserInterface" are
   //
   //      void DisplayMessage(UINode& GW, const string& message);
   //      bool Confirm(UINode& GW, const string& prompt);
   //      bool GetInputFileName(UINode& GW, string& filename);
   //      bool GetInputFileName(UINode& GW, string& filename, 
   //                            const string& prompt);
   //      bool GetOutputFileName(UINode& GW, string& filename);
   //      bool GetOutputFileName(UINode& GW, string& filename, 
   //                                   const string& prompt);
   //
   //  A window with help information is created (and run) using the private
   //  static member:
   //
   //      UserInterface::CreateShow_HelpDisplay(HelpHandler *handler) 
   //
   //  Generally, this is called by an object of class "HelpHandler" when
   //  the user requests that help for a given menu is displayed.
   //
   // ********************************************************************
  
class UserInterface 
{

   Glib::RefPtr<Gtk::Application> m_app;
   std::list<UINode*> m_windows;            // crucial list of window pointers

 public:

   UserInterface();

   ~UserInterface();                       // destructor
   void clear();

/*
   static void Initiate(int argc, char *argv[], const string& helpdir);

   static void CreateRun_InfoMenu(const InfoList& inflist,
                                  bool& proceed,
                                  UINode *outerwindow=0,
                                  const string& infofile=""); 

   static void CreateRun_TaskMenu(const TaskList& tasklist,
                                  LogHandler& outlog); 

   static void Run_DisplayProgress(UINode *outerwindow,
                                   const ComputeEngine& calculator,
                                   bool& completed,
                                   const string& title,
                                   bool PercentageMode=true);

   static void DeleteWindow(UINode *w);

           //  Miscellaneous helper routines

   static void DisplayMessage(UINode& GW, const string& message);

   static bool Confirm(UINode& GW, const string& prompt);

   static bool GetInputFileName(UINode& GW, string& filename);

   static bool GetInputFileName(UINode& GW, string& filename, 
                                const string& prompt);

   static bool GetOutputFileName(UINode& GW, string& filename);

   static bool GetOutputFileName(UINode& GW, string& filename, 
                                 const string& prompt);

*/
 private:

   UserInterface(const UserInterface&);            // copy prevented
   UserInterface& operator=(const UserInterface&); // copy assignment prevent

 public:


/*
 private:

   static list<UINode*> m_windows;      // crucial list of window pointers
   static bool m_initiated;

   static HelpDisplay* CreateShow_HelpDisplay(HelpHandler *handler);
*/
   friend class UINode;
//   friend class HelpHandler;
  
};

   // *************************************************************
   // 
   //  Monte Carlo Data Analysis window base class.  Each such window
   //  is contained in one "outer" window (which could be null) and can 
   //  contain several "inner" windows.  An object of class "UINode"
   //  stores a pointer to its outer window, and a list of pointers to 
   //  its inner windows.  When an object of class UINode wants to 
   //  create a new inner window, the signal is connected to UserInterface, 
   //  which does the actual creation.  Quit/close signals are also  
   //  connected to UserInterface, which actually hides and deletes the 
   //  window and all of its inner windows.  Its constructor and destructor
   //  are not public: only UserInterface can create and destroy such 
   //  windows.
   //
   //  We want the quit button and the close window to both transfer
   //  control to UserInterface to do the destruction of the window.
   //  But clicking the quit button produces a signal taking no parameters
   //  and returning void, whereas clicking the "close window" produces
   //  a signal which returns a bool and takes a GdkEventAny* parameter.
   //  Hence, suitable signal handlers "on_quit_button_clicked" and
   //  "on_close_window" are defined in this base class.  Both of these
   //  functions merely emit a "kill_me" signal, which is connected to
   //  the UserInterface::DeleteWindow function.  In this way, we don't
   //  have to do any awkward signal re-bindings.
   //
   // *************************************************************

class UINode : public Gtk::Window
{

 protected:

   UINode(UINode *outerwd, bool deactivate_outer);

   UINode *outer;
   std::list<UINode*> inner;

  // sigc::signal<void,UINode*> kill_me;

   bool m_active;
   bool m_deactivate_outer;

 private:

   UINode(const UINode&);
   UINode& operator=(const UINode&);

 public:

   virtual void show();
   virtual void hide();
   virtual ~UINode();

   void on_quit_button_confirm_clicked();
   void on_quit_button_clicked();
   virtual bool on_close_window(GdkEventAny* event);
   virtual bool on_delete_event(GdkEventAny* event);
   void on_help_about();
   virtual void Deactivate();
   virtual void Activate();

   friend class UserInterface;

};
/*
   // *************************************************************
   // 
   //  HelpDisplay and HelpHandler: Any window that one wishes to
   //  provide help for, such as a "TaskMenu" or "InfoMenu", should
   //  include a "HelpHandler" object in its definition.  An object
   //  of class "HelpHandler" stores a pointer to the UINode it is
   //  associated with, the name of the file containing the help
   //  information, a title to display on the help window, and the
   //  width of the help window, in pixels.  If no help file is given,
   //  the handler assumes that no help is available.  The use of a
   //  HelpHandler is given below:
   //
   //    HelpHandler helper(&outerwindow,helpfile,hsize,title);
   //    helper.Available();  // is help available?
   //    helper.Display();    // displays a help window (if available)
   //    helper.UnDisplay();  // stops help window display
   //    helper.ToggleDisplay();
   //
   //  The help window displayed is actually an object of class
   //  "HelpDisplay" which is derived from UINode.  The user cannot
   //  directly access such an object; it is manipulated only through
   //  a "HelpHandler" object.  The use of the class HelpHandler
   //  ensures that more than one help window cannot be displayed for
   //  a given menu, and the help window is killed when the original
   //  menu is closed.  Also, the display of the help window does NOT
   //  deactivate the original menu.
   //
   // *************************************************************

class HelpHandler
{

    bool m_available;
    bool m_display_status;
    string m_helpfile;
    int m_hsize;
    string m_title;
    UINode *m_outerwd;
    HelpDisplay *m_display;

    static string m_HelpDirectory;
    void Setup(const string& helpfile, int hsize, const string& title);

 public:

    HelpHandler(UINode *outer, const string& helpfile, int hsize, const string& title);
    void Reset(const string& helpfile, int hsize, const string& title);
    virtual ~HelpHandler(){}

    virtual void Deactivate();
    virtual void Activate();
    void Display();
    void UnDisplay();
    void ToggleDisplay();
    bool Available() const { return m_available; }

    virtual void show();
    virtual void hide();

    friend class UserInterface;
    friend class HelpDisplay;

 private:

   HelpHandler(const HelpHandler&);
   HelpHandler& operator=(const HelpHandler&);

};


class HelpDisplay : public UINode
{

 private:

   HelpDisplay(HelpHandler *handler);

 protected:
  
  //Child widgets:
   Gtk::Frame m_frame;
   Gtk::VBox m_VBox;

   Gtk::ScrolledWindow m_ScrolledWindow;
   Gtk::TextView m_TextView;
  
   Gtk::HButtonBox m_ButtonBox;
   Gtk::Button m_hide;

   HelpHandler *m_handler;

 private:

   HelpDisplay(const HelpDisplay&);
   HelpDisplay& operator=(const HelpDisplay&);

 public:

   virtual ~HelpDisplay();
   virtual void show();
   virtual void hide();

   friend class HelpHandler;
   friend class UserInterface;

}; 


   // *************************************************************
   // 
   //  Log Display: Objects of the class "TaskMenu" need to show output
   //  from the various tasks that are performed.  Such output is sent
   //  to an object of class "LogHandler" which stores the information
   //  in a log file.  Output from one given computation is stored in
   //  the file as a "Chapter".  The LogHandler keeps appending chapters
   //  to the log file as more and more computations are done.  However,
   //  the LogHandler can delete past chapters.  An object of class
   //  "LogDisplay" handles the display of the log file to the user.
   //  It displays only one chapter at a time.  It allows the user to
   //  browse all chapters in the file, however.  There are buttons to
   //  go show the next chapter, the previous chapter, or to jump to
   //  an arbitary chapter.  It also has a button to delete a given
   //  chapter, or to delete all chapters.
   //
   // *************************************************************

class LogDisplay : public Gtk::Frame
{

 private:

   LogDisplay(UINode *parent, LogHandler& outlog);

 protected:

   LogHandler& m_log;       // reference to LogHandler
   UINode *m_parent;        // parent window
   int m_current_chapter;

  //Child widgets:
   Gtk::VBox m_VBox;

   Gtk::ScrolledWindow m_ScrolledWindow;
   Gtk::TextView m_TextView;
  
   Gtk::HButtonBox m_ButtonBox;
   Gtk::Button m_jumpto;
   Gtk::Button m_back;
   Gtk::Button m_forward;
   Gtk::Button m_delete;
   Gtk::Button m_clear;

   Gtk::Menu m_Jumpto_Popup;  //  jump-to popup menu

   void on_jumpto_clicked();
   void on_back_clicked();
   void on_forward_clicked();
   void on_delete_clicked();
   void on_clear_clicked();
   void on_jumpto_select(int i);

   void on_new_chapter_done();

 private:

   LogDisplay(const LogDisplay&);
   LogDisplay& operator=(const LogDisplay&);

 public:

   virtual void show();
   virtual void hide();
   virtual ~LogDisplay();

   void DisplayChapter(int index);
   void DisplayLastChapter();

   friend class UINode;
   friend class TaskMenu;
   friend class UserInterface;

};

   // *************************************************************
   // 
   //  ProgressInterrupter: When a long computation is done, one should
   //  pop up an object of class "ProgressInterrupter" to show the user
   //  the progress that is being made, and this also gives the user
   //  the opportunity to pause or cancel the calculation.  To do this,
   //  the long task must be split up into pieces.  Each piece is done
   //  by a function which takes no arguments and returns a real number:
   //  the fraction of the calculation completed during this piece.
   //  The typedef "ComputeEngine" refers to a sigc::slot<Real> which
   //  refers to the function that computes the pieces of the task.
   //
   // *************************************************************

class ProgressInterrupter : public UINode
{

 private:

   ProgressInterrupter(UINode *outerwindow, const ComputeEngine& calculator,
                       bool& completed, const string& title, 
                       bool percentage_mode=true);

   ProgressInterrupter(const ProgressInterrupter&);
   ProgressInterrupter& operator=(const ProgressInterrupter&);

 protected:

        // Signal Handlers:
   bool on_idle();
   void on_cancel_clicked();
   void on_pause_clicked();
   void on_resume_clicked();
   bool on_finished();

      // Member data:
   Gtk::Label m_Title;
   Gtk::VBox m_Box;
   Gtk::Button m_Cancel;
   Gtk::Button m_Pause;
   Gtk::Button m_Resume;
   Gtk::HButtonBox m_hbox;
   Gtk::ProgressBar m_ProgressBar;
   sigc::connection m_engine_connection;
   sigc::connection m_cancel_connection;
   ComputeEngine m_calculator;

   bool m_percentage_mode;
   bool *m_completed;

   virtual ~ProgressInterrupter();
   virtual void show();
   virtual void hide();
   friend class UserInterface;

};

   // *************************************************************
   // 
   //  An object of class "TaskMenu" provides a menu bar, a vertical
   //  list of buttons which, when pressed, perform various tasks or
   //  produce a new list of buttons, and a log displayer.  The
   //  constructor requires an object of class "TaskList" (for details
   //  see "TaskList.h") and an object of class "LogHandler".
   //  This is the main menu for doing the tasks required for data
   //  analysis.
   //
   // *************************************************************

class TaskMenu : public UINode
{

 private:

   TaskMenu(const TaskList& tasklist, LogHandler& outlog);

   void ReDraw(TaskList* tasklist_ptr);

   TaskMenu(const TaskMenu& menu);              // undefined and private
   TaskMenu& operator=(const TaskMenu& menu);   //  to prevent copy

   void MakeButtons(const TaskList& tasklist);
   void ClearButtons();

 protected:

    // Member widgets:
   Gtk::VBox m_outerbox;
   Gtk::Frame m_frame;
   Gtk::VBox m_menubox;
   Gtk::HPaned m_paned;
   Gtk::VButtonBox m_vbox;
   Gtk::Alignment m_align;
   Gtk::Button m_quit;
   sigc::connection m_quit_connection;
   LogDisplay m_logdisplay;
   list<Gtk::Button*> m_tasks;
   HelpHandler m_helper;

   Glib::RefPtr<Gtk::UIManager> m_refUIManager;
   Glib::RefPtr<Gtk::ActionGroup> m_refActionGroup;


 public:

   virtual ~TaskMenu();  
   void Clear();
   virtual void show();
   virtual void hide();
   virtual void Activate();
   virtual void Deactivate();

   friend class UserInterface;

};

   // *****************************************************************
   //
   //  Objects of class "InfoMenu" interact with the end user to obtain
   //  information needed to perform some task.  The constructor requires
   //      an "InfoList" object, 
   //      a boolean variable, 
   //      and a pointer to an outer window (if omitted, the null pointer
   //        is the default value).  
   //  The purpose of an "InfoMenu" object is to assign values to the "InfoItem"
   //  objects pointed to in the "InfoList".  If the "infofile" parameter
   //  is given, this file is read to load initial values of all of the
   //  items in the InfoList.  An "InfoMenu" window is modal,
   //  that is, it disallows interactions with other windows until the quit 
   //  or submit button is pressed.  If the window is closed by the "quit"
   //  or "close" buttons, the boolean parameter "proceed" is set to "false";
   //  if the "submit" button is pressed to close the window, the "proceed"
   //  parameter is set to "true".  However, valid entries must be assigned
   //  or the submission fails and the window persists.  Extra checks can be
   //  performed when the submit button is clicked as well.  This must be
   //  specified in the "InfoList".  Successful submission requires that
   //  any extra checks must succeed. (See "InfoList.h" for more information.)
   //
   //  The constructor is private so that objects of class "InfoMenu" can 
   //  only be created using UserInterface.  In order to create different 
   //  types of menu items, it is necessary to define a base class 
   //  "InfoMenuItem" and several classes, such as "StringMenuItem", etc., 
   //  derived from "InfoMenuItem".  These are described below.
   //
   // ******************************************************************


class InfoMenu : public UINode
{

 private:

   InfoMenu(const InfoList& inflist, bool& proceed, UINode *outerwd,
            const string& infofile="");

   InfoMenu(const InfoMenu& menu);              // undefined and private
   InfoMenu& operator=(const InfoMenu& menu);   //  to prevent copy

 protected:

    // Member widgets:
   Gtk::VBox m_outerbox;
   Gtk::HButtonBox m_bottombox;
   Gtk::Table m_table;
   Gtk::Button m_quit;
   Gtk::Button m_submit;

   bool *m_proceed;
   bool m_extra_on_submit;
   InfoSubmitSlot m_extra_on_submit_slot;

   list<Gtk::Label*> m_labels;
   list<InfoMenuItem*> m_items;
   HelpHandler m_helper;

   Glib::RefPtr<Gtk::UIManager> m_refUIManager;
   Glib::RefPtr<Gtk::ActionGroup> m_refActionGroup;

   void get_info(const string& infofile, bool warn=true);

 public:

   virtual ~InfoMenu();  
   void Clear();
   virtual void show();
   virtual void hide();
   virtual void Activate();
   virtual void Deactivate();
   friend class UserInterface;

   void AddMenuItem(StringInfo& info);
   void AddMenuItem(OutputFileInfo& info);
   void AddMenuItem(InputFileInfo& info);
   void AddMenuItem(TextSelectionSetInfo& info);
   void AddMenuItem(RealInfo& info);

     // callback routines
   void on_submit_clicked();
   void on_save();
   void on_load();

};

   // **************************************************************
   // 
   //  Displaying and interacting with different types of menu items in
   //  "InfoMenu" necessitates polymorphism.  Hence, a base abstract
   //  class "InfoMenuItem" is declared below.  The base class is
   //  derived from Gtk::HBox; hence, each menu item is a horizontal
   //  box (an "HBox").  The base class declares three purely virtual 
   //  members:
   //
   //     bool InfoMenuItem::submit()
   //     bool InfoMenuItem::WriteInfoToStream(ostream& stream) const
   //     bool InfoMenuItem::ReadInfoFromStream(istream& stream)
   //
   //  so each derived class must define such functions.  The "submit"
   //  function must do the following: assign its corresponding "InfoItem"
   //  a valid value, return true if successful, false otherwise.
   //  The base class also defined a "show" member, which calls the
   //  "show_all" command of the HBox so each derived class need not
   //  define a "show" member.  The "WriteInfoToStream" member just
   //  needs to call the "WriteToStream" member of the associated
   //  "InfoItem", and "ReadInfoFromStream" must first call the
   //  "ReadFromStream" of the associated "InfoItem", then put the
   //  assignment into the menu.
   //
   //  Each class derived from "InfoMenuItem" contains (a) a pointer to
   //  an appropriate "InfoItem"; (b) the graphical elements needed to 
   //  query the user for the information to assign to the corresponding
   //  "InfoItem"; and (c) an appropriate "submit" member; (d) an
   //  appropriate "WriteInfoToStream" member; (e) an appropriate
   //  "ReadInfoFromStream" member.
   //
   // *****************************************************************

class InfoMenuItem : public Gtk::HBox
{

 protected:

   InfoMenuItem(InfoMenu *infmenu) : m_menu(infmenu){}
   virtual ~InfoMenuItem(){}
   virtual bool submit()=0;   // purely virtual
   void show();

   InfoMenu *m_menu;
   virtual bool WriteInfoToStream(ostream& stream) const=0;
   virtual bool ReadInfoFromStream(istream& stream)=0;

 private:

   InfoMenuItem(const InfoMenuItem& inf);             // undefined and private
   InfoMenuItem& operator=(const InfoMenuItem& inf);  //  to prevent copy

   friend class InfoMenu;

};



class StringMenuItem : public InfoMenuItem
{

 protected:

   Gtk::Entry m_stringval;
   Gtk::Button m_clear;
   StringInfo *m_infoPtr;

   StringMenuItem(StringInfo *inf, InfoMenu* infmenu);

   virtual bool submit();
   void on_clear_clicked();

   virtual bool WriteInfoToStream(ostream& stream) const;
   virtual bool ReadInfoFromStream(istream& stream);

 private:

   StringMenuItem(const StringMenuItem& inf);     // undefined to prevent copy
   StringMenuItem& operator=(const StringMenuItem& inf);

 public:

   virtual ~StringMenuItem(){}
   friend class InfoMenu;

};


class OutputFileMenuItem : public InfoMenuItem
{

 protected:

   Gtk::Entry m_stringval;
   Gtk::Button m_select;
   OutputFileInfo *m_infoPtr;

   OutputFileMenuItem(OutputFileInfo *inf, InfoMenu* infmenu);

   virtual bool submit();
   void on_select_clicked();

   virtual bool WriteInfoToStream(ostream& stream) const;
   virtual bool ReadInfoFromStream(istream& stream);

 private:

   OutputFileMenuItem(const OutputFileMenuItem& inf);     // undefined to prevent copy
   OutputFileMenuItem& operator=(const OutputFileMenuItem& inf);

 public:

   virtual ~OutputFileMenuItem(){}
   friend class InfoMenu;

};

class InputFileMenuItem : public InfoMenuItem
{

 protected:

   Gtk::Entry m_stringval;
   Gtk::Button m_select;
   InputFileInfo *m_infoPtr;

   InputFileMenuItem(InputFileInfo *inf, InfoMenu* infmenu);

   virtual bool submit();
   void on_select_clicked();

   virtual bool WriteInfoToStream(ostream& stream) const;
   virtual bool ReadInfoFromStream(istream& stream);

 private:

   InputFileMenuItem(const InputFileMenuItem& inf);     // undefined to prevent copy
   InputFileMenuItem& operator=(const InputFileMenuItem& inf);

 public:

   virtual ~InputFileMenuItem(){}
   friend class InfoMenu;

};

class TextSelectionSetMenuItem : public InfoMenuItem
{

 protected:

   Gtk::ComboBoxText m_choices;
   TextSelectionSetInfo *m_infoPtr;

   TextSelectionSetMenuItem(TextSelectionSetInfo *inf, InfoMenu* infmenu);

   virtual bool submit();

   virtual bool WriteInfoToStream(ostream& stream) const;
   virtual bool ReadInfoFromStream(istream& stream);

 private:

   TextSelectionSetMenuItem(const TextSelectionSetMenuItem& inf);     // undefined to prevent copy
   TextSelectionSetMenuItem& operator=(const TextSelectionSetMenuItem& inf);

 public:

   virtual ~TextSelectionSetMenuItem(){}
   friend class InfoMenu;

};

class RealNumberMenuItem : public InfoMenuItem
{

 protected:

   Gtk::SpinButton m_realvalue;
   Gtk::Button m_enter, m_limits;
   //Gtk::Tooltips m_tips;
   RealInfo *m_infoPtr;

   RealNumberMenuItem(RealInfo *inf, InfoMenu* infmenu);

   void on_value_changed();
   void on_enter_clicked();
   void on_limits_clicked();
   virtual bool submit();

   virtual bool WriteInfoToStream(ostream& stream) const;
   virtual bool ReadInfoFromStream(istream& stream);

 private:

   RealNumberMenuItem(const RealNumberMenuItem& inf);     // undefined to prevent copy
   RealNumberMenuItem& operator=(const RealNumberMenuItem& inf);

 public:

   virtual ~RealNumberMenuItem(){}
   friend class InfoMenu;

};
*/
  // *************************************************************

#endif 
