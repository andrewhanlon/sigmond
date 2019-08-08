#include "userinterface_gui.h"
#include <iostream>
using namespace std;

 // ***************************************************************
 //
 //                        UserInterface
 //
 // ***************************************************************

UserInterface::UserInterface()
{
 m_app = Gtk::Application::create("org.gtkmm.sigmond.main");


// static UserInterface singleton;   // don't actually need the singleton 
// return singleton;                 // object

// HelpHandler::m_HelpDirectory=helpdir;
}


void UserInterface::clear()
{
 list<UINode*>::iterator wt;
 for (wt=m_windows.begin();wt!=m_windows.end();wt++) delete (*wt);
 m_windows.clear();
}


UserInterface::~UserInterface()
{
 clear();
}

   // static data members must be declared at file scope
/*
list<UINode*> UserInterface::m_windows;
bool UserInterface::m_initiated=false;


void UserInterface::CreateRun_InfoMenu(const InfoList& inflist,
                                       bool& proceed,
                                       UINode *outerwindow,
                                       const string& infofile) 
{
 if (!m_initiated){
    cerr << "Cannot create any windows until UserInterface::Initiate(argc,argv) called"<<endl;
    exit(1);}
     // create the new window allocating memory
 InfoMenu *pw=new InfoMenu(inflist,proceed,outerwindow,infofile); 
 if (!pw){
    cout << "FATAL ERROR: could not create a window"<<endl; 
    exit(1);}
     // update outer window's list of inner windows
 if (outerwindow!=0){
    outerwindow->inner.push_back(pw);}
     // update handler's list of windows
 m_windows.push_back(pw);
     // show the window
 pw->show();
 if (outerwindow!=0) outerwindow->Deactivate();
 Gtk::Main::run(*pw);
 if (outerwindow!=0) outerwindow->Activate();
}


void UserInterface::CreateRun_TaskMenu(const TaskList& tasklist,
                                       LogHandler& outlog) 
{
 if (!m_initiated){
    cerr << "Cannot create any windows until UserInterface::Initiate(argc,argv) called"<<endl;
    exit(1);}
     // create the new window allocating memory
 TaskMenu *pw=new TaskMenu(tasklist,outlog); 
 if (!pw){
    cout << "FATAL ERROR: could not create a window"<<endl; 
    exit(1);}
     // update handler's list of windows
 m_windows.push_back(pw);
     // show the window
 pw->show();
 Gtk::Main::run(*pw);   // returns when TaskMenu deleted
}


HelpDisplay* UserInterface::CreateShow_HelpDisplay(HelpHandler *handler) 
{
 if (!m_initiated){
    cerr << "Cannot create any windows until UserInterface::Initiate(argc,argv) called"<<endl;
    exit(1);}
 if (handler==0) return 0;
 if (!handler->Available()) return 0;
     // create the new window allocating memory
 HelpDisplay *pw=new HelpDisplay(handler); 
 if (!pw){
    cout << "FATAL ERROR: could not create a window"<<endl; 
    exit(1);}
     // update handler's list of windows
 m_windows.push_back(pw);
     // show the window
 pw->show();
 return pw;
}


void UserInterface::Run_DisplayProgress(UINode *outerwindow,
                                        const ComputeEngine& calculator,
                                        bool& completed,
                                        const string& title,
                                        bool PercentageMode)
{
 if (!m_initiated){
    cerr << "Cannot create any windows until UserInterface::Initiate(argc,argv) called"<<endl;
    exit(1);}
     // create the new window allocating memory
 ProgressInterrupter *pw=new ProgressInterrupter(outerwindow,calculator,
                                                 completed,title,PercentageMode); 
 if (!pw){
    cout << "FATAL ERROR: could not create a window"<<endl; 
    exit(1);}
     // update outer window's list of inner windows
 if (outerwindow!=0){
    outerwindow->inner.push_back(pw);}
     // update handler's list of windows
 m_windows.push_back(pw);
     // show the window
 pw->show();
 if (outerwindow!=0) outerwindow->Deactivate();
 Gtk::Main::run(*pw);  // returns when ProgressInterrupter deleted
 if (outerwindow!=0) outerwindow->Activate();
}


void UserInterface::DeleteWindow(UINode *w)
{
     // destroy inner windows first
 while (w->inner.size()>0) DeleteWindow(w->inner.front());
     // remove from outer window's knowledge
 if (w->outer!=0){
    if (w->m_deactivate_outer) (w->outer)->Activate();
    (w->outer)->inner.remove(w);}
     // remove from handlers list of windows
 m_windows.remove(w);
     // now hide and delete
 w->hide();
 delete w;
 if (!w){
    cerr<<"FATAL error: could not delete a window"<<endl;
    exit(1);}
     // check if all windows are gone and quit if so
 if (m_windows.size()==0){
    Gtk::Main::quit();
    } 
}


void UserInterface::DisplayMessage(UINode& GW, const string& message)
{
 GW.Deactivate();
 Gtk::MessageDialog dialog(GW,message);
 dialog.run();
 GW.Activate();
}

bool UserInterface::Confirm(UINode& GW, const string& prompt)
{
 GW.Deactivate();
 Gtk::MessageDialog dialog(GW,prompt, false, Gtk::MESSAGE_QUESTION, 
                           Gtk::BUTTONS_YES_NO);
 int result = dialog.run();
 GW.Activate();
 if (result==Gtk::RESPONSE_YES) return true;
 else return false;
}

bool UserInterface::GetInputFileName(UINode& GW, string& filename)
{
 return GetInputFileName(GW,filename,"Choose input file");
}

bool UserInterface::GetInputFileName(UINode& GW, string& filename, 
                                     const string& prompt)
{
 GW.Deactivate();
 Gtk::FileChooserDialog dialog(prompt, Gtk::FILE_CHOOSER_ACTION_OPEN);
 dialog.set_transient_for(GW);
 dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
 dialog.add_button("Select", Gtk::RESPONSE_OK);

 int result = dialog.run();
 GW.Activate();
 if (result!=Gtk::RESPONSE_OK) return false;
 filename=dialog.get_filename();
 return true;
}


bool UserInterface::GetOutputFileName(UINode& GW, string& filename)
{
 return UserInterface::GetOutputFileName(GW,filename,"Choose name of output file");
}


bool UserInterface::GetOutputFileName(UINode& GW, string& filename, 
                                      const string& prompt)
{
 GW.Deactivate();
 Gtk::FileChooserDialog dialog(prompt, Gtk::FILE_CHOOSER_ACTION_SAVE);
 dialog.set_transient_for(GW);
 dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
 dialog.add_button("Select", Gtk::RESPONSE_OK);
  
 int result = dialog.run();
 GW.Activate();
 if (result!=Gtk::RESPONSE_OK) return false;
 filename=dialog.get_filename();
 return true;
}


 // ***************************************************************
 //
 //                           UINode
 //
 // ***************************************************************

UINode::UINode(UINode *outerwd, bool deactivate_outer)
          : outer(outerwd), 
            m_active(true), 
            m_deactivate_outer(deactivate_outer)
{
 kill_me.connect(sigc::ptr_fun(UserInterface::DeleteWindow));
}

void UINode::show()
{
 Gtk::Window::show();
}

void UINode::hide()
{
 Gtk::Window::hide();
}

UINode::~UINode(){}

void UINode::on_quit_button_confirm_clicked()
{
 if (UserInterface::Confirm(*this,"Are you sure you want to quit?")) kill_me(this);
}

void UINode::on_quit_button_clicked()
{
 kill_me(this);
}

bool UINode::on_close_window(GdkEventAny* event)
{
 if (m_active) kill_me(this);
 return true;
}

  // When the window is given the "delete_event" signal (this is given by the 
  // window manager, usually by the "close" option, or on the titlebar), the 
  // on_delete_event() slot is called.  If you return false, a "destroy" signal
  // is emitted. Returning true means you don't want the window to be destroyed.

bool UINode::on_delete_event(GdkEventAny* event)
{
 if (m_active) kill_me(this);
 return true;
}

void UINode::on_help_about()
{
 Deactivate();
 Gtk::AboutDialog about;
 about.set_name("SigMonD");
 about.set_border_width(25);
 Gtk::Label line1("Version: 0.1 (August 2007)");
 Gtk::Label line2("Author: Colin Morningstar");
 Gtk::Label line3; 
 line3.set_markup("<span foreground=\"blue\">SIG</span>nal extraction from <span foreground=\"blue\">MON</span>te carlo <span foreground=\"blue\">D</span>ata:");
 Gtk::Label line4("Software for the analysis, fitting, and\nvisualization of data from Monte Carlo\nestimates of correlation functions\nin lattice field theory");
 line4.set_alignment(0.0,0.5);
 about.get_vbox()->add(line1);
 about.get_vbox()->add(line2);
 about.get_vbox()->add(line3);
 about.get_vbox()->add(line4);
 line1.show();
 line2.show();
 line3.show();
 line4.show();
 about.run();
 Activate();
}

void UINode::Activate()
{
 if (m_active) return;
 set_sensitive(true);
 m_active=true;
}

void UINode::Deactivate()
{
 if (!m_active) return;
 set_sensitive(false);
 m_active=false;
}

 // ***************************************************************
 //
 //                 HelpHandler and HelpDisplay
 //
 // ***************************************************************

HelpHandler::HelpHandler(UINode *outer, const string& helpfile, int hsize,
                         const string& title)
{
 m_outerwd=outer;
 Setup(helpfile,hsize,title);
}

void HelpHandler::Setup(const string& helpfile, int hsize,
                        const string& title)
{
 m_title=title;
 m_helpfile=m_HelpDirectory+Trim(helpfile);
 m_hsize=hsize;
    // make sure help file can be opened and read; if not,
    // set helpfile name to null string
 if (!InputFileReadable(m_helpfile)){
    m_available=false;
    m_helpfile.clear();}
 else
    m_available=true;

 m_display_status=false;
 m_display=0;
}

void HelpHandler::Reset(const string& helpfile, int hsize,
                        const string& title)
{
 UnDisplay();
 Setup(helpfile,hsize,title);
}

void HelpHandler::Deactivate()
{
 if (m_display_status) m_display->Deactivate();
}

void HelpHandler::Activate()
{
 if (m_display_status) m_display->Activate();
}

void HelpHandler::Display()
{
 if (m_display_status) return;
 if (!m_available){
    UserInterface::DisplayMessage(*m_outerwd,"Sorry, no help for this menu.");
    return;}
 m_display=UserInterface::CreateShow_HelpDisplay(this); 
 if (m_display!=0) m_display_status=true;
}

void HelpHandler::UnDisplay()
{
 if (!m_display_status) return;
 if (m_display!=0) UserInterface::DeleteWindow(m_display);
 m_display_status=false;
}

void HelpHandler::ToggleDisplay()
{
 if (m_display_status) UnDisplay();
 else Display();
}

void HelpHandler::show()
{
 if (!(m_available && m_display_status)) return;
 m_display->show();
}

void HelpHandler::hide()
{
 if (!m_display_status) return;
 m_display->hide();
}


string HelpHandler::m_HelpDirectory;


    //  This routine facilitates the use of Pango mark up
    //  text in the Help textbuffer.  The string "markup"
    //  is the string to put into the buffer, which may contain
    //  Pango markup.  It is parsed, and then put in the buffer
    //  with appropriate tag information.  Not all of the Markup
    //  language is implemented, but the most useful ones, such
    //  as text color, italics, bold, are implemented.
     

void gtk_text_buffer_set_markup(Gtk::TextView& tview,
                                const Glib::ustring& markup)
{

 Glib::ustring text;
 gunichar accel_char;
 Pango::AttrList attrlist(markup,0,text,accel_char);
 if (!attrlist){
    cerr << "problem"<<endl;
    return;}

 tview.get_buffer()->set_text(text);
 Pango::AttrIter paiter=attrlist.get_iter();
 Gtk::TextBuffer::iterator bit=tview.get_buffer()->begin();
 Gtk::TextBuffer::iterator sit,eit;

 do {
    Pango::Attribute attr;
    int start,end;

    paiter.get_range(start,end);
    if (end>int(text.length())) end=-1;

    sit=bit; sit.forward_chars(start);
    eit=bit; eit.forward_chars(end);

    Glib::RefPtr<Gtk::TextTag> tag=tview.get_buffer()->create_tag();
    tview.get_buffer()->apply_tag(tag,sit,eit);

    attr=paiter.get_attribute(Pango::ATTR_LANGUAGE);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_language()=paiter.get_language().get_string();}

    attr=paiter.get_attribute(Pango::ATTR_FAMILY);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_family()=paiter.get_font_desc().get_family();} 

    attr=paiter.get_attribute(Pango::ATTR_STYLE);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_style()=paiter.get_font_desc().get_style();}

    attr=paiter.get_attribute(Pango::ATTR_WEIGHT);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_weight()=paiter.get_font_desc().get_weight();} 

    attr=paiter.get_attribute(Pango::ATTR_VARIANT);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_variant()=paiter.get_font_desc().get_variant();} 

    attr=paiter.get_attribute(Pango::ATTR_STRETCH);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_stretch()=paiter.get_font_desc().get_stretch();} 

    attr=paiter.get_attribute(Pango::ATTR_SIZE);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_size()=paiter.get_font_desc().get_size();} 

    attr=paiter.get_attribute(Pango::ATTR_FONT_DESC);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_font_desc()=paiter.get_font_desc();} 

    attr=paiter.get_attribute(Pango::ATTR_FOREGROUND);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       Pango::AttrColor pcolor(((PangoAttrColor*)attr.gobj()),true);
       Pango::Color col=pcolor.get_color();
       Gdk::Color gdkcol;
       gdkcol.set_rgb(col.get_red(),col.get_green(),col.get_blue());
       tag->property_foreground_gdk()=gdkcol;}

    attr=paiter.get_attribute(Pango::ATTR_BACKGROUND);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       Pango::AttrColor pcolor(((PangoAttrColor*)attr.gobj()),true);
       Pango::Color col=pcolor.get_color();
       Gdk::Color gdkcol;
       gdkcol.set_rgb(col.get_red(),col.get_green(),col.get_blue());
       tag->property_background_gdk()=gdkcol;}

    attr=paiter.get_attribute(Pango::ATTR_UNDERLINE);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_underline()=Pango::UNDERLINE_LOW;}

    attr=paiter.get_attribute(Pango::ATTR_STRIKETHROUGH);
    if (attr.get_type()!=Pango::ATTR_INVALID){
       tag->property_strikethrough()=true;}

    }
 while (paiter.next());

}
 
// **********************************************************************

HelpDisplay::HelpDisplay(HelpHandler *handler)
           : UINode(handler->m_outerwd,false),
             m_hide("Hide"),
             m_handler(handler)
{

 set_title("Help: "+m_handler->m_title);
 set_border_width(5);
 const int vsize=500;
 set_size_request(m_handler->m_hsize,vsize);

 m_ScrolledWindow.add(m_TextView);
 m_ScrolledWindow.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
 m_frame.set_border_width(5);
 m_frame.set_label("Help");
 m_frame.add(m_ScrolledWindow);
 m_VBox.pack_start(m_frame);

 m_ButtonBox.pack_start(m_hide, Gtk::PACK_SHRINK);
 m_ButtonBox.set_border_width(5);
 m_ButtonBox.set_spacing(5);
 m_ButtonBox.set_layout(Gtk::BUTTONBOX_END);
 m_VBox.pack_start(m_ButtonBox, Gtk::PACK_SHRINK);
 
 stringstream buffer;
 char c;

 ifstream fin(m_handler->m_helpfile.c_str());
 if (!fin){
    buffer.str("Sorry, there is no help for this menu.");
    }
 else{
    while (fin.get(c)) buffer.put(c);   // character by character copy
    }
 fin.close();

 gtk_text_buffer_set_markup(m_TextView, buffer.str());
 buffer.str("");  // erases the buffer
 m_TextView.set_justification(Gtk::JUSTIFY_LEFT);
 m_TextView.set_wrap_mode(Gtk::WRAP_WORD);
 m_TextView.set_editable(false);
 add(m_VBox);

 m_hide.signal_clicked().connect(
        sigc::mem_fun(*m_handler,&HelpHandler::UnDisplay));
}

void HelpDisplay::show()
{
 m_VBox.show();
 m_ScrolledWindow.show();
 m_TextView.show(); 
 m_ButtonBox.show();
 m_hide.show();
 m_frame.show();
 UINode::show(); 
}

void HelpDisplay::hide()
{
 m_VBox.hide();
 m_ScrolledWindow.hide();
 m_TextView.hide(); 
 m_ButtonBox.hide();
 m_hide.hide();
 m_frame.hide();
 UINode::hide();
}

HelpDisplay::~HelpDisplay() 
{
 m_handler->m_display_status=false;
 m_handler->m_display=0;
}


 // ***************************************************************
 //
 //                        LogDisplay
 //
 // ***************************************************************

LogDisplay::LogDisplay(UINode *parent, LogHandler& outlog)
           : m_log(outlog),
             m_parent(parent),
             m_current_chapter(0),
             m_jumpto(Gtk::Stock::JUMP_TO),
             m_back(Gtk::Stock::GO_BACK),
             m_forward(Gtk::Stock::GO_FORWARD),
             m_delete(Gtk::Stock::DELETE),
             m_clear(Gtk::Stock::CLEAR)
{
 set_label(m_log.LogTitle());
 set_border_width(5);
// set_size_request(500, 400);

 add(m_VBox);
 m_TextView.set_editable(false);
 m_ScrolledWindow.add(m_TextView);
 m_ScrolledWindow.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
 m_VBox.pack_start(m_ScrolledWindow);

 m_ButtonBox.pack_start(m_jumpto, Gtk::PACK_SHRINK);
 m_ButtonBox.pack_start(m_back, Gtk::PACK_SHRINK);
 m_ButtonBox.pack_start(m_forward, Gtk::PACK_SHRINK);
 m_ButtonBox.pack_start(m_delete, Gtk::PACK_SHRINK);
 m_ButtonBox.pack_start(m_clear, Gtk::PACK_SHRINK);
 m_ButtonBox.set_border_width(5);
 m_ButtonBox.set_spacing(5);
 m_ButtonBox.set_layout(Gtk::BUTTONBOX_END);
 m_VBox.pack_start(m_ButtonBox, Gtk::PACK_SHRINK);
 
 m_jumpto.signal_clicked().connect(
     sigc::mem_fun(*this,&LogDisplay::on_jumpto_clicked));
 m_back.signal_clicked().connect(
     sigc::mem_fun(*this,&LogDisplay::on_back_clicked));
 m_forward.signal_clicked().connect(
     sigc::mem_fun(*this,&LogDisplay::on_forward_clicked));
 m_delete.signal_clicked().connect(
     sigc::mem_fun(*this,&LogDisplay::on_delete_clicked));
 m_clear.signal_clicked().connect(
     sigc::mem_fun(*this,&LogDisplay::on_clear_clicked));

 m_log.signal_new_chapter_done().connect(
     sigc::mem_fun(*this,&LogDisplay::DisplayLastChapter));
}

void LogDisplay::show()
{
 m_VBox.show();
 m_ScrolledWindow.show();
 m_TextView.show(); 
 m_ButtonBox.show();
 m_jumpto.show();
 m_back.show();
 m_forward.show();
 m_delete.show();
 m_clear.show();
 Gtk::Frame::show();
}

void LogDisplay::hide()
{
 m_VBox.hide();
 m_ScrolledWindow.hide();
 m_TextView.hide(); 
 m_ButtonBox.hide();
 m_jumpto.hide();
 m_back.hide();
 m_forward.hide();
 m_delete.hide();
 m_clear.hide();
 Gtk::Frame::hide();
}

LogDisplay::~LogDisplay()
{
}

void LogDisplay::DisplayChapter(int index)
{
 stringstream buffer;
 if (m_log.ReadChapter(index,buffer)){
    m_current_chapter=index;
    set_label(m_log.LogTitle()+":  Chapter "+IntToString(index)
              +" of "+IntToString(m_log.NumberOfChapters()));}
 else{
    set_label(m_log.LogTitle());
    m_current_chapter=0;}
 m_TextView.get_buffer()->set_text(buffer.str());
 buffer.str("");  // erase temporary buffer
}

void LogDisplay::DisplayLastChapter()
{
 DisplayChapter(m_log.NumberOfChapters());
}

void LogDisplay::on_jumpto_clicked()
{
   //  create the pop up menu
 Gtk::Menu::MenuList& menulist = m_Jumpto_Popup.items();
 menulist.clear();
 for (int i=1;i<=m_log.NumberOfChapters();i++){
    menulist.push_back( Gtk::Menu_Helpers::MenuElem(
         "Chapter "+IntToString(i)+": "+m_log.ChapterTitle(i),
          sigc::bind(sigc::mem_fun(*this, 
          &LogDisplay::on_jumpto_select) ,i)) );}
 m_Jumpto_Popup.popup(1,0);
}

void LogDisplay::on_jumpto_select(int i)
{
 m_Jumpto_Popup.items().clear();
 DisplayChapter(i);
}

void LogDisplay::on_back_clicked()
{
 if (m_current_chapter>1) 
    DisplayChapter(m_current_chapter-1);
}

void LogDisplay::on_forward_clicked()
{
 if (m_current_chapter<m_log.NumberOfChapters()) 
    DisplayChapter(m_current_chapter+1);
}

void LogDisplay::on_delete_clicked()
{
 if (m_log.NumberOfChapters()>0){
    if (UserInterface::Confirm(*m_parent,"Delete current chapter?")){
       m_log.EraseChapter(m_current_chapter);
       if (m_current_chapter>m_log.NumberOfChapters())
          m_current_chapter--;
       DisplayChapter(m_current_chapter);}}
}

void LogDisplay::on_clear_clicked()
{
 if (m_log.NumberOfChapters()>0){
    if (UserInterface::Confirm(*m_parent,"Delete all chapters?")){
       m_log.EraseAllChapters();
       DisplayChapter(0);
       m_current_chapter=0;}}
}


   // *************************************************************
   // 
   //                      ProgressInterrupter
   //
   // *************************************************************

ProgressInterrupter::ProgressInterrupter(UINode *outerwd,
                                         const ComputeEngine& calculator,
                                         bool& completed,
                                         const string& title, 
                                         bool percentage_mode)
                    :  UINode(outerwd,true),           
                       m_Title(title),
                       m_Box(false, 15),
                       m_Cancel(Gtk::Stock::CANCEL),
                       m_Pause("Pause"),
                       m_Resume("Resume"),
                       m_hbox(Gtk::BUTTONBOX_SPREAD,8),
                       m_calculator(calculator),
                       m_percentage_mode(percentage_mode),
                       m_completed(&completed)
{
 set_border_width(15);
 set_title("Progress report");

 m_Box.pack_start(m_Title);
 m_Box.pack_start(m_ProgressBar);
 add(m_Box);
 m_ProgressBar.set_text("Running");
 if (!m_percentage_mode) m_ProgressBar.set_pulse_step(0.025);

 m_hbox.add(m_Cancel);
 m_hbox.add(m_Pause);
 m_hbox.add(m_Resume);
 m_Resume.set_sensitive(false);
 m_Box.pack_start(m_hbox);

 // Connect the signal handlers:
 m_cancel_connection=m_Cancel.signal_clicked().connect(sigc::mem_fun(*this,
                             &ProgressInterrupter::on_cancel_clicked));
 m_Pause.signal_clicked().connect( sigc::mem_fun(*this,
             &ProgressInterrupter::on_pause_clicked) );
 m_Resume.signal_clicked().connect( sigc::mem_fun(*this,
             &ProgressInterrupter::on_resume_clicked) );

     // idle signal handler - called as quickly as possible
 m_engine_connection=Glib::signal_idle().connect( 
         sigc::mem_fun(*this, &ProgressInterrupter::on_idle) );

 *m_completed=false;
}

void ProgressInterrupter::on_cancel_clicked()
{
 bool running=m_Pause.is_sensitive();
 if (running) on_pause_clicked();
 if (UserInterface::Confirm(*this,"Are you sure you want to quit?")) kill_me(this);
 else if (running) on_resume_clicked();
}

void ProgressInterrupter::on_pause_clicked()
{
 m_Pause.set_sensitive(false);
 m_Resume.set_sensitive(true);
 m_engine_connection.disconnect();
 m_ProgressBar.set_text("Paused");
}

void ProgressInterrupter::on_resume_clicked()
{
 m_Pause.set_sensitive(true);
 m_Resume.set_sensitive(false);
 m_engine_connection=Glib::signal_idle().connect( 
         sigc::mem_fun(*this, &ProgressInterrupter::on_idle) );
 m_ProgressBar.set_text("Running");
}

      // This idle callback function is executed as often as possible, 
      // hence it is ideal for processing intensive tasks.
      // The "ComputeEngine" function must take no arguments and return
      // a real number.  If the returned number >= 1.0, then the task
      // is done.  

bool ProgressInterrupter::on_idle()
{
 Real done=m_calculator();    // do a piece of the long computation
 done=(done>0.0)?done:0.0;
 if (m_percentage_mode){
    done+=m_ProgressBar.get_fraction();
    m_ProgressBar.set_fraction((done<1.0)?done:1.0);}
 else{
    m_ProgressBar.pulse();}
 if (done>=1.0){              // computation is done!!
    *m_completed=true;
    m_ProgressBar.set_text("Done");
    m_hbox.remove(m_Pause);
    m_hbox.remove(m_Resume);
    m_Cancel.set_label("Done");
    m_cancel_connection.disconnect();
    m_engine_connection.disconnect();
       //  sit idle to timeout_value milliseconds, then kill progress bar
    int timeout_value=3000;
    m_engine_connection=Glib::signal_timeout().connect(
        sigc::mem_fun(*this,&ProgressInterrupter::on_finished),timeout_value);
    return false;}      // return false when done
 return true;           // return true when not yet done
}

bool ProgressInterrupter::on_finished()
{
 m_engine_connection.disconnect();
 kill_me(this);
 return false;
}

void ProgressInterrupter::show()
{
 m_Title.show();
 m_Box.show();
 m_Cancel.show();
 m_Pause.show();
 m_Resume.show();
 m_hbox.show();
 m_ProgressBar.show();
 UINode::show();
}

void ProgressInterrupter::hide()
{
 m_Title.hide();
 m_Box.hide();
 m_Cancel.hide();
 m_Pause.hide();
 m_Resume.hide();
 m_hbox.hide();
 m_ProgressBar.hide();
 UINode::hide();
}

ProgressInterrupter::~ProgressInterrupter()
{
 hide();
}

 // ***************************************************************
 //
 //                            TaskMenu
 //
 // ***************************************************************

TaskMenu::TaskMenu(const TaskList& tasklist,
                   LogHandler& outlog)
              : UINode(0,false),
                m_frame(tasklist.Title()),
                m_vbox(Gtk::BUTTONBOX_START,15),
                m_align(0.5,1.0,0,0),
                m_quit("Quit this menu"),
                m_logdisplay(this,outlog),
                m_helper(this,tasklist.HelpFile(),
                         tasklist.HelpDisplayWidth(),
                         tasklist.Title())
{
 set_title("SigMonD");
// set_border_width(8);
 set_default_size(1024,600);
// set_resizable(false);

    //  Make the task menu buttons and connect quit menu signal
 m_vbox.set_child_min_width(256);
 m_vbox.set_child_min_height(36);
 m_quit.set_size_request(256,36);
 m_align.add(m_quit);
 m_menubox.set_homogeneous(false);
 MakeButtons(tasklist);
    //  connect close window
 signal_delete_event().connect(sigc::mem_fun(*this,&UINode::on_close_window));

        //Create actions and menus
 m_refActionGroup = Gtk::ActionGroup::create();
 Glib::ustring ui_info = "<ui>"
                         "  <menubar name='MenuBar'>";

 m_refActionGroup->add( Gtk::Action::create("MenuFile", "_File") );
// m_refActionGroup->add( Gtk::Action::create("Open", Gtk::Stock::OPEN),
//                        sigc::mem_fun(*this,&TaskMenu::on_open_new_project));
 m_refActionGroup->add( Gtk::Action::create("Quit", Gtk::Stock::QUIT),
                        sigc::mem_fun(*this,&UINode::on_quit_button_clicked));
 ui_info+="    <menu action='MenuFile'>"
//          "      <menuitem action='Open'/>"
          "      <menuitem action='Quit'/>"
          "    </menu>";

 m_refActionGroup->add( Gtk::Action::create("HelpMenu", "_Help") );
 if (m_helper.Available()){
    m_refActionGroup->add( Gtk::Action::create("Helper", Gtk::Stock::HELP),
                   sigc::mem_fun(m_helper, &HelpHandler::ToggleDisplay));}
 m_refActionGroup->add( Gtk::Action::create("HelpAbout", Gtk::Stock::ABOUT),
                sigc::mem_fun(*this,&UINode::on_help_about));
 ui_info+="    <menu action='HelpMenu'>";
 if (m_helper.Available()){
    ui_info+="      <menuitem action='Helper'/>";}
 ui_info+="      <menuitem action='HelpAbout'/>"
          "    </menu>";
 ui_info+="  </menubar>"
          "</ui>";

 m_refUIManager = Gtk::UIManager::create();
 m_refUIManager->insert_action_group(m_refActionGroup);
 add_accel_group(m_refUIManager->get_accel_group());
   
 m_refUIManager->add_ui_from_string(ui_info);
 Gtk::Widget* pMenuBar = m_refUIManager->get_widget("/MenuBar");
   
      //  Pack menu bar into m_outerbox   
 m_outerbox.pack_start(*pMenuBar, Gtk::PACK_SHRINK);

 m_frame.add(m_menubox);
 m_paned.add1(m_frame);
 m_paned.add2(m_logdisplay);
 m_paned.set_position(300);
 m_outerbox.pack_start(m_paned);
 add(m_outerbox);
}

void TaskMenu::MakeButtons(const TaskList& tasklist)
{
 Gtk::Button *bptr;
 TaskList::const_iterator it;
 for (it=tasklist.begin();it!=tasklist.end();it++){
    bptr=new Gtk::Button("    "+(*it)->Label());
    bptr->set_alignment(0.0,0.5);
    SingleTask *taskptr=GetTask(it);
    if (taskptr!=0){
       bptr->signal_clicked().connect(sigc::bind(taskptr->TaskReference(),this));
       }
    else{
       TaskList *tasklist_ptr=GetTaskList(it);
       if (tasklist_ptr!=0){
          bptr->signal_clicked().connect(sigc::bind(sigc::mem_fun(*this,&TaskMenu::ReDraw),
                                         tasklist_ptr));}}
    m_tasks.push_back(bptr);
    m_vbox.add(*bptr);
    }
 if (tasklist.Up()==0){
    m_quit_connection=m_quit.signal_clicked().connect(
          sigc::mem_fun(*this,&UINode::on_quit_button_clicked));
}
 else{
    m_quit_connection=m_quit.signal_clicked().connect(
          sigc::bind(sigc::mem_fun(*this,&TaskMenu::ReDraw), tasklist.Up()));
    }
 m_menubox.pack_start(m_vbox,Gtk::PACK_SHRINK,15);
 m_menubox.pack_start(m_align,Gtk::PACK_EXPAND_WIDGET,15);
}

void TaskMenu::ClearButtons()
{
 m_menubox.remove(m_vbox);
 m_menubox.remove(m_align);
 list<Gtk::Button*>::const_iterator it;
 for (it=m_tasks.begin();it!=m_tasks.end();it++){
    m_vbox.remove(*(*it));
    delete (*it);}
 m_tasks.clear();
 m_quit_connection.disconnect();
}

void TaskMenu::ReDraw(TaskList* tasklist)
{
 m_frame.set_label(tasklist->Title());
 ClearButtons();
 m_helper.Reset(tasklist->HelpFile(),
                tasklist->HelpDisplayWidth(),
                tasklist->Title());
 MakeButtons(*tasklist);
 show();
}


TaskMenu::~TaskMenu()
{
 hide();
 Clear();
} 

void TaskMenu::Clear()  
{
 ClearButtons();
}


void TaskMenu::show()
{
 list<Gtk::Button*>::const_iterator it;
 for (it=m_tasks.begin();it!=m_tasks.end();it++)
    (*it)->show();
 m_quit.show();
 m_vbox.show();
 m_frame.show();
 m_menubox.show();
 m_align.show();
 m_paned.show();
 m_logdisplay.show();
 m_outerbox.show();
 m_helper.show();
 UINode::show();
}

void TaskMenu::hide()
{
 list<Gtk::Button*>::const_iterator it;
 for (it=m_tasks.begin();it!=m_tasks.end();it++)
    (*it)->hide();
 m_quit.hide();
 m_vbox.hide();
 m_frame.hide();
 m_menubox.hide();
 m_align.hide();
 m_paned.hide();
 m_logdisplay.hide();
 m_outerbox.hide();
 m_helper.hide();
 UINode::hide();
}

void TaskMenu::Deactivate()
{
 UINode::Deactivate();
 m_helper.Deactivate();
}

void TaskMenu::Activate()
{
 UINode::Activate();
 m_helper.Activate();
}


 // ***************************************************************
 //
 //                            InfoMenu
 //
 // ***************************************************************

InfoMenu::InfoMenu(const InfoList& inflist, bool& proceed,
                   UINode *outerwd, const string& infofile)
              : UINode(outerwd,true),
                m_bottombox(Gtk::BUTTONBOX_SPREAD,20),
                m_table(inflist.size(),2,false),
                m_quit("Quit"),
                m_submit("Submit"),
                m_proceed(&proceed),
                m_extra_on_submit(inflist.ExtraOnSubmit()),
                m_extra_on_submit_slot(inflist.ExtraOnSubmitSlot()),
                m_helper(this,inflist.HelpFile(),
                         inflist.HelpDisplayWidth(),
                         inflist.Title())
{
 proceed=false;
 set_title(inflist.Title());
// set_border_width(25);
// set_default_size(400,400);
// set_resizable(false);
// set_modal(true);
// if (outerwd!=0) set_transient_for(*outerwd);

    //  Make the menu labels

 Gtk::Label *lptr;
 list<InfoItem*>::const_iterator it;
 list<Gtk::Label*>::const_iterator lt;
 list<InfoMenuItem*>::iterator mt;

 for (it=inflist.begin();it!=inflist.end();it++){
    lptr=new Gtk::Label((*it)->Label());
    lptr->set_alignment(1.0,0.5);
    m_labels.push_back(lptr);}

   //  Now make the actual menu items

 for (it=inflist.begin();it!=inflist.end();it++){
    (*it)->SetupMenuItem(*this);}


   // Connect signal handlers for bottom buttons

 if (inflist.QuitConfirm())
    m_quit.signal_clicked().connect(sigc::mem_fun(*this,&UINode::on_quit_button_confirm_clicked));
 else
    m_quit.signal_clicked().connect(sigc::mem_fun(*this,&UINode::on_quit_button_clicked));
 signal_delete_event().connect(sigc::mem_fun(*this,&UINode::on_close_window));
 m_submit.signal_clicked().connect(sigc::mem_fun(*this,&InfoMenu::on_submit_clicked));

   // Now put the appropriate widgets into the table (labels, info items)

 m_table.set_col_spacings(8);
 m_table.set_row_spacings(8);
 guint xpadding=5,ypadding=5;
 int hor=1,vert=1;
 for (lt=m_labels.begin();lt!=m_labels.end();lt++){
    m_table.attach(**lt,hor-1,hor,vert-1,vert,Gtk::FILL,Gtk::FILL,
                   xpadding,ypadding);
    vert++;}
 hor=2; vert=1;
 for (mt=m_items.begin();mt!=m_items.end();mt++){
    m_table.attach(**mt,hor-1,hor,vert-1,vert,Gtk::FILL,Gtk::FILL,
                   xpadding,ypadding);
    vert++;}
 m_table.set_border_width(25);

   //  Make the bottom box containing the 
   //  quit, help, save, load, submit  buttons.

 m_bottombox.add(m_quit);
 m_bottombox.add(m_submit);
// m_bottombox.set_border_width(5);
// m_bottombox.set_spacing(5);
// m_bottombox.set_layout(Gtk::BUTTONBOX_END);

 bool use_help=m_helper.Available();
 bool use_file=inflist.LoadSaveEnabled();
 bool use_menu=use_help || use_file;

 if (use_menu){

        //Create actions and menus
    m_refActionGroup = Gtk::ActionGroup::create();
    Glib::ustring ui_info = "<ui>"
                            "  <menubar name='MenuBar'>";

    if (use_file){
       m_refActionGroup->add( Gtk::Action::create("MenuFile", "_File") );
       m_refActionGroup->add( Gtk::Action::create("Open", Gtk::Stock::OPEN),
                sigc::mem_fun(*this,&InfoMenu::on_load));
       m_refActionGroup->add( Gtk::Action::create("Save", Gtk::Stock::SAVE_AS),
                sigc::mem_fun(*this,&InfoMenu::on_save));
       if (inflist.QuitConfirm())
          m_refActionGroup->add( Gtk::Action::create("Quit", Gtk::Stock::QUIT),
                   sigc::mem_fun(*this,&UINode::on_quit_button_confirm_clicked));
       else
          m_refActionGroup->add( Gtk::Action::create("Quit", Gtk::Stock::QUIT),
                   sigc::mem_fun(*this,&UINode::on_quit_button_clicked));
       ui_info+="    <menu action='MenuFile'>"
                "      <menuitem action='Open'/>"
                "      <menuitem action='Save'/>"
                "      <menuitem action='Quit'/>"
                "    </menu>";
       }

    if (use_help){
       m_refActionGroup->add( Gtk::Action::create("HelpMenu", "_Help") );
       m_refActionGroup->add( Gtk::Action::create("Helper", Gtk::Stock::HELP),
                sigc::mem_fun(m_helper, &HelpHandler::ToggleDisplay));
       m_refActionGroup->add( Gtk::Action::create("HelpAbout", Gtk::Stock::ABOUT),
                sigc::mem_fun(*this,&UINode::on_help_about));
       ui_info+="    <menu action='HelpMenu'>"
                "      <menuitem action='Helper'/>"
                "      <menuitem action='HelpAbout'/>"
                "    </menu>";
       }

    ui_info+="  </menubar>"
             "</ui>";

    m_refUIManager = Gtk::UIManager::create();
    m_refUIManager->insert_action_group(m_refActionGroup);
    add_accel_group(m_refUIManager->get_accel_group());
   
    m_refUIManager->add_ui_from_string(ui_info);
    Gtk::Widget* pMenuBar = m_refUIManager->get_widget("/MenuBar");
   
      //  Pack menu bar into m_outerbox   
    m_outerbox.pack_start(*pMenuBar, Gtk::PACK_SHRINK);
    }

      //  Pack table and bottom box into m_outerbox
      //  then put into the window with the "add" command.

 m_outerbox.pack_start(m_table);
 m_outerbox.pack_start(m_bottombox,Gtk::PACK_EXPAND_PADDING,0,20);
// m_paned.pack1(m_outerbox,false,false);
// add(m_paned);
 add(m_outerbox);

 if (infofile.length()>0) get_info(infofile,false);
}

InfoMenu::~InfoMenu()
{
 hide();
 Clear();
} 

void InfoMenu::Clear()  
{
 list<Gtk::Label*>::const_iterator it;
 for (it=m_labels.begin();it!=m_labels.end();it++)
    delete (*it);
 list<InfoMenuItem*>::const_iterator mt;
 for (mt=m_items.begin();mt!=m_items.end();mt++)
    delete (*mt);
 m_labels.clear();
 m_items.clear();
// m_helper.Clear();
}

void InfoMenu::Deactivate()
{
 UINode::Deactivate();
 m_helper.Deactivate();
}

void InfoMenu::Activate()
{
 UINode::Activate();
 m_helper.Activate();
}


void InfoMenu::AddMenuItem(StringInfo& info_item)
{
 StringMenuItem *m_get=new StringMenuItem(&info_item,this);
 m_items.push_back(m_get);
}

void InfoMenu::AddMenuItem(OutputFileInfo& info_item)
{
 OutputFileMenuItem *m_get=new OutputFileMenuItem(&info_item,this);
 m_items.push_back(m_get);
}

void InfoMenu::AddMenuItem(InputFileInfo& info_item)
{
 InputFileMenuItem *m_get=new InputFileMenuItem(&info_item,this);
 m_items.push_back(m_get);
}

void InfoMenu::AddMenuItem(TextSelectionSetInfo& info_item)
{
 TextSelectionSetMenuItem *m_get=new TextSelectionSetMenuItem(&info_item,this);
 m_items.push_back(m_get);
}

void InfoMenu::AddMenuItem(RealInfo& info_item)
{
 RealNumberMenuItem *m_get=new RealNumberMenuItem(&info_item,this);
 m_items.push_back(m_get);
}

void InfoMenu::show()
{
 list<Gtk::Label*>::const_iterator lt;
 list<InfoMenuItem*>::iterator mt;
 for (lt=m_labels.begin();lt!=m_labels.end();lt++)
    (*lt)->show();
 for (mt=m_items.begin();mt!=m_items.end();mt++)
    (*mt)->show();

 m_quit.show();
 m_submit.show();
 m_bottombox.show();
 m_table.show();
 m_outerbox.show();
 m_helper.show();
 UINode::show();
}

void InfoMenu::hide()
{
 list<Gtk::Label*>::const_iterator lt;
 list<InfoMenuItem*>::iterator mt;
 for (lt=m_labels.begin();lt!=m_labels.end();lt++)
    (*lt)->hide();
 for (mt=m_items.begin();mt!=m_items.end();mt++)
    (*mt)->hide();

 m_quit.hide();
 m_submit.hide();
 m_bottombox.hide();
 m_table.hide();
 m_outerbox.hide();
 m_helper.hide();
 UINode::hide();
}

void InfoMenu::on_submit_clicked()
{
 string stext;
 bool flag=true;
 list<InfoMenuItem*>::iterator mt;
 list<Gtk::Label*>::const_iterator lt=m_labels.begin();
 for (mt=m_items.begin();mt!=m_items.end();lt++,mt++){
    bool flag2=(*mt)->submit();
    if (!flag2){
       stext+=(*lt)->get_label()+"\n";
       flag=false;}
    }
 if (!flag){
    Gtk::MessageDialog dialog(*this, "Invalid entries:");
    dialog.set_secondary_text(stext);
    dialog.run();
    return;}
 if (m_extra_on_submit)
    flag=m_extra_on_submit_slot(this);
 if (flag){
    *m_proceed=true;
    kill_me(this);}
}

void InfoMenu::on_load()
{
 string info_file;
 if (!UserInterface::GetInputFileName(*this,
          info_file,"Select the info file")) return;
 get_info(info_file);
}

void InfoMenu::get_info(const string& info_file, bool warn)
{
 if (!InputFileReadable(info_file)){
    if (warn)
       UserInterface::DisplayMessage(*this,
             "Either no read permission or non-existent file");
    return;}
 
 ifstream fin(info_file.c_str());

 string str;
 getline(fin,str);
 bool flag=SubstringFind("InfoMenu Title: "+get_title(),str);
 getline(fin,str); 
 flag=flag && SubstringFind("Number of items: "+IntToString(m_items.size()),str);

 list<InfoMenuItem*>::const_iterator mt;
 for (mt=m_items.begin();mt!=m_items.end();mt++){
    getline(fin,str);
    flag=flag && (*mt)->ReadInfoFromStream(fin);
    }
 flag=flag && IstreamCheck(fin);
 fin.close();
 if (!flag)
    UserInterface::DisplayMessage(*this,"Problem encountered; read failed");
}


void InfoMenu::on_save()
{
 string info_file;
 if (!UserInterface::GetOutputFileName(*this,info_file,
        "Select a name for the info file")) return;

 bool exists,flag;
 flag=OutputFileWritable(info_file,exists);
 if (exists){
    if (flag){
       Gtk::MessageDialog dialog(*this,"File exists; overwrite?",false, 
                                 Gtk::MESSAGE_QUESTION,Gtk::BUTTONS_YES_NO);
       int result=dialog.run();
       if (result!=Gtk::RESPONSE_YES) return;
       }
    else{
       UserInterface::DisplayMessage(*this,"File exists and is write protected");
       return;}
    }

 ofstream fout(info_file.c_str());
 fout << "InfoMenu Title: "<< get_title() <<endl;
 fout << "Number of items: "<< m_items.size()<<endl;
 flag=OstreamCheck(fout);

 list<InfoMenuItem*>::const_iterator mt;
 for (mt=m_items.begin();mt!=m_items.end();mt++){
    fout << endl;
    (*mt)->submit();
    flag=flag && (*mt)->WriteInfoToStream(fout);
    }
 flag=flag && OstreamCheck(fout);
 fout.close();
 if (flag)
    UserInterface::DisplayMessage(*this,"InfoMenu items saved");
 else
    UserInterface::DisplayMessage(*this,"Problem encountered; could not save");
}

 // ***************************************************************
 //
 //                        InfoMenuItem
 //
 // ***************************************************************

void InfoMenuItem::show()
{
 HBox::show_all();
}

 // ***************************************************************
 //
 //                        StringMenuItem
 //
 // ***************************************************************

StringMenuItem::StringMenuItem(StringInfo *inf, InfoMenu* infmenu)
               : InfoMenuItem(infmenu),
                 m_stringval(),
                 m_clear("Clear"),
                 m_infoPtr(inf)
{
 m_stringval.set_width_chars(32);
 if (inf->IsAssigned()) m_stringval.set_text(inf->Value());
 m_clear.signal_clicked().connect( sigc::mem_fun(*this,&StringMenuItem::on_clear_clicked));
 pack_start(m_stringval,Gtk::PACK_EXPAND_WIDGET,6);
 pack_start(m_clear,Gtk::PACK_EXPAND_WIDGET,6);
}

bool StringMenuItem::WriteInfoToStream(ostream& stream) const
{
 return m_infoPtr->WriteValueToStream(stream);
}

bool StringMenuItem::ReadInfoFromStream(istream& stream)
{
 bool flag=m_infoPtr->ReadValueFromStream(stream);
 if (m_infoPtr->IsAssigned()) m_stringval.set_text(m_infoPtr->Value());
 return flag;
}


void StringMenuItem::on_clear_clicked()
{
 m_stringval.set_text("");
}

bool StringMenuItem::submit()
{
 m_infoPtr->Assign(m_stringval.get_text());
 return true;
}


 // ***************************************************************
 //
 //                        OutputFileMenuItem
 //
 // ***************************************************************

OutputFileMenuItem::OutputFileMenuItem(OutputFileInfo *inf, InfoMenu* infmenu)
               : InfoMenuItem(infmenu),
                 m_stringval(),
                 m_select("Select"),
                 m_infoPtr(inf)
{
 m_stringval.set_width_chars(32);
 if (inf->IsAssigned()) m_stringval.set_text(inf->Value());
// m_stringval.set_editable(false);
 m_select.signal_clicked().connect( 
      sigc::mem_fun(*this,&OutputFileMenuItem::on_select_clicked));

 pack_start(m_stringval,Gtk::PACK_EXPAND_WIDGET,6);
 pack_start(m_select,Gtk::PACK_EXPAND_WIDGET,6);
}

bool OutputFileMenuItem::WriteInfoToStream(ostream& stream) const
{
 return m_infoPtr->WriteValueToStream(stream);
}

bool OutputFileMenuItem::ReadInfoFromStream(istream& stream)
{
 bool flag=m_infoPtr->ReadValueFromStream(stream);
 if (m_infoPtr->IsAssigned()) m_stringval.set_text(m_infoPtr->Value());
 return flag;
}


void OutputFileMenuItem::on_select_clicked()
{
 string filename;
 if (UserInterface::GetOutputFileName(*m_menu,filename))
    m_stringval.set_text(filename);
}

bool OutputFileMenuItem::submit()
{
 m_infoPtr->Assign(m_stringval.get_text());
 return m_infoPtr->IsAssigned();
}





 // ***************************************************************
 //
 //                        InputFileMenuItem
 //
 // ***************************************************************

InputFileMenuItem::InputFileMenuItem(InputFileInfo *inf, InfoMenu* infmenu)
               : InfoMenuItem(infmenu),
                 m_stringval(),
                 m_select("Select"),
                 m_infoPtr(inf)
{
 m_stringval.set_width_chars(32);
 if (inf->IsAssigned()) m_stringval.set_text(inf->Value());
// m_stringval.set_editable(false);
 m_select.signal_clicked().connect( 
     sigc::mem_fun(*this,&InputFileMenuItem::on_select_clicked));

 pack_start(m_stringval,Gtk::PACK_EXPAND_WIDGET,6);
 pack_start(m_select,Gtk::PACK_EXPAND_WIDGET,6);
}

bool InputFileMenuItem::WriteInfoToStream(ostream& stream) const
{
 return m_infoPtr->WriteValueToStream(stream);
}

bool InputFileMenuItem::ReadInfoFromStream(istream& stream)
{
 bool flag=m_infoPtr->ReadValueFromStream(stream);
 if (m_infoPtr->IsAssigned()) m_stringval.set_text(m_infoPtr->Value());
 return flag;
}


void InputFileMenuItem::on_select_clicked()
{
 string filename;
 if (UserInterface::GetInputFileName(*m_menu,filename))
    m_stringval.set_text(filename);
}

bool InputFileMenuItem::submit()
{
 m_infoPtr->Assign(m_stringval.get_text());
 return m_infoPtr->IsAssigned();
}




 // ***************************************************************
 //
 //                        TextSelectionSetMenuItem
 //
 // ***************************************************************

TextSelectionSetMenuItem::TextSelectionSetMenuItem(TextSelectionSetInfo *inf, 
                                                   InfoMenu* infmenu)
                         : InfoMenuItem(infmenu),
                           m_infoPtr(inf)
{
     // Create the combobox
 TextSelectionSetInfo::const_iterator ct;
 for (ct=m_infoPtr->begin();ct!=m_infoPtr->end();ct++)
    m_choices.append_text(*ct);
     // Set active choice if assigned
 if (inf->IsAssigned()) m_choices.set_active_text(inf->Value());
     // Put combobox in HBox
 pack_start(m_choices,Gtk::PACK_SHRINK,6);
}

bool TextSelectionSetMenuItem::WriteInfoToStream(ostream& stream) const
{
 return m_infoPtr->WriteValueToStream(stream);
}

bool TextSelectionSetMenuItem::ReadInfoFromStream(istream& stream)
{
 bool flag=m_infoPtr->ReadValueFromStream(stream);
 if (m_infoPtr->IsAssigned()) m_choices.set_active_text(m_infoPtr->Value());
 return flag;
}

bool TextSelectionSetMenuItem::submit()
{
 m_infoPtr->Assign(m_choices.get_active_text());
 return m_infoPtr->IsAssigned();
}

 // ***************************************************************
 //
 //                      RealNumberMenuItem
 //
 // ***************************************************************

RealNumberMenuItem::RealNumberMenuItem(RealInfo *inf, InfoMenu* infmenu)
               : InfoMenuItem(infmenu),
                 m_realvalue(1.0,inf->Precision()),
                 m_enter("Enter"),
                 m_limits("Limits"),
                 m_infoPtr(inf)
{
// m_realvalue.set_numeric(true);
 m_realvalue.set_size_request(128, -1);
 m_realvalue.set_range(-10000.0,10000);
 Real increment=1.0;
 for (uint i=0;i<inf->Precision();i++) increment/=10.0;
 m_realvalue.set_increments(increment,10.0*increment);
// m_realvalue.set_update_policy(Gtk::UPDATE_ALWAYS);

 m_realvalue.set_text(" ");

 m_realvalue.signal_value_changed().connect(sigc::mem_fun(*this,
             &RealNumberMenuItem::on_value_changed));
 m_enter.signal_clicked().connect(sigc::mem_fun(*this,
             &RealNumberMenuItem::on_enter_clicked));
 m_limits.signal_clicked().connect(sigc::mem_fun(*this,
             &RealNumberMenuItem::on_limits_clicked));
       
// m_realvalue.set_digits(inf->Precision());
// m_stringval.set_width_chars(32);
// if (inf->IsAssigned()) m_stringval.set_text(inf->Value());
// m_clear.signal_clicked().connect( sigc::mem_fun(*this,&RealNumberMenuItem::on_clear_clicked));
// pack_start(m_stringval,Gtk::PACK_EXPAND_WIDGET,6);
 pack_start(m_realvalue,Gtk::PACK_SHRINK,6);
 pack_start(m_enter,Gtk::PACK_SHRINK,6);
 pack_start(m_limits,Gtk::PACK_SHRINK,6);
}

bool RealNumberMenuItem::WriteInfoToStream(ostream& stream) const
{
 return m_infoPtr->WriteValueToStream(stream);
}

bool RealNumberMenuItem::ReadInfoFromStream(istream& stream)
{
 bool flag=m_infoPtr->ReadValueFromStream(stream);
// if (m_infoPtr->IsAssigned()) m_stringval.set_text(m_infoPtr->Value());
 return flag;
}


void RealNumberMenuItem::on_value_changed()
{
 cout << "value changed"<<endl;
 m_realvalue.update();
 double vmin,vmax,val;
 m_realvalue.get_range(vmin,vmax);
 cout << "vmin="<<vmin<<"  vmax="<<vmax<<endl;
 val=m_realvalue.get_value();
 cout << "val = "<<val<<endl;
 if (val>vmax){
    m_realvalue.set_value(vmax);m_realvalue.set_text("too large");
    m_realvalue.update();
   }
 else if (val<vmin){
    m_realvalue.set_value(vmin);
    m_realvalue.update();}
}

void RealNumberMenuItem::on_enter_clicked()
{
 cout << "enter clicked"<<endl;
 m_realvalue.update();
}

void RealNumberMenuItem::on_limits_clicked()
{
 cout << "limits clicked"<<endl;
}

bool RealNumberMenuItem::submit()
{
// m_infoPtr->Assign(m_stringval.get_text());
 return true;
}
*/
// ************************************************************************

