#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H


  // In batch mode, this user interface is not used, so this
  // is just a dummy implementation.  It is meant for interactive
  // mode.


class UserInterface
{

 public:

   UserInterface() {}
   UserInterface& operator=(const UserInterface& ui)
    {return *this;}

   void pressEnterToContinue();

};


// **************************************************************
#endif
