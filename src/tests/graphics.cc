#include "fx.h"
#include "glviewWindow.h"
#include <stdio.h>
#include <iostream>
#include <fstream>


int main(int argc,char *argv[])
{
  cout << "I am main for frontend!\n";  
  // Make application
  FXApp* application=new FXApp("GLViewer","Test");
  
  // Open the display
  application->init(argc,argv);

  // Make window
  new GLViewWindow(application);
  
  // Create the application's windows
  application->create();
  
  // Run the application
  application->run();
  
  return 0;
}
