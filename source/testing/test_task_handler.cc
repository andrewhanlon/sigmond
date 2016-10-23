#include "task_handler.h"
#include <cstdio>
#include <ctime>
#include <map>

using namespace std;



void testTaskHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestTaskHandler")==0)
 return;

 cout << endl << "Starting test_task_handler"<<endl;

 TaskHandler TH(xml_in); 

}
