add_library(tasks_libs SHARED single_pivot.cc       
   rolling_pivot.cc       
   stopwatch.cc          
   task_handler.cc       
   task_check.cc         
   task_fit.cc           
   task_get_from_pivot.cc 
   task_obsfunc.cc       
   task_plot.cc          
   task_print.cc         
   task_rotate_corrs.cc  
   task_utils.cc
   xml_handler.cc
   scalar_defs.cc)
target_precompile_headers( tasks_libs PUBLIC args_handler.h 
   log_helper.h     
   multi_compare.h  
   scalar_defs.h    
   #pivoter.h   
   single_pivot.h   
   rolling_pivot.h   
   stopwatch.h      
   task_handler.h   
   task_utils.h     
   xml_handler.h)
