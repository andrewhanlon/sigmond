
add_library(data_handling_libs bins_info.cc 
	byte_handler.cc 
	corr_data_handler.cc 
	ensemble_info.cc 
	filelist_info.cc 
	io_handler_fstream.cc 
	io_handler_hdf5.cc 
	obs_get_handler.cc 
	vev_data_handler.cc)
target_precompile_headers(data_handling_libs PUBLIC array.h    
   bins_handler.h         
   bins_info.h            
   byte_handler.h         
   corr_data_handler.h    
   data_io_handler.h      
   ensemble_info.h        
   filelist_info.h        
   io_handler_fstream.h   
   io_handler_hdf5.h      
   io_map.h               
   io_map_base.h          
   io_map_fstream.h       
   io_map_hdf5.h          
   laph_data_io_handler.h 
   obs_get_handler.h      
   samplings_handler.h    
   vev_data_handler.h)
