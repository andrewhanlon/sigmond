add_library(analysis_libs SHARED bootstrapper.cc histogram.cc matrix.cc mc_estimate.cc mcobs_handler.cc sampling_info.cc)
target_precompile_headers(analysis_libs PUBLIC bootstrapper.h       
   histogram.h          
   matrix.h             
   mc_estimate.h        
   mcobs_handler.h      
   sampling_info.h)
target_compile_definitions(analysis_libs PUBLIC analysis_libs)
