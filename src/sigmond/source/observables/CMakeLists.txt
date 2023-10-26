add_library(observables_libs  basic_laph_operator_info.cc 
   diag_corr_set.cc 
   correlator_info.cc 
   correlator_matrix_info.cc 
   encoder.cc  
   gen_irrep_operator_info.cc 
   mcobs_info.cc 
   momenta.cc 
   operator_info.cc)
target_precompile_headers(observables_libs PUBLIC basic_laph_operator_info.h  
   diag_corr_set.h             
   correlator_info.h           
   correlator_matrix_info.h    
   encoder.h                   
   gen_irrep_operator_info.h   
   mcobs_info.h                
   momenta.h                   
   operator_info.h)
