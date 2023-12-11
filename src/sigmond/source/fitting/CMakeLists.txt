add_library(fitting SHARED chisq_anisotropy.cc 
   chisq_base.cc 
   chisq_disp.cc     
   chisq_fit.cc           
   chisq_tcorr.cc 
   chisq_logtcorr.cc 
   lmder.cc       
   minimizer.cc 
   minpack.cc 
   model_logtcorr.cc 
   model_tcorr.cc 
   nl2sol.cc )
target_precompile_headers(fitting PUBLIC chisq_anisotropy.h    
   chisq_base.h            
   chisq_disp.h      
   chisq_fit.h             
   chisq_tcorr.h           
   minimizer.h             
   minpack.h               
   model_logtcorr.h        
   model_tcorr.h)
