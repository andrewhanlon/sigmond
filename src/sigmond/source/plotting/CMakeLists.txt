add_library(plotting SHARED create_plots.cc grace_plot.cc)
target_precompile_headers(plotting PUBLIC create_plots.h 
   grace_plot.h)
