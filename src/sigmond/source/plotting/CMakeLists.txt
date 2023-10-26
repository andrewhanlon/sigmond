add_library(plotting_libs create_plots.cc grace_plot.cc)
target_precompile_headers(plotting_libs PUBLIC create_plots.h 
   grace_plot.h)
