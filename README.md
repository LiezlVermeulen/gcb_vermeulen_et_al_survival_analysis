# gcb_vermeulen_et_al_survival_analysis
Code accompanying the paper:  L.M. Vermeulen, B. Verbist, K. Van Meerbeek, J. Slingsby, P.N. Bernardino, B. Somers. 2024. Wetness severity increases abrupt shifts in ecosystem functioning in arid savannas.

 calc_climate.R: script for calculating the annual drought and wetness severity per pixel
This script allows calculating annual drought and wetness severity, to be used in the surivvial analysis model. The main input is a raster stack of rainfall data e.g. CHIRPS. 

survival_analysis.R: script for performing survivial analysis 
This script allows performing survivial analysis.  Static explanatory variables (e.g. ecosystem characteristics) should be provided in raster format. Temporal explanatory variables (e.g. stack of annual drought or wetness severity) should be provided as raster stacks. Binary event response variable (e.g. occurrence of breakpoint or no-breakpoint) should be provided in both raster and shapefile format, with a binary attribute/value indicating whether the event has happened (1) or not (0). 

extract_breakpoints.R: script for detecting and characterising breakpoints
This script allows detecting and characterising breakpoints in a rain us efficiency (RaUE/RUE) time series using BFAST01  for the purpose of finding abrupt shifts in ecosystem functioning. The main input is a raster stack of RaUE/RUE time series.  

