###########################################
# MIMICS repeat run function
###########################################

MIMrepeat <- function(forcing_df, rparams, output_type = "summary") {
  
  # Set global model parameters
  .GlobalEnv$Vslope = Vslope_default * rparams$Vslope_x[1]
  .GlobalEnv$Vint = Vint_default * rparams$Vint_x[1]
  .GlobalEnv$Kslope = Kslope_default * rparams$Kslope_x[1]
  .GlobalEnv$Kint = Kint_default * rparams$Kint_x[1]
  # .GlobalEnv$Tau_MULT = Tau_MULT_default * rparams$Tau_x[1]
  # .GlobalEnv$CUE = CUE_default * rparams$CUE_x[1]
  # .GlobalEnv$desorb_MULT = desorb_MULT_default * rparams$desorb_x[1]
  # .GlobalEnv$fPHYS_MULT = fPHYS_MULT_default * rparams$fPHYS_x[1]
  .GlobalEnv$avGWC = rparams$avGWC[1]
  .GlobalEnv$akGWC = rparams$akGWC[1]
  
  #full run of forcing data csv
  MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% map(MIMICS_INC_HOUR) %>% bind_rows() 

  #Optional combine MIMout with forcing data
  MIMrun <- forcing_df %>% cbind(MIMrun %>% select(-Site, -TSOI))
  
  #add run number
  MIMrun$run_num <- rparams$run_num[1]
  
  # return MIMrun
  return(MIMrun)
}

