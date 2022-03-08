# Set working directory
setwd("C:/github/MIMICS_MSBio")

# Bring in MIMICS RXEQ function
source("MIMICS_ftns/RXEQ_ftn.R")

# Data table input to the MIMICS_24 ftn requires the following data 
# ..with matching column names

# "SITE"  Identifier for the sampling location
# "MAT"   Mean annual temperature (or temp of incubation) in deg C
# "CLAY"  Soil clay content in percent
# "LIG"   Lignin content of the litter 
# "N"     N content of the litter
# "CN"    C:N ratio of the litter 

###########################################
# MIMICS single point function
###########################################
MIMICS_INC_DAILY <- function(df, ndays){
  
  #Convert parameter rates from hourly to daily
  Vslope <- Vslope
  Vint <- Vint
  tau_r <- tau_r
  tau_K <- tau_K
  desorb_MULT <- desorb_MULT 

  ###Bring in forcing data
  ANPP       <- df$ANPP/2 # Estimate SOM inputs to be 50% of ANPP
  fCLAY      <- df$CLAY/100  # Convert clay from percent to fraction
  TSOI       <- df$MAT # Use MAT as proxy for soil temperature
  lig_N      <- df$lig_N 
  
  ### Set fMET (=partioning coefficient for metabolic vs. structural litter pools)
  # Defualt fMET equation using lig:N values
  fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
  
  #Debug fMET
  .GlobalEnv$fMET<- fMET #store as global to check value after run
  
  ############################################################
  # MIMICS MODEL PARAMETERIZATION STARTS HERE
  ############################################################
  
  # Calc litter input rate
  EST_LIT <- ANPP/(365*24) # Input C per day
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3]
    
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))    
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT    
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))               
  
  desorb <- desorb * desorb_MULT  
  fPHYS <- fPHYS * fPHYS_MULT
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  #MC scalars on VMAX and Km
  VMAX <- VMAX * VMAX_MULT
  KM <- KM * KM_MULT
  
  
  ############################################################
  # MIMICS MODEL RUN STARTS HERE
  ############################################################
 
  # Set run length
  nday   <- ndays # number of days
   
  # Open matrices to store model output
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- rep(NA, dim=6)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  # Create dataframe to store model output
  MIMout <- data.frame(SITE = rep(df$SITE, nday),
                       DAY = rep(NA, nday),
                       LITm = rep(NA, nday),
                       LITs = rep(NA, nday),
                       MICr = rep(NA, nday),
                       MICK = rep(NA, nday),
                       SOMp = rep(NA, nday),
                       SOMc = rep(NA, nday),
                       SOMa = rep(NA, nday),
                       CO2_MICr = rep(NA, nday),
                       CO2_MICK = rep(NA, nday))
  
  # Initialize model pools and fluxes (same as sandbox)
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  lit     <- I   
  mic     <- I  
  som     <- rep(NA, 3) 
  som[1]  <- I[1]
  som[2]  <- I[2]
  som[3]  <- I[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  ## Set global parameters to pass to RXEQ function
  .GlobalEnv$VMAX <- VMAX 
  .GlobalEnv$KM <- KM 
  .GlobalEnv$fPHYS <- fPHYS
  .GlobalEnv$fCHEM <- fCHEM
  .GlobalEnv$fAVAI <- fAVAI
  .GlobalEnv$I <- I
  .GlobalEnv$tau <- tau 
  .GlobalEnv$LITmin <- LITmin
  .GlobalEnv$SOMmin <- SOMmin
  .GlobalEnv$MICtrn <- MICtrn
  .GlobalEnv$desorb <- desorb 
  .GlobalEnv$DEsorb <- DEsorb
  .GlobalEnv$OXIDAT <- OXIDAT
  
  # Create vector of parameter values
  tpars <- c(I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
             fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
             tau   = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
             desorb= desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
  
  # Set initial pool values
  LIT_1  <- lit[1]
  LIT_2  <- lit[2]
  MIC_1  <- mic[1]
  MIC_2  <- mic[2]
  SOM_1  <- som[1]
  SOM_2  <- som[2]
  SOM_3  <- som[3]
  CO2_1  <- 0 
  CO2_2  <- 0
  
  ### BEGIN MODEL LOOP ###
  # Loop over 24 hours for specified number of days 
  for (d in 1:nday)  {
      
      # Get model output from RXEQ ftn
      update <- RXEQ(y = c( LIT_1 = LIT_1, LIT_2 = LIT_2, 
                            MIC_1 = MIC_1 , MIC_2 = MIC_2, 
                            SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3),
                     pars = tpars)
      
      # Update C pools
      LIT_1  <- LIT_1 + update[[1]][1]
      LIT_2  <- LIT_2 + update[[1]][2]
      MIC_1  <- MIC_1 + update[[1]][3]
      MIC_2  <- MIC_2 + update[[1]][4]
      SOM_1  <- SOM_1 + update[[1]][5]
      SOM_2  <- SOM_2 + update[[1]][6]
      SOM_3  <- SOM_3 + update[[1]][7]
      CO2_1  <- CO2_1 + update[[1]][8] #daily rate, not cumulative
      CO2_2  <- CO2_2 + update[[1]][9] #daily rate, not cumulative
      #remove(UPpars, UPy, update)
      
      # Store delta
      MIMout$DAY[d] <- d
      MIMout$dLITm[d] <- update[[1]][1]
      MIMout$dLITs[d] <- update[[1]][2]
      MIMout$dMICr[d] <- update[[1]][3]
      MIMout$dMICK[d] <- update[[1]][4]
      MIMout$dSOMp[d] <- update[[1]][5]
      MIMout$dSOMc[d] <- update[[1]][6]
      MIMout$dSOMa[d] <- update[[1]][7]
      MIMout$dCO2_MICr[d] <- update[[1]][8]
      MIMout$dCO2_MICK[d] <- update[[1]][9]
      
      # Store daily pool size
      MIMout$LITm[d] <- LIT_1
      MIMout$LITs[d] <- LIT_2
      MIMout$MICr[d] <- MIC_1
      MIMout$MICK[d] <- MIC_2
      MIMout$SOMp[d] <- SOM_1
      MIMout$SOMc[d] <- SOM_2
      MIMout$SOMa[d] <- SOM_3
      MIMout$CO2_MICr[d] <- CO2_1
      MIMout$CO2_MICK[d] <- CO2_2
      
      #Print tracker
      print(d)
      
  }	#close daily loop	

  # Return daily model output as datatable
  return(MIMout)
  
} 

#####################
# Example use of 
#####################

# ##############################################
# #single point run
# ##############################################
df <- data.frame(SITE = 'HARV',
                   ANPP = 744,
                   MAT = 7.1,#7.1,
                   CLAY = 15,
                   LIG = 21,
                   N = 0.75,
                   lig_N = 21/0.75,
                   CN = 49.01960784)

# Load model parameters
source("MIMICS_set_parameters.R")

# Run daily model function
MIMout <- MIMICS_INC_DAILY(df[1,], 30000)


# ##############################################
# # Full forcing dataset run
# ##############################################
# data <- data <- read.csv("RCrk_Modelling_Data/LTER_SITE_1.csv", as.is=T)

# MIMrun <- data %>% split(1:nrow(data)) %>% map(~ MIMICS1(df=.)) %>% bind_rows()
# MIMrun <- data %>% cbind(MIMrun %>% select(-Site, -TSOI))


############
# Plots
###########
library(ggplot2)
library(ggpubr)

# Litter mass
plot_LIT <- ggplot(MIMout, aes(y=LITs, x=DAY, color="Structural")) + geom_line(size=1) +
  geom_line(aes(y=LITm, x=DAY, color="Metabolic"), size=1) +
  theme_bw() +
  ylab("Litter C") +
  xlab("Time (days)") +
  labs(color = "Litter Pool")

# SOM & MIC pools
plot_SOM <- ggplot(MIMout, aes(y=SOMc, x=DAY, color="SOMc")) + geom_line(size=1) +
  geom_line(aes(y=SOMa, x=DAY, color="SOMa"), size=1) +
  geom_line(aes(y=SOMp, x=DAY, color="SOMp"), size=1) +
  theme_bw() +
  ylab("SOM") +
  xlab("Time (days)") +
  labs(color = "SOC Pool")

# SOM & MIC pools
plot_MIC <- ggplot(MIMout, aes(y=MICK, x=DAY, color="MIC-K")) + geom_line(size=1) +
  geom_line(aes(y=MICr, x=DAY, color="MIC-r"), size=1) +
  theme_bw() +
  ylab("Microbial C") +
  xlab("Time (days)") +
  labs(color = "MIC-C Pool")

# CO2 fraction
plot_CO2 <- ggplot(MIMout, aes(y=CO2_MICK,
                   x=DAY, color="dMIC-K")) + geom_line(size=1) +
  geom_line(aes(y=CO2_MICr, x=DAY, color="dMIC-r"), size=1) +
  theme_bw() +
  ylab("Respiration Total (CO2-C)") +
  xlab("Time (days)") +
  labs(color = "MIC")



# Build a panel plot
ggarrange(plot_LIT, plot_SOM, plot_MIC, plot_CO2,
          nrow=4,
          ncol=1)

####### Rate plots ################
# Litter mass
plot_dLIT <- ggplot(MIMout, aes(y=dLITs, x=DAY, color="Structural")) + geom_line(size=1) +
  geom_line(aes(y=dLITm, x=DAY, color="Metabolic"), size=1) +
  theme_bw() +
  ylab("Litter C") +
  xlab("Time (days)") +
  labs(color = "Litter Pool")

# SOM & MIC pools
plot_dSOM <- ggplot(MIMout, aes(y=dSOMc, x=DAY, color="SOMc")) + geom_line(size=1) +
  geom_line(aes(y=dSOMa, x=DAY, color="SOMa"), size=1) +
  geom_line(aes(y=dSOMp, x=DAY, color="SOMp"), size=1) +
  theme_bw() +
  ylab("SOM") +
  xlab("Time (days)") +
  labs(color = "SOC Pool")

# SOM & MIC pools
plot_dMIC <- ggplot(MIMout, aes(y=dMICK, x=DAY, color="dMIC-K")) + geom_line(size=1) +
  geom_line(aes(y=dMICr, x=DAY, color="dMIC-r"), size=1) +
  theme_bw() +
  ylab("Microbial C") +
  xlab("Time (days)") +
  labs(color = "MIC-C Pool")

# CO2 fraction
plot_dCO2 <- ggplot(MIMout, aes(y=dCO2_MICK,
                               x=DAY, color="MIC-K")) + geom_line(size=1) +
  geom_line(aes(y=dCO2_MICr, x=DAY, color="MIC-r"), size=1) +
  theme_bw() +
  ylab("Respiration (CO2-C d-1)") +
  xlab("Time (days)") +
  labs(color = "MIC")



# Build a panel plot
ggarrange(plot_LIT, plot_dLIT, 
          plot_SOM, plot_dSOM,
          plot_MIC, plot_dMIC, 
          plot_CO2, plot_dCO2,
          nrow=4,
          ncol=2)

