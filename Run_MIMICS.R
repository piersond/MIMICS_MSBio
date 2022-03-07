## Set working drive to repo parent
setwd("C:/github/MIMICS_MSBio")

# Load packages
library(dplyr) # used for pipes (%>%)
library(purrr) # used for map() function

# Load MIMICS model ftn script
source("MIMICS_ftns/MIMICS_INC_hourly.R")

# Load parameters for MIMICS
source("MIMICS_set_parameters.R")


#########################
# Example MIMICS runs
#########################

### Single-point run
df <- data.frame(SITE = 'HARV',
                   ANPP = 750,
                   MAT = 25,
                   CLAY = 15,
                   LIG = 20,
                   N = 1,
                   CN = 49)
# Run MIMICS and store output as MIMout 
MIMout <- MIMICS_INC_HOUR(df[1,], ndays=200) #hourly step 

### Multi-point run (hourly)
df50 <- df[rep(seq_len(nrow(df)), each = 50),] # Create example dataframe with 10 replicate rows
df50$SITE <-  sprintf("SITE%d",seq(1:50)) # Set unique SITE names

# Run MIMICS using each row of SITE data and row bind the model output into one dataframe
MIMout_50 <- df50 %>% split(1:nrow(df50)) %>% map(~ MIMICS_INC_HOUR(df=., ndays=200)) %>% bind_rows()




