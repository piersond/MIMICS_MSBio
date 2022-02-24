## Set working drive to repo parent
setwd("C:/github/MIMICS_MSBio")

# Load packages
library(dplyr) # used for pipes (%>%)
library(purrr) # used for map() function

# Load MIMICS model ftn script
source("MIMICS_ftns/MIMICS_24hr.R")

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
MIMout <- MIMICS_24(df[1,])


### Multi-point run
df10 <- df[rep(seq_len(nrow(df)), each = 10),] # Create example dataframe with 10 replicate rows
df10$SITE <-  sprintf("SITE%d",seq(1:10)) # Set unique SITE names

# Run MIMICS using each row of SITE data and row bind the model output into one dataframe
MIMout2 <- df10 %>% split(1:nrow(df10)) %>% map(~ MIMICS_24(df=.)) %>% bind_rows() 



