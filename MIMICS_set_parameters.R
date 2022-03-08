########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.0000000125, 6)
Kslope  <- rep(0.02, 6)
Kint    <- rep(3.19, 6)
aK      <- rep(0.15625, 6)
vMOD    <- c(2, 0.4, 2, 0.6, 0.6, 0.4) * 5
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.5, 0.25, 0.7, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.6, 1.3, 3.5)
Tau_MULT <- 1
fPHYS_r <- c(0.15, 1.3)
fPHYS_K <- c(0.1, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.00002, -4.5)
PHYS_scalar <- c(3, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(0.5, 0.85, 0.013)
depth <- 20 # Set soil depth

#Set required multipliers to 1 (e.g. use default)
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1
VMAX_MULT = 1
KM_MULT = 1

########################################
# Apply parameter multipliers
########################################
# Vslope = Vslope * 1
# Vint = Vint * 1
# Kslope = Kslope * 1
# Kint = Kint * 1
# CUE = CUE * 1
# Tau_MULT = 1
# desorb_MULT = 1
# fPHYS_MULT = 1
# VMAX_MULT = 1
# KM_MULT = 1
