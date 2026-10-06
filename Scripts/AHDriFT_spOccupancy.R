# Title: AHDriFT spOccupancy
# Wildlife Insights Data Processing workflow, Script 2 of _
# Author: Will Harrod
# Date Created: 2026-09-29
################################################################################

################################################################################
# 1) Prep ######################################################################
################################################################################

# Add packages
library(tidyverse)
library(fs)
library(spOccupancy)
library(MCMCvis)

# Clear environments
rm(list = ls())

# Set Working Directory
wd <- "/home/will/NCSU/R_Code/Wildlife_Insights_Data_Management/Data"

# View file names
basename(dir_ls(wd))

# Read the AHDriFT data back in
ahdrift_dat_raw <- read_csv(path(wd, "ahdrift_data_cleaned.csv"))
# View
glimpse(ahdrift_dat_raw)

# Clean the data for spOccupancy
ahdrift_dat_spOcc <-  ahdrift_dat_raw |> 
  # Clean the Spcies Names
  mutate(
    Common.Name = str_replace_all(Common.Name, " ", "."),
    Common.Name = str_replace_all(Common.Name, "-", "."),
    Common.Name = str_replace_all(Common.Name, "'", ".")
  ) |> 
  # Convert categorical variable to factors and scale numeric covs
  mutate(Plot.Type = str_sub(Plot.ID, start = 1, end = 2)) |> 
  arrange(Plot.Type, Plot.ID, Date) |> 
  mutate(Plot.ID.fct = as.numeric(factor(Plot.ID)),
         Plot.Type.fct = as.numeric(factor(Plot.Type)),
         Date.num = parse_number(as.character(Date - ymd("2026-01-01"))),
         Date.num = Date.num - min(Date.num) + 1,
         Date.scl = scale(Date.num)[,1]) |> 
  # Select Only necessary columns 
  select(Common.Name, Plot.ID, Plot.ID.fct, Plot.Type, Plot.Type.fct, 
         Date, Date.num, Date.scl, Huts.Active) 


################################################################################
# 2) Prepare the data for spOccupancy ##########################################
################################################################################
glimpse(ahdrift_dat_spOcc)

# 2.1) Prepare the detection matrix --------------------------------------------

# List of species 
sp_ls <-  ahdrift_dat_spOcc |>
  filter(Common.Name != "No.Detections") |> 
  distinct(Common.Name) |> 
  arrange(Common.Name) |> 
  pull(Common.Name)
sp_ls

# List of plots 
plt_ls <-  ahdrift_dat_spOcc |> distinct(Plot.ID.fct) |> 
  arrange(Plot.ID.fct) |> 
  pull(Plot.ID.fct)
plt_ls

# List of sampling occasions 
vst_ls <-  ahdrift_dat_spOcc |> 
  distinct(Date.num) |> 
  arrange(Date.num) |> 
  pull(Date.num)
vst_ls

# Matrix dimensions
n_sp <-  length(sp_ls)
n_plt <- length(plt_ls)
n_vst <-  length(vst_ls)
c(n_sp, n_plt, n_vst)

# Create the blank detection maxtrix (as an array)
dct_mtx <- array(data = NA, dim = c(n_sp, n_plt, n_vst))
str(dct_mtx)
dct_mtx[1, ,]

# Create matrixes for detection covariates 
day_mtx <- matrix(data = NA, nrow = n_plt, ncol = n_vst) # Survey Dates
str(day_mtx)
eft_mtx <- matrix(data = NA, nrow = n_plt, ncol = n_vst) # Number of active buckets 
str(eft_mtx)

# List of visits by plot
plt_vst_ls <-  ahdrift_dat_spOcc |> 
  distinct(Plot.ID, Plot.ID.fct, Date, Date.num, Date.scl) |> 
  arrange(Plot.ID.fct, Date.num)
glimpse(plt_vst_ls)

# Fill in the Detection Matrix 
for(s in 1:n_sp){

  # Define a species 
  sp <- sp_ls[s]
  
  # Filter detections of that species
  sp_dct <- ahdrift_dat_spOcc |> 
    filter(Common.Name == sp) |> 
    select(Plot.ID.fct, Date.num) |> 
    mutate(Detected = 1) |> 
    distinct() |> 
    right_join(plt_vst_ls, by = c("Plot.ID.fct", "Date.num")) |> 
    mutate(Detected = replace_na(Detected, 0)) |> 
    arrange(Date, Plot.ID)
  
  # Loop over the plots to fill the matrix 
  for(j in 1:n_plt){
   
    # Filter detentions to that plot
    plt_dct <- sp_dct |> 
      filter(Plot.ID.fct == j) |>
      pull(Detected)
    
    # Pull out detection information
    plt_det_cov_dat <- ahdrift_dat_spOcc |> 
      filter(Plot.ID.fct == j) |> 
      select(Plot.ID.fct, Date.scl, Huts.Active) |> 
      distinct()
    
    # Number of surveys at that plot
    n_vst_plt <-  length(plt_dct)
    
    # Assign those values the appropriate part of the detection matrix or covariate matrix
    dct_mtx[s, j, 1:n_vst_plt] <- plt_dct
    day_mtx[j, 1:n_vst_plt] <- plt_det_cov_dat$Date.scl
    eft_mtx[j, 1:n_vst_plt] <- plt_det_cov_dat$Huts.Active
  }
  
  # Message
  message("Created detection matrix for ", sp, " Species ", s, " out of ", n_sp)
  
}

# View
str(dct_mtx)
dct_mtx[6, ,]
str(day_mtx)
day_mtx
str(eft_mtx)
eft_mtx

# Collapse detection histories for occupancy priors
occupied <- apply(dct_mtx, c(1, 2), max, na.rm = TRUE)
occupied

# 2.2) Prepare detection covariates --------------------------------------------

# Combine detection covs as a list
det_covs <- list(
  date = day_mtx,
  daily.effort = eft_mtx
)
det_covs

# 2.3) Prepare occupancy covariates --------------------------------------------

# View the full dataset
glimpse(ahdrift_dat_spOcc)

# Calculate how many survey occations there were per plot
tot_eft <- ahdrift_dat_spOcc |> 
  distinct(Plot.ID, Date.num) |> 
  count(Plot.ID)
#
tot_eft

# Filter the occupancy level covariates to be distinct at each plot
occ_covs <- ahdrift_dat_spOcc |> 
  distinct(Plot.ID, Plot.ID.fct, Plot.Type.fct) |> 
  arrange(Plot.ID) |>
  left_join(tot_eft, by = "Plot.ID") |> 
  # Change covariate names
  mutate(
    plot.type = Plot.Type.fct,
    total.effort = scale(n)[,1]
  ) |> 
  # Here is where you can change what goes in the occupancy covariates
  select(
    plot.type,
    total.effort
    ) |> 
  # Convert too matrix
  as.matrix() 
# View
occ_covs

################################################################################
# 3) Run spOccupancy ###########################################################
################################################################################

# 3.1) Priors, inits, and data bundle ------------------------------------------

# Bundle the data 
dat_lst <- list(
  y = dct_mtx,
  occ.covs = occ_covs,
  det.covs = det_covs
)
dat_lst

# Define inits 
inits <- list(
  alpha.comm = 0,
  beta.comm = 0,
  alpha = 0,
  beta = 0, 
  tau.sq.beta = 1, 
  tau.sq.alpha = 1,
  z = occupied
)
inits

# Define priors
priors <- list(
  beta.comm.normal = list(mean = 0, var = 2.72),
  alpha.comm.normal = list(mean = 0, var = 2.72), 
  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1)
)

# MCMC parameters
n_sample <- 20000
n_rprt <- n_sample/2
n_burn <- n_sample/2
n_thin <-  125
n_chains <- 3

# How many samples per chain?
((n_sample - n_burn) / n_thin)

# Model formulas 
occ_formu <- ~ factor(plot.type) + total.effort
det_formu <- ~ factor(daily.effort) + date + I(date^2)

# Directory for the model output
mcmc_dir <- "/home/will/NCSU/Model_Outputs"

# 3.2) Run the model -----------------------------------------------------------

# Run the model
ahdrifft_occ_mod1 <- msPGOcc(
  occ.formula = occ_formu,  # Occupancy Formula
  det.formula = det_formu,  # Detection formula
  data = dat_lst,           # Data
  inits = inits,            # Initial Values
  priors = priors,          # Priors
  verbose = TRUE,           # Display messages
  n.samples = n_sample,
  n.report = n_rprt,
  n.burn = n_burn,
  n.thin = n_thin,
  n.chains = n_chains,
  n.omp.threads = 4
)

# Make sure all is good
summary(ahdrifft_occ_mod1)

# Save the Model summary
saveRDS(ahdrifft_occ_mod1, path(mcmc_dir, "ahdrifft_occ_mod1.rds"))

################################################################################
# 4) Model Diagnostics #########################################################
################################################################################

# Directory for the model output
mcmc_dir <- "/home/will/NCSU/Model_Outputs"

# Directory for figures
fig_dir <- "/home/will/NCSU/Figures"

# Load the output back in
ahdrifft_occ_mod1 <- readRDS(path(mcmc_dir, "ahdrifft_occ_mod1.rds"))

# View MCMC summary
summary(ahdrifft_occ_mod1)

# Traceplots 
plot(ahdrifft_occ_mod1, 'beta', density = FALSE) 


# Are R-hat values good?
summary(ahdrifft_occ_mod1)

# View MCMC plot
MCMCplot(
  object = ahdrifft_occ_mod1$samples,
  # excl = c("fit_pa", "fit_pa_new", "fit_pd", "fit_pd_new"),
  guide_lines = TRUE,
  params = sobs_params
         )

