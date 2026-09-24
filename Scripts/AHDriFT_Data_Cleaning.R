# ################################################################################
# Title: AHDriFT Data Cleaning
# Wildlife Insights Data Processing workflow, Script _ of _
# Author: Will Harrod
# Date Created: 2026-09-22
################################################################################

################################################################################
# 1) Prep ######################################################################
################################################################################

# Add packages
library(tidyverse)
library(fs)

# Clear environments
rm(list = ls())

# Set Working Directory
wd <- "/home/will/NCSU/R_Code/Wildlife_Insights_Data_Management/Data"

# View file names
basename(dir_ls(wd))

# Add the Raw data
seqs_raw <-  read_csv(path(wd, "sequences.csv"))
dpys_raw <- read_csv(path(wd, "deployments.csv"))

# View the data
glimpse(seqs_raw)
glimpse(dpys_raw)

################################################################################
# 2) Clean the sequence data ###################################################
################################################################################

# 2.1) CLeaning ----------------------------------------------------------------

# Remove unnecessary columns
seqs <- seqs_raw |> 
  # Select columns
  select(deployment_id, class, order, family, genus, species, common_name, start_time, end_time, group_size) |> 
  # Rename
  rename(
    Deployment = deployment_id,
    Class = class,
    Order = order,
    Family = family,
    Genus = genus,
    Species = species,
    Common.Name = common_name,
    Time.Start = start_time,
    Time.End = end_time,
    Group.Size = group_size
  ) |> 
  # Remove invertebrates
  filter(Class %in% c("Mammalia", "Aves", "Amphibia","Reptilia")) |> 
  # Add site information
  mutate(
    Plot.ID = str_sub(Deployment, start = 1, end = 4),
    Hut.ID = str_sub(Deployment, start = 11, end = 14),
    Sensor.Type = str_sub(Deployment, start = 6, end = 7),
    Deployment.ID = str_sub(Deployment, start = 21, end = 22),
    Settings = str_sub(Deployment, start = 24, end = 25)
  ) |> 
  # Select only one type of images
  filter(Sensor.Type %in% c("AD") & Settings %in% c("MC")) |> 
  select(-Sensor.Type, - Deployment)

# View
glimpse(seqs)

# 2.2) Sequence Summaries ------------------------------------------------------

# List of species 
seqs |> 
  count(Common.Name, Genus, Species) |> 
  arrange(-n) |> 
  print(n = Inf)

# 2.3) Remove species ----------------------------------------------------------

# Lump some rare etections by genus 
seqs_lumped <- seqs |> 
  mutate(Common.Name = case_when(
    Common.Name == "Broadhead Skink" ~ "Plestiodon Species",
    Common.Name == "Golden Mouse" ~ "Peromyscus or Ochrotomys Species",
    Common.Name == "Peromyscus Species" ~ "Peromyscus or Ochrotomys Species",
    Genus == "Dryophytes" ~ "Tree Frog Species",
    TRUE ~ Common.Name
    ))


# Cutoff for the fewest number of detections to include
sp_cutoff <-  8

# List of species with few detections
rare_sp <- seqs_lumped |> 
  count(Common.Name) |> 
  filter(n < sp_cutoff) |> 
  pull(Common.Name)
rare_sp
  
# Filter out those rare species
seqs_common <- seqs_lumped |> 
  filter(!Common.Name %in% rare_sp)

# List of all remaining species
seqs_common |> 
  count(Common.Name) |> 
  arrange(-n) |> 
  print(n = Inf)

# List of species to remove
drop_sp <- c("Carolina Wren", "Woodrat or Rat or Mouse Species")

# Filter those out
seqs_model <- seqs_common |> 
  filter(!Common.Name %in% drop_sp)

# View emaining species
seqs_model |> 
  count(Common.Name) |> 
  arrange(-n) |> 
  print(n = Inf)
glimpse(seqs_model)
  
################################################################################
# 3) Clean the deployment ######################################################
################################################################################

# 3.1) Initial cleaning --------------------------------------------------------

# View the deployments 
glimpse(dpys_raw)

# clean
dpys <- dpys_raw |> 
  # Rename
  rename(
    Deployment = deployment_id,
    Start.Date = start_date,
    End.Date = end_date,
    Camera.Model = camera_name,
    Notes = remarks,
    x = longitude,
    y = latitude
  ) |> 
  # Break the deployment ID into useful parts
  mutate(
    Plot.Type = str_sub(Deployment, start = 1, end = 2),     
    Plot.ID = str_sub(Deployment, start = 1, end = 4),
    Hut.ID = str_sub(Deployment, start = 11, end = 13),
    Sensor.Type = str_sub(Deployment, start = 6, end = 7),
    Deployment.ID = str_sub(Deployment, start = 21, end = 22),
    Settings = str_sub(Deployment, start = 24, end = 25)
  ) |> 
  mutate(
    Plot.Name = case_when(
    Plot.Type == "IF" ~ "Interior Forest",
    Plot.Type == "RE" ~ "Reference Edge",
    Plot.Type == "TO" ~ "Turbine Opening",
    Plot.Type == "TE" ~ "Turbine Edge"
  )) |> 
  # Select only one type of images
  filter(Sensor.Type %in% c("AD") & Settings %in% c("MC")) |> 
  # select useful columns
  select(
    Plot.ID, 
    Plot.Name,
    Hut.ID,
    Plot.Type, 
    Start.Date,
    End.Date, 
    Deployment.ID, 
    Camera.Model,
    Settings,
    Deployment,
    Notes, 
    x, 
    y
    ) 

# View after cleaning
glimpse(dpys)

# 3.2 Clean dates --------------------------------------------------------------

# View all the start dates
dpys |> distinct(Start.Date) |> arrange(Start.Date) |> print(n = Inf)
# View all the end dates
dpys |> distinct(End.Date) |> arrange(End.Date) |> print(n = Inf)
# Wow! No errors

# 3.3) Clean camera model info -------------------------------------------------

# View all of camera models
dpys |> 
  count(Camera.Model) |> 
  arrange(-n) |> 
  print(n = Inf)

# View the NAs 
dpys |> filter(is.na(Camera.Model) | Camera.Model == "?") |> select(Deployment, Camera.Model) |> print(n = Inf)

# List of browning camera models
hp5 <- c("Browning Recon Force HP5", "Browning Recon Force Elite")
dfndr <- c("Browning Defender Pro Scout Max Extreme HD")
reconyx <- c("Reconyx HC500 Hyperfire", "Reconyx Hyperfire PC500", "Reconyx Hyperfire PC800", "Reconyx Hyperfire HC500", "Reconyx PC900")

# Clean the camera models
dpys_cammod <- dpys |> 
  mutate(
    Camera.Model = case_when(
      Camera.Model %in% hp5 ~ "Browning HP5",
      Camera.Model %in% dfndr ~ "Browning Defender",
      Camera.Model %in% reconyx ~ "Reconyx Hyperfire 1",
      TRUE ~ Camera.Model
    ))

# View all of camera models again
dpys_cammod |> 
  count(Camera.Model) |> 
  arrange(-n) |> 
  print(n = Inf)
