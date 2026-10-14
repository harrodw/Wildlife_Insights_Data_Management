# ################################################################################
# Title: AHDriFT Data Cleaning
# Wildlife Insights Data Processing workflow, Script 1 of _
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
  # Pull out only the camera traps
  filter(str_detect(Deployment, "CT")) |> 
  # Add site information
  mutate(
    Plot.ID = str_sub(Deployment, start = 1, end = 4),
    Camera.ID = str_sub(Deployment, start = 6, end = 9),
    Deployment.ID = str_sub(Deployment, start = 16, end = 17)
  ) 

# View
glimpse(seqs)

# List of Plots
seqs |> 
  count(Plot.ID) |> 
  arrange(-n) |> 
  print(n = Inf)

# List of sensors by plot
seqs |> 
  distinct(Plot.ID, Camera.ID) |> 
  count(Plot.ID) |> 
  arrange(-n) |> 
  print(n = Inf)

# List of species 
seqs |> 
  count(Common.Name, Class, Genus, Species) |> 
  arrange(-n) |> 
  print(n = Inf)

# Now with just the species
seqs |> 
  count(Common.Name, Class, Genus, Species) |> 
  arrange(-n) |> 
  pull(Common.Name)

# 2.2) Remove species ----------------------------------------------------------

  
# Filter to only the things that I want
seqs_common <- seqs |>  
  filter(Common.Name %in%c(
    "White-tailed Deer",
    "Vehicle", 
    "Human-Camera Trapper",
    "Wild Turkey",
    "Domestic Dog",
    "Marsh Rabbit",
    "Northern Raccoon",
    "Coyote",
    "American Black Bear",
    "Eastern Gray Squirrel", 
    "Virginia Opossum", 
    "Human", 
    "Eastern Cottontail",
    "Bobcat",
    "North American River Otter"
  ))

# List of all remaining species
seqs_common |> 
  count(Common.Name) |> 
  arrange(-n) |> 
  print(n = Inf)

# 2.3) Combine detections at the same time -------------------------------------

# Prep the data to combine detections within 60 minutes
seq_model <- seqs_common |> 
  #Sort chronologically by location, species, and timestamp
  arrange(Plot.ID, Camera.ID, Class, Genus, Species, Common.Name, Time.Start) %>%
  # Group by location and taxonomy/species metadata
  group_by(Plot.ID, Common.Name) %>%
  # Identify continuous events separated by gaps > 60 minutes
  mutate(
    # Time difference from previous detection in minutes
    time_gap = as.numeric(difftime(Time.Start, lag(Time.Start), units = "mins")),
    # Increment event ID whenever gap > 60 mins (or on the first row where lag is NA)
    event_id = cumsum(coalesce(time_gap > 60, TRUE))
  ) %>%
  # Collapse each detection event into a single record
  group_by(Common.Name, Plot.ID, Camera.ID, Class, Genus, Species, event_id) %>%
  summarise(
    Time.Start = min(Time.Start),
    Time.End = max(Time.End),
    n.Seq = n(),
    Max.Group.Size = max(Group.Size, na.rm = TRUE),
    # Hut.ID = 
    .groups= "drop"
  ) |> 
  mutate(Date = as.Date(Time.Start),
         Plot.Date = paste(Plot.ID, Camera.ID, Date, sep = "-")) |> 
  relocate(Date, .before = Time.Start) |> 
  # Remove columns that are no longer needed
  select(-event_id, -Plot.ID, -Camera.ID, -Date)

# View
glimpse(seq_model)
seq_model |> 
  count(Common.Name) |> 
  arrange(-n) |> 
  print(n = Inf)

################################################################################
# 3) Clean the deployment data #################################################
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
  # Pull out only the camera traps
  filter(str_detect(Deployment, "CT")) |> 
  # Add site information
  mutate(
    Plot.ID = str_sub(Deployment, start = 1, end = 4),
    Plot.Type = str_sub(Deployment, start = 1, end = 2),
    Camera.ID = str_sub(Deployment, start = 6, end = 9),
    Deployment.ID = str_sub(Deployment, start = 16, end = 17)
  ) |> 
  mutate(
    Plot.Name = case_when(
    Plot.Type == "IF" ~ "Interior Forest",
    Plot.Type == "RE" ~ "Reference Edge",
    Plot.Type == "TO" ~ "Turbine Opening",
    Plot.Type == "TE" ~ "Turbine Edge"),
  Start.Date = date(Start.Date),
  End.Date = date(End.Date)
  ) |> 
  # select useful columns
  select(
    Plot.ID, 
    Camera.ID,
    Plot.Type,
    Plot.Name,
    Plot.Type, 
    Start.Date,
    End.Date, 
    Deployment.ID, 
    Camera.ID,
    Camera.Model,
    Deployment,
    Notes, 
    x, 
    y
    ) 

# View after cleaning
glimpse(dpys)
dpys |> count(Plot.ID) |>  print(n = Inf)
dpys |> count(Plot.Name)

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
hp5 <- c("Browning Recon Force HP5", "MM251")
dfndr <- c("Browning Defender Pro Scout Max Extreme HD", "Browning Strike Force Pro", "Browning Recon Force Elite")
reconyx <- c("Reconyx HC500 Hyperfire", "Reconyx Hyperfire HC500", "Reconyx Hyperfire 4K Ultra HD")

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
glimpse(dpys_cammod)
dpys_cammod |> 
  count(Camera.Model) |> 
  arrange(-n) |> 
  print(n = Inf)

# 3.3) Make detection data longer ----------------------------------------------

# Find each survey occation for the cameras
surveys <- dpys_cammod |> 
  distinct(Plot.ID, Plot.Type, Plot.Name, Camera.ID, Camera.Model, Deployment.ID, Start.Date, End.Date, x, y) |> 
  arrange(Plot.ID, Camera.ID, Camera.Model, Start.Date, End.Date)
# View
glimpse(surveys)

# Go through the deployments to make a row for every date
survey_dates <- surveys %>%
  mutate(
    Start.Date = ymd(Start.Date),
    End.Date   = ymd(End.Date)
  ) %>%
  # Expand date sequences for each deployment row
  rowwise() %>%
  mutate(Date = list(seq(Start.Date, End.Date, by = "day"))) %>%
  unnest(Date) %>%
  mutate(Plot.Date = paste(Plot.ID, Camera.ID, Date, sep = "-")) |> 
  # Keep core columns and remove duplicate dates
  select(Plot.Date, Plot.ID, Camera.ID, Plot.Type, Plot.Name, Deployment.ID, Date,  Camera.Model, x, y) %>%
  distinct()

# View
glimpse(survey_dates)
survey_dates |> count(Plot.Date) |> filter(n > 1)

# 3.4) Combine deployment and detection data -----------------------------------

# View
glimpse(seq_model)
glimpse(survey_dates)

# Find the dates where no species were found
blank_dates <- survey_dates |> 
  filter(!Plot.Date %in% seq_model$Plot.Date)
glimpse(blank_dates)

# Combine the data and add the missing dates
seq_dpy <- seq_model |> 
  left_join(survey_dates, by = c("Plot.Date")) |> 
  bind_rows(blank_dates) |> 
  # select(-Plot.Date) |> 
  arrange(Plot.ID, Date, Time.Start, Common.Name) |> 
  mutate(Common.Name = replace_na(Common.Name, "No Detections")) |> 
  mutate(across(
    .cols = c("Class", "Genus", "Species"),
    .fns = ~ifelse(Common.Name == "No Detections",
                   yes = replace_na(., "No Detections"),
                   no = .
                   ))) |> 
  mutate(across(.cols = c("n.Seq", "Max.Group.Size"), .fns = ~ replace_na(0)))
# View
glimpse(seq_dpy)

# How many of the trap days had detections?
seq_dpy |> count(Common.Name) |> arrange(-n) |> print(n = Inf)

# See which detection evenets are missing information
# These should be fixed manually in Wildlife Insights
na_seq <- seq_dpy |> 
  filter(is.na(Plot.ID)) |> 
  arrange(Plot.Date) |> 
  select(Common.Name, Plot.Date, Plot.ID, Plot.Type, Plot.Name, Camera.ID, Date) 
# View
na_seq |> print(n = Inf)

# Remove the ones that could not be fixed manually 
seq_true <- seq_dpy |> 
  filter(!is.na(Plot.ID))

# View one last time
glimpse(seq_true)

# Save
write_csv(seq_true, path(wd, "camera_trap_data_cleaned.csv"))
