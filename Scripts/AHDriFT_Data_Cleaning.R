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

# List of species 
seqs |> 
  count(Common.Name, Genus, Species) |> 
  arrange(-n) |> 
  print(n = Inf)

# 2.2) Remove species ----------------------------------------------------------

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
sp_cutoff <-  3

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
drop_sp <- c(
  "Carolina Wren", 
  "Woodrat or Rat or Mouse Species", 
  "Frogs", 
  "Small Mammal", 
  "Eulipotyphla Order", 
  "American Black Bear"
             )

# Filter those out
seqs_sp_lst <- seqs_common |> 
  filter(!Common.Name %in% drop_sp)

# View remaining species
seqs_sp_lst  |> 
  count(Common.Name) |> 
  arrange(-n) |> 
  print(n = Inf)
glimpse(seqs_sp_lst)

# 2.3) Combine detections at the same time -------------------------------------

# Prep the data to combine detections within 60 minutes
seq_model <- seqs_sp_lst %>%
  #Sort chronologically by location, species, and timestamp
  arrange(Plot.ID, Class, Genus, Species, Common.Name, Time.Start) %>%
  # Group by location and taxonomy/species metadata
  group_by(Plot.ID, Class, Genus, Species, Common.Name) %>%
  # Identify continuous events separated by gaps > 60 minutes
  mutate(
    # Time difference from previous detection in minutes
    time_gap = as.numeric(difftime(Time.Start, lag(Time.Start), units = "mins")),
    # Increment event ID whenever gap > 60 mins (or on the first row where lag is NA)
    event_id = cumsum(coalesce(time_gap > 60, TRUE))
  ) %>%
  # Collapse each detection event into a single record
  group_by(Common.Name, Plot.ID, Class, Genus, Species, event_id) %>%
  summarise(
    Time.Start = min(Time.Start),
    Time.End = max(Time.End),
    n.Seq = n(),
    Max.Group.Size = max(Group.Size, na.rm = TRUE),
    # Hut.ID = Hut.ID,
    .groups= "drop"
  ) |> 
  # Remove columns that are no longer needed
  select(-event_id) |> 
  mutate(Date = as.Date(Time.Start),
         Plot.Date = paste(Plot.ID, Date, sep = "-")) |> 
  relocate(Date, .before = Time.Start) 
# View
glimpse(seq_model)

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
  # Break the deployment ID into useful parts
  mutate(
    Plot.Type = str_sub(Deployment, start = 1, end = 2),     
    Plot.ID = str_sub(Deployment, start = 1, end = 4),
    Hut.ID = str_sub(Deployment, start = 11, end = 14),
    Sensor.Type = str_sub(Deployment, start = 6, end = 7),
    Deployment.ID = str_sub(Deployment, start = 21, end = 22),
    Settings = str_sub(Deployment, start = 24, end = 25)
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
dpys |>  count(Hut.ID)
dpys |> count(Plot.ID) |>  print(n = Inf)
dpys |>  filter(!Hut.ID %in% c("CamA", "CamB")) |> select(Plot.ID, Deployment.ID, Deployment, Hut.ID)

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

# 3.3) Make detection data longer ----------------------------------------------

# Find each survey occation for the cameras
surveys <- dpys_cammod |> 
  distinct(Plot.ID, Hut.ID, Start.Date, End.Date) |> 
  arrange(Plot.ID, Hut.ID, Start.Date, End.Date)
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
  # Keep core columns and remove duplicate dates
  select(Plot.ID, Hut.ID, Date) %>%
  distinct()

# View
glimpse(survey_dates)

# Pivot wider to combine huts and estimate effort
survey_dates_efft <- survey_dates |> 
  group_by(Plot.ID, Date) |> 
  summarise(
    Huts.Active = n(),
    .groups = "drop"
  ) |> 
  mutate(Plot.Date = paste(Plot.ID, Date, sep = "-"))
# View
glimpse(survey_dates_efft)
survey_dates_efft |> count(Huts.Active)
survey_dates_efft |> filter(Huts.Active > 2)

# 3.4) Combine deployment and detection data -----------------------------------

# View
glimpse(seq_model)
glimpse(survey_dates_efft)

# Find the dates where no species were found
blank_dates <- survey_dates_efft |> 
  filter(!Plot.Date %in% seq_model$Plot.Date)
glimpse(blank_dates)

# Combine the data and add the missing dates
seq_dpy <- seq_model |> 
  left_join(survey_dates_efft, by = c("Plot.ID", "Date", "Plot.Date")) |> 
  bind_rows(blank_dates) |> 
  select(-Plot.Date) |> 
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

# How many of the hut-nights had detections?
seq_dpy |> count(Common.Name) |> print(n = Inf)

# View the detections that do not match an official deployment date
seq_dpy |> filter(is.na(Huts.Active) & Date > ymd("2026-01-01")) |> 
  select(Common.Name, Plot.ID, Date, n.Seq, Huts.Active) |> 
  print(n = Inf)

# Remove the detections with no date if they look like true post-deployment detections
ahdrift_dat <- seq_dpy |> filter(!is.na(Huts.Active))
# View one last time
glimpse(ahdrift_dat)

# Save
write_csv(ahdrift_dat, path(wd, "ahdrift_data_cleaned.csv"))
