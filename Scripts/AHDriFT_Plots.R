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
library(suncalc)

# Clear environments
rm(list = ls())

# Set Working Directory
wd <- "/home/will/NCSU/R_Code/Wildlife_Insights_Data_Management/Data"

# Directory for figures
fig_dir <- "/home/will/NCSU/Figures"

# View file names
basename(dir_ls(wd))

# Add the data
ahdrift_dat_raw <- read_csv(path(wd, "ahdrift_data_cleaned.csv"))
# View
glimpse(ahdrift_dat_raw)

# Clean the data 
ahdrift_dat <- ahdrift_dat_raw |> 
  mutate(Plot.Type = str_sub(Plot.ID, start = 1, end = 2)) |> 
  mutate(Plot.Type.Name = case_when(
    Plot.Type == "IF" ~ "Interior Forest",
    Plot.Type == "RE" ~ "Reference Edge",
    Plot.Type == "TE" ~ "Turbine Edge",
    Plot.Type == "TO" ~ "Turbine Opening" 
  ))

# View
glimpse(ahdrift_dat)
ahdrift_dat |> count(Common.Name, Plot.Type.Name)

# List of species 
sp_ls <- ahdrift_dat |> 
  filter(!Common.Name %in% c("No Detections", "Human", "Human-Camera Trapper")) |> 
  distinct(Common.Name) |> 
  arrange(Common.Name) |> 
  pull(Common.Name)
# View
sp_ls

################################################################################
# 2) Plots ################################################################
################################################################################

# Palette
wind_pal <- c(
  "Interior Forest" = "forestgreen", 
  "Reference Edge" = "goldenrod3", 
  "Turbine Edge" = "violetred4",
  "Turbine Opening" = "darkorchid4"
)

# 2.1) Species by Plot Type Scaled for effort----------------------------------------------------

# Caculate effort
ahdrift_efft <-  ahdrift_dat |>
  distinct(Plot.ID, Plot.Type.Name, Date) |> 
  count(Plot.Type.Name) |> 
  rename(Effort = n)
# View
glimpse(ahdrift_efft)

# Prep the data 
ahdrift_dat_efft <- ahdrift_dat |> 
  filter(Common.Name != "No Detections") |> 
  count(Common.Name, Plot.Type.Name) |> 
  rename(Count = n) |> 
  distinct() |> 
  left_join(ahdrift_efft,  by = "Plot.Type.Name") |> 
  mutate(Effort.Count = Count / Effort) |> 
  mutate(Common.Name = fct_reorder(Common.Name, Effort.Count, .fun = sum)) 
# View
glimpse(ahdrift_dat_efft)
ahdrift_dat_efft |> select(Common.Name, Plot.Type.Name, Count) |> print(n = Inf)

# Make the plot
ahdrift_dat_efft_fig <- ahdrift_dat_efft |> 
  ggplot(aes(x = Common.Name, y = Effort.Count, fill = Plot.Type.Name)) +
  geom_col() +
  coord_flip() + 
  scale_fill_manual(values = wind_pal) +
  scale_color_manual(values = wind_pal) +
  scale_y_continuous() +
  theme_classic() +
  labs(y = "Number of Detections per Trap Night") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_blank()) 
# View the plot
ahdrift_dat_efft_fig

# Save the plot as a png
ggsave(plot = ahdrift_dat_efft_fig,
       path(fig_dir, "ahdrift_effort_raw_count.png"),
       width = 200,
       height = 175,
       units = "mm",
       dpi = 300)

# 2.2) Activity Curves by Species ---------------------------------------------

# View the data again
glimpse(ahdrift_dat)

# View the species with only a few detection
sp_count <- ahdrift_dat |> 
  filter(Common.Name %in% sp_ls) |> 
  count(Common.Name, Class) |> 
  arrange(-n)
sp_count |> print(n = Inf)

# Isolate the species with few detections
min_dct <- 24
common_sp <- sp_count |>
  filter(n >= min_dct) |>
  pull(Common.Name)

# Define the classes
mams <- sp_count |> filter(Class == "Mammalia") |> pull(Common.Name)
reps <- sp_count |> filter(Class == "Reptilia") |> pull(Common.Name)
amphs <- sp_count |> filter(Class == "Amphibia") |> pull(Common.Name)

# Constrain to species over 24 detections
ahdrift_dat_common <- ahdrift_dat |> 
  filter(Common.Name %in% common_sp) 

# List of plot types
plt_typ <- ahdrift_dat |> 
  distinct(Plot.Type.Name) |> 
  pull(Plot.Type.Name)

# List of dates
dates <- ahdrift_dat |> 
  distinct(Date) |> 
  arrange(Date) |> 
  pull(Date)

# Center of the study area 
sr_center <- ahdrift_dat |> 
  distinct(x, y) |> 
  mutate(x = mean(x), y = mean(y)) |> 
  distinct()  

# Calculate Sunrises and sunsets
sunr_tm <- getSunlightTimes(
  date = dates,
  lat = sr_center$y,
  lon = sr_center$x,
  keep = c("sunrise", "sunset")
) |> 
  mutate(
    # Convert timestamps from UTC to local time
    sunrise = with_tz(sunrise, tzone = "America/New_York"),
    sunset  = with_tz(sunset, tzone = "America/New_York"),
    # Recalculate minutes past midnight in local time
    Rise.Time = hour(sunrise) * 60 + minute(sunrise),
    Set.Time = hour(sunset) * 60 + minute(sunset),
  ) 

# View
glimpse(sunr_tm)

# Take summary stats (In hours)
rise_mode <- mean(sunr_tm$Rise.Time) / 60
rise_min <- quantile(sunr_tm$Rise.Time, 0.025) / 60
rise_max <- quantile(sunr_tm$Rise.Time, 0.975) / 60
set_mode <- mean(sunr_tm$Set.Time) / 60
set_min <- quantile(sunr_tm$Set.Time, 0.025) / 60
set_max <- quantile(sunr_tm$Set.Time, 0.975) / 60

# Prep the Detection Data 
acvty_dat <- ahdrift_dat_common |> 
  mutate(Hour = hour(Time.Start)) |> 
  count(Common.Name, 
        # Plot.Type.Name, 
        Hour, name = "Hour.Detections") |>
  complete(
    Common.Name = common_sp, 
    # Plot.Type.Name = plt_typ,
    Hour = 0:23, 
    fill = list(Hour.Detections = 0)
  ) |> 
  mutate(Hour.Detections = replace_na(Hour.Detections, 0)) |> 
  group_by(
    # Plot.Type.Name,
    Common.Name
  ) |> 
  mutate(
    Total.Detections = sum(Hour.Detections),
    Prop.Detections = if_else(Total.Detections > 0, Hour.Detections / Total.Detections, 0)
  ) |> 
  ungroup() |> 
  arrange(
    Common.Name, 
    # Plot.Type.Name, 
    Hour
  ) |> 
  # Rename Classes
  mutate(
    Class = case_when(
      Common.Name %in% mams ~ "Mammals",
      Common.Name %in% reps ~ "Reptiles",
      Common.Name %in% amphs ~ "Amphibians"
    ))

# View again
glimpse(acvty_dat)
count(acvty_dat, Common.Name, Class)

# Make the Plot
acvty_plt <- acvty_dat |> 
  mutate(Common.Name = case_when(
    Common.Name == "Peromyscus or Ochrotomys Species" ~ "Cricetidae Mice",
    Common.Name == "Southern Short-tailed Shrew" ~ "Short-tailed Shrew",
    Common.Name == "North American Least Shrew" ~ "Least Shrew",
    TRUE ~ Common.Name
    ),
    Class = factor(Class, levels = c("Mammals", "Reptiles", "Amphibians"))
    ) |> 
  arrange(Class, Common.Name) |> 
  mutate(Common.Name = fct_inorder(Common.Name)) |>
  ggplot(aes(x = Hour, y = Prop.Detections)) +
  # Detection Data
  geom_smooth(aes(color = Class), span = 0.35, se = FALSE, linewidth = 1) +
  # geom_point(aes(color = "Relitive Animal Activity"),size = 1.5) +
  # Sunrise Lines
  geom_vline(aes(xintercept = rise_mode, color = "Sunrise"), 
             linewidth = 0.8, alpha = 0.75, linetype = "solid") +
  geom_vline(aes(xintercept = rise_min, color = "Sunrise"), 
             linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  geom_vline(aes(xintercept = rise_max, color = "Sunrise"), 
             linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  # Sunset Lines
  geom_vline(aes(xintercept = set_mode, color = "Sunset"), 
             linewidth = 0.8, alpha = 0.75, linetype = "solid") +
  geom_vline(aes(xintercept = set_min, color = "Sunset"), 
             linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  geom_vline(aes(xintercept = set_max, color = "Sunset"), 
             linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  # Legend
  scale_color_manual(
    name = "", 
    values = c(
      "Sunrise" = "goldenrod", "Sunset" = "orangered4", 
               "Mammals" = "tan4", "Reptiles" = "springgreen4", "Amphibians" = "blue2" 
  )) +
  # Other Formating
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 23, by = 4)) +
  labs(x = "Hour", y = "Proportion of Detection Events") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom",          
    legend.box = "horizontal"            
  ) +
  facet_wrap(~Common.Name)

# View the Plot
acvty_plt

# Save the plot as a png
ggsave(plot = acvty_plt,
       path(fig_dir, "ahdrift_activity_curves.png"),
       width = 200,
       height = 175,
       units = "mm",
       dpi = 300)

