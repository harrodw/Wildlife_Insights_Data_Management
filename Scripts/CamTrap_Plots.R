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
camtrap_dat_raw <- read_csv(path(wd, "camera_trap_data_cleaned.csv"))
# View
glimpse(camtrap_dat_raw)

# Clean the data 
camtrap_dat <- camtrap_dat_raw |> 
  mutate(Plot.Type = str_sub(Plot.ID, start = 1, end = 2)) 

# View
glimpse(camtrap_dat)
camtrap_dat |> count(Common.Name, Plot.Type.Name)

# List of species 
sp_ls <- camtrap_dat |> 
  filter(!Common.Name %in% c(
    "No Detections",
    "Domestic Dog",
    "Human",
    "Vehicle",
    "Human-Camera Trapper"
                             )) |> 
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
camtrap_efft <-  camtrap_dat |>
  distinct(Plot.ID, Plot.Type.Name, Date) |> 
  count(Plot.Type.Name) |> 
  rename(Effort = n)
# View
camtrap_efft |> print(n = Inf)

# Prep the data 
camtrap_dat_efft <- camtrap_dat |> 
  filter(Common.Name %in% sp_ls) |> 
  count(Common.Name, Plot.Type.Name) |> 
  rename(Count = n) |> 
  distinct() |> 
  left_join(camtrap_efft,  by = "Plot.Type.Name") |> 
  mutate(Effort.Count = 100*Count / Effort) |> 
  mutate(Common.Name = fct_reorder(Common.Name, Effort.Count, .fun = sum)) 
# View
glimpse(camtrap_dat_efft)
camtrap_dat_efft |> select(Common.Name, Plot.Name, Count) |> print(n = Inf)

# Make the plot
camtrap_dat_efft_fig <- camtrap_dat_efft |> 
  ggplot(aes(x = Common.Name, y = Effort.Count, fill = Plot.Name)) +
  geom_col() +
  coord_flip() + 
  scale_fill_manual(values = wind_pal) +
  scale_color_manual(values = wind_pal) +
  scale_y_continuous() +
  theme_classic() +
  labs(y = "Number of Detections per 100 Trap Nights") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_blank()) 
# View the plot
camtrap_dat_efft_fig

# Save the plot as a png
ggsave(plot = camtrap_dat_efft_fig,
       path(fig_dir, "camtrap_effort_raw_count.png"),
       width = 200,
       height = 175,
       units = "mm",
       dpi = 300)

# 2.2) Activity Curves by Species ---------------------------------------------

# View the data again
glimpse(camtrap_dat)

# View the species with only a few detection
sp_count <- camtrap_dat |> 
  filter(Common.Name %in% sp_ls) |> 
  count(Common.Name) |> 
  arrange(-n)
sp_count

# Isolate the species with few detections
min_dct <- 24
common_sp <- sp_count |> 
  filter(n >= min_dct) |> 
  pull(Common.Name)
  
# Constrain to species over 24 detections
camtrap_dat_common <- camtrap_dat |> 
  filter(Common.Name %in% common_sp) 

# List of plot types
plt_typ <- camtrap_dat |> 
  distinct(Plot.Name) |> 
  pull(Plot.Name)

# List of dates
dates <- camtrap_dat |> 
  distinct(Date) |> 
  arrange(Date) |> 
  pull(Date)

# Center of the study area 
sr_center <- camtrap_dat |> 
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
sunr_tm |> 
  ggplot(aes(x = Rise.Time)) +
  geom_histogram()

# Take summary stats (In hours)
rise_mode <- mean(sunr_tm$Rise.Time) / 60
rise_min <- quantile(sunr_tm$Rise.Time, 0.025) / 60
rise_max <- quantile(sunr_tm$Rise.Time, 0.975) / 60
set_mode <- mean(sunr_tm$Set.Time) / 60
set_min <- quantile(sunr_tm$Set.Time, 0.025) / 60
set_max <- quantile(sunr_tm$Set.Time, 0.975) / 60

# Prep the Detection Data 
acvty_dat <- camtrap_dat_common |> 
  mutate(Hour = hour(Time.Start)) |> 
  count(Common.Name, 
        # Plot.Name, 
        Hour, name = "Hour.Detections") |>
  complete(
    Common.Name = common_sp, 
    # Plot.Name = plt_typ,
    Hour = 0:23, 
    fill = list(Hour.Detections = 0)
    ) |> 
  mutate(Hour.Detections = replace_na(Hour.Detections, 0)) |> 
  group_by(
    # Plot.Name,
    Common.Name
    ) |> 
  mutate(
    Total.Detections = sum(Hour.Detections),
    Prop.Detections = if_else(Total.Detections > 0, Hour.Detections / Total.Detections, 0)
    ) |> 
  ungroup() |> 
  arrange(
    Common.Name, 
    # Plot.Name, 
    Hour
    )


# View again
glimpse(acvty_dat)
count(acvty_dat, Common.Name)

# Make the Plot
acvty_plt <- acvty_dat |> 
  ggplot(aes(x = Hour, y = Prop.Detections)) +
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
  # Detection Data
  geom_smooth(aes(color = "Relitive Animal Activity"), span = 0.35, se = FALSE, linewidth = 1) +
  # geom_point(aes(color = "Relitive Animal Activity"),size = 1.5) +
  # Legend
  scale_color_manual(
    name = "", 
    values = c("Sunrise" = "goldenrod", "Sunset" = "orangered4", "Relitive Animal Activity" = "black")
  ) +
  # Other Formating
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 23, by = 4)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 1)) +
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
       path(fig_dir, "mammal_activity_curves.png"),
       width = 200,
       height = 175,
       units = "mm",
       dpi = 300)

