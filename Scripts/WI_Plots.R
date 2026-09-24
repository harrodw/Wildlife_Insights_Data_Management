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

# Add the data

################################################################################
# 2) Plots ################################################################
################################################################################

# How many plot types
unique(seqs$Plot.Type)
unique(seqs$Plot.Name)

# List of species
seqs |> distinct(Common.Name) |> arrange(Common.Name) |>  print(n = Inf)

# Palette
wind_pal <- c(
  "Interior Forest" = "forestgreen", 
  "Reference Edge" = "goldenrod3", 
  "Turbine Opening" = "darkorchid4",
  "Turbine Edge" = "violetred4"
)

# Species to ignore
ignore_sp <- c(
  "Woodrat or Rat or Mouse Species",
  "Bird",
  "Rodent",
  "Eulipotyphla Order",
  "Eulipotyphla Order",
  "Frogs",
  "Mud Snake",
  "Small Mammal",
  "Rabbit and Hare Family"
)

# 2.1) Species by Plot Type ----------------------------------------------------
seqs |> 
  group_by(Common.Name, Plot.Name) |>
  reframe(Common.Name, Plot.Name, Count = n()) |> 
  distinct() |> 
  filter(!Common.Name %in% ignore_sp) |> 
  mutate(Common.Name = fct_reorder(Common.Name, Count, .fun = sum)) |>
  filter(Count > 5) |> 
  ggplot(aes(x = Common.Name, y = Count, fill = Plot.Name)) +
  geom_col() +
  coord_flip() + 
  scale_fill_manual(values = wind_pal) +
  scale_color_manual(values = wind_pal) +
  scale_y_continuous() +
  theme_classic() +
  labs(title = "Number of AHDriFT Detections") +
  labs(y = "Number of Detections") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_blank()) 

# 2.2) Proportion of sites occupied ---------------------------------------------

# number of plots by type.convert
n_plots <- seqs |> 
  distinct(Plot.Name, Plot.ID) |> 
  count(Plot.Name) |> 
  rename(n.Plots = n)
n_plots 

# Make the plot
seqs |> 
  filter(!Common.Name %in% ignore_sp) |> 
  group_by(Common.Name, Plot.ID) |> 
  reframe(Common.Name, Plot.ID, Plot.Type, Plot.Name, n.Detections = n()) |> 
  distinct() |> 
  mutate(Present = case_when(n.Detections > 0 ~ 1, TRUE ~ 0)) |> 
  group_by(Common.Name, Plot.Type) |> 
  reframe(Common.Name, Plot.Type, Plot.Name, n.Sites.Present = sum(Present)) |> 
  distinct() |> 
  left_join(n_plots, by = "Plot.Name") |> 
  mutate(Prop.Present = n.Sites.Present / n.Plots) |> 
  # Remove rare species
  group_by(Common.Name) |> 
  mutate(Total.Prop = sum(Prop.Present)) |> 
  ungroup() |> 
  filter(Total.Prop >= 0.12) |> 
  ggplot(aes(x = Plot.Type, y = Prop.Present, fill = Plot.Name)) +
  geom_col() +
  scale_fill_manual(values = wind_pal) +
  scale_color_manual(values = wind_pal) +
  theme_classic() +
  labs(title = "Proportion of AHDriFT arrays with at least one detection") +
  facet_wrap(~Common.Name) +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 8),
        legend.title = element_blank()) 
