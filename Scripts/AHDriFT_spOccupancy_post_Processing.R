# ##############################################################################
# Title: Occupancy modeling using spOccupancy for AHDriFT photo data
# Author: Will Harrod
# Date Created: 2026-10-06
################################################################################

################################################################################
# 1) Prep ######################################################################
################################################################################

# 1.1) Add packages and data ---------------------------------------------------

# Clear environments
rm(list = ls())

# Load Packages
library(tidyverse)
library(spOccupancy)
library(fs)

# Set Seed
set.seed(27606)

# 1.2) File directories (Change these for your device) -------------------------

# Set Working Directory for csv's 
wd <- "/home/will/NCSU/R_Code/Wildlife_Insights_Data_Management/Data"

# Directory for the model output
mcmc_dir <- "/home/will/NCSU/Model_Outputs"

# Directory for figures
fig_dir <- "/home/will/NCSU/Figures"

# Start Date for the bird surveys
start_date <- ymd("2026-05-01")

# 1.3) Add the filtered BirdNET detections and other data ----------------------

# View file names
basename(dir_ls(wd))

# Read the AHDriFT data back in
ahdrift_dat_raw <- read_csv(path(wd, "ahdrift_data_cleaned.csv"))
# View
glimpse(ahdrift_dat_raw)

# Clean the data for spOccupancy
ahdrift_dat_spOcc <-  ahdrift_dat_raw |> 
  # Clean the Species Names
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
# View
glimpse(ahdrift_dat_spOcc)

# List of species 
sp_ls <- ahdrift_dat_spOcc |> 
  distinct(Common.Name) |> 
  arrange(Common.Name) |> 
  pull(Common.Name)
# View
sp_ls

# Load the spOccupancy output back in
ahdrift_occ_mod <- readRDS(path(mcmc_dir, "ahdrifft_occ_mod1.rds"))

# View MCMC summary
names(ahdrift_occ_mod)
summary(ahdrift_occ_mod)

################################################################################
# 2) Model Diagnostics #########################################################
################################################################################

# 2.1) Preliminary diagnostics -------------------------------------------------

# Traceplots 
# plot(ahdrift_occ_mod, 'beta', density = FALSE)
# plot(ahdrift_occ_mod, "alpha", density = FALSE)

# 2.2) Goodness of fit stats by plot -------------------------------------------
ahdrift_occ_out_plt <- ppcOcc(ahdrift_occ_mod, fit.stat = "freeman-tukey", group = 1) 

# View
summary(ahdrift_occ_out_plt)
str(ahdrift_occ_out_plt)

# Convert to a data frame
ahdrift_occ_out_plt_tbl <- data.frame(
  fit = ahdrift_occ_out_plt$fit.y,
  fit.rep = ahdrift_occ_out_plt$fit.y.rep
) |> 
  tibble()

# View
glimpse(ahdrift_occ_out_plt_tbl)

# Pivot Longer
ahdrift_occ_out_plt_lng <- ahdrift_occ_out_plt_tbl |> 
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "species"),
    names_pattern = "^(fit\\.rep|fit)\\.(\\d+)$"
  ) %>%
  mutate(
    species = as.numeric(species),
    rep.greater = fit.rep > fit
  )
# View
glimpse(ahdrift_occ_out_plt_lng)

# Visualize
ahdrift_occ_out_plt_lng |> 
  ggplot(aes(x = fit, y = fit.rep)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.8) +
  geom_point(aes(fill = rep.greater), shape = 21, color = "black", size = 2.5, alpha = 0.8) +
  scale_fill_manual(
    values = c("FALSE" = "lightskyblue1", "TRUE" = "lightsalmon"),
    guide = "none" 
  ) +
  labs(
    x = "True",
    y = "Fit",
    title = "Posterior Predictive Check by Plot"
  ) +
  theme_classic() 

# Which plots contribute to the poor goodness of fit? 
ahdrift_occ_out_fit_plt <- ahdrift_occ_out_plt$fit.y.rep.group.quants[3, , ] - ahdrift_occ_out_plt$fit.y.group.quants[3, , ]
plot(ahdrift_occ_out_fit_plt, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')

# 2.3 Repeat at the visit level ----------------------------------------------------
ahdrift_occ_out_vst <- ppcOcc(ahdrift_occ_mod, fit.stat = "freeman-tukey", group = 2) 

# View
summary(ahdrift_occ_out_vst)
str(ahdrift_occ_out_vst)

# Convert to a data frame
ahdrift_occ_out_vst_tbl <- data.frame(
  fit = ahdrift_occ_out_vst$fit.y,
  fit.rep = ahdrift_occ_out_vst$fit.y.rep
) |> 
  tibble()

# View
glimpse(ahdrift_occ_out_vst_tbl)

# Pivot Longer
ahdrift_occ_out_vst_lng <- ahdrift_occ_out_vst_tbl |> 
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "species"),
    names_pattern = "^(fit\\.rep|fit)\\.(\\d+)$"
  ) %>%
  mutate(
    species = as.numeric(species),
    rep.greater = fit.rep > fit
  )
# View
glimpse(ahdrift_occ_out_vst_lng)

# Visualize
ahdrift_occ_out_vst_lng |> 
  ggplot(aes(x = fit, y = fit.rep)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.8) +
  geom_point(aes(fill = rep.greater), shape = 21, color = "black", size = 2.5, alpha = 0.8) +
  scale_fill_manual(
    values = c("FALSE" = "lightskyblue1", "TRUE" = "lightsalmon"),
    guide = "none" 
  ) +
  labs(
    x = "True",
    y = "Fit",
    title = "Posterior Predictive Check by replicate"
  ) +
  theme_classic() 

################################################################################
# 3) Predictions ###############################################################
################################################################################

# View the output again
summary(ahdrift_occ_mod)

# Define a credible interval width
ci_max <- 0.875
ci_min <- 1 - ci_max

# Model formulas 
occ_formu <- ~ factor(plot.type) + total.effort
det_formu <- ~ date + I(date^2) + factor(daily.effort)

# Generate sample plots for predictions
pred_df <- data.frame(
  plot.type = factor(c("IF", "RE", "TE", "TO"), levels = c("IF", "RE", "TE", "TO"))
) |> 
  mutate(total.effort = rep(0, nrow(pred_df)))

# Generate design matrix using your original occurrence formula
X.0 <- model.matrix(occ_formu, data = pred_df)
# View
X.0

# Generate average occupancy predictions at those locations
ahdrift_occ_pred <- predict(ahdrift_occ_mod, X.0)

# Separate the occupancy probabilities
ahdrift_occ_prob_pred <- ahdrift_occ_pred$psi.0.samples
colnames(ahdrift_occ_prob_pred) <- sp_ls
# View
ahdrift_occ_prob_pred
str(ahdrift_occ_prob_pred)

# Summarize by plot type
ahdrift_occ_prob_pred_if <- data.frame(ahdrift_occ_prob_pred[, , 1]) |> 
  tibble() |> 
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Plot.Type = "Interior Forest")
ahdrift_occ_prob_pred_re <- data.frame(ahdrift_occ_prob_pred[, , 2]) |> 
  tibble() |>
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Plot.Type = "Reference Edge")
ahdrift_occ_prob_pred_te <- data.frame(ahdrift_occ_prob_pred[, , 3]) |> 
  tibble() |>
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Plot.Type = "Turbine Edge")
ahdrift_occ_prob_pred_to <- data.frame(ahdrift_occ_prob_pred[, , 4]) |> 
  tibble() |>
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Plot.Type = "Turbine Opening")

# Combine
ahdrift_occ_prob_pred_full <- bind_rows(
  ahdrift_occ_prob_pred_if,
  ahdrift_occ_prob_pred_re,
  ahdrift_occ_prob_pred_te,
  ahdrift_occ_prob_pred_to
)

# View
glimpse(ahdrift_occ_prob_pred_full)
count(ahdrift_occ_prob_pred_full, Plot.Type)

# Summarize all samples by species
ahdrift_occ_prob_pred_sum <- ahdrift_occ_prob_pred_full |> 
  group_by(Plot.Type, Species) |> 
  reframe(Mean = mean(Psi),
          CI.Low = quantile(Psi, probs = ci_min),
          CI.High = quantile(Psi, probs = ci_max)
  ) 
# View
glimpse(ahdrift_occ_prob_pred_sum)

################################################################################
# 4) Posterior Plots ###########################################################
################################################################################

# View the data again
glimpse(ahdrift_occ_prob_pred_sum)

# 3.1) Occupancy Probability by plot type --------------------------------------

# Palette
wind_pal <- c(
  "Interior Forest" = "forestgreen",
  "Reference Edge" = "goldenrod3",
  "Turbine Edge" = "violetred3",
  "Turbine Opening" = "darkorchid4"
  )

# Make the plot
ahdrift_occ_prob_pred_fig <- ahdrift_occ_prob_pred_sum |>
  mutate(Species = str_replace_all(Species, fixed("."), " ")) |>
  mutate(Species = factor(Species)) |> 
  mutate(Species = fct_reorder(Species, desc(Species))) |> 
  ggplot(aes(y = Species)) +
  # Add points at the mean values for each parameter
  geom_point(
    aes(
      x = Mean, 
      colour = Plot.Type
    ), 
    shape = 15, 
    size = 1.2, 
    alpha = 0.8
  ) +
  # Add whiskers for Credible intervals
  geom_linerange(
    aes(
      xmin = CI.Low,
      xmax = CI.High,
      # linetype = Supported,
      colour = Plot.Type
    ),
    linewidth = 0.4,
    alpha = 0.7
  ) +
  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = wind_pal) +
  # Change the Labels
  labs(x = "Parameter Estimate", y = "") + 
  # Simple theme
  theme_bw() +
  # Edit theme
  theme(
    legend.position = "top", 
    legend.text = element_text(size = 10), 
    legend.title = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10)
  ) +
  facet_wrap(~Plot.Type)

# View the plot 
ahdrift_occ_prob_pred_fig

# Save the plot as a png
ggsave(plot = ahdrift_occ_prob_pred_fig,
       path(fig_dir, "ahdrift_occupancy_prob_plot_type_whiskers.png"),
       width = 200,
       height = 175,
       units = "mm",
       dpi = 300)


# 3.2) Difference between turbines and other plots -----------------------------------

# Palette
signif_pal <- c(
  "Supported Difference" = "darkblue",
  "No Supported Difference" = "gray30"
)

# Extract differences between turbines and the other plot types
ahdrift_occ_diff_te_if <- data.frame(ahdrift_occ_prob_pred[, , 3] - ahdrift_occ_prob_pred[, , 1]) |> 
  tibble() |> 
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Psi.Diff = "Turbine Edge vs Interior Forest")
ahdrift_occ_diff_te_re <- data.frame(ahdrift_occ_prob_pred[, , 3] - ahdrift_occ_prob_pred[, , 2]) |> 
  tibble() |> 
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Psi.Diff = "Turbine Edge vs Reference Edge")
ahdrift_occ_diff_to_if <- data.frame(ahdrift_occ_prob_pred[, , 4] - ahdrift_occ_prob_pred[, , 1]) |> 
  tibble() |> 
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Psi.Diff = "Turbine Opening vs Interior Forest")
ahdrift_occ_diff_to_re <- data.frame(ahdrift_occ_prob_pred[, , 4] - ahdrift_occ_prob_pred[, , 2]) |> 
  tibble() |> 
  pivot_longer(names_to = "Species", values_to = "Psi", cols = everything()) |> 
  mutate(Psi.Diff = "Turbine Opening vs Reference Edge")
# View
glimpse(ahdrift_occ_diff_te_if)
glimpse(ahdrift_occ_diff_te_re)
glimpse(ahdrift_occ_diff_to_if)
glimpse(ahdrift_occ_diff_to_re)

# Combine
ahdrift_occ_diff_full <- bind_rows(
  ahdrift_occ_diff_te_if,
  ahdrift_occ_diff_te_re,
  ahdrift_occ_diff_to_if,
  ahdrift_occ_diff_to_re
) |> 
  group_by(Species, Psi.Diff) |> 
  reframe(Mean = mean(Psi),
          CI.Low = quantile(Psi, probs = ci_min),
          CI.High = quantile(Psi, probs = ci_max)
  )  |> 
  distinct() |> 
  mutate(Supported = case_when(CI.High*CI.Low > 0 ~ "Supported Difference",
                               CI.High*CI.Low <= 0 ~ "No Supported Difference"
  )) 
# View
glimpse(ahdrift_occ_diff_full)

# Make the plot
ahdrift_occ_diff_plot_whisker <- ahdrift_occ_diff_full |>
  mutate(Species = str_replace_all(Species, fixed("."), " ")) |>
  mutate(Species = factor(Species)) |> 
  mutate(Species = fct_reorder(Species, desc(Species))) |> 
  ggplot(aes(y = Species)) +
  # Add points at the mean values for each parameter
  geom_point(
    aes(x = Mean, color = Supported), 
    shape = 15, 
    size = 1.2, 
    alpha = 0.8
  ) +
  # Add whiskers for Credible intervals
  geom_linerange(
    aes(
      xmin = CI.Low,
      xmax = CI.High,
      colour = Supported
    ),
    linewidth = 0.4,
    alpha = 0.7
  ) +
  # Change colors
  scale_color_manual(values = signif_pal) +
  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  # Change the Labels
  labs(x = "Parameter Estimate", y = "") + 
  # Simple theme
  theme_bw() +
  # Edit theme
  theme(
    legend.position = "top", 
    legend.text = element_text(size = 10), 
    legend.title = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10)
  ) +
  facet_wrap(~Psi.Diff)

# View the plot 
ahdrift_occ_diff_plot_whisker

# Save the plot as a png
ggsave(plot = ahdrift_occ_diff_plot_whisker,
       path(fig_dir, "ahdrift_occupancy_diff_plot_type_whiskers.png"),
       width = 200,
       height = 175,
       units = "mm",
       dpi = 300)

