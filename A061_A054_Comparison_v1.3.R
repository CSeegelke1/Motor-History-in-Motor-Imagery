# Title: Distinct contributions of motor imagery and execution to history-dependent biases in reaching; Experiment 1 (A061) & Experiment 2 (A054)
# Authors: Seegelke, Heed 
# Script: Christian Seegelke 01/09/2026
# INPUT:  A061_data.csv, A054_data.csv (preprocessed data), A061_SubjInfo.xlsx, A054_SubjInfo.xlsx (Demographics & MI Questionnaire data)
# OUTPUT: STATS OF EXP 1 & EXP 2; FIGURE 2, FIGURE 4, SUPPLEMENTARY FIGURES S2, S3, S4, S5
#========================================================================================================================
#===================================================================================================================


# ==== REPRODUCABILITY =====
set.seed(1234)
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  scipen = 999,            # fewer scientific notation surprises
  warn = 1                 # show warnings as they occur
)

# ==== INSTALL PACKAGES IF MISSING =====
pkgs <- c(
  "tidyverse","rio","psych","afex","emmeans","cowplot","patchwork","ggpubr","ggpattern",
  "sdamr","brms","GGally","see","bayestestR","tidybayes","BayesFactor",
  "flextable","officer","conflicted","renv","scales","rlang"
)
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install)

# ==== LOAD PACKAGES =====
suppressPackageStartupMessages({
  library(tidyverse)
  library(rio)                      # for data I/O
  library(psych)
  library(afex)
  library(emmeans)
  library(cowplot)
  library(patchwork)               # for arranging plots
  library(ggpubr)                  # for arranging plots
  library(ggpattern)               # for pattern in plots
  library(sdamr)                   # for position_jitternudge in plotting
  library(brms)                    # Bürkner, P. (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80, 1
  library(GGally)                  # for plot checks of multicollinearity
  library(see)
  library(bayestestR)              # Makowski, D., Ben-Shachar, M. S., &  Lüdecke, D.(2019). bayestestR: Describing Effects and their Uncertainty, Existence and Significance within the Bayesian Framework. The Journal of Open Source Software, 4(40), 1541.
  library(tidybayes)
  library(BayesFactor)             # Richard D. Morey, Jeffrey N. Rouder 
  library(flextable)               # for creating nice tables
  library(officer)                 # for exporting tables to word
  #library(conflicted)              # tidyverse plus MASS/plyr/psych/see/ggpubr often create masked objects. Use conflicted to make masking explicit and prefer the intended functions.
  library(scales)
  library(rlang)
  library(purrr)
})


# ==== Resolve common conflicts explicitly =====
# Prefer dplyr verbs over stats/MASS
#conflict_prefer("filter", "dplyr")
#conflict_prefer("select", "dplyr")
#conflict_prefer("lag", "dplyr")
#conflict_prefer("summarise", "dplyr")
#conflict_prefer("rename", "dplyr")
#conflict_prefer("align_plots", "patchwork")
#conflict_prefer("alpha", "scales")
#conflict_prefer("ar", "brms")
#conflict_prefer("cs", "brms")
#conflict_prefer("dstudent_t", "brms")
#conflict_prefer("border", "ggpubr")
#conflict_prefer("col_factor", "scales")
#conflict_prefer("compose", "purrr")
#conflict_prefer("discard", "purrr")
# Prefer broom tidiers if used alongside others
#if ("broom" %in% installed.packages()[,"Package"]) {
#  library(broom)
#  conflict_prefer("tidy", "broom")
#  conflict_prefer("glance", "broom")
#}
## Prefer emmeans::emmeans explicitly
#conflict_prefer("emmeans", "emmeans")
## psych::describe vs Hmisc::describe (if Hmisc is loaded elsewhere)
#conflict_prefer("describe", "psych")

# ==== afex / contrasts / emmeans defaults =====
# Set Type-III SS compatible contrasts for factorial designs
options(contrasts = c("contr.sum", "contr.poly"))
afex::afex_options(
  type = 3,
  es_aov = "ges",
  return_aov = "afex_aov",
  emmeans_model = "multivariate"
)


# ==== brms / Stan options =====
# Use all cores and a reasonable default backend
options(mc.cores = parallel::detectCores())
# Faster compilation for development; consider O3 for final runs
Sys.setenv(LOCAL_CPPFLAGS = "-O2")
# If you use cmdstanr instead of rstan, set: brms.backend = "cmdstanr"
# options(brms.backend = "cmdstanr")





#==== A061: IMPORT SINGLE TRIAL DATA  ====
# clear workspace
# rm(list = ls(all = TRUE))  # Generally avoid in shared scripts

# ==== Paths (portable) =====
basePath   <- "C:/Experiments/A061_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")

# ==== Load single-trial data =====
# Note: N = 32, 704 trials (incl. 64 practise trials) per subject, subject#07 only 640 trial = 22464 trials
# Expecting column 'subID' present
d <- import(file.path(dataPath, "A061_data.csv")) %>%
  as_tibble()






# ==== Load demographics (sheet 2) =====
# Ensure there is a subject identifier column to join on.
# If sheet lacks an explicit subID column but rows correspond to subjects 1..32,
# create it via row_number().
d_SubInfo_raw <- import(file.path(docPath, "A061_SubjInfo.xlsx"), which = 2) %>%
  as_tibble()

if (!"subID" %in% names(d_SubInfo_raw)) {
  d_SubInfo <- d_SubInfo_raw %>%
    mutate(subID = row_number())
} else {
  d_SubInfo <- d_SubInfo_raw
}

# Keep only needed columns and standardize names
d_SubInfo <- d_SubInfo %>%
  rename(
    Age = Age,
    Gender = Gender,
    Handedness = Handedness,
    EHI = EHI
  ) %>%
  select(subID, Age, Gender, Handedness, EHI)

# ==== Merge demographics to trial data =====
d <- d %>%
  left_join(d_SubInfo, by = "subID")



# ==== Quick checks =====
# 1) Remove training block and subject 07 BEFORE converting subID to factor
#    (or convert with character comparison)
d <- d %>%
  filter(blocks_thisN != 0)


# If subID is numeric here, remove by numeric id, otherwise coerce to character
if (is.numeric(d$subID)) {
  d <- filter(d, subID != 7)
} else {
  d <- filter(d, as.character(subID) != "7")
}

# 2) Recode to factor in a vectorized way
d <- d %>%
  mutate(
    across(
      c(stop_signal_prime, obstacle_prime, obstacle_probe,
        target_location_probe, target_location_prime,
        deg_targets_prime, deg_targets_probe,
        Error_stop_signal, Error_too_slow_prime, Error_obstacle_prime,
        Error_too_slow_probe, Error_obstacle_probe,
        Error_MovOnsetDetection_Prime, Error_MovOnsetDetection_Probe,
        Error_Any, Gender, Handedness),
      ~ factor(.)
    ),
    subID = factor(subID)  # after filtering
  )



# 3) Compute TRT_Prime without a loop
#    MI trials (stop) use MIT_Prime; Execution trials (go) sum RT + MT + MT_back
#    Make sure these components are numeric
num_cols <- c("MIT_Prime","RT_Prime","MT_Prime","MT_back_Prime")
d[num_cols] <- lapply(d[num_cols], function(x) suppressWarnings(as.numeric(x)))

d <- d %>%
  mutate(
    TRT_Prime = ifelse(
      stop_signal_prime == "stop",
      MIT_Prime,
      RT_Prime + MT_Prime + MT_back_Prime
    )
  )


# 4) Drop unused levels after filtering
d <- d %>%
  mutate(across(where(is.factor), fct_drop))


# 5) Quick checks
# Check subject ids remaining
print(levels(d$subID))
# Verify TRT has no unintended NAs beyond missing inputs
sum(is.na(d$TRT_Prime))



# Adjust signed Reach difference
# signed ReachDiff should be between -180 and 180 degrees
plot(d$ReachDiff_sig_Prime)
plot(d$ReachDiff2_sig_Prime)
plot(d$ReachDiff_sig_Probe)
plot(d$ReachDiff2_sig_Probe)

wrap180 <- function(x) ((x + 180) %% 360) - 180
abs180  <- function(x) abs(wrap180(x))   # minimal absolute angular deviation

d$ReachDiff_sig_Prime  <- wrap180(d$ReachDiff_sig_Prime)
d$ReachDiff2_sig_Prime <- wrap180(d$ReachDiff2_sig_Prime)
d$ReachDiff_sig_Probe  <- wrap180(d$ReachDiff_sig_Probe)
d$ReachDiff2_sig_Probe <- wrap180(d$ReachDiff2_sig_Probe)


# Calculate Final Reach Error
circ_diff_deg <- function(a, b) wrap180(a - b)
d$EndPointError_sig_Prime <- circ_diff_deg(
  d$V_HandMovOnsetTargetCenter_Prime_Ori,
  d$V_HandMovOnsetHandTargetReached_Prime_Ori
)
d$EndPointError_abs_Prime <- abs(d$EndPointError_sig_Prime)

d$EndPointError_sig_Probe <- circ_diff_deg(
  d$V_HandMovOnsetTargetCenter_Probe_Ori,
  d$V_HandMovOnsetHandTargetReached_Probe_Ori
)
d$EndPointError_abs_Probe <- abs(d$EndPointError_sig_Probe)

plot(d$EndPointError_sig_Prime);  plot(d$EndPointError_abs_Prime)
plot(d$EndPointError_sig_Probe);  plot(d$EndPointError_abs_Probe)


# TIME BETWEEN MOVEMENT ONSET AND CURSOR OUT OF START
# movement onset should not be later than cursor out of start => values should lie under diagonal
# Prime movements
plot((d$frameNum_prime_time_cursor_out_of_start)-d$frameNum_prime_target_onset,d$mov_onset_Prime)
abline(coef = c(0,1))
#Probe movements
plot(d$frameNum_probe_time_cursor_out_of_start,d$mov_onset_Probe)
abline(coef = c(0,1))

# difference between time cursor out of start and movement onset (should be positive and close to zero; converted to ms)
d$time_diff_mov_onset_cursor_out_of_start_Prime <- (( (d$frameNum_prime_time_cursor_out_of_start) - d$frameNum_prime_target_onset) - d$mov_onset_Prime) *1000/d$frameRate
plot(d$time_diff_mov_onset_cursor_out_of_start_Prime)
d$time_diff_mov_onset_cursor_out_of_start_Probe <-    (d$frameNum_probe_time_cursor_out_of_start    - d$mov_onset_Probe) *1000/d$frameRate
plot(d$time_diff_mov_onset_cursor_out_of_start_Probe)


# ==== REMOVE ERROR TRIALS ====
d_raw <- d

# 1) Normalize error flags to logicals
# Works if columns are logical, numeric 0/1, or factor/character "0"/"1"
error_cols <- c("Error_stop_signal","Error_too_slow_prime","Error_obstacle_prime",
                "Error_too_slow_probe","Error_obstacle_probe",
                "Error_MovOnsetDetection_Prime","Error_MovOnsetDetection_Probe",
                "Error_Any")

to_logical_01 <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  # factors/characters: coerce to character then numeric
  xv <- suppressWarnings(as.numeric(as.character(x)))
  return(xv == 1)
}

d <- d %>%
  mutate(across(all_of(error_cols), to_logical_01))

# 2) Compute denominators
n_total <- nrow(d)
n_stop  <- sum(d$stop_signal_prime == "stop", na.rm = TRUE)
n_go    <- sum(d$stop_signal_prime == "go",  na.rm = TRUE)

# 3) Helper to print counts and percents
report <- function(n, denom, label_n, label_pct) {
  cat(paste0(label_n, ": ", n, "\n"))
  if (!is.na(denom) && denom > 0) {
    cat(paste0(label_pct, ": ", round(100 * n / denom, 1), "%\n\n"))
  } else {
    cat("(denominator is 0 or NA)\n\n")
  }
}

cat(paste0("Total number of trials: ", n_total, "\n\n"))

# 4) Reports
n_err_stop   <- sum(d$Error_stop_signal, na.rm = TRUE)
report(n_err_stop, n_stop,
       "Number of MI error (stop-signal) trials",
       "Percent of stop-signal trials with MI error")

n_err_slow_pr <- sum(d$Error_too_slow_prime, na.rm = TRUE)
report(n_err_slow_pr, n_go,
       "Number of Too Slow Prime trials",
       "Percent of go (prime) trials too slow")

n_err_obst_pr <- sum(d$Error_obstacle_prime, na.rm = TRUE)
report(n_err_obst_pr, n_total,
       "Number of Obstacle Hit Prime trials",
       "Percent of all trials (Obstacle Hit Prime)")

n_err_slow_pb <- sum(d$Error_too_slow_probe, na.rm = TRUE)
report(n_err_slow_pb, n_total,
       "Number of Too Slow Probe trials",
       "Percent of all trials (Too Slow Probe)")

n_err_obst_pb <- sum(d$Error_obstacle_probe, na.rm = TRUE)
report(n_err_obst_pb, n_total,
       "Number of Obstacle Hit Probe trials",
       "Percent of all trials (Obstacle Hit Probe)")

n_err_mo_pr <- sum(d$Error_MovOnsetDetection_Prime, na.rm = TRUE)
report(n_err_mo_pr, n_total,
       "Number of Error Mov Onset Prime trials",
       "Percent of all trials (Mov Onset Prime error)")

n_err_mo_pb <- sum(d$Error_MovOnsetDetection_Probe, na.rm = TRUE)
report(n_err_mo_pb, n_total,
       "Number of Error Mov Onset Probe trials",
       "Percent of all trials (Mov Onset Probe error)")

n_err_any <- sum(d$Error_Any, na.rm = TRUE)
report(n_err_any, n_total,
       "Number of All Error trials",
       "Percent of all trials with any error")

# 5) Remove error trials
d <- d %>% filter(!Error_Any)


# Keep a copy before non-smooth exclusion
d_nonsmooth <- d

# Ensure numeric columns for peak velocity counts
to_num <- function(x) suppressWarnings(as.numeric(x))
d <- d %>%
  mutate(
    num_peak_vel_Probe = to_num(num_peak_vel_Probe),
    num_peak_vel_Prime = to_num(num_peak_vel_Prime)
  )

# Filter: probe must have <= 2 peaks; prime must have <= 2 peaks or missing (NA/NaN)
d <- d %>%
  filter(
    num_peak_vel_Probe <= 2 | is.na(num_peak_vel_Probe),
    num_peak_vel_Prime <= 2 | is.na(num_peak_vel_Prime)
  )


n_before_nonsmooth <- nrow(d_nonsmooth)
n_after_nonsmooth  <- nrow(d)
n_excl_nonsmooth   <- n_before_nonsmooth - n_after_nonsmooth

cat("Number of trials before non-smooth exclude: ", n_before_nonsmooth, "\n")
cat("Number of trials after  non-smooth exclude: ", n_after_nonsmooth,  "\n")
cat("Number of trials excluded (non-smooth): ", n_excl_nonsmooth, "\n")
cat("Percent non-smooth (of total trials): ",
    round(100 * n_excl_nonsmooth / n_total, 1), "%\n\n")


# Total excluded relative to original d_raw
n_total_excl <- nrow(d_raw) - n_after_nonsmooth
cat("Total trials excluded since d_raw: ", n_total_excl, "\n")
cat("Percent of trials excluded since d_raw: ",
    round(100 * n_total_excl / nrow(d_raw), 1), "%\n")










# ==== OUTLIER REMOVAL ====
to_num <- function(x) suppressWarnings(as.numeric(x))

# Make a cleaned numeric copy
d_num <- d %>%
  mutate(
    RT_Probe   = to_num(RT_Probe),
    MT_Probe   = to_num(MT_Probe),
    RT_Prime   = to_num(RT_Prime),
    MT_Prime   = to_num(MT_Prime)
  )

n_all <- nrow(d_num)

# Define keep rules
keep_probe <- with(d_num, RT_Probe > 100 & MT_Probe <= 1000)

# For prime-phase thresholds, assume:
# - apply thresholds only to go trials,
# - let stop trials pass regardless of RT_Prime/MT_Prime (which may be NA)
keep_prime <- with(d_num,
                   ifelse(stop_signal_prime == "go",
                          RT_Prime > 100 & MT_Prime <= 1000,
                          TRUE)
)

dClean <- d_num %>% filter(keep_probe & keep_prime)

# Reports
cat("Anzahl Trials gesamt: ", n_all, "\n")
cat("Anzahl Trials nach Cleaning: ", nrow(dClean), "\n")
cat("Number of outlier trials: ", n_all - nrow(dClean), "\n")
cat("Prozent eliminiert: ", round(100 * (n_all - nrow(dClean)) / n_all, 2), "%\n")

dClean.A061 <- dClean # backup

# ==== DEMOGRAPHICS OF FINAL SAMPLE ====
# One row per subject
subjects <- dClean %>%
  distinct(subID, Handedness, Age, Gender, EHI)

# Quick check: each subject appears once
stopifnot(!any(duplicated(subjects$subID)))

# Make sure Age/EHI are numeric (if read as character)
subjects <- subjects %>%
  mutate(
    Age = suppressWarnings(as.numeric(Age)),
    EHI = suppressWarnings(as.numeric(EHI)),
    Gender = as.character(Gender),
    Handedness = as.character(Handedness)
  )

# Overall demographics
demo_overall <- subjects %>%
  summarise(
    n_subjects = n(),
    female  = sum(Gender %in% c("f","F","female","Female"), na.rm = TRUE),
    male    = sum(Gender %in% c("m","M","male","Male"), na.rm = TRUE),
    diverse = sum(Gender %in% c("d","D","diverse","nonbinary","non-binary"), na.rm = TRUE),
    meanAge = mean(Age, na.rm = TRUE),
    sdAge   = sd(Age, na.rm = TRUE),
    minAge  = min(Age, na.rm = TRUE),
    maxAge  = max(Age, na.rm = TRUE),
    meanEHI = mean(EHI, na.rm = TRUE),
    sdEHI   = sd(EHI, na.rm = TRUE),
    minEHI  = min(EHI, na.rm = TRUE),
    maxEHI  = max(EHI, na.rm = TRUE)
  )

# By handedness
subj_by_hand <- subjects %>%
  group_by(Handedness) %>%
  summarise(
    n_subjects = n(),
    female  = sum(Gender %in% c("f","F","female","Female"), na.rm = TRUE),
    male    = sum(Gender %in% c("m","M","male","Male"), na.rm = TRUE),
    diverse = sum(Gender %in% c("d","D","diverse","nonbinary","non-binary"), na.rm = TRUE),
    meanAge = mean(Age, na.rm = TRUE),
    sdAge   = sd(Age, na.rm = TRUE),
    meanEHI = mean(EHI, na.rm = TRUE),
    sdEHI   = sd(EHI, na.rm = TRUE),
    .groups = "drop"
  )

# Show demographics
demo_overall
subj_by_hand






#==== SETTINGS FOR PLOTTING =======
#define some graphical params like themes
custom_plot_theme <- theme(strip.background =element_rect(fill="white", linewidth = 2),
                           strip.text = element_text(size = rel(1), margin = margin(1,5,5,0, "pt")), #in ggplot2 clockwise starting from top: trbl
                           plot.title = element_text(size = rel(1.5)),
                           panel.background = element_blank())

color_exe1     <- "#008000"
color_mi1      <- "#000080"
color_Neut     <- "#9C9C9C"
color_exe_mi1  <- "#008080"
color_None     <- "#808000"
color_exe_None <- "#804000"
color1_exe_mi  <- "#800080"

# For prime actions
color_exe2     <- "#00c000"
color_mi2      <- "#0000c0"
color_exe_mi2  <- "#00c0c0"
color_Neut2    <- "#BDB8AD"





#==== A061: FINAL REACH ERROR AS A FUNCTION OF INITIAL REACH ERROR ====
q_all <- dClean %>%
  filter(stop_signal_prime %in% c("go", "stop")) %>%
  select(subID,
         stop_signal_prime,
         obstacle_prime,
         obstacle_probe,
         EndPointError_abs_Probe,
         ReachDiff2_Probe) %>%
  
  # drop NA RTs before binning
  filter(!is.na(ReachDiff2_Probe)) %>%
  # create subject-wise RT deciles within each condition
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  # optionally enforce a minimum number of trials per subject-condition
  # filter(dplyr::n() >= 20) %>%
  mutate(Mydeciles = ntile(ReachDiff2_Probe, 10)) %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe, Mydeciles) %>%
  summarise(
    IREBin    = median(ReachDiff2_Probe, na.rm = TRUE),
    FREEffect = median(EndPointError_abs_Probe, na.rm = TRUE),
    nTrials  = dplyr::n(),              # for diagnostics
    .groups  = "drop"
  ) %>%
  # add human-readable labels
  mutate(
    trial_type  = if_else(stop_signal_prime == "go", "Execution", "Motor Imagery"),
    prime       = if_else(obstacle_prime == "yes", "with obstacle", "without obstacle"),
    probe       = if_else(obstacle_probe == "yes", "with obstacle", "without obstacle"),
  )


# group-level aggregation
q_all_agg <- q_all %>%
  group_by(trial_type, prime, probe, Mydeciles) %>%
  summarise(
    meanBin      = mean(IREBin, na.rm = TRUE),
    medianBin    = median(IREBin, na.rm = TRUE),
    meanEffect   = mean(FREEffect, na.rm = TRUE),
    medianEffect = median(FREEffect, na.rm = TRUE),
    IQREffect    = IQR(FREEffect, na.rm = TRUE),
    seEffect     = sd(FREEffect, na.rm = TRUE) / sqrt(dplyr::n()),
    n_subj       = dplyr::n(),  # number of subject-bins contributing
    .groups      = "drop"
  )




#==== A061: FINAL REACH ERROR AS A FUNCTION OF INITIAL REACH ERROR: PLOTTING (SUPPLEMENTARY FIGURE S2 ABEF) ====
colorValues <- c(color_exe1, color_Neut)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY <- ""

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_GoNoObs.A061 <- q_all_agg %>%
  dplyr::filter(trial_type == "Execution", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 5, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("without obstacle")

g.FRE_IRE_GoNoObs.A061


colorValues <- c(color_exe1, color_Neut)  # adjust if needed
fontSize    <- 9
titleX      <- ""
titleY      <- "Final Reach Error (°)"

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_GoObs.A061 <- q_all_agg %>%
  dplyr::filter(trial_type == "Execution", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 5, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("with obstacle")

g.FRE_IRE_GoObs.A061



colorValues <- c(color_mi1,color_Neut)
fontSize    <- 9
titleX      <- "Initial Reach Error Deciles"
titleY      <- ""

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_StopNoObs.A061 <- q_all_agg %>%
  dplyr::filter(trial_type == "Motor Imagery", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("with obstacle")

g.FRE_IRE_StopNoObs.A061



colorValues <- c(color_mi1,color_Neut)
fontSize    <- 9
titleX      <- "Initial Reach Error Deciles"
titleY      <- "Final Reach Error (°)"

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_StopObs.A061 <- q_all_agg %>%
  dplyr::filter(trial_type == "Motor Imagery", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = TRUE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("with obstacle")

g.FRE_IRE_StopObs.A061












######################################################################################################
#========================================= A061: BAYESIAN REGRESSION MODELS ==========================
######################################################################################################
#==== A061: INITIAL REACH ERROR PROBE ====
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe
# ==== MODEL FITTING: INITIAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (UB=180) LOGNORMAL DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A061 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_noobs$ReachDiff2_Probe)
sd(log(dClean_noobs$ReachDiff2_Probe))
# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = lognormal(),
          data   = dClean_noobs)

# define priors
# Intercept: median exp(2.4) = 11; 1-sigma range exp(2.4 ± 0.5) ≈ [6.7, 18.0]
# Fixed effects: sd = 0.3 ⇒ 1-sigma multiplicative factor exp(±0.3) ≈ ×[0.74, 1.35]
# Group-level SD: exponential(rate = 2) ⇒ mean = 0.5, median = log(2)/2 ≈ 0.347
# Residual log-SD: exponential(rate = 0.8) ⇒ mean = 1/0.8 = 1.25 (matches sd(log X))
# Correlations among random effects
prior_IRE_Probe_NoObs <- c(prior(normal(    2.4,    0.5  ),     class = "Intercept"),  
                           prior(normal(       0,    0.3  ),     class = "b"),          
                           prior(exponential(2),                 class = "sd"),         # mean 0.5 
                           prior(exponential(0.8),               class = "sigma"),      
                           prior(lkj(2),                         class = "cor"))      


# we model the mu here, we include random effects for participants
# runs about 45 min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.IRE_Probe_NoObs_OPRIxSS <- brm( 
  bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_noobs,              # data
  family = lognormal(),               # distribution of the response variable
  prior  = prior_IRE_Probe_NoObs,     # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control   = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.IRE_Probe_NoObs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A061_fit.IRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(IRE_Probe_NoObs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
IRE_Probe_NoObs.posteriors_long <- IRE_Probe_NoObs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_NoObs.EMM <- IRE_Probe_NoObs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_NoObs.EMM



# Conditions
dClean_noobs$subID <- droplevels(dClean_noobs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_noobs$subID),
                                  obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.IRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
IRE_NoObs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)





# ==== MODEL FITTING: INITIAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=180) GAUSSIAN DISTRIBUTION ====
dClean_obs <- dClean %>% filter(obstacle_probe=="yes")

dClean_obs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A061 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_obs$ReachDiff2_Probe)

# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_obs)

# define priors
# Intercept around the observed mean, with moderate width
# Fixed effects (main effects and interaction), mildly weakly informative
# Group-level SDs: with scale ~10 on each random effect under subID
# Residual SD: half-Student-t or half-normal with scale ~15
# Correlations among random effects
prior_IRE_Probe_Obs <-  c(prior(normal(     54,   10  ), class = "Intercept"), 
                          prior(normal(      0,   7.5 ), class = "b"),
                          prior(student_t(3, 0,  10   ), class = "sd"), 
                          prior(lkj(2)                 , class = "cor"),
                          prior(student_t(3, 0,  15   ), class = "sigma"))


# we model the mu here. We include random effects for participants
# runs about 1.5 hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.IRE_Probe_Obs_OPRIxSS <- brm( 
  bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_obs,                # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_IRE_Probe_Obs,       # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.IRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A061_fit.IRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(IRE_Probe_Obs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
IRE_Probe_Obs.posteriors_long <- IRE_Probe_Obs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_Obs.EMM <- IRE_Probe_Obs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_Obs.EMM



# Conditions
dClean_obs$subID <- droplevels(dClean_obs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_obs$subID),
                                  obstacle_prime     = levels(dClean_obs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_obs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.IRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
IRE_Obs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)








# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 2AB) ====
# Prepare and tidy data
IRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
IRE.EMM.A061 <- rbind(IRE_NoObs.EMM,IRE_Obs.EMM)

IRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
IRE.subj.EMM.A061 <- rbind(IRE_NoObs.subj.EMM,IRE_Obs.subj.EMM)

# rename factors
IRE.EMM.A061 <- IRE.EMM.A061  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

IRE.EMM.A061

# rename factors
IRE.subj.EMM.A061 <- IRE.subj.EMM.A061  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.IRE_NoObs.EMM.A061 <-
  IRE.subj.EMM.A061 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(IRE.EMM.A061, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(0, 80),
    breaks = seq(0, 80, 15)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.IRE_NoObs.EMM.A061



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_Obs.EMM.A061 <-
  IRE.subj.EMM.A061 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(IRE.EMM.A061, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(0, 80),
    breaks = seq(0, 80, 15)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

g.IRE_Obs.EMM.A061





# ==== INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
IRE.Contrasts.A061           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(IRE.Contrasts.A061) <- c("Execution", "Motor Imagery", "EX minus MI", 
                             "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                             "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials
IRE.Contrasts.A061$`Execution`           <- (IRE_Probe_NoObs.posteriors$`yes_go` + IRE_Probe_Obs.posteriors$`yes_go`)/2 - ( IRE_Probe_NoObs.posteriors$`no_go` + IRE_Probe_Obs.posteriors$`no_go`)/2
# Motor Imagery trials
IRE.Contrasts.A061$`Motor Imagery`       <- (IRE_Probe_NoObs.posteriors$`yes_stop` + IRE_Probe_Obs.posteriors$`yes_stop`)/2 - ( IRE_Probe_NoObs.posteriors$`no_stop` + IRE_Probe_Obs.posteriors$`no_stop`)/2
# Execution vs Motor Imagery
IRE.Contrasts.A061$`EX minus MI`         <- (IRE.Contrasts.A061$`Execution` - IRE.Contrasts.A061$`Motor Imagery`)
# Execution trials Probe without Obstacle
IRE.Contrasts.A061$`Execution NoObs`     <- (IRE_Probe_NoObs.posteriors$`yes_go` - IRE_Probe_NoObs.posteriors$`no_go`)
# Execution trials with Obstacle
IRE.Contrasts.A061$`Execution Obs`       <- (IRE_Probe_Obs.posteriors$`yes_go` - IRE_Probe_Obs.posteriors$`no_go`)
# Motor Imagery trials Probe without Obstacle
IRE.Contrasts.A061$`Motor Imagery NoObs` <- (IRE_Probe_NoObs.posteriors$`yes_stop` - IRE_Probe_NoObs.posteriors$`no_stop`)
# Motor Imagery trials Probe with Obstacle
IRE.Contrasts.A061$`Motor Imagery Obs`   <- (IRE_Probe_Obs.posteriors$`yes_stop` - IRE_Probe_Obs.posteriors$`no_stop`)
# Execution vs Motor Imagery without Obstacle
IRE.Contrasts.A061$`EX minus MI NoObs`   <- IRE.Contrasts.A061$`Execution NoObs` - IRE.Contrasts.A061$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
IRE.Contrasts.A061$`EX minus MI Obs`     <- IRE.Contrasts.A061$`Execution Obs` - IRE.Contrasts.A061$`Motor Imagery Obs`


IRE.Contrasts_long.A061 <- pivot_longer(IRE.Contrasts.A061, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "Motor Imagery Obs",
  "Execution NoObs", "Motor Imagery NoObs", 
  "EX minus MI Obs", "EX minus MI NoObs", 
  "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
IRE.Contrasts_long.A061$contrast <- factor(IRE.Contrasts_long.A061$contrast, levels = contrast_order)

# Check
head(IRE.Contrasts_long.A061)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
IRE.contrasts.summary.A061 <-
  IRE.Contrasts_long.A061 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

IRE.contrasts.summary.A061$pd <- format(IRE.contrasts.summary.A061$pd, nsmall = 4)  
#print(IRE.contrasts.summary.A061, n = Inf, width = Inf)


#### Calculate ROPE
### Calculate and back-transform ROPE, this is a bit tricky
## this does not take into account that values are on log scale so don't use it
## apparently, back-transformation is now implemented :)
#RR_NoObs <- rope_range(A061_fit.IRE_Probe_NoObs_OPRIxSS) 
# compute a ROPE on the response scale for a lognormal outcome
## calculate manually: # IMPORTANT:SD(x) != exp(SD(log(x)))
dClean$logReachDiff2_Probe <- log(dClean$ReachDiff2_Probe) #log transform raw data
dClean_noobs               <- dClean %>% filter(obstacle_probe=="no")
# geometric mean of raw values ≈ exp(mean(log(x)))
#the geometric mean multiplied by the log-standard deviation. This should approximate the "natural" standard deviation pretty well.
geom_mean <- exp(mean(dClean_noobs$logReachDiff2_Probe))
# SD of log-transformed data
sd_log <- sd(dClean_noobs$logReachDiff2_Probe)
# ROPE = 0.1 * "natural" SD
rope_value_RE_log <- 0.1 * geom_mean * sd_log
rope_value_RE_log

( RR_IRE.NoObs.A061 <- c(-rope_value_RE_log, rope_value_RE_log) )
( RR_IRE.Obs.A061   <- rope_range(A061_fit.IRE_Probe_Obs_OPRIxSS) )
( RR.IRE.A061       <- (RR_IRE.NoObs.A061 + RR_IRE.Obs.A061) /2 )



# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE.A061                       <- as.data.frame(IRE.contrasts.summary.A061)
IRE.contrast_in_ROPE.A061$lowerROPE             <- NA
IRE.contrast_in_ROPE.A061$lowerROPE[c(1, 2, 5)] <- RR_IRE.Obs.A061[1]
IRE.contrast_in_ROPE.A061$lowerROPE[c(3, 4, 6)] <- RR_IRE.NoObs.A061[1]
IRE.contrast_in_ROPE.A061$lowerROPE[c(7:9)]     <- RR.IRE.A061[1]
IRE.contrast_in_ROPE.A061$upperROPE             <- NA
IRE.contrast_in_ROPE.A061$upperROPE[c(1, 2, 5)] <- RR_IRE.Obs.A061[2]
IRE.contrast_in_ROPE.A061$upperROPE[c(3, 4, 6)] <- RR_IRE.NoObs.A061[2]
IRE.contrast_in_ROPE.A061$upperROPE[c(7:9)]     <- RR.IRE.A061[2]
IRE.contrast_in_ROPE.A061$CI_range              <- IRE.contrast_in_ROPE.A061$upper - IRE.contrast_in_ROPE.A061$lower
IRE.contrast_in_ROPE.A061$minUpper              <- IRE.contrast_in_ROPE.A061 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
IRE.contrast_in_ROPE.A061$maxLower              <- IRE.contrast_in_ROPE.A061 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
IRE.contrast_in_ROPE.A061$DiffminUppermaxLower  <- IRE.contrast_in_ROPE.A061$minUpper  - IRE.contrast_in_ROPE.A061$maxLower 
IRE.contrast_in_ROPE.A061$Zeros                 <- rep(0,nrow(IRE.contrast_in_ROPE.A061))
IRE.contrast_in_ROPE.A061$Overlap               <- IRE.contrast_in_ROPE.A061  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
IRE.contrast_in_ROPE.A061$perc_in_ROPE          <- (IRE.contrast_in_ROPE.A061$Overlap*100)/IRE.contrast_in_ROPE.A061$CI_range
IRE.contrast_in_ROPE.A061[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
IRE.subj.contrast.A061 <- IRE.subj.EMM.A061 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "Motor Imagery" & probe == "without obstacle" ~ "Motor Imagery NoObs",
           trial_type == "Motor Imagery" & probe == "with obstacle"    ~ "Motor Imagery Obs")
  )

print(IRE.subj.contrast.A061)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
IRE.subj.contrast_pooled.A061 <- IRE.subj.EMM.A061 %>%
  group_by(subID, trial_type, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution"  ~ "Execution",
           trial_type == "Execution"  ~ "Execution",
           trial_type == "Motor Imagery" ~ "Motor Imagery",
           trial_type == "Motor Imagery" ~ "Motor Imagery"))

print(IRE.subj.contrast_pooled.A061)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
IRE.subj.diffContrast.A061 <- IRE.subj.contrast.A061 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus MI NoObs",
           probe == "with obstacle"  ~ "EX minus MI Obs")
  )

print(IRE.subj.diffContrast.A061)


IRE.subj.diffContrast_pooled.A061 <- IRE.subj.contrast_pooled.A061 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

print(IRE.subj.diffContrast_pooled.A061)
















#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS (FIGURE 2EFIJ) ====
# Probe without obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_NoObs.Contrast.A061 <- IRE.Contrasts_long.A061 %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "Motor Imagery NoObs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A061, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.NoObs.A061[1], xmax = RR_IRE.NoObs.A061[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast.A061, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,22),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs.Contrast.A061


# Probe with obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_Obs.Contrast.A061 <- IRE.Contrasts_long.A061 %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "Motor Imagery Obs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",     # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,       # No border on the stripes
    trim = TRUE,            # Trim the ends of the distributions
    linewidth = 0,          # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A061, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.Obs.A061[1], xmax = RR_IRE.Obs.A061[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast.A061, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,12),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") + #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs.Contrast.A061



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_NoObs_EXvsMI.Contrast.A061 <- IRE.Contrasts_long.A061 %>% filter(contrast=="EX minus MI NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI NoObs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A061, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.NoObs.A061[1], xmax = RR_IRE.NoObs.A061[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast.A061, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,22),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs_EXvsMI.Contrast.A061




# Probe with obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_Obs_EXvsMI.Contrast.A061 <- IRE.Contrasts_long.A061 %>% filter(contrast=="EX minus MI Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI Obs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A061, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.Obs.A061[1], xmax = RR_IRE.Obs.A061[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast.A061, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,12),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs_EXvsMI.Contrast.A061



















#==== A061: REACTION TIME (RT) PROBE ====
# ==== MODEL FITTING: REACTION TIME: 3-FACTOR INTERACTION MODEL WITH SHIFTED LOG-NORMAL DISTRIBUTION ====
dClean %>%
  ggplot( aes(x = RT_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 10, boundary = 0, position = 'identity' ) +
  labs( title = "A061 - Histogram", subtitle = "RT Probe" ) +
  scale_x_continuous( name = 'RT [ms]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean$RT_Probe)
sd(log(dClean$RT_Probe))

# looking at prior values
get_prior(bf(RT_Probe ~ obstacle_prime * obstacle_probe * stop_signal_prime + ( obstacle_prime * obstacle_probe * stop_signal_prime | subID )),
          family = shifted_lognormal(),
          data   = dClean)


# define priors
# the linear predictor is on the log of the shifted RT, i.e., y=log(RT-ndt). SO the Intercept corresponds to the typical value of (RT-ndt)
# Intercept: median exp(5.35), given median RT = 300 and ndt = 85, Intercept should be ~300-85 ~215 ~ 5.35
# Fixed effects: sd = 0.31 ⇒ 1-sigma multiplicative factor exp(±0.1) ≈ ×[0.90, 1.11]
# Group-level SD: reasonable weakly-informative choice.
# Residual log-SD: somewhat smaller than the observed sd(log RT) ≈ 0.254
# Shift: Centering near 90 is pragmatic given min RT ≈ 100.
# Correlations among random effects
prior_RT.Probe      <- c(prior(normal( 5.35,  0.5  ), class = "Intercept"), 
                         prior(normal( 0,     0.1  ), class = "b"),
                         prior(normal( 0,     0.1),   class = "sd"),
                         prior(normal( 0,     0.1),   class = "sigma"),
                         prior(normal(85,    15),     class = "ndt",  lb = 0),  
                         prior(lkj(3),                class = "cor")
)

# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 2.5 hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.RT_Probe_OPRIxOPROxSS <- brm( 
  bf(RT_Probe ~ obstacle_prime * obstacle_probe * stop_signal_prime + ( obstacle_prime * obstacle_probe * stop_signal_prime | subID )),     # model specification
  data   = dClean,                    # data
  family = shifted_lognormal(),       # distribution of the response variable
  prior  = prior_RT.Probe,            # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.RT_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )










#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean$obstacle_prime),
                        obstacle_probe     = levels(dClean$obstacle_probe),
                        stop_signal_prime  = levels(dClean$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
RT_Probe.posteriors <- as.data.frame(fitted(
  A061_fit.RT_Probe_OPRIxOPROxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(obstacle_prime, obstacle_probe, stop_signal_prime, sep = "_", drop = TRUE))
colnames(RT_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
RT_Probe.posteriors_long <- RT_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("obstacle_prime", "obstacle_probe", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
RT.EMM <- RT_Probe.posteriors_long %>%
  group_by(obstacle_prime, obstacle_probe, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
RT.EMM



# Summaries: median and 95% HDI for Execution and Motor Imagery
RT_Probe.posteriors_long %>%
  group_by(stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
dClean$subID    <- droplevels(dClean$subID)
exp.cond.subj   <- expand.grid(subID              = levels(dClean$subID),
                               obstacle_prime     = levels(dClean$obstacle_prime),
                               obstacle_probe     = levels(dClean_obs$obstacle_probe),
                               stop_signal_prime  = levels(dClean$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.RT_Probe_OPRIxOPROxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
RT.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, obstacle_probe, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) ====
# Prepare and tidy data
# rename factors
RT.EMM.A061 <- RT.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

RT.EMM.A061

# rename factors
RT.subj.EMM.A061 <- RT.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.RT_NoObs.EMM.A061 <-
  RT.subj.EMM.A061 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(RT.EMM.A061, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 500),
    breaks = seq(200, 500, 100)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.RT_NoObs.EMM.A061



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Reaction time (ms)"

g.RT_Obs.EMM.A061 <-
  RT.subj.EMM.A061 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(RT.EMM.A061, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 500),
    breaks = seq(200, 500, 100)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.RT_Obs.EMM.A061

















#==== A061: FINAL REACH ERROR PROBE ====
# measured as angular difference btw cursor mov_onset target center and cursor mov_onset cursor at target hit (EndPointError_abs_Probe)
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe
# ==== MODEL FITTING: FINAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=4) GAUSSIAN DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A061- Histogram", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_noobs$EndPointError_abs_Probe)
sd(log(dClean_noobs$EndPointError_abs_Probe))


get_prior(bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_noobs)


# define priors
# Intercept: median exp(0.3) = 1.35; 1-sigma range exp(0.3 ± 0.5) ≈ [0.8, 2.23]
# Fixed effects: sd = 0.2 ⇒ 1-sigma multiplicative effect for a unit change in a predictor ≈ exp(±0.2) = ×[0.82, 1.22]
# Group-level SD: Half-normal(0, 0.2) on SDs of random intercepts/slopes (log scale), Implied per-subject multiplicative spread ≈ exp(±SD) ≈ ×[0.85, 1.17]
# Residual log-SD: Half-normal(0, 0.2) on sigma (log scale)
# Correlations among random effects
prior_FRE          <- c(prior(normal(  1.4,     0.5   ),  class = "Intercept"), 
                        prior(normal(    0,     0.3   ),  class = "b"),
                        prior(normal(0, 0.3),    class = "sd",    lb = 0),
                        prior(normal(0, 0.3),    class = "sigma", lb = 0),
                        prior(lkj(2), class = "cor"))


# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 1 hour
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.FRE_Probe_NoObs_OPRIxSS <- brm( 
  bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_noobs,              # data
  family = gaussian(),               # distribution of the response variable
  prior  = prior_FRE,                 # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.FRE_Probe_NoObs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )







#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A061_fit.FRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_NoObs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_NoObs.posteriors_long <- FRE_Probe_NoObs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_NoObs.EMM <- FRE_Probe_NoObs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_NoObs.EMM



# Conditions
dClean_noobs$subID <- droplevels(dClean_noobs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_noobs$subID),
                                  obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.FRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_NoObs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)









# ==== MODEL FITTING: FINAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH SKEWED GAUSSIAN DISTRIBUTION ====
dClean_obs <- dClean %>% filter(obstacle_probe=="yes" & EndPointError_abs_Probe<=4)

dClean_obs %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A061- Histogram", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_obs$EndPointError_abs_Probe)
sd(log(dClean_obs$EndPointError_abs_Probe))

# looking at prior values
get_prior(bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = skew_normal(),
          data   = dClean_obs)



# define priors
prior_FRE_obs <-      c(prior(normal(  2.2,     0.5   ),  class = "Intercept"), 
                        prior(normal(    0,     0.3   ),  class = "b"),
                        prior(normal(    0,     0.3   ),  class = "sd",    lb = 0),
                        prior(normal(    0,     2     ),  class = "alpha"),
                        prior(normal(    0,     0.5   ),  class = "sigma",lb = 0),
                        prior(lkj(2), class = "cor"))



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 30min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.FRE_Probe_Obs_OPRIxSS <- brm( 
  bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_obs,                # data
  family = skew_normal(),             # distribution of the response variable
  prior  = prior_FRE_obs,             # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.FRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A061_fit.FRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_Obs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_Obs.posteriors_long <- FRE_Probe_Obs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_Obs.EMM <- FRE_Probe_Obs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_Obs.EMM



# Conditions
dClean_obs$subID <- droplevels(dClean_obs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_obs$subID),
                                  obstacle_prime     = levels(dClean_obs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_obs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.FRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_Obs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)






# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 4AB) ====
# Prepare and tidy data
FRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
FRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
FRE.EMM.A061 <- rbind(FRE_NoObs.EMM,FRE_Obs.EMM)

FRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
FRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
FRE.subj.EMM.A061 <- rbind(FRE_NoObs.subj.EMM,FRE_Obs.subj.EMM)

# rename factors
FRE.EMM.A061 <- FRE.EMM.A061  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

FRE.EMM.A061

# rename factors
FRE.subj.EMM.A061 <- FRE.subj.EMM.A061  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Final Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.FRE_NoObs.EMM.A061 <-
  FRE.subj.EMM.A061 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.A061, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(0, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.FRE_NoObs.EMM.A061



#plot Final Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Final Reach Error (°)"

g.FRE_Obs.EMM.A061 <-
  FRE.subj.EMM.A061 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.A061, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(1, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

g.FRE_Obs.EMM.A061













# ==== FINAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
# IMPORTANT!!!!: FOR POOLED CONSTRASTS: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
FRE.Contrasts.A061           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(FRE.Contrasts.A061) <- c("Execution", "Motor Imagery", "EX minus MI", 
                             "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                             "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials
FRE.Contrasts.A061$`Execution`           <- (FRE_Probe_NoObs.posteriors$`yes_go` + FRE_Probe_Obs.posteriors$`no_go`)/2 - ( FRE_Probe_NoObs.posteriors$`no_go` + FRE_Probe_Obs.posteriors$`yes_go`)/2
# Motor Imagery trials
FRE.Contrasts.A061$`Motor Imagery`       <- (FRE_Probe_NoObs.posteriors$`yes_stop` + FRE_Probe_Obs.posteriors$`no_stop`)/2 - ( FRE_Probe_NoObs.posteriors$`no_stop` + FRE_Probe_Obs.posteriors$`yes_stop`)/2
# Execution vs Motor Imagery
FRE.Contrasts.A061$`EX minus MI`         <- (FRE.Contrasts.A061$`Execution` - FRE.Contrasts.A061$`Motor Imagery`)
# Execution trials Probe without Obstacle
FRE.Contrasts.A061$`Execution NoObs`     <- (FRE_Probe_NoObs.posteriors$`yes_go` - FRE_Probe_NoObs.posteriors$`no_go`)
# Execution trials with Obstacle
FRE.Contrasts.A061$`Execution Obs`       <- (FRE_Probe_Obs.posteriors$`yes_go` - FRE_Probe_Obs.posteriors$`no_go`)
# Motor Imagery trials Probe without Obstacle
FRE.Contrasts.A061$`Motor Imagery NoObs` <- (FRE_Probe_NoObs.posteriors$`yes_stop` - FRE_Probe_NoObs.posteriors$`no_stop`)
# Motor Imagery trials Probe with Obstacle
FRE.Contrasts.A061$`Motor Imagery Obs`   <- (FRE_Probe_Obs.posteriors$`yes_stop` - FRE_Probe_Obs.posteriors$`no_stop`)
# Execution vs Motor Imagery without Obstacle
FRE.Contrasts.A061$`EX minus MI NoObs`   <- FRE.Contrasts.A061$`Execution NoObs` - FRE.Contrasts.A061$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
FRE.Contrasts.A061$`EX minus MI Obs`     <- FRE.Contrasts.A061$`Execution Obs` - FRE.Contrasts.A061$`Motor Imagery Obs`


FRE.Contrasts_long.A061 <- pivot_longer(FRE.Contrasts.A061, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "Motor Imagery Obs",
  "Execution NoObs", "Motor Imagery NoObs", 
  "EX minus MI Obs", "EX minus MI NoObs", 
  "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
FRE.Contrasts_long.A061$contrast <- factor(FRE.Contrasts_long.A061$contrast, levels = contrast_order)

# Check
head(FRE.Contrasts_long.A061)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
FRE.contrasts.summary.A061 <-
  FRE.Contrasts_long.A061 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

FRE.contrasts.summary.A061$pd <- format(FRE.contrasts.summary.A061$pd, nsmall = 4)  
#print(FRE.contrasts.summary.A061, n = Inf, width = Inf)


#### Calculate ROPE
( RR_FRE.NoObs.A061 <- rope_range(A061_fit.FRE_Probe_NoObs_OPRIxSS) )
( RR_FRE.Obs.A061   <- rope_range(A061_fit.FRE_Probe_Obs_OPRIxSS) )
( RR.FRE.A061       <- ( RR_FRE.NoObs.A061 + RR_FRE.Obs.A061)/2 )

# Calculate percent in ROPE for contrasts
options(digits=3)
FRE.contrast_in_ROPE.A061                       <- as.data.frame(FRE.contrasts.summary.A061)
FRE.contrast_in_ROPE.A061$lowerROPE             <- NA
FRE.contrast_in_ROPE.A061$lowerROPE[c(1, 2, 5)] <- RR_FRE.Obs.A061[1]
FRE.contrast_in_ROPE.A061$lowerROPE[c(3, 4, 6)] <- RR_FRE.NoObs.A061[1]
FRE.contrast_in_ROPE.A061$lowerROPE[c(7:9)]     <- RR.FRE.A061[1]
FRE.contrast_in_ROPE.A061$upperROPE             <- NA
FRE.contrast_in_ROPE.A061$upperROPE[c(1, 2, 5)] <- RR_FRE.Obs.A061[2]
FRE.contrast_in_ROPE.A061$upperROPE[c(3, 4, 6)] <- RR_FRE.NoObs.A061[2]
FRE.contrast_in_ROPE.A061$upperROPE[c(7:9)]     <- RR.FRE.A061[2]
FRE.contrast_in_ROPE.A061$CI_range              <- FRE.contrast_in_ROPE.A061$upper - FRE.contrast_in_ROPE.A061$lower
FRE.contrast_in_ROPE.A061$minUpper              <- FRE.contrast_in_ROPE.A061 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
FRE.contrast_in_ROPE.A061$maxLower              <-  FRE.contrast_in_ROPE.A061 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
FRE.contrast_in_ROPE.A061$DiffminUppermaxLower  <- FRE.contrast_in_ROPE.A061$minUpper  - FRE.contrast_in_ROPE.A061$maxLower 
FRE.contrast_in_ROPE.A061$Zeros                 <- rep(0,nrow(FRE.contrast_in_ROPE.A061))
FRE.contrast_in_ROPE.A061$Overlap               <- FRE.contrast_in_ROPE.A061  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
FRE.contrast_in_ROPE.A061$perc_in_ROPE          <- (FRE.contrast_in_ROPE.A061$Overlap*100)/FRE.contrast_in_ROPE.A061$CI_range
FRE.contrast_in_ROPE.A061[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
FRE.subj.contrast.A061 <- FRE.subj.EMM.A061 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "Motor Imagery" & probe == "without obstacle" ~ "Motor Imagery NoObs",
           trial_type == "Motor Imagery" & probe == "with obstacle"    ~ "Motor Imagery Obs")
  )

#print(FRE.subj.contrast.A061)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
FRE.subj.contrast_pooled.A061 <- FRE.subj.EMM.A061 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(FRE.subj.contrast_pooled.A061)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
FRE.subj.diffContrast.A061 <- FRE.subj.contrast.A061 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus MI NoObs",
           probe == "with obstacle"  ~ "EX minus MI Obs")
  )

#print(FRE.subj.diffContrast.A061)


FRE.subj.diffContrast_pooled.A061 <- FRE.subj.contrast_pooled.A061 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(FRE.subj.diffContrast_pooled.A061)









#  ==== FINAL REACH ERROR: PLOTTING CONTRASTS  (FIGURE 4EFIJ) ====
# Probe without obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_NoObs.Contrast.A061 <- FRE.Contrasts_long.A061 %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "Motor Imagery NoObs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A061, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.NoObs.A061[1], xmax = RR_FRE.NoObs.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.A061, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs.Contrast.A061


# Probe with obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_Obs.Contrast.A061 <- FRE.Contrasts_long.A061 %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "Motor Imagery Obs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A061, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.Obs.A061[1], xmax = RR_FRE.Obs.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.A061, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs.Contrast.A061



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_NoObs_EXvsMI.Contrast.A061 <- FRE.Contrasts_long.A061 %>% filter(contrast=="EX minus MI NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI NoObs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A061, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.NoObs.A061[1], xmax = RR_FRE.NoObs.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.A061, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs_EXvsMI.Contrast.A061




# Probe with obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_Obs_EXvsMI.Contrast.A061 <- FRE.Contrasts_long.A061 %>% filter(contrast=="EX minus MI Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI Obs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A061, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.Obs.A061[1], xmax = RR_FRE.Obs.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.A061, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs_EXvsMI.Contrast.A061













#==== A061: ACCURACY: END POINT ERROR 150ms PROBE====
# measured as absolute distance of cursor 150ms after target hit to target center (EndPointError_150_Probe)
# ==== MODEL FITTING: ACCURACY: 3-FACTOR INTERACTION MODEL WITH TRUNCATED GAUSSIAN DISTRIBUTION ====
dClean_EPE150 <- dClean %>% filter(EndPointError_150_Probe <=0.5 & Vel_HandTargetReached_150_Probe<=5 )

dClean_EPE150 %>%
  ggplot( aes(x = EndPointError_150_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.02, boundary = 0, position = 'identity' ) +
  labs( title = "A061- Histogram", subtitle = "End Point Error 150ms after Target Hit Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [cm]') +      
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_EPE150$EndPointError_150_Probe)

# looking at prior values
get_prior(bf(EndPointError_150_Probe | trunc(lb = 0, ub = 0.5) ~ obstacle_prime * obstacle_probe * stop_signal_prime + ( obstacle_prime * obstacle_probe * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_EPE150)

# define priors
prior_absEPE150    <- c(prior(normal(      0.2,   0.3   ),  class = "Intercept"), 
                        prior(normal(      0,     0.2   ),  class = "b"),
                        prior(student_t(3, 0,     0.2),  class = "sd"),
                        prior(student_t(3, 0,     0.2),  class = "sigma"),
                        prior(lkj(2),                     class = "cor"))



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 3.25 hours, introducing truncation in prior predicitve checks causes NA's
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.EPE150_Probe_OPRIxOPROxSS <- brm( 
  bf(EndPointError_150_Probe | trunc(lb = 0, ub = 0.5) ~ obstacle_prime * obstacle_probe * stop_signal_prime + ( obstacle_prime * obstacle_probe * stop_signal_prime | subID )),     # model specification
  data   = dClean_EPE150,             # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_absEPE150,           # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.EPE150_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )







#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean$obstacle_prime),
                        obstacle_probe     = levels(dClean$obstacle_probe),
                        stop_signal_prime  = levels(dClean$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
EPE150_Probe.posteriors <- as.data.frame(fitted(
  A061_fit.EPE150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(obstacle_prime, obstacle_probe, stop_signal_prime, sep = "_", drop = TRUE))
colnames(EPE150_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
EPE150_Probe.posteriors_long <- EPE150_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("obstacle_prime", "obstacle_probe", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
EPE150.EMM <- EPE150_Probe.posteriors_long %>%
  group_by(obstacle_prime, obstacle_probe, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
EPE150.EMM



# Summaries: median and 95% HDI for Execution and Motor Imagery
EPE150_Probe.posteriors_long %>%
  group_by(stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
dClean$subID    <- droplevels(dClean$subID)
exp.cond.subj   <- expand.grid(subID              = levels(dClean$subID),
                               obstacle_prime     = levels(dClean$obstacle_prime),
                               obstacle_probe     = levels(dClean$obstacle_probe),
                               stop_signal_prime  = levels(dClean$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.EPE150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
EPE150.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, obstacle_probe, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)





# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S5AB) ====
# Prepare and tidy data
# rename factors
EPE150.EMM.A061 <- EPE150.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    .value      = .value*10, .lower = .lower*10 , .upper = .upper*10  # convert from cm to mm
    # convert from cm to mm
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns
EPE150.EMM.A061

# rename factors
EPE150.subj.EMM.A061 <- EPE150.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    .value      = .value*10, .lower = .lower*10 , .upper = .upper*10  # convert from cm to mm
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.EPE150_NoObs.EMM.A061 <-
  EPE150.subj.EMM.A061 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPE150.EMM.A061, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPE150.EMM.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(2, 3),
    breaks = seq(1, 5, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.EPE150_NoObs.EMM.A061



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Endpoint Error (mm)"

g.EPE150_Obs.EMM.A061 <-
  EPE150.subj.EMM.A061 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPE150.EMM.A061, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPE150.EMM.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(2, 3),
    breaks = seq(1, 5, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

g.EPE150_Obs.EMM.A061















# ==== ACCURACY: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
EPE150.Contrasts.A061           <- data.frame(matrix(ncol = 10, nrow = 40000))
colnames(EPE150.Contrasts.A061) <- c("Execution", "Motor Imagery", "EX minus MI", 
                                "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                                "EX minus MI NoObs", "EX minus MI Obs",
                                "EX vs MI")

# FOR POOLED DATA: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
# Execution trials (calculated as different minus same movement context)
EPE150.Contrasts.A061$`Execution`           <- ((EPE150_Probe.posteriors$`yes_no_go` - EPE150_Probe.posteriors$`no_no_go`) + ( EPE150_Probe.posteriors$`no_yes_go` - EPE150_Probe.posteriors$`yes_yes_go`) )/2
# Motor Imagery trials (calculated as different minus same movement context)
EPE150.Contrasts.A061$`Motor Imagery`       <- ((EPE150_Probe.posteriors$`yes_no_stop` - EPE150_Probe.posteriors$`no_no_stop`) + ( EPE150_Probe.posteriors$`no_yes_stop` - EPE150_Probe.posteriors$`yes_yes_stop`) )/2
# Difference of Execution vs Motor Imagery
EPE150.Contrasts.A061$`EX minus MI`         <- (EPE150.Contrasts.A061$`Execution` - EPE150.Contrasts.A061$`Motor Imagery`)
# Execution trials Probe without Obstacle
EPE150.Contrasts.A061$`Execution NoObs`     <- (EPE150_Probe.posteriors$`yes_no_go` - EPE150_Probe.posteriors$`no_no_go`)
# Execution trials with Obstacle
EPE150.Contrasts.A061$`Execution Obs`       <- (EPE150_Probe.posteriors$`yes_yes_go` - EPE150_Probe.posteriors$`no_yes_go`)
# Motor Imagery trials Probe without Obstacle
EPE150.Contrasts.A061$`Motor Imagery NoObs` <- (EPE150_Probe.posteriors$`yes_no_stop` - EPE150_Probe.posteriors$`no_no_stop`)
# Motor Imagery trials Probe with Obstacle
EPE150.Contrasts.A061$`Motor Imagery Obs`   <- (EPE150_Probe.posteriors$`yes_yes_stop` - EPE150_Probe.posteriors$`no_yes_stop`)
# Execution vs Motor Imagery without Obstacle
EPE150.Contrasts.A061$`EX minus MI NoObs`   <- EPE150.Contrasts.A061$`Execution NoObs` - EPE150.Contrasts.A061$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
EPE150.Contrasts.A061$`EX minus MI Obs`     <- EPE150.Contrasts.A061$`Execution Obs` - EPE150.Contrasts.A061$`Motor Imagery Obs`
# Execution vs Motor Imagery (Overall difference)
EPE150.Contrasts.A061$`EX vs MI`            <- ((EPE150_Probe.posteriors$`yes_yes_go` + EPE150_Probe.posteriors$`yes_no_go` + EPE150_Probe.posteriors$`no_yes_go` + EPE150_Probe.posteriors$`no_no_go` ) /4 ) - ((EPE150_Probe.posteriors$`yes_yes_stop` + EPE150_Probe.posteriors$`yes_no_stop` + EPE150_Probe.posteriors$`no_yes_stop` + EPE150_Probe.posteriors$`no_no_stop` ) /4 )

EPE150.Contrasts_long.A061 <- pivot_longer(EPE150.Contrasts.A061*10, cols = everything(), # convert from cm to mm
                                      names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c("EX vs MI",
                    "Execution Obs", "Motor Imagery Obs",
                    "Execution NoObs", "Motor Imagery NoObs", 
                    "EX minus MI Obs", "EX minus MI NoObs", 
                    "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
EPE150.Contrasts_long.A061$contrast <- factor(EPE150.Contrasts_long.A061$contrast, levels = contrast_order)

# Check
head(EPE150.Contrasts_long.A061)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
EPE150.contrasts.summary.A061 <-
  EPE150.Contrasts_long.A061 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

EPE150.contrasts.summary.A061$pd <- format(EPE150.contrasts.summary.A061$pd, nsmall = 1)  
print(EPE150.contrasts.summary.A061, n = Inf, width = Inf)


#### Calculate ROPE
( RR.A061       <- rope_range(A061_fit.EPE150_Probe_OPRIxOPROxSS)*10 ) # convert from cm to mm


# Calculate percent in ROPE for contrasts
options(digits=3)
EPE150.contrast_in_ROPE.A061                       <- as.data.frame(EPE150.contrasts.summary.A061)
EPE150.contrast_in_ROPE.A061$lowerROPE             <- RR.A061[1]
EPE150.contrast_in_ROPE.A061$upperROPE             <- RR.A061[2]
EPE150.contrast_in_ROPE.A061$CI_range              <- EPE150.contrast_in_ROPE.A061$upper - EPE150.contrast_in_ROPE.A061$lower
EPE150.contrast_in_ROPE.A061$minUpper              <- EPE150.contrast_in_ROPE.A061 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
EPE150.contrast_in_ROPE.A061$maxLower              <- EPE150.contrast_in_ROPE.A061 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
EPE150.contrast_in_ROPE.A061$DiffminUppermaxLower  <- EPE150.contrast_in_ROPE.A061$minUpper  - EPE150.contrast_in_ROPE.A061$maxLower 
EPE150.contrast_in_ROPE.A061$Zeros                 <- rep(0,nrow(EPE150.contrast_in_ROPE.A061))
EPE150.contrast_in_ROPE.A061$Overlap               <- EPE150.contrast_in_ROPE.A061  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
EPE150.contrast_in_ROPE.A061$perc_in_ROPE          <- (EPE150.contrast_in_ROPE.A061$Overlap*100)/EPE150.contrast_in_ROPE.A061$CI_range
EPE150.contrast_in_ROPE.A061[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
EPE150.ExvsMI.subj.contrast.A061 <- EPE150.subj.EMM.A061 %>%
  group_by(subID, trial_type) %>%
  summarise( .value = mean(.value)) %>%
  tidyr::pivot_wider(
    names_from = trial_type,
    values_from = .value
  ) %>%
  mutate(diff = `Execution` - `Motor Imagery`,
         contrast = "EX vs MI"
  )

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
EPE150.subj.contrast.A061 <- EPE150.subj.EMM.A061 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "Motor Imagery" & probe == "without obstacle" ~ "Motor Imagery NoObs",
           trial_type == "Motor Imagery" & probe == "with obstacle"    ~ "Motor Imagery Obs")
  )

#print(EPE150.subj.contrast.A061)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
EPE150.subj.contrast_pooled.A061 <- EPE150.subj.EMM.A061 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(EPE150.subj.contrast_pooled.A061)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
EPE150.subj.diffContrast.A061 <- EPE150.subj.contrast.A061 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus MI NoObs",
           probe == "with obstacle"  ~ "EX minus MI Obs")
  )

#print(EPE150.subj.diffContrast.A061)


EPE150.subj.diffContrast_pooled.A061 <- EPE150.subj.contrast_pooled.A061 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(EPE150.subj.diffContrast_pooled.A061)






#  ==== ACCURACY: PLOTTING CONTRASTS OVERALL EXECUTION VS MOTOR IMAGERY (SUPPLEMENTARY FIGURE S5E) ====
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- "EX minus MI"
g.EPE150_EXvsMI.Contrast.A061 <- EPE150.Contrasts_long.A061 %>% filter(contrast=="EX vs MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX vs MI" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_mi1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A061, contrast == "EX vs MI"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.ExvsMI.subj.contrast.A061, contrast == "EX vs MI"),
             aes(y=0, x=diff,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.22,0.12),breaks=c(-0.2,-0.1,0,0.1,0.2)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_EXvsMI.Contrast.A061



#  ==== ACCURACY:: PLOTTING CONTRASTS  (SUPPLEMENTARY FIGURE S5GHKL) ====
# Probe without obstacle
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150_NoObs.Contrast.A061 <- EPE150.Contrasts_long.A061 %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "Motor Imagery NoObs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A061, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast.A061, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_NoObs.Contrast.A061


# Probe with obstacle
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150_Obs.Contrast.A061 <- EPE150.Contrasts_long.A061 %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "Motor Imagery Obs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A061, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast.A061, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_Obs.Contrast.A061



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~EPE~Difference(mm))
titleY <- ""
g.EPE150_NoObs_EXvsMI.Contrast.A061 <- EPE150.Contrasts_long.A061 %>% filter(contrast=="EX minus MI NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI NoObs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A061, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast.A061, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_NoObs_EXvsMI.Contrast.A061




# Probe with obstacle
titleX <- expression(Delta~EPE~Difference(mm))
titleY <- ""
g.EPE150_Obs_EXvsMI.Contrast.A061 <- EPE150.Contrasts_long.A061 %>% filter(contrast=="EX minus MI Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI Obs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A061, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast.A061, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_Obs_EXvsMI.Contrast.A061

















#==== A061: PRECISION: END POINT ERROR 150ms PROBE ====
# ==== PRECISION: END POINT ERROR 150ms PROBE DATA PREPARATION ====
# Endpoint error cannot be larger 0.5cm, also we only take values if velocity is below 5cm/s
dClean_EPP150 <- dClean %>%
  filter(
    !is.na(x_pos_HandTargetReached_150_Probe),
    x_pos_HandTargetReached_150_Probe >= -0.5 & 
      x_pos_HandTargetReached_150_Probe <= 0.5,
    y_pos_HandTargetReached_150_Probe >= 7.5 &
      y_pos_HandTargetReached_150_Probe <= 8.5,
    Vel_HandTargetReached_150_Probe<=5)

summary(dClean_EPP150$x_pos_HandTargetReached_150_Probe)
summary(dClean_EPP150$y_pos_HandTargetReached_150_Probe)
summary(dClean_EPP150$Vel_HandTargetReached_150_Probe)

# Compute condition-wise centers (convert from cm to mm)
centers <- dClean_EPP150 %>%
  group_by(stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  summarise(
    mean_x = mean(x_pos_HandTargetReached_150_Probe*10, na.rm = TRUE),
    mean_y = mean(y_pos_HandTargetReached_150_Probe*10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns


# Plot
dClean_EPP150 %>% 
  ggplot(aes(
    x = x_pos_HandTargetReached_150_Probe*10, #(convert from cm to mm)
    y = y_pos_HandTargetReached_150_Probe*10  #(convert from cm to mm)
  )) +
  geom_point(alpha = 0.6, size = 2) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_point(
    data = centers,
    aes(x = mean_x, y = mean_y),
    color = "red",
    size = 4
  ) +
  facet_grid(
    rows = vars(stop_signal_prime),
    cols = vars(obstacle_prime, obstacle_probe)
  ) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "X Position",
    y = "Y Position",
    title = "Reach Endpoint Positions (150 ms Post Target Hit)",
    subtitle = "Red dots = condition means (ellipse centers)"
  )


#calculate are for 95% ellipses of reach end points
chi2_95 <- qchisq(0.95, df = 2)

data_ellipse_area <- dClean_EPP150 %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  summarise(
    # covariance matrix of x and y
    cov_mat = list(cov(cbind(
      x_pos_HandTargetReached_150_Probe*10,
      y_pos_HandTargetReached_150_Probe*10
    ))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # eigenvalues of covariance matrix
    eig = list(eigen(cov_mat)$values),
    lambda1 = eig[[1]],
    lambda2 = eig[[2]],
    
    # ellipse semi-axes
    a = sqrt(lambda1 * chi2_95),
    b = sqrt(lambda2 * chi2_95),
    
    # area of 95% ellipse
    ellipse_area = pi * a * b
  ) %>%
  select(
    subID, stop_signal_prime, obstacle_prime, obstacle_probe,
    ellipse_area
  ) %>%
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns















# ==== MODEL FITTING: PRECISION: 3-FACTOR INTERACTION MODEL WITH TRUNCATED GAUSSIAN DISTRIBUTION ====
data_ellipse_area %>%
  ggplot( aes(x = ellipse_area, fill = prime) ) +
  facet_wrap(~probe:trial_type) +
  geom_histogram( alpha = 0.3, binwidth = 5, boundary = 0, position = 'identity' ) +
  labs( title = "A061- Histogram", subtitle = "End Point Precision 150ms after Target Hit Probe" ) +
  scale_x_continuous( name = 'Area Ellipse [mm²]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(data_ellipse_area$ellipse_area)

# looking at prior values
get_prior(bf(ellipse_area | trunc(lb = 0) ~ prime * probe * trial_type + ( prime * probe * trial_type | subID )),
          family = gaussian(),
          data   = data_ellipse_area)

# define priors
prior_EPP150    <- c(prior(normal(        70,     5   ),  class = "Intercept"), 
                     prior(normal(      0,     5   ),  class = "b"),
                     prior(student_t(3, 0,     5),  class = "sd"),
                     prior(student_t(3, 0,     5),  class = "sigma"),
                     prior(lkj(2),                     class = "cor"))




# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about x min, introducing truncation in prior predicitve checks causes NA's
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.EPP150_Probe_OPRIxOPROxSS <- brm( 
  bf(ellipse_area | trunc(lb = 0)  ~ prime * probe * trial_type + ( prime * probe * trial_type | subID )),     # model specification
  data   = data_ellipse_area,         # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_EPP150,              # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.EPP150_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )






#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(prime       = levels(data_ellipse_area$prime),
                        probe       = levels(data_ellipse_area$probe),
                        trial_type  = levels(data_ellipse_area$trial_type))


# Posterior draws of expected values (population-level, no random effects)
EPP150_Probe.posteriors <- as.data.frame(fitted(
  A061_fit.EPP150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(prime, probe, trial_type, sep = "_", drop = TRUE))
colnames(EPP150_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
EPP150_Probe.posteriors_long <- EPP150_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("prime", "probe", "trial_type"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
EPP150.EMM <- EPP150_Probe.posteriors_long %>%
  group_by(prime, probe, trial_type) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
EPP150.EMM



# Summaries: median and 95% HDI for Execution and Motor Imagery
EPP150_Probe.posteriors_long %>%
  group_by(trial_type) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
data_ellipse_area$subID    <- droplevels(data_ellipse_area$subID)
exp.cond.subj   <- expand.grid(subID       = levels(data_ellipse_area$subID),
                               prime       = levels(data_ellipse_area$prime),
                               probe       = levels(data_ellipse_area$probe),
                               trial_type  = levels(data_ellipse_area$trial_type))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.EPP150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
EPP150.subj.EMM <- draws_subj_df |>
  group_by(subID, prime, probe, trial_type) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)





# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S4AB) ====
# Prepare and tidy data
# rename factors
EPP150.EMM.A061 <- EPP150.EMM  %>% 
  mutate(
    trial_type  = factor(trial_type),
    prime       = factor(prime),
    probe       = factor(probe),
  )   
EPP150.EMM.A061

# rename factors
EPP150.subj.EMM.A061 <- EPP150.subj.EMM  %>% 
  mutate(
    trial_type  = factor(trial_type),
    prime       = factor(prime),
    probe       = factor(probe),
  ) 



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.EPP150_NoObs.EMM.A061 <-
  EPP150.subj.EMM.A061 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPP150.EMM.A061, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPP150.EMM.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(40, 100),
    breaks = seq(0, 100, 10)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.EPP150_NoObs.EMM.A061



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Area Ellipse (mm²)"

g.EPP150_Obs.EMM.A061 <-
  EPP150.subj.EMM.A061 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPP150.EMM.A061, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPP150.EMM.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(40, 100),
    breaks = seq(0, 100, 10)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

g.EPP150_Obs.EMM.A061


# ==== PRECISION: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
EPP150.Contrasts.A061           <- data.frame(matrix(ncol = 10, nrow = 40000))
colnames(EPP150.Contrasts.A061) <- c("Execution", "Motor Imagery", "EX minus MI", 
                                "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                                "EX minus MI NoObs", "EX minus MI Obs",
                                "EX vs MI")

# FOR POOLED DATA: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
# Execution trials (calculated as different minus same movement context)
EPP150.Contrasts.A061$`Execution`           <- ((EPP150_Probe.posteriors$`with obstacle_without obstacle_Execution` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Execution`) + ( EPP150_Probe.posteriors$`without obstacle_with obstacle_Execution` - EPP150_Probe.posteriors$`with obstacle_with obstacle_Execution`) )/2
# Motor Imagery trials (calculated as different minus same movement context)
EPP150.Contrasts.A061$`Motor Imagery`       <- ((EPP150_Probe.posteriors$`with obstacle_without obstacle_Motor Imagery` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Motor Imagery`) + ( EPP150_Probe.posteriors$`without obstacle_with obstacle_Motor Imagery` - EPP150_Probe.posteriors$`with obstacle_with obstacle_Motor Imagery`) )/2
# Difference of Execution vs Motor Imagery
EPP150.Contrasts.A061$`EX minus MI`         <- (EPP150.Contrasts.A061$`Execution` - EPP150.Contrasts.A061$`Motor Imagery`)
# Execution trials Probe without Obstacle
EPP150.Contrasts.A061$`Execution NoObs`     <- (EPP150_Probe.posteriors$`with obstacle_without obstacle_Execution` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Execution`)
# Execution trials with Obstacle
EPP150.Contrasts.A061$`Execution Obs`       <- (EPP150_Probe.posteriors$`with obstacle_with obstacle_Execution` - EPP150_Probe.posteriors$`without obstacle_with obstacle_Execution`)
# Motor Imagery trials Probe without Obstacle
EPP150.Contrasts.A061$`Motor Imagery NoObs` <- (EPP150_Probe.posteriors$`with obstacle_without obstacle_Motor Imagery` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Motor Imagery`)
# Motor Imagery trials Probe with Obstacle
EPP150.Contrasts.A061$`Motor Imagery Obs`   <- (EPP150_Probe.posteriors$`with obstacle_with obstacle_Motor Imagery` - EPP150_Probe.posteriors$`without obstacle_with obstacle_Motor Imagery`)
# Execution vs Motor Imagery without Obstacle
EPP150.Contrasts.A061$`EX minus MI NoObs`   <- EPP150.Contrasts.A061$`Execution NoObs` - EPP150.Contrasts.A061$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
EPP150.Contrasts.A061$`EX minus MI Obs`     <- EPP150.Contrasts.A061$`Execution Obs` - EPP150.Contrasts.A061$`Motor Imagery Obs`
# Execution vs Motor Imagery (Overall difference)
EPP150.Contrasts.A061$`EX vs MI`            <- ((EPP150_Probe.posteriors$`with obstacle_with obstacle_Execution` + EPP150_Probe.posteriors$`with obstacle_without obstacle_Execution` + EPP150_Probe.posteriors$`without obstacle_with obstacle_Execution` + EPP150_Probe.posteriors$`without obstacle_without obstacle_Execution` ) /4 ) - ((EPP150_Probe.posteriors$`with obstacle_with obstacle_Motor Imagery` + EPP150_Probe.posteriors$`with obstacle_without obstacle_Motor Imagery` + EPP150_Probe.posteriors$`without obstacle_with obstacle_Motor Imagery` + EPP150_Probe.posteriors$`without obstacle_without obstacle_Motor Imagery` ) /4 )

EPP150.Contrasts_long.A061 <- pivot_longer(EPP150.Contrasts.A061, cols = everything(),
                                      names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c("EX vs MI",
                    "Execution Obs", "Motor Imagery Obs",
                    "Execution NoObs", "Motor Imagery NoObs", 
                    "EX minus MI Obs", "EX minus MI NoObs", 
                    "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
EPP150.Contrasts_long.A061$contrast <- factor(EPP150.Contrasts_long.A061$contrast, levels = contrast_order)

# Check
head(EPP150.Contrasts_long.A061)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
EPP150.contrasts.summary.A061 <-
  EPP150.Contrasts_long.A061 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

EPP150.contrasts.summary.A061$pd <- format(EPP150.contrasts.summary.A061$pd, nsmall = 1)  
print(EPP150.contrasts.summary.A061, n = Inf, width = Inf)


#### Calculate ROPE
( RR.A061       <- rope_range(A061_fit.EPP150_Probe_OPRIxOPROxSS) ) 


# Calculate percent in ROPE for contrasts
options(digits=3)
EPP150.contrast_in_ROPE.A061                       <- as.data.frame(EPP150.contrasts.summary.A061)
EPP150.contrast_in_ROPE.A061$lowerROPE             <- RR.A061[1]
EPP150.contrast_in_ROPE.A061$upperROPE             <- RR.A061[2]
EPP150.contrast_in_ROPE.A061$CI_range              <- EPP150.contrast_in_ROPE.A061$upper - EPP150.contrast_in_ROPE.A061$lower
EPP150.contrast_in_ROPE.A061$minUpper              <- EPP150.contrast_in_ROPE.A061 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
EPP150.contrast_in_ROPE.A061$maxLower              <- EPP150.contrast_in_ROPE.A061 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
EPP150.contrast_in_ROPE.A061$DiffminUppermaxLower  <- EPP150.contrast_in_ROPE.A061$minUpper  - EPP150.contrast_in_ROPE.A061$maxLower 
EPP150.contrast_in_ROPE.A061$Zeros                 <- rep(0,nrow(EPP150.contrast_in_ROPE.A061))
EPP150.contrast_in_ROPE.A061$Overlap               <- EPP150.contrast_in_ROPE.A061  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
EPP150.contrast_in_ROPE.A061$perc_in_ROPE          <- (EPP150.contrast_in_ROPE.A061$Overlap*100)/EPP150.contrast_in_ROPE.A061$CI_range
EPP150.contrast_in_ROPE.A061[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
EPP150.ExvsMI.subj.contrast.A061 <- EPP150.subj.EMM.A061 %>%
  group_by(subID, trial_type) %>%
  summarise( .value = mean(.value)) %>%
  tidyr::pivot_wider(
    names_from = trial_type,
    values_from = .value
  ) %>%
  mutate(diff = `Execution` - `Motor Imagery`,
         contrast = "EX vs MI"
  )

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
EPP150.subj.contrast.A061 <- EPP150.subj.EMM.A061 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "Motor Imagery" & probe == "without obstacle" ~ "Motor Imagery NoObs",
           trial_type == "Motor Imagery" & probe == "with obstacle"    ~ "Motor Imagery Obs")
  )

#print(EPP150.subj.contrast.A061)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
EPP150.subj.contrast_pooled.A061 <- EPP150.subj.EMM.A061 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(EPP150.subj.contrast_pooled.A061)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
EPP150.subj.diffContrast.A061 <- EPP150.subj.contrast.A061 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus MI NoObs",
           probe == "with obstacle"  ~ "EX minus MI Obs")
  )

#print(EPP150.subj.diffContrast.A061)


EPP150.subj.diffContrast_pooled.A061 <- EPP150.subj.contrast_pooled.A061 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(EPP150.subj.diffContrast_pooled.A061)











#  ==== PRECISION: PLOTTING CONTRASTS OVERALL EXECUTION VS MOTOR IMAGERY (SUPPLEMENTARY FIGURE S4E) ====
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- "EX minus MI"
g.EPP150_EXvsMI.Contrast.A061 <- EPP150.Contrasts_long.A061 %>% filter(contrast=="EX vs MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX vs MI" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_mi1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A061, contrast == "EX vs MI"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -0.05, ymax = 0.5, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.ExvsMI.subj.contrast.A061, contrast == "EX vs MI"),
             aes(y=0, x=diff,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-11,5),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.5)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_EXvsMI.Contrast.A061



#  ==== PRECISION: PLOTTING CONTRASTS  (SUPPLEMENTARY FIGURE S4GHKL) ====
# Probe without obstacle
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150_NoObs.Contrast.A061 <- EPP150.Contrasts_long.A061 %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "Motor Imagery NoObs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A061, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast.A061, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_NoObs.Contrast.A061


# Probe with obstacle
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150_Obs.Contrast.A061 <- EPP150.Contrasts_long.A061 %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "Motor Imagery Obs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A061, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast.A061, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_Obs.Contrast.A061



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~EA~Difference~(mm^2))
titleY <- ""
g.EPP150_NoObs_EXvsMI.Contrast.A061 <- EPP150.Contrasts_long.A061 %>% filter(contrast=="EX minus MI NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI NoObs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A061, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast.A061, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_NoObs_EXvsMI.Contrast.A061




# Probe with obstacle
titleX <- expression(Delta~EA~Difference~(mm^2))
titleY <- ""
g.EPP150_Obs_EXvsMI.Contrast.A061 <- EPP150.Contrasts_long.A061 %>% filter(contrast=="EX minus MI Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI Obs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A061, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A061[1], xmax = RR.A061[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast.A061, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_Obs_EXvsMI.Contrast.A061















########################################## TRIMMED DATA ##############################################
#==== A061: TRIM DATA SUCH THAT INITIAL REACH ANGLES ARE SIMILAR BETWEEN OBSTACLE PRIME AND NO OBSTACLE PRIME (SEPARATE FOR EX/MI AND PROBE OBS/NO OBS) ====
d.Go_OPRO.A061 <- dClean.A061 %>% filter(stop_signal_prime=="go" & obstacle_probe=="yes")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

d.Go_NOPRO.A061 <- dClean.A061 %>% filter(stop_signal_prime=="go" & obstacle_probe=="no")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

d.Stop_OPRO.A061 <- dClean.A061 %>% filter(stop_signal_prime=="stop" & obstacle_probe=="yes")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

d.Stop_NOPRO.A061 <- dClean.A061 %>% filter(stop_signal_prime=="stop" & obstacle_probe=="no")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

#==== TRIMMING FUNCTION ====
balance_subject_means <- function(df,
                                  value_col = "ReachDiff2_Probe",
                                  cond_col  = "obstacle_prime",
                                  tol = 0.05,
                                  min_n_per_cond = 8) {
  dat <- df
  
  # Helper to extract condition means safely (expects levels "yes" and "no")
  get_means <- function(d) {
    s <- d %>%
      group_by(.data[[cond_col]]) %>%
      summarize(m = mean(.data[[value_col]], na.rm = TRUE),
                n = n(),
                .groups = "drop")
    means <- setNames(rep(NA_real_, 2), c("no", "yes"))
    ns    <- setNames(rep(0L,       2), c("no", "yes"))
    if (nrow(s)) {
      means[s[[cond_col]]] <- s$m
      ns[s[[cond_col]]]    <- s$n
    }
    list(means = means, ns = ns)
  }
  
  # If fewer than 2 conditions present, return as-is
  levs_present <- unique(dat[[cond_col]])
  if (length(levs_present) < 2) return(dat)
  
  repeat {
    stats <- get_means(dat)
    m_yes <- stats$means["yes"]
    m_no  <- stats$means["no"]
    n_yes <- stats$ns["yes"]
    n_no  <- stats$ns["no"]
    
    # Stop if any condition is at/below min_n or if a mean is missing
    if (is.na(m_yes) || is.na(m_no)) break
    if (any(c(n_yes, n_no) <= min_n_per_cond)) break
    
    diff_abs <- abs(m_yes - m_no)
    if (diff_abs <= tol) break
    
    target <- as.numeric((m_yes + m_no) / 2)
    
    # Choose which condition to trim (the one farther from target)
    dev_yes <- abs(m_yes - target)
    dev_no  <- abs(m_no  - target)
    cond_to_trim <- if (dev_yes >= dev_no) "yes" else "no"
    
    # Ensure we won't violate the minimum count
    n_trim <- if (cond_to_trim == "yes") n_yes else n_no
    if (n_trim <= min_n_per_cond) break
    
    # Find the single row to drop: farthest from target within cond_to_trim
    idx_to_drop <- dat %>%
      mutate(.row_id_tmp = row_number()) %>%
      filter(.data[[cond_col]] == cond_to_trim) %>%
      mutate(.dist = abs(.data[[value_col]] - target)) %>%
      arrange(desc(.dist)) %>%
      slice(1) %>%
      pull(.row_id_tmp)
    
    # If nothing to drop (shouldn't happen), stop
    if (length(idx_to_drop) == 0) break
    
    # Drop it
    dat <- dat[-idx_to_drop, , drop = FALSE]
  }
  
  dat
}


#==== APPLY TRIMMING AND CHECK DATA ====
# Apply within each subject
d_balanced.Go_OPRO.A061 <- d.Go_OPRO.A061 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))

d_balanced.Go_NOPRO.A061 <- d.Go_NOPRO.A061 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))

d_balanced.Stop_OPRO.A061 <- d.Stop_OPRO.A061 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))

d_balanced.Stop_NOPRO.A061 <- d.Stop_NOPRO.A061 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))


dClean_trimmed.A061       <- rbind(d_balanced.Go_OPRO.A061,d_balanced.Go_NOPRO.A061,d_balanced.Stop_OPRO.A061,d_balanced.Stop_NOPRO.A061)

# number of excluded trials after trimming
# Reports
cat("Anzahl Trials gesamt: ", nrow(dClean.A061), "\n")
cat("Anzahl Trials nach Trimming: ", nrow(dClean_trimmed.A061), "\n")
cat("Number of eliminated trials: ", nrow(dClean.A061) - nrow(dClean_trimmed.A061), "\n")
cat("Prozent eliminiert: ", round(100 * (nrow(dClean.A061) - nrow(dClean_trimmed.A061)) / nrow(dClean.A061), 2), "%\n")

#==== TRIMMED DATA: FINAL REACH ERROR PROBE ====
# measured as angular difference btw cursor mov_onset target center and cursor mov_onset cursor at target hit (EndPointError_abs_Probe)
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe



# ==== MODEL FITTING: FINAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=4) GAUSSIAN DISTRIBUTION ====
#merge trimmed dataframes
dClean_trimmed.A061       <- rbind(d_balanced.Go_OPRO.A061,d_balanced.Go_NOPRO.A061,d_balanced.Stop_OPRO.A061,d_balanced.Stop_NOPRO.A061)
dClean_trimmed_noobs.A061 <- dClean_trimmed.A061 %>% filter(obstacle_probe=="no")

dClean_trimmed_noobs.A061 %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A061 - TRIMMED DATA", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_trimmed_noobs.A061$EndPointError_abs_Probe)
sd(log(dClean_trimmed_noobs.A061$EndPointError_abs_Probe))


get_prior(bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_trimmed_noobs.A061)


# define priors
# Intercept: median exp(0.3) = 1.35; 1-sigma range exp(0.3 ± 0.5) ≈ [0.8, 2.23]
# Fixed effects: sd = 0.2 ⇒ 1-sigma multiplicative effect for a unit change in a predictor ≈ exp(±0.2) = ×[0.82, 1.22]
# Group-level SD: Half-normal(0, 0.2) on SDs of random intercepts/slopes (log scale), Implied per-subject multiplicative spread ≈ exp(±SD) ≈ ×[0.85, 1.17]
# Residual log-SD: Half-normal(0, 0.2) on sigma (log scale)
# Correlations among random effects
prior_FRE_noobs.trim   <- c(prior(normal(  1.4,     0.5   ),  class = "Intercept"), 
                            prior(normal(    0,     0.3   ),  class = "b"),
                            prior(normal(0, 0.3),    class = "sd",    lb = 0),
                            prior(normal(0, 0.3),    class = "sigma", lb = 0),
                            prior(lkj(2), class = "cor"))


# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 1 hour
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED <- brm( 
  bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_trimmed_noobs.A061, # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_FRE_noobs.trim,      # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED'), #save model
)
endTime <- Sys.time()
( endTime- startTime )


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_trimmed_noobs.A061$obstacle_prime),
                        stop_signal_prime  = levels(dClean_trimmed_noobs.A061$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_NoObs.posteriors.TRIMMED <- as.data.frame(fitted(
  A061_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_NoObs.posteriors.TRIMMED) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_NoObs.posteriors.TRIMMED_long <- FRE_Probe_NoObs.posteriors.TRIMMED %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_NoObs.EMM.TRIMMED <- FRE_Probe_NoObs.posteriors.TRIMMED_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_NoObs.EMM.TRIMMED



# Conditions
dClean_trimmed_noobs.A061$subID <- droplevels(dClean_trimmed_noobs.A061$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_trimmed_noobs.A061$subID),
                                  obstacle_prime     = levels(dClean_trimmed_noobs.A061$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_trimmed_noobs.A061$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_NoObs.subj.EMM.TRIMMED <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)






# ==== MODEL FITTING: FINAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH SKEWED GAUSSIAN DISTRIBUTION ====
#merge trimmed dataframes
dClean_trimmed.A061       <- rbind(d_balanced.Go_OPRO.A061,d_balanced.Go_NOPRO.A061,d_balanced.Stop_OPRO.A061,d_balanced.Stop_NOPRO.A061)
dClean_trimmed_obs.A061   <- dClean_trimmed.A061 %>% filter(obstacle_probe=="yes" & EndPointError_abs_Probe<=4)


dClean_trimmed_obs.A061 %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A061 - Histogram", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_trimmed_obs.A061$EndPointError_abs_Probe)
sd(log(dClean_trimmed_obs.A061$EndPointError_abs_Probe))

# looking at prior values
get_prior(bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = skew_normal(),
          data   = dClean_trimmed_obs.A061)



# define priors
prior_FRE_obs.trim <- c(prior(normal(  2.2,     0.5   ),  class = "Intercept"), 
                        prior(normal(    0,     0.3   ),  class = "b"),
                        prior(normal(    0,     0.3   ),  class = "sd",    lb = 0),
                        prior(normal(    0,     2     ),  class = "alpha"),
                        prior(normal(    0,     0.5   ),  class = "sigma",lb = 0),
                        prior(lkj(2), class = "cor"))



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 30min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A061_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED <- brm( 
  bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_trimmed_obs.A061,    # data
  family = skew_normal(),             # distribution of the response variable
  prior  = prior_FRE_obs.trim,        # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A061_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED'), #save model
)
endTime <- Sys.time()
( endTime- startTime )




#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_trimmed_obs.A061$obstacle_prime),
                        stop_signal_prime  = levels(dClean_trimmed_obs.A061$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_Obs.posteriors.TRIMMED <- as.data.frame(fitted(
  A061_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_Obs.posteriors.TRIMMED) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_Obs.posteriors.TRIMMED_long <- FRE_Probe_Obs.posteriors.TRIMMED %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_Obs.EMM.TRIMMED <- FRE_Probe_Obs.posteriors.TRIMMED_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_Obs.EMM.TRIMMED



# Conditions
dClean_trimmed_obs.A061$subID <- droplevels(dClean_trimmed_obs.A061$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_trimmed_obs.A061$subID),
                                  obstacle_prime     = levels(dClean_trimmed_obs.A061$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_trimmed_obs.A061$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A061_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_Obs.subj.EMM.TRIMMED <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)







# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S3AB) ====
# Prepare and tidy data
FRE_NoObs.EMM.TRIMMED$obstacle_probe   <- as.factor("no") 
FRE_Obs.EMM.TRIMMED$obstacle_probe     <- as.factor("yes") 
FRE.EMM.TRIMMED <- rbind(FRE_NoObs.EMM.TRIMMED,FRE_Obs.EMM.TRIMMED)

FRE_NoObs.subj.EMM.TRIMMED$obstacle_probe   <- as.factor("no") 
FRE_Obs.subj.EMM.TRIMMED$obstacle_probe     <- as.factor("yes") 
FRE.subj.EMM.TRIMMED <- rbind(FRE_NoObs.subj.EMM.TRIMMED,FRE_Obs.subj.EMM.TRIMMED)

# rename factors
FRE.EMM.TRIMMED.A061 <- FRE.EMM.TRIMMED  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

FRE.EMM.TRIMMED.A061

# rename factors
FRE.subj.EMM.TRIMMED.A061 <- FRE.subj.EMM.TRIMMED  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.FRE_NoObs.EMM.TRIMMED.A061 <-
  FRE.subj.EMM.TRIMMED.A061 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.TRIMMED.A061, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.TRIMMED.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(0, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.FRE_NoObs.EMM.TRIMMED.A061



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Final Reach Error (°)"

g.FRE_Obs.EMM.TRIMMED.A061 <-
  FRE.subj.EMM.TRIMMED.A061 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.TRIMMED.A061, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.TRIMMED.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(1, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

g.FRE_Obs.EMM.TRIMMED.A061



# ==== FINAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
# IMPORTANT!!!!: FOR POOLED CONSTRASTS: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
FRE.Contrasts.TRIMMED.A061           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(FRE.Contrasts.TRIMMED.A061) <- c("Execution", "Motor Imagery", "EX minus MI", 
                                     "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                                     "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials
FRE.Contrasts.TRIMMED.A061$`Execution`           <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_go` + FRE_Probe_Obs.posteriors.TRIMMED$`no_go`)/2 - ( FRE_Probe_NoObs.posteriors.TRIMMED$`no_go` + FRE_Probe_Obs.posteriors.TRIMMED$`yes_go`)/2
# Motor Imagery trials
FRE.Contrasts.TRIMMED.A061$`Motor Imagery`       <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_stop` + FRE_Probe_Obs.posteriors.TRIMMED$`no_stop`)/2 - ( FRE_Probe_NoObs.posteriors.TRIMMED$`no_stop` + FRE_Probe_Obs.posteriors.TRIMMED$`yes_stop`)/2
# Execution vs Motor Imagery
FRE.Contrasts.TRIMMED.A061$`EX minus MI`         <- (FRE.Contrasts.TRIMMED.A061$`Execution` - FRE.Contrasts.TRIMMED.A061$`Motor Imagery`)
# Execution trials Probe without Obstacle
FRE.Contrasts.TRIMMED.A061$`Execution NoObs`     <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_go` - FRE_Probe_NoObs.posteriors.TRIMMED$`no_go`)
# Execution trials with Obstacle
FRE.Contrasts.TRIMMED.A061$`Execution Obs`       <- (FRE_Probe_Obs.posteriors.TRIMMED$`yes_go` - FRE_Probe_Obs.posteriors.TRIMMED$`no_go`)
# Motor Imagery trials Probe without Obstacle
FRE.Contrasts.TRIMMED.A061$`Motor Imagery NoObs` <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_stop` - FRE_Probe_NoObs.posteriors.TRIMMED$`no_stop`)
# Motor Imagery trials Probe with Obstacle
FRE.Contrasts.TRIMMED.A061$`Motor Imagery Obs`   <- (FRE_Probe_Obs.posteriors.TRIMMED$`yes_stop` - FRE_Probe_Obs.posteriors.TRIMMED$`no_stop`)
# Execution vs Motor Imagery without Obstacle
FRE.Contrasts.TRIMMED.A061$`EX minus MI NoObs`   <- FRE.Contrasts.TRIMMED.A061$`Execution NoObs` - FRE.Contrasts.TRIMMED.A061$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
FRE.Contrasts.TRIMMED.A061$`EX minus MI Obs`     <- FRE.Contrasts.TRIMMED.A061$`Execution Obs` - FRE.Contrasts.TRIMMED.A061$`Motor Imagery Obs`


FRE.Contrasts_long.TRIMMED.A061 <- pivot_longer(FRE.Contrasts.TRIMMED.A061, cols = everything(),
                                           names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "Motor Imagery Obs",
  "Execution NoObs", "Motor Imagery NoObs", 
  "EX minus MI Obs", "EX minus MI NoObs", 
  "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
FRE.Contrasts_long.TRIMMED.A061$contrast <- factor(FRE.Contrasts_long.TRIMMED.A061$contrast, levels = contrast_order)

# Check
head(FRE.Contrasts_long.TRIMMED.A061)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
FRE.contrasts.summary.TRIMMED.A061 <-
  FRE.Contrasts_long.TRIMMED.A061 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

FRE.contrasts.summary.TRIMMED.A061$pd <- format(FRE.contrasts.summary.TRIMMED.A061$pd, nsmall = 4)  
#print(FRE.contrasts.summary.TRIMMED.A061, n = Inf, width = Inf)


#### Calculate ROPE
( RR_NoObs_trim.A061 <- rope_range(A061_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED) )
( RR_Obs_trim.A061   <- rope_range(A061_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED) )
( RR_trim.A061       <- ( RR_NoObs_trim.A061 + RR_Obs_trim.A061)/2 )

# Calculate percent in ROPE for contrasts
options(digits=3)
FRE.contrast_in_ROPE.TRIMMED.A061                       <- as.data.frame(FRE.contrasts.summary.TRIMMED.A061)
FRE.contrast_in_ROPE.TRIMMED.A061$lowerROPE             <- NA
FRE.contrast_in_ROPE.TRIMMED.A061$lowerROPE[c(1, 2, 5)] <- RR_Obs_trim.A061[1]
FRE.contrast_in_ROPE.TRIMMED.A061$lowerROPE[c(3, 4, 6)] <- RR_NoObs_trim.A061[1]
FRE.contrast_in_ROPE.TRIMMED.A061$lowerROPE[c(7:9)]     <- RR_trim.A061[1]
FRE.contrast_in_ROPE.TRIMMED.A061$upperROPE             <- NA
FRE.contrast_in_ROPE.TRIMMED.A061$upperROPE[c(1, 2, 5)] <- RR_Obs_trim.A061[2]
FRE.contrast_in_ROPE.TRIMMED.A061$upperROPE[c(3, 4, 6)] <- RR_NoObs_trim.A061[2]
FRE.contrast_in_ROPE.TRIMMED.A061$upperROPE[c(7:9)]     <- RR_trim.A061[2]
FRE.contrast_in_ROPE.TRIMMED.A061$CI_range              <- FRE.contrast_in_ROPE.TRIMMED.A061$upper - FRE.contrast_in_ROPE.TRIMMED.A061$lower
FRE.contrast_in_ROPE.TRIMMED.A061$minUpper              <- FRE.contrast_in_ROPE.TRIMMED.A061 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
FRE.contrast_in_ROPE.TRIMMED.A061$maxLower              <-  FRE.contrast_in_ROPE.TRIMMED.A061 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
FRE.contrast_in_ROPE.TRIMMED.A061$DiffminUppermaxLower  <- FRE.contrast_in_ROPE.TRIMMED.A061$minUpper  - FRE.contrast_in_ROPE.TRIMMED.A061$maxLower 
FRE.contrast_in_ROPE.TRIMMED.A061$Zeros                 <- rep(0,nrow(FRE.contrast_in_ROPE.TRIMMED.A061))
FRE.contrast_in_ROPE.TRIMMED.A061$Overlap               <- FRE.contrast_in_ROPE.TRIMMED.A061  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
FRE.contrast_in_ROPE.TRIMMED.A061$perc_in_ROPE          <- (FRE.contrast_in_ROPE.TRIMMED.A061$Overlap*100)/FRE.contrast_in_ROPE.TRIMMED.A061$CI_range
FRE.contrast_in_ROPE.TRIMMED.A061[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
FRE.subj.contrast.TRIMMED.A061 <- FRE.subj.EMM.TRIMMED.A061 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "Motor Imagery" & probe == "without obstacle" ~ "Motor Imagery NoObs",
           trial_type == "Motor Imagery" & probe == "with obstacle"    ~ "Motor Imagery Obs")
  )

#print(FRE.subj.contrast.TRIMMED.A061)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
FRE.subj.contrast_pooled.TRIMMED.A061 <- FRE.subj.EMM.TRIMMED.A061 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(FRE.subj.contrast_pooled.TRIMMED.A061)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
FRE.subj.diffContrast.TRIMMED.A061 <- FRE.subj.contrast.TRIMMED.A061 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus MI NoObs",
           probe == "with obstacle"  ~ "EX minus MI Obs")
  )

#print(FRE.subj.diffContrast.TRIMMED.A061)


FRE.subj.diffContrast_pooled.TRIMMED.A061 <- FRE.subj.contrast_pooled.TRIMMED.A061 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(FRE.subj.diffContrast_pooled.TRIMMED.A061)











#  ==== FINAL REACH ERROR: PLOTTING CONTRASTS  (SUPPLEMENTARY FIGURE S3EFIJ) ====
# Probe without obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_NoObs.Contrast.TRIMMED.A061 <- FRE.Contrasts_long.TRIMMED.A061 %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "Motor Imagery NoObs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A061, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs_trim.A061[1], xmax = RR_NoObs_trim.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.TRIMMED.A061, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs.Contrast.TRIMMED.A061


# Probe with obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_Obs.Contrast.TRIMMED.A061 <- FRE.Contrasts_long.TRIMMED.A061 %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "Motor Imagery Obs" = "Motor Imagery"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A061, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs_trim.A061[1], xmax = RR_Obs_trim.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.TRIMMED.A061, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs.Contrast.TRIMMED.A061



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_NoObs_EXvsMI.Contrast.TRIMMED.A061 <- FRE.Contrasts_long.TRIMMED.A061 %>% filter(contrast=="EX minus MI NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI NoObs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A061, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs_trim.A061[1], xmax = RR_NoObs_trim.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.TRIMMED.A061, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs_EXvsMI.Contrast.TRIMMED.A061




# Probe with obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_Obs_EXvsMI.Contrast.TRIMMED.A061 <- FRE.Contrasts_long.TRIMMED.A061 %>% filter(contrast=="EX minus MI Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI Obs" = "EX minus MI"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_mi1)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A061, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs_trim.A061[1], xmax = RR_Obs_trim.A061[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.TRIMMED.A061, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs_EXvsMI.Contrast.TRIMMED.A061










######################################################################################################
#========================================= A054 =====================================================
######################################################################################################
#==== A054: IMPORT SINGLE TRIAL DATA  ====
# clear workspace
# rm(list = ls(all = TRUE))  # Generally avoid in shared scripts
# ===== Paths (portable) =====
basePath     <- "C:/Experiments/A054_Hand Path Priming Obstacle Avoidance"
dataPath     <- file.path(basePath, "2_data", "4_clean")
docPath      <- file.path(basePath, "2_data", "2_data_documentation")
figurePath   <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")

# ===== Load single-trial data =====
# Note: N = 31, 704 trials (incl. 64 practise trials) per subject = 21824 trials
# Expecting column 'subID' present
d <- import(file.path(dataPath, "A054_data.csv")) %>%
  as_tibble()

#recode S999 to S31 for demographics
for( n in 1:nrow(d))
{
  if (d$subID[[n]] == 999)
  {
    d$subID[[n]] <- 31
  }
}


# ===== Load demographics (sheet 2) =====
# Ensure there is a subject identifier column to join on.
# If sheet lacks an explicit subID column but rows correspond to subjects 1..32,
# create it via row_number().
d_SubInfo_raw <- import(file.path(docPath, "A054_SubjInfo.xlsx"), which = 2) %>%
  as_tibble()

if (!"subID" %in% names(d_SubInfo_raw)) {
  d_SubInfo <- d_SubInfo_raw %>%
    mutate(subID = row_number())
} else {
  d_SubInfo <- d_SubInfo_raw
}

# Keep only needed columns and standardize names
d_SubInfo <- d_SubInfo %>%
  rename(
    Age = Age,
    Gender = Gender,
    Handedness = Handedness,
    EHI = EHI
  ) %>%
  select(subID, Age, Gender, Handedness, EHI)


# ===== Merge demographics to trial data =====
d <- d %>%
  left_join(d_SubInfo, by = "subID")

#recode S31 back to S999 
for( n in 1:nrow(d))
{
  if (d$subID[[n]] == 31)
  {
    d$subID[[n]] <- 999
  }
}

# ===== Quick checks =====
glimpse(d)
summary(select(d, Age, EHI))

# Example sanity plots
# library(ggplot2)
ggplot(d, aes(subID, Age)) + geom_point(alpha = 0.2, position = position_jitter(width = 0.15))
ggplot(d, aes(subID, EHI)) + geom_point(alpha = 0.2, position = position_jitter(width = 0.15))



# 1) Remove training block and subject 28 BEFORE converting subID to factor
#    (or convert with character comparison)
d <- d %>%
  filter(blocks_thisN != 0)


# If subID is numeric here, remove by numeric id, otherwise coerce to character
if (is.numeric(d$subID)) {
  d <- filter(d, subID != 28)
} else {
  d <- filter(d, as.character(subID) != "28")
}

# 2) Recode to factor in a vectorized way
d <- d %>%
  mutate(
    across(
      c(stop_signal_prime, obstacle_prime, obstacle_probe,
        target_location_probe, target_location_prime,
        deg_targets_prime, deg_targets_probe,
        Error_stop_signal, Error_too_slow_prime, Error_obstacle_prime,
        Error_too_slow_probe, Error_obstacle_probe,
        Error_MovOnsetDetection_Prime, Error_MovOnsetDetection_Probe,
        Error_Any, Gender, Handedness),
      ~ factor(.)
    ),
    subID = factor(subID)  # after filtering
  )



# 4) Drop unused levels after filtering
d <- d %>%
  mutate(across(where(is.factor), fct_drop))


# 5) Quick checks
# Check subject ids remaining
print(levels(d$subID))

# Adjust signed Reach difference
# signed ReachDiff should be between -180 and 180 degrees
plot(d$ReachDiff_sig_Prime)
plot(d$ReachDiff2_sig_Prime)
plot(d$ReachDiff_sig_Probe)
plot(d$ReachDiff2_sig_Probe)

wrap180 <- function(x) ((x + 180) %% 360) - 180
abs180  <- function(x) abs(wrap180(x))   # minimal absolute angular deviation

d$ReachDiff_sig_Prime  <- wrap180(d$ReachDiff_sig_Prime)
d$ReachDiff2_sig_Prime <- wrap180(d$ReachDiff2_sig_Prime)
d$ReachDiff_sig_Probe  <- wrap180(d$ReachDiff_sig_Probe)
d$ReachDiff2_sig_Probe <- wrap180(d$ReachDiff2_sig_Probe)


# Calculate Final Reach Error
circ_diff_deg <- function(a, b) wrap180(a - b)
d$EndPointError_sig_Prime <- circ_diff_deg(
  d$V_HandMovOnsetTargetCenter_Prime_Ori,
  d$V_HandMovOnsetHandTargetReached_Prime_Ori
)
d$EndPointError_abs_Prime <- abs(d$EndPointError_sig_Prime)

d$EndPointError_sig_Probe <- circ_diff_deg(
  d$V_HandMovOnsetTargetCenter_Probe_Ori,
  d$V_HandMovOnsetHandTargetReached_Probe_Ori
)
d$EndPointError_abs_Probe <- abs(d$EndPointError_sig_Probe)

plot(d$EndPointError_sig_Prime);  plot(d$EndPointError_abs_Prime)
plot(d$EndPointError_sig_Probe);  plot(d$EndPointError_abs_Probe)


# TIME BETWEEN MOVEMENT ONSET AND CURSOR OUT OF START
# movement onset should not be later than cursor out of start => values should lie under diagonal
# Prime movements
plot((d$frameNum_prime_time_cursor_out_of_start)-d$frameNum_prime_target_onset,d$mov_onset_Prime)
abline(coef = c(0,1))
#Probe movements
plot(d$frameNum_probe_time_cursor_out_of_start,d$mov_onset_Probe)
abline(coef = c(0,1))

# difference between time cursor out of start and movement onset (should be positive and close to zero; converted to ms)
d$time_diff_mov_onset_cursor_out_of_start_Prime <- (( (d$frameNum_prime_time_cursor_out_of_start) - d$frameNum_prime_target_onset) - d$mov_onset_Prime) *1000/d$frameRate
plot(d$time_diff_mov_onset_cursor_out_of_start_Prime)
d$time_diff_mov_onset_cursor_out_of_start_Probe <-    (d$frameNum_probe_time_cursor_out_of_start    - d$mov_onset_Probe) *1000/d$frameRate
plot(d$time_diff_mov_onset_cursor_out_of_start_Probe)


#==== REMOVE ERROR TRIALS ====
d_raw <- d

# 1) Normalize error flags to logicals
# Works if columns are logical, numeric 0/1, or factor/character "0"/"1"
error_cols <- c("Error_stop_signal","Error_too_slow_prime","Error_obstacle_prime",
                "Error_too_slow_probe","Error_obstacle_probe",
                "Error_MovOnsetDetection_Prime","Error_MovOnsetDetection_Probe",
                "Error_Any")

to_logical_01 <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  # factors/characters: coerce to character then numeric
  xv <- suppressWarnings(as.numeric(as.character(x)))
  return(xv == 1)
}

d <- d %>%
  mutate(across(all_of(error_cols), to_logical_01))

# 2) Compute denominators
n_total <- nrow(d)
n_stop  <- sum(d$stop_signal_prime == "stop", na.rm = TRUE)
n_go    <- sum(d$stop_signal_prime == "go",  na.rm = TRUE)

# 3) Helper to print counts and percents
report <- function(n, denom, label_n, label_pct) {
  cat(paste0(label_n, ": ", n, "\n"))
  if (!is.na(denom) && denom > 0) {
    cat(paste0(label_pct, ": ", round(100 * n / denom, 1), "%\n\n"))
  } else {
    cat("(denominator is 0 or NA)\n\n")
  }
}

cat(paste0("Total number of trials: ", n_total, "\n\n"))

# 4) Reports
n_err_stop   <- sum(d$Error_stop_signal, na.rm = TRUE)
report(n_err_stop, n_stop,
       "Number of MI error (stop-signal) trials",
       "Percent of stop-signal trials with MI error")

n_err_slow_pr <- sum(d$Error_too_slow_prime, na.rm = TRUE)
report(n_err_slow_pr, n_go,
       "Number of Too Slow Prime trials",
       "Percent of go (prime) trials too slow")

n_err_obst_pr <- sum(d$Error_obstacle_prime, na.rm = TRUE)
report(n_err_obst_pr, n_total,
       "Number of Obstacle Hit Prime trials",
       "Percent of all trials (Obstacle Hit Prime)")

n_err_slow_pb <- sum(d$Error_too_slow_probe, na.rm = TRUE)
report(n_err_slow_pb, n_total,
       "Number of Too Slow Probe trials",
       "Percent of all trials (Too Slow Probe)")

n_err_obst_pb <- sum(d$Error_obstacle_probe, na.rm = TRUE)
report(n_err_obst_pb, n_total,
       "Number of Obstacle Hit Probe trials",
       "Percent of all trials (Obstacle Hit Probe)")

n_err_mo_pr <- sum(d$Error_MovOnsetDetection_Prime, na.rm = TRUE)
report(n_err_mo_pr, n_total,
       "Number of Error Mov Onset Prime trials",
       "Percent of all trials (Mov Onset Prime error)")

n_err_mo_pb <- sum(d$Error_MovOnsetDetection_Probe, na.rm = TRUE)
report(n_err_mo_pb, n_total,
       "Number of Error Mov Onset Probe trials",
       "Percent of all trials (Mov Onset Probe error)")

n_err_any <- sum(d$Error_Any, na.rm = TRUE)
report(n_err_any, n_total,
       "Number of All Error trials",
       "Percent of all trials with any error")

# 5) Remove error trials
d <- d %>% filter(!Error_Any)


# Keep a copy before non-smooth exclusion
d_nonsmooth <- d

# Ensure numeric columns for peak velocity counts
to_num <- function(x) suppressWarnings(as.numeric(x))
d <- d %>%
  mutate(
    num_peak_vel_Probe = to_num(num_peak_vel_Probe),
    num_peak_vel_Prime = to_num(num_peak_vel_Prime)
  )

# Filter: probe must have <= 2 peaks; prime must have <= 2 peaks or missing (NA/NaN)
d <- d %>%
  filter(
    num_peak_vel_Probe <= 2 | is.na(num_peak_vel_Probe),
    num_peak_vel_Prime <= 2 | is.na(num_peak_vel_Prime)
  )


n_before_nonsmooth <- nrow(d_nonsmooth)
n_after_nonsmooth  <- nrow(d)
n_excl_nonsmooth   <- n_before_nonsmooth - n_after_nonsmooth

cat("Number of trials before non-smooth exclude: ", n_before_nonsmooth, "\n")
cat("Number of trials after  non-smooth exclude: ", n_after_nonsmooth,  "\n")
cat("Number of trials excluded (non-smooth): ", n_excl_nonsmooth, "\n")
cat("Percent non-smooth (of total trials): ",
    round(100 * n_excl_nonsmooth / n_total, 1), "%\n\n")


# Total excluded relative to original d_raw
n_total_excl <- nrow(d_raw) - n_after_nonsmooth
cat("Total trials excluded since d_raw: ", n_total_excl, "\n")
cat("Percent of trials excluded since d_raw: ",
    round(100 * n_total_excl / nrow(d_raw), 1), "%\n")


#==== OUTLIER REMOVAL ====
to_num <- function(x) suppressWarnings(as.numeric(x))

# Make a cleaned numeric copy
d_num <- d %>%
  mutate(
    RT_Probe   = to_num(RT_Probe),
    MT_Probe   = to_num(MT_Probe),
    RT_Prime   = to_num(RT_Prime),
    MT_Prime   = to_num(MT_Prime)
  )

n_all <- nrow(d_num)

# Define keep rules
keep_probe <- with(d_num, RT_Probe > 100 & MT_Probe <= 1000)

# For prime-phase thresholds, assume:
# - apply thresholds only to go trials,
# - let stop trials pass regardless of RT_Prime/MT_Prime (which may be NA)
keep_prime <- with(d_num,
                   ifelse(stop_signal_prime == "go",
                          RT_Prime > 100 & MT_Prime <= 1000,
                          TRUE)
)

dClean <- d_num %>% filter(keep_probe & keep_prime)

# Reports
cat("Anzahl Trials gesamt: ", n_all, "\n")
cat("Anzahl Trials nach Cleaning: ", nrow(dClean), "\n")
cat("Number of outlier trials: ", n_all - nrow(dClean), "\n")
cat("Prozent eliminiert: ", round(100 * (n_all - nrow(dClean)) / n_all, 2), "%\n")

dClean.A054 <- dClean

#==== DEMOGRAPHICS OF FINAL SAMPLE ====
# One row per subject
subjects <- dClean %>%
  distinct(subID, Handedness, Age, Gender, EHI)

# Quick check: each subject appears once
stopifnot(!any(duplicated(subjects$subID)))

# Make sure Age/EHI are numeric (if read as character)
subjects <- subjects %>%
  mutate(
    Age = suppressWarnings(as.numeric(Age)),
    EHI = suppressWarnings(as.numeric(EHI)),
    Gender = as.character(Gender),
    Handedness = as.character(Handedness)
  )

# Overall demographics
demo_overall <- subjects %>%
  summarise(
    n_subjects = n(),
    female  = sum(Gender %in% c("f","F","female","Female"), na.rm = TRUE),
    male    = sum(Gender %in% c("m","M","male","Male"), na.rm = TRUE),
    diverse = sum(Gender %in% c("d","D","diverse","nonbinary","non-binary"), na.rm = TRUE),
    meanAge = mean(Age, na.rm = TRUE),
    sdAge   = sd(Age, na.rm = TRUE),
    minAge  = min(Age, na.rm = TRUE),
    maxAge  = max(Age, na.rm = TRUE),
    meanEHI = mean(EHI, na.rm = TRUE),
    sdEHI   = sd(EHI, na.rm = TRUE),
    minEHI  = min(EHI, na.rm = TRUE),
    maxEHI  = max(EHI, na.rm = TRUE)
  )

# By handedness
subj_by_hand <- subjects %>%
  group_by(Handedness) %>%
  summarise(
    n_subjects = n(),
    female  = sum(Gender %in% c("f","F","female","Female"), na.rm = TRUE),
    male    = sum(Gender %in% c("m","M","male","Male"), na.rm = TRUE),
    diverse = sum(Gender %in% c("d","D","diverse","nonbinary","non-binary"), na.rm = TRUE),
    meanAge = mean(Age, na.rm = TRUE),
    sdAge   = sd(Age, na.rm = TRUE),
    meanEHI = mean(EHI, na.rm = TRUE),
    sdEHI   = sd(EHI, na.rm = TRUE),
    .groups = "drop"
  )

# Show demographics
demo_overall
subj_by_hand






#==== A054: FINAL REACH ERROR AS A FUNCTION OF INITIAL REACH ERROR ====
q_all <- dClean %>%
  filter(stop_signal_prime %in% c("go", "stop")) %>%
  select(subID,
         stop_signal_prime,
         obstacle_prime,
         obstacle_probe,
         EndPointError_abs_Probe,
         ReachDiff2_Probe) %>%
  
  # drop NA RTs before binning
  filter(!is.na(ReachDiff2_Probe)) %>%
  # create subject-wise RT deciles within each condition
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  # optionally enforce a minimum number of trials per subject-condition
  # filter(dplyr::n() >= 20) %>%
  mutate(Mydeciles = ntile(ReachDiff2_Probe, 10)) %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe, Mydeciles) %>%
  summarise(
    IREBin    = median(ReachDiff2_Probe, na.rm = TRUE),
    FREEffect = median(EndPointError_abs_Probe, na.rm = TRUE),
    nTrials  = dplyr::n(),              # for diagnostics
    .groups  = "drop"
  ) %>%
  # add human-readable labels
  mutate(
    trial_type  = if_else(stop_signal_prime == "go", "Execution", "No Movement"),
    prime       = if_else(obstacle_prime == "yes", "with obstacle", "without obstacle"),
    probe       = if_else(obstacle_probe == "yes", "with obstacle", "without obstacle"),
  )


# group-level aggregation
q_all_agg <- q_all %>%
  group_by(trial_type, prime, probe, Mydeciles) %>%
  summarise(
    meanBin      = mean(IREBin, na.rm = TRUE),
    medianBin    = median(IREBin, na.rm = TRUE),
    meanEffect   = mean(FREEffect, na.rm = TRUE),
    medianEffect = median(FREEffect, na.rm = TRUE),
    IQREffect    = IQR(FREEffect, na.rm = TRUE),
    seEffect     = sd(FREEffect, na.rm = TRUE) / sqrt(dplyr::n()),
    n_subj       = dplyr::n(),  # number of subject-bins contributing
    .groups      = "drop"
  )




#==== A054: FINAL REACH ERROR AS A FUNCTION OF INITIAL REACH ERROR: PLOTTING (SUPPLEMENTARY FIGURE S2 CDGH) ====
colorValues <- c(color_exe1, color_Neut)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY <- ""

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_GoNoObs.A054 <- q_all_agg %>%
  dplyr::filter(trial_type == "Execution", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 5, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("without obstacle")

g.FRE_IRE_GoNoObs.A054


colorValues <- c(color_exe1, color_Neut)  # adjust if needed
fontSize    <- 9
titleX      <- ""
titleY      <- ""

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_GoObs.A054 <- q_all_agg %>%
  dplyr::filter(trial_type == "Execution", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 5, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("with obstacle")

g.FRE_IRE_GoObs.A054



colorValues <- c(color_None,color_Neut)
fontSize    <- 9
titleX      <- "Initial Reach Error Deciles"
titleY      <- ""

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_StopNoObs.A054 <- q_all_agg %>%
  dplyr::filter(trial_type == "No Movement", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("with obstacle")

g.FRE_IRE_StopNoObs.A054



colorValues <- c(color_None,color_Neut)
fontSize    <- 9
titleX      <- "Initial Reach Error Deciles"
titleY      <- ""

pos_dodge <- position_dodge(width = 0.5)

g.FRE_IRE_StopObs.A054 <- q_all_agg %>%
  dplyr::filter(trial_type == "No Movement", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  geom_line(aes(group = prime), linewidth = 1, show.legend = FALSE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(1.0, 2.5),
                     breaks = c(1.0,1.5,2.0,2.5)) +
  scale_x_continuous(name = titleX, limits = c(-5, 95),
                     breaks = c(0, 30, 60, 90)) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  ggtitle("with obstacle")

g.FRE_IRE_StopObs.A054













######################################################################################################
#========================================= A054: BAYESIAN REGRESSION MODELS ==========================
######################################################################################################
#==== A054: INITIAL REACH ERROR PROBE ====
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe




# ==== MODEL FITTING: INITIAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (UB=180) LOGNORMAL DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_noobs$ReachDiff2_Probe)
sd(log(dClean_noobs$ReachDiff2_Probe))
# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = lognormal(),
          data   = dClean_noobs)

# define priors
# Intercept: median exp(2.3) = 10; 1-sigma range exp(2.3 ± 0.5) ≈ [6.05, 16.4]
# Fixed effects: sd = 0.3 ⇒ 1-sigma multiplicative factor exp(±0.3) ≈ ×[0.74, 1.35]
# Group-level SD: exponential(rate = 2) ⇒ mean = 0.5, median = log(2)/2 ≈ 0.347
# Residual log-SD: exponential(rate = 0.8) ⇒ mean = 1/0.8 = 1.25 (matches sd(log X))
# Correlations among random effects
prior_IRE_Probe_NoObs <- c(prior(normal(    2.3,    0.5  ),     class = "Intercept"),  
                           prior(normal(       0,    0.3  ),     class = "b"),          
                           prior(exponential(2),                 class = "sd"),         # mean 0.5 
                           prior(exponential(0.8),               class = "sigma"),      
                           prior(lkj(2),                         class = "cor"))      


# we model the mu here, we include random effects for participants
# runs about 45 min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.IRE_Probe_NoObs_OPRIxSS <- brm( 
  bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_noobs,              # data
  family = lognormal(),               # distribution of the response variable
  prior  = prior_IRE_Probe_NoObs,     # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control   = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.IRE_Probe_NoObs_OPRIxSS') #save model
)
endTime <- Sys.time()
( endTime- startTime )



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A054_fit.IRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(IRE_Probe_NoObs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
IRE_Probe_NoObs.posteriors_long <- IRE_Probe_NoObs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_NoObs.EMM <- IRE_Probe_NoObs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_NoObs.EMM



# Conditions
dClean_noobs$subID <- droplevels(dClean_noobs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_noobs$subID),
                                  obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.IRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
IRE_NoObs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)




# ==== MODEL FITTING: INITIAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=180) GAUSSIAN DISTRIBUTION ====
dClean_obs <- dClean %>% filter(obstacle_probe=="yes")

dClean_obs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_obs$ReachDiff2_Probe)

# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_obs)

# define priors
# Intercept around the observed mean, with moderate width
# Fixed effects (main effects and interaction), mildly weakly informative
# Group-level SDs: with scale ~10 on each random effect under subID
# Residual SD: half-Student-t or half-normal with scale ~15
# Correlations among random effects
prior_IRE_Probe_Obs <-  c(prior(normal(     54,   10  ), class = "Intercept"), 
                          prior(normal(      0,   7.5 ), class = "b"),
                          prior(student_t(3, 0,  10   ), class = "sd"), 
                          prior(lkj(2)                 , class = "cor"),
                          prior(student_t(3, 0,  15   ), class = "sigma"))


# we model the mu here. We include random effects for participants
# runs about 1.2 hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.IRE_Probe_Obs_OPRIxSS <- brm( 
  bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_obs,                # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_IRE_Probe_Obs,       # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.IRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )












#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A054_fit.IRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(IRE_Probe_Obs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
IRE_Probe_Obs.posteriors_long <- IRE_Probe_Obs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_Obs.EMM <- IRE_Probe_Obs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_Obs.EMM



# Conditions
dClean_obs$subID <- droplevels(dClean_obs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_obs$subID),
                                  obstacle_prime     = levels(dClean_obs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_obs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.IRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
IRE_Obs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)









# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 2CD) ====
# Prepare and tidy data
IRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
IRE.EMM.A054 <- rbind(IRE_NoObs.EMM,IRE_Obs.EMM)
IRE.EMM.A054

IRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
IRE.subj.EMM.A054 <- rbind(IRE_NoObs.subj.EMM,IRE_Obs.subj.EMM)

# rename factors
IRE.EMM.A054 <- IRE.EMM.A054  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

IRE.EMM.A054

# rename factors
IRE.subj.EMM.A054 <- IRE.subj.EMM.A054  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.IRE_NoObs.EMM.A054 <-
  IRE.subj.EMM.A054 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(IRE.EMM.A054, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM.A054, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  

  scale_y_continuous(
    name = titleY,
    limits = c(0, 80),
    breaks = seq(0, 80, 15)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.IRE_NoObs.EMM.A054



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.IRE_Obs.EMM.A054 <-
  IRE.subj.EMM.A054 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(IRE.EMM.A054, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM.A054, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(0, 80),
    breaks = seq(0, 80, 15)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.IRE_Obs.EMM.A054




# ==== INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
IRE.Contrasts.A054           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(IRE.Contrasts.A054) <- c("Execution", "No Movement", "EX minus NM", 
                                  "Execution NoObs", "Execution Obs", "No Movement NoObs", "No Movement Obs",
                                  "EX minus NM NoObs", "EX minus NM Obs")

# Execution trials
IRE.Contrasts.A054$`Execution`           <- (IRE_Probe_NoObs.posteriors$`yes_go` + IRE_Probe_Obs.posteriors$`yes_go`)/2 - ( IRE_Probe_NoObs.posteriors$`no_go` + IRE_Probe_Obs.posteriors$`no_go`)/2
# No Movement trials
IRE.Contrasts.A054$`No Movement`        <- (IRE_Probe_NoObs.posteriors$`yes_stop` + IRE_Probe_Obs.posteriors$`yes_stop`)/2 - ( IRE_Probe_NoObs.posteriors$`no_stop` + IRE_Probe_Obs.posteriors$`no_stop`)/2
# Execution vs No Movement
IRE.Contrasts.A054$`EX minus NM`         <- (IRE.Contrasts.A054$`Execution` - IRE.Contrasts.A054$`No Movement`)
# Execution trials Probe without Obstacle
IRE.Contrasts.A054$`Execution NoObs`     <- (IRE_Probe_NoObs.posteriors$`yes_go` - IRE_Probe_NoObs.posteriors$`no_go`)
# Execution trials with Obstacle
IRE.Contrasts.A054$`Execution Obs`       <- (IRE_Probe_Obs.posteriors$`yes_go` - IRE_Probe_Obs.posteriors$`no_go`)
# No Movement trials Probe without Obstacle
IRE.Contrasts.A054$`No Movement NoObs` <- (IRE_Probe_NoObs.posteriors$`yes_stop` - IRE_Probe_NoObs.posteriors$`no_stop`)
# No Movement trials Probe with Obstacle
IRE.Contrasts.A054$`No Movement Obs`   <- (IRE_Probe_Obs.posteriors$`yes_stop` - IRE_Probe_Obs.posteriors$`no_stop`)
# Execution vs No Movement without Obstacle
IRE.Contrasts.A054$`EX minus NM NoObs`   <- IRE.Contrasts.A054$`Execution NoObs` - IRE.Contrasts.A054$`No Movement NoObs` 
# Execution vs No Movement trials with Obstacle
IRE.Contrasts.A054$`EX minus NM Obs`     <- IRE.Contrasts.A054$`Execution Obs` - IRE.Contrasts.A054$`No Movement Obs`


IRE.Contrasts_long.A054 <- pivot_longer(IRE.Contrasts.A054, cols = everything(),
                                        names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "No Movement Obs",
  "Execution NoObs", "No Movement NoObs", 
  "EX minus NM Obs", "EX minus NM NoObs", 
  "Execution", "No Movement", "EX minus NM"
)

# Convert 'contrast' to a factor with this order
IRE.Contrasts_long.A054$contrast <- factor(IRE.Contrasts_long.A054$contrast, levels = contrast_order)

# Check
head(IRE.Contrasts_long.A054)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
IRE.contrasts.summary.A054 <-
  IRE.Contrasts_long.A054 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

IRE.contrasts.summary.A054$pd <- format(IRE.contrasts.summary.A054$pd, nsmall = 4)  
#print(IRE.contrasts.summary.A054, n = Inf, width = Inf)


#### Calculate ROPE
### Calculate and back-transform ROPE, this is a bit tricky
## this does not take into account that values are on log scale so don't use it
## apparently, back-transformation is now implemented :)
#RR_NoObs <- rope_range(A054_fit.IRE_Probe_NoObs_OPRIxSS) 
# compute a ROPE on the response scale for a lognormal outcome
## calculate manually: # IMPORTANT:SD(x) != exp(SD(log(x)))
dClean$logReachDiff2_Probe <- log(dClean$ReachDiff2_Probe) #log transform raw data
dClean_noobs               <- dClean %>% filter(obstacle_probe=="no")
# geometric mean of raw values ≈ exp(mean(log(x)))
#the geometric mean multiplied by the log-standard deviation. This should approximate the "natural" standard deviation pretty well.
geom_mean <- exp(mean(dClean_noobs$logReachDiff2_Probe))
# SD of log-transformed data
sd_log <- sd(dClean_noobs$logReachDiff2_Probe)
# ROPE = 0.1 * "natural" SD
rope_value_RE_log <- 0.1 * geom_mean * sd_log
rope_value_RE_log

( RR_IRE.NoObs.A054 <- c(-rope_value_RE_log, rope_value_RE_log) )
( RR_IRE.Obs.A054   <- rope_range(A054_fit.IRE_Probe_Obs_OPRIxSS) )
( RR.IRE.A054       <- (RR_IRE.NoObs.A054 + RR_IRE.Obs.A054) /2 )



# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE.A054                       <- as.data.frame(IRE.contrasts.summary.A054)
IRE.contrast_in_ROPE.A054$lowerROPE             <- NA
IRE.contrast_in_ROPE.A054$lowerROPE[c(1, 2, 5)] <- RR_IRE.Obs.A054[1]
IRE.contrast_in_ROPE.A054$lowerROPE[c(3, 4, 6)] <- RR_IRE.NoObs.A054[1]
IRE.contrast_in_ROPE.A054$lowerROPE[c(7:9)]     <- RR.IRE.A054[1]
IRE.contrast_in_ROPE.A054$upperROPE             <- NA
IRE.contrast_in_ROPE.A054$upperROPE[c(1, 2, 5)] <- RR_IRE.Obs.A054[2]
IRE.contrast_in_ROPE.A054$upperROPE[c(3, 4, 6)] <- RR_IRE.NoObs.A054[2]
IRE.contrast_in_ROPE.A054$upperROPE[c(7:9)]     <- RR.IRE.A054[2]
IRE.contrast_in_ROPE.A054$CI_range              <- IRE.contrast_in_ROPE.A054$upper - IRE.contrast_in_ROPE.A054$lower
IRE.contrast_in_ROPE.A054$minUpper              <- IRE.contrast_in_ROPE.A054 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
IRE.contrast_in_ROPE.A054$maxLower              <- IRE.contrast_in_ROPE.A054 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
IRE.contrast_in_ROPE.A054$DiffminUppermaxLower  <- IRE.contrast_in_ROPE.A054$minUpper  - IRE.contrast_in_ROPE.A054$maxLower 
IRE.contrast_in_ROPE.A054$Zeros                 <- rep(0,nrow(IRE.contrast_in_ROPE.A054))
IRE.contrast_in_ROPE.A054$Overlap               <- IRE.contrast_in_ROPE.A054  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
IRE.contrast_in_ROPE.A054$perc_in_ROPE          <- (IRE.contrast_in_ROPE.A054$Overlap*100)/IRE.contrast_in_ROPE.A054$CI_range
IRE.contrast_in_ROPE.A054[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
IRE.subj.contrast.A054 <- IRE.subj.EMM.A054 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "No Movement" & probe == "without obstacle" ~ "No Movement NoObs",
           trial_type == "No Movement" & probe == "with obstacle"    ~ "No Movement Obs")
  )

print(IRE.subj.contrast.A054)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
IRE.subj.contrast_pooled.A054 <- IRE.subj.EMM.A054 %>%
  group_by(subID, trial_type, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution"  ~ "Execution",
           trial_type == "Execution"  ~ "Execution",
           trial_type == "No Movement" ~ "No Movement",
           trial_type == "No Movement" ~ "No Movement"))

print(IRE.subj.contrast_pooled.A054)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus No Movement (Prime with obstacle minus Prime without obstacle) 
IRE.subj.diffContrast.A054 <- IRE.subj.contrast.A054 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_NM = Execution - `No Movement`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus NM NoObs",
           probe == "with obstacle"  ~ "EX minus NM Obs")
  )

print(IRE.subj.diffContrast.A054)


IRE.subj.diffContrast_pooled.A054 <- IRE.subj.contrast_pooled.A054 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_NM = Execution - `No Movement`,
         contrast = "EX minus MI"
  )

print(IRE.subj.diffContrast_pooled.A054)
















#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS  (FIGURE 2GHKL) ====
# Probe without obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_NoObs.Contrast.A054 <- IRE.Contrasts_long.A054 %>% filter(contrast=="Execution NoObs" | contrast=="No Movement NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "No Movement NoObs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A054, contrast == "Execution NoObs" | contrast=="No Movement NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.NoObs.A054[1], xmax = RR_IRE.NoObs.A054[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast.A054, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-3,22),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs.Contrast.A054


# Probe with obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_Obs.Contrast.A054 <- IRE.Contrasts_long.A054 %>% filter(contrast=="Execution Obs" | contrast=="No Movement Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "No Movement Obs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",     # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,       # No border on the stripes
    trim = TRUE,            # Trim the ends of the distributions
    linewidth = 0,          # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A054, contrast == "Execution Obs" | contrast=="No Movement Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.Obs.A054[1], xmax = RR_IRE.Obs.A054[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast.A054, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-3,12),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") + #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs.Contrast.A054



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_NoObs_EXvsNM.Contrast.A054 <- IRE.Contrasts_long.A054 %>% filter(contrast=="EX minus NM NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM NoObs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A054, contrast == "EX minus NM NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.NoObs.A054[1], xmax = RR_IRE.NoObs.A054[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast.A054, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_NM,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-3,22),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs_EXvsNM.Contrast.A054




# Probe with obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_Obs_EXvsNM.Contrast.A054 <- IRE.Contrasts_long.A054 %>% filter(contrast=="EX minus NM Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM Obs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A054, contrast == "EX minus NM Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.Obs.A054[1], xmax = RR_IRE.Obs.A054[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast.A054, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_NM,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-3,12),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs_EXvsNM.Contrast.A054















#==== A054: FINAL REACH ERROR PROBE ====
# measured as angular difference btw cursor mov_onset target center and cursor mov_onset cursor at target hit (EndPointError_abs_Probe)
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe
# ==== MODEL FITTING: FINAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=4) GAUSSIAN DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_noobs$EndPointError_abs_Probe)
sd(log(dClean_noobs$EndPointError_abs_Probe))


get_prior(bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_noobs)


# define priors
# Intercept: median exp(0.3) = 1.35; 1-sigma range exp(0.3 ± 0.5) ≈ [0.8, 2.23]
# Fixed effects: sd = 0.2 ⇒ 1-sigma multiplicative effect for a unit change in a predictor ≈ exp(±0.2) = ×[0.82, 1.22]
# Group-level SD: Half-normal(0, 0.2) on SDs of random intercepts/slopes (log scale), Implied per-subject multiplicative spread ≈ exp(±SD) ≈ ×[0.85, 1.17]
# Residual log-SD: Half-normal(0, 0.2) on sigma (log scale)
# Correlations among random effects
prior_FRE          <- c(prior(normal(  1.4,     0.5   ),  class = "Intercept"), 
                        prior(normal(    0,     0.3   ),  class = "b"),
                        prior(normal(0, 0.3),    class = "sd",    lb = 0),
                        prior(normal(0, 0.3),    class = "sigma", lb = 0),
                        prior(lkj(2), class = "cor"))


# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 1 hour
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.FRE_Probe_NoObs_OPRIxSS <- brm( 
  bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_noobs,              # data
  family = gaussian(),               # distribution of the response variable
  prior  = prior_FRE,                 # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.FRE_Probe_NoObs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A054_fit.FRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_NoObs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_NoObs.posteriors_long <- FRE_Probe_NoObs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_NoObs.EMM <- FRE_Probe_NoObs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_NoObs.EMM



# Conditions
dClean_noobs$subID <- droplevels(dClean_noobs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_noobs$subID),
                                  obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.FRE_Probe_NoObs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_NoObs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)





# ==== MODEL FITTING: FINAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH SKEWED GAUSSIAN DISTRIBUTION ====
dClean_obs <- dClean %>% filter(obstacle_probe=="yes" & EndPointError_abs_Probe<=4)

dClean_obs %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_obs$EndPointError_abs_Probe)
sd(log(dClean_obs$EndPointError_abs_Probe))

# looking at prior values
get_prior(bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = skew_normal(),
          data   = dClean_obs)



# define priors
prior_FRE_obs <-      c(prior(normal(  2.2,     0.5   ),  class = "Intercept"), 
                        prior(normal(    0,     0.3   ),  class = "b"),
                        prior(normal(    0,     0.3   ),  class = "sd",    lb = 0),
                        prior(normal(    0,     2     ),  class = "alpha"),
                        prior(normal(    0,     0.5   ),  class = "sigma",lb = 0),
                        prior(lkj(2), class = "cor"))



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 30min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.FRE_Probe_Obs_OPRIxSS <- brm( 
  bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_obs,                # data
  family = skew_normal(),             # distribution of the response variable
  prior  = prior_FRE_obs,             # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.FRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )









#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A054_fit.FRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_Obs.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_Obs.posteriors_long <- FRE_Probe_Obs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_Obs.EMM <- FRE_Probe_Obs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_Obs.EMM



# Conditions
dClean_obs$subID <- droplevels(dClean_obs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_obs$subID),
                                  obstacle_prime     = levels(dClean_obs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_obs$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.FRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_Obs.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)






# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 4CD) ====
# Prepare and tidy data
FRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
FRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
FRE.EMM.A054 <- rbind(FRE_NoObs.EMM,FRE_Obs.EMM)

FRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
FRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
FRE.subj.EMM.A054 <- rbind(FRE_NoObs.subj.EMM,FRE_Obs.subj.EMM)

# rename factors
FRE.EMM.A054 <- FRE.EMM.A054  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

FRE.EMM.A054

# rename factors
FRE.subj.EMM.A054 <- FRE.subj.EMM.A054  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.FRE_NoObs.EMM.A054 <-
  FRE.subj.EMM.A054 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.A054, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.A054, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  

  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(0, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.FRE_NoObs.EMM.A054



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.FRE_Obs.EMM.A054 <-
  FRE.subj.EMM.A054 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.A054, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.A054, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  

  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(1, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.FRE_Obs.EMM.A054
















# ==== FINAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
# IMPORTANT!!!!: FOR POOLED CONSTRASTS: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
FRE.Contrasts.A054           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(FRE.Contrasts.A054) <- c("Execution", "No Movement", "EX minus NM", 
                                  "Execution NoObs", "Execution Obs", "No Movement NoObs", "No Movement Obs",
                                  "EX minus NM NoObs", "EX minus NM Obs")

# Execution trials
FRE.Contrasts.A054$`Execution`           <- (FRE_Probe_NoObs.posteriors$`yes_go` + FRE_Probe_Obs.posteriors$`no_go`)/2 - ( FRE_Probe_NoObs.posteriors$`no_go` + FRE_Probe_Obs.posteriors$`yes_go`)/2
# Motor Imagery trials
FRE.Contrasts.A054$`No Movement`       <- (FRE_Probe_NoObs.posteriors$`yes_stop` + FRE_Probe_Obs.posteriors$`no_stop`)/2 - ( FRE_Probe_NoObs.posteriors$`no_stop` + FRE_Probe_Obs.posteriors$`yes_stop`)/2
# Execution vs Motor Imagery
FRE.Contrasts.A054$`EX minus NM`         <- (FRE.Contrasts.A054$`Execution` - FRE.Contrasts.A054$`No Movement`)
# Execution trials Probe without Obstacle
FRE.Contrasts.A054$`Execution NoObs`     <- (FRE_Probe_NoObs.posteriors$`yes_go` - FRE_Probe_NoObs.posteriors$`no_go`)
# Execution trials with Obstacle
FRE.Contrasts.A054$`Execution Obs`       <- (FRE_Probe_Obs.posteriors$`yes_go` - FRE_Probe_Obs.posteriors$`no_go`)
# Motor Imagery trials Probe without Obstacle
FRE.Contrasts.A054$`No Movement NoObs` <- (FRE_Probe_NoObs.posteriors$`yes_stop` - FRE_Probe_NoObs.posteriors$`no_stop`)
# Motor Imagery trials Probe with Obstacle
FRE.Contrasts.A054$`No Movement Obs`   <- (FRE_Probe_Obs.posteriors$`yes_stop` - FRE_Probe_Obs.posteriors$`no_stop`)
# Execution vs Motor Imagery without Obstacle
FRE.Contrasts.A054$`EX minus NM NoObs`   <- FRE.Contrasts.A054$`Execution NoObs` - FRE.Contrasts.A054$`No Movement NoObs` 
# Execution vs Motor Imagery trials with Obstacle
FRE.Contrasts.A054$`EX minus NM Obs`     <- FRE.Contrasts.A054$`Execution Obs` - FRE.Contrasts.A054$`No Movement Obs`


FRE.Contrasts_long.A054 <- pivot_longer(FRE.Contrasts.A054, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "No Movement Obs",
  "Execution NoObs", "No Movement NoObs", 
  "EX minus NM Obs", "EX minus NM NoObs", 
  "Execution", "No Movement", "EX minus NM"
)

# Convert 'contrast' to a factor with this order
FRE.Contrasts_long.A054$contrast <- factor(FRE.Contrasts_long.A054$contrast, levels = contrast_order)

# Check
head(FRE.Contrasts_long.A054)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
FRE.contrasts.summary.A054 <-
  FRE.Contrasts_long.A054 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

FRE.contrasts.summary.A054$pd <- format(FRE.contrasts.summary.A054$pd, nsmall = 4)  
#print(FRE.contrasts.summary.A054, n = Inf, width = Inf)


#### Calculate ROPE
( RR_FRE.NoObs.A054 <- rope_range(A054_fit.FRE_Probe_NoObs_OPRIxSS) )
( RR_FRE.Obs.A054   <- rope_range(A054_fit.FRE_Probe_Obs_OPRIxSS) )
( RR.FRE.A054       <- ( RR_FRE.NoObs.A054 + RR_FRE.Obs.A054)/2 )

# Calculate percent in ROPE for contrasts
options(digits=3)
FRE.contrast_in_ROPE.A054                       <- as.data.frame(FRE.contrasts.summary.A054)
FRE.contrast_in_ROPE.A054$lowerROPE             <- NA
FRE.contrast_in_ROPE.A054$lowerROPE[c(1, 2, 5)] <- RR_FRE.Obs.A054[1]
FRE.contrast_in_ROPE.A054$lowerROPE[c(3, 4, 6)] <- RR_FRE.NoObs.A054[1]
FRE.contrast_in_ROPE.A054$lowerROPE[c(7:9)]     <- RR.FRE.A054[1]
FRE.contrast_in_ROPE.A054$upperROPE             <- NA
FRE.contrast_in_ROPE.A054$upperROPE[c(1, 2, 5)] <- RR_FRE.Obs.A054[2]
FRE.contrast_in_ROPE.A054$upperROPE[c(3, 4, 6)] <- RR_FRE.NoObs.A054[2]
FRE.contrast_in_ROPE.A054$upperROPE[c(7:9)]     <- RR.FRE.A054[2]
FRE.contrast_in_ROPE.A054$CI_range              <- FRE.contrast_in_ROPE.A054$upper - FRE.contrast_in_ROPE.A054$lower
FRE.contrast_in_ROPE.A054$minUpper              <- FRE.contrast_in_ROPE.A054 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
FRE.contrast_in_ROPE.A054$maxLower              <-  FRE.contrast_in_ROPE.A054 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
FRE.contrast_in_ROPE.A054$DiffminUppermaxLower  <- FRE.contrast_in_ROPE.A054$minUpper  - FRE.contrast_in_ROPE.A054$maxLower 
FRE.contrast_in_ROPE.A054$Zeros                 <- rep(0,nrow(FRE.contrast_in_ROPE.A054))
FRE.contrast_in_ROPE.A054$Overlap               <- FRE.contrast_in_ROPE.A054  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
FRE.contrast_in_ROPE.A054$perc_in_ROPE          <- (FRE.contrast_in_ROPE.A054$Overlap*100)/FRE.contrast_in_ROPE.A054$CI_range
FRE.contrast_in_ROPE.A054[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
FRE.subj.contrast.A054 <- FRE.subj.EMM.A054 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "No Movement" & probe == "without obstacle" ~ "No Movement NoObs",
           trial_type == "No Movement" & probe == "with obstacle"    ~ "No Movement Obs")
  )

#print(FRE.subj.contrast.A054)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
FRE.subj.contrast_pooled.A054 <- FRE.subj.EMM.A054 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(FRE.subj.contrast_pooled.A054)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
FRE.subj.diffContrast.A054 <- FRE.subj.contrast.A054 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(go_vs_stop = Execution - `No Movement`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus NM NoObs",
           probe == "with obstacle"  ~ "EX minus NM Obs")
  )

#print(FRE.subj.diffContrast.A054)


FRE.subj.diffContrast_pooled.A054 <- FRE.subj.contrast_pooled.A054 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(go_vs_stop = Execution - `No Movement`,
         contrast = "EX minus NM"
  )

#print(FRE.subj.diffContrast_pooled.A054)









#  ==== FINAL REACH ERROR: PLOTTING CONTRASTS  (FIGURE 4GHKL) ====
# Probe without obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_NoObs.Contrast.A054 <- FRE.Contrasts_long.A054 %>% filter(contrast=="Execution NoObs" | contrast=="No Movement NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "No Movement NoObs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A054, contrast == "Execution NoObs" | contrast=="No Movement NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.NoObs.A054[1], xmax = RR_FRE.NoObs.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.A054, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs.Contrast.A054


# Probe with obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_Obs.Contrast.A054 <- FRE.Contrasts_long.A054 %>% filter(contrast=="Execution Obs" | contrast=="No Movement Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "No Movement Obs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A054, contrast == "Execution Obs" | contrast=="No Movement Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.Obs.A054[1], xmax = RR_FRE.Obs.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.A054, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs.Contrast.A054



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_NoObs_EXvsNM.Contrast.A054 <- FRE.Contrasts_long.A054 %>% filter(contrast=="EX minus NM NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM NoObs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A054, contrast == "EX minus NM NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.NoObs.A054[1], xmax = RR_FRE.NoObs.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.A054, probe == "without obstacle"),
             aes(y=0, x=go_vs_stop,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs_EXvsNM.Contrast.A054




# Probe with obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_Obs_EXvsNM.Contrast.A054 <- FRE.Contrasts_long.A054 %>% filter(contrast=="EX minus NM Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM Obs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.A054, contrast == "EX minus NM Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_FRE.Obs.A054[1], xmax = RR_FRE.Obs.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.A054, probe == "with obstacle"),
             aes(y=0, x=go_vs_stop,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs_EXvsNM.Contrast.A054












#==== ACCURACY: END POINT ERROR 150ms PROBE====
# measured as absolute distance of cursor 150ms after target hit to target center (EndPointError_150_Probe)
# ==== MODEL FITTING: ACCURACY: 3-FACTOR INTERACTION MODEL WITH TRUNCATED GAUSSIAN DISTRIBUTION ====
dClean_EPE150 <- dClean %>% filter(EndPointError_150_Probe <=0.5 & Vel_HandTargetReached_150_Probe<=5 )

dClean_EPE150 %>%
  ggplot( aes(x = EndPointError_150_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.02, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "End Point Error 150ms after Target Hit Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [cm]') +      
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_EPE150$EndPointError_150_Probe)

# looking at prior values
get_prior(bf(EndPointError_150_Probe | trunc(lb = 0, ub = 0.5) ~ obstacle_prime * obstacle_probe * stop_signal_prime + ( obstacle_prime * obstacle_probe * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_EPE150)

# define priors
prior_absEPE150    <- c(prior(normal(      0.3,   0.3   ),  class = "Intercept"), 
                        prior(normal(      0,     0.2   ),  class = "b"),
                        prior(student_t(3, 0,     0.2),  class = "sd"),
                        prior(student_t(3, 0,     0.2),  class = "sigma"),
                        prior(lkj(2),                     class = "cor"))



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 6 hours, introducing truncation in prior predicitve checks causes NA's
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.EPE150_Probe_OPRIxOPROxSS <- brm( 
  bf(EndPointError_150_Probe | trunc(lb = 0, ub = 0.5) ~ obstacle_prime * obstacle_probe * stop_signal_prime + ( obstacle_prime * obstacle_probe * stop_signal_prime | subID )),     # model specification
  data   = dClean_EPE150,             # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_absEPE150,           # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  init   = 0,
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 15),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/',  'A054_fit.EPE150_Probe_OPRIxOPROxSS') #save model
)
  
endTime <- Sys.time()
( endTime- startTime )











#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_EPE150$obstacle_prime),
                        obstacle_probe     = levels(dClean_EPE150$obstacle_probe),
                        stop_signal_prime  = levels(dClean_EPE150$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
EPE150_Probe.posteriors <- as.data.frame(fitted(
  A054_fit.EPE150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(obstacle_prime, obstacle_probe, stop_signal_prime, sep = "_", drop = TRUE))
colnames(EPE150_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
EPE150_Probe.posteriors_long <- EPE150_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("obstacle_prime", "obstacle_probe", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
EPE150.EMM <- EPE150_Probe.posteriors_long %>%
  group_by(obstacle_prime, obstacle_probe, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
EPE150.EMM



# Summaries: median and 95% HDI for Go and Stop
EPE150_Probe.posteriors_long %>%
  group_by(stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
dClean_EPE150$subID    <- droplevels(dClean_EPE150$subID)
exp.cond.subj   <- expand.grid(subID              = levels(dClean_EPE150$subID),
                               obstacle_prime     = levels(dClean_EPE150$obstacle_prime),
                               obstacle_probe     = levels(dClean_EPE150$obstacle_probe),
                               stop_signal_prime  = levels(dClean_EPE150$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.EPE150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
EPE150.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, obstacle_probe, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)



# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S5CD) ====
# Prepare and tidy data
# rename factors
EPE150.EMM.A054 <- EPE150.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    .value      = .value*10, .lower = .lower*10 , .upper = .upper*10  # convert from cm to mm
    # convert from cm to mm
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns
EPE150.EMM.A054

# rename factors
EPE150.subj.EMM.A054 <- EPE150.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    .value      = .value*10, .lower = .lower*10 , .upper = .upper*10  # convert from cm to mm
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.EPE150_NoObs.EMM.A054 <-
  EPE150.subj.EMM.A054 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPE150.EMM.A054, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPE150.EMM.A054, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(2, 3),
    breaks = seq(1, 5, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.EPE150_NoObs.EMM.A054



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.EPE150_Obs.EMM.A054 <-
  EPE150.subj.EMM.A054 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPE150.EMM.A054, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPE150.EMM.A054, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(2, 3),
    breaks = seq(1, 5, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.EPE150_Obs.EMM.A054















# ==== ACCURACY: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
EPE150.Contrasts.A054           <- data.frame(matrix(ncol = 10, nrow = 40000))
colnames(EPE150.Contrasts.A054) <- c("Execution", "No Movement", "EX minus NM", 
                                     "Execution NoObs", "Execution Obs", "No Movement NoObs", "No Movement Obs",
                                     "EX minus NM NoObs", "EX minus NM Obs",
                                     "EX vs NM")

# FOR POOLED DATA: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
# Execution trials (calculated as different minus same movement context)
EPE150.Contrasts.A054$`Execution`           <- ((EPE150_Probe.posteriors$`yes_no_go` - EPE150_Probe.posteriors$`no_no_go`) + ( EPE150_Probe.posteriors$`no_yes_go` - EPE150_Probe.posteriors$`yes_yes_go`) )/2
# No Movement trials (calculated as different minus same movement context)
EPE150.Contrasts.A054$`No Movement`       <- ((EPE150_Probe.posteriors$`yes_no_stop` - EPE150_Probe.posteriors$`no_no_stop`) + ( EPE150_Probe.posteriors$`no_yes_stop` - EPE150_Probe.posteriors$`yes_yes_stop`) )/2
# Difference of Execution vs No Movement
EPE150.Contrasts.A054$`EX minus NM`         <- (EPE150.Contrasts.A054$`Execution` - EPE150.Contrasts.A054$`No Movement`)
# Execution trials Probe without Obstacle
EPE150.Contrasts.A054$`Execution NoObs`     <- (EPE150_Probe.posteriors$`yes_no_go` - EPE150_Probe.posteriors$`no_no_go`)
# Execution trials with Obstacle
EPE150.Contrasts.A054$`Execution Obs`       <- (EPE150_Probe.posteriors$`yes_yes_go` - EPE150_Probe.posteriors$`no_yes_go`)
# No Movement trials Probe without Obstacle
EPE150.Contrasts.A054$`No Movement NoObs` <- (EPE150_Probe.posteriors$`yes_no_stop` - EPE150_Probe.posteriors$`no_no_stop`)
# No Movement trials Probe with Obstacle
EPE150.Contrasts.A054$`No Movement Obs`   <- (EPE150_Probe.posteriors$`yes_yes_stop` - EPE150_Probe.posteriors$`no_yes_stop`)
# Execution vs No Movement without Obstacle
EPE150.Contrasts.A054$`EX minus NM NoObs`   <- EPE150.Contrasts.A054$`Execution NoObs` - EPE150.Contrasts.A054$`No Movement NoObs` 
# Execution vs No Movement trials with Obstacle
EPE150.Contrasts.A054$`EX minus NM Obs`     <- EPE150.Contrasts.A054$`Execution Obs` - EPE150.Contrasts.A054$`No Movement Obs`
# Execution vs No Movement (Overall difference)
EPE150.Contrasts.A054$`EX vs NM`            <- ((EPE150_Probe.posteriors$`yes_yes_go` + EPE150_Probe.posteriors$`yes_no_go` + EPE150_Probe.posteriors$`no_yes_go` + EPE150_Probe.posteriors$`no_no_go` ) /4 ) - ((EPE150_Probe.posteriors$`yes_yes_stop` + EPE150_Probe.posteriors$`yes_no_stop` + EPE150_Probe.posteriors$`no_yes_stop` + EPE150_Probe.posteriors$`no_no_stop` ) /4 )

EPE150.Contrasts_long.A054 <- pivot_longer(EPE150.Contrasts.A054*10, cols = everything(), # convert from cm to mm
                                           names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c("EX vs NM",
                    "Execution Obs", "No Movement Obs",
                    "Execution NoObs", "No Movement NoObs", 
                    "EX minus NM Obs", "EX minus NM NoObs", 
                    "Execution", "No Movement", "EX minus NM"
)

# Convert 'contrast' to a factor with this order
EPE150.Contrasts_long.A054$contrast <- factor(EPE150.Contrasts_long.A054$contrast, levels = contrast_order)

# Check
head(EPE150.Contrasts_long.A054)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
EPE150.contrasts.summary.A054 <-
  EPE150.Contrasts_long.A054 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

EPE150.contrasts.summary.A054$pd <- format(EPE150.contrasts.summary.A054$pd, nsmall = 1)  
print(EPE150.contrasts.summary.A054, n = Inf, width = Inf)


#### Calculate ROPE
( RR.A054       <- rope_range(A054_fit.EPE150_Probe_OPRIxOPROxSS)*10 ) # convert from cm to mm


# Calculate percent in ROPE for contrasts
options(digits=3)
EPE150.contrast_in_ROPE.A054                       <- as.data.frame(EPE150.contrasts.summary.A054)
EPE150.contrast_in_ROPE.A054$lowerROPE             <- RR.A054[1]
EPE150.contrast_in_ROPE.A054$upperROPE             <- RR.A054[2]
EPE150.contrast_in_ROPE.A054$CI_range              <- EPE150.contrast_in_ROPE.A054$upper - EPE150.contrast_in_ROPE.A054$lower
EPE150.contrast_in_ROPE.A054$minUpper              <- EPE150.contrast_in_ROPE.A054 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
EPE150.contrast_in_ROPE.A054$maxLower              <- EPE150.contrast_in_ROPE.A054 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
EPE150.contrast_in_ROPE.A054$DiffminUppermaxLower  <- EPE150.contrast_in_ROPE.A054$minUpper  - EPE150.contrast_in_ROPE.A054$maxLower 
EPE150.contrast_in_ROPE.A054$Zeros                 <- rep(0,nrow(EPE150.contrast_in_ROPE.A054))
EPE150.contrast_in_ROPE.A054$Overlap               <- EPE150.contrast_in_ROPE.A054  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
EPE150.contrast_in_ROPE.A054$perc_in_ROPE          <- (EPE150.contrast_in_ROPE.A054$Overlap*100)/EPE150.contrast_in_ROPE.A054$CI_range
EPE150.contrast_in_ROPE.A054[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
EPE150.ExvsMI.subj.contrast.A054 <- EPE150.subj.EMM.A054 %>%
  group_by(subID, trial_type) %>%
  summarise( .value = mean(.value)) %>%
  tidyr::pivot_wider(
    names_from = trial_type,
    values_from = .value
  ) %>%
  mutate(diff = `Execution` - `No Movement`,
         contrast = "EX vs NM"
  )

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
EPE150.subj.contrast.A054 <- EPE150.subj.EMM.A054 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "No Movement" & probe == "without obstacle" ~ "No Movement NoObs",
           trial_type == "No Movement" & probe == "with obstacle"    ~ "No Movement Obs")
  )

#print(EPE150.subj.contrast.A054)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
EPE150.subj.contrast_pooled.A054 <- EPE150.subj.EMM.A054 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(EPE150.subj.contrast_pooled.A054)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
EPE150.subj.diffContrast.A054 <- EPE150.subj.contrast.A054 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_NM = Execution - `No Movement`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus NM NoObs",
           probe == "with obstacle"  ~ "EX minus NM Obs")
  )

#print(EPE150.subj.diffContrast.A054)


EPE150.subj.diffContrast_pooled.A054 <- EPE150.subj.contrast_pooled.A054 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_NM = Execution - `No Movement`,
         contrast = "EX minus NM"
  )

#print(EPE150.subj.diffContrast_pooled.A054)






#  ==== ACCURACY: PLOTTING CONTRASTS OVERALL EXECUTION VS NO MOVEMENT  (SUPPLEMENTARY FIGURE S5F) ====
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- "EX minus NM"
g.EPE150_EXvsNM.Contrast.A054 <- EPE150.Contrasts_long.A054 %>% filter(contrast=="EX vs NM") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX vs NM" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_None,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A054, contrast == "EX vs NM"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -2, ymax = 22, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.ExvsMI.subj.contrast.A054, contrast == "EX vs NM"),
             aes(y=0, x=diff,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.22,0.12),breaks=c(-0.2,-0.1,0,0.1,0.2)) + 
  scale_y_continuous(limits = c(-2, 22)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_EXvsNM.Contrast.A054



#  ==== ACCURACY:: PLOTTING CONTRASTS  (SUPPLEMENTARY FIGURE SIJMN) ====
# Probe without obstacle
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150_NoObs.Contrast.A054 <- EPE150.Contrasts_long.A054 %>% filter(contrast=="Execution NoObs" | contrast=="No Movement NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "No Movement NoObs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A054, contrast == "Execution NoObs" | contrast=="No Movement NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast.A054, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_NoObs.Contrast.A054


# Probe with obstacle
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150_Obs.Contrast.A054 <- EPE150.Contrasts_long.A054 %>% filter(contrast=="Execution Obs" | contrast=="No Movement Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "No Movement Obs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A054, contrast == "Execution Obs" | contrast=="No Movement Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast.A054, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_Obs.Contrast.A054



### Execution vs No Movement ###
# Probe without obstacle
titleX <- expression(Delta~EPE~Difference(mm))
titleY <- ""
g.EPE150_NoObs_EXvsNM.Contrast.A054 <- EPE150.Contrasts_long.A054 %>% filter(contrast=="EX minus NM NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM NoObs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A054, contrast == "EX minus NM NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast.A054, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_NM,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_NoObs_EXvsNM.Contrast.A054




# Probe with obstacle
titleX <- expression(Delta~EPE~Difference(mm))
titleY <- ""
g.EPE150_Obs_EXvsNM.Contrast.A054 <- EPE150.Contrasts_long.A054 %>% filter(contrast=="EX minus NM Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM Obs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary.A054, contrast == "EX minus NM Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast.A054, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_NM,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.25,0.3),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 12)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_Obs_EXvsNM.Contrast.A054









#==== PRECISION: END POINT ERROR 150ms PROBE ====
# ==== PRECISION: END POINT ERROR 150ms PROBE DATA PREPARATION ====
# Endpoint error cannot be larger 0.5cm, also we only take values if velocity is below 5cm/s
dClean_EPP150 <- dClean %>%
  filter(
    !is.na(x_pos_HandTargetReached_150_Probe),
    x_pos_HandTargetReached_150_Probe >= -0.5 & 
      x_pos_HandTargetReached_150_Probe <= 0.5,
    y_pos_HandTargetReached_150_Probe >= 7.5 &
      y_pos_HandTargetReached_150_Probe <= 8.5,
    Vel_HandTargetReached_150_Probe<=5)

summary(dClean_EPP150$x_pos_HandTargetReached_150_Probe)
summary(dClean_EPP150$y_pos_HandTargetReached_150_Probe)
summary(dClean_EPP150$Vel_HandTargetReached_150_Probe)

# Compute condition-wise centers (convert from cm to mm)
centers <- dClean_EPP150 %>%
  group_by(stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  summarise(
    mean_x = mean(x_pos_HandTargetReached_150_Probe*10, na.rm = TRUE),
    mean_y = mean(y_pos_HandTargetReached_150_Probe*10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "go", "stop")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns


# Plot
dClean_EPP150 %>% 
  ggplot(aes(
    x = x_pos_HandTargetReached_150_Probe*10, #(convert from cm to mm)
    y = y_pos_HandTargetReached_150_Probe*10  #(convert from cm to mm)
  )) +
  geom_point(alpha = 0.6, size = 2) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_point(
    data = centers,
    aes(x = mean_x, y = mean_y),
    color = "red",
    size = 4
  ) +
  facet_grid(
    rows = vars(stop_signal_prime),
    cols = vars(obstacle_prime, obstacle_probe)
  ) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "X Position",
    y = "Y Position",
    title = "Reach Endpoint Positions (150 ms Post Target Hit)",
    subtitle = "Red dots = condition means (ellipse centers)"
  )


#calculate are for 95% ellipses of reach end points
chi2_95 <- qchisq(0.95, df = 2)

data_ellipse_area <- dClean_EPP150 %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  summarise(
    # covariance matrix of x and y
    cov_mat = list(cov(cbind(
      x_pos_HandTargetReached_150_Probe*10,
      y_pos_HandTargetReached_150_Probe*10
    ))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # eigenvalues of covariance matrix
    eig = list(eigen(cov_mat)$values),
    lambda1 = eig[[1]],
    lambda2 = eig[[2]],
    
    # ellipse semi-axes
    a = sqrt(lambda1 * chi2_95),
    b = sqrt(lambda2 * chi2_95),
    
    # area of 95% ellipse
    ellipse_area = pi * a * b
  ) %>%
  select(
    subID, stop_signal_prime, obstacle_prime, obstacle_probe,
    ellipse_area
  ) %>%
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "go", "stop")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



# ==== MODEL FITTING: PRECISION: 3-FACTOR INTERACTION MODEL WITH TRUNCATED GAUSSIAN DISTRIBUTION ====
data_ellipse_area %>%
  ggplot( aes(x = ellipse_area, fill = prime) ) +
  facet_wrap(~probe:trial_type) +
  geom_histogram( alpha = 0.3, binwidth = 5, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "End Point Precision 150ms after Target Hit Probe" ) +
  scale_x_continuous( name = 'Area Ellipse [mm²]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(data_ellipse_area$ellipse_area)

# looking at prior values
get_prior(bf(ellipse_area | trunc(lb = 0) ~ prime * probe * trial_type + ( prime * probe * trial_type | subID )),
          family = gaussian(),
          data   = data_ellipse_area)

# define priors
prior_EPP150    <- c(prior(normal(        70,     5   ),  class = "Intercept"), 
                     prior(normal(      0,     5   ),  class = "b"),
                     prior(student_t(3, 0,     5),  class = "sd"),
                     prior(student_t(3, 0,     5),  class = "sigma"),
                     prior(lkj(2),                     class = "cor"))




# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about x min, introducing truncation in prior predicitve checks causes NA's
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.EPP150_Probe_OPRIxOPROxSS <- brm( 
  bf(ellipse_area | trunc(lb = 0)  ~ prime * probe * trial_type + ( prime * probe * trial_type | subID )),     # model specification
  data   = data_ellipse_area,         # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_EPP150,              # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.EPP150_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )








#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(prime       = levels(data_ellipse_area$prime),
                        probe       = levels(data_ellipse_area$probe),
                        trial_type  = levels(data_ellipse_area$trial_type))


# Posterior draws of expected values (population-level, no random effects)
EPP150_Probe.posteriors <- as.data.frame(fitted(
  A054_fit.EPP150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(prime, probe, trial_type, sep = "_", drop = TRUE))
colnames(EPP150_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
EPP150_Probe.posteriors_long <- EPP150_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("prime", "probe", "trial_type"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
EPP150.EMM <- EPP150_Probe.posteriors_long %>%
  group_by(prime, probe, trial_type) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
EPP150.EMM



# Summaries: median and 95% HDI for Execution and Motor Imagery
EPP150_Probe.posteriors_long %>%
  group_by(trial_type) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
data_ellipse_area$subID    <- droplevels(data_ellipse_area$subID)
exp.cond.subj   <- expand.grid(subID       = levels(data_ellipse_area$subID),
                               prime       = levels(data_ellipse_area$prime),
                               probe       = levels(data_ellipse_area$probe),
                               trial_type  = levels(data_ellipse_area$trial_type))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.EPP150_Probe_OPRIxOPROxSS,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
EPP150.subj.EMM <- draws_subj_df |>
  group_by(subID, prime, probe, trial_type) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)





# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S4CD) ====
# Prepare and tidy data
# rename factors
EPP150.EMM.A054 <- EPP150.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(trial_type == "go", "Execution", "No Movement")),
    prime       = factor(prime),
    probe       = factor(probe),
  )   
EPP150.EMM.A054

# rename factors
EPP150.subj.EMM.A054 <- EPP150.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(trial_type == "go", "Execution", "No Movement")),
    prime       = factor(prime),
    probe       = factor(probe),
  ) 



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.EPP150_NoObs.EMM.A054 <-
  EPP150.subj.EMM.A054 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPP150.EMM.A054, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPP150.EMM.A054, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(40, 100),
    breaks = seq(0, 100, 10)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.EPP150_NoObs.EMM.A054



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.EPP150_Obs.EMM.A054 <-
  EPP150.subj.EMM.A054 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(EPP150.EMM.A054, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPP150.EMM.A054, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(40, 100),
    breaks = seq(0, 100, 10)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.EPP150_Obs.EMM.A054


# ==== PRECISION: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
EPP150.Contrasts.A054           <- data.frame(matrix(ncol = 10, nrow = 40000))
colnames(EPP150.Contrasts.A054) <- c("Execution", "No Movement", "EX minus NM", 
                                     "Execution NoObs", "Execution Obs", "No Movement NoObs", "No Movement Obs",
                                     "EX minus NM NoObs", "EX minus NM Obs",
                                     "EX vs NM")

# FOR POOLED DATA: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
# Execution trials (calculated as different minus same movement context)
EPP150.Contrasts.A054$`Execution`           <- ((EPP150_Probe.posteriors$`with obstacle_without obstacle_go` - EPP150_Probe.posteriors$`without obstacle_without obstacle_go`) + ( EPP150_Probe.posteriors$`without obstacle_with obstacle_go` - EPP150_Probe.posteriors$`with obstacle_with obstacle_go`) )/2
# No Movement trials (calculated as different minus same movement context)
EPP150.Contrasts.A054$`No Movement`       <- ((EPP150_Probe.posteriors$`with obstacle_without obstacle_stop` - EPP150_Probe.posteriors$`without obstacle_without obstacle_stop`) + ( EPP150_Probe.posteriors$`without obstacle_with obstacle_stop` - EPP150_Probe.posteriors$`with obstacle_with obstacle_stop`) )/2
# Difference of Execution vs No Movement
EPP150.Contrasts.A054$`EX minus NM`         <- (EPP150.Contrasts.A054$`Execution` - EPP150.Contrasts.A054$`No Movement`)
# Execution trials Probe without Obstacle
EPP150.Contrasts.A054$`Execution NoObs`     <- (EPP150_Probe.posteriors$`with obstacle_without obstacle_go` - EPP150_Probe.posteriors$`without obstacle_without obstacle_go`)
# Execution trials with Obstacle
EPP150.Contrasts.A054$`Execution Obs`       <- (EPP150_Probe.posteriors$`with obstacle_with obstacle_go` - EPP150_Probe.posteriors$`without obstacle_with obstacle_go`)
# No Movement trials Probe without Obstacle
EPP150.Contrasts.A054$`No Movement NoObs` <- (EPP150_Probe.posteriors$`with obstacle_without obstacle_stop` - EPP150_Probe.posteriors$`without obstacle_without obstacle_stop`)
# No Movement trials Probe with Obstacle
EPP150.Contrasts.A054$`No Movement Obs`   <- (EPP150_Probe.posteriors$`with obstacle_with obstacle_stop` - EPP150_Probe.posteriors$`without obstacle_with obstacle_stop`)
# Execution vs No Movement without Obstacle
EPP150.Contrasts.A054$`EX minus NM NoObs`   <- EPP150.Contrasts.A054$`Execution NoObs` - EPP150.Contrasts.A054$`No Movement NoObs` 
# Execution vs No Movement trials with Obstacle
EPP150.Contrasts.A054$`EX minus NM Obs`     <- EPP150.Contrasts.A054$`Execution Obs` - EPP150.Contrasts.A054$`No Movement Obs`
# Execution vs No Movement (Overall difference)
EPP150.Contrasts.A054$`EX vs NM`            <- ((EPP150_Probe.posteriors$`with obstacle_with obstacle_go` + EPP150_Probe.posteriors$`with obstacle_without obstacle_go` + EPP150_Probe.posteriors$`without obstacle_with obstacle_go` + EPP150_Probe.posteriors$`without obstacle_without obstacle_go` ) /4 ) - ((EPP150_Probe.posteriors$`with obstacle_with obstacle_stop` + EPP150_Probe.posteriors$`with obstacle_without obstacle_stop` + EPP150_Probe.posteriors$`without obstacle_with obstacle_stop` + EPP150_Probe.posteriors$`without obstacle_without obstacle_stop` ) /4 )

EPP150.Contrasts_long.A054 <- pivot_longer(EPP150.Contrasts.A054, cols = everything(),
                                           names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c("EX vs NM",
                    "Execution Obs", "No Movement Obs",
                    "Execution NoObs", "No Movement NoObs", 
                    "EX minus NM Obs", "EX minus NM NoObs", 
                    "Execution", "No Movement", "EX minus NM"
)


# Convert 'contrast' to a factor with this order
EPP150.Contrasts_long.A054$contrast <- factor(EPP150.Contrasts_long.A054$contrast, levels = contrast_order)

# Check
head(EPP150.Contrasts_long.A054)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
EPP150.contrasts.summary.A054 <-
  EPP150.Contrasts_long.A054 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

EPP150.contrasts.summary.A054$pd <- format(EPP150.contrasts.summary.A054$pd, nsmall = 1)  
print(EPP150.contrasts.summary.A054, n = Inf, width = Inf)


#### Calculate ROPE
( RR.A054       <- rope_range(A054_fit.EPP150_Probe_OPRIxOPROxSS) ) 


# Calculate percent in ROPE for contrasts
options(digits=3)
EPP150.contrast_in_ROPE.A054                       <- as.data.frame(EPP150.contrasts.summary.A054)
EPP150.contrast_in_ROPE.A054$lowerROPE             <- RR.A054[1]
EPP150.contrast_in_ROPE.A054$upperROPE             <- RR.A054[2]
EPP150.contrast_in_ROPE.A054$CI_range              <- EPP150.contrast_in_ROPE.A054$upper - EPP150.contrast_in_ROPE.A054$lower
EPP150.contrast_in_ROPE.A054$minUpper              <- EPP150.contrast_in_ROPE.A054 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
EPP150.contrast_in_ROPE.A054$maxLower              <- EPP150.contrast_in_ROPE.A054 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
EPP150.contrast_in_ROPE.A054$DiffminUppermaxLower  <- EPP150.contrast_in_ROPE.A054$minUpper  - EPP150.contrast_in_ROPE.A054$maxLower 
EPP150.contrast_in_ROPE.A054$Zeros                 <- rep(0,nrow(EPP150.contrast_in_ROPE.A054))
EPP150.contrast_in_ROPE.A054$Overlap               <- EPP150.contrast_in_ROPE.A054  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
EPP150.contrast_in_ROPE.A054$perc_in_ROPE          <- (EPP150.contrast_in_ROPE.A054$Overlap*100)/EPP150.contrast_in_ROPE.A054$CI_range
EPP150.contrast_in_ROPE.A054[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
EPP150.ExvsNM.subj.contrast.A054 <- EPP150.subj.EMM.A054 %>%
  group_by(subID, trial_type) %>%
  summarise( .value = mean(.value)) %>%
  tidyr::pivot_wider(
    names_from = trial_type,
    values_from = .value
  ) %>%
  mutate(diff = `Execution` - `No Movement`,
         contrast = "EX vs NM"
  )

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
EPP150.subj.contrast.A054 <- EPP150.subj.EMM.A054 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "No Movement" & probe == "without obstacle" ~ "No Movement NoObs",
           trial_type == "No Movement" & probe == "with obstacle"    ~ "No Movement Obs")
  )

#print(EPP150.subj.contrast.A054)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
EPP150.subj.contrast_pooled.A054 <- EPP150.subj.EMM.A054 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(EPP150.subj.contrast_pooled.A054)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
EPP150.subj.diffContrast.A054 <- EPP150.subj.contrast.A054 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_NM = Execution - `No Movement`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus NM NoObs",
           probe == "with obstacle"  ~ "EX minus NM Obs")
  )

#print(EPP150.subj.diffContrast.A054)


EPP150.subj.diffContrast_pooled.A054 <- EPP150.subj.contrast_pooled.A054 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `No Movement`,
         contrast = "EX minus NM"
  )

#print(EPP150.subj.diffContrast_pooled.A054)











#  ==== PRECISION: PLOTTING CONTRASTS OVERALL EXECUTION VS MOTOR IMAGERY  (SUPPLEMENTARY FIGURE S4F) ====
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- "EX minus NM"
g.EPP150_EXvsNM.Contrast.A054 <- EPP150.Contrasts_long.A054 %>% filter(contrast=="EX vs NM") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX vs NM" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_None,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A054, contrast == "EX vs NM"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -0.05, ymax = 0.5, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.ExvsNM.subj.contrast.A054, contrast == "EX vs NM"),
             aes(y=0, x=diff,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-11,5),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.5)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_EXvsNM.Contrast.A054



#  ==== PRECISION: PLOTTING CONTRASTS (SUPPLEMENTARY FIGURE S4IJMN) ====
# Probe without obstacle
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150_NoObs.Contrast.A054 <- EPP150.Contrasts_long.A054 %>% filter(contrast=="Execution NoObs" | contrast=="No Movement NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "No Movement NoObs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A054, contrast == "Execution NoObs" | contrast=="No Movement NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast.A054, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_NoObs.Contrast.A054


# Probe with obstacle
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150_Obs.Contrast.A054 <- EPP150.Contrasts_long.A054 %>% filter(contrast=="Execution Obs" | contrast=="No Movement Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "No Movement Obs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A054, contrast == "Execution Obs" | contrast=="No Movement Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast.A054, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_Obs.Contrast.A054



### Execution vs No Movement ###
# Probe without obstacle
titleX <- expression(Delta~EA~Difference~(mm^2))
titleY <- ""
g.EPP150_NoObs_EXvsNM.Contrast.A054 <- EPP150.Contrasts_long.A054 %>% filter(contrast=="EX minus NM NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM NoObs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A054, contrast == "EX minus NM NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast.A054, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_NM,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_NoObs_EXvsNM.Contrast.A054




# Probe with obstacle
titleX <- expression(Delta~EA~Difference~(mm^2))
titleY <- ""
g.EPP150_Obs_EXvsNM.Contrast.A054 <- EPP150.Contrasts_long.A054 %>% filter(contrast=="EX minus NM Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM Obs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary.A054, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR.A054[1], xmax = RR.A054[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast.A054, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_NM,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-10,12),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_Obs_EXvsNM.Contrast.A054








########################################## TRIMMED DATA ##############################################
#==== A054: TRIM DATA SUCH THAT INITIAL REACH ANGLES ARE SIMILAR BETWEEN OBSTACLE PRIME AND NO OBSTACLE PRIME (SEPARATE FOR EX/MI AND PROBE OBS/NO OBS) ====
d.Go_OPRO.A054 <- dClean.A054 %>% filter(stop_signal_prime=="go" & obstacle_probe=="yes")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

d.Go_NOPRO.A054 <- dClean.A054 %>% filter(stop_signal_prime=="go" & obstacle_probe=="no")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

d.Stop_OPRO.A054 <- dClean.A054 %>% filter(stop_signal_prime=="stop" & obstacle_probe=="yes")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

d.Stop_NOPRO.A054 <- dClean.A054 %>% filter(stop_signal_prime=="stop" & obstacle_probe=="no")  %>% 
  select(c("subID", "obstacle_prime","obstacle_probe","stop_signal_prime","ReachDiff2_Probe","EndPointError_abs_Probe"))

#==== TRIMMING FUNCTION ====
balance_subject_means <- function(df,
                                  value_col = "ReachDiff2_Probe",
                                  cond_col  = "obstacle_prime",
                                  tol = 0.05,
                                  min_n_per_cond = 8) {
  dat <- df
  
  # Helper to extract condition means safely (expects levels "yes" and "no")
  get_means <- function(d) {
    s <- d %>%
      group_by(.data[[cond_col]]) %>%
      summarize(m = mean(.data[[value_col]], na.rm = TRUE),
                n = n(),
                .groups = "drop")
    means <- setNames(rep(NA_real_, 2), c("no", "yes"))
    ns    <- setNames(rep(0L,       2), c("no", "yes"))
    if (nrow(s)) {
      means[s[[cond_col]]] <- s$m
      ns[s[[cond_col]]]    <- s$n
    }
    list(means = means, ns = ns)
  }
  
  # If fewer than 2 conditions present, return as-is
  levs_present <- unique(dat[[cond_col]])
  if (length(levs_present) < 2) return(dat)
  
  repeat {
    stats <- get_means(dat)
    m_yes <- stats$means["yes"]
    m_no  <- stats$means["no"]
    n_yes <- stats$ns["yes"]
    n_no  <- stats$ns["no"]
    
    # Stop if any condition is at/below min_n or if a mean is missing
    if (is.na(m_yes) || is.na(m_no)) break
    if (any(c(n_yes, n_no) <= min_n_per_cond)) break
    
    diff_abs <- abs(m_yes - m_no)
    if (diff_abs <= tol) break
    
    target <- as.numeric((m_yes + m_no) / 2)
    
    # Choose which condition to trim (the one farther from target)
    dev_yes <- abs(m_yes - target)
    dev_no  <- abs(m_no  - target)
    cond_to_trim <- if (dev_yes >= dev_no) "yes" else "no"
    
    # Ensure we won't violate the minimum count
    n_trim <- if (cond_to_trim == "yes") n_yes else n_no
    if (n_trim <= min_n_per_cond) break
    
    # Find the single row to drop: farthest from target within cond_to_trim
    idx_to_drop <- dat %>%
      mutate(.row_id_tmp = row_number()) %>%
      filter(.data[[cond_col]] == cond_to_trim) %>%
      mutate(.dist = abs(.data[[value_col]] - target)) %>%
      arrange(desc(.dist)) %>%
      slice(1) %>%
      pull(.row_id_tmp)
    
    # If nothing to drop (shouldn't happen), stop
    if (length(idx_to_drop) == 0) break
    
    # Drop it
    dat <- dat[-idx_to_drop, , drop = FALSE]
  }
  
  dat
}


#==== APPLY TRIMMING AND CHECK DATA ====
d_balanced.Go_OPRO.A054 <- d.Go_OPRO.A054 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))

d_balanced.Go_NOPRO.A054 <- d.Go_NOPRO.A054 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))

d_balanced.Stop_OPRO.A054 <- d.Stop_OPRO.A054 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))

d_balanced.Stop_NOPRO.A054 <- d.Stop_NOPRO.A054 %>%
  group_by(subID) %>%
  group_split() %>%
  map_dfr(~ balance_subject_means(.x, tol = 0.05, min_n_per_cond = 40))


dClean_trimmed.A054       <- rbind(d_balanced.Go_OPRO.A054,d_balanced.Go_NOPRO.A054,d_balanced.Stop_OPRO.A054,d_balanced.Stop_NOPRO.A054)

# number of excluded trials after trimming
# Reports
cat("Anzahl Trials gesamt: ", nrow(dClean.A054), "\n")
cat("Anzahl Trials nach Trimming: ", nrow(dClean_trimmed.A054), "\n")
cat("Number of eliminated trials: ", nrow(dClean.A054) - nrow(dClean_trimmed.A054), "\n")
cat("Prozent eliminiert: ", round(100 * (nrow(dClean.A054) - nrow(dClean_trimmed.A054)) / nrow(dClean.A054), 2), "%\n")

#==== TRIMMED DATA: FINAL REACH ERROR PROBE ====
# measured as angular difference btw cursor mov_onset target center and cursor mov_onset cursor at target hit (EndPointError_abs_Probe)
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe

# ==== MODEL FITTING: FINAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=4) GAUSSIAN DISTRIBUTION ====
#merge trimmed dataframes
dClean_trimmed.A054       <- rbind(d_balanced.Go_OPRO.A054,d_balanced.Go_NOPRO.A054,d_balanced.Stop_OPRO.A054,d_balanced.Stop_NOPRO.A054)
dClean_trimmed_noobs.A054 <- dClean_trimmed.A054 %>% filter(obstacle_probe=="no")

dClean_trimmed_noobs.A054 %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - TRIMMED DATA", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_trimmed_noobs.A054$EndPointError_abs_Probe)
sd(log(dClean_trimmed_noobs.A054$EndPointError_abs_Probe))


get_prior(bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = gaussian(),
          data   = dClean_trimmed_noobs.A054)


# define priors
# Intercept: median exp(0.3) = 1.35; 1-sigma range exp(0.3 ± 0.5) ≈ [0.8, 2.23]
# Fixed effects: sd = 0.2 ⇒ 1-sigma multiplicative effect for a unit change in a predictor ≈ exp(±0.2) = ×[0.82, 1.22]
# Group-level SD: Half-normal(0, 0.2) on SDs of random intercepts/slopes (log scale), Implied per-subject multiplicative spread ≈ exp(±SD) ≈ ×[0.85, 1.17]
# Residual log-SD: Half-normal(0, 0.2) on sigma (log scale)
# Correlations among random effects
prior_FRE_noobs.trim   <- c(prior(normal(  1.3,     0.5   ),  class = "Intercept"), 
                            prior(normal(    0,     0.3   ),  class = "b"),
                            prior(normal(0, 0.3),    class = "sd",    lb = 0),
                            prior(normal(0, 0.3),    class = "sigma", lb = 0),
                            prior(lkj(2), class = "cor"))


# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 1 hour
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED <- brm( 
  bf(EndPointError_abs_Probe | trunc(lb = 0, ub = 4) ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_trimmed_noobs.A054, # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_FRE_noobs.trim,      # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED'), #save model
)
endTime <- Sys.time()
( endTime- startTime )





#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_trimmed_noobs.A054$obstacle_prime),
                        stop_signal_prime  = levels(dClean_trimmed_noobs.A054$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_NoObs.posteriors.TRIMMED <- as.data.frame(fitted(
  A054_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_NoObs.posteriors.TRIMMED) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_NoObs.posteriors.TRIMMED_long <- FRE_Probe_NoObs.posteriors.TRIMMED %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_NoObs.EMM.TRIMMED <- FRE_Probe_NoObs.posteriors.TRIMMED_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_NoObs.EMM.TRIMMED



# Conditions
dClean_trimmed_noobs.A054$subID <- droplevels(dClean_trimmed_noobs.A054$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_trimmed_noobs.A054$subID),
                                  obstacle_prime     = levels(dClean_trimmed_noobs.A054$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_trimmed_noobs.A054$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_NoObs.subj.EMM.TRIMMED <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)








# ==== MODEL FITTING: FINAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH SKEWED GAUSSIAN DISTRIBUTION ====
#merge trimmed dataframes
dClean_trimmed.A054       <- rbind(d_balanced.Go_OPRO.A054,d_balanced.Go_NOPRO.A054,d_balanced.Stop_OPRO.A054,d_balanced.Stop_NOPRO.A054)
dClean_trimmed_obs.A054   <- dClean_trimmed.A054 %>% filter(obstacle_probe=="yes" & EndPointError_abs_Probe<=4)


dClean_trimmed_obs.A054 %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A054- Histogram", subtitle = "Absolute End Point Error Probe" ) +
  scale_x_continuous( name = 'Absolute End Point Error [°]') +
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_trimmed_obs.A054$EndPointError_abs_Probe)
sd(log(dClean_trimmed_obs.A054$EndPointError_abs_Probe))

# looking at prior values
get_prior(bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = skew_normal(),
          data   = dClean_trimmed_obs.A054)



# define priors
prior_FRE_obs.trim <- c(prior(normal(  2.2,     0.5   ),  class = "Intercept"), 
                        prior(normal(    0,     0.3   ),  class = "b"),
                        prior(normal(    0,     0.3   ),  class = "sd",    lb = 0),
                        prior(normal(    0,     2     ),  class = "alpha"),
                        prior(normal(    0,     0.5   ),  class = "sigma",lb = 0),
                        prior(lkj(2), class = "cor"))



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 30min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED <- brm( 
  bf(EndPointError_abs_Probe ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean_trimmed_obs.A054,   # data
  family = skew_normal(),             # distribution of the response variable
  prior  = prior_FRE_obs.trim,        # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED'), #save model
)
endTime <- Sys.time()
( endTime- startTime )





#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_trimmed_obs.A054$obstacle_prime),
                        stop_signal_prime  = levels(dClean_trimmed_obs.A054$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_Obs.posteriors.TRIMMED <- as.data.frame(fitted(
  A054_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(FRE_Probe_Obs.posteriors.TRIMMED) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
FRE_Probe_Obs.posteriors.TRIMMED_long <- FRE_Probe_Obs.posteriors.TRIMMED %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
FRE_Obs.EMM.TRIMMED <- FRE_Probe_Obs.posteriors.TRIMMED_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
FRE_Obs.EMM.TRIMMED



# Conditions
dClean_trimmed_obs.A054$subID <- droplevels(dClean_trimmed_obs.A054$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_trimmed_obs.A054$subID),
                                  obstacle_prime     = levels(dClean_trimmed_obs.A054$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_trimmed_obs.A054$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A054_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED,
  newdata    = exp.cond.subj,
  re_formula = NULL,
  summary    = FALSE
)
# Long format with condition labels
draws_subj_df <- as.data.frame(draws_subj) |>
  mutate(.draw = row_number()) |>
  pivot_longer(cols = - .draw, names_to = "row", values_to = "value") |>
  mutate(row = as.integer(gsub("^V", "", row))) |>
  left_join(exp.cond.subj |> mutate(row = row_number()), by = "row")

# Median and 95% HDI per subject-condition
FRE_Obs.subj.EMM.TRIMMED <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)







# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S3CD) ====
# Prepare and tidy data
FRE_NoObs.EMM.TRIMMED$obstacle_probe   <- as.factor("no") 
FRE_Obs.EMM.TRIMMED$obstacle_probe     <- as.factor("yes") 
FRE.EMM.TRIMMED <- rbind(FRE_NoObs.EMM.TRIMMED,FRE_Obs.EMM.TRIMMED)

FRE_NoObs.subj.EMM.TRIMMED$obstacle_probe   <- as.factor("no") 
FRE_Obs.subj.EMM.TRIMMED$obstacle_probe     <- as.factor("yes") 
FRE.subj.EMM.TRIMMED <- rbind(FRE_NoObs.subj.EMM.TRIMMED,FRE_Obs.subj.EMM.TRIMMED)

# rename factors
FRE.EMM.TRIMMED.A054 <- FRE.EMM.TRIMMED  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

FRE.EMM.TRIMMED.A054

# rename factors
FRE.subj.EMM.TRIMMED.A054 <- FRE.subj.EMM.TRIMMED  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.FRE_NoObs.EMM.TRIMMED.A054 <-
  FRE.subj.EMM.TRIMMED.A054 %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.TRIMMED.A054, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.TRIMMED.A054, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(0, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.FRE_NoObs.EMM.TRIMMED.A054



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.FRE_Obs.EMM.TRIMMED.A054 <-
  FRE.subj.EMM.TRIMMED.A054 %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_wrap(~probe, strip.position = "top") + 
  
  # Subject-level points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = .8
    ),
    size = 1,
    alpha = .5,
    show.legend = FALSE
  ) +
  
  # Group-level HDI
  geom_errorbar(
    data = dplyr::filter(FRE.EMM.TRIMMED.A054, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM.TRIMMED.A054, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(1, 2.5),
    breaks = seq(1, 3, 0.5)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.FRE_Obs.EMM.TRIMMED.A054


# ==== FINAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
# IMPORTANT!!!!: FOR POOLED CONSTRASTS: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
FRE.Contrasts.TRIMMED.A054           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(FRE.Contrasts.TRIMMED.A054) <- c("Execution", "No Movement", "EX minus NM", 
                                          "Execution NoObs", "Execution Obs", "No Movement NoObs", "No Movement Obs",
                                          "EX minus NM NoObs", "EX minus NM Obs")

# Execution trials
FRE.Contrasts.TRIMMED.A054$`Execution`           <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_go` + FRE_Probe_Obs.posteriors.TRIMMED$`no_go`)/2 - ( FRE_Probe_NoObs.posteriors.TRIMMED$`no_go` + FRE_Probe_Obs.posteriors.TRIMMED$`yes_go`)/2
# Motor Imagery trials
FRE.Contrasts.TRIMMED.A054$`No Movement`       <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_stop` + FRE_Probe_Obs.posteriors.TRIMMED$`no_stop`)/2 - ( FRE_Probe_NoObs.posteriors.TRIMMED$`no_stop` + FRE_Probe_Obs.posteriors.TRIMMED$`yes_stop`)/2
# Execution vs Motor Imagery
FRE.Contrasts.TRIMMED.A054$`EX minus NM`         <- (FRE.Contrasts.TRIMMED.A054$`Execution` - FRE.Contrasts.TRIMMED.A054$`No Movement`)
# Execution trials Probe without Obstacle
FRE.Contrasts.TRIMMED.A054$`Execution NoObs`     <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_go` - FRE_Probe_NoObs.posteriors.TRIMMED$`no_go`)
# Execution trials with Obstacle
FRE.Contrasts.TRIMMED.A054$`Execution Obs`       <- (FRE_Probe_Obs.posteriors.TRIMMED$`yes_go` - FRE_Probe_Obs.posteriors.TRIMMED$`no_go`)
# Motor Imagery trials Probe without Obstacle
FRE.Contrasts.TRIMMED.A054$`No Movement NoObs` <- (FRE_Probe_NoObs.posteriors.TRIMMED$`yes_stop` - FRE_Probe_NoObs.posteriors.TRIMMED$`no_stop`)
# Motor Imagery trials Probe with Obstacle
FRE.Contrasts.TRIMMED.A054$`No Movement Obs`   <- (FRE_Probe_Obs.posteriors.TRIMMED$`yes_stop` - FRE_Probe_Obs.posteriors.TRIMMED$`no_stop`)
# Execution vs Motor Imagery without Obstacle
FRE.Contrasts.TRIMMED.A054$`EX minus NM NoObs`   <- FRE.Contrasts.TRIMMED.A054$`Execution NoObs` - FRE.Contrasts.TRIMMED.A054$`No Movement NoObs` 
# Execution vs Motor Imagery trials with Obstacle
FRE.Contrasts.TRIMMED.A054$`EX minus NM Obs`     <- FRE.Contrasts.TRIMMED.A054$`Execution Obs` - FRE.Contrasts.TRIMMED.A054$`No Movement Obs`


FRE.Contrasts_long.TRIMMED.A054 <- pivot_longer(FRE.Contrasts.TRIMMED.A054, cols = everything(),
                                           names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "No Movement Obs",
  "Execution NoObs", "No Movement NoObs", 
  "EX minus NM Obs", "EX minus NM NoObs", 
  "Execution", "No Movement", "EX minus NM"
)

# Convert 'contrast' to a factor with this order
FRE.Contrasts_long.TRIMMED.A054$contrast <- factor(FRE.Contrasts_long.TRIMMED.A054$contrast, levels = contrast_order)

# Check
head(FRE.Contrasts_long.TRIMMED.A054)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
FRE.contrasts.summary.TRIMMED.A054 <-
  FRE.Contrasts_long.TRIMMED.A054 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

FRE.contrasts.summary.TRIMMED.A054$pd <- format(FRE.contrasts.summary.TRIMMED.A054$pd, nsmall = 4)  
print(FRE.contrasts.summary.TRIMMED.A054, n = Inf, width = Inf)


#### Calculate ROPE
( RR_NoObs_trim.A054 <- rope_range(A054_fit.FRE_Probe_NoObs_OPRIxSS_TRIMMED) )
( RR_Obs_trim.A054   <- rope_range(A054_fit.FRE_Probe_Obs_OPRIxSS_TRIMMED) )
( RR_trim.A054       <- ( RR_NoObs_trim.A054 + RR_Obs_trim.A054)/2 )

# Calculate percent in ROPE for contrasts
options(digits=3)
FRE.contrast_in_ROPE.TRIMMED.A054                       <- as.data.frame(FRE.contrasts.summary.TRIMMED.A054)
FRE.contrast_in_ROPE.TRIMMED.A054$lowerROPE             <- NA
FRE.contrast_in_ROPE.TRIMMED.A054$lowerROPE[c(1, 2, 5)] <- RR_Obs_trim.A054[1]
FRE.contrast_in_ROPE.TRIMMED.A054$lowerROPE[c(3, 4, 6)] <- RR_NoObs_trim.A054[1]
FRE.contrast_in_ROPE.TRIMMED.A054$lowerROPE[c(7:9)]     <- RR_trim.A054[1]
FRE.contrast_in_ROPE.TRIMMED.A054$upperROPE             <- NA
FRE.contrast_in_ROPE.TRIMMED.A054$upperROPE[c(1, 2, 5)] <- RR_Obs_trim.A054[2]
FRE.contrast_in_ROPE.TRIMMED.A054$upperROPE[c(3, 4, 6)] <- RR_NoObs_trim.A054[2]
FRE.contrast_in_ROPE.TRIMMED.A054$upperROPE[c(7:9)]     <- RR_trim.A054[2]
FRE.contrast_in_ROPE.TRIMMED.A054$CI_range              <- FRE.contrast_in_ROPE.TRIMMED.A054$upper - FRE.contrast_in_ROPE.TRIMMED.A054$lower
FRE.contrast_in_ROPE.TRIMMED.A054$minUpper              <- FRE.contrast_in_ROPE.TRIMMED.A054 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
FRE.contrast_in_ROPE.TRIMMED.A054$maxLower              <-  FRE.contrast_in_ROPE.TRIMMED.A054 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
FRE.contrast_in_ROPE.TRIMMED.A054$DiffminUppermaxLower  <- FRE.contrast_in_ROPE.TRIMMED.A054$minUpper  - FRE.contrast_in_ROPE.TRIMMED.A054$maxLower 
FRE.contrast_in_ROPE.TRIMMED.A054$Zeros                 <- rep(0,nrow(FRE.contrast_in_ROPE.TRIMMED.A054))
FRE.contrast_in_ROPE.TRIMMED.A054$Overlap               <- FRE.contrast_in_ROPE.TRIMMED.A054  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
FRE.contrast_in_ROPE.TRIMMED.A054$perc_in_ROPE          <- (FRE.contrast_in_ROPE.TRIMMED.A054$Overlap*100)/FRE.contrast_in_ROPE.TRIMMED.A054$CI_range
FRE.contrast_in_ROPE.TRIMMED.A054[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
FRE.subj.contrast.TRIMMED.A054 <- FRE.subj.EMM.TRIMMED.A054 %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution" & probe == "without obstacle" ~ "Execution NoObs",
           trial_type == "Execution" & probe == "with obstacle"    ~ "Execution Obs",
           trial_type == "No Movement" & probe == "without obstacle" ~ "No Movement NoObs",
           trial_type == "No Movement" & probe == "with obstacle"    ~ "No Movement Obs")
  )

#print(FRE.subj.contrast.TRIMMED.A054)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
FRE.subj.contrast_pooled.TRIMMED.A054 <- FRE.subj.EMM.TRIMMED.A054 %>%
  group_by(subID, trial_type, probe, prime) %>%
  summarise( .value = mean(.value))  %>%
  select(subID, trial_type, probe, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = case_when(
    probe == "with obstacle" ~ `without obstacle` - `with obstacle`,
    probe == "without obstacle" ~ `with obstacle` - `without obstacle`,
    TRUE ~ NA_real_)
  ) %>%
  group_by(subID, trial_type) %>%      # ⬅ remove probe from grouping
  summarise(
    diff = mean(diff, na.rm = TRUE),  # ⬅ average across probes
    .groups = "drop"
  ) %>%
  mutate(contrast = trial_type)

#print(FRE.subj.contrast_pooled.TRIMMED.A054)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
FRE.subj.diffContrast.TRIMMED.A054 <- FRE.subj.contrast.TRIMMED.A054 %>%
  select(subID, trial_type, probe, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(go_vs_stop = Execution - `No Movement`,
         contrast = case_when(
           probe == "without obstacle"  ~ "EX minus NM NoObs",
           probe == "with obstacle"  ~ "EX minus NM Obs")
  )

#print(FRE.subj.diffContrast.TRIMMED.A054)


FRE.subj.diffContrast_pooled.TRIMMED.A054 <- FRE.subj.contrast_pooled.TRIMMED.A054 %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(go_vs_stop = Execution - `No Movement`,
         contrast = "EX minus NM"
  )

#print(FRE.subj.diffContrast_pooled.TRIMMED.A054)











#  ==== FINAL REACH ERROR: PLOTTING CONTRASTS  (SUPPLEMENTARY FIGURE S3GHKL) ====
# Probe without obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_NoObs.Contrast.TRIMMED.A054 <- FRE.Contrasts_long.TRIMMED.A054 %>% filter(contrast=="Execution NoObs" | contrast=="No Movement NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution NoObs" = "Execution",
    "No Movement NoObs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A054, contrast == "Execution NoObs" | contrast=="No Movement NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs_trim.A054[1], xmax = RR_NoObs_trim.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.TRIMMED.A054, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs.Contrast.TRIMMED.A054


# Probe with obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_Obs.Contrast.TRIMMED.A054 <- FRE.Contrasts_long.TRIMMED.A054 %>% filter(contrast=="Execution Obs" | contrast=="No Movement Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution Obs" = "Execution",
    "No Movement Obs" = "No Movement"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe1,color_None)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A054, contrast == "Execution Obs" | contrast=="No Movement Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs_trim.A054[1], xmax = RR_Obs_trim.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast.TRIMMED.A054, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs.Contrast.TRIMMED.A054



### Execution vs No Movement ###
# Probe without obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_NoObs_EXvsNM.Contrast.TRIMMED.A054 <- FRE.Contrasts_long.TRIMMED.A054 %>% filter(contrast=="EX minus NM NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM NoObs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A054, contrast == "EX minus NM NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs_trim.A054[1], xmax = RR_NoObs_trim.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.TRIMMED.A054, probe == "without obstacle"),
             aes(y=0, x=go_vs_stop,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs_EXvsNM.Contrast.TRIMMED.A054




# Probe with obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_Obs_EXvsNM.Contrast.TRIMMED.A054 <- FRE.Contrasts_long.TRIMMED.A054 %>% filter(contrast=="EX minus NM Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus NM Obs" = "EX minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color_exe_None)) + 
  scale_pattern_fill_manual(values = c(color_exe_None)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary.TRIMMED.A054, contrast == "EX minus NM Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs_trim.A054[1], xmax = RR_Obs_trim.A054[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast.TRIMMED.A054, probe == "with obstacle"),
             aes(y=0, x=go_vs_stop,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_None)) +
  scale_x_continuous(name= titleX, limits=c(-0.3,0.3),breaks=c(-0.25,0,0.25)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs_EXvsNM.Contrast.TRIMMED.A054





######################################################################################################
#========================================= A061 vs A054 ==============================================
######################################################################################################
# ==== Reset Paths (portable) =====
basePath   <- "C:/Experiments/A061_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")

#==== A061 vs A054: INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
IRE.Contrasts.A061_A054           <- data.frame(matrix(ncol = 2, nrow = 40000))
colnames(IRE.Contrasts.A061_A054) <- c("Difference Effect NoObs","Difference Effect Obs")

# Difference A061 minus A054) in Reach Bias Difference (EX minus MI/NM): no obstacle probe
IRE.Contrasts.A061_A054$`Difference Effect NoObs`  <-  IRE.Contrasts.A061$`EX minus MI NoObs` - IRE.Contrasts.A054$`EX minus NM NoObs`
# Difference A061 minus A054) in Reach Bias Difference (EX minus MI/NM): obstacle probe
IRE.Contrasts.A061_A054$`Difference Effect Obs`    <-  IRE.Contrasts.A061$`EX minus MI Obs` - IRE.Contrasts.A054$`EX minus NM Obs` 


IRE.Contrasts_long.A061_A054 <- pivot_longer(IRE.Contrasts.A061_A054, cols = everything(),
                                        names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Difference Effect NoObs", "Difference Effect Obs"
)

# Convert 'contrast' to a factor with this order
IRE.Contrasts_long.A061_A054$contrast <- factor(IRE.Contrasts_long.A061_A054$contrast, levels = contrast_order)

# Check
head(IRE.Contrasts_long.A061_A054)


#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
IRE.contrasts.summary.A061_A054 <-
  IRE.Contrasts_long.A061_A054 %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

IRE.contrasts.summary.A061_A054$pd <- format(IRE.contrasts.summary.A061_A054$pd, nsmall = 4)  
print(IRE.contrasts.summary.A061_A054, n = Inf, width = Inf)



### Calculate and back-transform ROPE, this is a bit tricky
## this does not take into account that values are on log scale so don't use it
## apparently, back-transformation is now implemented :)
( RR_IRE.NoObs <- (RR_IRE.NoObs.A061 + RR_IRE.NoObs.A054)/2 )
( RR_IRE.Obs   <- (RR_IRE.Obs.A061   + RR_IRE.Obs.A054)/2 )

# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE.A061_A054                       <- as.data.frame(IRE.contrasts.summary.A061_A054)
IRE.contrast_in_ROPE.A061_A054$lowerROPE             <- NA
IRE.contrast_in_ROPE.A061_A054$lowerROPE[c(1)]       <- RR_IRE.NoObs[1]
IRE.contrast_in_ROPE.A061_A054$lowerROPE[c(2)]       <- RR_IRE.Obs[1]
IRE.contrast_in_ROPE.A061_A054$upperROPE             <- NA
IRE.contrast_in_ROPE.A061_A054$upperROPE[c(1)]       <- RR_IRE.NoObs[2]
IRE.contrast_in_ROPE.A061_A054$upperROPE[c(2)]       <- RR_IRE.Obs[2]
IRE.contrast_in_ROPE.A061_A054$CI_range              <- IRE.contrast_in_ROPE.A061_A054$upper - IRE.contrast_in_ROPE.A061_A054$lower
IRE.contrast_in_ROPE.A061_A054$minUpper              <- IRE.contrast_in_ROPE.A061_A054 %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
IRE.contrast_in_ROPE.A061_A054$maxLower              <- IRE.contrast_in_ROPE.A061_A054 %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
IRE.contrast_in_ROPE.A061_A054$DiffminUppermaxLower  <- IRE.contrast_in_ROPE.A061_A054$minUpper  - IRE.contrast_in_ROPE.A061_A054$maxLower 
IRE.contrast_in_ROPE.A061_A054$Zeros                 <- rep(0,nrow(IRE.contrast_in_ROPE.A061_A054))
IRE.contrast_in_ROPE.A061_A054$Overlap               <- IRE.contrast_in_ROPE.A061_A054  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
IRE.contrast_in_ROPE.A061_A054$perc_in_ROPE          <- (IRE.contrast_in_ROPE.A061_A054$Overlap*100)/IRE.contrast_in_ROPE.A061_A054$CI_range
IRE.contrast_in_ROPE.A061_A054[,c(1:7,14)]


#==== INITIAL REACH ERROR: PLOTTING CONTRASTS (FIGURE 2 MN) ====
### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Reach~Bias~Reduction)
titleY <- ""
g.IRE_NoObs_MIvsNM.Contrast.A061_A054 <- IRE.Contrasts_long.A061_A054 %>% filter(contrast=="Difference Effect NoObs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Difference Effect NoObs" = "MI minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color1_exe_mi)) + 
  scale_pattern_fill_manual(values = c(color1_exe_mi)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A061_A054, contrast == "Difference Effect NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.NoObs[1], xmax = RR_IRE.NoObs[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  #geom_point(data=dplyr::filter(IRE.subj.diffContrast.A054, probe == "without obstacle"),
  #           aes(y=0, x=Execution_vs_NM,color=contrast),
  #           position = position_jitternudge(jitter.width =  0,
  #                                           jitter.height = 0.02,
  #                                           nudge.y = -0.05,
  #                                           nudge.x = 0,
  #                                           seed = NA),       
  #           size = 2,
  #           show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color1_exe_mi)) +
  scale_x_continuous(name= titleX, limits=c(-15,5),breaks=c(-15,-10,-5,0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= "MI minus NM") + guides(y = "none") + #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs_MIvsNM.Contrast.A061_A054




# Probe with obstacle
titleX <- expression(Delta~Reach~Bias~Reduction)
titleY <- ""
g.IRE_Obs_MIvsNM.Contrast.A061_A054 <- IRE.Contrasts_long.A061_A054 %>% filter(contrast=="Difference Effect Obs") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Difference Effect Obs" = "MI minus NM"
  ))) + 
  geom_density_pattern(
    linetype = 0,
    pattern = "stripe",  # Stripes
    pattern_density = 0.5,  # Take up 50% of the pattern (i.e. stripes equally sized)
    pattern_spacing = 0.1,  # Thicker stripes
    pattern_linetype = 0,
    pattern_size = 0,  # No border on the stripes
    trim = TRUE,  # Trim the ends of the distributions
    linewidth = 0,  # No border on the distributions
    position = position_nudge(y=0)
  ) +
  scale_fill_manual(values = c(color1_exe_mi)) + 
  scale_pattern_fill_manual(values = c(color1_exe_mi)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary.A061_A054, contrast == "Difference Effect Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_IRE.Obs[1], xmax = RR_IRE.Obs[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  #geom_point(data=dplyr::filter(IRE.subj.diffContrast.A054, probe == "without obstacle"),
  #           aes(y=0, x=Execution_vs_NM,color=contrast),
  #           position = position_jitternudge(jitter.width =  0,
  #                                           jitter.height = 0.02,
  #                                           nudge.y = -0.05,
  #                                           nudge.x = 0,
  #                                           seed = NA),       
  #           size = 2,
  #           show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color1_exe_mi)) +
  scale_x_continuous(name= titleX, limits=c(-15,5),breaks=c(-15,-10,-5,0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= "MI minus NM") + guides(y = "none") + #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs_MIvsNM.Contrast.A061_A054






#==== A061 vs A054: INITIAL REACH ERROR: COMBINE PLOTS (FIGURE 2) =====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 6, r = 5), # p1
  area(t = 1, l = 6, b = 6, r = 10), # p2
  area(t = 1, l = 11, b = 6, r = 15),  # p3
  area(t = 1, l = 16, b = 6, r = 20),   #p4
  area(t = 7, l = 1, b = 12, r = 4),   #p5
  area(t = 7, l = 5, b = 12, r = 8),   #p6
  area(t = 7, l = 9, b = 12, r = 12),   #p7
  area(t = 7, l = 13, b = 12, r = 16),   #p8
  area(t = 13, l = 1, b = 15, r = 4),   #p9
  area(t = 13, l = 5, b = 15, r = 8),   #p10
  area(t = 13, l = 9, b = 15, r = 12),   #p11
  area(t = 13, l = 13, b = 15, r = 16),  #p12
  area(t = 7, l = 17, b = 11, r = 20),  #p13
  area(t = 12, l = 17, b = 15, r = 20)  #p14
)



#merge into a 1x2 plot grid (library patchwork)
g.IRE.Contrast_EMM.A061_A054 <-
  g.IRE_Obs.EMM.A061 + g.IRE_NoObs.EMM.A061 + g.IRE_Obs.EMM.A054 + g.IRE_NoObs.EMM.A054 +
     g.IRE_Obs.Contrast.A061 + g.IRE_NoObs.Contrast.A061 + g.IRE_Obs.Contrast.A054 + g.IRE_NoObs.Contrast.A054 +
     g.IRE_Obs_EXvsMI.Contrast.A061 + g.IRE_NoObs_EXvsMI.Contrast.A061 + g.IRE_Obs_EXvsNM.Contrast.A054 + g.IRE_NoObs_EXvsNM.Contrast.A054 +
     g.IRE_Obs_MIvsNM.Contrast.A061_A054 + g.IRE_NoObs_MIvsNM.Contrast.A061_A054 +
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
  #theme(legend.position = "bottom")

g.IRE.Contrast_EMM.A061_A054


# Build a valid path
outfile <- file.path(figurePath, "Fig2_A061_A054_InitialReachError_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.IRE.Contrast_EMM.A061_A054,
  width = 36, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "Fig2_A061_A054_InitialReachError_EMM_Contrasts_ALL.svg")

ggsave(
  filename = outfile,
  plot = g.IRE.Contrast_EMM.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "Fig2_A061_A054_InitialReachError_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.IRE.Contrast_EMM.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)







#==== A061 vs A054: FINAL REACH ERROR: COMBINE PLOTS (FIGURE 4) ====
layout <- c(
  area(t = 1, l = 1, b = 6, r = 5), # p1
  area(t = 1, l = 6, b = 6, r = 10), # p2
  area(t = 1, l = 11, b = 6, r = 15),  # p3
  area(t = 1, l = 16, b = 6, r = 20),   #p4
  area(t = 7, l = 1, b = 12, r = 5),   #p5
  area(t = 7, l = 6, b = 12, r = 10),   #p6
  area(t = 7, l = 11, b = 12, r = 15),   #p7
  area(t = 7, l = 16, b = 12, r = 20),   #p8
  area(t = 13, l = 1, b = 15, r = 5),   #p9
  area(t = 13, l = 6, b = 15, r = 10),   #p10
  area(t = 13, l = 11, b = 15, r = 15),   #p11
  area(t = 13, l = 16, b = 15, r = 20) #p12
)


#merge into a 1x2 plot grid (library patchwork)
g.FRE.Contrast_EMM.A061_A054 <-
  g.FRE_Obs.EMM.A061 + g.FRE_NoObs.EMM.A061 + g.FRE_Obs.EMM.A054 + g.FRE_NoObs.EMM.A054 +
  g.FRE_Obs.Contrast.A061 + g.FRE_NoObs.Contrast.A061 + g.FRE_Obs.Contrast.A054 + g.FRE_NoObs.Contrast.A054 +
  g.FRE_Obs_EXvsMI.Contrast.A061 + g.FRE_NoObs_EXvsMI.Contrast.A061 + g.FRE_Obs_EXvsNM.Contrast.A054 + g.FRE_NoObs_EXvsNM.Contrast.A054 +
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
#theme(legend.position = "bottom")

g.FRE.Contrast_EMM.A061_A054



# Build a valid path
outfile <- file.path(figurePath, "Fig4_A061_A054_FinalReachError_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.FRE.Contrast_EMM.A061_A054,
  width = 36, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)


# Build a valid path
outfile <- file.path(figurePath, "Fig4_A061_A054_FinalReachError_EMM_Contrasts_ALL.svg")

ggsave(
  filename = outfile,
  plot = g.FRE.Contrast_EMM.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "Fig4_A061_A054_FinalReachError_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.FRE.Contrast_EMM.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)



#==== A061 vs A054: FINAL REACH ERROR AS A FUNCTION OF INITIAL REACH ERROR: COMBINE PLOTS (SUPPLEMENTARY FIGURE S2) =====
#merge into a 2x2 plot grid (library patchwork)

g.FRE_IRE.A061_A054 <- g.FRE_IRE_GoObs.A061 + g.FRE_IRE_GoNoObs.A061 + g.FRE_IRE_GoObs.A054 + g.FRE_IRE_GoNoObs.A054 + 
  g.FRE_IRE_StopObs.A061 + g.FRE_IRE_StopNoObs.A061 + g.FRE_IRE_StopObs.A054 + g.FRE_IRE_StopNoObs.A054 +
  plot_layout(widths = c(1, 1, 1, 1), guides='collect') +   # relation of x-axes
  plot_annotation(tag_levels = "A") & # adding panel labels
  theme(legend.position='bottom')  # only 1 legend
g.FRE_IRE.A061_A054




# Build a valid path
outfile <- file.path(figurePath, "FigS2_A061_A054_FinalReachError_InitialReachError.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.FRE_IRE.A061_A054,
  width = 36, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS2_A061_A054_FinalReachError_InitialReachError.svg")

ggsave(
  filename = outfile,
  plot = g.FRE_IRE.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "FigS2_A061_A054_FinalReachError_InitialReachError.png")

ggsave(
  filename = outfile,
  plot = g.FRE_IRE.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)












#==== A061 vs A054: ACCURACY END POINT ERROR: COMBINE PLOTS (SUPPLEMENTARY FIGURE S5) ====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 6, r = 5), # p1
  area(t = 1, l = 6, b = 6, r = 10), # p2
  area(t = 1, l = 11, b = 6, r = 15),  # p3
  area(t = 1, l = 16, b = 6, r = 20),   #p4
  area(t = 7, l = 1, b = 11, r = 4),   #p5
  area(t = 12,l = 1, b = 15, r = 4),   #p6
  
  area(t = 7, l = 5, b = 12, r = 8),   #p7
  area(t = 7, l = 9, b = 12, r = 12),   #p8
  area(t = 7, l = 13, b = 12, r = 16),   #p9
  area(t = 7, l = 17, b = 12, r = 20),   #p10
  
  area(t = 13, l = 5, b = 15, r = 8),   #p11
  area(t = 13, l = 9, b = 15, r = 12),  #p12
  area(t = 13, l = 13, b = 15, r = 16),  #p13
  area(t = 13, l = 17, b = 15, r = 20)  #p14
)


#merge into a 1x2 plot grid (library patchwork)
g.EPE150.Contrast.A061_A054 <-
  g.EPE150_Obs.EMM.A061 + g.EPE150_NoObs.EMM.A061 + g.EPE150_Obs.EMM.A054 + g.EPE150_NoObs.EMM.A054 +
  g.EPE150_EXvsMI.Contrast.A061 + g.EPE150_EXvsNM.Contrast.A054 + 
  g.EPE150_Obs.Contrast.A061 + g.EPE150_NoObs.Contrast.A061 + g.EPE150_Obs.Contrast.A054 + g.EPE150_NoObs.Contrast.A054 +
  g.EPE150_Obs_EXvsMI.Contrast.A061 + g.EPE150_NoObs_EXvsMI.Contrast.A061 + g.EPE150_Obs_EXvsNM.Contrast.A054 + g.EPE150_NoObs_EXvsNM.Contrast.A054 +
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
g.EPE150.Contrast.A061_A054



# Build a valid path
outfile <- file.path(figurePath, "FigS5_A061_A054_EPE150_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.EPE150.Contrast.A061_A054,
  width = 45, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS5_A061_A054_EPE150_EMM_Contrasts_ALL.svg")

ggsave(
  filename = outfile,
  plot = g.EPE150.Contrast.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "FigS5_A061_A054_EPE150_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.EPE150.Contrast.A061_A054,
  width = 45, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)













#==== A061 vs A054: PRECISION END POINT ERROR: COMBINE PLOTS (SUPPLEMENTARY FIGURE S4) ====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 6, r = 5), # p1
  area(t = 1, l = 6, b = 6, r = 10), # p2
  area(t = 1, l = 11, b = 6, r = 15),  # p3
  area(t = 1, l = 16, b = 6, r = 20),   #p4
  area(t = 7, l = 1, b = 11, r = 4),   #p5
  area(t = 12,l = 1, b = 15, r = 4),   #p6
  
  area(t = 7, l = 5, b = 12, r = 8),   #p7
  area(t = 7, l = 9, b = 12, r = 12),   #p8
  area(t = 7, l = 13, b = 12, r = 16),   #p9
  area(t = 7, l = 17, b = 12, r = 20),   #p10
  
  area(t = 13, l = 5, b = 15, r = 8),   #p11
  area(t = 13, l = 9, b = 15, r = 12),  #p12
  area(t = 13, l = 13, b = 15, r = 16),  #p13
  area(t = 13, l = 17, b = 15, r = 20)  #p14
)


#merge into a 1x2 plot grid (library patchwork)
g.EPP150.Contrast.A061_A054 <-
  g.EPP150_Obs.EMM.A061 + g.EPP150_NoObs.EMM.A061 + g.EPP150_Obs.EMM.A054 + g.EPP150_NoObs.EMM.A054 +
  g.EPP150_EXvsMI.Contrast.A061 + g.EPP150_EXvsNM.Contrast.A054 + 
  g.EPP150_Obs.Contrast.A061 + g.EPP150_NoObs.Contrast.A061 + g.EPP150_Obs.Contrast.A054 + g.EPP150_NoObs.Contrast.A054 +
  g.EPP150_Obs_EXvsMI.Contrast.A061 + g.EPP150_NoObs_EXvsMI.Contrast.A061 + g.EPP150_Obs_EXvsNM.Contrast.A054 + g.EPP150_NoObs_EXvsNM.Contrast.A054 +
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
g.EPP150.Contrast.A061_A054



# Build a valid path
outfile <- file.path(figurePath, "FigS4_A061_A054_EPP150_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.EPP150.Contrast.A061_A054,
  width = 45, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS4_A061_A054_EPP150_EMM_Contrasts_ALL.svg")

ggsave(
  filename = outfile,
  plot = g.EPP150.Contrast.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "FigS4_A061_A054_EPP150_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.EPP150.Contrast.A061_A054,
  width = 45, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)













#==== A061 vs A054: TRIMMED: FINAL REACH ERROR: COMBINE PLOTS (SUPPLEMENTARY FIGURE S3) ====
layout <- c(
  area(t = 1, l = 1, b = 6, r = 5), # p1
  area(t = 1, l = 6, b = 6, r = 10), # p2
  area(t = 1, l = 11, b = 6, r = 15),  # p3
  area(t = 1, l = 16, b = 6, r = 20),   #p4
  area(t = 7, l = 1, b = 12, r = 5),   #p5
  area(t = 7, l = 6, b = 12, r = 10),   #p6
  area(t = 7, l = 11, b = 12, r = 15),   #p7
  area(t = 7, l = 16, b = 12, r = 20),   #p8
  area(t = 13, l = 1, b = 15, r = 5),   #p9
  area(t = 13, l = 6, b = 15, r = 10),   #p10
  area(t = 13, l = 11, b = 15, r = 15),   #p11
  area(t = 13, l = 16, b = 15, r = 20) #p12
)


#merge into a 1x2 plot grid (library patchwork)
g.FRE.Contrast_EMM.TRIMMED.A061_A054 <-
  g.FRE_Obs.EMM.TRIMMED.A061 + g.FRE_NoObs.EMM.TRIMMED.A061 + g.FRE_Obs.EMM.TRIMMED.A054 + g.FRE_NoObs.EMM.TRIMMED.A054 +
  g.FRE_Obs.Contrast.TRIMMED.A061 + g.FRE_NoObs.Contrast.TRIMMED.A061 + g.FRE_Obs.Contrast.TRIMMED.A054 + g.FRE_NoObs.Contrast.TRIMMED.A054 +
  g.FRE_Obs_EXvsMI.Contrast.TRIMMED.A061 + g.FRE_NoObs_EXvsMI.Contrast.TRIMMED.A061 + g.FRE_Obs_EXvsNM.Contrast.TRIMMED.A054 + g.FRE_NoObs_EXvsNM.Contrast.TRIMMED.A054 +
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
#theme(legend.position = "bottom")

g.FRE.Contrast_EMM.TRIMMED.A061_A054



# Build a valid path
outfile <- file.path(figurePath, "FigS3_A061_A054_FinalReachError_EMM_Contrasts_TRIMMED_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.FRE.Contrast_EMM.TRIMMED.A061_A054,
  width = 36, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)


# Build a valid path
outfile <- file.path(figurePath, "FigS3_A061_A054_FinalReachError_EMM_Contrasts_TRIMMED_ALL.svg")

ggsave(
  filename = outfile,
  plot = g.FRE.Contrast_EMM.TRIMMED.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "FigS3_A061_A054_FinalReachError_EMM_Contrasts_TRIMMED_ALL.png")

ggsave(
  filename = outfile,
  plot = g.FRE.Contrast_EMM.TRIMMED.A061_A054,
  width = 36, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)





