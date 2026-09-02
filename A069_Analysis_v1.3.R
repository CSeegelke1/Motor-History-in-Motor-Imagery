# Title: Distinct contributions of motor imagery and execution to history-dependent biases in reaching (Experiment 3, A069)
# Authors: Seegelke, Heed 
# Script: Christian Seegelke 01/09/2026
# INPUT:  A069_data.csv (preprocessed data), A069_SubjInfo.xlsx (Demographics & MI Questionnaire data)
# OUTPUT: REPORTED STATS OF EXPERIMENT 3, SUPPLEMENTARY FIGURES S6, S7, S8, S9, S10, S11, S12, S13 Supplementary Table S2, S3, S4
#========================================================================================================================
#===================================================================================================================

# ===== REPRODUCABILITY =====
set.seed(1234)
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  scipen = 999,            # fewer scientific notation surprises
  warn = 1                 # show warnings as they occur
)

# ===== INSTALL PACKAGES IF MISSING =====
pkgs <- c(
  "tidyverse","rio","psych","afex","emmeans","cowplot","patchwork","ggpubr","ggpattern",
  "sdamr","brms","GGally","see","bayestestR","tidybayes","BayesFactor",
  "flextable","officer","conflicted","renv","scales","rlang"
)
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install)

# ===== LOAD PACKAGES =====
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
})


# ===== Resolve common conflicts explicitly =====
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

# ===== afex / contrasts / emmeans defaults =====
# Set Type-III SS compatible contrasts for factorial designs
options(contrasts = c("contr.sum", "contr.poly"))
afex::afex_options(
  type = 3,
  es_aov = "ges",
  return_aov = "afex_aov",
  emmeans_model = "multivariate"
)


# ===== brms / Stan options =====
# Use all cores and a reasonable default backend
options(mc.cores = parallel::detectCores())
# Faster compilation for development; consider O3 for final runs
Sys.setenv(LOCAL_CPPFLAGS = "-O2")
# If you use cmdstanr instead of rstan, set: brms.backend = "cmdstanr"
# options(brms.backend = "cmdstanr")



#==== IMPORT SINGLE TRIAL DATA  ====
# clear workspace
# rm(list = ls(all = TRUE))  # Generally avoid in shared scripts

# ===== Paths (portable) =====
basePath   <- "C:/Experiments/A069_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")





# ===== Load single-trial data =====
# Note: N = 30, 704 trials (incl. 64 practise trials) per subject
# Expecting column 'subID' present
d <- import(file.path(dataPath, "A069_data.csv")) %>%
  as_tibble()


# ===== Load demographics (sheet 2) =====
# Ensure there is a subject identifier column to join on.
# If sheet lacks an explicit subID column but rows correspond to subjects 1..32,
# create it via row_number().
d_SubInfo_raw <- import(file.path(docPath, "A069_SubjInfo.xlsx"), which = 2) %>%
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




# ===== Load MIQ-RS scores (sheet 3) and compute subscale/overall =====
d_MIQ_raw <- import(file.path(docPath, "A069_SubjInfo.xlsx"), which = 3) %>%
  as_tibble()

if (!"subID" %in% names(d_MIQ_raw)) {
  d_MIQ <- d_MIQ_raw %>%
    mutate(subID = row_number())
} else {
  d_MIQ <- d_MIQ_raw
}

# Identify K and V items robustly
k_items <- names(d_MIQ)[str_detect(names(d_MIQ), "K$")]
v_items <- names(d_MIQ)[str_detect(names(d_MIQ), "V$")]

# Coerce to numeric and compute means rowwise (ignoring NA)
d_MIQ <- d_MIQ %>%
  mutate(across(all_of(c(k_items, v_items)), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    MIQ_Score_Kin = rowMeans(across(all_of(k_items)), na.rm = TRUE),
    MIQ_Score_Vis = rowMeans(across(all_of(v_items)), na.rm = TRUE),
    MIQ_Score_Overall = (MIQ_Score_Kin + MIQ_Score_Vis) / 2
  ) %>%
  select(subID, MIQ_Score_Kin, MIQ_Score_Vis, MIQ_Score_Overall)


# ===== Load MI “Ease” per block (sheet 4), compute overall (mean across non-missing blocks) =====
d_Ease_MI_raw <- import(file.path(docPath, "A069_SubjInfo.xlsx"), which = 4) %>%
  as_tibble()

if (!"subID" %in% names(d_Ease_MI_raw)) {
  d_Ease_MI <- d_Ease_MI_raw %>%
    mutate(subID = row_number())
} else {
  d_Ease_MI <- d_Ease_MI_raw
}

block_cols <- names(d_Ease_MI)[str_detect(names(d_Ease_MI), "^Block\\d{2}$")]

d_Ease_MI <- d_Ease_MI %>%
  mutate(across(all_of(block_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    nmb_col_val = rowSums(!is.na(across(all_of(block_cols)))),
    MI_Ease_Overall = rowMeans(across(all_of(block_cols)), na.rm = TRUE)
  ) %>%
  select(subID, MI_Ease_Overall)


# ===== Load MI “Count” per block (sheet 5), compute overall =====
d_Count_MI_raw <- import(file.path(docPath, "A069_SubjInfo.xlsx"), which = 5) %>%
  as_tibble()

if (!"subID" %in% names(d_Count_MI_raw)) {
  d_Count_MI <- d_Count_MI_raw %>%
    mutate(subID = row_number())
} else {
  d_Count_MI <- d_Count_MI_raw
}

count_block_cols <- names(d_Count_MI)[str_detect(names(d_Count_MI), "^Block\\d{2}$")]

d_Count_MI <- d_Count_MI %>%
  mutate(across(all_of(count_block_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    nmb_col_val = rowSums(!is.na(across(all_of(count_block_cols)))),
    # FIXED: divide by Count's own valid column count (or use rowMeans)
    MI_Count_Overall = rowMeans(across(all_of(count_block_cols)), na.rm = TRUE)
  ) %>%
  select(subID, MI_Count_Overall)

# ===== Merge MIQ-RS/ MI ease/count to trial data =====
d <- d %>%
  left_join(d_MIQ,      by = "subID") %>%
  left_join(d_Ease_MI,  by = "subID") %>%
  left_join(d_Count_MI, by = "subID")

# ===== Quick checks =====
glimpse(d)
summary(select(d, Age, EHI, starts_with("MIQ_"), MI_Ease_Overall, MI_Count_Overall))

# Example sanity plots
# library(ggplot2)
ggplot(d, aes(subID, Age)) + geom_point(alpha = 0.2, position = position_jitter(width = 0.15))
ggplot(d, aes(subID, MIQ_Score_Overall)) + geom_point(alpha = 0.2, position = position_jitter(width = 0.15))
ggplot(d, aes(subID, MI_Ease_Overall)) + geom_point(alpha = 0.2, position = position_jitter(width = 0.15))
ggplot(d, aes(subID, MI_Count_Overall)) + geom_point(alpha = 0.2, position = position_jitter(width = 0.15))



# 1) Remove training block BEFORE converting subID to factor
#    (or convert with character comparison)
d <- d %>%
  filter(blocks_thisN != 0)


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


dClean.A069 <- dClean









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




#==== SETTINGS FOR PLOTTING ====
#define some graphical params like themes
custom_plot_theme <- theme(strip.background =element_rect(fill="white", linewidth = 2),
                           strip.text = element_text(size = rel(1), margin = margin(1,5,5,0, "pt")), #in ggplot2 clockwise starting from top: trbl
                           plot.title = element_text(size = rel(1.5)),
                           panel.background = element_blank())

# for probe actions
color_exe1    <- "#008000"
color_mi1     <- "#000080"
color_Neut    <- "#9C9C9C"
color_exe_mi1 <- "#008080"
color_None    <- "#808000"

# For prime actions
color_exe2     <- "#00c000"
color_mi2      <- "#0000c0"
color_exe_mi2  <- "#00c0c0"
color_Neut2    <- "#BDB8AD"


#==== END OF PREPROCESSING =========================================================================================




######################################################################################################
#========================================= BAYESIAN REGRESSION MODELS ================================
######################################################################################################
#==== INITIAL REACH ERROR PROBE ====
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe
# ==== MODEL FITTING: INITIAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (UB=180) LOGNORMAL DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "Reach Direction Probe" ) +
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
A069_fit.IRE_Probe_NoObs_OPRIxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.IRE_Probe_NoObs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )






#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSS,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSS,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSS,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.IRE_Probe_NoObs_OPRIxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.IRE_Probe_NoObs_OPRIxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A069_fit.IRE_Probe_NoObs_OPRIxSS,
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
  A069_fit.IRE_Probe_NoObs_OPRIxSS,
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
  labs( title = "A069 - Histogram", subtitle = "Reach Direction Probe" ) +
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
A069_fit.IRE_Probe_Obs_OPRIxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.IRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )






#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.IRE_Probe_Obs_OPRIxSS,   type = "dens_overlay", ndraws = 100)  # good
pp_check(A069_fit.IRE_Probe_Obs_OPRIxSS,   type = "hist", ndraws = 10)           # 
pp_check(A069_fit.IRE_Probe_Obs_OPRIxSS,   type = "boxplot", ndraws = 10)        # 


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.IRE_Probe_Obs_OPRIxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.IRE_Probe_Obs_OPRIxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)








#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A069_fit.IRE_Probe_Obs_OPRIxSS,
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
  A069_fit.IRE_Probe_Obs_OPRIxSS,
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





# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY S6AB) ====
# Prepare and tidy data
IRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
IRE.EMM <- rbind(IRE_NoObs.EMM,IRE_Obs.EMM)

IRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
IRE.subj.EMM <- rbind(IRE_NoObs.subj.EMM,IRE_Obs.subj.EMM)

# rename factors
IRE.EMM <- IRE.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

IRE.EMM

# rename factors
IRE.subj.EMM <- IRE.subj.EMM  %>% 
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

g.IRE_NoObs.EMM.A069 <-
  IRE.subj.EMM %>%
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
    data = dplyr::filter(IRE.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM, probe == "without obstacle"),
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

g.IRE_NoObs.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_Obs.EMM.A069 <-
  IRE.subj.EMM %>%
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
    data = dplyr::filter(IRE.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM, probe == "with obstacle"),
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

g.IRE_Obs.EMM.A069



# ==== INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
IRE.Contrasts           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(IRE.Contrasts) <- c("Execution", "Motor Imagery", "EX minus MI", 
                             "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                             "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials
IRE.Contrasts$`Execution`           <- (IRE_Probe_NoObs.posteriors$`yes_go` + IRE_Probe_Obs.posteriors$`yes_go`)/2 - ( IRE_Probe_NoObs.posteriors$`no_go` + IRE_Probe_Obs.posteriors$`no_go`)/2
# Motor Imagery trials
IRE.Contrasts$`Motor Imagery`       <- (IRE_Probe_NoObs.posteriors$`yes_stop` + IRE_Probe_Obs.posteriors$`yes_stop`)/2 - ( IRE_Probe_NoObs.posteriors$`no_stop` + IRE_Probe_Obs.posteriors$`no_stop`)/2
# Execution vs Motor Imagery
IRE.Contrasts$`EX minus MI`         <- (IRE.Contrasts$`Execution` - IRE.Contrasts$`Motor Imagery`)
# Execution trials Probe without Obstacle
IRE.Contrasts$`Execution NoObs`     <- (IRE_Probe_NoObs.posteriors$`yes_go` - IRE_Probe_NoObs.posteriors$`no_go`)
# Execution trials with Obstacle
IRE.Contrasts$`Execution Obs`       <- (IRE_Probe_Obs.posteriors$`yes_go` - IRE_Probe_Obs.posteriors$`no_go`)
# Motor Imagery trials Probe without Obstacle
IRE.Contrasts$`Motor Imagery NoObs` <- (IRE_Probe_NoObs.posteriors$`yes_stop` - IRE_Probe_NoObs.posteriors$`no_stop`)
# Motor Imagery trials Probe with Obstacle
IRE.Contrasts$`Motor Imagery Obs`   <- (IRE_Probe_Obs.posteriors$`yes_stop` - IRE_Probe_Obs.posteriors$`no_stop`)
# Execution vs Motor Imagery without Obstacle
IRE.Contrasts$`EX minus MI NoObs`   <- IRE.Contrasts$`Execution NoObs` - IRE.Contrasts$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
IRE.Contrasts$`EX minus MI Obs`     <- IRE.Contrasts$`Execution Obs` - IRE.Contrasts$`Motor Imagery Obs`


IRE.Contrasts_long <- pivot_longer(IRE.Contrasts, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "Motor Imagery Obs",
  "Execution NoObs", "Motor Imagery NoObs", 
  "EX minus MI Obs", "EX minus MI NoObs", 
  "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
IRE.Contrasts_long$contrast <- factor(IRE.Contrasts_long$contrast, levels = contrast_order)

# Check
head(IRE.Contrasts_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
IRE.contrasts.summary <-
  IRE.Contrasts_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

IRE.contrasts.summary$pd <- format(IRE.contrasts.summary$pd, nsmall = 4)  
#print(IRE.contrasts.summary, n = Inf, width = Inf)


#### Calculate ROPE
### Calculate and back-transform ROPE, this is a bit tricky
## this does not take into account that values are on log scale so don't use it
## apparently, back-transformation is now implemented :)
#RR_NoObs <- rope_range(A069_fit.IRE_Probe_NoObs_OPRIxSS) 
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

( RR_NoObs <- c(-rope_value_RE_log, rope_value_RE_log) )
( RR_Obs   <- rope_range(A069_fit.IRE_Probe_Obs_OPRIxSS) )
( RR       <- (RR_NoObs + RR_Obs) /2 )



# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE                       <- as.data.frame(IRE.contrasts.summary)
IRE.contrast_in_ROPE$lowerROPE             <- NA
IRE.contrast_in_ROPE$lowerROPE[c(1, 2, 5)] <- RR_Obs[1]
IRE.contrast_in_ROPE$lowerROPE[c(3, 4, 6)] <- RR_NoObs[1]
IRE.contrast_in_ROPE$lowerROPE[c(7:9)]     <- RR[1]
IRE.contrast_in_ROPE$upperROPE             <- NA
IRE.contrast_in_ROPE$upperROPE[c(1, 2, 5)] <- RR_Obs[2]
IRE.contrast_in_ROPE$upperROPE[c(3, 4, 6)] <- RR_NoObs[2]
IRE.contrast_in_ROPE$upperROPE[c(7:9)]     <- RR[2]
IRE.contrast_in_ROPE$CI_range              <- IRE.contrast_in_ROPE$upper - IRE.contrast_in_ROPE$lower
IRE.contrast_in_ROPE$minUpper              <- IRE.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
IRE.contrast_in_ROPE$maxLower              <-  IRE.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
IRE.contrast_in_ROPE$DiffminUppermaxLower  <- IRE.contrast_in_ROPE$minUpper  - IRE.contrast_in_ROPE$maxLower 
IRE.contrast_in_ROPE$Zeros                 <- rep(0,nrow(IRE.contrast_in_ROPE))
IRE.contrast_in_ROPE$Overlap               <- IRE.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
IRE.contrast_in_ROPE$perc_in_ROPE          <- (IRE.contrast_in_ROPE$Overlap*100)/IRE.contrast_in_ROPE$CI_range
IRE.contrast_in_ROPE[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
IRE.subj.contrast <- IRE.subj.EMM %>%
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

print(IRE.subj.contrast)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
IRE.subj.contrast_pooled <- IRE.subj.EMM %>%
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

print(IRE.subj.contrast_pooled)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
IRE.subj.diffContrast <- IRE.subj.contrast %>%
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

print(IRE.subj.diffContrast)


IRE.subj.diffContrast_pooled <- IRE.subj.contrast_pooled %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

print(IRE.subj.diffContrast_pooled)







#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS (SUPPLEMENTARY FIGURE S6CDEF) ====
# Probe without obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_NoObs.Contrast.A069 <- IRE.Contrasts_long %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs[1], xmax = RR_NoObs[2], ymin = -0.1, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast, probe == "without obstacle"),
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
  scale_y_continuous(limits = c(-0.1, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs.Contrast.A069


# Probe with obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_Obs.Contrast.A069 <- IRE.Contrasts_long %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs[1], xmax = RR_Obs[2], ymin = -0.1, ymax = 0.4, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast, probe == "with obstacle"),
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
  scale_y_continuous(limits = c(-0.1, 0.4)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") + #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs.Contrast.A069



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_NoObs_EXvsMI.Contrast.A069 <- IRE.Contrasts_long %>% filter(contrast=="EX minus MI NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs[1], xmax = RR_NoObs[2], ymin = -0.1, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast, probe == "without obstacle"),
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
  scale_y_continuous(limits = c(-0.1, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs_EXvsMI.Contrast.A069




# Probe with obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_Obs_EXvsMI.Contrast.A069 <- IRE.Contrasts_long %>% filter(contrast=="EX minus MI Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs[1], xmax = RR_Obs[2], ymin = -0.1, ymax = 0.4, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast, probe == "with obstacle"),
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
  scale_y_continuous(limits = c(-0.1, 0.4)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs_EXvsMI.Contrast.A069







#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS POOLED OVER OBSTACLE PROBE ====
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast2.A069 <- IRE.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution" = "Execution",
    "Motor Imagery" = "Motor Imagery"
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.contrast_pooled),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,15),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast2.A069



### Execution vs Motor Imagery ###
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_EXvsMI.Contrast2.A069 <- IRE.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI" = "EX minus MI"
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.1, ymax = 0.7, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(IRE.subj.diffContrast_pooled),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.02,
                                             nudge.y = -0.05,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,15),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_EXvsMI.Contrast2.A069









#  ==== INITIAL REACH ERROR: COMBINE PLOTS (SUPPLEMENTARY FIGURE S6) ====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 2, r = 1), # p1
  area(t = 1, l = 2, b = 2, r = 2), # p2
  area(t = 3, l = 1, b = 4, r = 1),  # p3
  area(t = 3, l = 2, b = 4, r = 2),   #p4
  area(t = 5, l = 1, b = 5, r = 1),   #p5
  area(t = 5, l = 2, b = 5, r = 2)  #p6
  
)


g.RE.A069 <-  g.IRE_Obs.EMM.A069 + g.IRE_NoObs.EMM.A069 + 
  g.IRE_Obs.Contrast.A069 + g.IRE_NoObs.Contrast.A069 + 
  g.IRE_Obs_EXvsMI.Contrast.A069 + g.IRE_NoObs_EXvsMI.Contrast.A069 + 
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
#theme(legend.position = "bottom")
g.RE.A069



# Build a valid path
outfile <- file.path(figurePath, "FigS6_A069_InitialReachError_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.RE.A069,
  width = 24, height = 28, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)


# Build a valid path
outfile <- file.path(figurePath, "FigS6_A069_InitialReachError_EMM_Contrasts_ALL.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.RE.A069,
  width = 24, height = 28, units = "cm",
  device = "svg"
)



outfile <- file.path(figurePath, "FigS6_A069_InitialReachError_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.RE.A069,
  width = 24, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)




#==== INITIAL REACH ERROR PROBE WITH TARGET LOCATIONS ====
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe
# ==== MODEL FITTING: INITIAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (UB=180) LOGNORMAL DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime:deg_targets_probe) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_noobs$ReachDiff2_Probe)
sd(log(dClean_noobs$ReachDiff2_Probe))

# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime * deg_targets_probe + ( 1 | subID )),
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
                           prior(exponential(0.8),               class = "sigma")      
                           )      


# we model the mu here, we include random effects for participants
# runs about 45 min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A069_fit.IRE_Probe_NoObs_OPRIxSSxTL <- brm( 
  bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime * deg_targets_probe + ( 1 | subID )),     # model specification
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
  file = paste0(analysisPath,'/', 'A069_fit.IRE_Probe_NoObs_OPRIxSSxTL'), #save model
)
endTime <- Sys.time()
( endTime- startTime )









#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSSxTL,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSSxTL,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSSxTL,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.IRE_Probe_NoObs_OPRIxSSxTL)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.IRE_Probe_NoObs_OPRIxSSxTL,          
     combo = c("dens", "trace"),
     ask   = FALSE
)


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime),
                        deg_targets_probe   = levels(dClean_noobs$deg_targets_probe))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A069_fit.IRE_Probe_NoObs_OPRIxSSxTL,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS_TL = interaction(obstacle_prime, stop_signal_prime, deg_targets_probe, sep = "_", drop = TRUE))
colnames(IRE_Probe_NoObs.posteriors) <- exp.cond_posterior$OPRI_SS_TL


# Long format (use pivot_longer)
IRE_Probe_NoObs.posteriors_long <- IRE_Probe_NoObs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS_TL", values_to = "value") %>%
  separate(OPRI_SS_TL, into = c("obstacle_prime", "stop_signal_prime", "deg_targets_probe"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_NoObs.EMM <- IRE_Probe_NoObs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime, deg_targets_probe) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_NoObs.EMM



# Conditions
dClean_noobs$subID <- droplevels(dClean_noobs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_noobs$subID),
                                  obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_noobs$stop_signal_prime),
                                  deg_targets_probe   = levels(dClean_noobs$deg_targets_probe))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A069_fit.IRE_Probe_NoObs_OPRIxSSxTL,
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
  group_by(subID, obstacle_prime, stop_signal_prime, deg_targets_probe) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)





# ==== MODEL FITTING: INITIAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=180) GAUSSIAN DISTRIBUTION ====
dClean_obs <- dClean %>% filter(obstacle_probe=="yes")

dClean_obs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime:deg_targets_probe) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_obs$ReachDiff2_Probe)

# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime * deg_targets_probe + ( 1 | subID )),
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
                          prior(student_t(3, 0,  15   ), class = "sigma"))


# we model the mu here. We include random effects for participants
# runs about 1.5 hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A069_fit.IRE_Probe_Obs_OPRIxSSxTL <- brm( 
  bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime * deg_targets_probe + ( 1 | subID )),     # model specification
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
  file = paste0(analysisPath,'/', 'A069_fit.IRE_Probe_Obs_OPRIxSSxTL'), #save model
)
endTime <- Sys.time()
( endTime- startTime )









#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.IRE_Probe_Obs_OPRIxSSxTL,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.IRE_Probe_Obs_OPRIxSSxTL,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.IRE_Probe_Obs_OPRIxSSxTL,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.IRE_Probe_Obs_OPRIxSSxTL)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.IRE_Probe_Obs_OPRIxSSxTL,          
     combo = c("dens", "trace"),
     ask   = FALSE
)


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime),
                        deg_targets_probe   = levels(dClean_obs$deg_targets_probe))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A069_fit.IRE_Probe_Obs_OPRIxSSxTL,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS_TL = interaction(obstacle_prime, stop_signal_prime, deg_targets_probe, sep = "_", drop = TRUE))
colnames(IRE_Probe_Obs.posteriors) <- exp.cond_posterior$OPRI_SS_TL


# Long format (use pivot_longer)
IRE_Probe_Obs.posteriors_long <- IRE_Probe_Obs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS_TL", values_to = "value") %>%
  separate(OPRI_SS_TL, into = c("obstacle_prime", "stop_signal_prime", "deg_targets_probe"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_Obs.EMM <- IRE_Probe_Obs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime, deg_targets_probe) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_Obs.EMM



# Conditions
dClean_obs$subID <- droplevels(dClean_obs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_obs$subID),
                                  obstacle_prime     = levels(dClean_obs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_obs$stop_signal_prime),
                                  deg_targets_probe   = levels(dClean_obs$deg_targets_probe))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A069_fit.IRE_Probe_Obs_OPRIxSSxTL,
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
  group_by(subID, obstacle_prime, stop_signal_prime, deg_targets_probe) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)






# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) ====
# Prepare and tidy data
IRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
IRE.EMM <- rbind(IRE_NoObs.EMM,IRE_Obs.EMM)

IRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
IRE.subj.EMM <- rbind(IRE_NoObs.subj.EMM,IRE_Obs.subj.EMM)

# rename factors
IRE.EMM <- IRE.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "EX", "MI")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    target      = factor(deg_targets_probe)
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe,deg_targets_probe)) # remove old columns

print(IRE.EMM, n = Inf)

# rename factors
IRE.subj.EMM <- IRE.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "EX", "MI")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    target      = factor(deg_targets_probe)
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe,deg_targets_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_NoObs_Targets.EMM.A069 <-
  IRE.subj.EMM %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_grid(~target) + 
  
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
    data = dplyr::filter(IRE.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM, probe == "without obstacle"),
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
  ) +
  ggtitle("without obstacle") 

g.IRE_NoObs_Targets.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_Obs_Targets.EMM.A069 <-
  IRE.subj.EMM %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_grid(~target) + 
  
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
    data = dplyr::filter(IRE.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM, probe == "with obstacle"),
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
    plot.title = element_text(size = 12, face = "bold")
  ) +
  ggtitle("with obstacle") 

g.IRE_Obs_Targets.EMM.A069


#merge into a 1x2 plot grid (library patchwork)
g.IRE_Targets.EMM.A069 <-
  (g.IRE_Obs_Targets.EMM.A069 + g.IRE_NoObs_Targets.EMM.A069) +
  plot_layout(
    heights = c(1,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.IRE_Targets.EMM.A069




# ==== INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# we compare between Movement No obstacle prime and obstacle prime, separately for stop-signal (EX, MI)
IRE.Contrasts_TL           <- data.frame(matrix(ncol = 16, nrow = 40000))
colnames(IRE.Contrasts_TL) <- c( "EX 22.5",  "MI 22.5",
                                 "EX 67.5",  "MI 67.5",
                                 "EX 112.5", "MI 112.5", 
                                 "EX 157.5", "MI 157.5",
                                 "EX 202.5", "MI 202.5",
                                 "EX 247.5", "MI 247.5",
                                 "EX 292.5", "MI 292.5",
                                 "EX 337.5", "MI 337.5")

IRE.Contrasts_TL$`EX 22.5`      <- (IRE_Probe_NoObs.posteriors$`yes_go_22.5`    + IRE_Probe_Obs.posteriors$`yes_go_22.5`)/2     - (IRE_Probe_NoObs.posteriors$`no_go_22.5`    + IRE_Probe_Obs.posteriors$`no_go_22.5`)/2
IRE.Contrasts_TL$`MI 22.5`      <- (IRE_Probe_NoObs.posteriors$`yes_stop_22.5`  + IRE_Probe_Obs.posteriors$`yes_stop_22.5`)/2   - (IRE_Probe_NoObs.posteriors$`no_stop_22.5`  + IRE_Probe_Obs.posteriors$`no_stop_22.5`)/2
IRE.Contrasts_TL$`EX 67.5`      <- (IRE_Probe_NoObs.posteriors$`yes_go_67.5`    + IRE_Probe_Obs.posteriors$`yes_go_67.5`)/2     - (IRE_Probe_NoObs.posteriors$`no_go_67.5`    + IRE_Probe_Obs.posteriors$`no_go_67.5`)/2
IRE.Contrasts_TL$`MI 67.5`      <- (IRE_Probe_NoObs.posteriors$`yes_stop_67.5`  + IRE_Probe_Obs.posteriors$`yes_stop_67.5`)/2   - (IRE_Probe_NoObs.posteriors$`no_stop_67.5`  + IRE_Probe_Obs.posteriors$`no_stop_67.5`)/2
IRE.Contrasts_TL$`EX 112.5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_112.5`   + IRE_Probe_Obs.posteriors$`yes_go_112.5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_112.5`   + IRE_Probe_Obs.posteriors$`no_go_112.5`)/2
IRE.Contrasts_TL$`MI 112.5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_112.5` + IRE_Probe_Obs.posteriors$`yes_stop_112.5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_112.5` + IRE_Probe_Obs.posteriors$`no_stop_112.5`)/2
IRE.Contrasts_TL$`EX 157.5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_157.5`   + IRE_Probe_Obs.posteriors$`yes_go_157.5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_157.5`   + IRE_Probe_Obs.posteriors$`no_go_157.5`)/2
IRE.Contrasts_TL$`MI 157.5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_157.5` + IRE_Probe_Obs.posteriors$`yes_stop_157.5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_157.5` + IRE_Probe_Obs.posteriors$`no_stop_157.5`)/2
IRE.Contrasts_TL$`EX 202.5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_202.5`   + IRE_Probe_Obs.posteriors$`yes_go_202.5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_202.5`   + IRE_Probe_Obs.posteriors$`no_go_202.5`)/2
IRE.Contrasts_TL$`MI 202.5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_202.5` + IRE_Probe_Obs.posteriors$`yes_stop_202.5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_202.5` + IRE_Probe_Obs.posteriors$`no_stop_202.5`)/2
IRE.Contrasts_TL$`EX 247.5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_247.5`   + IRE_Probe_Obs.posteriors$`yes_go_247.5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_247.5`   + IRE_Probe_Obs.posteriors$`no_go_247.5`)/2
IRE.Contrasts_TL$`MI 247.5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_247.5` + IRE_Probe_Obs.posteriors$`yes_stop_247.5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_247.5` + IRE_Probe_Obs.posteriors$`no_stop_247.5`)/2
IRE.Contrasts_TL$`EX 292.5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_292.5`   + IRE_Probe_Obs.posteriors$`yes_go_292.5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_292.5`   + IRE_Probe_Obs.posteriors$`no_go_292.5`)/2
IRE.Contrasts_TL$`MI 292.5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_292.5` + IRE_Probe_Obs.posteriors$`yes_stop_292.5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_292.5` + IRE_Probe_Obs.posteriors$`no_stop_292.5`)/2
IRE.Contrasts_TL$`EX 337.5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_337.5`   + IRE_Probe_Obs.posteriors$`yes_go_337.5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_337.5`   + IRE_Probe_Obs.posteriors$`no_go_337.5`)/2
IRE.Contrasts_TL$`MI 337.5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_337.5` + IRE_Probe_Obs.posteriors$`yes_stop_337.5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_337.5` + IRE_Probe_Obs.posteriors$`no_stop_337.5`)/2


IRE.Contrasts_TL_long <- pivot_longer(IRE.Contrasts_TL, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c( "EX 22.5",  "MI 22.5",
                     "EX 67.5",  "MI 67.5",
                     "EX 112.5", "MI 112.5", 
                     "EX 157.5", "MI 157.5",
                     "EX 202.5", "MI 202.5",
                     "EX 247.5", "MI 247.5",
                     "EX 292.5", "MI 292.5",
                     "EX 337.5", "MI 337.5")

# Convert 'contrast' to a factor with this order
IRE.Contrasts_TL_long$contrast <- factor(IRE.Contrasts_TL_long$contrast, levels = contrast_order)

# Check
head(IRE.Contrasts_TL_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
IRE.contrasts_TL.summary <-
  IRE.Contrasts_TL_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

IRE.contrasts_TL.summary$pd <- format(IRE.contrasts_TL.summary$pd, nsmall = 4)  
print(IRE.contrasts_TL.summary, n = Inf, width = Inf)


#### Calculate ROPE
### Calculate and back-transform ROPE, this is a bit tricky
## this does not take into account that values are on log scale so don't use it
## apparently, back-transformation is now implemented :)
#RR_NoObs <- rope_range(A069_fit.IRE_Probe_NoObs_OPRIxSS) 
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

( RR_NoObs <- c(-rope_value_RE_log, rope_value_RE_log) )
( RR_Obs   <- rope_range(A069_fit.IRE_Probe_Obs_OPRIxSSxTL) )
( RR       <- (RR_NoObs + RR_Obs) /2 )



# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE                       <- as.data.frame(IRE.contrasts_TL.summary)
IRE.contrast_in_ROPE$lowerROPE             <- RR[1]
IRE.contrast_in_ROPE$upperROPE             <- RR[2]
IRE.contrast_in_ROPE$CI_range              <- IRE.contrast_in_ROPE$upper - IRE.contrast_in_ROPE$lower
IRE.contrast_in_ROPE$minUpper              <- IRE.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
IRE.contrast_in_ROPE$maxLower              <-  IRE.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
IRE.contrast_in_ROPE$DiffminUppermaxLower  <- IRE.contrast_in_ROPE$minUpper  - IRE.contrast_in_ROPE$maxLower 
IRE.contrast_in_ROPE$Zeros                 <- rep(0,nrow(IRE.contrast_in_ROPE))
IRE.contrast_in_ROPE$Overlap               <- IRE.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
IRE.contrast_in_ROPE$perc_in_ROPE          <- (IRE.contrast_in_ROPE$Overlap*100)/IRE.contrast_in_ROPE$CI_range
IRE.contrast_in_ROPE[,c(1:7,14)]


#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS FOR EACH TARGET POOLED OVER OBSTACLE PROBE (SUPPLEMENTARY FIGURE S7) ====
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL22.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 22.5" | contrast=="MI 22.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 22.5" = "",
    "MI 22.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 22.5" | contrast=="MI 22.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("22.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL22.5.A069


titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL67.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 67.5" | contrast=="MI 67.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 67.5" = "",
    "MI 67.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 67.5" | contrast=="MI 67.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("67.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL67.5.A069


titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL112.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 112.5" | contrast=="MI 112.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 112.5" = "",
    "MI 112.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 112.5" | contrast=="MI 112.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("112.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL112.5.A069



titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL157.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 157.5" | contrast=="MI 157.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 157.5" = "",
    "MI 157.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 157.5" | contrast=="MI 157.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("157.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL157.5.A069



titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL202.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 202.5" | contrast=="MI 202.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 202.5" = "",
    "MI 202.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 202.5" | contrast=="MI 202.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("202.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL202.5.A069



titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL247.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 247.5" | contrast=="MI 247.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 247.5" = "",
    "MI 247.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 247.5" | contrast=="MI 247.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("247.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL247.5.A069



titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL292.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 292.5" | contrast=="MI 292.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 292.5" = "",
    "MI 292.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 292.5" | contrast=="MI 292.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("292.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL292.5.A069



titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_TL337.5.A069 <- IRE.Contrasts_TL_long %>% filter(contrast=="EX 337.5" | contrast=="MI 337.5") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 337.5" = "",
    "MI 337.5" = ""
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
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_TL.summary, contrast == "EX 337.5" | contrast=="MI 337.5" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("337.5°") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_TL337.5.A069





## arrange in circle
plots <- c(g.IRE.Contrast_TL22.5.A069,g.IRE.Contrast_TL67.5.A069,g.IRE.Contrast_TL112.5.A069,g.IRE.Contrast_TL157.5.A069,g.IRE.Contrast_TL202.5.A069,g.IRE.Contrast_TL247.5.A069,g.IRE.Contrast_TL292.5.A069,g.IRE.Contrast_TL337.5.A069)

#### Parameters for the circular layout 
n <- length(plots)
radius <- 0.35     # distance from center (0..0.5 recommended)
pw <- 0.28         # plot width (tweak as needed)
ph <- 0.28         # plot height (tweak as needed)
start_angle <-(pi/180)*22.5# pi/2  # start at 22.5°; change to rotate all positions

# compute positions around the circle (centered in [0,1] x [0,1])
angles <- seq(0, 2*pi, length.out = n + 1)[- (n+1)] + start_angle
xs <- 0.5 + radius * cos(angles)
ys <- 0.5 + radius * sin(angles)

#Compose the final canvas and draw each plot at its (x,y
g.IRE.Contrast_TL <- ggdraw()  # blank canvas with coordinates 0..1
for (i in seq_len(n)) {
  # draw_plot places plots with x,y = lower-left corner, so shift by half width/height
  g.IRE.Contrast_TL <- g.IRE.Contrast_TL + draw_plot(plots[[i]],
                                     x = xs[i] - pw/2,
                                     y = ys[i] - ph/2,
                                     width = pw,
                                     height = ph)
}
g.IRE.Contrast_TL


# create task setup plot
# parameters
center <- c(x = 0, y = 0)
start_radius  <- 0.1        # radius of central "start" circle
target_radius <- 0.1       # radius of each target circle
ring_radius <- 0.8          # distance from center to center of target circles
first_angle_deg <- 22.5     # first target angle in degrees
n_targets <- 8
angles_deg <- first_angle_deg + seq(0, by = 45, length.out = n_targets)
angles_rad <- angles_deg * pi / 180

# build data frame for targets
targets <- data.frame(
  id = seq_len(n_targets),
  angle_deg = angles_deg,
  x = center["x"] + ring_radius * cos(angles_rad),
  y = center["y"] + ring_radius * sin(angles_rad)
)

# center circle data (ggforce::geom_circle needs center + r)
center_df <- data.frame(x0 = center["x"], y0 = center["y"], r = start_radius)

# targets as circles
targets_df <- data.frame(x0 = targets$x, y0 = targets$y, r = target_radius, id = targets$id)

# Plot
g.tasksetup <- ggplot() +
  # target circles
  ggforce::geom_circle(data = targets_df, aes(x0 = x0, y0 = y0, r = r, group = id),
                       fill = "white", alpha = 0.8, color = "white", size = 0.5) +
  # center/start circle
  ggforce::geom_circle(data = center_df, aes(x0 = x0, y0 = y0, r = r),
                       fill = "white", alpha = 0.9, color = "white", size = 0.5) +
  # optional: center point and target labels
  #geom_point(aes(x = center["x"], y = center["y"]), color = "black", size = 0.5) +
  geom_text(data = targets, aes(x = x, y = y, label = c("22.5","67.5","112.5","157.5","202.5","247.5","292.5","337.5")), vjust = -1.1, size = 3) +
  # equal aspect so circles are not ellipses
  coord_equal() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "grey75", color = NA),
    panel.background = element_rect(fill = "grey75", color = NA))+ 
  #ggtitle("Center start circle with 8 targets (first at 22.5°)") +
  xlim(- (ring_radius + target_radius + 0.2), ring_radius + target_radius + 0.2) +
  ylim(- (ring_radius + target_radius + 0.2), ring_radius + target_radius + 0.2)

g.tasksetup

# add plot to center
g.IRE.Contrast_TL.A069 <- g.IRE.Contrast_TL + draw_plot(g.tasksetup, x = 0.5 - pw/2, y = 0.5 - ph/2,
                                   width = pw, height = ph)
g.IRE.Contrast_TL.A069


# Build a valid path
outfile <- file.path(figurePath, "FigS7_A069_InitialReachError_Contrast_TargetsAll_PooledOverProbeObs_TargetLayout.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.IRE.Contrast_TL.A069,
  width = 24, height = 24, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS7_A069_InitialReachError_Contrast_TargetsAll_PooledOverProbeObs_TargetLayout.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.IRE.Contrast_TL.A069,
  width = 24, height = 24, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "FigS7_A069_InitialReachError_Contrast_TargetsAll_PooledOverProbeObs_TargetLayout.png")

ggsave(
  filename = outfile,
  plot = g.IRE.Contrast_TL.A069,
  width = 28, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)

()




#  ==== SUMMARY TABLE (SUPPLEMENTARY TABLE S2) ====
RT.TL2_pd   <- p_direction(IRE.Contrasts_TL)
table.RE.TL <- IRE.contrast_in_ROPE[,c(1:4,14)]
table.RE.TL <- cbind(table.RE.TL,RT.TL2_pd[1])

# Reorder rows to match the Parameter column
table.RE.TL    <- table.RE.TL[match(table.RE.TL$Parameter, table.RE.TL$contrast), ]
table.RE.TL    <- cbind(table.RE.TL,RT.TL2_pd[2]) #add pd values
table.RE.TL    <- table.RE.TL %>% select(-c(Parameter)) # remove Parameter column
table.RE.TL$pd <- table.RE.TL$pd *100
table.RE.TL <- table.RE.TL %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
table.RE.TL <- table.RE.TL %>% rename("Bias" = median, "% in ROPE" = perc_in_ROPE, "pd (%)" = pd) %>% unite("HDI",lower:upper)
table.RE.TL <- table.RE.TL  %>% mutate(HDI = str_c("[", HDI, "]"))
table.RE.TL <- table.RE.TL %>%
  mutate(HDI = str_replace(HDI, "_", ","))
table.RE.TL <- table.RE.TL %>%
  separate(col = contrast, into = c("Trial Type", "Target Location (°)"), sep = " ")
table.RE.TL   <- table.RE.TL %>% mutate( `Trial Type`    = recode_factor( as.factor(`Trial Type` ) , "EX" = "Execution", "MI" = "Motor Imagery"))


ftable.IRE.TL <- flextable(table.RE.TL)  %>%  
  fontsize(size = 12, part = "header") %>% 
  bold(part = "header") %>% 
  align(align = "center", part = "header") %>% 
  align(j = c("Target Location (°)", "Bias", "HDI", "% in ROPE","pd (%)"), align = "center", part = "all") %>% 
  #bg(j="r",bg = scales::col_numeric(palette = "viridis", domain = c(-0.3, 0.3))) %>% 
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 2),part = "header" ) %>% 
  hline_top(border = fp_border(color = "black", width = 2) ) %>% 
  hline_bottom(border = fp_border(color = "black", width = 2) ) %>% 
  hline(i = c(2,4,6,8,10,12,14), border = fp_border(color = "grey85", width = 1.5)) %>% 
  set_caption("Table S2: Reach Bias as a function of trial type and target location in Experiment 3") 

ftable.IRE.TL <- autofit(ftable.IRE.TL)
ftable.IRE.TL

# export to doc
doc <- read_docx()
doc <- body_add_flextable(doc, ftable.IRE.TL)
print(doc, target = "TableS2_A069_Table_ReachBiasTargetLocation.docx")




#==== INITIAL REACH ERROR PROBE WITH BLOCKS ====
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe



# ==== MODEL FITTING: INITIAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (UB=180) LOGNORMAL DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")
dClean_noobs$blocks_thisN <- as.factor(dClean_noobs$blocks_thisN)

dClean_noobs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime:blocks_thisN) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_noobs$ReachDiff2_Probe)
sd(log(dClean_noobs$ReachDiff2_Probe))

# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime * blocks_thisN + ( 1 | subID )),
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
                           prior(exponential(0.8),               class = "sigma")      
)      


# we model the mu here, we include random effects for participants
# runs about 45 min
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK <- brm( 
  bf(ReachDiff2_Probe | trunc(ub = 180) ~ obstacle_prime * stop_signal_prime * blocks_thisN + ( 1 | subID )),     # model specification
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
  file = paste0(analysisPath,'/', 'A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK'), #save model
)
endTime <- Sys.time()
( endTime- startTime )










#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK,          
     combo = c("dens", "trace"),
     ask   = FALSE
)


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime),
                        blocks_thisN       = levels(dClean_noobs$blocks_thisN))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS_BL = interaction(obstacle_prime, stop_signal_prime, blocks_thisN, sep = "_", drop = TRUE))
colnames(IRE_Probe_NoObs.posteriors) <- exp.cond_posterior$OPRI_SS_BL


# Long format (use pivot_longer)
IRE_Probe_NoObs.posteriors_long <- IRE_Probe_NoObs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS_BL", values_to = "value") %>%
  separate(OPRI_SS_BL, into = c("obstacle_prime", "stop_signal_prime", "blocks_thisN"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_NoObs.EMM <- IRE_Probe_NoObs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime, blocks_thisN) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_NoObs.EMM



# Conditions
dClean_noobs$subID <- droplevels(dClean_noobs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_noobs$subID),
                                  obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_noobs$stop_signal_prime),
                                  blocks_thisN       = levels(dClean_noobs$blocks_thisN))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A069_fit.IRE_Probe_NoObs_OPRIxSSxBLOCK,
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
  group_by(subID, obstacle_prime, stop_signal_prime, blocks_thisN) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)




# ==== MODEL FITTING: INITIAL REACH ERROR OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=180) GAUSSIAN DISTRIBUTION ====
dClean_obs <- dClean %>% filter(obstacle_probe=="yes")
dClean_obs$blocks_thisN <- as.factor(dClean_obs$blocks_thisN)

dClean_obs %>%
  ggplot( aes(x = ReachDiff2_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime:blocks_thisN) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(dClean_obs$ReachDiff2_Probe)

# looking at prior values
get_prior(bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime * blocks_thisN + ( 1 | subID )),
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
                          prior(student_t(3, 0,  15   ), class = "sigma"))


# we model the mu here. We include random effects for participants
# runs about 1.5 hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A069_fit.IRE_Probe_Obs_OPRIxSSxBLOCK <- brm( 
  bf(ReachDiff2_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * stop_signal_prime * blocks_thisN + ( 1 | subID )),     # model specification
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
  file = paste0(analysisPath,'/', 'A069_fit.IRE_Probe_Obs_OPRIxSSxBLOCK'), #save model
)
endTime <- Sys.time()
( endTime- startTime )





#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime),
                        blocks_thisN       = levels(dClean_obs$blocks_thisN))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A069_fit.IRE_Probe_Obs_OPRIxSSxBLOCK,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS_BL = interaction(obstacle_prime, stop_signal_prime, blocks_thisN, sep = "_", drop = TRUE))
colnames(IRE_Probe_Obs.posteriors) <- exp.cond_posterior$OPRI_SS_BL


# Long format (use pivot_longer)
IRE_Probe_Obs.posteriors_long <- IRE_Probe_Obs.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS_BL", values_to = "value") %>%
  separate(OPRI_SS_BL, into = c("obstacle_prime", "stop_signal_prime", "blocks_thisN"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE_Obs.EMM <- IRE_Probe_Obs.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime, blocks_thisN) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE_Obs.EMM



# Conditions
dClean_obs$subID <- droplevels(dClean_obs$subID)
exp.cond.subj      <- expand.grid(subID              = levels(dClean_obs$subID),
                                  obstacle_prime     = levels(dClean_obs$obstacle_prime),
                                  stop_signal_prime  = levels(dClean_obs$stop_signal_prime),
                                  blocks_thisN   = levels(dClean_obs$blocks_thisN))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A069_fit.IRE_Probe_Obs_OPRIxSSxBLOCK,
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
  group_by(subID, obstacle_prime, stop_signal_prime, blocks_thisN) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)












# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) ====
# Prepare and tidy data
IRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
IRE.EMM <- rbind(IRE_NoObs.EMM,IRE_Obs.EMM)

IRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
IRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
IRE.subj.EMM <- rbind(IRE_NoObs.subj.EMM,IRE_Obs.subj.EMM)

# rename factors
IRE.EMM <- IRE.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "EX", "MI")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    block       = factor(blocks_thisN)
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe,blocks_thisN)) # remove old columns

print(IRE.EMM, n = Inf)

# rename factors
IRE.subj.EMM <- IRE.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "EX", "MI")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    block       = factor(blocks_thisN)
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe,blocks_thisN)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_NoObs_Blocks.EMM.A069 <-
  IRE.subj.EMM %>%
  dplyr::filter(probe == "without obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_grid(~block) + 
  
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
    data = dplyr::filter(IRE.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM, probe == "without obstacle"),
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
  ) +
  ggtitle("without obstacle") 

g.IRE_NoObs_Blocks.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_Obs_Blocks.EMM.A069 <-
  IRE.subj.EMM %>%
  dplyr::filter(probe == "with obstacle") %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
  facet_grid(~block) + 
  
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
    data = dplyr::filter(IRE.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(IRE.EMM, probe == "with obstacle"),
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
    plot.title = element_text(size = 12, face = "bold")
  ) +
  ggtitle("with obstacle") 

g.IRE_Obs_Blocks.EMM.A069


#merge into a 1x2 plot grid (library patchwork)
g.IRE_Blocks.EMM.A069 <-
  (g.IRE_Obs_Blocks.EMM.A069 + g.IRE_NoObs_Blocks.EMM.A069) +
  plot_layout(
    heights = c(1,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.IRE_Blocks.EMM.A069




# ==== INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# we compare between Movement No obstacle prime and obstacle prime, separately for stop-signal (EX, MI)
IRE.Contrasts_BL           <- data.frame(matrix(ncol = 20, nrow = 40000))
colnames(IRE.Contrasts_BL) <- c( "EX 1", "MI 1",
                                 "EX 2", "MI 2",
                                 "EX 3", "MI 3", 
                                 "EX 4", "MI 4",
                                 "EX 5", "MI 5",
                                 "EX 6", "MI 6",
                                 "EX 7", "MI 7",
                                 "EX 8", "MI 8",
                                 "EX 9", "MI 9",
                                 "EX 10","MI 10")

IRE.Contrasts_BL$`EX 1`      <- (IRE_Probe_NoObs.posteriors$`yes_go_1`    + IRE_Probe_Obs.posteriors$`yes_go_1`)/2     - (IRE_Probe_NoObs.posteriors$`no_go_1`    + IRE_Probe_Obs.posteriors$`no_go_1`)/2
IRE.Contrasts_BL$`MI 1`      <- (IRE_Probe_NoObs.posteriors$`yes_stop_1`  + IRE_Probe_Obs.posteriors$`yes_stop_1`)/2   - (IRE_Probe_NoObs.posteriors$`no_stop_1`  + IRE_Probe_Obs.posteriors$`no_stop_1`)/2
IRE.Contrasts_BL$`EX 2`      <- (IRE_Probe_NoObs.posteriors$`yes_go_2`    + IRE_Probe_Obs.posteriors$`yes_go_2`)/2     - (IRE_Probe_NoObs.posteriors$`no_go_2`    + IRE_Probe_Obs.posteriors$`no_go_2`)/2
IRE.Contrasts_BL$`MI 2`      <- (IRE_Probe_NoObs.posteriors$`yes_stop_2`  + IRE_Probe_Obs.posteriors$`yes_stop_2`)/2   - (IRE_Probe_NoObs.posteriors$`no_stop_2`  + IRE_Probe_Obs.posteriors$`no_stop_2`)/2
IRE.Contrasts_BL$`EX 3`     <- (IRE_Probe_NoObs.posteriors$`yes_go_3`   + IRE_Probe_Obs.posteriors$`yes_go_3`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_3`   + IRE_Probe_Obs.posteriors$`no_go_3`)/2
IRE.Contrasts_BL$`MI 3`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_3` + IRE_Probe_Obs.posteriors$`yes_stop_3`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_3` + IRE_Probe_Obs.posteriors$`no_stop_3`)/2
IRE.Contrasts_BL$`EX 4`     <- (IRE_Probe_NoObs.posteriors$`yes_go_4`   + IRE_Probe_Obs.posteriors$`yes_go_4`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_4`   + IRE_Probe_Obs.posteriors$`no_go_4`)/2
IRE.Contrasts_BL$`MI 4`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_4` + IRE_Probe_Obs.posteriors$`yes_stop_4`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_4` + IRE_Probe_Obs.posteriors$`no_stop_4`)/2
IRE.Contrasts_BL$`EX 5`     <- (IRE_Probe_NoObs.posteriors$`yes_go_5`   + IRE_Probe_Obs.posteriors$`yes_go_5`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_5`   + IRE_Probe_Obs.posteriors$`no_go_5`)/2
IRE.Contrasts_BL$`MI 5`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_5` + IRE_Probe_Obs.posteriors$`yes_stop_5`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_5` + IRE_Probe_Obs.posteriors$`no_stop_5`)/2
IRE.Contrasts_BL$`EX 6`     <- (IRE_Probe_NoObs.posteriors$`yes_go_6`   + IRE_Probe_Obs.posteriors$`yes_go_6`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_6`   + IRE_Probe_Obs.posteriors$`no_go_6`)/2
IRE.Contrasts_BL$`MI 6`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_6` + IRE_Probe_Obs.posteriors$`yes_stop_6`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_6` + IRE_Probe_Obs.posteriors$`no_stop_6`)/2
IRE.Contrasts_BL$`EX 7`     <- (IRE_Probe_NoObs.posteriors$`yes_go_7`   + IRE_Probe_Obs.posteriors$`yes_go_7`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_7`   + IRE_Probe_Obs.posteriors$`no_go_7`)/2
IRE.Contrasts_BL$`MI 7`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_7` + IRE_Probe_Obs.posteriors$`yes_stop_7`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_7` + IRE_Probe_Obs.posteriors$`no_stop_7`)/2
IRE.Contrasts_BL$`EX 8`     <- (IRE_Probe_NoObs.posteriors$`yes_go_8`   + IRE_Probe_Obs.posteriors$`yes_go_8`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_8`   + IRE_Probe_Obs.posteriors$`no_go_8`)/2
IRE.Contrasts_BL$`MI 8`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_8` + IRE_Probe_Obs.posteriors$`yes_stop_8`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_8` + IRE_Probe_Obs.posteriors$`no_stop_8`)/2
IRE.Contrasts_BL$`EX 9`     <- (IRE_Probe_NoObs.posteriors$`yes_go_9`   + IRE_Probe_Obs.posteriors$`yes_go_9`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_9`   + IRE_Probe_Obs.posteriors$`no_go_9`)/2
IRE.Contrasts_BL$`MI 9`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_9` + IRE_Probe_Obs.posteriors$`yes_stop_9`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_9` + IRE_Probe_Obs.posteriors$`no_stop_9`)/2
IRE.Contrasts_BL$`EX 10`     <- (IRE_Probe_NoObs.posteriors$`yes_go_10`   + IRE_Probe_Obs.posteriors$`yes_go_10`)/2    - (IRE_Probe_NoObs.posteriors$`no_go_10`   + IRE_Probe_Obs.posteriors$`no_go_10`)/2
IRE.Contrasts_BL$`MI 10`     <- (IRE_Probe_NoObs.posteriors$`yes_stop_10` + IRE_Probe_Obs.posteriors$`yes_stop_10`)/2  - (IRE_Probe_NoObs.posteriors$`no_stop_10` + IRE_Probe_Obs.posteriors$`no_stop_10`)/2

IRE.Contrasts_BL_long <- pivot_longer(IRE.Contrasts_BL, cols = everything(),
                                      names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c( "EX 1", "MI 1",
                     "EX 2", "MI 2",
                     "EX 3", "MI 3", 
                     "EX 4", "MI 4",
                     "EX 5", "MI 5",
                     "EX 6", "MI 6",
                     "EX 7", "MI 7",
                     "EX 8", "MI 8",
                     "EX 9", "MI 9",
                     "EX 10","MI 10")

# Convert 'contrast' to a factor with this order
IRE.Contrasts_BL_long$contrast <- factor(IRE.Contrasts_BL_long$contrast, levels = contrast_order)

# Check
head(IRE.Contrasts_BL_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
IRE.contrasts_BL.summary <-
  IRE.Contrasts_BL_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

IRE.contrasts_BL.summary$pd <- format(IRE.contrasts_BL.summary$pd, nsmall = 4)  
print(IRE.contrasts_BL.summary, n = Inf, width = Inf)


#### Calculate ROPE
### Calculate and back-transform ROPE, this is a bit tricky
## this does not take into account that values are on log scale so don't use it
## apparently, back-transformation is now implemented :)
#RR_NoObs <- rope_range(A069_fit.IRE_Probe_NoObs_OPRIxSS) 
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

( RR_NoObs <- c(-rope_value_RE_log, rope_value_RE_log) )
( RR_Obs   <- rope_range(A069_fit.IRE_Probe_Obs_OPRIxSSxBLOCK) )
( RR       <- (RR_NoObs + RR_Obs) /2 )



# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE                       <- as.data.frame(IRE.contrasts_BL.summary)
IRE.contrast_in_ROPE$lowerROPE             <- RR[1]
IRE.contrast_in_ROPE$upperROPE             <- RR[2]
IRE.contrast_in_ROPE$CI_range              <- IRE.contrast_in_ROPE$upper - IRE.contrast_in_ROPE$lower
IRE.contrast_in_ROPE$minUpper              <- IRE.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
IRE.contrast_in_ROPE$maxLower              <-  IRE.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
IRE.contrast_in_ROPE$DiffminUppermaxLower  <- IRE.contrast_in_ROPE$minUpper  - IRE.contrast_in_ROPE$maxLower 
IRE.contrast_in_ROPE$Zeros                 <- rep(0,nrow(IRE.contrast_in_ROPE))
IRE.contrast_in_ROPE$Overlap               <- IRE.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
IRE.contrast_in_ROPE$perc_in_ROPE          <- (IRE.contrast_in_ROPE$Overlap*100)/IRE.contrast_in_ROPE$CI_range
IRE.contrast_in_ROPE[,c(1:7,14)]


#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS FOR EACH TARGET POOLED OVER OBSTACLE PROBE (SUPPLEMENTARY FIGURE S8) ====
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_BL_EX.A069 <- IRE.Contrasts_BL_long %>% filter(contrast=="EX 1" | contrast=="EX 2" | contrast=="EX 3" | contrast=="EX 4" | contrast=="EX 5" | contrast=="EX 6" | contrast=="EX 7" | contrast=="EX 8" | contrast=="EX 9" | contrast=="EX 10") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX 1" = "Block 01", "EX 2" = "Block 02", "EX 3" = "Block 03", "EX 4" = "Block 04", "EX 5" = "Block 05", "EX 6" = "Block 06", "EX 7" = "Block 07", "EX 8" = "Block 08", "EX 9" = "Block 09", "EX 10" = "Block 10"
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
  scale_fill_manual(values = c(color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1)) + 
  scale_pattern_fill_manual(values = c(color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_BL.summary, contrast=="EX 1" | contrast=="EX 2" | contrast=="EX 3" | contrast=="EX 4" | contrast=="EX 5" | contrast=="EX 6" | contrast=="EX 7" | contrast=="EX 8" | contrast=="EX 9" | contrast=="EX 10" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1,color_exe1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("Execution") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_BL_EX.A069


titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast_BL_MI.A069 <- IRE.Contrasts_BL_long %>% filter(contrast=="MI 1" | contrast=="MI 2" | contrast=="MI 3" | contrast=="MI 4" | contrast=="MI 5" | contrast=="MI 6" | contrast=="MI 7" | contrast=="MI 8" | contrast=="MI 9" | contrast=="MI 10") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "MI 1" = "Block 01", "MI 2" = "Block 02", "MI 3" = "Block 03", "MI 4" = "Block 04", "MI 5" = "Block 05", "MI 6" = "Block 06", "MI 7" = "Block 07", "MI 8" = "Block 08", "MI 9" = "Block 09", "MI 10" = "Block 10"
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
  scale_fill_manual(values = c(color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1)) + 
  scale_pattern_fill_manual(values = c(color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts_BL.summary, contrast=="MI 1" | contrast=="MI 2" | contrast=="MI 3" | contrast=="MI 4" | contrast=="MI 5" | contrast=="MI 6" | contrast=="MI 7" | contrast=="MI 8" | contrast=="MI 9" | contrast=="MI 10" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-3,20),breaks=c(0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("Motor Imagery") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast_BL_MI.A069


#merge into a 1x2 plot grid (library patchwork)
g.IRE.Blocks.EMM.A069 <-
  (g.IRE.Contrast_BL_EX.A069 + g.IRE.Contrast_BL_MI.A069) +
  plot_layout(
    widths = c(1, 1.2),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.IRE.Blocks.EMM.A069




# Build a valid path
outfile <- file.path(figurePath, "FigS8_A069_InitialReachError_Contrast_BlocksAll_PooledOverProbeObs.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.IRE.Blocks.EMM.A069,
  width = 18, height = 28, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS8_A069_InitialReachError_Contrast_BlocksAll_PooledOverProbeObs.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.IRE.Blocks.EMM.A069,
  width = 18, height = 28, units = "cm",
  device = "svg"
)

outfile <- file.path(figurePath, "FigS8_A069_InitialReachError_Contrast_BlocksAll_PooledOverProbeObs.png")

ggsave(
  filename = outfile,
  plot = g.IRE.Blocks.EMM.A069,
  width = 18, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)








#  ==== SUMMARY TABLE (SUPPLEMENTARY TABLE S3) ====
IRE.BL_pd   <- p_direction(IRE.Contrasts_BL)
table.IRE.BL <- IRE.contrast_in_ROPE[,c(1:4,14)]
table.IRE.BL <- cbind(table.IRE.BL,IRE.BL_pd[1])

# Reorder rows to match the Parameter column
table.IRE.BL    <- table.IRE.BL[match(table.IRE.BL$Parameter, table.IRE.BL$contrast), ]
table.IRE.BL    <- cbind(table.IRE.BL,IRE.BL_pd[2]) #add pd values
table.IRE.BL    <- table.IRE.BL %>% select(-c(Parameter)) # remove Parameter column
table.IRE.BL$pd <- table.IRE.BL$pd *100
table.IRE.BL <- table.IRE.BL %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
table.IRE.BL <- table.IRE.BL %>% rename("Bias" = median, "% in ROPE" = perc_in_ROPE, "pd (%)" = pd) %>% unite("HDI",lower:upper)
table.IRE.BL <- table.IRE.BL  %>% mutate(HDI = str_c("[", HDI, "]"))
table.IRE.BL <- table.IRE.BL %>%
  mutate(HDI = str_replace(HDI, "_", ","))
table.IRE.BL <- table.IRE.BL %>%
  separate(col = contrast, into = c("Trial Type", "Block Number"), sep = " ")
table.IRE.BL   <- table.IRE.BL %>% mutate( `Trial Type`    = recode_factor( as.factor(`Trial Type` ) , "EX" = "Execution", "MI" = "Motor Imagery"))


ftable.IRE.BL <- flextable(table.IRE.BL)  %>%  
  fontsize(size = 12, part = "header") %>% 
  bold(part = "header") %>% 
  align(align = "center", part = "header") %>% 
  align(j = c("Block Number", "Bias", "HDI", "% in ROPE","pd (%)"), align = "center", part = "all") %>% 
  #bg(j="r",bg = scales::col_numeric(palette = "viridis", domain = c(-0.3, 0.3))) %>% 
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 2),part = "header" ) %>% 
  hline_top(border = fp_border(color = "black", width = 2) ) %>% 
  hline_bottom(border = fp_border(color = "black", width = 2) ) %>% 
  hline(i = c(2,4,6,8,10,12,14), border = fp_border(color = "grey85", width = 1.5)) %>% 
  set_caption("Table S3. Reach Bias as a function of trial type and block number in Experiment 3") 

ftable.IRE.BL <- autofit(ftable.IRE.BL)
ftable.IRE.BL

# export to doc
doc <- read_docx()
doc <- body_add_flextable(doc, ftable.IRE.BL)
print(doc, target = "TableS3_A069_Table_ReachBiasBlockNumber.docx")






#==== REACTION TIME (RT) PROBE ====
# ==== MODEL FITTING: REACTION TIME: 3-FACTOR INTERACTION MODEL WITH SHIFTED LOG-NORMAL DISTRIBUTION ====
dClean %>%
  ggplot( aes(x = RT_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 10, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "RT Probe" ) +
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
A069_fit.RT_Probe_OPRIxOPROxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.RT_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )





#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.RT_Probe_OPRIxOPROxSS,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.RT_Probe_OPRIxOPROxSS,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.RT_Probe_OPRIxOPROxSS,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.RT_Probe_OPRIxOPROxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.RT_Probe_OPRIxOPROxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)




#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean$obstacle_prime),
                        obstacle_probe     = levels(dClean$obstacle_probe),
                        stop_signal_prime  = levels(dClean$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
RT_Probe.posteriors <- as.data.frame(fitted(
  A069_fit.RT_Probe_OPRIxOPROxSS,
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
                               obstacle_probe     = levels(dClean$obstacle_probe),
                               stop_signal_prime  = levels(dClean$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A069_fit.RT_Probe_OPRIxOPROxSS,
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
RT.EMM <- RT.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

RT.EMM

# rename factors
RT.subj.EMM <- RT.subj.EMM  %>% 
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

g.RT_NoObs.EMM.A069 <-
  RT.subj.EMM %>%
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
    data = dplyr::filter(RT.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM, probe == "without obstacle"),
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

g.RT_NoObs.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Reaction time (ms)"

g.RT_Obs.EMM.A069 <-
  RT.subj.EMM %>%
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
    data = dplyr::filter(RT.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM, probe == "with obstacle"),
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

g.RT_Obs.EMM.A069

#merge into a 1x2 plot grid (library patchwork)
g.RT.EMM.A069 <-
  (g.RT_Obs.EMM.A069 + g.RT_NoObs.EMM.A069) +
  plot_layout(
    widths = c(1, 1.2),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.RT.EMM.A069












# ==== REACTION TIME: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
RT.Contrasts           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(RT.Contrasts) <- c("Execution", "Motor Imagery", "EX minus MI", 
                            "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                            "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials (calculated as different minus same movement context)
RT.Contrasts$`Execution`           <- ((RT_Probe.posteriors$`yes_no_go` - RT_Probe.posteriors$`no_no_go`) + ( RT_Probe.posteriors$`no_yes_go` - RT_Probe.posteriors$`yes_yes_go`) )/2
# Motor Imagery trials (calculated as different minus same movement context)
RT.Contrasts$`Motor Imagery`       <- ((RT_Probe.posteriors$`yes_no_stop` - RT_Probe.posteriors$`no_no_stop`) + ( RT_Probe.posteriors$`no_yes_stop` - RT_Probe.posteriors$`yes_yes_stop`) )/2
# Execution vs Motor Imagery
RT.Contrasts$`EX minus MI`         <- (RT.Contrasts$`Execution` - RT.Contrasts$`Motor Imagery`)
# Execution trials Probe without Obstacle
RT.Contrasts$`Execution NoObs`     <- (RT_Probe.posteriors$`yes_no_go` - RT_Probe.posteriors$`no_no_go`)
# Execution trials with Obstacle
RT.Contrasts$`Execution Obs`       <- (RT_Probe.posteriors$`yes_yes_go` - RT_Probe.posteriors$`no_yes_go`)
# Motor Imagery trials Probe without Obstacle
RT.Contrasts$`Motor Imagery NoObs` <- (RT_Probe.posteriors$`yes_no_stop` - RT_Probe.posteriors$`no_no_stop`)
# Motor Imagery trials Probe with Obstacle
RT.Contrasts$`Motor Imagery Obs`   <- (RT_Probe.posteriors$`yes_yes_stop` - RT_Probe.posteriors$`no_yes_stop`)
# Execution vs Motor Imagery without Obstacle
RT.Contrasts$`EX minus MI NoObs`   <- RT.Contrasts$`Execution NoObs` - RT.Contrasts$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
RT.Contrasts$`EX minus MI Obs`     <- RT.Contrasts$`Execution Obs` - RT.Contrasts$`Motor Imagery Obs`


RT.Contrasts_long <- pivot_longer(RT.Contrasts, cols = everything(),
                                  names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "Motor Imagery Obs",
  "Execution NoObs", "Motor Imagery NoObs", 
  "EX minus MI Obs", "EX minus MI NoObs", 
  "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
RT.Contrasts_long$contrast <- factor(RT.Contrasts_long$contrast, levels = contrast_order)

# Check
head(RT.Contrasts_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
RT.contrasts.summary <-
  RT.Contrasts_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

RT.contrasts.summary$pd <- format(RT.contrasts.summary$pd, nsmall = 1)  
#print(RT.contrasts.summary, n = Inf, width = Inf)


#### Calculate ROPE
( RR       <- rope_range(A069_fit.RT_Probe_OPRIxOPROxSS) )


# Calculate percent in ROPE for contrasts
options(digits=3)
RT.contrast_in_ROPE                       <- as.data.frame(RT.contrasts.summary)
RT.contrast_in_ROPE$lowerROPE             <- RR[1]
RT.contrast_in_ROPE$upperROPE             <- RR[2]
RT.contrast_in_ROPE$CI_range              <- RT.contrast_in_ROPE$upper - RT.contrast_in_ROPE$lower
RT.contrast_in_ROPE$minUpper              <- RT.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
RT.contrast_in_ROPE$maxLower              <- RT.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
RT.contrast_in_ROPE$DiffminUppermaxLower  <- RT.contrast_in_ROPE$minUpper  - RT.contrast_in_ROPE$maxLower 
RT.contrast_in_ROPE$Zeros                 <- rep(0,nrow(RT.contrast_in_ROPE))
RT.contrast_in_ROPE$Overlap               <- RT.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
RT.contrast_in_ROPE$perc_in_ROPE          <- (RT.contrast_in_ROPE$Overlap*100)/RT.contrast_in_ROPE$CI_range
RT.contrast_in_ROPE[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
RT.subj.contrast <- RT.subj.EMM %>%
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

#print(RT.subj.contrast)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
RT.subj.contrast_pooled <- RT.subj.EMM %>%
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

#print(RT.subj.contrast_pooled)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
RT.subj.diffContrast <- RT.subj.contrast %>%
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

#print(RT.subj.diffContrast)


RT.subj.diffContrast_pooled <- RT.subj.contrast_pooled %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(RT.subj.diffContrast_pooled)














#  ==== REACTION TIME: PLOTTING CONTRASTS  ====
# Probe without obstacle
titleX <- "RT difference (°)"
titleY <- ""
g.RT_NoObs.Contrast.A069 <- RT.Contrasts_long %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(RT.contrasts.summary, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.02, ymax = 0.15, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(RT.subj.contrast, probe == "without obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.005,
                                             nudge.y = -0.01,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-50,50),breaks=c(-50,0,50)) + 
  scale_y_continuous(limits = c(-0.02, 0.15)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.RT_NoObs.Contrast.A069


# Probe with obstacle
titleX <- "RT difference (°)"
titleY <- ""
g.RT_Obs.Contrast.A069 <- RT.Contrasts_long %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(RT.contrasts.summary, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.02, ymax = 0.15, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(RT.subj.contrast, probe == "with obstacle"),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.005,
                                             nudge.y = -0.01,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-100,50),breaks=c(-100,-50,0,50)) + 
  scale_y_continuous(limits = c(-0.02, 0.15)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.RT_Obs.Contrast.A069



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~RT~difference~(ms))
titleY <- ""
g.RT_NoObs_EXvsMI.Contrast.A069 <- RT.Contrasts_long %>% filter(contrast=="EX minus MI NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(RT.contrasts.summary, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.02, ymax = 0.15, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(RT.subj.diffContrast, probe == "without obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.005,
                                             nudge.y = -0.01,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-50,50),breaks=c(-50,0,50)) + 
  scale_y_continuous(limits = c(-0.02, 0.15)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.RT_NoObs_EXvsMI.Contrast.A069




# Probe with obstacle
titleX <- expression(Delta~RT~difference~(ms))
titleY <- ""
g.RT_Obs_EXvsMI.Contrast.A069 <- RT.Contrasts_long %>% filter(contrast=="EX minus MI Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(RT.contrasts.summary, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.02, ymax = 0.15, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(RT.subj.diffContrast, probe == "with obstacle"),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.005,
                                             nudge.y = -0.01,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-100,50),breaks=c(-100,-50,0,50)) + 
  scale_y_continuous(limits = c(-0.02, 0.15)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.RT_Obs_EXvsMI.Contrast.A069



#merge into a 1x2 plot grid (library patchwork)
g.RT.Contrast.A069 <-
  (g.RT_Obs.Contrast.A069 + g.RT_NoObs.Contrast.A069 + g.RT_Obs_EXvsMI.Contrast.A069 + g.RT_NoObs_EXvsMI.Contrast.A069) +
  plot_layout(
    widths = c(1, 110/100), heights = c(2,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.RT.Contrast.A069





#  ==== REACTION TIME: PLOTTING CONTRASTS POOLED OVER OBSTACLE PROBE ====
titleX <- "RT difference (°)"
titleY <- ""
g.RT.Contrast2.A069 <- RT.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution" = "Execution",
    "Motor Imagery" = "Motor Imagery"
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
  geom_pointinterval(data = dplyr::filter(RT.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.02, ymax = 0.2, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(RT.subj.contrast_pooled),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.005,
                                             nudge.y = -0.01,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-50,75),breaks=c(-50,-25,0,25,50,75)) + 
  scale_y_continuous(limits = c(-0.02, 0.2)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.RT.Contrast2.A069



### Execution vs Motor Imagery ###
titleX <- expression(Delta~RT~difference~(ms))
titleY <- ""
g.RT_EXvsMI.Contrast2.A069 <- RT.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI" = "EX minus MI"
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
  geom_pointinterval(data = dplyr::filter(RT.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.02, ymax = 0.15, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(RT.subj.diffContrast_pooled),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.005,
                                             nudge.y = -0.01,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-50,75),breaks=c(-50,-25,0,25,50,75)) + 
  scale_y_continuous(limits = c(-0.02, 0.15)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.RT_EXvsMI.Contrast2.A069



#merge into a 1x2 plot grid (library patchwork)
g.RT.Contrast_Pooled.A069 <-
  (g.RT.Contrast2.A069 + g.RT_EXvsMI.Contrast2.A069) +
  plot_layout(
    widths = c(1), heights = c(2,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.RT.Contrast_Pooled.A069




#==== FINAL REACH ERROR PROBE ====
# measured as angular difference btw cursor mov_onset target center and cursor mov_onset cursor at target hit (EndPointError_abs_Probe)
# Since obstacle probe & no obstacle probe trials differ wrt to their underlying distribution, we model them separately
# log_normal for no obstacle probe, gaussian for obstacle probe


# ==== MODEL FITTING: FINAL REACH ERROR NO OBSTACLE PROBE: 2-FACTOR INTERACTION MODEL WITH TRUNCATED (LB=0, UB=4) GAUSSIAN DISTRIBUTION ====
dClean_noobs <- dClean %>% filter(obstacle_probe=="no")

dClean_noobs %>%
  ggplot( aes(x = EndPointError_abs_Probe, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.3, boundary = 0, position = 'identity' ) +
  labs( title = "A069- Histogram", subtitle = "Absolute End Point Error Probe" ) +
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
A069_fit.FRE_Probe_NoObs_OPRIxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.FRE_Probe_NoObs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )

#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.FRE_Probe_NoObs_OPRIxSS,   type = "dens_overlay", ndraws = 100)  # 
pp_check(A069_fit.FRE_Probe_NoObs_OPRIxSS,   type = "hist", ndraws = 10)           # 
pp_check(A069_fit.FRE_Probe_NoObs_OPRIxSS,   type = "boxplot", ndraws = 10)  # 


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.FRE_Probe_NoObs_OPRIxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.FRE_Probe_NoObs_OPRIxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_noobs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_noobs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_NoObs.posteriors <- as.data.frame(fitted(
  A069_fit.FRE_Probe_NoObs_OPRIxSS,
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
  A069_fit.FRE_Probe_NoObs_OPRIxSS,
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
  labs( title = "A069- Histogram", subtitle = "Absolute End Point Error Probe" ) +
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
prior_FRE_obs <-      c(prior(normal(  2.1,     0.5   ),  class = "Intercept"), 
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
A069_fit.FRE_Probe_Obs_OPRIxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.FRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )







#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.FRE_Probe_Obs_OPRIxSS,   type = "dens_overlay", ndraws = 100)  # looks ok
pp_check(A069_fit.FRE_Probe_Obs_OPRIxSS,   type = "hist", ndraws = 10)           # looks ok
pp_check(A069_fit.FRE_Probe_Obs_OPRIxSS,   type = "boxplot", ndraws = 10)        # looks ok


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.FRE_Probe_Obs_OPRIxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.FRE_Probe_Obs_OPRIxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean_obs$obstacle_prime),
                        stop_signal_prime  = levels(dClean_obs$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
FRE_Probe_Obs.posteriors <- as.data.frame(fitted(
  A069_fit.FRE_Probe_Obs_OPRIxSS,
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
  A069_fit.FRE_Probe_Obs_OPRIxSS,
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






# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S11AB) ====
# Prepare and tidy data
FRE_NoObs.EMM$obstacle_probe   <- as.factor("no") 
FRE_Obs.EMM$obstacle_probe     <- as.factor("yes") 
FRE.EMM <- rbind(FRE_NoObs.EMM,FRE_Obs.EMM)

FRE_NoObs.subj.EMM$obstacle_probe   <- as.factor("no") 
FRE_Obs.subj.EMM$obstacle_probe     <- as.factor("yes") 
FRE.subj.EMM <- rbind(FRE_NoObs.subj.EMM,FRE_Obs.subj.EMM)

# rename factors
FRE.EMM <- FRE.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

FRE.EMM

# rename factors
FRE.subj.EMM <- FRE.subj.EMM  %>% 
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

g.FRE_NoObs.EMM.A069 <-
  FRE.subj.EMM %>%
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
    data = dplyr::filter(FRE.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM, probe == "without obstacle"),
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

g.FRE_NoObs.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Final Reach Error (°)"

g.FRE_Obs.EMM.A069 <-
  FRE.subj.EMM %>%
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
    data = dplyr::filter(FRE.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(FRE.EMM, probe == "with obstacle"),
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

g.FRE_Obs.EMM.A069









# ==== FINAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
# IMPORTANT!!!!: FOR POOLED CONSTRASTS: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
FRE.Contrasts           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(FRE.Contrasts) <- c("Execution", "Motor Imagery", "EX minus MI", 
                             "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                             "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials
FRE.Contrasts$`Execution`           <- (FRE_Probe_NoObs.posteriors$`yes_go` + FRE_Probe_Obs.posteriors$`no_go`)/2 - ( FRE_Probe_NoObs.posteriors$`no_go` + FRE_Probe_Obs.posteriors$`yes_go`)/2
# Motor Imagery trials
FRE.Contrasts$`Motor Imagery`       <- (FRE_Probe_NoObs.posteriors$`yes_stop` + FRE_Probe_Obs.posteriors$`no_stop`)/2 - ( FRE_Probe_NoObs.posteriors$`no_stop` + FRE_Probe_Obs.posteriors$`yes_stop`)/2
# Execution vs Motor Imagery
FRE.Contrasts$`EX minus MI`         <- (FRE.Contrasts$`Execution` - FRE.Contrasts$`Motor Imagery`)
# Execution trials Probe without Obstacle
FRE.Contrasts$`Execution NoObs`     <- (FRE_Probe_NoObs.posteriors$`yes_go` - FRE_Probe_NoObs.posteriors$`no_go`)
# Execution trials with Obstacle
FRE.Contrasts$`Execution Obs`       <- (FRE_Probe_Obs.posteriors$`yes_go` - FRE_Probe_Obs.posteriors$`no_go`)
# Motor Imagery trials Probe without Obstacle
FRE.Contrasts$`Motor Imagery NoObs` <- (FRE_Probe_NoObs.posteriors$`yes_stop` - FRE_Probe_NoObs.posteriors$`no_stop`)
# Motor Imagery trials Probe with Obstacle
FRE.Contrasts$`Motor Imagery Obs`   <- (FRE_Probe_Obs.posteriors$`yes_stop` - FRE_Probe_Obs.posteriors$`no_stop`)
# Execution vs Motor Imagery without Obstacle
FRE.Contrasts$`EX minus MI NoObs`   <- FRE.Contrasts$`Execution NoObs` - FRE.Contrasts$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
FRE.Contrasts$`EX minus MI Obs`     <- FRE.Contrasts$`Execution Obs` - FRE.Contrasts$`Motor Imagery Obs`


FRE.Contrasts_long <- pivot_longer(FRE.Contrasts, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution Obs", "Motor Imagery Obs",
  "Execution NoObs", "Motor Imagery NoObs", 
  "EX minus MI Obs", "EX minus MI NoObs", 
  "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
FRE.Contrasts_long$contrast <- factor(FRE.Contrasts_long$contrast, levels = contrast_order)

# Check
head(FRE.Contrasts_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
FRE.contrasts.summary <-
  FRE.Contrasts_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

FRE.contrasts.summary$pd <- format(FRE.contrasts.summary$pd, nsmall = 4)  
#print(FRE.contrasts.summary, n = Inf, width = Inf)


#### Calculate ROPE
( RR_NoObs <- rope_range(A069_fit.FRE_Probe_NoObs_OPRIxSS) )
( RR_Obs   <- rope_range(A069_fit.FRE_Probe_Obs_OPRIxSS) )
( RR       <- ( RR_NoObs + RR_Obs)/2 )

# Calculate percent in ROPE for contrasts
options(digits=3)
FRE.contrast_in_ROPE                       <- as.data.frame(FRE.contrasts.summary)
FRE.contrast_in_ROPE$lowerROPE             <- NA
FRE.contrast_in_ROPE$lowerROPE[c(1, 2, 5)] <- RR_Obs[1]
FRE.contrast_in_ROPE$lowerROPE[c(3, 4, 6)] <- RR_NoObs[1]
FRE.contrast_in_ROPE$lowerROPE[c(7:9)]     <- RR[1]
FRE.contrast_in_ROPE$upperROPE             <- NA
FRE.contrast_in_ROPE$upperROPE[c(1, 2, 5)] <- RR_Obs[2]
FRE.contrast_in_ROPE$upperROPE[c(3, 4, 6)] <- RR_NoObs[2]
FRE.contrast_in_ROPE$upperROPE[c(7:9)]     <- RR[2]
FRE.contrast_in_ROPE$CI_range              <- FRE.contrast_in_ROPE$upper - FRE.contrast_in_ROPE$lower
FRE.contrast_in_ROPE$minUpper              <- FRE.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
FRE.contrast_in_ROPE$maxLower              <-  FRE.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
FRE.contrast_in_ROPE$DiffminUppermaxLower  <- FRE.contrast_in_ROPE$minUpper  - FRE.contrast_in_ROPE$maxLower 
FRE.contrast_in_ROPE$Zeros                 <- rep(0,nrow(FRE.contrast_in_ROPE))
FRE.contrast_in_ROPE$Overlap               <- FRE.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
FRE.contrast_in_ROPE$perc_in_ROPE          <- (FRE.contrast_in_ROPE$Overlap*100)/FRE.contrast_in_ROPE$CI_range
FRE.contrast_in_ROPE[,c(1:7,14)]


# subject-wise contrasts: Prime with obstacle minus Prime without obstacle 
FRE.subj.contrast <- FRE.subj.EMM %>%
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

#print(FRE.subj.contrast)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
FRE.subj.contrast_pooled <- FRE.subj.EMM %>%
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

#print(FRE.subj.contrast_pooled)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
FRE.subj.diffContrast <- FRE.subj.contrast %>%
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

#print(FRE.subj.diffContrast)


FRE.subj.diffContrast_pooled <- FRE.subj.contrast_pooled %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(FRE.subj.diffContrast_pooled)









#  ==== FINAL REACH ERROR: PLOTTING CONTRASTS (SUPPLEMENTARY FIGURE S11CDEF) ====
# Probe without obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_NoObs.Contrast.A069 <- FRE.Contrasts_long %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs[1], xmax = RR_NoObs[2], ymin = -2, ymax = 15, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast, probe == "without obstacle"),
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
  scale_y_continuous(limits = c(-2, 15)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs.Contrast.A069


# Probe with obstacle
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE_Obs.Contrast.A069 <- FRE.Contrasts_long %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs[1], xmax = RR_Obs[2], ymin = -2, ymax = 16, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast, probe == "with obstacle"),
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
  scale_y_continuous(limits = c(-2, 16)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs.Contrast.A069



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_NoObs_EXvsMI.Contrast.A069 <- FRE.Contrasts_long %>% filter(contrast=="EX minus MI NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_NoObs[1], xmax = RR_NoObs[2], ymin = -2, ymax = 10, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast, probe == "without obstacle"),
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
  scale_y_continuous(limits = c(-2, 10)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_NoObs_EXvsMI.Contrast.A069




# Probe with obstacle
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_Obs_EXvsMI.Contrast.A069 <- FRE.Contrasts_long %>% filter(contrast=="EX minus MI Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR_Obs[1], xmax = RR_Obs[2], ymin = -2, ymax = 11, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast, probe == "with obstacle"),
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
  scale_y_continuous(limits = c(-2, 11)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_Obs_EXvsMI.Contrast.A069



#  ==== FINAL REACH ERROR: PLOTTING CONTRASTS POOLED OVER OBSTACLE PROBE ====
titleX <- "Final Reach Bias (°)"
titleY <- ""
g.FRE.Contrast2.A069 <- FRE.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution" = "Execution",
    "Motor Imagery" = "Motor Imagery"
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
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 22, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.contrast_pooled),
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
  scale_y_continuous(limits = c(-2, 22)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE.Contrast2.A069



### Execution vs Motor Imagery ###
titleX <- expression(Delta~Final~Reach~Bias~(degree))
titleY <- ""
g.FRE_EXvsMI.Contrast2.A069 <- FRE.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI" = "EX minus MI"
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
  geom_pointinterval(data = dplyr::filter(FRE.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 22, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(FRE.subj.diffContrast_pooled),
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
  scale_y_continuous(limits = c(-2, 22)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.FRE_EXvsMI.Contrast2.A069



#merge into a 1x2 plot grid (library patchwork)
g.FRE.Contrast_Pooled.A069 <-
  (g.FRE.Contrast2.A069 + g.FRE_EXvsMI.Contrast2.A069) +
  plot_layout(
    widths = c(1), heights = c(2,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.FRE.Contrast_Pooled.A069


#===== MAKE ONE Figure EMM + Constrasts Final Reach Error (SUPPLEMENTARY FIGURE S11) ====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 2, r = 1), # p1
  area(t = 1, l = 2, b = 2, r = 2), # p2
  area(t = 3, l = 1, b = 4, r = 1),  # p3
  area(t = 3, l = 2, b = 4, r = 2),   #p4
  area(t = 5, l = 1, b = 5, r = 1),   #p5
  area(t = 5, l = 2, b = 5, r = 2)  #p6
  
)


g.FRE.A069 <-  g.FRE_Obs.EMM.A069 + g.FRE_NoObs.EMM.A069 + g.FRE_Obs.Contrast.A069 + g.FRE_NoObs.Contrast.A069 + 
  g.FRE_Obs_EXvsMI.Contrast.A069 + g.FRE_NoObs_EXvsMI.Contrast.A069 + 
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A")
g.FRE.A069



# Build a valid path
outfile <- file.path(figurePath, "FigS11_A069_FinalReachError_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.FRE.A069,
  width = 24, height = 28, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)



# Build a valid path
outfile <- file.path(figurePath, "FigS11_A069_FinalReachError_EMM_Contrasts_ALL.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.FRE.A069,
  width = 24, height = 28, units = "cm",
  device = "svg"
)

outfile <- file.path(figurePath, "FigS11_A069_FinalReachError_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.FRE.A069,
  width = 24, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)




#==== ACCURACY: END POINT ERROR 150ms PROBE====
# measured as absolute distance of cursor 150ms after target hit to target center (EndPointError_150_Probe)
# ==== MODEL FITTING: ACCURACY: 3-FACTOR INTERACTION MODEL WITH TRUNCATED GAUSSIAN DISTRIBUTION ====
dClean_EPE150 <- dClean %>% filter(EndPointError_150_Probe <=0.5 & Vel_HandTargetReached_150_Probe<=5 )

dClean_EPE150 %>%
  ggplot( aes(x = EndPointError_150_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 0.02, boundary = 0, position = 'identity' ) +
  labs( title = "A069- Histogram", subtitle = "End Point Error 150ms after Target Hit Probe" ) +
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
A069_fit.EPE150_Probe_OPRIxOPROxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.EPE150_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )





















#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.EPE150_Probe_OPRIxOPROxSS,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.EPE150_Probe_OPRIxOPROxSS,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.EPE150_Probe_OPRIxOPROxSS,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.EPE150_Probe_OPRIxOPROxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.EPE150_Probe_OPRIxOPROxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean$obstacle_prime),
                        obstacle_probe     = levels(dClean$obstacle_probe),
                        stop_signal_prime  = levels(dClean$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
EPE150_Probe.posteriors <- as.data.frame(fitted(
  A069_fit.EPE150_Probe_OPRIxOPROxSS,
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
  A069_fit.EPE150_Probe_OPRIxOPROxSS,
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





# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY S13AB) ====
# Prepare and tidy data
# rename factors
EPE150.EMM <- EPE150.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle")),
    .value      = .value*10, .lower = .lower*10 , .upper = .upper*10  # convert from cm to mm
    # convert from cm to mm
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns
EPE150.EMM

# rename factors
EPE150.subj.EMM <- EPE150.subj.EMM  %>% 
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

g.EPE150_NoObs.EMM.A069 <-
  EPE150.subj.EMM %>%
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
    data = dplyr::filter(EPE150.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPE150.EMM, probe == "without obstacle"),
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

g.EPE150_NoObs.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Endpoint Error (mm)"

g.EPE150_Obs.EMM.A069 <-
  EPE150.subj.EMM %>%
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
    data = dplyr::filter(EPE150.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPE150.EMM, probe == "with obstacle"),
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

g.EPE150_Obs.EMM.A069







# ==== ACCURACY: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
EPE150.Contrasts           <- data.frame(matrix(ncol = 10, nrow = 40000))
colnames(EPE150.Contrasts) <- c("Execution", "Motor Imagery", "EX minus MI", 
                                "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                                "EX minus MI NoObs", "EX minus MI Obs",
                                "EX vs MI")

# FOR POOLED DATA: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
# Execution trials (calculated as different minus same movement context)
EPE150.Contrasts$`Execution`           <- ((EPE150_Probe.posteriors$`yes_no_go` - EPE150_Probe.posteriors$`no_no_go`) + ( EPE150_Probe.posteriors$`no_yes_go` - EPE150_Probe.posteriors$`yes_yes_go`) )/2
# Motor Imagery trials (calculated as different minus same movement context)
EPE150.Contrasts$`Motor Imagery`       <- ((EPE150_Probe.posteriors$`yes_no_stop` - EPE150_Probe.posteriors$`no_no_stop`) + ( EPE150_Probe.posteriors$`no_yes_stop` - EPE150_Probe.posteriors$`yes_yes_stop`) )/2
# Difference of Execution vs Motor Imagery
EPE150.Contrasts$`EX minus MI`         <- (EPE150.Contrasts$`Execution` - EPE150.Contrasts$`Motor Imagery`)
# Execution trials Probe without Obstacle
EPE150.Contrasts$`Execution NoObs`     <- (EPE150_Probe.posteriors$`yes_no_go` - EPE150_Probe.posteriors$`no_no_go`)
# Execution trials with Obstacle
EPE150.Contrasts$`Execution Obs`       <- (EPE150_Probe.posteriors$`yes_yes_go` - EPE150_Probe.posteriors$`no_yes_go`)
# Motor Imagery trials Probe without Obstacle
EPE150.Contrasts$`Motor Imagery NoObs` <- (EPE150_Probe.posteriors$`yes_no_stop` - EPE150_Probe.posteriors$`no_no_stop`)
# Motor Imagery trials Probe with Obstacle
EPE150.Contrasts$`Motor Imagery Obs`   <- (EPE150_Probe.posteriors$`yes_yes_stop` - EPE150_Probe.posteriors$`no_yes_stop`)
# Execution vs Motor Imagery without Obstacle
EPE150.Contrasts$`EX minus MI NoObs`   <- EPE150.Contrasts$`Execution NoObs` - EPE150.Contrasts$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
EPE150.Contrasts$`EX minus MI Obs`     <- EPE150.Contrasts$`Execution Obs` - EPE150.Contrasts$`Motor Imagery Obs`
# Execution vs Motor Imagery (Overall difference)
EPE150.Contrasts$`EX vs MI`            <- ((EPE150_Probe.posteriors$`yes_yes_go` + EPE150_Probe.posteriors$`yes_no_go` + EPE150_Probe.posteriors$`no_yes_go` + EPE150_Probe.posteriors$`no_no_go` ) /4 ) - ((EPE150_Probe.posteriors$`yes_yes_stop` + EPE150_Probe.posteriors$`yes_no_stop` + EPE150_Probe.posteriors$`no_yes_stop` + EPE150_Probe.posteriors$`no_no_stop` ) /4 )

EPE150.Contrasts_long <- pivot_longer(EPE150.Contrasts*10, cols = everything(), # convert from cm to mm
                                      names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c("EX vs MI",
                    "Execution Obs", "Motor Imagery Obs",
                    "Execution NoObs", "Motor Imagery NoObs", 
                    "EX minus MI Obs", "EX minus MI NoObs", 
                    "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
EPE150.Contrasts_long$contrast <- factor(EPE150.Contrasts_long$contrast, levels = contrast_order)

# Check
head(EPE150.Contrasts_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
EPE150.contrasts.summary <-
  EPE150.Contrasts_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

EPE150.contrasts.summary$pd <- format(EPE150.contrasts.summary$pd, nsmall = 1)  
print(EPE150.contrasts.summary, n = Inf, width = Inf)


#### Calculate ROPE
( RR       <- rope_range(A069_fit.EPE150_Probe_OPRIxOPROxSS)*10 ) # convert from cm to mm


# Calculate percent in ROPE for contrasts
options(digits=3)
EPE150.contrast_in_ROPE                       <- as.data.frame(EPE150.contrasts.summary)
EPE150.contrast_in_ROPE$lowerROPE             <- RR[1]
EPE150.contrast_in_ROPE$upperROPE             <- RR[2]
EPE150.contrast_in_ROPE$CI_range              <- EPE150.contrast_in_ROPE$upper - EPE150.contrast_in_ROPE$lower
EPE150.contrast_in_ROPE$minUpper              <- EPE150.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
EPE150.contrast_in_ROPE$maxLower              <- EPE150.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
EPE150.contrast_in_ROPE$DiffminUppermaxLower  <- EPE150.contrast_in_ROPE$minUpper  - EPE150.contrast_in_ROPE$maxLower 
EPE150.contrast_in_ROPE$Zeros                 <- rep(0,nrow(EPE150.contrast_in_ROPE))
EPE150.contrast_in_ROPE$Overlap               <- EPE150.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
EPE150.contrast_in_ROPE$perc_in_ROPE          <- (EPE150.contrast_in_ROPE$Overlap*100)/EPE150.contrast_in_ROPE$CI_range
EPE150.contrast_in_ROPE[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
EPE150.ExvsMI.subj.contrast <- EPE150.subj.EMM %>%
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
EPE150.subj.contrast <- EPE150.subj.EMM %>%
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

#print(EPE150.subj.contrast)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
EPE150.subj.contrast_pooled <- EPE150.subj.EMM %>%
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

#print(EPE150.subj.contrast_pooled)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
EPE150.subj.diffContrast <- EPE150.subj.contrast %>%
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

#print(EPE150.subj.diffContrast)


EPE150.subj.diffContrast_pooled <- EPE150.subj.contrast_pooled %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(EPE150.subj.diffContrast_pooled)






#  ==== ACCURACY: PLOTTING CONTRASTS OVERALL EXECUTION VS MOTOR IMAGERY (SUPPLEMENTARY S13C)  ====
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- "EX minus MI"
g.EPE150_EXvsMI.Contrast.A069 <- EPE150.Contrasts_long %>% filter(contrast=="EX vs MI") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "EX vs MI"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.ExvsMI.subj.contrast, contrast == "EX vs MI"),
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

g.EPE150_EXvsMI.Contrast.A069



#  ==== ACCURACY:: PLOTTING CONTRASTS (SUPPLEMENTARY S13DEFG)  ====
# Probe without obstacle
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150_NoObs.Contrast.A069 <- EPE150.Contrasts_long %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast, probe == "without obstacle"),
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

g.EPE150_NoObs.Contrast.A069


# Probe with obstacle
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150_Obs.Contrast.A069 <- EPE150.Contrasts_long %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 12, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast, probe == "with obstacle"),
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

g.EPE150_Obs.Contrast.A069



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~EPE~Difference(mm))
titleY <- ""
g.EPE150_NoObs_EXvsMI.Contrast.A069 <- EPE150.Contrasts_long %>% filter(contrast=="EX minus MI NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 10, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast, probe == "without obstacle"),
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
  scale_y_continuous(limits = c(-2, 10)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_NoObs_EXvsMI.Contrast.A069




# Probe with obstacle
titleX <- expression(Delta~EPE~Difference(mm))
titleY <- ""
g.EPE150_Obs_EXvsMI.Contrast.A069 <- EPE150.Contrasts_long %>% filter(contrast=="EX minus MI Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 10, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast, probe == "with obstacle"),
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
  scale_y_continuous(limits = c(-2, 10)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_Obs_EXvsMI.Contrast.A069





#  ==== ACCURACY: PLOTTING CONTRASTS POOLED OVER OBSTACLE PROBE ====
titleX <- expression(Delta~Endpoint~Error~(mm))
titleY <- ""
g.EPE150.Contrast2.A069 <- EPE150.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution" = "Execution",
    "Motor Imagery" = "Motor Imagery"
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.contrast_pooled),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.2,0.2),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150.Contrast2.A069



### Execution vs Motor Imagery ###
titleX <- expression(Delta~Endpoint~Error~Difference(mm))
titleY <- ""
g.EPE150_EXvsMI.Contrast2.A069 <- EPE150.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI" = "EX minus MI"
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
  geom_pointinterval(data = dplyr::filter(EPE150.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -2, ymax = 20, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPE150.subj.diffContrast_pooled),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.5,
                                             nudge.y = -1,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-0.2,0.2),breaks=c(-0.2,0,0.2)) + 
  scale_y_continuous(limits = c(-2, 20)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPE150_EXvsMI.Contrast2.A069



#merge into a 1x2 plot grid (library patchwork)
g.EPE150.Contrast_Pooled.A069 <-
  (g.EPE150.Contrast2.A069 + g.EPE150_EXvsMI.Contrast2.A069) +
  plot_layout(
    widths = c(1), heights = c(2,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.EPE150.Contrast_Pooled.A069




#  ==== ACCURACY: COMBINE PLOTS (SUPPLEMENTARY FIGURE S13) ====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 2, r = 3), # p1
  area(t = 1, l = 4, b = 2, r = 6), # p2
  area(t = 1, l = 7, b = 2, r = 8), # p3
  
  area(t = 3, l = 1, b = 4, r = 4),  # p4
  area(t = 3, l = 5, b = 4, r = 8),   #p5
  area(t = 5, l = 1, b = 5, r = 4),   #p6
  area(t = 5, l = 5, b = 5, r = 8)  #p7
  
)


g.EPE150.A069 <-  g.EPE150_Obs.EMM.A069 + g.EPE150_NoObs.EMM.A069 + g.EPE150_EXvsMI.Contrast.A069 + 
  g.EPE150_Obs.Contrast.A069 + g.EPE150_NoObs.Contrast.A069 + 
  g.EPE150_Obs_EXvsMI.Contrast.A069 + g.EPE150_NoObs_EXvsMI.Contrast.A069 + 
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A")
g.EPE150.A069




# Build a valid path
outfile <- file.path(figurePath, "FigS13_A069_EPE150_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.EPE150.A069,
  width = 24, height = 28, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS13_A069_EPE150_EMM_Contrasts_ALL.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.EPE150.A069,
  width = 24, height = 28, units = "cm",
  device = "svg"
)



outfile <- file.path(figurePath, "FigS13_A069_EPE150_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.EPE150.A069,
  width = 24, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)


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


# Run repeated-measures ANOVA
anova_ellipse <- aov_ez(
  id = "subID",
  dv = "ellipse_area",
  within = c("trial_type", "prime", "probe"),
  data = data_ellipse_area,
  anova_table = list(correction = "GG", es = "pes")  # Greenhouse-Geisser correction
)

knitr::kable(nice(anova_ellipse))

#marginal means for stop-signal
anova_ellipse.SS <- emmeans(anova_ellipse, ~trial_type)
anova_ellipse.SS

#marginal means for Interaction obstacle_prime x obstacle_probe x stop-signal
anova_ellipse.OPRIxOPROxSS <- emmeans(anova_ellipse, ~prime:probe:trial_type)
anova_ellipse.OPRIxOPROxSS

# Quick plot Interaction: obstacle_prime x obstacle_probe
ggplot(as.data.frame(anova_ellipse.OPRIxOPROxSS), aes(x=probe, y=emmean, color=prime)) + 
  facet_wrap(~trial_type) +
  geom_point(position = position_dodge(width=0.5), size=3) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position = position_dodge(width=0.5), width=0.2) + 
  ylab(label = "Ellipse Area (cm²)") + 
  ggtitle(label = "End Point Precision 150ms")

#comparison with fdr corrected p-values
update(pairs(anova_ellipse.OPRIxOPROxSS, by = c("probe","trial_type"), adjust="fdr"))







# ==== MODEL FITTING: PRECISION: 3-FACTOR INTERACTION MODEL WITH TRUNCATED GAUSSIAN DISTRIBUTION ====
data_ellipse_area %>%
  ggplot( aes(x = ellipse_area, fill = prime) ) +
  facet_wrap(~probe:trial_type) +
  geom_histogram( alpha = 0.3, binwidth = 5, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "End Point Precision 150ms after Target Hit Probe" ) +
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
# runs about 6 min, introducing truncation in prior predicitve checks causes NA's
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A069_fit.EPP150_Probe_OPRIxOPROxSS <- brm( 
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
  file = paste0(analysisPath,'/', 'A069_fit.EPP150_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )






#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.EPP150_Probe_OPRIxOPROxSS,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.EPP150_Probe_OPRIxOPROxSS,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.EPP150_Probe_OPRIxOPROxSS,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.EPP150_Probe_OPRIxOPROxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.EPP150_Probe_OPRIxOPROxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)



#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(prime       = levels(data_ellipse_area$prime),
                        probe       = levels(data_ellipse_area$probe),
                        trial_type  = levels(data_ellipse_area$trial_type))


# Posterior draws of expected values (population-level, no random effects)
EPP150_Probe.posteriors <- as.data.frame(fitted(
  A069_fit.EPP150_Probe_OPRIxOPROxSS,
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
  A069_fit.EPP150_Probe_OPRIxOPROxSS,
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





# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S12AB) ====
# Prepare and tidy data
# rename factors
EPP150.EMM <- EPP150.EMM  %>% 
  mutate(
    trial_type  = factor(trial_type),
    prime       = factor(prime),
    probe       = factor(probe),
  )   
EPP150.EMM

# rename factors
EPP150.subj.EMM <- EPP150.subj.EMM  %>% 
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

g.EPP150_NoObs.EMM.A069 <-
  EPP150.subj.EMM %>%
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
    data = dplyr::filter(EPP150.EMM, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPP150.EMM, probe == "without obstacle"),
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

g.EPP150_NoObs.EMM.A069



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_mi1,color_exe1,color_mi1)
fontSize    <- 9
titleX      <- ""
titleY      <- "Area Ellipse (mm²)"

g.EPP150_Obs.EMM.A069 <-
  EPP150.subj.EMM %>%
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
    data = dplyr::filter(EPP150.EMM, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(EPP150.EMM, probe == "with obstacle"),
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

g.EPP150_Obs.EMM.A069



# ==== PRECISION: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
EPP150.Contrasts           <- data.frame(matrix(ncol = 10, nrow = 40000))
colnames(EPP150.Contrasts) <- c("Execution", "Motor Imagery", "EX minus MI", 
                                "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                                "EX minus MI NoObs", "EX minus MI Obs",
                                "EX vs MI")

# FOR POOLED DATA: REVERSE DIFFERENCE FOR PROBE TRIALS WITH OBSTACLE SUCH THAT SIMILARITY TO PREVIOUS MOVEMENT MATCHES WITH TRIALS WITHOUT PROBE OBSTACLE
# Execution trials (calculated as different minus same movement context)
EPP150.Contrasts$`Execution`           <- ((EPP150_Probe.posteriors$`with obstacle_without obstacle_Execution` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Execution`) + ( EPP150_Probe.posteriors$`without obstacle_with obstacle_Execution` - EPP150_Probe.posteriors$`with obstacle_with obstacle_Execution`) )/2
# Motor Imagery trials (calculated as different minus same movement context)
EPP150.Contrasts$`Motor Imagery`       <- ((EPP150_Probe.posteriors$`with obstacle_without obstacle_Motor Imagery` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Motor Imagery`) + ( EPP150_Probe.posteriors$`without obstacle_with obstacle_Motor Imagery` - EPP150_Probe.posteriors$`with obstacle_with obstacle_Motor Imagery`) )/2
# Difference of Execution vs Motor Imagery
EPP150.Contrasts$`EX minus MI`         <- (EPP150.Contrasts$`Execution` - EPP150.Contrasts$`Motor Imagery`)
# Execution trials Probe without Obstacle
EPP150.Contrasts$`Execution NoObs`     <- (EPP150_Probe.posteriors$`with obstacle_without obstacle_Execution` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Execution`)
# Execution trials with Obstacle
EPP150.Contrasts$`Execution Obs`       <- (EPP150_Probe.posteriors$`with obstacle_with obstacle_Execution` - EPP150_Probe.posteriors$`without obstacle_with obstacle_Execution`)
# Motor Imagery trials Probe without Obstacle
EPP150.Contrasts$`Motor Imagery NoObs` <- (EPP150_Probe.posteriors$`with obstacle_without obstacle_Motor Imagery` - EPP150_Probe.posteriors$`without obstacle_without obstacle_Motor Imagery`)
# Motor Imagery trials Probe with Obstacle
EPP150.Contrasts$`Motor Imagery Obs`   <- (EPP150_Probe.posteriors$`with obstacle_with obstacle_Motor Imagery` - EPP150_Probe.posteriors$`without obstacle_with obstacle_Motor Imagery`)
# Execution vs Motor Imagery without Obstacle
EPP150.Contrasts$`EX minus MI NoObs`   <- EPP150.Contrasts$`Execution NoObs` - EPP150.Contrasts$`Motor Imagery NoObs` 
# Execution vs Motor Imagery trials with Obstacle
EPP150.Contrasts$`EX minus MI Obs`     <- EPP150.Contrasts$`Execution Obs` - EPP150.Contrasts$`Motor Imagery Obs`
# Execution vs Motor Imagery (Overall difference)
EPP150.Contrasts$`EX vs MI`            <- ((EPP150_Probe.posteriors$`with obstacle_with obstacle_Execution` + EPP150_Probe.posteriors$`with obstacle_without obstacle_Execution` + EPP150_Probe.posteriors$`without obstacle_with obstacle_Execution` + EPP150_Probe.posteriors$`without obstacle_without obstacle_Execution` ) /4 ) - ((EPP150_Probe.posteriors$`with obstacle_with obstacle_Motor Imagery` + EPP150_Probe.posteriors$`with obstacle_without obstacle_Motor Imagery` + EPP150_Probe.posteriors$`without obstacle_with obstacle_Motor Imagery` + EPP150_Probe.posteriors$`without obstacle_without obstacle_Motor Imagery` ) /4 )

EPP150.Contrasts_long <- pivot_longer(EPP150.Contrasts, cols = everything(),
                                      names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c("EX vs MI",
                    "Execution Obs", "Motor Imagery Obs",
                    "Execution NoObs", "Motor Imagery NoObs", 
                    "EX minus MI Obs", "EX minus MI NoObs", 
                    "Execution", "Motor Imagery", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
EPP150.Contrasts_long$contrast <- factor(EPP150.Contrasts_long$contrast, levels = contrast_order)

# Check
head(EPP150.Contrasts_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
EPP150.contrasts.summary <-
  EPP150.Contrasts_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

EPP150.contrasts.summary$pd <- format(EPP150.contrasts.summary$pd, nsmall = 1)  
print(EPP150.contrasts.summary, n = Inf, width = Inf)


#### Calculate ROPE
( RR       <- rope_range(A069_fit.EPP150_Probe_OPRIxOPROxSS) ) 


# Calculate percent in ROPE for contrasts
options(digits=3)
EPP150.contrast_in_ROPE                       <- as.data.frame(EPP150.contrasts.summary)
EPP150.contrast_in_ROPE$lowerROPE             <- RR[1]
EPP150.contrast_in_ROPE$upperROPE             <- RR[2]
EPP150.contrast_in_ROPE$CI_range              <- EPP150.contrast_in_ROPE$upper - EPP150.contrast_in_ROPE$lower
EPP150.contrast_in_ROPE$minUpper              <- EPP150.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
EPP150.contrast_in_ROPE$maxLower              <- EPP150.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
EPP150.contrast_in_ROPE$DiffminUppermaxLower  <- EPP150.contrast_in_ROPE$minUpper  - EPP150.contrast_in_ROPE$maxLower 
EPP150.contrast_in_ROPE$Zeros                 <- rep(0,nrow(EPP150.contrast_in_ROPE))
EPP150.contrast_in_ROPE$Overlap               <- EPP150.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
EPP150.contrast_in_ROPE$perc_in_ROPE          <- (EPP150.contrast_in_ROPE$Overlap*100)/EPP150.contrast_in_ROPE$CI_range
EPP150.contrast_in_ROPE[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
EPP150.ExvsMI.subj.contrast <- EPP150.subj.EMM %>%
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
EPP150.subj.contrast <- EPP150.subj.EMM %>%
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

#print(EPP150.subj.contrast)

# subject-wise contrasts: Prime with obstacle minus Prime without obstacle (pooled across Probe)
EPP150.subj.contrast_pooled <- EPP150.subj.EMM %>%
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

#print(EPP150.subj.contrast_pooled)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
EPP150.subj.diffContrast <- EPP150.subj.contrast %>%
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

#print(EPP150.subj.diffContrast)


EPP150.subj.diffContrast_pooled <- EPP150.subj.contrast_pooled %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(EPP150.subj.diffContrast_pooled)






#  ==== PRECISION: PLOTTING CONTRASTS OVERALL EXECUTION VS MOTOR IMAGERY (SUPPLEMENTARY FIGURE S12C) ====
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- "EX minus MI"
g.EPP150_EXvsMI.Contrast.A069 <- EPP150.Contrasts_long %>% filter(contrast=="EX vs MI") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "EX vs MI"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.35, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.ExvsMI.subj.contrast, contrast == "EX vs MI"),
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
  scale_y_continuous(limits = c(-0.05, 0.35)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_EXvsMI.Contrast.A069



#  ==== PRECISION: PLOTTING CONTRASTS  (SUPPLEMENTARY FIGURE S12DEFG) ====
# Probe without obstacle
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150_NoObs.Contrast.A069 <- EPP150.Contrasts_long %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast, probe == "without obstacle"),
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

g.EPP150_NoObs.Contrast.A069


# Probe with obstacle
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150_Obs.Contrast.A069 <- EPP150.Contrasts_long %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast, probe == "with obstacle"),
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

g.EPP150_Obs.Contrast.A069



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~EA~Difference~(mm^2))
titleY <- ""
g.EPP150_NoObs_EXvsMI.Contrast.A069 <- EPP150.Contrasts_long %>% filter(contrast=="EX minus MI NoObs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.2, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast, probe == "without obstacle"),
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
  scale_y_continuous(limits = c(-0.05, 0.2)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_NoObs_EXvsMI.Contrast.A069




# Probe with obstacle
titleX <- expression(Delta~EA~Difference~(mm^2))
titleY <- ""
g.EPP150_Obs_EXvsMI.Contrast.A069 <- EPP150.Contrasts_long %>% filter(contrast=="EX minus MI Obs") %>% 
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.2, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast, probe == "with obstacle"),
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
  scale_y_continuous(limits = c(-0.05, 0.2)) +
  theme_cowplot() + custom_plot_theme +
  # theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_Obs_EXvsMI.Contrast.A069







#  ==== PRECISION: PLOTTING CONTRASTS POOLED OVER OBSTACLE PROBE ====
titleX <- expression(Delta~Ellipse~Area~(mm^2))
titleY <- ""
g.EPP150.Contrast2.A069 <- EPP150.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution" = "Execution",
    "Motor Imagery" = "Motor Imagery"
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.5, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.contrast_pooled),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-5,10),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.5)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150.Contrast2.A069



### Execution vs Motor Imagery ###
titleX <- expression(Delta~Ellipse~Area~Difference~(mm^2))
titleY <- ""
g.EPP150_EXvsMI.Contrast2.A069 <- EPP150.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI" = "EX minus MI"
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
  geom_pointinterval(data = dplyr::filter(EPP150.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.05, ymax = 0.5, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(EPP150.subj.diffContrast_pooled),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.01,
                                             nudge.y = -0.025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-5,10),breaks=c(-15,-10,-5,0,5,10,15)) + 
  scale_y_continuous(limits = c(-0.05, 0.5)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.EPP150_EXvsMI.Contrast2.A069



#merge into a 1x2 plot grid (library patchwork)
g.EPP150.Contrast_Pooled.A069 <-
  (g.EPP150.Contrast2.A069 + g.EPP150_EXvsMI.Contrast2.A069) +
  plot_layout(
    widths = c(1), heights = c(2,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.EPP150.Contrast_Pooled.A069



#  ==== PRECISION: COMBINE PLOTS (SUPPLEMENTARY FIGURE S12) ====
### ALL in ONE
# merge into a 5x4 plot grid (library patchwork)
layout <- c(
  area(t = 1, l = 1, b = 2, r = 3), # p1
  area(t = 1, l = 4, b = 2, r = 6), # p2
  area(t = 1, l = 7, b = 2, r = 8), # p3
  
  area(t = 3, l = 1, b = 4, r = 4),  # p4
  area(t = 3, l = 5, b = 4, r = 8),   #p5
  area(t = 5, l = 1, b = 5, r = 4),   #p6
  area(t = 5, l = 5, b = 5, r = 8)  #p7
  
)


g.EPP150.A069 <-  g.EPP150_Obs.EMM.A069 + g.EPP150_NoObs.EMM.A069 + g.EPP150_EXvsMI.Contrast.A069 + 
  g.EPP150_Obs.Contrast.A069 + g.EPP150_NoObs.Contrast.A069 + 
  g.EPP150_Obs_EXvsMI.Contrast.A069 + g.EPP150_NoObs_EXvsMI.Contrast.A069 + 
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A")
g.EPP150.A069



# Build a valid path
outfile <- file.path(figurePath, "FigS12_A069_EPP150_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.EPP150.A069,
  width = 24, height = 28, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)


# Build a valid path
outfile <- file.path(figurePath, "FigS12_A069_EPP150_EMM_Contrasts_ALL.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.EPP150.A069,
  width = 24, height = 28, units = "cm",
  device = "svg"
)

outfile <- file.path(figurePath, "FigS12_A069_EPP150_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.EPP150.A069,
  width = 24, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)




#==== CORRELATION ANALYSES ====
# ==== CORRELATION ANALYSES: DATA PREPARATION ====
d.corr <- dClean %>%
  group_by(subID, obstacle_prime, obstacle_probe, stop_signal_prime) %>%
  summarise(
    medianReachDiff   = median(ReachDiff2_Probe, na.rm = TRUE),
    MIQ_Vis           = mean(MIQ_Score_Vis,     na.rm = TRUE),
    MIQ_Kin           = mean(MIQ_Score_Kin,     na.rm = TRUE),
    MIQ_Overall       = mean(MIQ_Score_Overall, na.rm = TRUE),
    MI_Ease_Overall   = mean(MI_Ease_Overall,   na.rm = TRUE),
    MI_Count_Overall  = mean(MI_Count_Overall,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = obstacle_prime,   # "yes"/"no"
    values_from = medianReachDiff
  ) %>%
  mutate(
    ReachBias         = yes - no,
    probe    = recode(obstacle_probe, "no" = "without obstacle", "yes" = "with obstacle"),
    trial_type = recode(stop_signal_prime, "go" = "Execution", "stop" = "Motor Imagery")
  )  %>% 
  select(-c(stop_signal_prime,obstacle_probe)) # remove old columns


d.corr.noobs_go   <- d.corr  %>%  filter(probe=="without obstacle" & trial_type=="Execution")
d.corr.noobs_stop <- d.corr  %>%  filter(probe=="without obstacle" & trial_type=="Motor Imagery")
d.corr.obs_go     <- d.corr  %>%  filter(probe=="with obstacle" & trial_type=="Execution")
d.corr.obs_stop   <- d.corr  %>%  filter(probe=="with obstacle" & trial_type=="Motor Imagery")


# ==== CORRELATION OVERALL MIQ-SCORE with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
#Bayesian
corr.MIQ_Overall.noobs_go    <- correlationBF(d.corr.noobs_go$MIQ_Overall, d.corr.noobs_go$ReachBias, posterior = FALSE)
table.MIQ_Overall.noobs_go   <- describe_posterior(corr.MIQ_Overall.noobs_go, ci_method = "HDI")
table.MIQ_Overall.noobs_go   <- table.MIQ_Overall.noobs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Overall.noobs_stop  <- correlationBF(d.corr.noobs_stop$MIQ_Overall, d.corr.noobs_stop$ReachBias, posterior = FALSE)
table.MIQ_Overall.noobs_stop <- describe_posterior(corr.MIQ_Overall.noobs_stop, ci_method = "HDI")
table.MIQ_Overall.noobs_stop   <- table.MIQ_Overall.noobs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median)

corr.MIQ_Overall.obs_go      <- correlationBF(d.corr.obs_go$MIQ_Overall, d.corr.obs_go$ReachBias, posterior = FALSE)
table.MIQ_Overall.obs_go     <- describe_posterior(corr.MIQ_Overall.obs_go, ci_method = "HDI")
table.MIQ_Overall.obs_go     <- table.MIQ_Overall.obs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median)

corr.MIQ_Overall.obs_stop    <- correlationBF(d.corr.obs_stop$MIQ_Overall, d.corr.obs_stop$ReachBias, posterior = FALSE)
table.MIQ_Overall.obs_stop   <- describe_posterior(corr.MIQ_Overall.obs_stop, ci_method = "HDI")
table.MIQ_Overall.obs_stop   <- table.MIQ_Overall.obs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median)

#combine tables
table.MIQ_Overall <- bind_rows(
  table.MIQ_Overall.noobs_go,
  table.MIQ_Overall.noobs_stop,
  table.MIQ_Overall.obs_go,
  table.MIQ_Overall.obs_stop
) %>%
  as.data.frame() %>%
  add_column(Measure = "MIQ Score Overall") %>%
  relocate(Measure, .before = `Trial Type`) %>%
  mutate(pd = pd * 100) %>%
  # Optional: numeric rounding (keeps columns numeric)
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # Force two-decimal display for the columns you print
  mutate(
    Median  = sprintf("%.2f", Median),
    pd      = sprintf("%.2f", pd),
    CI_low  = sprintf("%.2f", CI_low),
    CI_high = sprintf("%.2f", CI_high)
  ) %>%
  rename(
    r = Median,
    `pd (%)` = pd
  ) %>%
  unite("HDI", CI_low:CI_high, sep = ", ", remove = TRUE) %>%
  mutate(HDI = str_c("[", HDI, "]"))

table.MIQ_Overall



# ==== CORRELATION VISUAL MIQ-SCORE with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
#Bayesian
corr.MIQ_Visual.noobs_go   <- correlationBF(d.corr.noobs_go$MIQ_Vis, d.corr.noobs_go$ReachBias, posterior = FALSE)
table.MIQ_Visual.noobs_go   <- describe_posterior(corr.MIQ_Visual.noobs_go, ci_method = "HDI")
table.MIQ_Visual.noobs_go   <- table.MIQ_Visual.noobs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Visual.noobs_stop   <- correlationBF(d.corr.noobs_stop$MIQ_Vis, d.corr.noobs_stop$ReachBias, posterior = FALSE)
table.MIQ_Visual.noobs_stop   <- describe_posterior(corr.MIQ_Visual.noobs_stop, ci_method = "HDI")
table.MIQ_Visual.noobs_stop   <- table.MIQ_Visual.noobs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Visual.obs_go   <- correlationBF(d.corr.obs_go$MIQ_Vis, d.corr.obs_go$ReachBias, posterior = FALSE)
table.MIQ_Visual.obs_go   <- describe_posterior(corr.MIQ_Visual.obs_go, ci_method = "HDI")
table.MIQ_Visual.obs_go   <- table.MIQ_Visual.obs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Visual.obs_stop   <- correlationBF(d.corr.obs_stop$MIQ_Vis, d.corr.obs_stop$ReachBias, posterior = FALSE)
table.MIQ_Visual.obs_stop   <- describe_posterior(corr.MIQ_Visual.obs_stop, ci_method = "HDI")
table.MIQ_Visual.obs_stop   <- table.MIQ_Visual.obs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

#combine tables
table.MIQ_Visual <- bind_rows(
  table.MIQ_Visual.noobs_go,
  table.MIQ_Visual.noobs_stop,
  table.MIQ_Visual.obs_go,
  table.MIQ_Visual.obs_stop
) %>%
  as.data.frame() %>%
  add_column(Measure = "MIQ Score Visual") %>%
  relocate(Measure, .before = `Trial Type`) %>%
  mutate(pd = pd * 100) %>%
  # Optional: numeric rounding (keeps columns numeric)
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # Force two-decimal display for the columns you print
  mutate(
    Median  = sprintf("%.2f", Median),
    pd      = sprintf("%.2f", pd),
    CI_low  = sprintf("%.2f", CI_low),
    CI_high = sprintf("%.2f", CI_high)
  ) %>%
  rename(
    r = Median,
    `pd (%)` = pd
  ) %>%
  unite("HDI", CI_low:CI_high, sep = ", ", remove = TRUE) %>%
  mutate(HDI = str_c("[", HDI, "]"))

table.MIQ_Visual


# ==== CORRELATION KINESTHETIC MIQ-SCORE with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
#Bayesian
corr.MIQ_Kin.noobs_go   <- correlationBF(d.corr.noobs_go$MIQ_Kin, d.corr.noobs_go$ReachBias, posterior = FALSE)
table.MIQ_Kin.noobs_go   <- describe_posterior(corr.MIQ_Kin.noobs_go, ci_method = "HDI")
table.MIQ_Kin.noobs_go   <- table.MIQ_Kin.noobs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Kin.noobs_stop   <- correlationBF(d.corr.noobs_stop$MIQ_Kin, d.corr.noobs_stop$ReachBias, posterior = FALSE)
table.MIQ_Kin.noobs_stop   <- describe_posterior(corr.MIQ_Kin.noobs_stop, ci_method = "HDI")
table.MIQ_Kin.noobs_stop   <- table.MIQ_Kin.noobs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Kin.obs_go   <- correlationBF(d.corr.obs_go$MIQ_Kin, d.corr.obs_go$ReachBias, posterior = FALSE)
table.MIQ_Kin.obs_go   <- describe_posterior(corr.MIQ_Kin.obs_go, ci_method = "HDI")
table.MIQ_Kin.obs_go   <- table.MIQ_Kin.obs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MIQ_Kin.obs_stop   <- correlationBF(d.corr.obs_stop$MIQ_Kin, d.corr.obs_stop$ReachBias, posterior = FALSE)
table.MIQ_Kin.obs_stop   <- describe_posterior(corr.MIQ_Kin.obs_stop, ci_method = "HDI")
table.MIQ_Kin.obs_stop   <- table.MIQ_Kin.obs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

#combine tables
table.MIQ_Kin <- bind_rows(
  table.MIQ_Kin.noobs_go,
  table.MIQ_Kin.noobs_stop,
  table.MIQ_Kin.obs_go,
  table.MIQ_Kin.obs_stop
) %>%
  as.data.frame() %>%
  add_column(Measure = "MIQ Score Kinesthetic") %>%
  relocate(Measure, .before = `Trial Type`) %>%
  mutate(pd = pd * 100) %>%
  # Optional: numeric rounding (keeps columns numeric)
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # Force two-decimal display for the columns you print
  mutate(
    Median  = sprintf("%.2f", Median),
    pd      = sprintf("%.2f", pd),
    CI_low  = sprintf("%.2f", CI_low),
    CI_high = sprintf("%.2f", CI_high)
  ) %>%
  rename(
    r = Median,
    `pd (%)` = pd
  ) %>%
  unite("HDI", CI_low:CI_high, sep = ", ", remove = TRUE) %>%
  mutate(HDI = str_c("[", HDI, "]"))

table.MIQ_Kin



# ==== CORRELATION EASE OF MI with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
#Bayesian
corr.MI_Ease.noobs_go <- correlationBF(d.corr.noobs_go$MI_Ease_Overall, d.corr.noobs_go$ReachBias, posterior = FALSE)
table.MI_Ease.noobs_go   <- describe_posterior(corr.MI_Ease.noobs_go, ci_method = "HDI")
table.MI_Ease.noobs_go   <- table.MI_Ease.noobs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MI_Ease.noobs_stop <- correlationBF(d.corr.noobs_stop$MI_Ease_Overall, d.corr.noobs_stop$ReachBias, posterior = FALSE)
table.MI_Ease.noobs_stop   <- describe_posterior(corr.MI_Ease.noobs_stop, ci_method = "HDI")
table.MI_Ease.noobs_stop   <- table.MI_Ease.noobs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MI_Ease.obs_go <- correlationBF(d.corr.obs_go$MI_Ease_Overall, d.corr.obs_go$ReachBias, posterior = FALSE)
table.MI_Ease.obs_go   <- describe_posterior(corr.MI_Ease.obs_go, ci_method = "HDI")
table.MI_Ease.obs_go   <- table.MI_Ease.obs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MI_Ease.obs_stop    <- correlationBF(d.corr.obs_stop$MI_Ease_Overall, d.corr.obs_stop$ReachBias, posterior = FALSE)
table.MI_Ease.obs_stop   <- describe_posterior(corr.MI_Ease.obs_stop, ci_method = "HDI")
table.MI_Ease.obs_stop   <- table.MI_Ease.obs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 


#combine tables
table.MI_Ease <- bind_rows(
  table.MI_Ease.noobs_go,
  table.MI_Ease.noobs_stop,
  table.MI_Ease.obs_go,
  table.MI_Ease.obs_stop
) %>%
  as.data.frame() %>%
  add_column(Measure = "Ease of Imagery") %>%
  relocate(Measure, .before = `Trial Type`) %>%
  mutate(pd = pd * 100) %>%
  # Optional: numeric rounding (keeps columns numeric)
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # Force two-decimal display for the columns you print
  mutate(
    Median  = sprintf("%.2f", Median),
    pd      = sprintf("%.2f", pd),
    CI_low  = sprintf("%.2f", CI_low),
    CI_high = sprintf("%.2f", CI_high)
  ) %>%
  rename(
    r = Median,
    `pd (%)` = pd
  ) %>%
  unite("HDI", CI_low:CI_high, sep = ", ", remove = TRUE) %>%
  mutate(HDI = str_c("[", HDI, "]"))

table.MI_Ease


# ==== CORRELATION COUNT OF MI with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
#Bayesian
corr.MI_Count.noobs_go <- correlationBF(d.corr.noobs_go$MI_Count_Overall, d.corr.noobs_go$ReachBias, posterior = FALSE)
table.MI_Count.noobs_go   <- describe_posterior(corr.MI_Count.noobs_go, ci_method = "HDI")
table.MI_Count.noobs_go   <- table.MI_Count.noobs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MI_Count.noobs_stop <- correlationBF(d.corr.noobs_stop$MI_Count_Overall, d.corr.noobs_stop$ReachBias, posterior = FALSE)
table.MI_Count.noobs_stop   <- describe_posterior(corr.MI_Count.noobs_stop, ci_method = "HDI")
table.MI_Count.noobs_stop   <- table.MI_Count.noobs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "without obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MI_Count.obs_go <- correlationBF(d.corr.obs_go$MI_Count_Overall, d.corr.obs_go$ReachBias, posterior = FALSE)
table.MI_Count.obs_go   <- describe_posterior(corr.MI_Count.obs_go, ci_method = "HDI")
table.MI_Count.obs_go   <- table.MI_Count.obs_go %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Execution", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 

corr.MI_Count.obs_stop    <- correlationBF(d.corr.obs_stop$MI_Count_Overall, d.corr.obs_stop$ReachBias, posterior = FALSE)
table.MI_Count.obs_stop   <- describe_posterior(corr.MI_Count.obs_stop, ci_method = "HDI")
table.MI_Count.obs_stop   <- table.MI_Count.obs_stop %>% select(c(Median,CI_low,CI_high,pd))  %>% 
  add_column(`Trial Type` = "Motor Imagery", Probe = "with obstacle")  %>% 
  relocate(`Trial Type`, .before = Median)  %>% 
  relocate(Probe, .before = Median) 


#combine tables
table.MI_Count <- bind_rows(
  table.MI_Count.noobs_go,
  table.MI_Count.noobs_stop,
  table.MI_Count.obs_go,
  table.MI_Count.obs_stop
) %>%
  as.data.frame() %>%
  add_column(Measure = "Count of Imagery") %>%
  relocate(Measure, .before = `Trial Type`) %>%
  mutate(pd = pd * 100) %>%
  # Optional: numeric rounding (keeps columns numeric)
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # Force two-decimal display for the columns you print
  mutate(
    Median  = sprintf("%.2f", Median),
    pd      = sprintf("%.2f", pd),
    CI_low  = sprintf("%.2f", CI_low),
    CI_high = sprintf("%.2f", CI_high)
  ) %>%
  rename(
    r = Median,
    `pd (%)` = pd
  ) %>%
  unite("HDI", CI_low:CI_high, sep = ", ", remove = TRUE) %>%
  mutate(HDI = str_c("[", HDI, "]"))

table.MI_Count





# ==== DESIGN AND EXPORT CORRELATION TABLE (SUPPLEMENTARY TABLE S4) ====
# combine dataframes
table.MI <- rbind(table.MIQ_Overall,table.MIQ_Visual,table.MIQ_Kin,table.MI_Ease,table.MI_Count)
table.MI

table.MI$r <- as.numeric(table.MI$r)

ftable.MI <- flextable(table.MI)  %>%  
  fontsize(size = 12, part = "header") %>% 
  bold(part = "header") %>% 
  align(align = "center", part = "header") %>% 
  align(j = c("r", "HDI", "pd (%)"), align = "center", part = "all") %>% 
  bg(j="r",bg = scales::col_numeric(palette = "viridis", domain = c(-0.3, 0.3))) %>% 
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 2),part = "header" ) %>% 
  hline_top(border = fp_border(color = "black", width = 2) ) %>% 
  hline_bottom(border = fp_border(color = "black", width = 2) ) %>% 
  hline(i = c(4,8,12,16), border = fp_border(color = "black", width = 1.5)) %>%
  set_caption("Table S4. Correlations between Reach Bias and Motor Imagery Questionnaire Scores (Experiment 3)") 

ftable.MI <- autofit(ftable.MI)
ftable.MI

# export to doc
doc <- read_docx()
doc <- body_add_flextable(doc, ftable.MI)
print(doc, target = "TableS4_A069_Table_MIQ_ReachBias_Correlations.docx")




#==== CORRELATION PLOTS ====
# ==== CORRELATION PLOT: OVERALL MIQ-SCORE with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
# probe without obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Overall MIQ Score"
titleY      <- ""

EXE_r_idx   <- formatC( table.MIQ_Overall$r[table.MIQ_Overall$`Trial Type`   =="Execution"     & table.MIQ_Overall$Probe=="without obstacle"], format = "f", digits = 2)
EXE_HDI_idx <- formatC( table.MIQ_Overall$HDI[table.MIQ_Overall$`Trial Type` =="Execution"     & table.MIQ_Overall$Probe=="without obstacle"], format = "f", digits = 2)
MI_r_idx    <- formatC( table.MIQ_Overall$r[table.MIQ_Overall$`Trial Type`   =="Motor Imagery" & table.MIQ_Overall$Probe=="without obstacle"], format = "f", digits = 2)
MI_HDI_idx  <- formatC( table.MIQ_Overall$HDI[table.MIQ_Overall$`Trial Type` =="Motor Imagery" & table.MIQ_Overall$Probe=="without obstacle"], format = "f", digits = 2)

g.corr_NoObs_MIQ <- ggplot(
  data = d.corr %>% filter(probe == "without obstacle"),
  aes(x = MIQ_Overall, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_NoObs_MIQ


# probe with obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Overall MIQ Score"
titleY      <- "Reach Bias (°)"

EXE_r_idx   <- formatC( table.MIQ_Overall$r[table.MIQ_Overall$`Trial Type`   =="Execution"     & table.MIQ_Overall$Probe=="with obstacle"], format = "f", digits = 2)
EXE_HDI_idx <- formatC( table.MIQ_Overall$HDI[table.MIQ_Overall$`Trial Type` =="Execution"     & table.MIQ_Overall$Probe=="with obstacle"], format = "f", digits = 2)
MI_r_idx    <- formatC( table.MIQ_Overall$r[table.MIQ_Overall$`Trial Type`   =="Motor Imagery" & table.MIQ_Overall$Probe=="with obstacle"], format = "f", digits = 2)
MI_HDI_idx  <- formatC( table.MIQ_Overall$HDI[table.MIQ_Overall$`Trial Type` =="Motor Imagery" & table.MIQ_Overall$Probe=="with obstacle"], format = "f", digits = 2)

g.corr_Obs_MIQ <- ggplot(
  data = d.corr %>% filter(probe == "with obstacle"),
  aes(x = MIQ_Overall, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_Obs_MIQ


# ==== CORRELATION PLOT: VISUAL MIQ-SCORE with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
# probe without obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Visual MIQ Score"
titleY      <- ""

EXE_r_idx   <- formatC( table.MIQ_Visual$r[table.MIQ_Visual$`Trial Type`  =="Execution" & table.MIQ_Visual$Probe=="without obstacle"], format = "f", digits = 2)
EXE_HDI_idx <- formatC( table.MIQ_Visual$HDI[table.MIQ_Visual$`Trial Type` =="Execution" & table.MIQ_Visual$Probe=="without obstacle"], format = "f", digits = 2)
MI_r_idx    <- formatC( table.MIQ_Visual$r[table.MIQ_Visual$`Trial Type`   =="Motor Imagery" & table.MIQ_Visual$Probe=="without obstacle"], format = "f", digits = 2)
MI_HDI_idx  <- formatC( table.MIQ_Visual$HDI[table.MIQ_Visual$`Trial Type` =="Motor Imagery" & table.MIQ_Visual$Probe=="without obstacle"], format = "f", digits = 2)


g.corr_NoObs_MIQ_vis <- ggplot(
  data = d.corr %>% filter(probe == "without obstacle"),
  aes(x = MIQ_Vis, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_NoObs_MIQ_vis



# probe with obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Visual MIQ Score"
titleY      <- "Reach Bias (°)"

EXE_r_idx   <- formatC( table.MIQ_Visual$r[table.MIQ_Visual$`Trial Type`   =="Execution"     & table.MIQ_Visual$Probe=="with obstacle"], format = "f", digits = 2)
EXE_HDI_idx <- formatC( table.MIQ_Visual$HDI[table.MIQ_Visual$`Trial Type` =="Execution"     & table.MIQ_Visual$Probe=="with obstacle"], format = "f", digits = 2)
MI_r_idx    <- formatC( table.MIQ_Visual$r[table.MIQ_Visual$`Trial Type`   =="Motor Imagery" & table.MIQ_Visual$Probe=="with obstacle"], format = "f", digits = 2)
MI_HDI_idx  <- formatC( table.MIQ_Visual$HDI[table.MIQ_Visual$`Trial Type` =="Motor Imagery" & table.MIQ_Visual$Probe=="with obstacle"], format = "f", digits = 2)


g.corr_Obs_MIQ_vis <- ggplot(
  data = d.corr %>% filter(probe == "with obstacle"),
  aes(x = MIQ_Vis, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 44, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_Obs_MIQ_vis


# ==== CORRELATION PLOT: KINESTHETIC MIQ-SCORE with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
# probe without obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Kinesthetic MIQ Score"
titleY      <- ""

EXE_r_idx   <- formatC( table.MIQ_Kin$r[table.MIQ_Kin$`Trial Type`   =="Execution"     & table.MIQ_Kin$Probe=="without obstacle"], format = "f", digits = 2)
EXE_HDI_idx <- formatC( table.MIQ_Kin$HDI[table.MIQ_Kin$`Trial Type` =="Execution"     & table.MIQ_Kin$Probe=="without obstacle"], format = "f", digits = 2)
MI_r_idx    <- formatC( table.MIQ_Kin$r[table.MIQ_Kin$`Trial Type`   =="Motor Imagery" & table.MIQ_Kin$Probe=="without obstacle"], format = "f", digits = 2)
MI_HDI_idx  <- formatC( table.MIQ_Kin$HDI[table.MIQ_Kin$`Trial Type` =="Motor Imagery" & table.MIQ_Kin$Probe=="without obstacle"], format = "f", digits = 2)


g.corr_NoObs_MIQ_kin <- ggplot(
  data = d.corr %>% filter(probe == "without obstacle"),
  aes(x = MIQ_Kin, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_NoObs_MIQ_kin


# probe with obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Kinesthetic MIQ Score"
titleY      <- "Reach Bias (°)"

EXE_r_idx   <- formatC( table.MIQ_Kin$r[table.MIQ_Kin$`Trial Type`   =="Execution"     & table.MIQ_Kin$Probe=="with obstacle"], format = "f", digits = 2)
EXE_HDI_idx <- formatC( table.MIQ_Kin$HDI[table.MIQ_Kin$`Trial Type` =="Execution"     & table.MIQ_Kin$Probe=="with obstacle"], format = "f", digits = 2)
MI_r_idx    <- formatC( table.MIQ_Kin$r[table.MIQ_Kin$`Trial Type`   =="Motor Imagery" & table.MIQ_Kin$Probe=="with obstacle"], format = "f", digits = 2)
MI_HDI_idx  <- formatC( table.MIQ_Kin$HDI[table.MIQ_Kin$`Trial Type` =="Motor Imagery" & table.MIQ_Kin$Probe=="with obstacle"], format = "f", digits = 2)


g.corr_Obs_MIQ_kin <- ggplot(
  data = d.corr %>% filter(probe == "with obstacle"),
  aes(x = MIQ_Kin, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_Obs_MIQ_kin



# ==== CORRELATION PLOT EASE OF MI with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
# probe without obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Ease of Motor Imagery"
titleY      <- ""

EXE_r_idx   <- formatC( table.MI_Ease$r[table.MI_Ease$`Trial Type`   =="Execution"     & table.MI_Ease$Probe=="without obstacle"], format ="f", digits=2)
EXE_HDI_idx <- formatC( table.MI_Ease$HDI[table.MI_Ease$`Trial Type` =="Execution"     & table.MI_Ease$Probe=="without obstacle"], format ="f", digits=2)
MI_r_idx    <- formatC( table.MI_Ease$r[table.MI_Ease$`Trial Type`   =="Motor Imagery" & table.MI_Ease$Probe=="without obstacle"], format ="f", digits=2)
MI_HDI_idx  <- formatC( table.MI_Ease$HDI[table.MI_Ease$`Trial Type` =="Motor Imagery" & table.MI_Ease$Probe=="without obstacle"], format ="f", digits=2)


g.corr_NoObs_MIQ_Ease <- ggplot(
  data = d.corr %>% filter(probe == "without obstacle"),
  aes(x = MI_Ease_Overall, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_NoObs_MIQ_Ease


# probe with obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Ease of Motor Imagery"
titleY      <- "Reach Bias (°)"

EXE_r_idx   <- formatC( table.MI_Ease$r[table.MI_Ease$`Trial Type`   =="Execution"     & table.MI_Ease$Probe=="with obstacle"], format ="f", digits=2)
EXE_HDI_idx <- formatC( table.MI_Ease$HDI[table.MI_Ease$`Trial Type` =="Execution"     & table.MI_Ease$Probe=="with obstacle"], format ="f", digits=2)
MI_r_idx    <- formatC( table.MI_Ease$r[table.MI_Ease$`Trial Type`   =="Motor Imagery" & table.MI_Ease$Probe=="with obstacle"], format ="f", digits=2)
MI_HDI_idx  <- formatC( table.MI_Ease$HDI[table.MI_Ease$`Trial Type` =="Motor Imagery" & table.MI_Ease$Probe=="with obstacle"], format ="f", digits=2)

g.corr_Obs_MIQ_Ease <- ggplot(
  data = d.corr %>% filter(probe == "with obstacle"),
  aes(x = MI_Ease_Overall, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 4, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 4, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(3, 7),
    breaks = 3:7
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_Obs_MIQ_Ease


# ==== CORRELATION PLOT COUNT OF MI with REACH BIAS (Obstacle Prime - No Obstacle Prime) ====
# probe without obstacle
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Count of Motor Imagery"
titleY      <- ""

EXE_r_idx   <- formatC( table.MI_Count$r[table.MI_Count$`Trial Type`   =="Execution"     & table.MI_Count$Probe=="without obstacle"], format ="f", digits=2)
EXE_HDI_idx <- formatC( table.MI_Count$HDI[table.MI_Count$`Trial Type` =="Execution"     & table.MI_Count$Probe=="without obstacle"], format ="f", digits=2)
MI_r_idx    <- formatC( table.MI_Count$r[table.MI_Count$`Trial Type`   =="Motor Imagery" & table.MI_Count$Probe=="without obstacle"], format ="f", digits=2)
MI_HDI_idx  <- formatC( table.MI_Count$HDI[table.MI_Count$`Trial Type` =="Motor Imagery" & table.MI_Count$Probe=="without obstacle"], format ="f", digits=2)

g.corr_NoObs_MIQ_Count <- ggplot(
  data = d.corr %>% filter(probe == "without obstacle"),
  aes(x = MI_Count_Overall, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 1.5, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 1.5, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(1, 3),
    breaks = 1:3
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_NoObs_MIQ_Count


# obs probe
colorValues <- c(color_exe1,color_mi1)
fontSize    <- 9
titleX      <- "Mean Count of Motor Imagery"
titleY      <- "Reach Bias (°)"

EXE_r_idx   <- formatC( table.MI_Count$r[table.MI_Count$`Trial Type`   =="Execution"     & table.MI_Count$Probe=="with obstacle"], format ="f", digits=2)
EXE_HDI_idx <- formatC( table.MI_Count$HDI[table.MI_Count$`Trial Type` =="Execution"     & table.MI_Count$Probe=="with obstacle"], format ="f", digits=2)
MI_r_idx    <- formatC( table.MI_Count$r[table.MI_Count$`Trial Type`   =="Motor Imagery" & table.MI_Count$Probe=="with obstacle"], format ="f", digits=2)
MI_HDI_idx  <- formatC( table.MI_Count$HDI[table.MI_Count$`Trial Type` =="Motor Imagery" & table.MI_Count$Probe=="with obstacle"], format ="f", digits=2)

g.corr_Obs_MIQ_Count <- ggplot(
  data = d.corr %>% filter(probe == "with obstacle"),
  aes(x = MI_Count_Overall, y = ReachBias, color = trial_type, fill = trial_type)
) +
  facet_wrap(~probe, strip.position = "top") + 
  geom_point(
    position = position_jitter(width = 0, height = 0), # no dodge/jitter needed
    size = 3, alpha = 0.7, show.legend = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
  annotate("text",
           label = paste("R =", EXE_r_idx, EXE_HDI_idx),
           x = 1.5, y = 45, size = 4, colour = color_exe1) +
  annotate("text",
           label = paste("R =", MI_r_idx, MI_HDI_idx),
           x = 1.5, y = 38, size = 4, colour = color_mi1) +
  scale_color_manual(values = colorValues) +
  scale_fill_manual(values = colorValues) +
  scale_x_continuous(
    name = titleX,
    limits = c(1, 3),
    breaks = 1:3
  ) +
  scale_y_continuous(
    name = titleY,
    limits = c(-10, 50),
    breaks = seq(-10, 50, by = 10)
  ) +
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

g.corr_Obs_MIQ_Count




# ==== MAKE PLOT (5x2 plot grid, library patchwork) (SUPPLEMENTARY FIGURE S10) ====
g.corr_MI <- (
  g.corr_Obs_MIQ + g.corr_NoObs_MIQ +
    g.corr_Obs_MIQ_vis + g.corr_NoObs_MIQ_vis +
    g.corr_Obs_MIQ_kin + g.corr_NoObs_MIQ_kin +
    g.corr_Obs_MIQ_Ease + g.corr_NoObs_MIQ_Ease +
    g.corr_Obs_MIQ_Count + g.corr_NoObs_MIQ_Count
) +
  plot_layout(
    ncol = 2,                 # arrange as 2 columns x 5 rows (10 plots)
    widths = c(1, 1),         # two columns -> two widths
    guides = "collect"        # single shared legend
  ) +
  plot_annotation(tag_levels = "A") &  # auto labels A, B, C, ...
  theme(legend.position = "bottom")

g.corr_MI


# Build a valid path
outfile <- file.path(figurePath, "FigS10_A069_Correlation_MIQ_ReachBias.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.corr_MI,
  width = 24, height = 40, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS10_A069_Correlation_MIQ_ReachBias.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.corr_MI,
  width = 24, height = 40, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "FigS10_A069_Correlation_MIQ_ReachBias.png")

ggsave(
  filename = outfile,
  plot = g.corr_MI,
  width = 24, height = 40, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)





#==== TOTAL RESPONSE TIME PRIME (TRT) ====
# Calculation of TRT:
# for MI trials equal to MIT = time from prime target onset until participants pressed mouse button to indicate finish of imagery
# for execution trials = RT + MT + MT_back

# ==== MODEL FITTING: PRIME TOTAL RESPONSE TIME: 2-FACTOR INTERACTION MODEL WITH SHIFTED LOG-NORMAL DISTRIBUTION ====
summary(dClean$TRT_Prime) #6 NA's present, total of 50 trials removed here
dClean.TRT       <- subset(dClean, (!is.na(TRT_Prime) & TRT_Prime > 200 & TRT_Prime < 5000))

dClean.TRT %>%
  ggplot( aes(x = TRT_Prime, fill = obstacle_prime) ) +
  facet_wrap(~stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 100, boundary = 0, position = 'identity' ) +
  labs( title = "A069 - Histogram", subtitle = "RT Probe" ) +
  scale_x_continuous( name = 'RT [ms]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )


summary(dClean.TRT$TRT_Prime)
sd(log(dClean.TRT$TRT_Prime))


# looking at prior values
get_prior(bf(TRT_Prime ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),
          family = shifted_lognormal(),
          data   = dClean.TRT)


# define priors
# the linear predictor is on the log of the shifted RT, i.e., y=log(RT-ndt). SO the Intercept corresponds to the typical value of (RT-ndt)
# Intercept: median exp(7), given median RT = 1380 and ndt = 250, Intercept should be ~1380-250 ~1130 ~ 7
# Fixed effects: normal(0, 0.2), class = "b") implies a 1 SD effect of about exp(+-0.2) ~ x[0.82,1.22]
# Group-level SD: reasonable weakly-informative choice.
# Residual log-SD: somewhat smaller than the observed sd(log RT) ≈ 0.4
# Shift: Centering near 200 is pragmatic given min RT ≈ 213
# Correlations among random effects
prior_TRT.Prime      <- c(prior(normal( 7,    0.3  ), class = "Intercept"), 
                          prior(normal( 0,     0.1  ), class = "b"),
                          prior(normal( 0,     0.1),   class = "sd"),
                          prior(normal( 0,     0.2),   class = "sigma"),
                          prior(normal(200,    30),    class = "ndt",  lb = 0),  
                          prior(lkj(3),                class = "cor")
)



# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 3 hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A069_fit.TRT_Prime_OPRIxSS <- brm( 
  bf(TRT_Prime ~ obstacle_prime * stop_signal_prime + ( obstacle_prime * stop_signal_prime | subID )),     # model specification
  data   = dClean.TRT,                # data
  family = shifted_lognormal(),       # distribution of the response variable
  prior  = prior_TRT.Prime,           # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A069_fit.TRT_Prime_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )













#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(A069_fit.TRT_Prime_OPRIxSS,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(A069_fit.TRT_Prime_OPRIxSS,   type = "hist", ndraws = 10)           # looks good
pp_check(A069_fit.TRT_Prime_OPRIxSS,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(A069_fit.TRT_Prime_OPRIxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(A069_fit.TRT_Prime_OPRIxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)




#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(dClean$obstacle_prime),
                        stop_signal_prime  = levels(dClean$stop_signal_prime))


# Posterior draws of expected values (population-level, no random effects)
TRT_Prime.posteriors <- as.data.frame(fitted(
  A069_fit.TRT_Prime_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_SS = interaction(obstacle_prime, stop_signal_prime, sep = "_", drop = TRUE))
colnames(TRT_Prime.posteriors) <- exp.cond_posterior$OPRI_SS


# Long format (use pivot_longer)
TRT_Prime.posteriors_long <- TRT_Prime.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_SS", values_to = "value") %>%
  separate(OPRI_SS, into = c("obstacle_prime", "stop_signal_prime"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
TRT.EMM <- TRT_Prime.posteriors_long %>%
  group_by(obstacle_prime, stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
TRT.EMM



# Summaries: median and 95% HDI for Execution and Motor Imagery
TRT_Prime.posteriors_long %>%
  group_by(stop_signal_prime) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
dClean$subID    <- droplevels(dClean$subID)
exp.cond.subj   <- expand.grid(subID              = levels(dClean$subID),
                               obstacle_prime     = levels(dClean$obstacle_prime),
                               stop_signal_prime  = levels(dClean$stop_signal_prime))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  A069_fit.TRT_Prime_OPRIxSS,
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
TRT.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, stop_signal_prime) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)



# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S9A) ====
# Prepare and tidy data
# rename factors
TRT.EMM <- TRT.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime)) # remove old columns
TRT.EMM

# rename factors
TRT.subj.EMM <- TRT.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime)) # remove old columns



###
colorValues <- c(color_exe2,color_mi2,color_exe2,color_mi2)
fontSize    <- 9
titleX      <- ""
titleY      <- "Prime duration (ms)"

g.TRT.EMM.A069 <-
  TRT.subj.EMM %>%
  ggplot(aes(x = trial_type,
             y = .value,
             color = prime:trial_type,
             shape = prime)) +
  
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
    data = dplyr::filter(TRT.EMM),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(TRT.EMM),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(400, 2500),
    breaks = seq(500, 2500, 500)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_Neut, color_Neut)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    legend.position = "top",
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.TRT.EMM.A069






# ==== TOTAL RESPONSE TIME: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
TRT.Contrasts           <- data.frame(matrix(ncol = 4, nrow = 40000))
colnames(TRT.Contrasts) <- c("Execution", "Motor Imagery", "EX vs MI", 
                             "EX minus MI")

# Execution trials (obs prime - noobs prime)
TRT.Contrasts$`Execution`             <- (TRT_Prime.posteriors$`yes_go` - TRT_Prime.posteriors$`no_go`)
# MI trials (obs prime - noobs prime)
TRT.Contrasts$`Motor Imagery`         <- (TRT_Prime.posteriors$`yes_stop` - TRT_Prime.posteriors$`no_stop`)
# EX minus MI (pooled over prime)
TRT.Contrasts$`EX vs MI`            <- (TRT_Prime.posteriors$`yes_go` + TRT_Prime.posteriors$`no_go`)/2 - (TRT_Prime.posteriors$`yes_stop` + TRT_Prime.posteriors$`no_stop`)/2 
# Diff in obs prime effect btw Ex and MI
TRT.Contrasts$`EX minus MI`     <- (TRT.Contrasts$`Execution` - TRT.Contrasts$`Motor Imagery`)


TRT.Contrasts_long <- pivot_longer(TRT.Contrasts, cols = everything(),
                                   names_to = "contrast", values_to = "value")
# Define the desired order
contrast_order <- c(
  "Execution", "Motor Imagery", "EX vs MI", "EX minus MI"
)

# Convert 'contrast' to a factor with this order
TRT.Contrasts_long$contrast <- factor(TRT.Contrasts_long$contrast, levels = contrast_order)

# Check
head(TRT.Contrasts_long)



#get point summary (median) and interval (95% hdci) and Probability of direction (pd) of contrasts
TRT.contrasts.summary <-
  TRT.Contrasts_long %>%
  group_by(contrast) %>%
  summarise(
    median = median(value),
    lower  = hdci(value, .width = 0.95)[1],
    upper  = hdci(value, .width = 0.95)[2],
    pd     = bayestestR::p_direction(value)$pd*100,
  ) 

TRT.contrasts.summary$pd <- format(TRT.contrasts.summary$pd, nsmall = 1)  
#print(TRT.contrasts.summary, n = Inf, width = Inf)


#### Calculate ROPE
( RR       <- rope_range(A069_fit.TRT_Prime_OPRIxSS) )


# Calculate percent in ROPE for contrasts
options(digits=3)
TRT.contrast_in_ROPE                       <- as.data.frame(TRT.contrasts.summary)
TRT.contrast_in_ROPE$lowerROPE             <- RR[1]
TRT.contrast_in_ROPE$upperROPE             <- RR[2]
TRT.contrast_in_ROPE$CI_range              <- TRT.contrast_in_ROPE$upper - TRT.contrast_in_ROPE$lower
TRT.contrast_in_ROPE$minUpper              <- TRT.contrast_in_ROPE %>% select(4, 7) %>%  apply(1, FUN=min)  # calculate min of upper limits
TRT.contrast_in_ROPE$maxLower              <- TRT.contrast_in_ROPE %>% select(3, 6) %>%  apply(1, FUN=max) # calculate max of lower limits
TRT.contrast_in_ROPE$DiffminUppermaxLower  <- TRT.contrast_in_ROPE$minUpper  - TRT.contrast_in_ROPE$maxLower 
TRT.contrast_in_ROPE$Zeros                 <- rep(0,nrow(TRT.contrast_in_ROPE))
TRT.contrast_in_ROPE$Overlap               <- TRT.contrast_in_ROPE  %>% select(11,12) %>%  apply(1, FUN=max) # calculate overlap
TRT.contrast_in_ROPE$perc_in_ROPE          <- (TRT.contrast_in_ROPE$Overlap*100)/TRT.contrast_in_ROPE$CI_range
TRT.contrast_in_ROPE[,c(1:7,14)]


# subject-wise contrasts: Overall difference Execution vs Motor Imagery
TRT.ExvsMI.subj.contrast <- TRT.subj.EMM %>%
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
TRT.subj.contrast <- TRT.subj.EMM %>%
  select(subID, trial_type, prime, .value) %>%
  tidyr::pivot_wider(
    names_from = prime,
    values_from = .value
  ) %>%
  mutate(diff = `with obstacle` - `without obstacle`,
         contrast = case_when(
           trial_type == "Execution"  ~ "Execution",
           trial_type == "Motor Imagery"  ~ "Motor Imagery")
  )

#print(TRT.subj.contrast)


# subject-wise contrasts: Execution (Prime with obstacle minus Prime without obstacle) minus Motor Imagery (Prime with obstacle minus Prime without obstacle) 
TRT.subj.diffContrast <- TRT.subj.contrast %>%
  select(subID, trial_type, diff) %>%
  pivot_wider(
    names_from = trial_type,
    values_from = diff
  ) %>%
  mutate(Execution_vs_MI = Execution - `Motor Imagery`,
         contrast = "EX minus MI"
  )

#print(TRT.subj.diffContrast)



#  ==== TOTAL RESPONSE TIME PRIME: PLOTTING CONTRASTS   (SUPPLEMENTARY FIGURE S9BCD) ====
# Contrasts: Prime with obstacle minus Prime without obstacle
titleX <- expression(Delta~Prime~duration~(ms))
titleY <- ""
g.TRT.Contrast.A069 <- TRT.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "Execution" = "Execution",
    "Motor Imagery" = "Motor Imagery"
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
  scale_fill_manual(values = c(color_exe2,color_mi2)) + 
  scale_pattern_fill_manual(values = c(color_exe2,color_mi2)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(TRT.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.005, ymax = 0.04, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(TRT.subj.contrast),
             aes(y=0, x=diff,color=trial_type),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.001,
                                             nudge.y = -0.0025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe2,color_mi2)) +
  scale_x_continuous(name= titleX, limits=c(-50,200),breaks=c(-50,0,50,100,150,200)) + 
  scale_y_continuous(limits = c(-0.005, 0.04)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.TRT.Contrast.A069


### Difference in Execution minus Difference Motor Imagery ###
titleX <- expression(Delta~Prime~duration~Difference(ms))
titleY <- "EX minus MI"
g.TRT_EXvsMI.Contrast.A069 <- TRT.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
  ggplot(aes(x = value, fill = contrast, pattern_fill = contrast)) +
  facet_wrap(~contrast, strip.position = "left", ncol = 1,labeller = as_labeller(c(
    "EX minus MI" = "EX minus MI"
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
  scale_fill_manual(values = c(color_exe_mi2)) + 
  scale_pattern_fill_manual(values = c(color_exe_mi2)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(TRT.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.005, ymax = 0.04, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(TRT.subj.diffContrast),
             aes(y=0, x=Execution_vs_MI,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.001,
                                             nudge.y = -0.0025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe_mi2)) +
  scale_x_continuous(name= titleX, limits=c(-100,150),breaks=c(-100,-50,0,50,100,150)) + 
  scale_y_continuous(limits = c(-0.005, 0.04)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.TRT_EXvsMI.Contrast.A069




### Execution minus Motor Imagery overall ###
titleX <- expression(Delta~Prime~duration~(ms))
titleY <- "EX minus MI"
g.TRT_EXvsMI_overall.Contrast.A069 <- TRT.Contrasts_long %>% filter(contrast=="EX vs MI") %>% 
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
  scale_fill_manual(values = c(color_exe2,color_mi2)) + 
  scale_pattern_fill_manual(values = c(color_mi2,color_mi2)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(TRT.contrasts.summary, contrast == "EX vs MI"),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.0005, ymax = 0.006, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  geom_point(data=dplyr::filter(TRT.ExvsMI.subj.contrast, contrast == "EX vs MI"),
             aes(y=0, x=diff,color=contrast),
             position = position_jitternudge(jitter.width =  0,
                                             jitter.height = 0.0001,
                                             nudge.y = -0.00025,
                                             nudge.x = 0,
                                             seed = NA),       
             size = 2,
             show.legend=FALSE) +
  scale_color_manual(labels = c(""),values = c(color_exe2,color_mi2)) +
  scale_x_continuous(name= titleX, limits=c(-500,1200),breaks=c(-500,0,500,1000)) + 
  scale_y_continuous(limits = c(-0.0005, 0.006)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.TRT_EXvsMI_overall.Contrast.A069






# ==== CORRELATION OF PRIME DURATION OF EXECUTION WITH MOTOR IMAGERY TRIALS  (SUPPLEMENTARY FIGURE S9E)====
# Execution vs MI Overall Prime duration
d.corr_TRT <- dClean  %>% 
  group_by(subID, stop_signal_prime) %>% 
  summarise(medianTRT   = median(TRT_Prime, na = T)
  )

#make wide format
d.corr_TRT            <- d.corr_TRT    %>% pivot_wider(names_from = c(stop_signal_prime), values_from = c(medianTRT))

#Bayesian
corr.EXEvsMI_Overall <- correlationBF(d.corr_TRT$go, d.corr_TRT$stop, posterior = FALSE)
describe_posterior(corr.EXEvsMI_Overall)

## Execution vs MI Prime duration, separate for obstacle and no obstacle 
d.corr_TRT2 <- dClean  %>% 
  group_by(subID, obstacle_prime, stop_signal_prime) %>% 
  summarise(medianTRT   = median(TRT_Prime, na = T)
  )

#make wide format
d.corr_TRT2            <- d.corr_TRT2    %>% pivot_wider(names_from = c(stop_signal_prime), values_from = c(medianTRT))

d.corr_TRT2_obs     <- d.corr_TRT2  %>%  filter(obstacle_prime=="yes")
d.corr_TRT2_noobs   <- d.corr_TRT2  %>%  filter(obstacle_prime=="no")

#Bayesian
#obs prime
corr.EXEvsMI_obs <- correlationBF(d.corr_TRT2_obs$go, d.corr_TRT2_obs$stop, posterior = FALSE)
describe_posterior(corr.EXEvsMI_obs)
#no obs prime
corr.EXEvsMI_noobs <- correlationBF(d.corr_TRT2_noobs$go, d.corr_TRT2_noobs$stop, posterior = FALSE)
describe_posterior(corr.EXEvsMI_noobs)



## Execution vs MI Difference in Prime duration 
d.corr_TRT3 <- dClean  %>% 
  group_by(subID, obstacle_prime, stop_signal_prime) %>% 
  summarise(medianTRT   = median(TRT_Prime, na = T)
  )

d.corr_TRT3           <- d.corr_TRT3    %>% pivot_wider(names_from = c(obstacle_prime), values_from = c(medianTRT))
d.corr_TRT3$TRT_Diff  <- d.corr_TRT3$yes - d.corr_TRT3$no #TRT_Diff = difference between prime with obstacle and without obstacle
d.corr_TRT3           <- d.corr_TRT3    %>%  select(-c("no","yes"))
d.corr_TRT3           <- d.corr_TRT3    %>% pivot_wider(names_from = c(stop_signal_prime), values_from = c(TRT_Diff))
#Bayesian
corr.EXEvsMI_Diff <- correlationBF(d.corr_TRT3$go, d.corr_TRT3$stop, posterior = FALSE)
( table.corr.EXEvsMI_Diff <- describe_posterior(corr.EXEvsMI_Diff) )


table.corr.EXEvsMI_Diff$pd <- table.corr.EXEvsMI_Diff$pd *100
table.corr.EXEvsMI_Diff <- table.corr.EXEvsMI_Diff %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  rename("r" = Median, "pd (%)" = pd) %>% unite("HDI",CI_low:CI_high) %>% 
  mutate(HDI = str_c("[", HDI, "]")) %>%
  mutate(HDI = str_replace(HDI, "_", ","))
table.corr.EXEvsMI_Diff

#
EXEvsMI_Diff_r_idx   <- formatC( table.corr.EXEvsMI_Diff$r, format = "f", digits = 2)
EXEvsMI_Diff_HDI_idx <- formatC( table.corr.EXEvsMI_Diff$HDI, format = "f", digits = 2)


#plots 
colorValues <- c(color_exe_mi1)
fontSize    <- 9
titleX      <- expression(Delta~Prime~duration~Execution~(ms))
titleY      <- expression(Delta~Prime~duration~Motor~Imagery~(ms))


g.trt.corr <- ggplot( data = d.corr_TRT3, aes(x=go , y=stop) ) + 
  geom_point(aes(colour = color_exe_mi1),position = position_jitterdodge(jitter.width = 0.0,
                                                                         jitter.height = 0,
                                                                         dodge.width = .0,
                                                                         seed = NA), size = 3, alpha = .7 , show.legend=FALSE) + 
  geom_smooth(aes(colour = color_exe_mi1),method = "lm", se = TRUE, show.legend=FALSE) +
  annotate("text", label = paste("R = ", EXEvsMI_Diff_r_idx, " ", EXEvsMI_Diff_HDI_idx),x = -20, y = 250, size = 4, colour = "black") + 
  scale_color_manual(values = colorValues ) + 
  scale_fill_manual(values = colorValues ) + 
  scale_x_continuous(name= titleX, limits=c(-100, 310),breaks=c(-100,0,100,200,300)) + 
  scale_y_continuous(name= titleY, limits=c(-100, 310),breaks=c(-100,0,100,200,300)) + 
  theme( axis.title.x = element_text(size=fontSize ) ) + 
  theme( axis.text = element_text( size = fontSize ) ) + 
  theme( axis.title.y = element_text(size=fontSize ) ) + 
  theme_cowplot() + custom_plot_theme + 
  theme(legend.title = element_blank()) +
  theme(
    plot.title = element_text(size = 12, face = "bold")
  )
g.trt.corr




# ==== PRIME TRT: MAKE ONE FIGURE  (SUPPLEMENTARY FIGURE S9) =====
design <- "ABCE
           ABDE"

p.trt.A069 <-
  (g.TRT.EMM.A069 +
     g.TRT.Contrast.A069 +
     g.TRT_EXvsMI.Contrast.A069 +
     g.TRT_EXvsMI_overall.Contrast.A069 +
     g.trt.corr) +
  plot_layout(widths = c(1.5, 1, 1, 1.5), design = design) +
  plot_annotation(tag_levels = "A")

p.trt.A069



# Build a valid path
outfile <- file.path(figurePath, "Figs9_A069_TRT_Prime_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = p.trt.A069,
  width = 48, height = 18, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "Figs9_A069_TRT_Prime_ALL.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = p.trt.A069,
  width = 48, height = 18, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "Figs9_A069_TRT_Prime_ALL.png")

ggsave(
  filename = outfile,
  plot = p.trt.A069,
  width = 48, height = 18, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)



