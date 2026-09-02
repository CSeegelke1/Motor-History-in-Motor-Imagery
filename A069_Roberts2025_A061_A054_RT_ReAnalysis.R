# Title: Distinct contributions of motor imagery and execution to history-dependent biases in reaching
# Re-Analysis of Data from Roberts, J.W., Wakefield, C.J., & Owen, R. (2025).
# Trajectory priming through obstacle avoidance in motor imagery – does motor imagery comprise the spatial characteristics of movement?
# Experimental Brain Research, 243:9
# Authors: Seegelke, Heed 
# Script: Christian Seegelke 01/09/2026
# INPUT:  Data_Roberts_2025_EBR.xlsx, A061_data.csv, A054_data.csv, A069_data.csv (preprocessed data)
# OUTPUT: FIGURE 5 and FIGURE 6
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



#==== SETTINGS FOR PLOTTING =======
#define some graphical params like themes
custom_plot_theme <- theme(strip.background =element_rect(fill="white", linewidth = 2),
                           strip.text = element_text(size = rel(1), margin = margin(1,5,5,0, "pt")), #in ggplot2 clockwise starting from top: trbl
                           plot.title = element_text(size = rel(1.5)),
                           panel.background = element_blank())

# for probe actions
color_A069        <- "grey80"
color_Neut        <- "#9C9C9C"
color_exeRob      <- "#006000"
color_miRob       <- "#000060"
color_NeutRob     <- "grey40"
color_exe_miRob   <- "#006060"

color_exe1    <- "#008000"
color_mi1     <- "#000080"
color_exe_mi1 <- "#008080"
color_None    <- "#808000"


#========================================= ROBERTS 2025 IMPORT SINGLE TRIAL DATA  ======================================
# ===== Paths (portable) =====
basePath   <- "C:/Experiments/A069_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")





# ===== Load single-trial data =====
# Note: N = 13, 40 trials per subject = 520 trials
# Expecting column 'subID' present
d <- import(file.path(dataPath, "Data_Roberts_2025_EBR.xlsx")) %>%
  as_tibble()


  d <- d %>%
  mutate(
    across(
      c(trial_type_str, obstacle_prime, obstacle_probe,
        Error),
      ~ factor(.)
    ),
    subID = factor(subID)  # after filtering
  )

# ==== REMOVE ERROR TRIALS ====
d_raw <-d

#Missing trials (N = 15, = 2.9% of all trials)
d.Error_Missing <- subset(d, Error ==1)
print( paste( "Total number of trials:", nrow( d ) ) )
print( paste( "Number of missing trials: ", nrow( d.Error_Missing ) ) )
print( paste( "Percent missing Trials:", round( nrow(d.Error_Missing ) / (nrow( d )) * 100, digits = 1 ) ) )


d$Error_ObsPrime <- NA
for( n in 1:nrow(d))
{
  if ( is.na(d$PeakHeight_Prime[[n]] ) )
      d$Error_ObsPrime[[n]] = NA
  
  else if ( d$obstacle_prime[[n]]  ==  "yes" & d$PeakHeight_Prime[[n]] <=20 )
    d$Error_ObsPrime[[n]] = 1
  
  else
    d$Error_ObsPrime[[n]] = 0
}


d$Error_ObsProbe <- NA
for( n in 1:nrow(d))
{
  if ( is.na(d$PeakHeight_Probe[[n]] ) )
    d$Error_ObsProbe[[n]] = NA
  
  else if ( d$obstacle_probe[[n]]  ==  "yes" & d$PeakHeight_Probe[[n]] <=20 )
    d$Error_ObsProbe[[n]] = 1
  
  else
    d$Error_ObsProbe[[n]] = 0
}


#failed obs avoidance, <20mm peak height when obstacle was present (N = 0, = 0% of all trials)
d.Error_Obs <- subset(d, Error_ObsPrime ==1 | Error_ObsProbe ==1)
print( paste( "Total number of trials:", nrow( d ) ) )
print( paste( "Number of obstacle hit trials: ", nrow( d.Error_Obs ) ) )
print( paste( "Percent obstacle hit Trials:", round( nrow(d.Error_Obs ) / (nrow( d )) * 100, digits = 1 ) ) )

#failed start, RT <100ms (N = 3, = 0.6% of all trials)
d.Error_Start <- subset(d, RT_Prime <100 | RT_Probe <100)
print( paste( "Total number of trials:", nrow( d ) ) )
print( paste( "Number of failed start trials: ", nrow( d.Error_Start ) ) )
print( paste( "Percent failed start Trials:", round( nrow(d.Error_Start ) / (nrow( d )) * 100, digits = 1 ) ) )

#remove errors, N = 18 (3.5%)
d       <- d[d$Error==0,]
d       <- d[d$RT_Prime >= 100  & d$RT_Probe >= 100,]


print( paste( "Total number of trials:", nrow( d_raw ) ) )
print( paste( "Number of removed trials: ", nrow( d_raw ) - nrow( d )) )
print( paste( "Percent removed Trials:", 100 - round( nrow(d) / (nrow( d_raw)) * 100, digits = 1 ) ) )




#==== PROBE RT: MODEL FITTING: 3-FACTOR INTERACTION MODEL WITH SHIFTED LOG-NORMAL DISTRIBUTION ====
d %>%
  ggplot( aes(x = RT_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:trial_type_str) +
  geom_histogram( alpha = 0.3, binwidth = 30, boundary = 0, position = 'identity' ) +
  labs( title = "Roberts 2025 - Histogram", subtitle = "RT Probe" ) +
  scale_x_continuous( name = 'RT [ms]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )

summary(d$RT_Probe)
sd(log(d$RT_Probe))


# looking at prior values
get_prior(bf(RT_Probe ~ obstacle_prime * obstacle_probe * trial_type_str + ( 1 | subID )), 
          family = shifted_lognormal(),
          data   = d)


prior_RT.Probe      <- c(prior(normal( 5.7,  0.5  ), class = "Intercept"), 
                         prior(normal( 0,     0.1  ), class = "b"),
                         prior(normal( 0,     0.1),   class = "sd"),
                         prior(normal( 0,     0.1),   class = "sigma"),
                         prior(normal(100,    15),     class = "ndt",  lb = 0)
)

# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about 2 mins
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
Roberts2025_fit.RT_Probe_OPRIxOPROxSS <- brm( 
  bf(RT_Probe ~ obstacle_prime * obstacle_probe * trial_type_str + ( 1 | subID )),    # model specification
  data   = d,                    # data
  family = shifted_lognormal(),       # distribution of the response variable
  prior  = prior_RT.Probe,       # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.99, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'Roberts2025_fit.RT_Probe_OPRIxOPROxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(d$obstacle_prime),
                        obstacle_probe     = levels(d$obstacle_probe),
                        trial_type_str  = levels(d$trial_type_str))


# Posterior draws of expected values (population-level, no random effects)
RT_Probe.posteriors <- as.data.frame(fitted(
  Roberts2025_fit.RT_Probe_OPRIxOPROxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(obstacle_prime, obstacle_probe, trial_type_str, sep = "_", drop = TRUE))
colnames(RT_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
RT_Probe.posteriors_long <- RT_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("obstacle_prime", "obstacle_probe", "trial_type_str"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
RT.EMM <- RT_Probe.posteriors_long %>%
  group_by(obstacle_prime, obstacle_probe, trial_type_str) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
RT.EMM



# Summaries: median and 95% HDI for Execution and Motor Imagery
RT_Probe.posteriors_long %>%
  group_by(trial_type_str) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)


# Conditions
d$subID    <- droplevels(d$subID)
exp.cond.subj   <- expand.grid(subID              = levels(d$subID),
                               obstacle_prime     = levels(d$obstacle_prime),
                               obstacle_probe     = levels(d$obstacle_probe),
                               trial_type_str  = levels(d$trial_type_str))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  Roberts2025_fit.RT_Probe_OPRIxOPROxSS,
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
  group_by(subID, obstacle_prime, obstacle_probe, trial_type_str) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 5DH) ====
# Prepare and tidy data
# rename factors
RT.EMM.Rob <- RT.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(trial_type_str == "Execution", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,trial_type_str,obstacle_probe)) # remove old columns

RT.EMM.Rob

# rename factors
RT.subj.EMM.Rob <- RT.subj.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(trial_type_str == "Execution", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,trial_type_str,obstacle_probe)) # remove old columns



#plot Initial Reach Error (EMM) for probe without obstacle
colorValues <- c(color_exeRob,color_miRob,color_exeRob,color_miRob)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.RT_NoObs.EMM.Rob <-
  RT.subj.EMM.Rob %>%
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
    data = dplyr::filter(RT.EMM.Rob, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.Rob, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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

g.RT_NoObs.EMM.Rob



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exeRob,color_miRob,color_exeRob,color_miRob)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.RT_Obs.EMM.Rob <-
  RT.subj.EMM.Rob %>%
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
    data = dplyr::filter(RT.EMM.Rob, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.Rob, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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

g.RT_Obs.EMM.Rob










#==== INITIAL REACH ERROR AS A FUNCTION OF PROBE RT ====
q_all.Rob <- d %>%
  filter(trial_type_str %in% c("Execution", "Imagery")) %>%
  select(subID,
         trial_type_str,
         obstacle_prime,
         obstacle_probe,
         RT_Probe,
         InitialReachError_Probe) %>%
  
  # drop NA RTs before binning
  filter(!is.na(RT_Probe)) %>%
  # create subject-wise RT deciles within each condition
  group_by(subID, trial_type_str, obstacle_prime, obstacle_probe) %>%
  # optionally enforce a minimum number of trials per subject-condition
  # filter(dplyr::n() >= 20) %>%
  mutate(Mydeciles = ntile(RT_Probe, 10)) %>%
  group_by(subID, trial_type_str, obstacle_prime, obstacle_probe, Mydeciles) %>%
  summarise(
    RTBin    = median(RT_Probe, na.rm = TRUE),
    REEffect = median(InitialReachError_Probe, na.rm = TRUE),
    nTrials  = dplyr::n(),              # for diagnostics
    .groups  = "drop"
  ) %>%
  # add human-readable labels
  mutate(
    trial_type  = if_else(trial_type_str == "Execution", "Execution", "Motor Imagery"),
    prime       = if_else(obstacle_prime == "yes", "with obstacle", "without obstacle"),
    probe       = if_else(obstacle_probe == "yes", "with obstacle", "without obstacle"),
  )


# group-level aggregation
q_all_agg.Rob <- q_all.Rob %>%
  group_by(trial_type, prime, probe, Mydeciles) %>%
  summarise(
    meanBin      = mean(RTBin, na.rm = TRUE),
    medianBin    = median(RTBin, na.rm = TRUE),
    meanEffect   = mean(REEffect, na.rm = TRUE),
    medianEffect = median(REEffect, na.rm = TRUE),
    IQREffect    = IQR(REEffect, na.rm = TRUE),
    seEffect     = sd(REEffect, na.rm = TRUE) / sqrt(dplyr::n()),
    n_subj       = dplyr::n(),  # number of subject-bins contributing
    .groups      = "drop"
  )



#========================================= A069: IMPORT SINGLE TRIAL DATA  ==============================================
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


dClean.A069 <- dClean











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




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 5CG) ====
# Prepare and tidy data
# rename factors
RT.EMM.A069 <- RT.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

RT.EMM.A069

# rename factors
RT.subj.EMM.A069 <- RT.subj.EMM  %>% 
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
  RT.subj.EMM.A069 %>%
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
    data = dplyr::filter(RT.EMM.A069, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A069, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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
titleY      <- ""

g.RT_Obs.EMM.A069 <-
  RT.subj.EMM.A069  %>%
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
    data = dplyr::filter(RT.EMM.A069, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A069, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 800), 
    breaks = seq(200, 800, 100)
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









#==== INITIAL REACH ERROR AS A FUNCTION OF PROBE RT ====
q_all.A069 <- dClean.A069 %>%
  filter(stop_signal_prime %in% c("go", "stop")) %>%
  select(subID,
         stop_signal_prime,
         obstacle_prime,
         obstacle_probe,
         RT_Probe,
         ReachDiff2_Probe) %>%
  
  # drop NA RTs before binning
  filter(!is.na(RT_Probe)) %>%
  # create subject-wise RT deciles within each condition
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  # optionally enforce a minimum number of trials per subject-condition
  # filter(dplyr::n() >= 20) %>%
  mutate(Mydeciles = ntile(RT_Probe, 10)) %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe, Mydeciles) %>%
  summarise(
    RTBin    = median(RT_Probe, na.rm = TRUE),
    REEffect = median(ReachDiff2_Probe, na.rm = TRUE),
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
q_all_agg.A069 <- q_all.A069 %>%
  group_by(trial_type, prime, probe, Mydeciles) %>%
  summarise(
    meanBin      = mean(RTBin, na.rm = TRUE),
    medianBin    = median(RTBin, na.rm = TRUE),
    meanEffect   = mean(REEffect, na.rm = TRUE),
    medianEffect = median(REEffect, na.rm = TRUE),
    IQREffect    = IQR(REEffect, na.rm = TRUE),
    seEffect     = sd(REEffect, na.rm = TRUE) / sqrt(dplyr::n()),
    n_subj       = dplyr::n(),  # number of subject-bins contributing
    .groups      = "drop"
  )






#========================================= A061: IMPORT SINGLE TRIAL DATA  ==============================================
# ===== Paths (portable) =====
basePath   <- "C:/Experiments/A061_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")

# ===== Load single-trial data =====
# Note: N = 32, 704 trials (incl. 64 practise trials) per subject, subject#07 only 640 trial = 22464 trials
# Expecting column 'subID' present
d <- import(file.path(dataPath, "A061_data.csv")) %>%
  as_tibble()


# ===== Load demographics (sheet 2) =====
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

# ===== Merge demographics to trial data =====
d <- d %>%
  left_join(d_SubInfo, by = "subID")



# ===== Load MIQ-RS scores (sheet 3) and compute subscale/overall =====
d_MIQ_raw <- import(file.path(docPath, "A061_SubjInfo.xlsx"), which = 3) %>%
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
d_Ease_MI_raw <- import(file.path(docPath, "A061_SubjInfo.xlsx"), which = 4) %>%
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
d_Count_MI_raw <- import(file.path(docPath, "A061_SubjInfo.xlsx"), which = 5) %>%
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


dClean.A061 <- dClean

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
                               obstacle_probe     = levels(dClean$obstacle_probe),
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




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 5AB)====
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
titleY      <- "Reaction time (ms)"

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
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A061, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A061, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  
  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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















#==== INITIAL REACH ERROR AS A FUNCTION OF PROBE RT ====
q_all.A061 <- dClean %>%
  filter(stop_signal_prime %in% c("go", "stop")) %>%
  select(subID,
         stop_signal_prime,
         obstacle_prime,
         obstacle_probe,
         RT_Probe,
         ReachDiff2_Probe) %>%
  
  # drop NA RTs before binning
  filter(!is.na(RT_Probe)) %>%
  # create subject-wise RT deciles within each condition
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  # optionally enforce a minimum number of trials per subject-condition
  # filter(dplyr::n() >= 20) %>%
  mutate(Mydeciles = ntile(RT_Probe, 10)) %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe, Mydeciles) %>%
  summarise(
    RTBin    = median(RT_Probe, na.rm = TRUE),
    REEffect = median(ReachDiff2_Probe, na.rm = TRUE),
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
q_all_agg.A061 <- q_all.A061 %>%
  group_by(trial_type, prime, probe, Mydeciles) %>%
  summarise(
    meanBin      = mean(RTBin, na.rm = TRUE),
    medianBin    = median(RTBin, na.rm = TRUE),
    meanEffect   = mean(REEffect, na.rm = TRUE),
    medianEffect = median(REEffect, na.rm = TRUE),
    IQREffect    = IQR(REEffect, na.rm = TRUE),
    seEffect     = sd(REEffect, na.rm = TRUE) / sqrt(dplyr::n()),
    n_subj       = dplyr::n(),  # number of subject-bins contributing
    .groups      = "drop"
  )





#========================================= A054: IMPORT SINGLE TRIAL DATA  ==============================================
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

# ==== MODEL FITTING: REACTION TIME: 3-FACTOR INTERACTION MODEL WITH SHIFTED LOG-NORMAL DISTRIBUTION ====
dClean %>%
  ggplot( aes(x = RT_Probe, fill = obstacle_prime) ) +
  facet_wrap(~obstacle_probe:stop_signal_prime) +
  geom_histogram( alpha = 0.3, binwidth = 10, boundary = 0, position = 'identity' ) +
  labs( title = "A054 - Histogram", subtitle = "RT Probe" ) +
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
# Intercept: median exp(5.4, given median RT = 320 and ndt = 85, Intercept should be ~320-85 ~235 ~ 5.4
# Fixed effects: sd = 0.31 ⇒ 1-sigma multiplicative factor exp(±0.1) ≈ ×[0.90, 1.11]
# Group-level SD: reasonable weakly-informative choice.
# Residual log-SD: somewhat smaller than the observed sd(log RT) ≈ 0.254
# Shift: Centering near 90 is pragmatic given min RT ≈ 100.
# Correlations among random effects
prior_RT.Probe      <- c(prior(normal( 5.4,  0.5  ), class = "Intercept"), 
                         prior(normal( 0,     0.1  ), class = "b"),
                         prior(normal( 0,     0.1),   class = "sd"),
                         prior(normal( 0,     0.1),   class = "sigma"),
                         prior(normal(85,    15),     class = "ndt",  lb = 0),  
                         prior(lkj(3),                class = "cor")
)

# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about XXX hours
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)
# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
A054_fit.RT_Probe_OPRIxOPROxSS <- brm( 
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
  control = list(adapt_delta = 0.9, max_treedepth = 15),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'A054_fit.RT_Probe_OPRIxOPROxSS'), #save model
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
  A054_fit.RT_Probe_OPRIxOPROxSS,
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



# Summaries: median and 95% HDI for Execution and No Movement
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
  A054_fit.RT_Probe_OPRIxOPROxSS,
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




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (FIGURE 5BF) ====
# Prepare and tidy data
# rename factors
RT.EMM.A054 <- RT.EMM  %>% 
  mutate(
    trial_type  = factor(if_else(stop_signal_prime == "go", "Execution", "No Movement")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,stop_signal_prime,obstacle_probe)) # remove old columns

RT.EMM.A054

# rename factors
RT.subj.EMM.A054 <- RT.subj.EMM  %>% 
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

g.RT_NoObs.EMM.A054<-
  RT.subj.EMM.A054 %>%
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
    data = dplyr::filter(RT.EMM.A054, probe == "without obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A054, probe == "without obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  

  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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

g.RT_NoObs.EMM.A054



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exe1,color_None,color_exe1,color_None)
fontSize    <- 9
titleX      <- ""
titleY      <- ""

g.RT_Obs.EMM.A054 <-
  RT.subj.EMM.A054 %>%
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
    data = dplyr::filter(RT.EMM.A054, probe == "with obstacle"),
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(.8),
    width = 0.2,
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Group-level points
  geom_point(
    data = dplyr::filter(RT.EMM.A054, probe == "with obstacle"),
    position = position_dodge(.8),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colorValues) +
  

  scale_y_continuous(
    name = titleY,
    limits = c(200, 800),
    breaks = seq(200, 800, 100)
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

g.RT_Obs.EMM.A054






#==== INITIAL REACH ERROR AS A FUNCTION OF PROBE RT ====
q_all.A054 <- dClean %>%
  filter(stop_signal_prime %in% c("go", "stop")) %>%
  select(subID,
         stop_signal_prime,
         obstacle_prime,
         obstacle_probe,
         RT_Probe,
         ReachDiff2_Probe) %>%
  
  # drop NA RTs before binning
  filter(!is.na(RT_Probe)) %>%
  # create subject-wise RT deciles within each condition
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe) %>%
  # optionally enforce a minimum number of trials per subject-condition
  # filter(dplyr::n() >= 20) %>%
  mutate(Mydeciles = ntile(RT_Probe, 10)) %>%
  group_by(subID, stop_signal_prime, obstacle_prime, obstacle_probe, Mydeciles) %>%
  summarise(
    RTBin    = median(RT_Probe, na.rm = TRUE),
    REEffect = median(ReachDiff2_Probe, na.rm = TRUE),
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
q_all_agg.A054 <- q_all.A054 %>%
  group_by(trial_type, prime, probe, Mydeciles) %>%
  summarise(
    meanBin      = mean(RTBin, na.rm = TRUE),
    medianBin    = median(RTBin, na.rm = TRUE),
    meanEffect   = mean(REEffect, na.rm = TRUE),
    medianEffect = median(REEffect, na.rm = TRUE),
    IQREffect    = IQR(REEffect, na.rm = TRUE),
    seEffect     = sd(REEffect, na.rm = TRUE) / sqrt(dplyr::n()),
    n_subj       = dplyr::n(),  # number of subject-bins contributing
    .groups      = "drop"
  )





##########################################################################################################
#######################################################################################################
#==== PLOT RT ACROSS ALL EXPERIMENTS  (FIGURE 5) ==================
# Paths (portable) 
basePath   <- "C:/Experiments/A069_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")



g.rt_probe_all <- g.RT_Obs.EMM.A061 + g.RT_Obs.EMM.A054 + g.RT_Obs.EMM.A069 + g.RT_Obs.EMM.Rob + 
  g.RT_NoObs.EMM.A061 + g.RT_NoObs.EMM.A054 + g.RT_NoObs.EMM.A069 +g.RT_NoObs.EMM.Rob + 
  plot_layout(widths = c(1, 1, 1, 1.2)) +   # relation of x-axes
  plot_annotation(tag_levels = "A") + # adding panel labels
  plot_layout(guides='collect')  & theme(legend.position='bottom')  # only 1 legend
g.rt_probe_all


# Build a valid path
outfile <- file.path(figurePath, "Fig5_A061_A054_A069_Roberts2025_RT_EMM.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.rt_probe_all,
  width = 36, height = 18, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "Fig5_A061_A054_A069_Roberts2025_RT_EMM.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.rt_probe_all,
  width = 36, height = 18, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "Fig5_A061_A054_A069_Roberts2025_RT_EMM.png")

ggsave(
  filename = outfile,
  plot = g.rt_probe_all,
  width = 36, height = 18, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)



#==== INITIAL REACH ERROR AS A FUNCTION OF PROBE RT =====
q_all_agg.Rob$Experiment    <- as.factor("Roberts2025")
q_all_agg.A061$Experiment   <- as.factor("A061")
q_all_agg.A054$Experiment   <- as.factor("A054")
q_all_agg.A069$Experiment   <- as.factor("A069")

#combine data frames
q_all_agg   <- rbind(q_all_agg.A061,q_all_agg.A054,q_all_agg.A069,q_all_agg.Rob)
q_all_agg$trial_type <- as.factor(q_all_agg$trial_type)
q_all_agg$prime <- as.factor(q_all_agg$prime)
q_all_agg$probe <- as.factor(q_all_agg$probe)

#==== INITIAL REACH ERROR AS A FUNCTION OF PROBE RT: PLOTTING (FIGURE 6) ====
# Paths (portable) 
basePath   <- "C:/Experiments/A069_Hand Path Priming Motor Imagery"
dataPath   <- file.path(basePath, "2_data", "4_clean")
docPath    <- file.path(basePath, "2_data", "2_data_documentation")
figurePath <- file.path(basePath, "4_figures")
analysisPath <- file.path(basePath, "3_code", "3_analysis")
colorValues <- c(color_A069,color_A069,color_exeRob,color_NeutRob)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY  <- "Initial Reach Error (°)"

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_GoNoObs <- q_all_agg %>% filter(Experiment=="A069" | Experiment=="Roberts2025") %>% 
  dplyr::filter(trial_type == "Execution", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = Experiment:prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(-5, 40),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_GoNoObs




colorValues <- c(color_A069,color_A069,color_exeRob,color_NeutRob)  # adjust if needed
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_GoObs <- q_all_agg %>% filter(Experiment=="A069" | Experiment=="Roberts2025") %>% 
  dplyr::filter(trial_type == "Execution", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = Experiment:prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(20, 70),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  annotate("text", label = paste("Roberts et al. 2025"), x = 750, y = 30, size = 4, colour = "grey40") +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  theme_cowplot() + custom_plot_theme +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text    = element_text(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) +
  ggtitle("with obstacle")

g.RE_RT_GoObs



colorValues <- c(color_A069,color_A069,color_miRob,color_NeutRob)  # adjust if needed
fontSize    <- 9
titleX      <- "RT Quantiles"
titleY      <- "Initial Reach Error (°)"

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_StopNoObs <- q_all_agg %>% filter(Experiment=="A069" | Experiment=="Roberts2025") %>% 
  dplyr::filter(trial_type == "Motor Imagery", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = Experiment:prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(-5, 40),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
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

g.RE_RT_StopNoObs



colorValues <- c(color_A069,color_A069,color_miRob,color_NeutRob)  # adjust if needed
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_StopObs <- q_all_agg %>% filter(Experiment=="A069" | Experiment=="Roberts2025") %>% 
  dplyr::filter(trial_type == "Motor Imagery", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = Experiment:prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(20, 70),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
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

g.RE_RT_StopObs


colorValues <- c(color_exe1,color_Neut)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY  <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_GoNoObs.A061 <- q_all_agg %>% filter(Experiment=="A061") %>% 
  dplyr::filter(trial_type == "Execution", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(-5, 40),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_GoNoObs.A061



colorValues <- c(color_exe1,color_Neut)  # adjust if needed
fontSize    <- 9
titleX      <- ""
titleY      <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_GoObs.A061 <- q_all_agg %>% filter(Experiment=="A061") %>% 
  dplyr::filter(trial_type == "Execution", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(20, 70),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  annotate("text", label = paste("Experiment 1"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_GoObs.A061




colorValues <- c(color_mi1,color_Neut)  # adjust if needed
fontSize <- 9
titleX <- "RT Quantiles"
titleY  <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_StopNoObs.A061 <- q_all_agg %>% filter(Experiment=="A061") %>% 
  dplyr::filter(trial_type == "Motor Imagery", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = TRUE) +
  guides(shape = guide_legend(override.aes = list(color = c(color_Neut, color_Neut)) ) ) +  # overrides legend color
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(-5, 40),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_StopNoObs.A061



colorValues <- c(color_mi1,color_Neut)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY  <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_StopObs.A061 <- q_all_agg %>% filter(Experiment=="A061") %>% 
  dplyr::filter(trial_type == "Motor Imagery", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(20, 70),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_StopObs.A061




colorValues <- c(color_exe1,color_Neut)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY  <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_GoNoObs.A054 <- q_all_agg %>% filter(Experiment=="A054") %>% 
  dplyr::filter(trial_type == "Execution", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(-5, 40),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_GoNoObs.A054



colorValues <- c(color_exe1,color_Neut)  # adjust if needed
fontSize    <- 9
titleX      <- ""
titleY      <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_GoObs.A054 <- q_all_agg %>% filter(Experiment=="A054") %>% 
  dplyr::filter(trial_type == "Execution", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(20, 70),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  annotate("text", label = paste("Experiment 2"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_GoObs.A054




colorValues <- c(color_None,color_Neut)  # adjust if needed
fontSize <- 9
titleX <- "RT Quantiles"
titleY  <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_StopNoObs.A054 <- q_all_agg %>% filter(Experiment=="A054") %>% 
  dplyr::filter(trial_type == "No Movement", probe == "without obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(-5, 40),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_StopNoObs.A054



colorValues <- c(color_None,color_Neut)  # adjust if needed
fontSize <- 9
titleX <- ""
titleY  <- ""

pos_dodge <- position_dodge(width = 0.5)

g.RE_RT_StopObs.A054 <- q_all_agg %>% filter(Experiment=="A054") %>% 
  dplyr::filter(trial_type == "No Movement", probe == "with obstacle") %>%
  ggplot(aes(x = meanBin, y = meanEffect, color = prime, shape = prime, group = prime)) +
  facet_wrap(~trial_type, strip.position = "bottom") + 
  geom_errorbar(
    aes(ymin = meanEffect - seEffect, ymax = meanEffect + seEffect),
    position = pos_dodge, width = 10, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(position = pos_dodge, size = 4, show.legend = FALSE) +
  scale_color_manual(values = colorValues) +
  scale_y_continuous(name = titleY, limits = c(20, 70),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
  scale_x_continuous(name = titleX, limits = c(150, 950),
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900)) +
  #annotate("text", label = paste("Experiment 3"), x = 700, y = 35, size = 4, colour = "grey80") +
  #annotate("text", label = paste("Roberts et al. 2025"), x = 745, y = 30, size = 4, colour = "grey40") +
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

g.RE_RT_StopObs.A054





g.re_rt_probe_all <- 
  g.RE_RT_GoObs   + g.RE_RT_GoObs.A061   + g.RE_RT_GoObs.A054   +
  g.RE_RT_GoNoObs + g.RE_RT_GoNoObs.A061 + g.RE_RT_GoNoObs.A054 +
  g.RE_RT_StopObs   + g.RE_RT_StopObs.A061   + g.RE_RT_StopObs.A054   +
  g.RE_RT_StopNoObs + g.RE_RT_StopNoObs.A061 + g.RE_RT_StopNoObs.A054 +
  plot_layout(ncol = 3, guides = "keep", widths = c(1, 1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.re_rt_probe_all


# Build a valid path
outfile <- file.path(figurePath, "Fig6_A061_A054_A069_Roberts2025_IRE_RT_EMM.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.re_rt_probe_all,
  width = 30, height = 30, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)


# Build a valid path
outfile <- file.path(figurePath, "Fig6_A061_A054_A069_Roberts2025_IRE_RT_EMM.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.re_rt_probe_all,
  width = 30, height = 30, units = "cm",
  device = "svg"
)


outfile <- file.path(figurePath, "Fig6_A061_A054_A069_Roberts2025_IRE_RT_EMM.png")

ggsave(
  filename = outfile,
  plot = g.re_rt_probe_all,
  width = 30, height = 30, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)


#==== MEAN MOVEMENT TIMES ========
describe(dClean.A061$MT_Probe)
describe(dClean.A054$MT_Probe)
describe(dClean.A069$MT_Probe)
