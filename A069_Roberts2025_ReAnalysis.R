# Title: Distinct contributions of motor imagery and execution to history-dependent biases in reaching
# Re-Analysis of Data from Roberts, J.W., Wakefield, C.J., & Owen, R. (2025).
# Trajectory priming through obstacle avoidance in motor imagery – does motor imagery comprise the spatial characteristics of movement?
# Experimental Brain Research, 243:9
# Authors: Seegelke, Heed 
# Script: Christian Seegelke 01/09/2026
# INPUT:  Data_Roberts_2025_EBR.xlsx
# OUTPUT: REPORTED STATS OF ROBERTS ET AL. 2025 (BAYESIAN), SUPPLEMENTARY FIGURE S14

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

#==== SANITY CHECK 1 ====
# Generic summariser: counts and within-subject proportions
summarise_by_subject <- function(df, var, id = subID) {
  var <- rlang::enquo(var)
  id  <- rlang::enquo(id)
  
  df %>%
    filter(!is.na(!!var)) %>%
    count(!!id, !!var, name = "Count") %>%
    group_by(!!id) %>%
    mutate(Freq = Count / sum(Count)) %>%
    ungroup()
}

# Plotters
plot_counts <- function(san, x = subID, fill) {
  fill <- rlang::enquo(fill)
  ggplot(san, aes(x = !!rlang::enquo(x), y = Count, fill = !!fill)) +
    geom_col(position = position_dodge(width = 0.9), color = "black") +
    labs(y = "Count")
}


plot_props <- function(san, x = subID, fill) {
  fill <- rlang::enquo(fill)
  ggplot(san, aes(x = !!rlang::enquo(x), y = Freq, fill = !!fill)) +
    geom_col(position = position_fill(reverse = FALSE), color = "black") +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Proportion")
}


# trial_type_str
san_stop <- summarise_by_subject(d, trial_type_str)
plot_counts(san_stop, fill = trial_type_str)
plot_props(san_stop,  fill = trial_type_str)

# obstacle_prime
san_obst_pr <- summarise_by_subject(d, obstacle_prime)
plot_counts(san_obst_pr, fill = obstacle_prime)
plot_props(san_obst_pr,  fill = obstacle_prime)

# obstacle_probe
san_obst_pb <- summarise_by_subject(d, obstacle_probe)
plot_counts(san_obst_pb, fill = obstacle_probe)
plot_props(san_obst_pb,  fill = obstacle_probe)


# For binary factors coded as "0/1" or "no/yes", proportions are more informative
san_err_stop <- summarise_by_subject(d, Error)
plot_props(san_err_stop, fill = Error)



#==== REMOVE ERROR TRIALS ====
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




#==== DATA EXPLORATION =====
# Preflight: ensure numeric, clean factors
num_cols <- c(
  "RT_Probe","MT_Probe",
  "InitialReachError_Probe","PeakHeight_Probe",
  "RT_Prime","MT_Prime",
  "InitialReachError_Prime",
  "PeakHeight_Prime"
)

to_num <- function(x) suppressWarnings(as.numeric(x))

d_plot <- d %>%
  mutate(across(all_of(num_cols), to_num),
         obstacle_prime = droplevels(obstacle_prime),
         obstacle_probe = droplevels(obstacle_probe),
         trial_type_str = droplevels(trial_type_str))


# Palette sized to obstacle_prime levels
pal_obst <- hue_pal()(nlevels(d$obstacle_prime))

# Legend placement that works on ggplot2 < 3.5 too
theme_legend_inside <- function() {
  # If "inside" is not supported, silently fall back to numeric coords
  tryCatch(
    theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85)),
    error = function(e) theme(legend.position = c(0.85, 0.85))
  )
}

plot_metric <- function(df, metric, bw, xlab, title_prefix,
                        normalize = TRUE,          # renamed from `density`
                        clip_quantiles = NULL,
                        angle_limits = NULL,
                        show_points = TRUE) {
  m <- ensym(metric)
  df2 <- df %>% filter(!is.na(!!m))
  
  # Choose y mapping for histogram outside aes()
  y_map <- if (normalize) {
    aes(y = after_stat(density))
  } else {
    aes(y = after_stat(count))
  }
  
  p_hist <- ggplot(df2, aes(x = !!m, fill = obstacle_prime, color = obstacle_prime)) +
    facet_grid(rows = vars(trial_type_str), cols = vars(obstacle_probe)) +
    geom_histogram(mapping = y_map, binwidth = bw, boundary = 0,
                   position = "identity", alpha = 0.25, na.rm = TRUE) +
    { if (normalize) geom_density(adjust = 1, linewidth = 0.8, alpha = 0.2, na.rm = TRUE) } +
    labs(
      title = paste(title_prefix, "distributions"),
      x = xlab,
      y = if (normalize) "Density" else "Count"
    ) +
    scale_fill_manual(values = pal_obst, name = "Prime obstacle") +
    scale_color_manual(values = pal_obst, guide = "none") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")) +
    theme_legend_inside()
  
  if (!is.null(clip_quantiles)) {
    lims <- quantile(df2[[as_string(m)]], clip_quantiles, na.rm = TRUE)
    p_hist <- p_hist + coord_cartesian(xlim = lims)
  }
  if (!is.null(angle_limits)) {
    p_hist <- p_hist + coord_cartesian(xlim = angle_limits)
  }
  
  p_box <- ggplot(df2, aes(x = obstacle_probe, y = !!m, fill = obstacle_prime)) +
    facet_wrap(~ trial_type_str) +
    geom_boxplot(width = 0.6, coef = 3, outlier.alpha = 0.25, na.rm = TRUE) +
    { if (show_points) geom_jitter(aes(color = obstacle_prime),
                                   width = 0.15, height = 0,
                                   alpha = 0.12, size = 0.6, show.legend = FALSE) } +
    labs(
      title = paste(title_prefix, "by obstacle conditions"),
      x = "Obstacle (probe phase)", y = xlab
    ) +
    scale_fill_manual(values = pal_obst, name = "Prime obstacle") +
    scale_color_manual(values = pal_obst, guide = "none") +
    theme_minimal(base_size = 11) +
    theme(axis.line = element_line(color = "black"))
  
  if (!is.null(clip_quantiles)) {
    lims <- quantile(df2[[as_string(m)]], clip_quantiles, na.rm = TRUE)
    p_box <- p_box + coord_cartesian(ylim = lims)
  }
  if (!is.null(angle_limits)) {
    p_box <- p_box + coord_cartesian(ylim = angle_limits)
  }
  
  list(hist = p_hist, box = p_box)
}

# ==== PROBE MOVEMENTS ==== 
# RT Probe
plots_rt <- plot_metric(
  d, metric = RT_Probe, bw = 10,
  xlab = "RT [ms]", title_prefix = "Roberts 2025 – RT Probe",
  normalize = TRUE, clip_quantiles = c(0.01, 0.99)
)
plots_rt$hist
plots_rt$box

# MT Probe
plots_mt <- plot_metric(
  d, metric = MT_Probe, bw = 20,
  xlab = "MT [ms]", title_prefix = "Roberts 2025 – MT Probe",
  normalize = TRUE, clip_quantiles = c(0.01, 0.99)
)
plots_mt$hist
plots_mt$box


# Initial Reach Error
plots_ire <- plot_metric(
  d, metric = InitialReachError_Probe, bw = 3,
  xlab = "Initial Absolute Reach Error [°]",
  title_prefix = "Roberts 2025 – Initial Absolute Reach Error Probe",
  normalize = TRUE, angle_limits = c(0, 180)
)
plots_ire$hist
plots_ire$box



# Peak Height
plots_ph<- plot_metric(
  d, metric = PeakHeight_Probe, bw = 2,
  xlab = "Peak Height",
  title_prefix = "Roberts 2025 – PeakHeight_Probe",
  normalize = TRUE, clip_quantiles = c(0.01, 0.99)
)
plots_ph$hist
plots_ph$box 


# ==== PRIME MOVEMENTS ==== 
# RT 
plots_rt <- plot_metric(
  d, metric = RT_Prime, bw = 10,
  xlab = "RT [ms]", title_prefix = "Roberts 2025 – RT Prime",
  normalize = TRUE, clip_quantiles = c(0.01, 0.99)
)
plots_rt$hist
plots_rt$box

# MT 
plots_mt <- plot_metric(
  d, metric = MT_Prime, bw = 20,
  xlab = "MT [ms]", title_prefix = "Roberts 2025 – MT Prime",
  normalize = TRUE, clip_quantiles = c(0.01, 0.99)
)
plots_mt$hist
plots_mt$box


# Initial Reach Error
plots_ire <- plot_metric(
  d, metric = InitialReachError_Prime, bw = 3,
  xlab = "Initial Absolute Reach Error [°]",
  title_prefix = "Roberts 2025 – Initial Absolute Reach Error Prime",
  normalize = TRUE, angle_limits = c(0, 180)
)
plots_ire$hist
plots_ire$box


# Peak Height
plots_ph<- plot_metric(
  d, metric = PeakHeight_Prime, bw = 2,
  xlab = "Peak Height",
  title_prefix = "Roberts 2025 – PeakHeight_Prime",
  normalize = TRUE, clip_quantiles = c(0.01, 0.99)
)
plots_ph$hist
plots_ph$box 




#==== AGGREGATE DATA FOR ANOVA: #2 trial_type (Execution, Motor Imagery) x 2 obstacle prime (yes, no) x 2 obstacle probe (yes,no) RM ANOVA ====
dm.subj <- d %>% 
  group_by(subID, obstacle_prime, obstacle_probe, trial_type_str) %>% 
  summarise(meanRT                 = mean(RT_Probe),
            medianRT               = median(RT_Probe),
            meanRT_Prime           = mean(RT_Prime),
            medianRT_Prime         = median(RT_Prime),
            meanMT                 = mean(MT_Probe),
            medianMT               = median(MT_Probe),
            meanMT_Prime           = mean(MT_Prime),
            medianMT_Prime         = median(MT_Prime),
            meanIRE                = mean(InitialReachError_Probe),
            medianIRE              = median(InitialReachError_Probe),
            meanPeakHeight         = mean(PeakHeight_Probe),
            medianPeakHeight       = median(PeakHeight_Probe),
            meanIRE_Prime          = mean(InitialReachError_Prime, na=T ),
            medianIRE_Prime        = median(InitialReachError_Prime , na=T),
            meanPeakHeight_Prime   = mean(PeakHeight_Prime    , na=T),
            medianPeakHeight_Prime = median(PeakHeight_Prime    , na=T)
  )


#==== ANOVAS ====
# ==== MEDIAN (ABSOLUTE) INITIAL REACH ERROR PROBE  ====
anova.ReachDiff_probe <-  afex::aov_ez(
  id     = "subID",
  dv     = "medianIRE",
  data   = dm.subj, 
  within = c("obstacle_prime", "obstacle_probe","trial_type_str"), 
  anova_table = list(es = "pes")
)

knitr::kable(nice(anova.ReachDiff_probe))


#marginal means for 3-way Interaction
anova.RE_probe.OPRIxOPROxSS <- emmeans(anova.ReachDiff_probe, ~obstacle_prime:obstacle_probe:trial_type_str)
anova.RE_probe.OPRIxOPROxSS

# Quick plot for 3-way Interaction
ggplot(as.data.frame(anova.RE_probe.OPRIxOPROxSS), aes(x=obstacle_probe, y=emmean, color=obstacle_prime)) + 
  facet_wrap(~trial_type_str) +
  geom_point(position = position_dodge(width=0.5), size=3) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position = position_dodge(width=0.5), width=0.2) + 
  ylab(label = "Median Initial Reach Error [°]") + 
  ggtitle(label = "Initial Reach Error Probe") +
  theme_classic() +
  theme(legend.position = "top")

#comparison with fdr corrected p-values
contrast(anova.RE_probe.OPRIxOPROxSS,
         method = "pairwise",
         by = c("obstacle_probe","trial_type_str"),
         adjust = "fdr")





# ==== MEDIAN PEAK HEIGHT PROBE  ====
anova.ph_probe <-  afex::aov_ez(
  id     = "subID",
  dv     = "medianPeakHeight",
  data   = dm.subj, 
  within = c("obstacle_prime", "obstacle_probe","trial_type_str"), 
  anova_table = list(es = "pes")
)

knitr::kable(nice(anova.ph_probe))


#marginal means for 3-way Interaction
anova.ph_probe.OPRIxOPROxSS <- emmeans(anova.ph_probe, ~obstacle_prime:obstacle_probe:trial_type_str)
anova.ph_probe.OPRIxOPROxSS

# Quick plot for 3-way Interaction
ggplot(as.data.frame(anova.ph_probe.OPRIxOPROxSS), aes(x=obstacle_probe, y=emmean, color=obstacle_prime)) + 
  facet_wrap(~trial_type_str) +
  geom_point(position = position_dodge(width=0.5), size=3) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position = position_dodge(width=0.5), width=0.2) + 
  ylab(label = "Median Peak Height") + 
  ggtitle(label = "Peak Height Probe") +
  theme_classic() +
  theme(legend.position = "top")

#comparison with fdr corrected p-values
contrast(anova.ph_probe.OPRIxOPROxSS,
         method = "pairwise",
         by = c("obstacle_probe","trial_type_str"),
         adjust = "fdr")






# ==== MEDIAN RT PROBE ====
anova.rt_probe <-  afex::aov_ez(
  id     = "subID",
  dv     = "medianRT",
  data   = dm.subj, 
  within = c("obstacle_prime", "obstacle_probe","trial_type_str"), 
  anova_table = list(es = "pes")
)

knitr::kable(nice(anova.rt_probe))


#marginal means for 3-way Interaction
anova.RT_probe.OPRIxOPROxSS <- emmeans(anova.rt_probe, ~obstacle_prime:obstacle_probe:trial_type_str)
anova.RT_probe.OPRIxOPROxSS

# Quick plot for 3-way Interaction
ggplot(as.data.frame(anova.RT_probe.OPRIxOPROxSS), aes(x=obstacle_probe, y=emmean, color=obstacle_prime)) + 
  facet_wrap(~trial_type_str) +
  geom_point(position = position_dodge(width=0.5), size=3) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position = position_dodge(width=0.5), width=0.2) + 
  ylab(label = "Median RT [ms]") + 
  ggtitle(label = "RT Probe") +
  theme_classic() +
  theme(legend.position = "top")

#comparison with fdr corrected p-values
contrast(anova.RT_probe.OPRIxOPROxSS,
         method = "pairwise",
         by = c("obstacle_probe","trial_type_str"),
         adjust = "fdr")





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

#####################################################################################################
#========================================= BAYESIAN REGRESSION MODELS ===========================================
######################################################################################################
#==== INITIAL REACH ERROR PROBE ====
#==== MODEL FITTING: INITIAL REACH ERROR PROBE: MODEL FITTING: 3-FACTOR INTERACTION MODEL WITH GAUSSIAN DISTRIBUTION ====
d %>%
  ggplot( aes(x = InitialReachError_Probe, fill = obstacle_prime) ) +
  facet_wrap(~trial_type_str:obstacle_probe) +
  geom_histogram( alpha = 0.3, binwidth = 3, boundary = 0, position = 'identity' ) +
  labs( title = "Roberts 2025 - Histogram", subtitle = "Reach Direction Probe" ) +
  scale_x_continuous( name = 'Reach Direction Difference [°]') +       
  scale_y_continuous( name = 'Count' )+
  theme( panel.background = element_blank(),
         axis.line = element_line( colour = "black" ),
         panel.grid.major = element_line( linewidth = 0.5, linetype = 'solid', colour = "grey90" ) )


# looking at prior values
get_prior(bf(InitialReachError_Probe | trunc(lb = 0, ub = 180) ~ obstacle_prime * obstacle_probe  * trial_type_str + ( 1 | subID )),
          family = gaussian(),
          data   = d)

# define priors
prior_RE.Probe       <- c(prior(normal(   30,     10  ),     class = "Intercept"), 
                          prior(normal(    0,     7.5  ),     class = "b"),
                          prior(student_t(3, 0,  10   ), class = "sd"), 
                          prior(student_t(3, 0,  15   ), class = "sigma"))


# model with orthogonal contrasts (intercepts corresponds to unweighted grand mean)
# we model the mu here, leaving the scale (sigma) fixed. We include random effects for participants
# runs about XX mins
# adjusting sampling behavior #adapt_delta (btw 0.8 (default) and 1) to decrease number of divergent transitions
# when the depth of the tree being evaluated in each iteration is exceeded => may also bias the posterior samples, increase max_treedepth (default = 10)

afex::set_sum_contrasts()
# prior predictive checks for reasonable priors
( startTime <- Sys.time() )
Roberts2025_fit.IRE_Probe_Obs_OPRIxSS <- brm( 
  bf(InitialReachError_Probe  ~ obstacle_prime * obstacle_probe * trial_type_str + ( 1 | subID )),    # model specification
  data   = d,                         # data
  family = gaussian(),                # distribution of the response variable
  prior  = prior_RE.Probe,            # prior distributions of model parameters
  sample_prior = TRUE,                # we also consider the data now and save samples from prior
  iter   = 11000,                     # number of iterations
  warmup = 1000,                      # number of iterations used as warmup
  chains = 4,                         # number of chains
  cores  = 4,                         # number of cores
  save_pars = save_pars(all = TRUE),  # set to true so model outcome can be saved
  control = list(adapt_delta = 0.9, max_treedepth = 10),
  seed =  1234,                       # set seed for reproducability
  file = paste0(analysisPath,'/', 'Roberts2025_fit.IRE_Probe_Obs_OPRIxSS'), #save model
)
endTime <- Sys.time()
( endTime- startTime )


#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(Roberts2025_fit.IRE_Probe_Obs_OPRIxSS,   type = "dens_overlay", ndraws = 100)  # good
pp_check(Roberts2025_fit.IRE_Probe_Obs_OPRIxSS,   type = "hist", ndraws = 10)           # 
pp_check(Roberts2025_fit.IRE_Probe_Obs_OPRIxSS,   type = "boxplot", ndraws = 10)        # 


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(Roberts2025_fit.IRE_Probe_Obs_OPRIxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(Roberts2025_fit.IRE_Probe_Obs_OPRIxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)


#  ==== EXTRACT ESTIMATED MARGINAL MEANS (EMM) ====
# using fitted function from brms
# Conditions
exp.cond <- expand.grid(obstacle_prime     = levels(d$obstacle_prime),
                        obstacle_probe     = levels(d$obstacle_probe),
                        trial_type_str  = levels(d$trial_type_str))


# Posterior draws of expected values (population-level, no random effects)
IRE_Probe.posteriors <- as.data.frame(fitted(
  Roberts2025_fit.IRE_Probe_Obs_OPRIxSS,
  newdata    = exp.cond,
  re_formula = NA,          # ignore random effects, set to NULL to include random effects
  summary    = FALSE,       # extract the full MCMC, set to TRUE to obtain point and interval estimates (specify those as well), see below
))


# Label columns with condition combinations
exp.cond_posterior <- exp.cond %>%
  mutate(OPRI_OPRO_SS = interaction(obstacle_prime, obstacle_probe, trial_type_str, sep = "_", drop = TRUE))
colnames(IRE_Probe.posteriors) <- exp.cond_posterior$OPRI_OPRO_SS


# Long format (use pivot_longer)
IRE_Probe.posteriors_long <- IRE_Probe.posteriors %>%
  pivot_longer(cols = everything(), names_to = "OPRI_OPRO_SS", values_to = "value") %>%
  separate(OPRI_OPRO_SS, into = c("obstacle_prime", "obstacle_probe", "trial_type_str"), sep = "_", remove = TRUE)


# Summaries: median and 95% HDI
IRE.EMM <- IRE_Probe.posteriors_long %>%
  group_by(obstacle_prime, obstacle_probe, trial_type_str) %>%
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)
IRE.EMM


# Conditions
d$subID    <- droplevels(d$subID)
exp.cond.subj   <- expand.grid(subID              = levels(d$subID),
                               obstacle_prime     = levels(d$obstacle_prime),
                               obstacle_probe     = levels(d$obstacle_probe),
                               trial_type_str    = levels(d$trial_type_str))

# Draws of expected values per subject-condition
draws_subj <- fitted(
  Roberts2025_fit.IRE_Probe_Obs_OPRIxSS,
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
IRE.subj.EMM <- draws_subj_df |>
  group_by(subID, obstacle_prime, obstacle_probe, trial_type_str) |>
  ggdist::point_interval(.value = value, .width = 0.95, .point = median, .interval = hdci)




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) (SUPPLEMENTARY FIGURE S14AB) ====
# Prepare and tidy data
# rename factors
IRE.EMM <- IRE.EMM  %>% #as.data.frame()  %>% 
  mutate(
    trial_type  = factor(if_else(trial_type_str    == "Execution", "Execution", "Motor Imagery")),
    prime       = factor(if_else(obstacle_prime    == "yes", "with obstacle", "without obstacle")),
    probe       = factor(if_else(obstacle_probe    == "yes", "with obstacle", "without obstacle"))
  )   %>% 
  select(-c(obstacle_prime,trial_type_str,obstacle_probe)) # remove old columns

IRE.EMM

# rename factors
IRE.subj.EMM <- IRE.subj.EMM  %>% 
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

g.IRE_NoObs.EMM.Rob<-
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
    limits = c(-10, 75),
    breaks = seq(0, 75, 15)
  ) +
  
  xlab(titleX) +
  
  theme_cowplot() +
  custom_plot_theme +
  
  guides(shape = guide_legend(
    override.aes = list(color = c(color_NeutRob, color_NeutRob)) ) ) + # overrides legend color
  
  guides(color = "none") +
  labs(shape = "Prime") +
  
  theme(
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text(size = 12, face = "bold")
  )

g.IRE_NoObs.EMM.Rob



#plot Initial Reach Error (EMM) for probe with obstacle
colorValues <- c(color_exeRob,color_miRob,color_exeRob,color_miRob)
fontSize    <- 9
titleX      <- ""
titleY      <- "Initial Reach Error (°)"

g.IRE_Obs.EMM.Rob <-
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
    limits = c(-10, 75),
    breaks = seq(0, 75, 15)
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

g.IRE_Obs.EMM.Rob
















# ==== INITIAL REACH ERROR: CALCULATE CONTRASTS, PD, AND ROPE ====
# contrasts are calculated as Prime with obstacle minus Prime without obstacle
IRE.Contrasts           <- data.frame(matrix(ncol = 9, nrow = 40000))
colnames(IRE.Contrasts) <- c("Execution", "Motor Imagery", "EX minus MI", 
                             "Execution NoObs", "Execution Obs", "Motor Imagery NoObs", "Motor Imagery Obs",
                             "EX minus MI NoObs", "EX minus MI Obs")

# Execution trials
IRE.Contrasts$`Execution`           <- (IRE_Probe.posteriors$`yes_yes_Execution` + IRE_Probe.posteriors$`yes_no_Execution`)/2 - ( IRE_Probe.posteriors$`no_yes_Execution` + IRE_Probe.posteriors$`no_no_Execution`)/2
# Motor Imagery trials
IRE.Contrasts$`Motor Imagery`       <- (IRE_Probe.posteriors$`yes_yes_Imagery` + IRE_Probe.posteriors$`yes_no_Imagery`)/2 - ( IRE_Probe.posteriors$`no_yes_Imagery` + IRE_Probe.posteriors$`no_no_Imagery`)/2
# Execution vs Motor Imagery
IRE.Contrasts$`EX minus MI`         <- (IRE.Contrasts$`Execution` - IRE.Contrasts$`Motor Imagery`)
# No obstacle probe go trials
IRE.Contrasts$`Execution NoObs`         <- (IRE_Probe.posteriors$`yes_no_Execution` - IRE_Probe.posteriors$`no_no_Execution`)
# obstacle probe go trials
IRE.Contrasts$`Execution Obs`            <- (IRE_Probe.posteriors$`yes_yes_Execution` - IRE_Probe.posteriors$`no_yes_Execution`)
# No obstacle probe stop trials
IRE.Contrasts$`Motor Imagery NoObs`       <- (IRE_Probe.posteriors$`yes_no_Imagery` - IRE_Probe.posteriors$`no_no_Imagery`)
# No obstacle probe stop trials
IRE.Contrasts$`Motor Imagery Obs`          <- (IRE_Probe.posteriors$`yes_yes_Imagery` - IRE_Probe.posteriors$`no_yes_Imagery`)
# No obstacle probe go vs stop trials
IRE.Contrasts$`EX minus MI NoObs`      <- IRE.Contrasts$`Execution NoObs` - IRE.Contrasts$`Motor Imagery NoObs` 
# Obstacle probe go vs stop trials
IRE.Contrasts$`EX minus MI Obs`       <- IRE.Contrasts$`Execution Obs` - IRE.Contrasts$`Motor Imagery Obs`


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
( RR <- rope_range(Roberts2025_fit.IRE_Probe_Obs_OPRIxSS) )




# Calculate percent in ROPE for contrasts
options(digits=3)
IRE.contrast_in_ROPE                       <- as.data.frame(IRE.contrasts.summary)
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



#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS  ====
# Probe without obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_NoObs.Contrast.Rob <- IRE.Contrasts_long %>% filter(contrast=="Execution NoObs" | contrast=="Motor Imagery NoObs") %>% 
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
  scale_fill_manual(values = c(color_exeRob,color_miRob)) + 
  scale_pattern_fill_manual(values = c(color_exeRob,color_miRob)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "Execution NoObs" | contrast=="Motor Imagery NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exeRob,color_miRob)) +
  scale_x_continuous(name= titleX, limits=c(-10,15),breaks=c(-10,-5,0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs.Contrast.Rob


# Probe with obstacle
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE_Obs.Contrast.Rob<- IRE.Contrasts_long %>% filter(contrast=="Execution Obs" | contrast=="Motor Imagery Obs") %>% 
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
  scale_fill_manual(values = c(color_exeRob,color_miRob)) + 
  scale_pattern_fill_manual(values = c(color_exeRob,color_miRob)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "Execution Obs" | contrast=="Motor Imagery Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.3, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe1,color_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,15),breaks=c(-10,-5,0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.3)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") + #remove y-axis
  ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs.Contrast.Rob



### Execution vs Motor Imagery ###
# Probe without obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_NoObs_EXvsMI.Contrast.Rob <- IRE.Contrasts_long %>% filter(contrast=="EX minus MI NoObs") %>% 
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
  scale_fill_manual(values = c(color_exe_miRob)) + 
  scale_pattern_fill_manual(values = c(color_exe_miRob)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "EX minus MI NoObs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.2, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,15),breaks=c(-10,-50,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.2)) +
  theme_cowplot() + custom_plot_theme +
  theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_NoObs_EXvsMI.Contrast.Rob




# Probe with obstacle
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_Obs_EXvsMI.Contrast.Rob <- IRE.Contrasts_long %>% filter(contrast=="EX minus MI Obs") %>% 
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
  scale_fill_manual(values = c(color_exe_miRob)) + 
  scale_pattern_fill_manual(values = c(color_exe_miRob)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "EX minus MI Obs" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.2, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe_mi1)) +
  scale_x_continuous(name= titleX, limits=c(-10,15),breaks=c(-10,-50,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.2)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("with obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_Obs_EXvsMI.Contrast.Rob




#  ==== INITIAL REACH ERROR: PLOTTING CONTRASTS POOLED OVER OBSTACLE PROBE ====
titleX <- "Reach Bias (°)"
titleY <- ""
g.IRE.Contrast2.Rob <- IRE.Contrasts_long %>% filter(contrast=="Execution" | contrast=="Motor Imagery") %>% 
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
  scale_fill_manual(values = c(color_exeRob,color_miRob)) + 
  scale_pattern_fill_manual(values = c(color_exeRob,color_miRob)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "Execution" | contrast=="Motor Imagery" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.4, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exeRob,color_miRob)) +
  scale_x_continuous(name= titleX, limits=c(-10,15),breaks=c(-10,-5,0,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.4)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE.Contrast2.Rob



### Execution vs Motor Imagery ###
titleX <- expression(Delta~Reach~Bias)
titleY <- ""
g.IRE_EXvsMI.Contrast2.Rob <- IRE.Contrasts_long %>% filter(contrast=="EX minus MI") %>% 
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
  scale_fill_manual(values = c(color_exe_miRob)) + 
  scale_pattern_fill_manual(values = c(color_exe_miRob)) +
  guides(fill = "none", pattern_fill = "none") +  # Turn off legends
  geom_pointinterval(data = dplyr::filter(IRE.contrasts.summary, contrast == "EX minus MI" ),
                     aes(x = median, xmin = lower, xmax = upper),
                     position = position_nudge(y=0),
                     size = 5) +
  
  annotate("rect", xmin = RR[1], xmax = RR[2], ymin = -0.01, ymax = 0.4, fill = color_Neut, alpha = .25)  +
  geom_vline(xintercept = 0, linewidth = 1, linetype = 'dashed') +
  scale_color_manual(labels = c(""),values = c(color_exe_miRob)) +
  scale_x_continuous(name= titleX, limits=c(-10,15),breaks=c(-10,-50,5,10,15,20,25)) + 
  scale_y_continuous(limits = c(-0.01, 0.4)) +
  theme_cowplot() + custom_plot_theme +
  #theme(strip.text.y = element_blank()) + #remove facet labels
  labs(y= NULL) + guides(y = "none") +    #remove y-axis
  #ggtitle("without obstacle") + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust =0.5)
  )

g.IRE_EXvsMI.Contrast2.Rob





#merge into a 1x2 plot grid (library patchwork)
g.IRE.Contrast_Pooled.Rob <-
  (g.IRE.Contrast2.Rob + g.IRE_EXvsMI.Contrast2.Rob) +
  plot_layout(
    widths = c(1), heights = c(2,1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.IRE.Contrast_Pooled.Rob







#  ==== INITIAL REACH ERROR: COMBINE PLOTS (SUPPLEMENTARY FIGURE S14) ====
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


g.RE.Rob <-  g.IRE_Obs.EMM.Rob + g.IRE_NoObs.EMM.Rob + 
  g.IRE_Obs.Contrast.Rob + g.IRE_NoObs.Contrast.Rob + 
  g.IRE_Obs_EXvsMI.Contrast.Rob + g.IRE_NoObs_EXvsMI.Contrast.Rob + 
  plot_layout(design = layout, guides='keep') +
  plot_annotation(tag_levels = "A") 
#theme(legend.position = "bottom")
g.RE.Rob



# Build a valid path
outfile <- file.path(figurePath, "FigS14_Roberts2025_InitialReachError_EMM_Contrasts_ALL.pdf")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.RE.Rob,
  width = 24, height = 28, units = "cm",
  device = grDevices::pdf,       # or simply rely on the .pdf extension
  dpi = 300,                     # ignored for pdf, but harmless
  limitsize = FALSE,
  useDingbats = FALSE            # avoids missing symbol issues in some viewers
)

# Build a valid path
outfile <- file.path(figurePath, "FigS14_Roberts2025_InitialReachError_EMM_Contrasts_ALL.svg")

# Save as PDF (base pdf device)
ggsave(
  filename = outfile,
  plot = g.RE.Rob,
  width = 24, height = 28, units = "cm",
  device = "svg"
)

outfile <- file.path(figurePath, "FigS14_Roberts2025_InitialReachError_EMM_Contrasts_ALL.png")

ggsave(
  filename = outfile,
  plot = g.RE.Rob,
  width = 24, height = 28, units = "cm",
  device = "png",          # or device = grDevices::png
  dpi = 300,
  bg = "white",
  type = "cairo",          # improves text rendering on Windows
  limitsize = FALSE
)







#==== REACTION TIME (RT) PROBE ====
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


#  ==== POSTERIOR PREDICTIVE CHECKS ====
pp_check(Roberts2025_fit.RT_Probe_OPRIxOPROxSS,   type = "dens_overlay", ndraws = 100)  # looks good
pp_check(Roberts2025_fit.RT_Probe_OPRIxOPROxSS,   type = "hist", ndraws = 10)           # looks good
pp_check(Roberts2025_fit.RT_Probe_OPRIxOPROxSS,   type = "boxplot", ndraws = 10)        # looks good


#  ==== CONVERGENCE AND HMC DIAGNOSTICS ====
summary(Roberts2025_fit.RT_Probe_OPRIxOPROxSS)

# Trace and density plots for MCMC Samples of all relevant parameters
# chains should mix and look like hairy caterpillars
plot(Roberts2025_fit.RT_Probe_OPRIxOPROxSS,          
     combo = c("dens", "trace"),
     ask   = FALSE
)




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




# ==== PLOTTING ESTIMATED MARGINAL MEANS (EMM) ====
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
titleY      <- "Reaction time (ms)"

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

#merge into a 1x2 plot grid (library patchwork)
g.RT.EMM.Rob <-
  (g.RT_Obs.EMM.Rob + g.RT_NoObs.EMM.Rob) +
  plot_layout(
    widths = c(1, 1.2),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

g.RT.EMM.Rob






#######################################################################################################