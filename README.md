# Title: Distinct contributions of motor imagery and execution to history-dependent biases in reaching
# DOI: https://doi.org/10.64898/2026.04.17.719269
# Authors: Seegelke, Heed
# Scripts by Christian Seegelke, 01/09/2026
# =========================================================================================================================================================================

# Contains the following Experiments:
# ================================================================================================

EXPERIMENT 1: A061
EXPERIMENT 2: A054
EXPERIMENT 3: A069
RE-ANALYSIS OF DATA FROM Roberts, J.W., Wakefield, C.J., & Owen, R. (2025). Trajectory priming through obstacle avoidance in motor imagery – does motor imagery comprise the spatial characteristics of movement? Experimental Brain Research, 243:9

# Data availability: The preprocessed datasets required to reproduce the analyses are available on Zenodo: DOI: 10.5281/zenodo.22251643.

# Scripts
# ================================================================================================

A061_Analysis_v1.3.R
- INPUT: A061_data.csv (preprocessed data), A061_SubjInfo.xlsx (Demographics & MI Questionnaire data)
- OUTPUT: Reported Stats of Experiment 1; Figure 3; Supplementary Table S1; Supplementary Figure S1


A061_A054_Comparison_v1.3.R
- INPUT: A061_data.csv, A054_data.csv (preprocessed data), A061_SubjInfo.xlsx (Demographics & MI Questionnaire data), A054_SubjInfo.xlsx (Demographics)
- OUTPUT: Reported Stats of Experiments 1 and 2; Figure 2, Figure 4, Supplementary Figures S2, S3, S4, S5


A069_Analysis_v1.3.R
- INPUT: A069_data.csv (preprocessed data), A069_SubjInfo.xlsx (Demographics & MI Questionnaire data)
- OUTPUT: Reported Stats of Experiment 3; Supplementary Figures S6, S7, S8, S9, S10, S11, S12, S13, Supplementary Table S2, S3, S4


A069_Roberts2025_ReAnalysis.R
- INPUT: Data_Roberts_2025_EBR.xlsx
- OUTPUT: Supplementary Figures S14


A069_Roberts2025_A061_A054_RT_ReAnalysis.R
- INPUT: Data_Roberts_2025_EBR.xlsx, A061_data.csv, A054_data.csv, A069_data.cs
- OUTPUT: Figure 5, Figure 6


A061_Plotting_Figure1.m
INPUT: A061_Data_complete_S18.mat
OUTPUT: Parts of Figure 1B


A054_Plotting_Figure1.m
INPUT: A054_Data_complete_S02.mat
OUTPUT: Parts of Figure 1B


# ================================================================================================
NOTE: PLEASE ADJUST DATA PATHS ACCORDINGLY

