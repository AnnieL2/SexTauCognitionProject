## Load libraries and data ----
library(dplyr)
library(nlme) 
library(sjPlot)
library(ggplot2)
library(emmeans)
library(clubSandwich)
library(gt)
library(gtsummary)
library(flextable)
library(patchwork)
library(ggtext)
library(interactions)

path <- '~/Desktop/sexcog/Figures/WRAP/'
load('~/Desktop/sexcog/Preprocessing/WRAP2025-06-15.RData')

# Demographics ----
WRAP_cross_sectional <-  WRAP %>%
  group_by(ID) %>%
  arrange(ID, pacc_time) %>%
  slice_tail()

demographics_WRAP <- WRAP_cross_sectional %>%
  ungroup() %>%
  dplyr::select(SEX, APOE, YrsEd, bl_tau_age, centiloid, pib_group, pacc_time, pacc_time_from_pacc,closest_PACC, PACC_tp, One_PACC_tp) %>%
  tbl_summary(
    by = SEX,  # Stratify by SEX
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",  # Mean (SD) for continuous variables
      all_categorical() ~ "{n} ({p}%)"    # n (%) for categorical variables
    ),
    digits = list(
      all_continuous() ~ 1  # Round continuous variables to 1 decimal place
    ),
    label = list(
      bl_tau_age ~ "Baseline tau-PET age (yrs)",
      APOE ~ "APOE ε4 carriers",
      YrsEd ~ "Education (yrs)",
      centiloid ~ "Baseline Centiloid (CL)",
      pacc_time ~ "Max PACC Time from baseline tau scan (yrs)",
      closest_PACC ~ "Closest PACC to tau scan",
      One_PACC_tp ~ "One PACC time point",
      PACC_tp ~ "PACC timepoints",
      pib_group ~ "Aβ Status"
    )
  ) %>%
  add_p() %>%
  modify_table_styling(
    columns = "p.value",
    rows = p.value < 0.05,
    text_format = "bold"
  ) %>%
  bold_labels()

demographics_WRAP

# Prospective and retrospective PACC stats
extract_pacc_stats <- function (data_frame, retro=TRUE) {
  if (retro) {
    updated_df <- subset(data_frame, pacc_time < 0)
  } else {
    updated_df <- subset(data_frame, pacc_time >= 0)
  }
  updated_df %>%
    dplyr::select(SEX, pacc_time) %>%
    tbl_summary(
      by = SEX,  # Group by sex
      statistic = list(pacc_time ~ "{mean} ({sd})"),  # Show mean (SD)
      digits = all_continuous() ~ 2  # Round numbers to 2 decimal places
    )  %>%
    add_p()
}
extract_pacc_stats(WRAP)
extract_pacc_stats(WRAP, FALSE)

## Main effect of sex -----
WRAP_cross_sectional$SEX <- relevel(WRAP_cross_sectional$SEX, ref = "M")
lm_SEX_Amyg <- lm(data = WRAP_cross_sectional, bl_Amyg ~ SEX +bl_tau_age)
summary(lm_SEX_Amyg)
tab_model(lm_SEX_Amyg, show.std = TRUE)

lm_SEX_Para <- lm(data = WRAP_cross_sectional, bl_Para ~ SEX +bl_tau_age)
summary(lm_SEX_Para)
tab_model(lm_SEX_Para, show.std = TRUE)

lm_SEX_EC <- lm(data = WRAP_cross_sectional, bl_EC ~ SEX +bl_tau_age)
summary(lm_SEX_EC)
tab_model(lm_SEX_EC, show.std = TRUE)

lm_SEX_Fus <- lm(data = WRAP_cross_sectional, bl_Fus ~ SEX +bl_tau_age)
summary(lm_SEX_Fus)
tab_model(lm_SEX_Fus, show.std = TRUE)

lm_SEX_ITG <- lm(data = WRAP_cross_sectional, bl_ITG ~ SEX +bl_tau_age)
summary(lm_SEX_ITG)
tab_model(lm_SEX_ITG, show.std = TRUE)

lm_SEX_MTG <- lm(data = WRAP_cross_sectional, bl_MTG ~ SEX +bl_tau_age)
summary(lm_SEX_MTG)
tab_model(lm_SEX_MTG, show.std = TRUE)
## PCA on tau ----
WRAP <- WRAP %>% ungroup()
tau.roi <- dplyr::select(WRAP, "bl_Amyg","bl_Para", "bl_ITG", "bl_EC", "bl_Fus", "bl_MTG")

#calculate principal components
results <- prcomp(tau.roi, scale = TRUE)

#display principal components
results$rotation
#view summary of the model
summary(results)

summ <- summary(results)

summ$importance
summ$importance[2,]

eigenvalues <- results$sdev^2
summ <- summary(results)

# Create a matrix with the values
summary_matrix <- rbind(
  Standard_Deviation = summ$importance[1, ],
  Proportion_Variance = summ$importance[2, ],
  Cumulative_Proportion = summ$importance[3, ],
  Eigenvalue = eigenvalues
)

# Round for readability
summary_matrix <- round(summary_matrix, 2)

# Convert to data frame for export
summary_df <- as.data.frame(summary_matrix)
## PC 1 captures 86% of variance

## Model 1: sex x regional tau-PET x time ----
WRAP$SEX <- relevel(WRAP$SEX, ref = "M")

lme_amyg <- lme(PACC~ SEX*bl_Amyg*pacc_time +
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP,
                
)
summary(lme_amyg)
tab_model(lme_amyg, show.std = TRUE)

lme_Para <- lme(PACC~ SEX*bl_Para*pacc_time +
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_Para)
tab_model(lme_Para, show.std = TRUE)

lme_EC <- lme(PACC~ SEX*bl_EC*pacc_time + 
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_EC)
tab_model(lme_EC, show.std = TRUE)

lme_Fus <- lme(PACC~ SEX*bl_Fus*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
               
)
summary(lme_Fus)
tab_model(lme_Fus, show.std = TRUE)

lme_ITG <- lme(PACC~ SEX*bl_ITG*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_ITG)
tab_model(lme_ITG, show.std = TRUE)


lme_MTG <- lme(PACC~ SEX*bl_MTG*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_MTG)
tab_model(lme_MTG, show.std = TRUE)

## Model 1: Visualization ----
# amygdala
outfile <- file.path(path, "model1_Amyg.png")
png(outfile, width = 17, height =7, units = "cm", res = 500)
plot_model(lme_amyg, type = "pred", terms = c("pacc_time [all]", "SEX", "bl_Amyg")) +
  ylim(-3, 1.5) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 17),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),  
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Amygdala tau", x = "PACC Time from baseline tau-PET (years)", y = "PACC score", color = "SEX")
dev.off() 

# parahippocampal
outfile <- file.path(path, "model1_Para.png")
png(outfile, width = 17, height =7, units = "cm", res = 500)
plot_model(lme_Para, type = "pred", terms = c("pacc_time [all]", "SEX", "bl_Para")) +
  ylim(-3, 1.5) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 17),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),  
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Parahippocampal tau", x = "PACC Time from baseline tau-PET (years)", y = "PACC score", color = "SEX")
dev.off()

# entorhinal
outfile <- file.path(path, "model1_EC.png")
png(outfile, width = 17, height =7, units = "cm", res = 500)
plot_model(lme_EC, type = "pred", terms = c("pacc_time [all]", "SEX", "bl_EC")) +
  ylim(-3, 1.5) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 17),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),  
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Entorhinal tau", x = "PACC Time from baseline tau-PET (years)", y = "PACC score", color = "SEX")
dev.off()

# fusiform
outfile <- file.path(path, "model1_Fus.png")
png(outfile, width = 17, height =7, units = "cm", res = 500)
plot_model(lme_Fus, type = "pred", terms = c("pacc_time [all]", "SEX", "bl_Fus")) +
  ylim(-2.5, 1.5) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 17),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),  
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Fusiform tau", x = "PACC Time from baseline tau-PET (years)", y = "PACC score", color = "SEX")
dev.off()

# ITG
outfile <- file.path(path, "model1_ITG.png")
png(outfile, width = 17, height =7, units = "cm", res = 500)
plot_model(lme_ITG, type = "pred", terms = c("pacc_time [all]", "SEX", "bl_ITG")) +
  ylim(-2.5, 1.5) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 17),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),  
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Inferior temporal tau", x = "PACC Time from baseline tau-PET (years)", y = "PACC score", color = "SEX")
dev.off()

# MTG
outfile <- file.path(path, "model1_MTG.png")
png(outfile, width = 17, height =7, units = "cm", res = 500)
plot_model(lme_MTG, type = "pred", terms = c("pacc_time [all]", "SEX", "bl_MTG")) +
  ylim(-2.5, 1.5) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 17),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),  
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Medial temporal tau", x = "PACC Time from baseline tau-PET (years)", y = "PACC score", color = "SEX")
dev.off()

## Model 1: Sensitivity analysis ----
# n = 424
lme_amyg_apoe <- lme(PACC~ SEX*bl_Amyg*pacc_time + 
                       bl_tau_age*pacc_time +
                       YrsEd*pacc_time + APOE*pacc_time, 
                     na.action = na.exclude, 
                     random = (~pacc_time | ID), 
                     data =WRAP, 
                     control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                     
)
summary(lme_amyg_apoe)
tab_model(lme_amyg_apoe, show.std = TRUE)

# n = 444
lme_amyg_amyloid <- lme(PACC~ SEX*bl_Amyg*pacc_time + 
                          bl_tau_age*pacc_time +
                          YrsEd*pacc_time + centiloid*pacc_time, 
                        na.action = na.exclude, 
                        random = (~pacc_time | ID), 
                        data =WRAP, 
                        control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                        
)
summary(lme_amyg_amyloid)
tab_model(lme_amyg_amyloid, show.std = TRUE)

lme_Para_apoe <- lme(PACC~ SEX*bl_Para*pacc_time + 
                       bl_tau_age*pacc_time +
                       YrsEd*pacc_time + APOE*pacc_time, 
                     na.action = na.exclude, 
                     random = (~pacc_time | ID), 
                     data =WRAP, 
                     control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                     
)
summary(lme_Para_apoe)
tab_model(lme_Para_apoe, show.std = TRUE)

lme_Para_amyloid <- lme(PACC~ SEX*bl_Para*pacc_time + 
                          bl_tau_age*pacc_time +
                          YrsEd*pacc_time + centiloid*pacc_time, 
                        na.action = na.exclude, 
                        random = (~pacc_time | ID), 
                        data =WRAP, 
                        control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                        
)
summary(lme_Para_amyloid)
tab_model(lme_Para_amyloid, show.std = TRUE)

lme_EC_apoe <- lme(PACC~ SEX*bl_EC*pacc_time + 
                     bl_tau_age*pacc_time +
                     YrsEd*pacc_time + APOE*pacc_time, 
                   na.action = na.exclude, 
                   random = (~pacc_time | ID), 
                   data =WRAP, 
                   control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                   
)
summary(lme_EC_apoe)
tab_model(lme_EC_apoe, show.std = TRUE)

lme_EC_amyloid <- lme(PACC~ SEX*bl_EC*pacc_time + 
                        bl_tau_age*pacc_time +
                        YrsEd*pacc_time + centiloid*pacc_time, 
                      na.action = na.exclude, 
                      random = (~pacc_time | ID), 
                      data =WRAP, 
                      control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                      
)
summary(lme_EC_amyloid)
tab_model(lme_EC_amyloid, show.std = TRUE)

lme_Fus_apoe <- lme(PACC~ SEX*bl_Fus*pacc_time + 
                      bl_tau_age*pacc_time +
                      YrsEd*pacc_time + APOE*pacc_time, 
                    na.action = na.exclude, 
                    random = (~pacc_time | ID), 
                    data =WRAP, 
                    control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                    
)
summary(lme_Fus_apoe)
tab_model(lme_Fus_apoe, show.std = TRUE)

lme_Fus_amyloid <- lme(PACC~ SEX*bl_Fus*pacc_time + 
                         bl_tau_age*pacc_time +
                         YrsEd*pacc_time + centiloid*pacc_time, 
                       na.action = na.exclude, 
                       random = (~pacc_time | ID), 
                       data =WRAP, 
                       control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                       
)
summary(lme_Fus_amyloid)
tab_model(lme_Fus_amyloid, show.std = TRUE)

lme_ITG_apoe <- lme(PACC~ SEX*bl_ITG*pacc_time + 
                      bl_tau_age*pacc_time +
                      YrsEd*pacc_time + APOE*pacc_time, 
                    na.action = na.exclude, 
                    random = (~pacc_time | ID), 
                    data =WRAP, 
                    control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                    
)
summary(lme_ITG_apoe)
tab_model(lme_ITG_apoe, show.std = TRUE)

lme_ITG_amyloid <- lme(PACC~ SEX*bl_ITG*pacc_time + 
                         bl_tau_age*pacc_time +
                         YrsEd*pacc_time + centiloid*pacc_time, 
                       na.action = na.exclude, 
                       random = (~pacc_time | ID), 
                       data =WRAP, 
                       control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                       
)
summary(lme_ITG_amyloid)
tab_model(lme_ITG_amyloid, show.std = TRUE)

lme_MTG_apoe <- lme(PACC~ SEX*bl_MTG*pacc_time + 
                      bl_tau_age*pacc_time +
                      YrsEd*pacc_time + APOE*pacc_time, 
                    na.action = na.exclude, 
                    random = (~pacc_time | ID), 
                    data =WRAP, 
                    control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                    
)
summary(lme_MTG_apoe)
tab_model(lme_MTG_apoe, show.std = TRUE)

lme_MTG_amyloid <- lme(PACC~ SEX*bl_MTG*pacc_time + 
                         bl_tau_age*pacc_time +
                         YrsEd*pacc_time + centiloid*pacc_time, 
                       na.action = na.exclude, 
                       random = (~pacc_time | ID), 
                       data =WRAP, 
                       control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                       
)
summary(lme_MTG_amyloid)
tab_model(lme_MTG_amyloid, show.std = TRUE)
## Model 1: Floodlight analysis ----
WRAP_PACC <- WRAP
pure_pacc <- lme(PACC ~ pacc_time, random = list(ID = pdDiag(~ pacc_time)), data = WRAP_PACC)
lme_slope <- coef(pure_pacc)$pacc_time
WRAP_cross_sectional$lme_pacc_slope <- lme_slope
WRAP_cross_sectional$SEX <- relevel(WRAP_cross_sectional$SEX, ref = "M")
WRAP_cross_sectional$SEX<-as.numeric(WRAP_cross_sectional$SEX)

## Johnson Neyman
mod<-lm(lme_pacc_slope ~ SEX*bl_Amyg, data=WRAP_cross_sectional, na.action = na.exclude); summary(mod)
jn_results_amyg <- johnson_neyman(model = mod, pred =SEX,  modx = bl_Amyg)
outfile <- file.path(path, "jn_amyg.png")
png(outfile, width = 16, height =14, units = "cm", res = 500)
jn_results_amyg$plot +
  labs(x = "Baseline Amygdala SUVR", 
       y = "Slope of SEX") +
  theme(legend.position = "bottom") 
dev.off() 

mod<-lm(lme_pacc_slope ~ SEX*bl_Para, data=WRAP_cross_sectional, na.action = na.exclude); summary(mod)
jn_results_para <- johnson_neyman(model = mod, pred =SEX,  modx = bl_Para)
outfile <- file.path(path, "jn_para.png")
png(outfile, width = 16, height =14, units = "cm", res = 500)
jn_results_para$plot +
  labs(x = "Baseline Parahippocampal SUVR", 
       y = "Slope of SEX") +
  theme(legend.position = "bottom") 
dev.off() 

## Model 1: sex x regional tau-PET x time (prospective only) ----
WRAP_prospective <- WRAP %>%
  filter(pacc_time >= 0)
WRAP_prospective$SEX <- relevel(WRAP_prospective$SEX, ref = "M")

lme_amyg <- lme(PACC~ SEX*bl_Amyg*pacc_time + 
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time , 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP_prospective, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_amyg)
tab_model(lme_amyg, show.std = TRUE)

lme_Para <- lme(PACC~ SEX*bl_Para*pacc_time + 
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP_prospective, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_Para)
tab_model(lme_Para, show.std = TRUE)

lme_EC <- lme(PACC~ SEX*bl_EC*pacc_time + 
                bl_tau_age*pacc_time +
                YrsEd*pacc_time, 
              na.action = na.exclude, 
              random = (~pacc_time | ID), 
              data =WRAP_prospective, 
              control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
              
)
summary(lme_EC)
tab_model(lme_EC, show.std = TRUE)

lme_Fus <- lme(PACC~ SEX*bl_Fus*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP_prospective, 
               control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
               
)
summary(lme_Fus)
tab_model(lme_Fus, show.std = TRUE)

lme_ITG <- lme(PACC~ SEX*bl_ITG*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time , 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP_prospective, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_ITG)
tab_model(lme_ITG, show.std = TRUE)


lme_MTG <- lme(PACC~ SEX*bl_MTG*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time , 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP_prospective, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_MTG)
tab_model(lme_MTG, show.std = TRUE)


## Model 1A: sex x regional tau-PET x time + sex x Aβ-CL x time ----
lme_amyg <- lme(PACC~ SEX*bl_Amyg*pacc_time + SEX*centiloid*pacc_time +
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_amyg)
tab_model(lme_amyg, show.std = TRUE)

lme_Para <- lme(PACC~ SEX*bl_Para*pacc_time + SEX*centiloid*pacc_time +
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_Para)
tab_model(lme_Para, show.std = TRUE)

lme_EC <- lme(PACC~ SEX*bl_EC*pacc_time + SEX*centiloid*pacc_time +
                bl_tau_age*pacc_time +
                YrsEd*pacc_time, 
              na.action = na.exclude, 
              random = (~pacc_time | ID), 
              data =WRAP, 
              control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
              
)
summary(lme_EC)
tab_model(lme_EC, show.std = TRUE)

lme_Fus <- lme(PACC~ SEX*bl_Fus*pacc_time + SEX*centiloid*pacc_time +
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
               
)
summary(lme_Fus)
tab_model(lme_Fus, show.std = TRUE)

lme_ITG <- lme(PACC~ SEX*bl_ITG*pacc_time + SEX*centiloid*pacc_time +
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_ITG)
tab_model(lme_ITG, show.std = TRUE)


lme_MTG <- lme(PACC~ SEX*bl_MTG*pacc_time + SEX*centiloid*pacc_time +
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_MTG)
tab_model(lme_MTG, show.std = TRUE)

## Model 1B: sex x regional tau-PET x time + sex x APOE x time ----
lme_amyg <- lme(PACC~ SEX*bl_Amyg*pacc_time + SEX*APOE*pacc_time +
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_amyg)
tab_model(lme_amyg, show.std = TRUE)

lme_Para <- lme(PACC~ SEX*bl_Para*pacc_time + SEX*APOE*pacc_time +
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_Para)
tab_model(lme_Para, show.std = TRUE)

lme_EC <- lme(PACC~ SEX*bl_EC*pacc_time + SEX*APOE*pacc_time +
                bl_tau_age*pacc_time +
                YrsEd*pacc_time, 
              na.action = na.exclude, 
              random = (~pacc_time | ID), 
              data =WRAP, 
              control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
              
)
summary(lme_EC)
tab_model(lme_EC, show.std = TRUE)

lme_Fus <- lme(PACC~ SEX*bl_Fus*pacc_time + SEX*APOE*pacc_time +
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
               
)
summary(lme_Fus)
tab_model(lme_Fus, show.std = TRUE)

lme_ITG <- lme(PACC~ SEX*bl_ITG*pacc_time + SEX*APOE*pacc_time +
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_ITG)
tab_model(lme_ITG, show.std = TRUE)


lme_MTG <- lme(PACC~ SEX*bl_MTG*pacc_time + SEX*APOE*pacc_time +
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_MTG)
tab_model(lme_MTG, show.std = TRUE)

## Model 2: sex x regional tau-PET x Aβ-CL x time ----
lme_amyg <- lme(PACC~ SEX*bl_Amyg*centiloid*pacc_time + 
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time , 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_amyg)
tab_model(lme_amyg, show.std = TRUE)

lme_Para <- lme(PACC~ SEX*bl_Para*centiloid*pacc_time + 
                  bl_tau_age*pacc_time +
                  YrsEd*pacc_time, 
                na.action = na.exclude, 
                random = (~pacc_time | ID), 
                data =WRAP, 
                control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
                
)
summary(lme_Para)
tab_model(lme_Para, show.std = TRUE)

lme_EC <- lme(PACC~ SEX*bl_EC*centiloid*pacc_time + 
                bl_tau_age*pacc_time +
                YrsEd*pacc_time, 
              na.action = na.exclude, 
              random = (~pacc_time | ID), 
              data =WRAP, 
              control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
              
)
summary(lme_EC)
tab_model(lme_EC, show.std = TRUE)

lme_Fus <- lme(PACC~ SEX*bl_Fus*centiloid*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time, 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim', msMaxIter = 80, msMaxEval = 80),
               
)
summary(lme_Fus)
tab_model(lme_Fus, show.std = TRUE)

lme_ITG <- lme(PACC~ SEX*bl_ITG*centiloid*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time , 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_ITG)
tab_model(lme_ITG, show.std = TRUE)


lme_MTG <- lme(PACC~ SEX*bl_MTG*centiloid*pacc_time + 
                 bl_tau_age*pacc_time +
                 YrsEd*pacc_time , 
               na.action = na.exclude, 
               random = (~pacc_time | ID), 
               data =WRAP, 
               control = lmeControl(opt = 'optim'),
               
)
summary(lme_MTG)
tab_model(lme_MTG, show.std = TRUE)
