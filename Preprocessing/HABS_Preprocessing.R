# Load librarys ----
library(dplyr)
library(nlme) 
library(sjPlot)
library(ggplot2)
library(emmeans)
library(clubSandwich)
library(gt)
library(gtsummary)
library(flextable)
library(flowchart)

load(file = "~/Desktop/sexcog/Raw_Data/HABS/HABS_DataCentral.RData")

## construct flowchart ----
data_bl <- data %>%
  arrange(SubjIDshort, NP_SessionDate) %>%
  group_by(SubjIDshort) %>%
  slice(1)
habs_fc <- data_bl %>%
  as_fc(label = "HABS Participants with Cognitive Data Available", text_pattern = "{label} \n N = {n}") %>%
  fc_filter(!is.na(TAU_Session_ID), label = "Eligible Participants", text_pattern = "N = {n} {label}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Without Tau-PET Data") %>%
  fc_filter(HABS_DX == "CN", label = "Eligible Participants", text_pattern = "N = {n} {label}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Missing Diagnosis At Baseline") %>%
  fc_filter(!is.na(YrsEd), label = "Eligible HABS Participants", text_pattern = "{label}\nN = {n}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Missing Years of Education") %>%
  fc_draw()


# set up data ----
## Remove unused columns and rename variables
habs_data <- data %>%
  select(-c("HollinScore", contains("j7"), "rn")) %>%
  dplyr::rename(APOE = E4_Status,
                ID = SubjIDshort,
                pib_group = PIB_FS_DVR_Group,
                bl_tau_age = TAU_Age,
                PACC = NP_PACC_PACC96,
                bl_Amyg = TAU_FS_SUVR_Amygdala_bh,
                bl_Para = TAU_FS_SUVR_parahippocampal_bh,
                bl_EC = TAU_FS_SUVR_entorhinal_bh,
                bl_Fus = TAU_FS_SUVR_fusiform_bh,
                bl_ITG = TAU_FS_SUVR_inferiortemporal_bh,
                bl_MTG = TAU_FS_SUVR_middletemporal_bh)

# make factors
habs_data$SEX<-as.factor(habs_data$SEX)
habs_data$APOE<-as.factor(habs_data$APOE)
habs_data$HABS_DX<-as.factor(habs_data$HABS_DX)
habs_data$pib_group<-as.factor(habs_data$pib_group)

## Filter rows w/o Tau
habs_data_filtered <- habs_data[!is.na(habs_data$TAU_Session_ID), ]

# create centiloid 
habs_data_filtered <- habs_data_filtered %>%
  mutate(centiloid = 143.06 * PIB_FS_DVR_FLR - 145.60)

# create baseline PACC age
habs_data_filtered <- habs_data_filtered %>%
  group_by(ID) %>%
  mutate(bl_pacc_age = min(NP_Age))

# Filter subjects with baseline PACC diagnosis = CN
baseline_data <- habs_data_filtered %>%
  arrange(ID, NP_SessionDate) %>%
  group_by(ID) %>%
  slice(1) %>%
  filter(HABS_DX == "CN")
baseline_data$ID <- as.character(baseline_data$ID)
cn_subjects <- baseline_data$ID

habs_data_filtered <- habs_data_filtered %>%
  filter(ID %in% cn_subjects)

habs_data_filtered <- habs_data_filtered %>%
  group_by(ID) %>% 
  mutate(first_NP_SessionDate = min(NP_SessionDate, na.rm = TRUE)) %>%
  mutate(pacc_time = as.numeric((NP_SessionDate - TAU_SessionDate) / 365.25),
         pacc_time_from_pacc = as.numeric((NP_SessionDate - first_NP_SessionDate) / 365.25))

## Filter NAs
habs_data_filtered %>%
  ungroup()

count_missing <- function(data, variables) {
  sapply(variables, function(var) sum(is.na(data[[var]])))
}

variables_to_check <- c(
  "PACC", "pacc_time", "bl_tau_age", "bl_Para", "bl_Amyg", "bl_EC", "bl_Fus",
  "bl_ITG", "bl_MTG", "centiloid", "YrsEd", "SEX", "APOE"
)

missing_counts <- count_missing(habs_data_filtered, variables_to_check)
missing_counts # 25 obs w/ APOE missing, 4 obs w/ YrsEd missing, 54 obs w/ centiloid missing


HABS <- habs_data_filtered %>% 
  filter(
    !is.na(YrsEd)
  )

length(unique(HABS$ID)) # n = 279

## create pacc time point
HABS$PACC_tp <- ave(as.character(HABS$ID), HABS$ID, FUN = length)
HABS$PACC_tp <- as.numeric(HABS$PACC_tp)

HABS <- HABS %>%
  mutate(One_PACC_tp = case_when(
    PACC_tp == 1 ~ T,  
    PACC_tp > 1 ~ F)
  )

HABS <- HABS %>%
  group_by(ID) %>%
  mutate(Avg_PACC_time_from_pacc = sum(pacc_time_from_pacc, na.rm = TRUE) / n(),
         Avg_PACC_time_from_TAU = sum(pacc_time, na.rm = TRUE)/ n(),
         Avg_pro_PACC_time_from_TAU = sum(pacc_time[pacc_time >= 0], na.rm = TRUE) / 
           sum(pacc_time >= 0, na.rm = TRUE),
         Avg_retro_PACC_time_from_TAU = sum(pacc_time[pacc_time < 0], na.rm = TRUE) / 
           sum(pacc_time < 0, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(prospective_vs_retro = case_when(
    pacc_time < 0 ~ "Retro N",
    pacc_time >= 0 ~ "Pros N"
  ))
HABS <- HABS %>%
  group_by(ID) %>%
  mutate(
    tau_pacc_age_diff = abs(NP_Age - bl_tau_age),
    closest_PACC = PACC[which.min(tau_pacc_age_diff)]
  ) %>%
  ungroup() %>%
  select(-c(contains("bh")))

save(HABS, file = paste0("~/Desktop/sexcog/Preprocessing/HABS", Sys.Date(), ".RData"))

