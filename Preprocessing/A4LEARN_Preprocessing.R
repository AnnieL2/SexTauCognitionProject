library(stringr)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(readxl)
library(flowchart)

adqs <- read.csv("~/Desktop/sexcog/Raw_Data/A4LEARN/ADQS.csv", na.strings = c("", " ", "NaN", "NA", NA))
amyloid <- read.csv("~/Desktop/sexcog/Raw_Data/A4LEARN/imaging_SUVR_amyloid.csv", na.strings = c("", " ", "NaN", "NA", NA))
tau <- read.csv("~/Desktop/sexcog/Raw_Data/A4LEARN/imaging_SUVR_tau.csv", na.strings = c("", " ", "NaN", "NA", NA))
subj <- read.csv("~/Desktop/sexcog/Raw_Data/A4LEARN/SUBJINFO.csv", na.strings = c("", " ", "NaN", "NA", NA))
visits <- read.csv("~/Desktop/sexcog/Raw_Data/A4LEARN/SV.csv", na.strings = c("", " ", "NaN", "NA", NA))
visit <- read.table(file = "~/Desktop/sexcog/Raw_Data/A4LEARN/visitinfo.csv", header = T, sep = ",", na.strings = c("NaN", NaN))
dose <- read_excel(path = "~/Desktop/sexcog/Raw_Data/A4LEARN/A4_DOSE_blind.xlsx")

## construct flowchart ----
eligible <- adqs %>%
  filter(QSTESTCD == "PACC") %>% 
  dplyr::rename(PACC = QSSTRESN) %>% 
  select(SUBSTUDY, BID, TX, AGEYR, SEX, EDCCNTU, RACE, ETHNIC, 
         APOE = APOEGNPRSNFLG, APOEGN, BMIBL,
         pacc_VISCODE = VISITCD, VISIT, EPOCH, QSVERSION, ASEQNCS, PACC,
         pacc_consent_time = QSDTC_DAYS_CONSENT) %>%
  group_by(BID) %>%
  arrange(BID, pacc_consent_time) %>%
  slice(1)

visits_fc <- select(visits, BID, VISCODE = VISITCD, SVSTDTC_DAYS_CONSENT)
visits_tau <- dplyr::rename(visits_fc, tau_VISCODE = VISCODE)

tau_roi <- c("Amygdala", "parahippocampal", "entorhinal", "fusiform", "inferiortemporal", "middletemporal")
tau_vars <- unique(filter(tau, str_detect(brain_region, paste(c(tau_roi), collapse = '|')))$brain_region)

tau_fc <- tau %>%
  filter(str_detect(brain_region, paste(c(tau_roi), collapse = '|'))) %>%
  spread(brain_region, suvr_crus) %>%
  group_by(BID, scan_number) %>%
  fill(tau_vars, .direction = "updown") %>%
  select(BID, tau_VISCODE = VISCODE, tau_scan_number = scan_number, scan_date_DAYS_CONSENT,
         contains(tau_roi)) %>%
  distinct() %>%
  left_join(visits_tau, by = c("BID", "tau_VISCODE")) %>%
  ungroup() %>%
  group_by(BID) %>%
  arrange(scan_date_DAYS_CONSENT, SVSTDTC_DAYS_CONSENT) %>%
  slice(1)

eligible <- eligible %>%
  left_join(tau_fc, by = "BID")

eligible$SUBSTUDY <- as.factor(eligible$SUBSTUDY)

A4LEARN_fc <- eligible %>%
  as_fc(label = "A4/LEARN Participants with Cognitive Data Available", text_pattern = "{label} \n N = {n}") %>%
  fc_filter(SUBSTUDY == "LEARN" | SUBSTUDY == "A4", label = "Eligible Participants", text_pattern = "N = {n} {label}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Screen Failed") %>%
  fc_filter(is.na(TX) | TX == "Placebo", label = "Eligible Participants", text_pattern = "N = {n} {label}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Assigned to Solanezumab") %>%
  fc_filter(!is.na(Amygdala_L), label = "Eligible A4/LEARN Participants", text_pattern = "{label} \n N = {n} (A4:196, LEARN:54)",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Without Tau-PET Data") %>%
  fc_draw()


# Visits ------------------------------------------------------------------

visits <- select(visits, BID, VISCODE = VISITCD, SVSTDTC_DAYS_CONSENT)


# Rescreens ---------------------------------------------------------------

rescreens <- subj %>% 
  filter(!is.na(PREVBID)) %>% 
  dplyr::rename(BID1 = PREVBID, BID2 = BID) %>% 
  select(BID1, BID2)


# PACC --------------------------------------------------------------------

data <- adqs %>%
  filter(QSTESTCD == "PACC") %>% #if you want cog tests other than pacc remove this line and convert data from long to wide format
  dplyr::rename(PACC = QSSTRESN) %>% 
  select(SUBSTUDY, BID, TX, AGEYR, SEX, EDCCNTU, RACE, ETHNIC, 
         APOE = APOEGNPRSNFLG, APOEGN, BMIBL,
         pacc_VISCODE = VISITCD, VISIT, EPOCH, QSVERSION, ASEQNCS, PACC,
         pacc_consent_time = QSDTC_DAYS_CONSENT) %>%
  group_by(BID) %>%
  dplyr::mutate(
    pacc_consent_time = pacc_consent_time / 365.25,
    pacc_time_old = pacc_consent_time - pacc_consent_time[1],
    pacc_age = AGEYR + pacc_time_old,
    bl_pacc_age = min(pacc_age, na.rm = TRUE),
    SEX = case_when(
      SEX == 1 ~ "Female",
      SEX == 2 ~ "Male", 
      TRUE ~ NA),
    ETHNIC = case_when(
      ETHNIC == 50 ~ "Hispanic or Latino",
      ETHNIC == 56 ~ "Not Hispanic or Latino",
      TRUE ~ NA),
    RACE = case_when(
      RACE == 1 ~ "White",
      RACE == 2 ~ "Black or African American",
      RACE == 58 ~ "Asian",
      RACE == 79 ~ "Native Hawaiian or Other Pacific Islander",
      RACE == 84 ~ "American Indian or Alaskan Native",
      RACE == 100 ~ "More than one race",
      TRUE ~ NA),
    APOE = case_when(
      APOE == 1 ~ "e4+",
      APOE == 0 ~ "e4-",
      TRUE ~ NA)
  ) %>%
  arrange(BID, pacc_consent_time) %>%
  select(-pacc_time_old)


# Amyloid -----------------------------------------------------------------

baseline_amyloid <- amyloid %>%
  filter(brain_region == "Composite_Summary") %>%
  left_join(visits, by = c("BID", "VISCODE")) %>%
  group_by(BID) %>%
  dplyr::mutate(ab_consent_time = ifelse(is.na(scan_date_DAYS_CONSENT), SVSTDTC_DAYS_CONSENT/365.25, scan_date_DAYS_CONSENT/365.25),
                overall_score = case_when(
                  SUBSTUDY == "A4" ~ "Ab+",
                  TRUE ~ "Ab-"
                )) %>%
  select(BID, ab_VISCODE = VISCODE, ab_consent_time, suvr_cer, overall_score) %>%
  arrange(BID, ab_consent_time) %>%
  mutate(
    centiloid = 183.07 * suvr_cer - 177.26
  ) %>%
  slice(1)


# Tau ---------------------------------------------------------------------

tau_roi <- c("Amygdala", "parahippocampal", "entorhinal", "fusiform", "inferiortemporal", "middletemporal")

tau_vars <- unique(filter(tau, str_detect(brain_region, paste(c(tau_roi), collapse = '|')))$brain_region)


visits_tau <- dplyr::rename(visits, tau_VISCODE = VISCODE)

tau_new <- tau %>%
  filter(str_detect(brain_region, paste(c(tau_roi), collapse = '|'))) %>%
  spread(brain_region, suvr_crus) %>%
  group_by(BID, scan_number) %>%
  fill(tau_vars, .direction = "updown") %>%
  select(BID, tau_VISCODE = VISCODE, tau_scan_number = scan_number, scan_date_DAYS_CONSENT,
         contains(tau_roi)) %>%
  distinct() %>%
  left_join(visits_tau, by = c("BID", "tau_VISCODE")) %>%
  group_by(BID) %>%
  dplyr::mutate(tau_consent_time = ifelse(is.na(scan_date_DAYS_CONSENT), SVSTDTC_DAYS_CONSENT/365.25, scan_date_DAYS_CONSENT/365.25)) %>%
  arrange(BID, tau_consent_time) 

bilateral_tau_calc <- function(data, data_name, roi){
  
  vector_names <- c(); left <- c(); right <- c()
  
  for (i in 1:length(tau_roi)){
    vector_names[i] <- paste0(tau_roi[i], "_bh")
    left[i] <- paste0(tau_roi[i], "_lh_VOI")
    right[i] <- paste0(tau_roi[i], "_rh_VOI")
    
    data[, vector_names[i]] <- (data[,left[i]] + data[,right[i]]) / 2
  }
   
  assign(data_name, data, .GlobalEnv)
}

bilateral_tau_calc(tau_new, "tau_new", tau_roi)

baseline_tau <- tau_new %>% 
  select(-c(contains("VOI"), SVSTDTC_DAYS_CONSENT, scan_date_DAYS_CONSENT)) %>%
  group_by(BID) %>%
  arrange(tau_consent_time) %>%
  slice(1) %>%
  select(-c(Amygdala_L, Amygdala_R))


# Merge -------------------------------------------------------------------

a4_data <- data %>%
  right_join(baseline_tau, by = "BID") %>%
  left_join(baseline_amyloid, by = "BID") %>%
  mutate(Cohort = case_when(
    SUBSTUDY == "SF" | SUBSTUDY == "LEARN" ~ "LEARN/SF",
    SUBSTUDY == "A4" & TX == "Placebo" ~ "A4 Placebo",
    SUBSTUDY == "A4" & TX == "Solanezumab" ~ "A4 Treated",
    TRUE ~ NA),
    ab_age = AGEYR + ab_consent_time) %>%
  select(SUBSTUDY, Cohort, BID, TX:BMIBL, VISIT, EPOCH, 
         pacc_VISCODE, QSVERSION:bl_pacc_age, tau_VISCODE:middletemporal_bh, ab_consent_time,
         overall_score, centiloid, Cohort) %>%
  group_by(BID) %>%
  mutate(pacc_time = pacc_consent_time - tau_consent_time[1],
         pacc_time_from_pacc = pacc_consent_time - pacc_consent_time[1])

length(unique(a4_data$BID)) # n = 449

# Create baseline tau age ------------------------------------------
a4_data <- a4_data %>%
  arrange(BID, pacc_consent_time) %>%
  group_by(BID) %>%
  dplyr::mutate(bl_tau_age = AGEYR + tau_consent_time)

# Calculates Cumulative Dose (Script from HK 1st Feb 2024) ------------------------------------------
visit <- visit %>%
  rename(subj = USUBJID,
         VISCODE = VISITCD) %>%
  dplyr::select(-c(VISIT, AVISIT))

dose <- rename(dose, subj = BID)

a4_data <- rename(a4_data, subj = BID)
a4_data <- rename(a4_data, VISCODE = pacc_VISCODE)
a4_data <- rename(a4_data, time_pacc = pacc_time)
a4_data <- rename(a4_data, data_epoch = EPOCH)


trt <- a4_data %>% 
  dplyr::select(subj, VISCODE, time_pacc, TX) %>%
  arrange(subj, time_pacc)

dose_trt <- merge(dose, trt, by = c("subj","VISCODE"), all.y = TRUE)
dose_trt <- merge(dose_trt, visit, by = c("subj", "VISCODE"), all.x = TRUE)

dose_trt <- dose_trt %>%
  mutate(Dose = case_when(
    TX == "Placebo" & EPOCH != "OPEN LABEL TREATMENT" ~ 0,
    EPOCH == "SCREENING" ~ 0,
    is.na(BLINDDOSE) ~ 0,
    TRUE ~ BLINDDOSE
  )) %>%
  group_by(subj) %>%
  arrange(subj, time_pacc) %>%
  mutate(Cumulative_Dose = cumsum(Dose)) %>%
  mutate(Max_Cumulative_Dose = max(Cumulative_Dose, na.rm = TRUE)) %>%
  dplyr::select(subj, VISCODE, EPOCH, Dose, Cumulative_Dose, Max_Cumulative_Dose)

new_a4_data <- merge(a4_data, dose_trt, by = c("subj", "VISCODE"), all.x = TRUE)

new_a4_data <- arrange(new_a4_data, subj, VISCODE)

## plots dose per substudy
ggplot(subset(new_a4_data, !is.na(time_pacc)), aes(x = time_pacc, y = Cumulative_Dose)) + 
  geom_line(aes(group = subj)) + 
  geom_point() + 
  facet_wrap(~TX)

# Scale cumulative dose - email from Mike D 9th Feb 2024
new_a4_data <- new_a4_data %>%
  mutate(Cumulative_Dose_Scaled = Cumulative_Dose / 10000)

# Comparing EPOCH from adqs and EPOCH from "~/Downloads/visitinfo.csv" (from Rory's folder)
# Using EPOCH from adqs instead (looks like it's more updated)
mismatched_rows <- new_a4_data %>%
  filter(is.na(data_epoch) != is.na(EPOCH) | data_epoch != EPOCH) %>%
  select(subj, SUBSTUDY, data_epoch, EPOCH, data_epoch, VISCODE)
new_a4_data <- new_a4_data %>%
  select(-EPOCH)

# Rename variables
new_a4_data <- rename(new_a4_data, BID = subj)
new_a4_data <- rename(new_a4_data, pacc_VISCODE = VISCODE)
new_a4_data <- rename(new_a4_data, pacc_time = time_pacc)
new_a4_data <- rename(new_a4_data, EPOCH = data_epoch)

length(unique(new_a4_data$BID)) # n = 449

## clean variables ----
new_a4_data$SEX<-as.factor(new_a4_data$SEX)
new_a4_data$TX<-as.factor(new_a4_data$TX)
new_a4_data$SUBSTUDY<-as.factor(new_a4_data$SUBSTUDY)
new_a4_data$RACE<-as.factor(new_a4_data$RACE)
new_a4_data$ETHNIC<-as.factor(new_a4_data$ETHNIC)
new_a4_data$APOE<-as.factor(new_a4_data$APOE)
new_a4_data$QSVERSION <- as.factor(new_a4_data$QSVERSION)
new_a4_data$overall_score <- as.factor(new_a4_data$overall_score)

## Filter out SF (no long PACC data)
new_a4_data <- new_a4_data[new_a4_data$SUBSTUDY != "SF", ]
new_a4_data$SUBSTUDY <- droplevels(new_a4_data$SUBSTUDY)

length(unique(new_a4_data$BID)) # n = 443

## Recreate SUBSTUDY variable including LEARN
new_a4_data %>%
  group_by(SUBSTUDY, TX) %>%
  dplyr::summarize(TX_count = n(), .groups = "drop")

new_a4_data <- new_a4_data %>%
  mutate(
    SUBSTUDY = case_when(
      is.na(TX) & SUBSTUDY == "LEARN" ~ "LEARN",
      TX == "Placebo" ~ "A4_PLACEBO",
      TX == "Solanezumab" ~ "A4_TREATMENT",
      TRUE ~ NA_character_  
    )
  )

new_a4_data$SUBSTUDY<-as.factor(new_a4_data$SUBSTUDY)

## Rename tau variables 
new_a4_data$bl_Amyg <- new_a4_data$Amygdala_bh
new_a4_data$bl_EC <- new_a4_data$entorhinal_bh
new_a4_data$bl_Para <- new_a4_data$parahippocampal_bh
new_a4_data$bl_Fus <- new_a4_data$fusiform_bh
new_a4_data$bl_ITG <- new_a4_data$inferiortemporal_bh
new_a4_data$bl_MTG <- new_a4_data$middletemporal_bh
hist(new_a4_data$bl_Amyg)
hist(new_a4_data$bl_EC)
hist(new_a4_data$bl_Para)
hist(new_a4_data$bl_Fus)
hist(new_a4_data$bl_ITG)
hist(new_a4_data$bl_MTG)

## Filter NAs
n_PACC_missing <- sum(is.na(new_a4_data$PACC)) # 85
n_pacc_time_missing <- sum(is.na(new_a4_data$pacc_time)) # 19
n_age_missing <- sum(is.na(new_a4_data$bl_tau_age)) # 0
n_bl_Para_missing <- sum(is.na(new_a4_data$bl_Para)) # 15
n_bl_Amyg_missing <- sum(is.na(new_a4_data$bl_Amyg)) # 15
n_bl_ITG_missing <- sum(is.na(new_a4_data$bl_ITG)) # 15
n_bl_MTG_missing <- sum(is.na(new_a4_data$bl_MTG)) # 15
n_bl_EC_missing <- sum(is.na(new_a4_data$bl_EC)) # 15
n_bl_Fus_missing <- sum(is.na(new_a4_data$bl_Fus)) # 15
n_centiloid_missing <- sum(is.na(new_a4_data$centiloid)) # 0
n_EDCCNTU_missing <- sum(is.na(new_a4_data$EDCCNTU)) # 0
n_APOE_missing <- sum(is.na(new_a4_data$APOE)) # 0
n_QSVERSION_missing <- sum(is.na(new_a4_data$QSVERSION)) # 75

## drop missing data
A4_filtered <- new_a4_data %>% 
  filter(
    !is.na(PACC),
    !is.na(pacc_time),
    !is.na(bl_Para),
    !is.na(bl_Amyg),
    !is.na(bl_ITG),
    !is.na(bl_MTG),
    !is.na(bl_EC),
    !is.na(bl_Fus),
    !is.na(bl_tau_age),
    !is.na(QSVERSION)
  )

length(unique(A4_filtered$BID)) # n = 442 

# Subset for placebo and LEARN only ----
A4_LEARN <- A4_filtered %>%
  subset(SUBSTUDY %in% c("LEARN", "A4_PLACEBO"))

A4_LEARN$SUBSTUDY <- droplevels(A4_LEARN$SUBSTUDY)

length(unique(A4_LEARN$BID)) # n = 250

A4_LEARN$PACC_tp <- ave(A4_LEARN$BID, A4_LEARN$BID, FUN = length)
A4_LEARN$PACC_tp <- as.numeric(as.character(A4_LEARN$PACC_tp))

A4_LEARN <- A4_LEARN %>%
  mutate(One_PACC_tp = case_when(
    PACC_tp == 1 ~ T,  
    PACC_tp > 1 ~ F)
  )

A4_LEARN <- A4_LEARN %>%
  group_by(BID) %>%
  mutate(Avg_PACC_time_from_PACC = sum(pacc_time_from_pacc, na.rm = TRUE) / n(),
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

A4_LEARN <- A4_LEARN %>%
  group_by(BID) %>%
  mutate(
    tau_pacc_age_diff = abs(pacc_age - bl_tau_age),
    closest_PACC = PACC[which.min(tau_pacc_age_diff)]
  ) %>%
  ungroup()

A4_LEARN <- A4_LEARN %>%
  dplyr::rename(YrsEd = EDCCNTU,
         ID = BID,
         pib_group = overall_score) %>%
  mutate(SEX = case_when(
    SEX == "Female" ~ "F",  
    SEX == "Male" ~ "M"),
    pib_group = case_when(
    pib_group == "Ab-" ~ "PIB-",
    pib_group == "Ab+" ~ "PIB+")
  ) %>%
  select(-c(Amygdala_bh, parahippocampal_bh, entorhinal_bh, fusiform_bh, inferiortemporal_bh, middletemporal_bh))

A4_LEARN$SEX<-as.factor(A4_LEARN$SEX)
A4_LEARN$pib_group<-as.factor(A4_LEARN$pib_group)
save(A4_LEARN, file = paste0("~/Desktop/sexcog/Preprocessing/A4LEARN", Sys.Date(), ".RData"))
