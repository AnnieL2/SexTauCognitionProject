# Load libraries ----
library(dplyr)
library(stringr)
library(data.table)

# Read in files ----
Tau <- read.csv('~/Desktop/sexcog/Raw_Data/WRAP/Tau_updated.csv', na.strings = c("", " ", NA, "NA"))
Amyloid <- read.csv('~/Desktop/sexcog/Raw_Data/WRAP/WRAP_AB_2024.csv', na.strings = c("", " ", NA, "NA"))
subject_id <- read.csv('~/Desktop/sexcog/Raw_Data/WRAP/subject_id.csv', na.strings = c("", " ", NA, "NA"))
demo <- read.csv('~/Desktop/sexcog/Raw_Data/WRAP/demo.csv', na.strings = c("", " ", NA, "NA"))
cognition <- read.csv('~/Desktop/sexcog/Raw_Data/WRAP/cognition.csv', na.strings = c("", " ", NA, "NA"))
visits <- read.csv("~/Desktop/sexcog/Raw_Data/WRAP/link.csv", na.strings = c("", " ", "NA", NA))
diagnosis <- read.csv('~/Desktop/sexcog/Raw_Data/WRAP/wrap_ConsensusConference_04302024.csv', na.strings = c("", " ", NA, "NA"))

# Rename variables across datasets ----
visits <- visits %>%
  dplyr::select(subject_id = WRAPNo, VisNo, age_at_appointment)

## rename the subject_id column
Tau$subject_id <- str_remove(Tau$subject_id, "sub-wrap")
Amyloid$subject_id <- str_remove(Amyloid$subject_id, "sub-wrap")
subject_id$subject_id <- str_remove(subject_id$subject_id, "sub-wrap")
visits$subject_id <- str_remove(visits$subject_id, "sub-wrap")
cognition <- cognition %>% dplyr::rename(subject_id = WRAPNo)
demo <- demo %>% dplyr::rename(subject_id = WRAPNo)

# Create unique age variables for each measure
Amyloid$age_at_appointment_amyloid <- Amyloid$age_at_appointment
Tau$age_at_appointment_tau <- Tau$age_at_appointment

# PACC --------------------------------------------------------------------

pacc <- cognition %>%
  dplyr::select(subject_id, VisNo, z_pacc4_ds_mmse) %>%
  filter(!is.na(z_pacc4_ds_mmse)) %>%
  group_by(subject_id) %>%
  arrange(subject_id, VisNo) %>%
  left_join(visits, by = c("subject_id", "VisNo")) %>%
  dplyr::mutate(pacc_age = age_at_appointment,
         bl_pacc_age = ifelse(is.na(pacc_age), NA, min(pacc_age, na.rm = TRUE))) %>%
  dplyr::select(-age_at_appointment)

# Amyloid ----
amyloid_new <- Amyloid %>%
  dplyr::select(subject_id, age_at_appointment, age_at_appointment_amyloid, pib_index) %>%
  group_by(subject_id) %>%
  arrange(subject_id, age_at_appointment) %>%
  mutate(pib_group = case_when(
           pib_index > 1.16 ~ "PIB+",
           pib_index < 1.16 ~ "PIB-",
           TRUE ~ NA)) %>%
  mutate(pib_index = case_when(
    pib_index == 0 ~ NA,  
    TRUE ~ pib_index
  )) %>%
  mutate(centiloid = pib_index * 148.33-154.96) 

# Clean tau data 
Tau$bl_Amyg <- (Tau$suvr_amygdala_l + Tau$suvr_amygdala_r)/2
Tau$bl_Para<- (Tau$suvr_parahippocampal_l+ Tau$suvr_parahippocampal_r)/2
Tau$bl_EC<- (Tau$suvr_parahippocampal_gy_ant_l + Tau$suvr_parahippocampal_gy_ant_r)/2
Tau$bl_Fus<- (Tau$suvr_fusiform_l+ Tau$suvr_fusiform_r)/2
Tau$bl_ITG<- (Tau$suvr_temporal_inf_l+ Tau$suvr_temporal_inf_r)/2
Tau$bl_MTG<- (Tau$suvr_temporal_mid_l+ Tau$suvr_temporal_mid_r)/2


Tau <- dplyr::select(Tau, subject_id, age_at_appointment, age_at_appointment_tau, 
                     bl_Amyg, bl_Para, bl_EC, bl_Fus, bl_ITG, bl_MTG)

## Baseline Tau
baseline_tau <- Tau %>%
  group_by(subject_id) %>%
  arrange(subject_id, age_at_appointment_tau) %>%
  slice(1)

# Demographics ----
subj_new <- subject_id %>%
  mutate(APOE = case_when(
    apoe_e1 == 4 | apoe_e2 == 4 ~ "e4+",
    is.na(apoe_e1) & is.na(apoe_e2) ~ NA,
    TRUE ~ "e4-")
  ) %>%
  dplyr::select(subject_id, APOE)

demo_new <- demo %>%
  mutate(sex = case_when(
    gender == 2 ~ "F",
    gender == 1 ~ "M"
  )) %>%
  dplyr::select(subject_id, EducYrs, sex)

## Join cross-sectional Pib closet to baseline tau
baseline_tau <- as.data.table(baseline_tau)
setkey(baseline_tau, subject_id, age_at_appointment)

amyloid_new <- as.data.table(amyloid_new)
setkey(amyloid_new, subject_id, age_at_appointment)

baseline_amyloid <- amyloid_new[baseline_tau, roll = "nearest"]
baseline_amyloid_tau <- dplyr::select(baseline_amyloid, -c(age_at_appointment))

# Merge --------
# tau and PACC
WRAP_data <- pacc %>%
  right_join(baseline_amyloid_tau, by = "subject_id") %>%
  left_join(demo_new, by = "subject_id") %>%
  left_join(subj_new, by = "subject_id") %>%
  mutate(Tau_Amyloid_Diff = as.numeric(age_at_appointment_tau - age_at_appointment_amyloid)) %>%
  group_by(subject_id) %>%
  mutate(Tau_PACC_Diff = as.numeric(bl_pacc_age - age_at_appointment_tau),
         pacc_time = pacc_age - age_at_appointment_tau,
         pacc_time_from_pacc = pacc_age - bl_pacc_age) %>%
  add_tally(name = "pacc_tp")

num <- unique(WRAP_data$subject_id) #499

# merge diagnosis
WRAP_data_final <- merge(
  WRAP_data,
  diagnosis[, c("WRAPNo", "VisNo", "Calculated_Consensus_dx")], 
  by.x = c("subject_id", "VisNo"), 
  by.y = c("WRAPNo", "VisNo"),
  all.x = TRUE
)

WRAP_data_final$Calculated_Consensus_dx<-as.factor(WRAP_data_final$Calculated_Consensus_dx)

# Construct flowchart ----
WRAP_data_all <- pacc %>%
  left_join(baseline_amyloid_tau, by = "subject_id") %>%
  left_join(demo_new, by = "subject_id") %>%
  left_join(subj_new, by = "subject_id")

WRAP_data_all <- merge(
  WRAP_data_all,
  diagnosis[, c("WRAPNo", "VisNo", "Calculated_Consensus_dx")], 
  by.x = c("subject_id", "VisNo"), 
  by.y = c("WRAPNo", "VisNo"),
  all.x = TRUE
)

WRAP_data_all_bl <- WRAP_data_all %>%
  group_by(subject_id) %>%
  arrange(subject_id, age_at_appointment_tau) %>%
  slice(1)

WRAP_fc <- WRAP_data_all_bl %>%
  as_fc(label = "WRAP Participants with Cognitive Data Available", text_pattern = "{label} \n N = {n}") %>%
  fc_filter(!is.na(age_at_appointment_tau), label = "Eligible Participants", text_pattern = "N = {n} {label}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Without Tau-PET Data") %>%
  fc_filter(Calculated_Consensus_dx == "Cog_Unimpaired_Stable" | Calculated_Consensus_dx == "Cog_Unimpaired_Declining", label = "Eligible WRAP Participants", text_pattern = "{label}\nN = {n}",
            show_exc = T, text_pattern_exc = "N = {n} {label}", label_exc = "Diagnosed With MCI At Baseline") %>%
  fc_draw()


## Clean data ----
WRAP_data <- WRAP_data_final %>% 
  dplyr::rename(PACC = z_pacc4_ds_mmse)

## make factor
WRAP_data$sex<-as.factor(WRAP_data$sex)
WRAP_data$APOE<-as.factor(WRAP_data$APOE)
WRAP_data$pib_group<-as.factor(WRAP_data$pib_group)
## Rename tau variables 
WRAP_data <- WRAP_data %>%
  dplyr::rename(
    bl_tau_age = age_at_appointment_tau,
    EDCCNTU = EducYrs,
    SEX = sex,
    ID = subject_id
  )

count_missing <- function(data, variables) {
  sapply(variables, function(var) sum(is.na(data[[var]])))
}
variables_to_check <- c(
  "PACC", "pacc_time", "bl_tau_age", "bl_Para", "bl_Amyg", "bl_EC", "bl_Fus",
  "bl_ITG", "bl_MTG", "centiloid", "YrsEd", "SEX", "APOE"
)

missing_counts <- count_missing(WRAP_data, variables_to_check)
missing_counts

WRAP <- WRAP_data %>% 
  filter(
    !is.na(PACC)
  )

n <- unique(WRAP$ID)


# Retain only cn subjects
baseline_data <- WRAP %>%
  arrange(ID, VisNo) %>%
  group_by(ID) %>%
  slice(1) %>%
  filter(Calculated_Consensus_dx == "Cog_Unimpaired_Stable" | Calculated_Consensus_dx == "Cog_Unimpaired_Declining")

WRAP$ID <- as.character(WRAP$ID)
baseline_data$ID <- as.character(baseline_data$ID)

cn_subjects <- baseline_data$ID

WRAP <- WRAP %>%
  filter(ID %in% cn_subjects)

WRAP$PACC_tp <- as.numeric(as.character(WRAP$pacc_tp))

WRAP <- WRAP %>%
  mutate(One_PACC_tp = case_when(
    PACC_tp == 1 ~ T,  
    PACC_tp > 1 ~ F)
  )

WRAP <- WRAP %>%
  group_by(ID) %>%
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

WRAP <- WRAP %>%
  group_by(ID) %>%
  mutate(
    tau_pacc_age_diff = abs(pacc_age - bl_tau_age),
    closest_PACC = PACC[which.min(tau_pacc_age_diff)]
  ) %>%
  ungroup()

n_sub <- unique(WRAP$ID)

save(WRAP, file = paste0("~/Desktop/sexcog/Preprocessing/WRAP", Sys.Date(), ".RData"))
