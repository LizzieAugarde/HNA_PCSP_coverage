################ HNA/PCSP coverage analysis - coverage by characteristics ###################

#Analysis for proposal aim 5 - proportion of people being offered an HNA/PCSP by 
#demographic and clinical characteristics 
#Used in Section 9 of initial results slide deck

#Created April 2024 by Lizzie Augarde 
##########################################################################################

library(scales)

patient_level_data <- patient_level_data |>
  mutate(hna_status = ifelse(hna_count != "0", "Has HNA", "No HNA"),
         pcsp_status = ifelse(pcsp_count != "0", "Has PCSP", "No PCSP"))

####by gender-------------------
gender_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(gender, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(gender) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(gender), hna_status == "Has HNA") |>
  mutate(gender = ifelse(gender == "1", "Male", "Female")) |>
  ungroup()

gender_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(gender, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(gender) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(gender), pcsp_status == "Has PCSP") |>
  mutate(gender = ifelse(gender == "1", "Male", "Female")) |>
  ungroup()

gender_hna_graph <- ggplot(gender_hna, aes(x = gender, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Gender", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

gender_pcsp_graph <- ggplot(gender_pcsp, aes(x = gender, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Gender", y = "Proportion of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#proportions with a HNA and no PCSP
gender_full_status <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE" & keepforhna == "INCLUDE") |>
  filter(hna_offered_code == "03") |>
  mutate(full_status = case_when(hna_status == "Has HNA" & pcsp_status == "No PCSP" ~ "HNA only", 
                                 hna_status == "Has HNA" & pcsp_status == "Has PCSP" ~ "Both", 
                                 TRUE ~"Other")) |>
  group_by(gender, full_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(gender) |>
  mutate(percent = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  mutate(gender = ifelse(gender == "1", "Male", "Female")) |>
  filter(!is.na(gender), full_status == "HNA only")



####by age at diag-------------------
patient_level_data <- patient_level_data |>
  mutate(age_group = case_when(age < 20 ~ "Under 20",
                               age < 30 & age > 19 ~ "20-29",
                               age < 40 & age > 29 ~ "30-39",
                               age < 50 & age > 39 ~ "40-49",
                               age < 60 & age > 49 ~ "50-59",
                               age < 70 & age > 59 ~ "60-69",
                               age < 80 & age > 69 ~ "70-79",
                               age > 79 ~ "80+", TRUE ~ NA)) |>
  mutate(age_group = factor(age_group, levels = c("Under 20", "20-29", "30-39", 
                                                  "40-49", "50-59", "60-69", 
                                                  "70-79", "80+")))

age_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(age_group, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(age_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(age_group), hna_status == "Has HNA") |>
  ungroup() 

age_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(age_group, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(age_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(age_group), pcsp_status == "Has PCSP") |>
  ungroup() 

age_hna_graph <- ggplot(age_hna, 
                        aes(x = age_group, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Age at diagnosis", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

age_pcsp_graph <- ggplot(age_pcsp,
                         aes(x = age_group, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Age at diagnosis", y = "Proportion of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#proportions with a HNA and no PCSP
age_full_status <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE" & keepforhna == "INCLUDE") |>
  filter(hna_offered_code == "03") |>
  mutate(full_status = case_when(hna_status == "Has HNA" & pcsp_status == "No PCSP" ~ "HNA only", 
                                 hna_status == "Has HNA" & pcsp_status == "Has PCSP" ~ "Both", 
                                 TRUE ~"Other")) |>
  group_by(age_group, full_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(age_group) |>
  mutate(percent = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(age_group), full_status == "HNA only")



####by ethnicity-------------------
patient_level_data <- patient_level_data |>
  mutate(ethnicity_group = case_when(ethnicity %in% c("A", "B", "C", "0") ~ "White",
                                     ethnicity %in% c("D", "E", "F", "G") ~ "Mixed",
                                     ethnicity %in% c("H", "J", "K", "L") ~ "Asian",
                                     ethnicity %in% c("M", "N", "P") ~ "Black",
                                     ethnicity == "R" ~ "Chinese", 
                                     ethnicity == "S" ~ "Other",
                                     ethnicity %in% c("Z", "X") ~ "Not known", TRUE ~ "Not known"))

ethnicity_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(ethnicity_group, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ethnicity_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ethnicity_group), hna_status == "Has HNA") |>
  ungroup()

ethnicity_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(ethnicity_group, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ethnicity_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ethnicity_group), pcsp_status == "Has PCSP") |>
  ungroup()

ethnicity_hna_graph <- ggplot(ethnicity_hna, 
                              aes(x = ethnicity_group, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Ethnicity", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ethnicity_pcsp_graph <- ggplot(filter(ethnicity_pcsp, ethnicity_pcsp$pcsp_status == "Has PCSP"), 
                               aes(x = ethnicity_group, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Ethnicity", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#proportions with a HNA and no PCSP
ethnicity_full_status <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE" & keepforhna == "INCLUDE") |>
  filter(hna_offered_code == "03") |>
  mutate(full_status = case_when(hna_status == "Has HNA" & pcsp_status == "No PCSP" ~ "HNA only", 
                                 hna_status == "Has HNA" & pcsp_status == "Has PCSP" ~ "Both", 
                                 TRUE ~"Other")) |>
  group_by(ethnicity_group, full_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ethnicity_group) |>
  mutate(percent = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(ethnicity_group), full_status == "HNA only")



####by deprivation-------------------
imd_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(imd19_decile_lsoas, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(imd19_decile_lsoas) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(imd19_decile_lsoas), hna_status == "Has HNA") |>
  ungroup() |>
  mutate(imd19_decile_lsoas = factor(imd19_decile_lsoas, levels = c(
    "1 - most deprived", "2", "3", "4", "5", 
    "6", "7", "8", "9", "10 - least deprived")))

imd_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(imd19_decile_lsoas, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(imd19_decile_lsoas) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(imd19_decile_lsoas), pcsp_status == "Has PCSP") |>
  ungroup() |>
  mutate(imd19_decile_lsoas = factor(imd19_decile_lsoas, levels = c(
    "1 - most deprived", "2", "3", "4", "5", 
    "6", "7", "8", "9", "10 - least deprived")))

imd_hna_graph <- ggplot(imd_hna, 
                        aes(x = imd19_decile_lsoas, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Deprivation decile", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

imd_pcsp_graph <- ggplot(imd_pcsp,
                         aes(x = imd19_decile_lsoas, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Deprivation decile", y = "Proportion of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#proportions with a HNA and no PCSP
imd_full_status <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE" & keepforhna == "INCLUDE") |>
  filter(hna_offered_code == "03") |>
  mutate(full_status = case_when(hna_status == "Has HNA" & pcsp_status == "No PCSP" ~ "HNA only", 
                                 hna_status == "Has HNA" & pcsp_status == "Has PCSP" ~ "Both", 
                                 TRUE ~"Other")) |>
  group_by(imd19_decile_lsoas, full_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(imd19_decile_lsoas) |>
  mutate(percent = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(imd19_decile_lsoas), full_status == "HNA only")



####by tumour type-----------
site_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE", ndrs_main != "to be grouped", ndrs_main != "See skin table") |>
  group_by(ndrs_main, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(suppress = ifelse(sum(number_patients)<1000, "Suppress", "Don't suppress")) 

site_suppression_hna <- site_hna |> select(ndrs_main, suppress) |> unique() #extracting suppression status for the HNA only table below

site_hna <- site_hna |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  mutate(number_patients_table = as.character(ifelse(suppress == "Suppress", "-", number_patients)),
         percent_table = ifelse(suppress == "Suppress", "-", percent_table),
         lower_table = ifelse(suppress == "Suppress", "-", lower_table),
         upper_table = ifelse(suppress == "Suppress", "-", upper_table)) |>
  filter(!is.na(ndrs_main), hna_status == "Has HNA") |>
  ungroup() 

site_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE", ndrs_main != "to be grouped", ndrs_main != "See skin table") |>
  group_by(ndrs_main, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(suppress = ifelse(sum(number_patients)<1000, "Suppress", "Don't suppress")) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  mutate(number_patients_table = as.character(ifelse(suppress == "Suppress", "-", number_patients)),
         percent_table = ifelse(suppress == "Suppress", "-", percent_table),
         lower_table = ifelse(suppress == "Suppress", "-", lower_table),
         upper_table = ifelse(suppress == "Suppress", "-", upper_table)) |>
  filter(!is.na(ndrs_main), pcsp_status == "Has PCSP") |>
  ungroup()

site_hna_graph <- ggplot(filter(site_hna, site_hna$suppress != "Suppress"), 
                        aes(x = ndrs_main, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Tumour type", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

site_pcsp_graph <- ggplot(filter(site_pcsp, site_pcsp$suppress != "Suppress"),  
                         aes(x = ndrs_main, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Tumour type", y = "Proportion of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#proportions with a HNA and no PCSP
site_full_status <- patient_level_data |>
  filter(ndrs_main != "to be grouped", ndrs_main != "See skin table") |>
  filter(keepforpcsp == "INCLUDE" & keepforhna == "INCLUDE") |>
  filter(hna_offered_code == "03") |>
  mutate(full_status = case_when(hna_status == "Has HNA" & pcsp_status == "No PCSP" ~ "HNA only", 
                                 hna_status == "Has HNA" & pcsp_status == "Has PCSP" ~ "Both", 
                                 TRUE ~"Other")) |>
  group_by(ndrs_main, full_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(percent = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  left_join(site_suppression_hna, by = "ndrs_main") |>
  mutate(number_patients = as.character(ifelse(suppress == "Suppress", "-", number_patients)),
         percent = ifelse(suppress == "Suppress", "-", percent)) |>
  select(-suppress) |>
  filter(!is.na(ndrs_main), full_status == "HNA only")



####by stage------------
stage_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(STAGE, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(STAGE) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(hna_status == "Has HNA", !is.na(STAGE), STAGE != "Error") |>
  ungroup() 

stage_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(STAGE, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(STAGE) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         proptest = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(proptest, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(proptest, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(pcsp_status == "Has PCSP", !is.na(STAGE), STAGE != "Error") |>
  ungroup()

stage_hna_graph <- ggplot(stage_hna, 
                         aes(x = STAGE, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Stage at diagnosis", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

stage_pcsp_graph <- ggplot(stage_pcsp, 
                          aes(x = STAGE, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Stage at diagnosis", y = "Proportion of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#proportions with a HNA and no PCSP
stage_full_status <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE" & keepforhna == "INCLUDE") |>
  filter(hna_offered_code == "03") |>
  mutate(full_status = case_when(hna_status == "Has HNA" & pcsp_status == "No PCSP" ~ "HNA only", 
                                 hna_status == "Has HNA" & pcsp_status == "Has PCSP" ~ "Both", 
                                 TRUE ~"Other")) |>
  group_by(STAGE, full_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(STAGE) |>
  mutate(percent = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(STAGE), full_status == "HNA only")

