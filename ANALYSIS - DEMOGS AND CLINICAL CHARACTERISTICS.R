################ HNA/PCSP coverage analysis - coverage by characteristics ###################

#Analysis for proposal aim XXXX - proportion of people being offered an HNA/PCSP by 
#demographic and clinical characteristics 

#Created April 2024 by Lizzie Augarde 
#Change log:
##########################################################################################


###### DEMOGRAPHICS (POST-DEDUPLICATION) ######
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
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(gender)) |>
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
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(gender)) |>
  mutate(gender = ifelse(gender == "1", "Male", "Female")) |>
  ungroup()

gender_hna_graph <- ggplot(filter(gender_hna, gender_hna$hna_status == "Has HNA"), aes(x = gender, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Gender", y = "Percentage of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

gender_pcsp_graph <- ggplot(filter(gender_pcsp, gender_pcsp$pcsp_status == "Has PCSP"), aes(x = gender, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Gender", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


####by age at diag-------------------
patient_level_data <- patient_level_data |>
  mutate(age_group = case_when(age < 20 ~ "Under 20",
                               age < 30 & age > 19 ~ "20-29",
                               age < 40 & age > 29 ~ "30-39",
                               age < 50 & age > 39 ~ "40-49",
                               age < 60 & age > 49 ~ "50-59",
                               age < 70 & age > 59 ~ "60-69",
                               age < 80 & age > 69 ~ "70-79",
                               age > 79 ~ "80+", TRUE ~ NA))

age_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(age_group, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(age_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(age_group)) |>
  ungroup() 

age_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(age_group, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(age_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(age_group)) |>
  ungroup() 

age_hna_graph <- ggplot(filter(age_hna, age_hna$hna_status == "Has HNA"), 
       aes(x = factor(age_group, levels = c("Under 20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")),
           y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Age at diagnosis", y = "Percentage of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

age_pcsp_graph <- ggplot(filter(age_pcsp, age_pcsp$pcsp_status == "Has PCSP"),
       aes(x = factor(age_group, levels = c("Under 20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")),
           y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Age at diagnosis", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


####by ethnicity-------------------
patient_level_data <- patient_level_data |>
  mutate(ethnicity_group = case_when(ethnicity %in% c("A", "B", "C") ~ "White",
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
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ethnicity_group)) |>
  ungroup()

ethnicity_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(ethnicity_group, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ethnicity_group) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ethnicity_group)) |>
  ungroup()

ethnicity_hna_graph <- ggplot(filter(ethnicity_hna, ethnicity_hna$hna_status == "Has HNA"), 
                              aes(x = ethnicity_group, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Ethnicity", y = "Percentage of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ethnicity_pcsp_graph <- ggplot(filter(ethnicity_pcsp, ethnicity_pcsp$pcsp_status == "Has PCSP"), 
                               aes(x = ethnicity_group, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Ethnicity", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


####by deprivation-------------------
imd_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(imd19_decile_lsoas, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(imd19_decile_lsoas) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(imd19_decile_lsoas)) |>
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
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(imd19_decile_lsoas)) |>
  ungroup() |>
  mutate(imd19_decile_lsoas = factor(imd19_decile_lsoas, levels = c(
    "1 - most deprived", "2", "3", "4", "5", 
    "6", "7", "8", "9", "10 - least deprived")))

imd_hna_graph <- ggplot(filter(imd_hna, imd_hna$hna_status == "Has HNA"), 
                              aes(x = imd19_decile_lsoas, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Deprivation decile", y = "Percentage of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

imd_pcsp_graph <- ggplot(filter(imd_pcsp, imd_pcsp$pcsp_status == "Has PCSP"), 
                               aes(x = imd19_decile_lsoas, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Deprivation decile", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))



####by tumour type-----------
site_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(ndrs_main, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ndrs_main)) |>
  ungroup() 

site_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(ndrs_main, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ndrs_main)) |>
  ungroup()

site_hna_graph <- ggplot(filter(site_hna, site_hna$hna_status == "Has HNA"), 
                        aes(x = ndrs_main, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Tumour type", y = "Percentage of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

site_pcsp_graph <- ggplot(filter(site_pcsp, site_pcsp$pcsp_status == "Has PCSP"), 
                         aes(x = ndrs_main, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Tumour type", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#by stage------------
site_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(ndrs_main, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ndrs_main)) |>
  ungroup() 

site_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(ndrs_main, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ndrs_main) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(!is.na(ndrs_main)) |>
  ungroup()

site_hna_graph <- ggplot(filter(site_hna, site_hna$hna_status == "Has HNA"), 
                         aes(x = ndrs_main, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Tumour type", y = "Percentage of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

site_pcsp_graph <- ggplot(filter(site_pcsp, site_pcsp$pcsp_status == "Has PCSP"), 
                          aes(x = ndrs_main, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Tumour type", y = "Percentage of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))