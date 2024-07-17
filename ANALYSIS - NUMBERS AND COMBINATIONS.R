################ HNA/PCSP coverage analysis - numbers and combos of HNAS/PCSPs ###################

#Analysis for proposal aim XXXX - people have more than 1 HNA/PCSP and people having both

#Created April 2024 by Lizzie Augarde 
#Change log:
#27/06/2024 reformatted percentage columns
##########################################################################################

library(scales)

###### MULTIPLE HNAS/PCSPs (POST-DEDUPLICATION) ######
patient_level_data <- patient_level_data %>%
    mutate(hna_count = as.character(hna_count),
           hna_count = ifelse(!hna_count %in% c("0", "1", "2"), "3 or more", hna_count),
           pcsp_count = as.character(pcsp_count),
           pcsp_count = ifelse(!pcsp_count %in% c("0", "1", "2"), "3 or more", pcsp_count))

hna_count_summary <- patient_level_data %>%
  group_by(hna_count) %>%
  summarise(number_patients_hna = n()) %>%
  ungroup() %>%
  mutate(percent_patients_hna = percent((number_patients_hna/sum(number_patients_hna)), accuracy = 0.1)) 

pcsp_count_summary <- patient_level_data %>%
  group_by(pcsp_count) %>%
  summarise(number_patients_pcsp = n()) %>%
  ungroup() %>%
  mutate(percent_patients_pcsp = percent((number_patients_pcsp/sum(number_patients_pcsp)), accuracy = 0.1)) 


