
library(tidyverse)

patient_level_data <- read.csv("N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level RCRD HNA and PCSP data 20240408.csv")

trust_hna <- patient_level_data %>%
  filter(keepforhna == "INCLUDE") %>%
  group_by(diag_trust, hna_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(diag_trust) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(diag_trust)) %>%
  ungroup()

trust_pcsp <- patient_level_data %>%
  filter(keepforpcsp == "INCLUDE") %>%
  group_by(diag_trust, pcsp_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(diag_trust) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(diag_trust)) %>%
  ungroup()

write.csv(trust_hna, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level HNA status by trust 20240508.csv")
write.csv(trust_pcsp, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level PCSP status by trust 20240508.csv")

