
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


lsoa_hna <- patient_level_data %>%
  filter(keepforhna == "INCLUDE") %>%
  group_by(lsoa11_code, hna_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(lsoa11_code) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(lsoa11_code)) %>%
  ungroup()

lsoa_pcsp <- patient_level_data %>%
  filter(keepforpcsp == "INCLUDE") %>%
  group_by(lsoa11_code, pcsp_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(lsoa11_code) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(lsoa11_code)) %>%
  ungroup()

write.csv(lsoa_hna, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level HNA status by LSOA 20240508.csv")
write.csv(lsoa_pcsp, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level PCSP status by LSOA 20240508.csv")


parl_hna <- patient_level_data %>%
  filter(keepforhna == "INCLUDE") %>%
  group_by(parl_con_name, hna_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(parl_con_name) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(parl_con_name)) %>%
  ungroup()

parl_pcsp <- patient_level_data %>%
  filter(keepforpcsp == "INCLUDE") %>%
  group_by(parl_con_name, pcsp_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(parl_con_name) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(parl_con_name)) %>%
  ungroup()

write.csv(parl_hna, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level HNA status by constituency 20240508.csv")
write.csv(parl_pcsp, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level PCSP status by constituency 20240508.csv")


utla_hna <- patient_level_data %>%
  filter(keepforhna == "INCLUDE") %>%
  group_by(utla_2021_name, hna_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(utla_2021_name) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(utla_2021_name)) %>%
  ungroup()

utla_pcsp <- patient_level_data %>%
  filter(keepforpcsp == "INCLUDE") %>%
  group_by(utla_2021_name, pcsp_status) %>%
  summarise(number_patients = n()) %>%
  ungroup() %>%
  group_by(utla_2021_name) %>%
  mutate(percent_patients = (number_patients / sum(number_patients)) * 100,
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = (sapply(lower, function(x) x$conf.int[2]))*100, 
         lower = (sapply(lower, function(x) x$conf.int[1]))*100) %>%
  filter(!is.na(utla_2021_name)) %>%
  ungroup()

write.csv(utla_hna, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level HNA status by constituency 20240508.csv")
write.csv(utla_pcsp, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level PCSP status by constituency 20240508.csv")


#boxplot 
parl_boxplot <- ggplot(filter(parl_hna, hna_status == "Has HNA"), aes(x = hna_status, y = percent_patients)) +
  geom_boxplot(outlier.size = 1.1) +
  ylab("Percentage of patients") +
  xlab("") +
  scale_y_continuous(limits = c(0,100))

parl_boxplot
