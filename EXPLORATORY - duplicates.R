################ HNA/PCSP coverage analysis - dupes ###################

#Investigating duplicate records

#Created March 2024 by Lizzie Augarde 
############################################################### 

library(xlsx) 

#individual patients review
pat_data <- dbGetQueryOracle(cas2401, "select * from analysisncr.at_rapid_pathway@cas2401 where patientid = '46753586'", rowlimit = NA)


#by Trust
raw_data <- dbGetQueryOracle(cas2401, "SELECT * FROM HNA_PCSP_2021DIAGS", rowlimit = NA)

raw_data <- raw_data %>% 
  clean_names() %>%
  filter(age > 17) %>%
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N"))

dups <- duplicated(raw_data) | duplicated(raw_data, fromLast = TRUE)
dups <- dups & !duplicated(raw_data)
raw_data$dup <- dups

raw_data_dups_hna <- raw_data %>%
  filter(hna == "Y") %>%
  mutate(dup = ifelse(dup == "FALSE", 1, 0)) %>%
  group_by(diag_trust) %>%
  summarise(dups_count = sum(dup)-1, total = n()) %>%
  ungroup() %>%
  mutate(prop_dups = (dups_count/total)*100)

raw_data_dups_pcsp <- raw_data %>%
  filter(pcsp == "Y") %>%
  mutate(dup = ifelse(dup == "FALSE", 1, 0)) %>%
  group_by(diag_trust) %>%
  summarise(dups_count = sum(dup)-1, total = n()) %>%
  ungroup() %>%
  mutate(prop_dups = (dups_count/total)*100)

write.xlsx(raw_data_dups_hna, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Results/Dups by trust 20240327.xlsx", sheetName = "HNA")
write.xlsx(raw_data_dups_pcsp, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Results/Dups by trust 20240327.xlsx", sheetName = "PCSP", append = TRUE)



##### Records within a week of each other at the same trust ######
hnas_within_week <- hna_pcsp_data |>
  select(patientid, hna, event_date, trust_code) |>
  filter(hna == "Y") |>
  arrange(patientid, trust_code, event_date) |>
  group_by(patientid, trust_code) |>
  mutate(days_between_hnas = difftime(event_date, lag(event_date), units = "days")) |>
  ungroup()

length(which(hnas_within_week$days_between_hnas<8)) ###6% 

pcsps_within_week <- hna_pcsp_data |>
  select(patientid, pcsp, event_date, trust_code) |>
  filter(pcsp == "Y") |>
  arrange(patientid, trust_code, event_date) |>
  group_by(patientid, trust_code) |>
  mutate(days_between_pcsps = difftime(event_date, lag(event_date), units = "days")) |>
  ungroup() 

length(which(pcsps_within_week$days_between_pcsps<8)) ###6% 
