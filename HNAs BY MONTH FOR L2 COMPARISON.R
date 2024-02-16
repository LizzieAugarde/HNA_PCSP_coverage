################ HNA/PCSP coverage analysis - HNAs by month/Trust ###################

#HNAs by month and Trust for comparison with the Level 2 COSD 
#extracts from CS2. Uses raw_data extracted in DATA PREPARATION script

#Created February 2024 by Lizzie Augarde 
#Change log:
############################################################### 

library(xlsx)

##### National-level #####
raw_data_l2_comp <- raw_data %>% 
  filter(event_type == 20)

hna_data_all <- raw_data_l2_comp %>% #not removing duplicates
  clean_names() %>%
  mutate(event_year = year(as.Date(event_date))) %>%
  mutate(event_month = month(as.Date(event_date)))

hnas_by_month_all <- hna_data_all %>%
  group_by(event_year, event_month) %>%
  summarise(count = n())

hna_data_unique <- raw_data_l2_comp %>% #unique HNA records
  unique() %>%
  clean_names() %>%
  mutate(event_year = year(as.Date(event_date))) %>%
  mutate(event_month = month(as.Date(event_date)))

hnas_by_month_unique <- hna_data_unique %>%
  group_by(event_year, event_month) %>%
  summarise(count = n())

write.xlsx(as.data.frame(hnas_by_month_all), "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Results/HNAs by month for L2 comparison 20240216.xlsx", 
           sheetName = "All records", row.names = FALSE)

write.xlsx(as.data.frame(hnas_by_month_unique), "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Results/HNAs by month for L2 comparison 20240216.xlsx", 
           sheetName = "Unique records", row.names = FALSE, append = TRUE)


##### By Trust #####
hnas_by_month_all_trust <- hna_data_all %>%
  group_by(event_year, event_month, diag_trust) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(event_ym = paste0(event_year, "-", event_month)) %>%
  select(-c(event_year, event_month)) %>%
  pivot_wider(names_from = event_ym, values_from = count) %>%
  mutate_if(is.numeric , replace_na, replace = 0)

hnas_by_month_unique_trust <- hna_data_unique %>%
  group_by(event_year, event_month, diag_trust) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(event_ym = paste0(event_year, "-", event_month)) %>%
  select(-c(event_year, event_month)) %>%
  pivot_wider(names_from = event_ym, values_from = count) %>%
  mutate_if(is.numeric , replace_na, replace = 0)

write.xlsx(as.data.frame(hnas_by_month_all_trust), "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Results/HNAs by month for L2 comparison 20240216.xlsx", 
           sheetName = "All records by Trust", row.names = FALSE, append = TRUE)

write.xlsx(as.data.frame(hnas_by_month_unique_trust), "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Results/HNAs by month for L2 comparison 20240216.xlsx", 
           sheetName = "Unique records by Trust", row.names = FALSE, append = TRUE)
