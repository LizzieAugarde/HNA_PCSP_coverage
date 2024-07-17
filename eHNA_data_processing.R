########## HNA 2021 coverage project - eHNA data processing ###########

#Script to process eHNA data into aggregate by trust and month
#for use in comparison to COSD for coverage project 2021

#Created July 2024 by Lizzie Augarde
#Change log:
######################################################################

library(tidyverse)
library(vroom)
library(janitor)

#read in
path <- "C:/Users/LAugarde/OneDrive - Macmillan Cancer Support/Projects/414 - HNA PCSP coverage analysis - NCRAS partnership project/eHNA data"
ehna_data_files <- list.files(path, full.names = TRUE)
ehna_data <- vroom(ehna_data_files)

#processing eHNA data
ehna_data <- ehna_data %>%
  unique() %>%
  filter(substr(DiagnosisCode, 1, 1) == "C" & DiagnosisCode != "C44") %>%
  select(c(AssessmentId, DateSetUpYear, DateSetUpMonth, Status, ODS, Organisation)) %>%
  filter(!is.na(AssessmentId)) %>%
  mutate(month_text = month.abb[DateSetUpMonth],
         monthyear_text = paste0(month_text, "-", (str_sub(as.character(DateSetUpYear), 3, 4)))) %>%
  filter(Status %in% c("locked", "submitted", "in progress")) %>% #matching modelled HNA process
  filter(str_sub(ODS, 1, 1) == "R") %>% #England trusts only 
  group_by(ODS, monthyear_text) %>%
  summarise(count = n()) %>%
  pivot_wider(., names_from = monthyear_text, values_from = count) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

#reorder columns 
ehna_data <- ehna_data[, c(1,12,11,13,2,7,6,5,3,10,9,8,4)]

#add England total by month and 2021 trust-level total
ehna_data <- ehna_data %>%
  adorn_totals(., where = "row", name = "England") %>%
  adorn_totals(., where = "col", name = "total_2021") 

#write out
write.csv(ehna_data, paste0(path, "/eHNA_data_2021.")) 

