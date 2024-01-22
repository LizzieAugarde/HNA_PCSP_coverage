################ HNA/PCSP coverage analysis ###################

#Data preparation and cleaning of COSD Level 3 registrations for 
#2021 linked to RCRD HNA and PCSP records for coverage analysis

#Created November 2023 by Lizzie Augarde 
#Change log:
############################################################### 

#prep
library(NDRSAfunctions)
library(tidyverse)
library(janitor)

casref01 <- createConnection()

#pulling from SQL
raw_data <- dbGetQueryOracle(casref01, "SELECT * FROM hna_pcsp_2021diags", rowlimit = NA)

raw_data <- clean_names(raw_data)

hna_pcsp_data <- raw_data %>% unique() 
#need to redo all this NUMBERS DON'T MAKE SENSE
hna_pcsp_data <- hna_pcsp_data %>%
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") %>%
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N"))

length(unique(hna_pcsp_data$patientid)) #311627 patients (329665 diagnoses in NS publication, remember DCO and C44)

#splitting into HNA/PCSP/neither
no_events <- hna_pcsp_data %>% filter(hna == "N" & pcsp == "N") %>% unique()
hna_data <- hna_pcsp_data %>% filter(hna == "Y") %>% unique()
pcsp_data <- hna_pcsp_data %>% filter(pcsp == "Y") %>% unique()

write.csv(no_events, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patients with no events raw data RCRD 20240116.csv")
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA events raw data RCRD 20240116.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSP events raw data RCRD 20240116.csv")

#adding IMD data 
query <- "select a.patientid,
                 g.lsoa11_code,
                 d.imd19_decile_lsoas,
                 d.imd19_quintile_lsoas
          from av2021.at_tumour_england@casref01 a 
          left join av2021.at_geography_england@casref01 g on a.tumourid = g.tumourid
          left join imd.imd2019_equal_lsoas@casref01 d on g.lsoa11_code = d.lsoa11_code
          where a.diagnosisyear = '2021'"#need this to select the most recent deprivation 

imd_data <- dbGetQueryOracle(casref01, query, rowlimit = NA)

imd_data <- imd_data %>% select(-LSOA11_CODE) %>% unique() 

hna_pcsp_data <- left_join(hna_pcsp_data, imd_data, by = c("patientid" = "PATIENTID"), relationship = "one-to-many")


#keeping earliest HNA/PCSP for each patient
need to rank across the hnas/pcsps for each person to pick the record on each date with the most info

hna_data <- hna_data %>%
  group_by(patientid) %>%
  slice_min(order_by = as.Date(event_date))

pcsp_data <- pcsp_data %>%
  group_by(patientid) %>%
  slice_min(order_by = as.Date(event_date))
