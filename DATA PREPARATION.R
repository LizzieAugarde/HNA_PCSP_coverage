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
cas2401 <- createConnection(port = 1525, sid = "cas2401")

#pulling from SQL 
raw_data <- dbGetQueryOracle(cas2401, "SELECT * FROM HNA_PCSP_2021DIAGS", rowlimit = NA)

raw_data <- clean_names(raw_data)

hna_pcsp_data <- raw_data %>% unique() #488976 unique rows 

hna_pcsp_data <- hna_pcsp_data %>%
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") %>%
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N"))

length(unique(hna_pcsp_data$patientid)) #311627 patients (329665 diagnoses in NS publication, remember DCO and C44)


#adding IMD data 
query <- "select a.patientid,
                 g.lsoa11_code,
                 d.imd19_decile_lsoas,
                 d.imd19_quintile_lsoas
          from analysiselizabethaugarde.hna_pcsp_2021diags@cas2401 a 
          left join av2021.at_geography_england@casref01 g on a.tumourid = g.tumourid
          left join imd.imd2019_equal_lsoas@casref01 d on g.lsoa11_code = d.lsoa11_code"

imd_data <- dbGetQueryOracle(casref01, query, rowlimit = NA)

imd_data <- imd_data %>% select(-LSOA11_CODE) %>% unique() 

hna_pcsp_data <- left_join(hna_pcsp_data, imd_data, by = c("patientid" = "PATIENTID"), relationship = "many-to-one")


#splitting into HNA/PCSP/neither
no_events <- hna_pcsp_data %>% filter(hna == "N" & pcsp == "N") %>% unique() #216669 patients with no events
hna_data <- hna_pcsp_data %>% filter(hna == "Y") %>% unique() #151216 HNA records
pcsp_data <- hna_pcsp_data %>% filter(pcsp == "Y") %>% unique() #121091 PCSP records

length(unique(no_events$patientid)) #check these patients are unique 
length(unique(hna_data$patientid)) #85798 unique patients with an HNA record
length(unique(pcsp_data$patientid)) #81216 unique patients with a PCSP record

write.csv(no_events, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patients with no events raw data RCRD 20240216.csv")
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA events raw data RCRD 20240216.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSP events raw data RCRD 20240216.csv")

