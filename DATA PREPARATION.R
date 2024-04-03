################ HNA/PCSP coverage analysis - data prep ###################

#Data preparation and cleaning of COSD Level 3 registrations for 
#2021 linked to RCRD HNA and PCSP records for coverage analysis

#Created November 2023 by Lizzie Augarde 
#Change log:
#15/03/2024 adjusted approach to counting patients
############################################################### 

#prep
library(NDRSAfunctions)
library(tidyverse)
library(janitor)

casref01 <- createConnection()
cas2401 <- createConnection(port = 1525, sid = "cas2401")

#pulling from SQL 
raw_data <- dbGetQueryOracle(cas2401, "SELECT * FROM HNA_PCSP_2021DIAGS", rowlimit = NA)

raw_data <- raw_data %>% 
  clean_names() %>%
  filter(age > 17) %>%
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N")) %>%
  
  #limiting to records within 2 years of diagnosis 
  mutate(time_diag_event = difftime(as.Date(event_date), as.Date(diagnosisdatebest), units = "days")) %>%
  filter((time_diag_event < 731 & time_diag_event >= 0) | is.na(time_diag_event))

sum(raw_data$hna == "Y" | raw_data$pcsp == "Y") #1752833 records of a HNA/PCSP

hna_pcsp_data <- raw_data %>%
  unique() %>% #removing duplicate rows
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") 

length(unique(hna_pcsp_data$patientid)) #309870 patients########not getting the same number of patients, need to work out why 
sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #272217 unique records of a HNA/PCSP
sum(hna_pcsp_data$hna == "Y") #151171 unique HNA records
sum(hna_pcsp_data$pcsp == "Y") #121046 unique PCSP records

#counting then excluding records with offered code 04 'Not offered'
sum(hna_pcsp_data$hna == "Y" & hna_pcsp_data$offered_code == "04") #479 not offered HNA records
sum(hna_pcsp_data$pcsp == "Y" & hna_pcsp_data$offered_code == "04") #9511 not offered PCSP records

hna_pcsp_data <- hna_pcsp_data %>%
  filter(offered_code != "04") #counting NA offered_code as offered for now

sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #262227 unique records of an OFFERED HNA/PCSP 
sum(hna_pcsp_data$hna == "Y") #150692 unique OFFERED HNA records
sum(hna_pcsp_data$pcsp == "Y") #111535 unique OFFERED PCSP records

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
no_events <- hna_pcsp_data %>% filter(hna == "N" & pcsp == "N") #214969 patients with no events
hna_data <- hna_pcsp_data %>% filter(hna == "Y") #151171 HNA records
pcsp_data <- hna_pcsp_data %>% filter(pcsp == "Y") #121046 PCSP records

#overall numbers of patients with HNAs/PCSPs/combination
length(unique(hna_data$patientid)) #85756 patients with at least one HNA record
length(unique(pcsp_data$patientid)) #81174 patients with at least one PCSP record
length(unique(rbind(hna_data, pcsp_data)$patientid))#94901 patients with either

ids_noevents <- unique(no_events$patientid)
ids_hna <- unique(hna_data$patientid)
ids_pcsp <- unique(pcsp_data$patientid)

only_noevents <- setdiff(ids_noevents, union(ids_hna, ids_pcsp)) #matches no_events row number, correct
in_noevents_and_hna_not_pcsp <- intersect(setdiff(ids_noevents, ids_pcsp), setdiff(ids_hna, ids_pcsp)) #empty, correct
in_noevents_and_pcsp_not_hna <- intersect(setdiff(ids_noevents, ids_hna), setdiff(ids_pcsp, ids_hna)) #empty, correct
only_hna <- setdiff(ids_hna, ids_pcsp) #13727 patients with only an HNA
only_pcsp <- setdiff(ids_pcsp, ids_hna) #9145 patients with only a PCSP
both <- intersect(ids_hna, ids_pcsp) #72029 patients with both 

#write out record level HNA and PCSP data, and unique patients with neither
write.csv(no_events, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patients with no events raw data RCRD 20240318.csv")
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA events raw data RCRD 20240318.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSP events raw data RCRD 20240318.csv")


###### Adding ranking of HNAs and PCSPs by date (earliest ranked 1) ######
rank_by_date <- function(data, pid_col, date_col) {
  data %>%
    group_by({{ pid_col }}) %>%
    arrange({{ date_col }}) %>%
    mutate(rank = row_number())
}

hna_data <- rank_by_date(hna_data, patientid, event_date)
pcsp_data <- rank_by_date(pcsp_data, patientid, event_date)

hna_pcsp_data <- rbind(hna_data, pcsp_data) #this dataset now contains all unique and ranked HNA and PCSP records

#write out record level ranked HNA and PCSP data
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA in RCRD ranked by patient 20240315.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSPs in RCRD ranked by patient 20240315.csv")
write.csv(hna_pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA and PCSPs in RCRD ranked by patient 20240315.csv")


####### Patient-level combined dataset - row for each patient with record of earliest HNA and PCSP for each person ######
#count of HNAs and PCSPs per patient
hna_count <- hna_data %>%
  group_by(patientid) %>%
  mutate(hna_count = n()) %>%
  select(patientid, hna_count) %>%
  unique()

pcsp_count <- pcsp_data %>%
  group_by(patientid) %>%
  mutate(pcsp_count = n()) %>%
  select(patientid, pcsp_count) %>%
  unique()

hna_data <- hna_data %>%
  filter(rank == 1) %>%
  rename("hna_event_type" = "event_type", "hna_point_of_pathway" = "point_of_pathway", "hna_staff_role" = "staff_role", "hna_date" = "event_date")

pcsp_data <- pcsp_data %>% 
  filter(rank == 1) %>%
  select(c(patientid, event_type, point_of_pathway, offered_code, staff_role, event_date)) %>%
  rename("pcsp_event_type" = "event_type", "pcsp_point_of_pathway" = "point_of_pathway", "pcsp_staff_role" = "staff_role", "pcsp_date" = "event_date")

patient_level_data <- inner_join(hna_data, pcsp_data, by = "patientid")

hna_only_patient_level <- anti_join(hna_data, patient_level_data, by = "patientid") %>%
  mutate(pcsp_event_type = NA, pcsp_point_of_pathway = NA, pcsp_staff_role = NA, pcsp_date = NA)

pcsp_only_patient_level <- anti_join(pcsp_data, patient_level_data, by = "patientid") %>%
  mutate(hna_event_type = NA, hna_point_of_pathway = NA, hna_staff_role = NA, hna_date = NA)

no_events <- no_events %>%
  mutate(hna_event_type = NA, hna_point_of_pathway = NA, hna_staff_role = NA, hna_date = NA,
         pcsp_event_type = NA, pcsp_point_of_pathway = NA, pcsp_staff_role = NA, pcsp_date = NA)

#combining all groups of patients into one patient-level data frame 
patient_level_data <- rbind(patient_level_data, hna_only_patient_level, pcsp_only_patient_level, no_events) 

#adding in HNA and PCSP counts for each patient
patient_level_data <- left_join(patient_level_data, hna_count, by = "patientid")
patient_level_data <- left_join(patient_level_data, pcsp_count, by = "patientid")

patient_level_data <- patient_level_data %>%
  mutate(hna_count = ifelse(is.na(hna_count), 0, hna_count),
         pcsp_count = ifelse(is.na(pcsp_count), 0, pcsp_count))
  

#write out patient-level HNA and PCSP data
write.csv(patient_level_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level RCRD HNA and PCSP data 20240403.csv")

