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
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N"))

sum(raw_data$hna == "Y" | raw_data$pcsp == "Y") #1918680 records of a HNA/PCSP

hna_pcsp_data <- raw_data %>%
  unique() %>% #removing duplicate rows
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") 

length(unique(hna_pcsp_data$patientid)) #311627 patients (329665 diagnoses in NS publication, remember DCO and C44)
sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #272307 unique records of a HNA/PCSP


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
no_events <- hna_pcsp_data %>% filter(hna == "N" & pcsp == "N") #216669 patients with no events
hna_data <- hna_pcsp_data %>% filter(hna == "Y") #151216 HNA records
pcsp_data <- hna_pcsp_data %>% filter(pcsp == "Y") #121091 PCSP records


#overall numbers of patients with HNAs/PCSPs/combination
length(unique(no_events$patientid)) #check these patients are unique 
length(unique(hna_data$patientid)) #85798 patients with at least one HNA record
length(unique(pcsp_data$patientid)) #81216 patients with at least one PCSP record

ids_noevents <- unique(no_events$patientid)
ids_hna <- unique(hna_data$patientid)
ids_pcsp <- unique(pcsp_data$patientid)

only_noevents <- setdiff(ids_noevents, union(ids_hna, ids_pcsp)) #matches no_events row number, correct
in_noevents_and_hna_not_pcsp <- intersect(setdiff(ids_noevents, ids_pcsp), setdiff(ids_hna, ids_pcsp)) #empty, correct
in_noevents_and_pcsp_not_hna <- intersect(setdiff(ids_noevents, ids_hna), setdiff(ids_pcsp, ids_hna)) #empty, correct
only_hna <- setdiff(ids_hna, ids_pcsp) #13742 patients with only an HNA
only_pcsp <- setdiff(ids_pcsp, ids_hna) #9160 patients with only a PCSP
both <- intersect(ids_hna, ids_pcsp) #72056 patients with both 


#write out record level HNA and PCSP data, and unique patients with neither
write.csv(no_events, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patients with no events raw data RCRD 20240315.csv")
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA events raw data RCRD 20240315.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSP events raw data RCRD 20240315.csv")


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
write.csv(no_events, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA and PCSPs in RCRD ranked by patient 20240315.csv")


####### Patient-level combined dataset - row for each patient with record of earliest HNA and PCSP for each person ######
hna_data <- hna_data %>%
  group_by(patientid) %>%
  mutate(hna_count = n())
  
  
include count of HNAs and PCSPs per patient 
