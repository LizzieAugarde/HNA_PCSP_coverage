################ HNA/PCSP coverage analysis - data prep ###################

#Data preparation and cleaning of COSD Level 3 registrations for 
#2021 linked to RCRD HNA and PCSP records for coverage analysis

#Created November 2023 by Lizzie Augarde 
#Change log:
#15/03/2024 adjusted approach to counting patients
#08/04/2024 adjusted query approach
############################################################### 

#prep
library(NDRSAfunctions)
library(tidyverse)
library(janitor)

cas2403 <- createConnection(port = 1525, sid = "cas2403")

###### EXTRACTING AND CLEANING RAW DATA #######
#patient cohort----
query <- "select * from (
  select pat.patientid,
  tum.tumourid,
  tum.diagnosisdatebest,
  tum.diagnosisyear,
  tum.site_icd10_o2_3char,
  tum.stage_best,
  tum.gender,
  tum.age,
  tum.ethnicity,
  tum.deathdatebest,
  tum.diag_trust,
  geo.lsoa11_code,
  geo.parl_con_name,
  geo.utla_2021_name,
  imd.imd19_decile_lsoas,
  rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join av2021.at_geography_england@casref01 geo on tum.tumourid = geo.tumourid
  left join imd.imd2019_equal_lsoas@casref01 imd on geo.lsoa11_code = imd.lsoa11_code
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  and tum.age >17
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44') 
where rank = 1"

patient_cohort <- dbGetQueryOracle(cas2403, query, rowlimit = NA)

patient_cohort <- patient_cohort %>% 
  clean_names() %>%
  unique()

#HNA and PCSP records------
query <- "select * from (
  select pat.patientid,
  rp.nhsnumber as rp_nhsnumber,
  rp.patientid as rp_patientid,       
  tum.tumourid,
  tum.diagnosisdatebest,
  tum.diagnosisyear,
  tum.diag_trust,
  rp.event_type,
  rp.event_property_1,
  rp.event_property_2,
  rp.event_property_3,
  rp.event_date,
  rp.event_end,
  rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join analysisncr.at_rapid_pathway@cas2403 rp on pat.nhsnumber = rp.nhsnumber 
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  and tum.age >17
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44'
  and rp.event_type in (20, 24))
where rank = 1"

raw_hna_pcsp_data <- dbGetQueryOracle(cas2403, query, rowlimit = NA)

raw_hna_pcsp_data <- raw_hna_pcsp_data %>% 
  clean_names() %>%

  #limiting to records within 2 years of diagnosis 
  mutate(time_diag_event = difftime(as.Date(event_date), as.Date(diagnosisdatebest), units = "days")) %>%
  filter((time_diag_event < 731 & time_diag_event >= 0)) %>%
  
  #marking HNAs and PCSPs
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N")) %>%
  
  #joining to patient cohort to get demographic info 
  left_join(., patient_cohort %>% select(patientid, site_icd10_o2_3char, stage_best, gender, age, ethnicity, deathdatebest, lsoa11_code, imd19_decile_lsoas, parl_con_name), 
            by = "patientid")

sum(raw_hna_pcsp_data$hna == "Y" | raw_hna_pcsp_data$pcsp == "Y") #1815300 records of a HNA/PCSP


###### REMOVING DUPLICATES, INITIAL NUMBERS #######
hna_pcsp_data <- raw_hna_pcsp_data %>%
  unique() %>% #removing duplicate rows
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") 

length(unique(patient_cohort$patientid)) #309870 patients in overall denominator
sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #244397 unique records of a HNA/PCSP
sum(hna_pcsp_data$hna == "Y") #133923 unique HNA records
sum(hna_pcsp_data$pcsp == "Y") #110474 unique PCSP records

#counting then excluding records with offered code 04 'Not offered'
sum(hna_pcsp_data$hna == "Y" & hna_pcsp_data$offered_code == "04") #465 not offered HNA records
sum(hna_pcsp_data$pcsp == "Y" & hna_pcsp_data$offered_code == "04") #8561 not offered PCSP records

hna_pcsp_data <- hna_pcsp_data %>%
  filter(offered_code != "04" | is.na(offered_code)) #counting NA or invalid offered_code as offered for now

sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #235371 unique records of an OFFERED HNA/PCSP 
sum(hna_pcsp_data$hna == "Y") #133458 unique OFFERED HNA records
sum(hna_pcsp_data$pcsp == "Y") #101913 unique OFFERED PCSP records


###### REMOVING NON-SUBMITTING TRUSTS #######
#number by trust
by_trust <- left_join((hna_pcsp_data %>%
                         mutate(hna = ifelse(hna == "Y", 1, 0), pcsp = ifelse(pcsp == "Y", 1, 0)) %>%
                         group_by(diag_trust) %>%
                         summarise(hna_count = sum(hna), pcsp_count = sum(pcsp))),
                      (patient_cohort %>%
                         select(c(diag_trust, patientid)) %>%
                         unique() %>%
                         group_by(diag_trust) %>%
                         summarise(pat_count = n())),
                      by = "diag_trust")

by_trust <- by_trust %>%
  mutate(keepforhna = ifelse((hna_count < 10 & pat_count <=100) | (hna_count < 20 & pat_count > 100), "EXCLUDE", "INCLUDE"),
         keepforpcsp = ifelse((pcsp_count < 10 & pat_count <=100) | (pcsp_count < 20 & pat_count > 100), "EXCLUDE", "INCLUDE")) %>%
  select(-c(hna_count, pcsp_count, pat_count))

hna_pcsp_data <- left_join(hna_pcsp_data, by_trust, by = "diag_trust", relationship = "many-to-one") 
hna_pcsp_data <- hna_pcsp_data %>% filter(keepforhna == "INCLUDE" | keepforpcsp == "INCLUDE")

#splitting into HNA/PCSP
hna_data <- hna_pcsp_data %>% filter(hna == "Y" & keepforhna == "INCLUDE") 
pcsp_data <- hna_pcsp_data %>% filter(pcsp == "Y" & keepforpcsp == "INCLUDE")

#patient denominator for HNA/PCSP after excluding trusts
patient_cohort <- left_join(patient_cohort, by_trust, by = "diag_trust", relationship = "many-to-one")
length(which(patient_cohort$keepforhna == "INCLUDE")) #304943 denominator for HNAs
length(which(patient_cohort$keepforpcsp == "INCLUDE")) #297165 denominator for PCSPs


###### PATIENT NUMBERS AND DATASETS #######
#overall numbers of patients with HNAs/PCSPs
length(unique(hna_data$patientid)) #79717 patients with at least one HNA record
length(unique(pcsp_data$patientid)) #71627 patients with at least one PCSP record

ids_hna <- unique(hna_data$patientid)
ids_pcsp <- unique(pcsp_data$patientid)
ids_either <- unique(c(ids_hna, ids_pcsp)) #86054 patients with at least one HNA or PCSP record

only_hna <- setdiff(ids_hna, ids_pcsp) #14427 patients with only an HNA 
only_pcsp <- setdiff(ids_pcsp, ids_hna) #6337 patients with only a PCSP 
both <- intersect(ids_hna, ids_pcsp) #65249 patients with both

#write out record level HNA and PCSP data
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA events raw data RCRD 20240508.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSP events raw data RCRD 20240508.csv")


###### RANKING BY DATE ######
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
write.csv(hna_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA in RCRD ranked by patient 20240408.csv")
write.csv(pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/PCSPs in RCRD ranked by patient 20240408.csv")
write.csv(hna_pcsp_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/HNA and PCSPs in RCRD ranked by patient 20240408.csv")


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
  select(-c(hna, pcsp, offered_code, event_end, keepforpcsp)) %>%
  rename("hna_event_type" = "event_type", "hna_point_of_pathway" = "point_of_pathway", "hna_staff_role" = "staff_role", "hna_date" = "event_date", 
         "hna_rank" = "rank", "hna_time_diag_event" = "time_diag_event")

pcsp_data <- pcsp_data %>% 
  filter(rank == 1) %>%
  select(c(patientid, event_type, point_of_pathway, offered_code, staff_role, event_date, rank, time_diag_event, keepforpcsp)) %>%
  rename("pcsp_event_type" = "event_type", "pcsp_point_of_pathway" = "point_of_pathway", "pcsp_staff_role" = "staff_role", "pcsp_date" = "event_date", 
         "pcsp_rank" = "rank", "pcsp_time_diag_event" = "time_diag_event")

both_patient_level <- inner_join(hna_data, pcsp_data, by = "patientid") 

hna_only_patient_level <- anti_join(hna_data, both_patient_level, by = "patientid") %>%
  mutate(pcsp_event_type = NA, pcsp_point_of_pathway = NA, pcsp_staff_role = NA, pcsp_date = NA)

pcsp_only_patient_level <- anti_join(pcsp_data, both_patient_level, by = "patientid") %>%
  mutate(hna_event_type = NA, hna_point_of_pathway = NA, hna_staff_role = NA, hna_date = NA)

no_events <- patient_cohort %>% 
  left_join(., hna_only_patient_level %>% select(c(patientid, hna_event_type)), by = "patientid") %>% 
  left_join(., pcsp_only_patient_level %>% select(c(patientid, pcsp_event_type)), by = "patientid") %>%
  left_join(., both_patient_level %>% select(c(patientid, pcsp_event_type)), by = "patientid") %>%
  filter(is.na(hna_event_type) & is.na(pcsp_event_type.x) & is.na(pcsp_event_type.y)) %>%
  mutate(hna_event_type = NA, hna_point_of_pathway = NA, hna_staff_role = NA, hna_date = NA, hna_rank = NA, hna_time_diag_event = NA,
         pcsp_event_type = NA, pcsp_point_of_pathway = NA, pcsp_staff_role = NA, pcsp_date = NA, pcsp_rank = NA, pcsp_time_diag_event = NA) %>%
  select(-c(pcsp_event_type.x, pcsp_event_type.y))

#combining all groups of patients into one patient-level data frame 
patient_level_data <- rbind(both_patient_level, hna_only_patient_level, pcsp_only_patient_level, no_events) 

#adding in HNA and PCSP counts for each patient
patient_level_data <- left_join(patient_level_data, hna_count, by = "patientid")
patient_level_data <- left_join(patient_level_data, pcsp_count, by = "patientid")

patient_level_data <- patient_level_data %>%
  mutate(hna_count = ifelse(is.na(hna_count), 0, hna_count),
         pcsp_count = ifelse(is.na(pcsp_count), 0, pcsp_count)) %>%
  select(-rank, hna_rank, pcsp_rank) %>%
  mutate(hna_status = ifelse(hna_count > 0, "Has HNA", "No HNA"),
         pcsp_status = ifelse(pcsp_count > 0, "Has PCSP", "No PCSP"))

#write out patient-level HNA and PCSP data
write.csv(patient_level_data, "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level RCRD HNA and PCSP data 20240508.csv")

