################ HNA/PCSP coverage analysis - data prep ###################

#Data preparation and cleaning of COSD Level 3 registrations for 
#2021 linked to RCRD HNA and PCSP records for coverage analysis

#Created November 2023 by Lizzie Augarde 
#Change log:
#15/03/2024 adjusted approach to counting patients
#08/04/2024 adjusted query approach
#25/07/2024 added staging functions and site from AT_SITE_ENGLAND
############################################################### 

#prep
library(NDRSAfunctions)
library(tidyverse)
library(janitor)

cas2406 <- createConnection(port = 1525, sid = "cas2406")

###### EXTRACTING AND CLEANING RAW DATA #######
#patient cohort----
query <- "select * from (
  select pat.patientid,
  tum.tumourid,
  tum.diagnosisdatebest,
  tum.site_icd10r4_o2_3char_from2013,
  tum.stage_pi_detail,
  tum.stage_best_system,
  tum.stage_best,
  tum.gender,
  tum.age,
  tum.ethnicity,
  tum.deathdatebest,
  tum.diag_trust,
  imd.imd19_decile_lsoas,
  site.ndrs_main,
  rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join analysispollyjeffrey.at_site_england@casref01 site on tum.tumourid = site.tumourid --added 25/07/2024
  left join av2021.at_geography_england@casref01 geo on tum.tumourid = geo.tumourid
  left join imd.imd2019_equal_lsoas@casref01 imd on geo.lsoa11_code = imd.lsoa11_code
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  --and tum.age >17 --removed 17/07/2024
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44') 
where rank = 1"

patient_cohort <- dbGetQueryOracle(cas2406, query, rowlimit = NA)

patient_cohort <- patient_cohort %>% 
  clean_names() %>%
  unique()


###### RECODING STAGE AND REMOVING STAGE 0 CASES ######
source("stage_functions_tidy.R")

#trans patient fix from 02_data_tidying_tidy.R received from Chloe Bright
patient_cohort <- patient_cohort %>% mutate(
  stage_pi_detail = if_else(
    condition = (
      (
        (gender == 2 & site_icd10r4_o2_3char_from2013 %in% c("C60", "C61", "C62", "C63")) |
          (gender == 1 & site_icd10r4_o2_3char_from2013 %in% c("C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58"))
      ) &
        !(stage_best == "?" | is.na(stage_best))
    ),
    true = "Y",
    false = stage_pi_detail
  )
)

patient_cohort <- stage_table(patient_cohort) 

patient_cohort <- patient_cohort |>
  select(-c(stage_best_system, stage_pi_detail, site_icd10r4_o2_3char_from2013, SUMMARY_STAGE, EARLY_ADVANCED_STAGE)) 

#HNA and PCSP records------
query <- "select * from (
  select pat.patientid,
  rp.nhsnumber as rp_nhsnumber,
  rp.patientid as rp_patientid,       
  tum.tumourid,
  tum.diagnosisdatebest,
  rp.event_type,
  rp.event_property_1,
  rp.event_property_2,
  rp.event_property_3,
  rp.event_date,
  rp.event_end,
  rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join analysisncr.at_rapid_pathway@cas2406 rp on pat.nhsnumber = rp.nhsnumber 
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  --and tum.age >17 ---removed 17/07/2024
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44'
  and rp.event_type in (20, 24))
where rank = 1"

raw_hna_pcsp_data <- dbGetQueryOracle(cas2406, query, rowlimit = NA)

raw_hna_pcsp_data <- raw_hna_pcsp_data %>% 
  clean_names() %>%

  #limiting to records within 2 years of diagnosis 
  mutate(time_diag_event = difftime(as.Date(event_date), as.Date(diagnosisdatebest), units = "days")) %>%
  filter((time_diag_event < 731 & time_diag_event >= 0)) %>%
  
  #marking HNAs and PCSPs
  mutate(hna = case_when(event_type == 20 ~ "Y", TRUE ~ "N"),
         pcsp = case_when(event_type == 24 ~ "Y", TRUE ~ "N")) %>%
  
  #joining to patient cohort to get demographic info 
  left_join(., patient_cohort, by = "patientid") |> 
  select(-c(tumourid.y, diagnosisdatebest.y, rank.x, rank.y)) |>
  rename("diagnosisdatebest" = "diagnosisdatebest.x", "tumourid" = "tumourid.x")

#count of all HNA and PCSP records 
all_records <- sum(raw_hna_pcsp_data$hna == "Y" | raw_hna_pcsp_data$pcsp == "Y") #1895742 records of a HNA/PCSP


###### REMOVING DUPLICATES, INITIAL NUMBERS #######
hna_pcsp_data <- raw_hna_pcsp_data %>%
  unique() %>% #removing duplicate rows
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") 

denom_size <- length(unique(patient_cohort$patientid)) #311627 patients in overall denominator
unique_hnas_pcsps <- sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #245592 unique records of a HNA/PCSP
unique_hnas <-sum(hna_pcsp_data$hna == "Y") #134332 unique HNA records
unique_pcsps <-sum(hna_pcsp_data$pcsp == "Y") #111260 unique PCSP records

#counting then excluding records with offered code 04 'Not offered'
hnas_notoff <- sum(hna_pcsp_data$hna == "Y" & hna_pcsp_data$offered_code == "04") #466 not offered HNA records
pcsps_notoff <- sum(hna_pcsp_data$pcsp == "Y" & hna_pcsp_data$offered_code == "04") #8743 not offered PCSP records

hna_pcsp_data <- hna_pcsp_data %>%
  filter(offered_code != "04" | is.na(offered_code)) #counting NA or invalid offered_code as offered for now

unique_off_hnas_pcsps <- sum(hna_pcsp_data$hna == "Y" | hna_pcsp_data$pcsp == "Y") #236383 unique records of an OFFERED HNA/PCSP 
unique_off_hnas <- sum(hna_pcsp_data$hna == "Y") #133866 unique OFFERED HNA records
unique_off_pcsps <- sum(hna_pcsp_data$pcsp == "Y") #102517 unique OFFERED PCSP records

pcsps_notreq <- sum(hna_pcsp_data$pcsp == "Y" & hna_pcsp_data$offered_code == "06") #40894 not required PCSP records


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
hnas_denom_size <- length(which(patient_cohort$keepforhna == "INCLUDE")) #306399 denominator for HNAs
pcsps_denom_size <- length(which(patient_cohort$keepforpcsp == "INCLUDE")) #300145 denominator for PCSPs


###### PATIENT NUMBERS AND DATASETS #######
#overall numbers of patients with HNAs/PCSPs
pts_with_hna <- length(unique(hna_data$patientid)) #79897 patients with at least one HNA record
pts_with_pcsp <- length(unique(pcsp_data$patientid)) #71920 patients with at least one PCSP record

ids_hna <- unique(hna_data$patientid)
ids_pcsp <- unique(pcsp_data$patientid)
ids_either <- length(unique(c(ids_hna, ids_pcsp))) #86424 patients with at least one HNA or PCSP record

only_hna <- length(setdiff(ids_hna, ids_pcsp)) #14504 patients with only an HNA 
only_pcsp <- length(setdiff(ids_pcsp, ids_hna)) #6527 patients with only a PCSP
both <- length(intersect(ids_hna, ids_pcsp)) #65393 patients with both


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
  select(-c(hna, pcsp, event_end, keepforpcsp)) %>%
  rename("hna_event_type" = "event_type", "hna_point_of_pathway" = "point_of_pathway", "hna_staff_role" = "staff_role", "hna_date" = "event_date", 
         "hna_rank" = "rank", "hna_time_diag_event" = "time_diag_event", "hna_offered_code" = "offered_code")

pcsp_data <- pcsp_data %>% 
  filter(rank == 1) %>%
  select(c(patientid, event_type, offered_code, point_of_pathway, staff_role, event_date, rank, time_diag_event, keepforpcsp, STAGE)) %>%
  rename("pcsp_event_type" = "event_type", "pcsp_point_of_pathway" = "point_of_pathway", "pcsp_staff_role" = "staff_role", "pcsp_date" = "event_date", 
         "pcsp_rank" = "rank", "pcsp_time_diag_event" = "time_diag_event", "pcsp_offered_code" = "offered_code")

both_patient_level <- inner_join(hna_data, pcsp_data, by = "patientid") |> select(-STAGE.y) |> rename("STAGE" = "STAGE.x")

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
#save(patient_level_data, file = "N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/Patient-level RCRD HNA and PCSP data 20240717.RData")

