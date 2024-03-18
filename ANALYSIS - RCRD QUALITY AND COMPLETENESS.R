################ HNA/PCSP coverage analysis - quality and completeness ###################

#Analysis for proposal aim 2 - the extent to which HNA and PCSP records in RCRD include all 
#the COSD data items, inaccurate records, PCSPs without previous HNA

#Created February 2024 by Lizzie Augarde 
#Change log:
#15/03/2024 removed section looking at PCSPs before HNAs due to previous diagnoses. Changed 
#approach to events before diagnosis
#18/03/2024 tidied date before diagnosis check section
##########################################################################################


###### Completeness of COSD data items ######
pathway_var <- hna_pcsp_data %>%
  select(-rank) %>%
  group_by(event_type, point_of_pathway) %>%
  summarise(count = n())

pathway_var %>% filter(point_of_pathway == "" | point_of_pathway == 99) #missing or invalid

staff_role_var <- hna_pcsp_data %>%
  select(-rank) %>%
  group_by(event_type, staff_role) %>%
  summarise(count = n())

staff_role_var %>% filter(staff_role == "" | staff_role == "EH") #missing or invalid


###### Records occurring after death date ######
death_check_hna <- hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 20) %>%
  mutate(death_event_diff = difftime(deathdatebest, event_date, unit = "days")) %>%
  filter(death_event_diff < 0)
  
death_check_pcsp <- hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 24) %>%
  mutate(death_event_diff = difftime(deathdatebest, event_date, unit = "days")) %>%
  filter(death_event_diff < 0)


###### Records occurring before diagnosis date for those with no previous diagnosis ######
prev_diags <- dbGetQueryOracle(casref01, "select p.patientid
                                                from analysiselizabethaugarde.hna_pcsp_2021diags@cas2401 p 
                                                left join av2021.at_tumour_england@casref01 t on p.patientid = t.patientid
                                                where t.diagnosisyear < 2021", rowlimit = NA)

prev_diags <- prev_diags %>% clean_names()

hna_pcsp_data_prevs <- left_join(hna_pcsp_data, prev_diags, by = "patientid", relationship = "many-to-many", keep = TRUE)

hna_pcsp_data_prevs <- hna_pcsp_data_prevs %>%
  filter(is.na(patientid.y))

sum(hna_pcsp_data_prevs$hna == "Y") #112110
sum(hna_pcsp_data_prevs$pcsp == "Y") #93207

diagdate_check_hna <- hna_pcsp_data_prevs %>%
  filter(!is.na(diagnosisdatebest), event_type == 20) %>%
  mutate(diag_event_diff = difftime(event_date, diagnosisdatebest, unit = "days")) %>%
  filter(diag_event_diff < 0) #4702

diagdate_check_pcsp <- hna_pcsp_data_prevs %>%
  filter(!is.na(diagnosisdatebest), event_type == 24) %>%
  mutate(diag_event_diff = difftime(event_date, diagnosisdatebest, unit = "days")) %>%
  filter(diag_event_diff < 0) #5572


