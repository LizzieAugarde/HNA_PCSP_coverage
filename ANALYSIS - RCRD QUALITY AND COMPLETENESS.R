################ HNA/PCSP coverage analysis - quality and completeness ###################

#Analysis for proposal aim 2 - the extent to which HNA and PCSP records in RCRD include all 
#the COSD data items, inaccurate records, PCSPs without previous HNA

#Created February 2024 by Lizzie Augarde 
#Change log:
#15/03/2024 removed section looking at PCSPs before HNAs due to previous diagnoses. Changed 
#approach to events before diagnosis
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


###### Records occuring after death date ######
hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 20) %>%
  summarise(count = sum(deathdatebest < event_date, na.rm = TRUE)) %>%
  filter(count == 1) %>%
  summarise(total_count = n())

hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 24) %>%
  summarise(count = sum(deathdatebest < event_date, na.rm = TRUE)) %>%
  filter(count == 1) %>%
  summarise(total_count = n())


###### Records occurring before diagnosis date ######
prev_diags <- dbGetQueryOracle(casref01, "select p.patientid
                                                from analysiselizabethaugarde.hna_pcsp_2021diags@cas2401 p 
                                                left join av2021.at_tumour_england@casref01 t on p.patientid = t.patientid
                                                where t.diagnosisyear < 2021", rowlimit = NA)

prev_diags <- prev_diags %>% clean_names()

hna_pcsp_data_prevs <- left_join(hna_pcsp_data, prev_diags, by = "patientid", relationship = "many-to-many", keep = TRUE)

hna_pcsp_data_prevs <- hna_pcsp_data_prevs %>%
  filter(is.na(patientid.y)) 
  
length(unique(hna_pcsp_data_prevs$patientid.x))

hna_pcsp_data_prevs %>%
  filter(!is.na(diagnosisdatebest), event_type == 20) %>%
  summarise(count = sum(event_date < diagnosisdatebest, na.rm = TRUE)) %>%
  filter(count == 1) %>%
  summarise(total_count = n()) #none

hna_pcsp_data_prevs %>%
  filter(!is.na(diagnosisdatebest), event_type == 24) %>%
  summarise(count = sum(event_date < diagnosisdatebest, na.rm = TRUE)) %>%
  filter(count == 1) %>%
  summarise(total_count = n()) #none





