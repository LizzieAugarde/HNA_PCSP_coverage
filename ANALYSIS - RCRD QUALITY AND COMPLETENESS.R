################ HNA/PCSP coverage analysis - quality and completeness ###################

#Analysis for proposal aim 2 - the extent to which HNA and PCSP records in RCRD include all 
#the COSD data items, inaccurate records, PCSPs without previous HNA

#Created February 2024 by Lizzie Augarde 
#Change log:
#15/03/2024 removed section looking at PCSPs before HNAs due to previous diagnoses. Changed 
#approach to events before diagnosis
#18/03/2024 removed date before diagnosis check section - excluding in data prep
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

