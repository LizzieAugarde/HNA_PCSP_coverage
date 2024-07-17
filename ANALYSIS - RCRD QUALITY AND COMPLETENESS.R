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
pathway_var_comp <- hna_pcsp_data %>%
  filter(point_of_pathway == "" | point_of_pathway == 99) %>% #missing or invalid
  group_by(event_type) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  mutate(measure = "Missing point of pathway code") %>%
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"))

staff_role_var_comp <- hna_pcsp_data %>%
  filter(staff_role == "" | staff_role == "EH") %>% #missing or invalid
  group_by(event_type) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  mutate(measure = "Missing staff role code") %>%
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"))


###### Records occurring after death date ######
death_check_hna <- hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 20) %>%
  mutate(death_event_diff = difftime(deathdatebest, event_date, unit = "days")) %>%
  filter(death_event_diff < 0)
  
death_check_pcsp <- hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 24) %>%
  mutate(death_event_diff = difftime(deathdatebest, event_date, unit = "days")) %>%
  filter(death_event_diff < 0)

death_check_hna <- length(death_check_hna$patientid)
death_check_pcsp <- length(death_check_pcsp$patientid)

death_check <- data.frame(event_type = c("HNA", "PCSP"), 
                          count = c(death_check_hna, death_check_pcsp)) %>%
  mutate(measure = "Event date after death date")


###### Completeness table ######
cosd_completeness <- rbind(pathway_var_comp, staff_role_var_comp, death_check)

cosd_completeness <- cosd_completeness %>%
  mutate(denom = ifelse(event_type == "HNA", unique_hnas, unique_pcsps),
         percent = percent((count/denom), accuracy = 0.1)) %>%
  select(-denom)

cosd_completeness <- cosd_completeness[,c(3,1,2,4)]
  
  
  
  
  
  
  
  
  
  
