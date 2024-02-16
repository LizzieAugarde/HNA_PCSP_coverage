################ HNA/PCSP coverage analysis - quality and completeness ###################

#Analysis for proposal aim 2 - the extent to which HNA and PCSP records in RCRD include all 
#the COSD data items, inaccurate records, PCSPs without previous HNA

#Created February 2024 by Lizzie Augarde 
#Change log:
##########################################################################################


###### Completeness of COSD data items ######
pathway_var <- hna_pcsp_data %>%
  select(-rank) %>%
  group_by(event_type, point_of_pathway) %>%
  summarise(count = n())

staff_role_var <- hna_pcsp_data %>%
  select(-rank) %>%
  group_by(event_type, staff_role) %>%
  summarise(count = n())


###### Records occuring after death date ######
hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 20) %>%
  summarise(count = sum(deathdatebest < event_date, na.rm = TRUE)) 

hna_pcsp_data %>%
  filter(!is.na(deathdatebest), event_type == 24) %>%
  summarise(count = sum(deathdatebest < event_date, na.rm = TRUE)) 

