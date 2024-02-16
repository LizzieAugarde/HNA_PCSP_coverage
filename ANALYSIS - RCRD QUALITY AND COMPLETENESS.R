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


###### PCSPs without a previous HNA record ######
pcsp_vs_hna_check <- hna_pcsp_data %>%
  filter(rank == 1) %>%
  rank_by_date(., patientid, event_date)

count <- pcsp_vs_hna_check %>% 
  filter(rank == 1 & pcsp == "Y") #total PCSP patients with no previous HNA or an HNA after PCSP

pcsp_vs_hna_check %>%
  group_by(patientid) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  summarise(count = sum(pcsp == "Y", na.rm = TRUE)) #total PCSP patients with no previous HNA

pcsp_vs_hna_check %>%
  group_by(patientid) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  summarise(count = sum(hna == "Y", na.rm = TRUE)) #total PCSP patients with a HNA after PCSP
  

