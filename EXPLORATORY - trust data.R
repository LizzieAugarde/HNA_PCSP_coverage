################ HNA/PCSP coverage analysis - trusts ###################

#Investigating records by trust for considering trust exclusions

#Created August 2024 by Lizzie Augarde 
############################################################### 

hna_pcsp_data_trust <- hna_pcsp_data |>
  mutate(event_year = year(event_date), 
         event_month = month(event_date),
         event_year_month = paste0(event_month, "-", event_year)) |>
  group_by(diag_trust, event_year_month, event_type) |>
  summarise(count = n())
