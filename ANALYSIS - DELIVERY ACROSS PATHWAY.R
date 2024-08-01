################ HNA/PCSP coverage analysis - delivery across the pathway ###################

#Analysis for proposal aim XXXX - delivery across the pathway in relation to key dates

#Created July 2024 by Lizzie Augarde 
#Change log:
##########################################################################################


####### DEATHS ######
patient_level_data <- patient_level_data |> 
  mutate(diag_death_diff = difftime(as.Date(deathdatebest), as.Date(diagnosisdatebest), units = "days")) |>
  mutate(death_status = ifelse(diag_death_diff <56, "Died within 8 weeks", "Survived at least 8 weeks"))

events_after_8weeks_deaths <- patient_level_data |>
  filter(death_status == "Died within 8 weeks") |>
  group_by(time_diag_event) |>
  summarise(count_patients = n())
