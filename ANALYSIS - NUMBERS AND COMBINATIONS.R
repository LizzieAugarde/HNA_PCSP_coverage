################ HNA/PCSP coverage analysis - numbers and combos of HNAS/PCSPs ###################

#Analysis for proposal aim XXXX - people have more than 1 HNA/PCSP and people having both

#Created April 2024 by Lizzie Augarde 
#Change log:
#27/06/2024 reformatted percentage columns
##########################################################################################

library(scales)

###### MULTIPLE HNAS/PCSPs (POST-DEDUPLICATION) ######
patient_level_data <- patient_level_data |>
    mutate(hna_count = as.character(hna_count),
           hna_count = ifelse(!hna_count %in% c("0", "1", "2"), "3 or more", hna_count),
           pcsp_count = as.character(pcsp_count),
           pcsp_count = ifelse(!pcsp_count %in% c("0", "1", "2"), "3 or more", pcsp_count))

hna_count_summary <- patient_level_data |>
  group_by(hna_count) |>
  summarise(number_patients_hna = n()) |>
  ungroup() |>
  mutate(percent_patients_hna = percent((number_patients_hna/sum(number_patients_hna)), accuracy = 0.1)) 

pcsp_count_summary <- patient_level_data |>
  group_by(pcsp_count) |>
  summarise(number_patients_pcsp = n()) |>
  ungroup() |>
  mutate(percent_patients_pcsp = percent((number_patients_pcsp/sum(number_patients_pcsp)), accuracy = 0.1)) 


###### TIME BETWEEN HNA AND PCSP ###### THIS NEEDS REDOING TO ACCOUNT FOR PREVIOUS DIAGNOSES 
time_between <- patient_level_data |> 
  ungroup() |>
  filter(hna_status == "Has HNA", pcsp_status == "Has PCSP") |>
  mutate(hna_pcsp_diff = difftime(as.Date(pcsp_date), as.Date(hna_date), units = "days")) |>
  mutate(hna_pcsp_diff = as.numeric(hna_pcsp_diff)) |>
  filter(hna_pcsp_diff > -1) |> 
  group_by(hna_pcsp_diff) |>
  summarise(count = n())

median_time_between <- median(time_between$hna_pcsp_diff)



###### TIME BETWEEN MULTIPLE HNAs ######
multi_hnas <- hna_pcsp_data |>
  filter(event_type == 20) |>
  group_by(patientid) |> 
  summarise(time_difference = as.numeric(diff(sort(event_date))),
            .groups = "drop")

median_multi_hnas <- median(multi_hnas$time_difference)
multi_hnas_within_week <- (length(which(multi_hnas$time_difference <8))) / nrow(multi_hnas) #HNAs within 1 week of each other

###### TIME BETWEEN MULTIPLE PCSPs ######
multi_pcsps <- hna_pcsp_data |>
  filter(event_type == 24) |>
  group_by(patientid) |> 
  summarise(time_difference = as.numeric(diff(sort(event_date))),
            .groups = "drop")

median_multi_pcsps <- median(multi_pcsps$time_difference)
multi_pcsps_within_week <- (length(which(multi_pcsps$time_difference <8))) / nrow(multi_pcsps) 



