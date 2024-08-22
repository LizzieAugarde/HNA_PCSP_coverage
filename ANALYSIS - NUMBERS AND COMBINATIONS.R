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

###### RELATIONSHIP BETWEEN HNA AND PCSP STATUS ######
offered_code_matrix <- patient_level_data |> 
  select(patientid, hna_offered_code, pcsp_offered_code) |>
  filter(!is.na(hna_offered_code), !is.na(pcsp_offered_code)) |>
  mutate(hna_offered_code_desc = case_when(hna_offered_code == "01" ~ "Offered and undecided", 
                                           hna_offered_code == "02" ~ "Offered and declined", 
                                           hna_offered_code == "03" ~ "Offered and accepted", 
                                           hna_offered_code == "05" ~ "Offered but patient unable to complete", 
                                           TRUE ~ "Not known"), 
         pcsp_offered_code_desc = case_when(pcsp_offered_code == "01" ~ "Offered and undecided", 
                                            pcsp_offered_code == "02" ~ "Offered and declined", 
                                            pcsp_offered_code == "03" ~ "Offered and accepted", 
                                            pcsp_offered_code == "05" ~ "Offered but patient unable to complete", 
                                            pcsp_offered_code == "06" ~ "Not required (no concerns from HNA)",
                                            TRUE ~ "Not known"))  |>
  count(hna_offered_code_desc, pcsp_offered_code_desc)  |>
  mutate(proportion = (n/sum(n))*100)

offered_code_matrix_plot <- ggplot(offered_code_matrix, aes(x = pcsp_offered_code_desc, 
                                                            y = hna_offered_code_desc, 
                                                            fill = proportion)) +
  geom_tile(color = "black") +
  geom_text(aes(label = paste0(round(proportion, digits = 1), "%"))) +
  scale_fill_gradient(low = "#d2f7dd", high = "#008A26", na.value = "white", 
                      name = "Proportion of patients") +
  labs(x = "PCSP status", y = "HNA status") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

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



