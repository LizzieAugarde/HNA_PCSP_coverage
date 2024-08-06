################ HNA/PCSP coverage analysis - delivery across the pathway ###################

#Analysis for proposal aim XXXX - delivery across the pathway in relation to key dates

#Created July 2024 by Lizzie Augarde 
#Change log:
##########################################################################################

library(scales)


####### DEATHS ######
patient_level_data <- patient_level_data |> 
  ungroup() |>
  mutate(diag_death_diff = difftime(as.Date(deathdatebest), as.Date(diagnosisdatebest), units = "days")) |>
  mutate(death_status = ifelse(diag_death_diff <56, "Died within 8 weeks", "Survived at least 8 weeks")) |>
  mutate(death_status = ifelse(is.na(diag_death_diff), "Survived at least 8 weeks", death_status))

deaths_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(death_status, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(death_status) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(hna_status == "Has HNA") |>
  ungroup()

deaths_hna_graph <- ggplot(deaths_hna, aes(x = death_status, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Survival status", y = "Proportion of patients offered a HNA") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() 

deaths_pcsp <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(death_status, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(death_status) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1),
         lower = lapply(number_patients, prop.test, n = sum(number_patients)), 
         upper = round((sapply(lower, function(x) x$conf.int[2]))*100, digits = 1), 
         lower = round((sapply(lower, function(x) x$conf.int[1]))*100, digits = 1),
         lower_table = paste0(lower,"%"), 
         upper_table = paste0(upper, "%")) |>
  filter(pcsp_status == "Has PCSP") |>
  ungroup()

deaths_pcsp_graph <- ggplot(deaths_pcsp, aes(x = death_status, y = percent)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008A26") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Survival status", y = "Proportion of patients offered a PCSP") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal()


####### RELATIONSHIP TO DIAGNOSIS DATE ######
patient_level_data <- patient_level_data |>
  mutate(hna_time_diag_event = as.numeric(hna_time_diag_event),
         pcsp_time_diag_event = as.numeric(pcsp_time_diag_event))

hna_hist <- ggplot(patient_level_data, aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black") +
  scale_y_continuous(limits = c(0,15000),
                     breaks = c(5000,10000,15000), 
                     labels = c("5,000", "10,000", "15,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

median_hna_diag_time <- median(patient_level_data$hna_time_diag_event, na.rm=TRUE)

pcsp_hist <- ggplot(patient_level_data, aes(x = pcsp_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black") +
  scale_y_continuous(limits = c(0,25000),
                     breaks = c(5000,10000,15000,20000,25000), 
                     labels = c("5,000", "10,000","15,000","20,000","25,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first PCSP", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

median_pcsp_diag_time <- median(patient_level_data$pcsp_time_diag_event, na.rm=TRUE)

