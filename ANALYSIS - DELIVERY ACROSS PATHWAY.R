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
  mutate(death_status_4wks = ifelse(diag_death_diff <28, "Died", "Survived"),
         death_status_6wks = ifelse(diag_death_diff <42, "Died", "Survived"),
         death_status_8wks = ifelse(diag_death_diff <56, "Died", "Survived"),
         death_status_12wks = ifelse(diag_death_diff <84, "Died", "Survived")) |>
  mutate(across(starts_with("death_status"), 
                function(x) ifelse(is.na(x), "Survived", x)))

deaths_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  pivot_longer(cols = starts_with("death_status"), names_to = "time_period", values_to = "death_status") |>
  group_by(time_period, death_status) |>
  summarize(total = n(), 
            number_with_hna = sum(hna_status == "Has HNA"),
            prop_with_hna = round(mean(hna_status == "Has HNA")*100, 1), 
            lower = round(prop.test(number_with_hna, total)$conf.int[1]*100, 1),
            upper = round(prop.test(number_with_hna, total)$conf.int[2]*100, 1)) |>
  mutate(time_period = case_when(time_period == "death_status_4wks" ~ "4 weeks",
                                 time_period == "death_status_6wks" ~ "6 weeks",
                                 time_period == "death_status_8wks" ~ "8 weeks",
                                 time_period == "death_status_12wks" ~ "12 weeks")) |>
  mutate(time_period = factor(time_period, levels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks")))

deaths_hna_graph <- ggplot(deaths_hna, aes(x = time_period, y = prop_with_hna, 
                                           fill = death_status)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, 
                position = position_dodge(width = 0.9)) +
  labs(x = "Time period", y = "Proportion of patients offered a HNA", 
       fill = "Survival status") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  theme_minimal() 

deaths_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  pivot_longer(cols = starts_with("death_status"), names_to = "time_period", values_to = "death_status") |>
  group_by(time_period, death_status) |>
  summarize(total = n(), 
            number_with_pcsp = sum(pcsp_status == "Has PCSP"),
            prop_with_pcsp = round(mean(pcsp_status == "Has PCSP")*100, 1), 
            lower = round(prop.test(number_with_pcsp, total)$conf.int[1]*100, 1),
            upper = round(prop.test(number_with_pcsp, total)$conf.int[2]*100, 1)) |>
  mutate(time_period = case_when(time_period == "death_status_4wks" ~ "4 weeks",
                                 time_period == "death_status_6wks" ~ "6 weeks",
                                 time_period == "death_status_8wks" ~ "8 weeks",
                                 time_period == "death_status_12wks" ~ "12 weeks")) |>
  mutate(time_period= factor(time_period, levels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks")))

deaths_pcsp_graph <- ggplot(deaths_pcsp, aes(x = time_period, y = prop_with_pcsp, 
                                           fill = death_status)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, 
                position = position_dodge(width = 0.9)) +
  labs(x = "Time period", y = "Proportion of patients offered a PCSP", 
       fill = "Survival status") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  theme_minimal() 


####### RELATIONSHIP TO DIAGNOSIS DATE ######
patient_level_data <- patient_level_data |>
  mutate(hna_time_diag_event = as.numeric(hna_time_diag_event),
         pcsp_time_diag_event = as.numeric(pcsp_time_diag_event))

hna_hist <- ggplot(patient_level_data, aes(x = hna_time_diag_event)) +
  geom_histogram(aes(y = cumsum(..count..)), 
                 binwidth = 10, fill = "#008A26", color = "black") +
  scale_y_continuous(limits = c(0,81000),
                     breaks = c(20000,40000, 60000, 80000), 
                     labels = label_comma()) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

pcsp_hist <- ggplot(patient_level_data, aes(x = pcsp_time_diag_event)) +
  geom_histogram(aes(y = cumsum(..count..)), 
                     binwidth = 10, fill = "#008A26", color = "black") +
  scale_y_continuous(limits = c(0,81000),
                     breaks = c(20000,40000, 60000, 80000), 
                     labels = label_comma()) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first PCSP", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#median number of days between diagnosis and first HNA/PCSP
median_hna_diag_time <- median(patient_level_data$hna_time_diag_event, na.rm=TRUE)
median_pcsp_diag_time <- median(patient_level_data$pcsp_time_diag_event, na.rm=TRUE)


#time in which 90% of patients with an offer have their HNA/PCSP 
hnas_90_percent <- patient_level_data |>
  filter(hna_status == "Has HNA", keepforhna == "INCLUDE") |>
  arrange(hna_time_diag_event) |> 
  mutate(cumulative_proportion = cumsum(rep(1/n(), n()))) |> 
  filter(cumulative_proportion >= 0.90) |> 
  slice(1) |> 
  pull(hna_time_diag_event)

pcsps_90_percent <- patient_level_data |>
  filter(pcsp_status == "Has PCSP", keepforpcsp == "INCLUDE") |>
  arrange(pcsp_time_diag_event) |> 
  mutate(cumulative_proportion = cumsum(rep(1/n(), n()))) |> 
  filter(cumulative_proportion >= 0.90) |> 
  slice(1) |> 
  pull(pcsp_time_diag_event)
