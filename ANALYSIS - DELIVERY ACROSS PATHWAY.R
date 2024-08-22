################ HNA/PCSP coverage analysis - delivery across the pathway ###################

#Analysis for proposal aim XXXX - delivery across the pathway in relation to key dates

#Created July 2024 by Lizzie Augarde 
#Change log:
##########################################################################################

library(scales)


####### SURVIVAL OVER TIME - HNAs ######
patient_level_data <- patient_level_data |> 
  ungroup() |>
  mutate(diag_death_diff = difftime(as.Date(deathdatebest), as.Date(diagnosisdatebest), units = "days")) |>
  mutate(surv_status_4wks = ifelse(diag_death_diff <28, "Died", "Survived"),
         surv_status_6wks = ifelse(diag_death_diff <42, "Died", "Survived"),
         surv_status_8wks = ifelse(diag_death_diff <56, "Died", "Survived"),
         surv_status_12wks = ifelse(diag_death_diff <84, "Died", "Survived")) |>
  mutate(across(starts_with("surv_status"), 
                function(x) ifelse(is.na(x), "Survived", x))) 

survival_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  pivot_longer(cols = starts_with("surv_status"), names_to = "time_period", values_to = "surv_status") |>
  group_by(time_period, surv_status) |>
  summarize(total = n(), 
            has_hna = sum(hna_status == "Has HNA")) 

deaths_hna <- survival_hna |>
  filter(surv_status == "Died") |>
  pivot_wider(names_from = "time_period", values_from = c("total", "has_hna")) |>
  mutate(total_surv_status_12wks = total_surv_status_12wks - total_surv_status_8wks,
         total_surv_status_8wks = total_surv_status_8wks - total_surv_status_6wks,
         total_surv_status_6wks = total_surv_status_6wks - total_surv_status_4wks,
         has_hna_surv_status_12wks = has_hna_surv_status_12wks - has_hna_surv_status_8wks,
         has_hna_surv_status_8wks = has_hna_surv_status_8wks - has_hna_surv_status_6wks,
         has_hna_surv_status_6wks = has_hna_surv_status_6wks - has_hna_surv_status_4wks) |>
  pivot_longer(cols = starts_with("total_"), names_to = "time_period", names_prefix = "total_", values_to = "total") |>
  pivot_longer(cols = starts_with("has_hna_"), names_to = "time_period_hna", names_prefix = "has_hna_", values_to = "has_hna") |>
  filter(time_period == time_period_hna) %>%
  select(surv_status, time_period, total, has_hna)

survival_hna <- survival_hna |> filter(surv_status == "Survived")

survival_hna <- rbind(survival_hna, deaths_hna) |>
  mutate(prop_with_hna = round((has_hna/total)*100, 1)) |>
  rowwise() |>
  mutate(ci = list(prop.test(has_hna, total)),
         lower = round(ci$conf.int[1] * 100, 1),
         upper = round(ci$conf.int[2] * 100, 1)) |>
  mutate(surv_status_table = case_when(grepl("4wks", time_period) & surv_status == "Survived" ~ "Survived at least 4 weeks",
                                 grepl("6wks", time_period) & surv_status == "Survived" ~ "Survived at least 6 weeks",
                                 grepl("8wks", time_period) & surv_status == "Survived" ~ "Survived at least 8 weeks",
                                 grepl("12wks", time_period) & surv_status == "Survived" ~ "Survived at least 12 weeks",
                                 grepl("4wks", time_period) & surv_status == "Died" ~ "Died within 4 weeks",
                                 grepl("6wks", time_period) & surv_status == "Died" ~ "Died between 4 and 6 weeks",
                                 grepl("8wks", time_period) & surv_status == "Died" ~ "Died between 6 and 8 weeks",
                                 grepl("12wks", time_period) & surv_status == "Died" ~ "Died between 8 and 12 weeks")) |>
  mutate(time_period = case_when(time_period == "surv_status_4wks" ~ "4 weeks",
                                 time_period == "surv_status_6wks" ~ "6 weeks",
                                 time_period == "surv_status_8wks" ~ "8 weeks",
                                 time_period == "surv_status_12wks" ~ "12 weeks")) |>
  mutate(time_period= factor(time_period, levels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks")))

surv_hna_graph <- ggplot(survival_hna, aes(x = time_period, y = prop_with_hna, 
                                           fill = surv_status)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, 
                position = position_dodge(width = 0.9)) +
  labs(x = "Time period", y = "Proportion of patients offered a HNA", 
       fill = "Survival status") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  theme_minimal() 


####### SURVIVAL OVER TIME - PCSPs ######
survival_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  pivot_longer(cols = starts_with("surv_status"), names_to = "time_period", values_to = "surv_status") |>
  group_by(time_period, surv_status) |>
  summarize(total = n(), 
            has_pcsp = sum(pcsp_status == "Has PCSP")) 

deaths_pcsp <- survival_pcsp |>
  filter(surv_status == "Died") |>
  pivot_wider(names_from = "time_period", values_from = c("total", "has_pcsp")) |>
  mutate(total_surv_status_12wks = total_surv_status_12wks - total_surv_status_8wks,
         total_surv_status_8wks = total_surv_status_8wks - total_surv_status_6wks,
         total_surv_status_6wks = total_surv_status_6wks - total_surv_status_4wks,
         has_pcsp_surv_status_12wks = has_pcsp_surv_status_12wks - has_pcsp_surv_status_8wks,
         has_pcsp_surv_status_8wks = has_pcsp_surv_status_8wks - has_pcsp_surv_status_6wks,
         has_pcsp_surv_status_6wks = has_pcsp_surv_status_6wks - has_pcsp_surv_status_4wks) |>
  pivot_longer(cols = starts_with("total_"), names_to = "time_period", names_prefix = "total_", values_to = "total") |>
  pivot_longer(cols = starts_with("has_pcsp_"), names_to = "time_period_pcsp", names_prefix = "has_pcsp_", values_to = "has_pcsp") |>
  filter(time_period == time_period_pcsp) %>%
  select(surv_status, time_period, total, has_pcsp)

survival_pcsp <- survival_pcsp |> filter(surv_status == "Survived")

survival_pcsp <- rbind(survival_pcsp, deaths_pcsp) |>
  mutate(prop_with_pcsp = round((has_pcsp/total)*100, 1)) |>
  rowwise() |>
  mutate(ci = list(prop.test(has_pcsp, total)),
         lower = round(ci$conf.int[1] * 100, 1),
         upper = round(ci$conf.int[2] * 100, 1)) |>
  mutate(surv_status_table = case_when(grepl("4wks", time_period) & surv_status == "Survived" ~ "Survived at least 4 weeks",
                                       grepl("6wks", time_period) & surv_status == "Survived" ~ "Survived at least 6 weeks",
                                       grepl("8wks", time_period) & surv_status == "Survived" ~ "Survived at least 8 weeks",
                                       grepl("12wks", time_period) & surv_status == "Survived" ~ "Survived at least 12 weeks",
                                       grepl("4wks", time_period) & surv_status == "Died" ~ "Died within 4 weeks",
                                       grepl("6wks", time_period) & surv_status == "Died" ~ "Died between 4 and 6 weeks",
                                       grepl("8wks", time_period) & surv_status == "Died" ~ "Died between 6 and 8 weeks",
                                       grepl("12wks", time_period) & surv_status == "Died" ~ "Died between 8 and 12 weeks")) |>
  mutate(time_period = case_when(time_period == "surv_status_4wks" ~ "4 weeks",
                                 time_period == "surv_status_6wks" ~ "6 weeks",
                                 time_period == "surv_status_8wks" ~ "8 weeks",
                                 time_period == "surv_status_12wks" ~ "12 weeks")) |>
  mutate(time_period= factor(time_period, levels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks")))

surv_pcsp_graph <- ggplot(survival_pcsp, aes(x = time_period, y = prop_with_pcsp, 
                                             fill = surv_status)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, 
                position = position_dodge(width = 0.9)) +
  labs(x = "Time period", y = "Proportion of patients offered a PCSP", 
       fill = "Survival status") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  theme_minimal() 


####### DAYS BETWEEN DIAGNOSIS AND FIRST HNA/PCSP ######
####### RELATIONSHIP TO DIAGNOSIS DATE ######
patient_level_data <- patient_level_data |>
  mutate(hna_time_diag_event = as.numeric(hna_time_diag_event),
         pcsp_time_diag_event = as.numeric(pcsp_time_diag_event))

hna_hist <- ggplot(patient_level_data, aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  scale_y_continuous(limits = c(0,80000),
                     breaks = c(20000,40000,60000,80000), 
                     labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

median_hna_diag_time <- median(patient_level_data$hna_time_diag_event, na.rm=TRUE)

pcsp_hist <- ggplot(patient_level_data, aes(x = pcsp_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  scale_y_continuous(limits = c(0,80000),
                     breaks = c(20000,40000,60000,80000), 
                     labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first PCSP", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

median_pcsp_diag_time <- median(patient_level_data$pcsp_time_diag_event, na.rm=TRUE)


####### TIME IN WHICH 90% OF PATIENTS WITH AN OFFER HAVE THEIR OFFER ######
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
