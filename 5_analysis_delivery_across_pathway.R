################ HNA/PCSP coverage analysis - delivery across the pathway ###################

#Analysis for proposal aim 4 - delivery across the pathway 
#Used in Section 7 of initial results slide deck

#Created July 2024 by Lizzie Augarde 
##########################################################################################

library(scales)


######### COSD POINT OF PATHWAY CODE #########
pathway_data <- hna_pcsp_data |>
  filter(offered_code != "04") |> #excluding records which are "Not offered"
  filter(point_of_pathway != "" & point_of_pathway != "09") |> #dealing with weird codes
  mutate(point_of_pathway = ifelse(point_of_pathway %in% c("98", "99"), "97", point_of_pathway)) |>
  group_by(event_type, point_of_pathway) |>
  summarise(count = n()) |>
  group_by(event_type) |>
  mutate(total_count = sum(count), #count of HNAs/PCSPs by point in pathway
         percent = percent((count/total_count), accuracy = 0.1),
         percent_graph = (count/total_count)*100) |>
  ungroup() |>
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"),
         point_of_pathway_desc = case_when(point_of_pathway == "01" ~ "Initial cancer diagnosis", 
                                           point_of_pathway == "02" ~ "Start of treatment", 
                                           point_of_pathway == "03" ~ "During treatment", 
                                           point_of_pathway == "04" ~ "End of treatment", 
                                           point_of_pathway == "05" ~ "Diagnosis of recurrence", 
                                           point_of_pathway == "06" ~ "Transition to palliative care",
                                           point_of_pathway == "07" ~ "Prehabilitation",
                                           TRUE ~ "Other")) |>
  mutate(point_of_pathway_desc = factor(point_of_pathway_desc, levels = c("Initial cancer diagnosis",
                                                                          "Start of treatment",
                                                                          "During treatment",
                                                                          "End of treatment",
                                                                          "Diagnosis of recurrence", 
                                                                          "Prehabilitation",
                                                                          "Transition to palliative care", "Other")))

#bar plot of HNAs/PCSPs by point in pathway code 
pathway_graph <- ggplot(pathway_data, aes(x = point_of_pathway_desc, y = percent_graph, group = event_type)) + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = event_type)) + 
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  labs(x = "Point of pathway", y = "Proportion of events", fill = "") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))


########## RELATIONSHIP BETWEEN DIAGNOSIS DATE AND FIRST HNA/PCSP #########
patient_level_data <- patient_level_data |>
  mutate(hna_time_diag_event = as.numeric(hna_time_diag_event), #time between diagnosis and first HNA/PCSP
         pcsp_time_diag_event = as.numeric(pcsp_time_diag_event))

#cumulative histogram of number of patients receiving first HNA by different times post-diagnosis
hna_hist <- ggplot(patient_level_data, aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  scale_y_continuous(limits = c(0,80000),
                     breaks = c(20000,40000,60000,80000), 
                     labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730), #adding in breaks for time periods of interest
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#cumulative histogram of number of patients receiving first PCSP by different times post-diagnosis
pcsp_hist <- ggplot(patient_level_data, aes(x = pcsp_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  scale_y_continuous(limits = c(0,80000),
                     breaks = c(20000,40000,60000,80000), 
                     labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730), #adding in breaks for time periods of interest
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first PCSP", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#median time between diagnosis and first HNA/PCSP
median_hna_diag_time <- median(patient_level_data$hna_time_diag_event, na.rm=TRUE)
median_pcsp_diag_time <- median(patient_level_data$pcsp_time_diag_event, na.rm=TRUE)


########## TIME IN WHICH 90% OF PATIENTS WITH AN OFFER HAVE THEIR OFFER #########
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


########## MULTIPLE HNAS/PCSPs########## 
patient_level_data <- patient_level_data |>
  mutate(hna_count = as.character(hna_count),
         hna_count = ifelse(!hna_count %in% c("0", "1", "2"), "3 or more", hna_count),
         pcsp_count = as.character(pcsp_count),
         pcsp_count = ifelse(!pcsp_count %in% c("0", "1", "2"), "3 or more", pcsp_count))

#number and % of patients with 0/1/2/3+ HNAs
hna_count_summary <- patient_level_data |>
  group_by(hna_count) |>
  summarise(number_patients_hna = n()) |>
  ungroup() |>
  mutate(percent_patients_hna = percent((number_patients_hna/sum(number_patients_hna)), accuracy = 0.1)) 

#number and % of patients with 0/1/2/3+ PCSPs
pcsp_count_summary <- patient_level_data |>
  group_by(pcsp_count) |>
  summarise(number_patients_pcsp = n()) |>
  ungroup() |>
  mutate(percent_patients_pcsp = percent((number_patients_pcsp/sum(number_patients_pcsp)), accuracy = 0.1)) 


###### TIME BETWEEN MULTIPLE HNAs ######
multi_hnas <- hna_pcsp_data |>
  filter(event_type == 20) |>
  group_by(patientid) |> 
  summarise(time_difference = as.numeric(diff(sort(event_date))),
            .groups = "drop")

#median time between HNAs where patients have more than 1
median_multi_hnas <- median(multi_hnas$time_difference)

#% of HNAs for the same patient which occur within 1 week of each other
multi_hnas_within_week <- (length(which(multi_hnas$time_difference <8))) / nrow(multi_hnas) 


###### TIME BETWEEN MULTIPLE PCSPs ######
multi_pcsps <- hna_pcsp_data |>
  filter(event_type == 24) |>
  group_by(patientid) |> 
  summarise(time_difference = as.numeric(diff(sort(event_date))),
            .groups = "drop")

#median time between PCSPs where patients have more than 1
median_multi_pcsps <- median(multi_pcsps$time_difference)

#% of PCSPs for the same patient which occur within 1 week of each other
multi_pcsps_within_week <- (length(which(multi_pcsps$time_difference <8))) / nrow(multi_pcsps) 


####### SURVIVAL OVER TIME - HNAs ######
patient_level_data <- patient_level_data |> 
  ungroup() |>
  mutate(diag_death_diff = difftime(as.Date(deathdatebest), as.Date(diagnosisdatebest), units = "days")) |>
  mutate(surv_status_4wks = ifelse(diag_death_diff <28, "Died", "Survived"), #creating survival variables
         surv_status_6wks = ifelse(diag_death_diff <42, "Died", "Survived"),
         surv_status_8wks = ifelse(diag_death_diff <56, "Died", "Survived"),
         surv_status_12wks = ifelse(diag_death_diff <84, "Died", "Survived")) |>
  mutate(across(starts_with("surv_status"), 
                function(x) ifelse(is.na(x), "Survived", x))) 

#number of patients surviving to each time period and number with a HNA
survival_hna <- patient_level_data |> 
  filter(keepforhna == "INCLUDE") |>
  pivot_longer(cols = starts_with("surv_status"), names_to = "time_period", values_to = "surv_status") |>
  group_by(time_period, surv_status) |>
  summarize(total = n(), 
            has_hna = sum(hna_status == "Has HNA")) 

#number of patients dying within each time period and number with a HNA
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

#combining into data frame of % surviving each time period/dying within each time period with a HNA
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

#bar graph of % with HNA by survival status in each time period
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
#number of patients surviving to each time period and number with a PCSP
survival_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  pivot_longer(cols = starts_with("surv_status"), names_to = "time_period", values_to = "surv_status") |>
  group_by(time_period, surv_status) |>
  summarize(total = n(), 
            has_pcsp = sum(pcsp_status == "Has PCSP")) 

#number of patients dying within each time period and number with a PCSP
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

#combining into data frame of % surviving each time period/dying within each time period with a PCSP
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

#bar graph of % with PCSP by survival status in each time period
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




