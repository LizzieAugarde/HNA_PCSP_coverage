################ HNA/PCSP coverage analysis - offered code ###################

#Analysis of offered codes for HNA and PCSP records 
#Used in Section 3 and Section 6 of initial results slide deck

#Created April 2024 by Lizzie Augarde 
##############################################################################

library(scales)

###### OFFERED CODES FOR HNAs AND PCSPs ######
offered_code_data <- hna_pcsp_data |>
  filter(offered_code != "04") |> #excluding records which are "Not offered"
  filter(offered_code != "" & offered_code != "08/03/2023" & offered_code != "16/06/2023") |> #removing weird codes
  mutate(offered_code = ifelse(offered_code == "001", "01", offered_code)) |> 
  group_by(event_type, offered_code) |>
  summarise(count = n()) |>
  group_by(event_type) |>
  mutate(total_count = sum(count), #count and % of HNAs/PCSPs by offered_code
         percent = percent((count/total_count), accuracy = 0.1),
         percent_graph = (count/total_count)*100) |>
  ungroup() |>
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"),
         offered_code_desc = case_when(offered_code == "01" ~ "Offered and undecided", 
                                       offered_code == "02" ~ "Offered and declined", 
                                       offered_code == "03" ~ "Offered and accepted", 
                                       offered_code == "05" ~ "Offered but patient unable to complete", 
                                       offered_code == "06" ~ "Not required (no concerns from HNA)",
                                       TRUE ~ "Not known"))


###### RELATIONSHIP BETWEEN FIRST HNA AND PCSP STATUS ######
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

#matrix plot of % of patients by status of 1st HNA and PCSP
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


