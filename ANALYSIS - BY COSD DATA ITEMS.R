################ HNA/PCSP coverage analysis - staff role and pathway data items ###################

#Analysis for proposal aim XXXX - the proportion of HNAs and PCSPs delivered by staff in 
#different roles and across the pathway based on the COSD pathway data item

#Created April 2024 by Lizzie Augarde 
##########################################################################################

library(scales)

###### Offered code ######
offered_code_data <- hna_pcsp_data |>
  filter(offered_code != "04") |> #excluding records which are "Not offered"
  filter(offered_code != "" & offered_code != "08/03/2023" & offered_code != "16/06/2023") |>
  mutate(offered_code = ifelse(offered_code == "001", "01", offered_code)) |>
  group_by(event_type, offered_code) |>
  summarise(count = n()) |>
  group_by(event_type) |>
  mutate(total_count = sum(count),
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
  

###### Staff role ######
staff_role_data <- hna_pcsp_data |>
  filter(offered_code != "04") |> #excluding records which are "Not offered"
  filter(staff_role != "" & staff_role != "EH") |>
  group_by(event_type, staff_role) |>
  summarise(count = n()) |>
  group_by(event_type) |>
  mutate(total_count = sum(count),
         percent = percent((count/total_count), accuracy = 0.1),
         percent_graph = (count/total_count)*100) |>
  ungroup() |>
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"),
         staff_role_desc = case_when(staff_role == "01" ~ "CNS", 
                                     staff_role == "02" ~ "Other nurse", 
                                     staff_role == "03" ~ "AHP", 
                                     staff_role == "04" ~ "Support worker", 
                                     staff_role == "05" ~ "Psychologist/MH professional", 
                                     staff_role == "06" ~ "Consultant/medical team",
                                     staff_role == "08" ~ "Other",
                                     TRUE ~ "Not known")) 


###### Point of pathway variable ######
pathway_data <- hna_pcsp_data |>
  filter(offered_code != "04") |> #excluding records which are "Not offered"
  filter(point_of_pathway != "" & point_of_pathway != "09") |>
  mutate(point_of_pathway = ifelse(point_of_pathway %in% c("98", "99"), "97", point_of_pathway)) |>
  group_by(event_type, point_of_pathway) |>
  summarise(count = n()) |>
  group_by(event_type) |>
  mutate(total_count = sum(count),
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

staff_role_graph <- ggplot(staff_role_data, aes(x = staff_role_desc, y = percent_graph, group = event_type)) + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = event_type)) + 
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  labs(x = "Staff role", y = "Proportion of events", fill = "") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))

pathway_graph <- ggplot(pathway_data, aes(x = point_of_pathway_desc, y = percent_graph, group = event_type)) + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = event_type)) + 
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  labs(x = "Point of pathway", y = "Proportion of events", fill = "") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))



