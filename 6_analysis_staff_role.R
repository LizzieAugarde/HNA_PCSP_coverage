################ HNA/PCSP coverage analysis - staff role ###################

#Analysis of staff role codes for HNA and PCSP records 
#Used in Section 8 of initial results slide deck

#Created April 2024 by Lizzie Augarde 
##############################################################################

library(scales)

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
         staff_role_desc = case_when(staff_role == "01" ~ "Cancer Nurse Specialist", 
                                     staff_role == "02" ~ "Other nurse", 
                                     staff_role == "03" ~ "Allied Health Professional", 
                                     staff_role == "04" ~ "Support worker", 
                                     staff_role == "05" ~ "Psychologist/Mental health professional", 
                                     staff_role == "06" ~ "Consultant/medical team",
                                     staff_role == "08" ~ "Other",
                                     TRUE ~ "Not known")) 

#bar graph of the % of HNAs/PCSPs by staff member delivering
staff_role_graph <- ggplot(staff_role_data, aes(x = staff_role_desc, y = percent_graph, group = event_type)) + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = event_type)) + 
  scale_fill_manual(values = c("#008A26", "#02D462")) +
  labs(x = "Staff role", y = "Proportion of events", fill = "") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))