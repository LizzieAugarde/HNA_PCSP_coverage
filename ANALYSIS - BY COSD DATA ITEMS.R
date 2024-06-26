################ HNA/PCSP coverage analysis - staff role and pathway data items ###################

#Analysis for proposal aim XXXX - the proportion of HNAs and PCSPs delivered by staff in 
#different roles and across the pathway based on the COSD pathway data item

#Created April 2024 by Lizzie Augarde 
#Change log:
##########################################################################################

library(plotly)

###### Staff role ######
staff_role_data <- hna_pcsp_data %>%
  filter(offered_code != "04") %>% #excluding records which are "Not offered"
  filter(staff_role != "" & staff_role != "EH") %>%
  group_by(event_type, staff_role) %>%
  summarise(count = n()) %>%
  group_by(event_type) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup() %>%
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"),
         staff_role = case_when(staff_role == "01" ~ "CNS", 
                                staff_role == "02" ~ "Other nurse", 
                                staff_role == "03" ~ "AHP", 
                                staff_role == "04" ~ "Support worker", 
                                staff_role == "05" ~ "Psychologist/MH professional", 
                                staff_role == "06" ~ "Consultant/medical team",
                                staff_role == "08" ~ "Other",
                                TRUE ~ "Not known"))


###### Point of pathway variable ######
pathway_data <- hna_pcsp_data %>%
  filter(offered_code != "04") %>% #excluding records which are "Not offered"
  filter(point_of_pathway != "" & point_of_pathway != "09") %>%
  mutate(point_of_pathway = ifelse(point_of_pathway %in% c("98", "99"), "97", point_of_pathway)) %>%
  group_by(event_type, point_of_pathway) %>%
  summarise(count = n()) %>%
  group_by(event_type) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup() %>%
  mutate(event_type = ifelse(event_type == 20, "HNA", "PCSP"),
         point_of_pathway = case_when(point_of_pathway == "01" ~ "Initial cancer diagnosis", 
                                      point_of_pathway == "02" ~ "Start of treatment", 
                                      point_of_pathway == "03" ~ "During treatment", 
                                      point_of_pathway == "04" ~ "End of treatment", 
                                      point_of_pathway == "05" ~ "Diagnosis of recurrence", 
                                      point_of_pathway == "06" ~ "Transition to palliative care",
                                      point_of_pathway == "07" ~ "Prehabilitation",
                                      TRUE ~ "Other"))


ggplot(staff_role_data, aes(x = staff_role, y = percentage, group = event_type)) + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = event_type)) + 
  labs(x = "Staff role", y = "Percentage", fill = "") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(pathway_data, aes(x = point_of_pathway, y = percentage, group = event_type)) + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = event_type)) + 
  labs(x = "Point of pathway", y = "Percentage", fill = "") + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


