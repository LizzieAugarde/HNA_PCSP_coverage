################ HNA/PCSP coverage analysis - geographic variation ###################

#Analysis for proposal aim XXXX - proportion of people being offered an HNA/PCSP by 
#diagnosis trust

#Created July 2024 by Lizzie Augarde 
#Change log:
#################################################################################

####HNA coverage categories-------------------
trust_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(diag_trust, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(diag_trust) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(diag_trust)) |>
  mutate(percent_group = case_when(percent <10 ~ "Less than 10%",
                                   percent >=10 & percent <20 ~ "10% - 20%",
                                   percent >=20 & percent <30 ~ "20% - 30%",
                                   percent >=30 & percent <40 ~ "30% - 40%",
                                   percent >=40 & percent <50 ~ "40% - 50%",
                                   percent >=50 & percent <60 ~ "50% - 60%",
                                   percent >=60 & percent <70 ~ "60% - 70%",
                                   percent >=70 & percent <80 ~ "70% - 80%",
                                   percent >=80 & percent <90 ~ "80% - 90%",
                                   percent >=90 & percent <100 ~ "90% - 100%"))

trust_coverage_hna_cats <- trust_hna |> 
  group_by(percent_group) |>
  summarise(trust_count = n()) |>
  mutate(percent_group = fct_relevel(percent_group, "Less than 10%"))

trust_hna_graph <- ggplot(trust_coverage_hna_cats, aes(x = percent_group, y = trust_count)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  labs(x = "Proportion of patients offered a HNA", y = "Number of trusts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))    


####PCSP coverage categories-------------------
trust_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(diag_trust, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(diag_trust) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(diag_trust)) |>
  mutate(percent_group = case_when(percent <10 ~ "Less than 10%",
                                   percent >=10 & percent <20 ~ "10% - 20%",
                                   percent >=20 & percent <30 ~ "20% - 30%",
                                   percent >=30 & percent <40 ~ "30% - 40%",
                                   percent >=40 & percent <50 ~ "40% - 50%",
                                   percent >=50 & percent <60 ~ "50% - 60%",
                                   percent >=60 & percent <70 ~ "60% - 70%",
                                   percent >=70 & percent <80 ~ "70% - 80%",
                                   percent >=80 & percent <90 ~ "80% - 90%",
                                   percent >=90 & percent <100 ~ "90% - 100%"))

trust_coverage_pcsp_cats <- trust_pcsp |> 
  group_by(percent_group) |>
  summarise(trust_count = n()) |>
  mutate(percent_group = fct_relevel(percent_group, "Less than 10%"))

trust_pcsp_graph <- ggplot(trust_coverage_pcsp_cats, aes(x = percent_group, y = trust_count)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") + 
  labs(x = "Proportion of patients offered a PCSP", y = "Number of trusts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))  


####HNA coverage map-------
library(sf)
library(tmap)

