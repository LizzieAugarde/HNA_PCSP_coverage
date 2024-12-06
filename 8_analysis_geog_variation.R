################ HNA/PCSP coverage analysis - geographic variation ###################

#Analysis of the proportion of patients being offered an HNA/PCSP by trust of diagnosis
#Used in Section 9g of initial results slide deck 

#Created July 2024 by Lizzie Augarde 
#################################################################################

############# HNA COVERAGE CATEGORISATION #############
trust_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(diag_trust, hna_status) |> #aggregating to trust level
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(diag_trust) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(diag_trust), hna_status == "Has HNA") |>
  mutate(percent_group = case_when(percent <10 ~ "Less than 10%", #grouping
                                   percent >=10 & percent <20 ~ "10% - 20%",
                                   percent >=20 & percent <30 ~ "20% - 30%",
                                   percent >=30 & percent <40 ~ "30% - 40%",
                                   percent >=40 & percent <50 ~ "40% - 50%",
                                   percent >=50 & percent <60 ~ "50% - 60%",
                                   percent >=60 & percent <70 ~ "60% - 70%",
                                   percent >=70 & percent <80 ~ "70% - 80%",
                                   percent >=80 & percent <90 ~ "80% - 90%",
                                   percent >=90 & percent <100 ~ "90% - 100%"))

#dummy table of all percentage groups 
percent_cats <- data.frame(percent_group = c("Less than 10%", "10% - 20%", "20% - 30%", "30% - 40%", "40% - 50%",
                                             "50% - 60%", "60% - 70%", "70% - 80%", "80% - 90%", "90% - 100%"),
                           trust_count_dummy = 0)

#number of trusts with each % coverage
trust_coverage_hna_cats <- trust_hna |> 
  group_by(percent_group) |>
  summarise(trust_count = n()) |>
  full_join(percent_cats, by = "percent_group") |>
  mutate(trust_count = ifelse(is.na(trust_count), trust_count_dummy, trust_count)) |>
  select(-trust_count_dummy) |>
  mutate(percent_group = fct_relevel(percent_group, "Less than 10%")) |>
  arrange(percent_group)

#bar graph of the number of trusts with each % coverage of HNAs
trust_hna_graph <- ggplot(trust_coverage_hna_cats, aes(x = percent_group, y = trust_count)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008D26") + 
  labs(x = "Proportion of patients offered a HNA", y = "Number of trusts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))    


############# PCSP COVERAGE CATEGORISATION #############
trust_pcsp <- patient_level_data |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(diag_trust, pcsp_status) |> #aggregating to trust level
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(diag_trust) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(diag_trust), pcsp_status == "Has PCSP") |>
  mutate(percent_group = case_when(percent <10 ~ "Less than 10%", #grouping
                                   percent >=10 & percent <20 ~ "10% - 20%",
                                   percent >=20 & percent <30 ~ "20% - 30%",
                                   percent >=30 & percent <40 ~ "30% - 40%",
                                   percent >=40 & percent <50 ~ "40% - 50%",
                                   percent >=50 & percent <60 ~ "50% - 60%",
                                   percent >=60 & percent <70 ~ "60% - 70%",
                                   percent >=70 & percent <80 ~ "70% - 80%",
                                   percent >=80 & percent <90 ~ "80% - 90%",
                                   percent >=90 & percent <100 ~ "90% - 100%"))

#number of trusts with each % coverage
trust_coverage_pcsp_cats <- trust_pcsp |> 
  group_by(percent_group) |>
  summarise(trust_count = n()) |>
  full_join(percent_cats, by = "percent_group") |>
  mutate(trust_count = ifelse(is.na(trust_count), trust_count_dummy, trust_count)) |>
  select(-trust_count_dummy) |>
  mutate(percent_group = fct_relevel(percent_group, "Less than 10%")) |>
  arrange(percent_group)

#bar graph of the number of trusts with each % coverage of PCSPs
trust_pcsp_graph <- ggplot(trust_coverage_pcsp_cats, aes(x = percent_group, y = trust_count)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008D26") + 
  labs(x = "Proportion of patients offered a PCSP", y = "Number of trusts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))  


############# COVERAGE BY ICB #############
#aggregating HNA coverage to ICB level
icb_coverage_hna <- patient_level_data |>
  filter(keepforhna == "INCLUDE") |>
  group_by(icb_2022_code, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(icb_2022_code) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(icb_2022_code), hna_status == "Has HNA") |>
  mutate(percent_group = case_when(percent <10 ~ "Less than 10%",
                                   percent >=10 & percent <20 ~ "10% - 20%",
                                   percent >=20 & percent <30 ~ "20% - 30%",
                                   percent >=30 & percent <40 ~ "30% - 40%",
                                   percent >=40 & percent <50 ~ "40% - 50%",
                                   percent >=50 & percent <60 ~ "50% - 60%",
                                   percent >=60 & percent <70 ~ "60% - 70%",
                                   percent >=70 & percent <80 ~ "70% - 80%",
                                   percent >=80 & percent <90 ~ "80% - 90%",
                                   percent >=90 & percent <100 ~ "90% - 100%")) |>
  mutate(percent_group = factor(percent_group, levels = c("Less than 10%",
                                                          "10% - 20%", "20% - 30%",
                                                          "30% - 40%", "40% - 50%",
                                                          "50% - 60%", "60% - 70%",
                                                          "70% - 80%", "80% - 90%",
                                                          "90% - 100%"))) |>
  ungroup()

#aggregating PCSP coverage to ICB level
icb_coverage_pcsp <- patients_geog |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(icb_2022_code, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(icb_2022_code) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(icb_2022_code), pcsp_status == "Has PCSP") |>
  mutate(percent_group = case_when(percent <10 ~ "Less than 10%",
                                   percent >=10 & percent <20 ~ "10% - 20%",
                                   percent >=20 & percent <30 ~ "20% - 30%",
                                   percent >=30 & percent <40 ~ "30% - 40%",
                                   percent >=40 & percent <50 ~ "40% - 50%",
                                   percent >=50 & percent <60 ~ "50% - 60%",
                                   percent >=60 & percent <70 ~ "60% - 70%",
                                   percent >=70 & percent <80 ~ "70% - 80%",
                                   percent >=80 & percent <90 ~ "80% - 90%",
                                   percent >=90 & percent <100 ~ "90% - 100%")) |>
  mutate(percent_group = factor(percent_group, levels = c("Less than 10%",
                                                          "10% - 20%", "20% - 30%",
                                                          "30% - 40%", "40% - 50%",
                                                          "50% - 60%", "60% - 70%",
                                                          "70% - 80%", "80% - 90%",
                                                          "90% - 100%"))) |>
  ungroup()


############# MAPS OF COVERAGE BY ICB #############
library(sf)
library(tmap)

#read in boundaries
icbs <- read_sf("N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/ICB_JUL_2022_EN_BFC_V2.shp")

#joining boundary to data 
icbs_hna <- left_join(icbs, icb_coverage_hna, by = c("ICB22CD" = "icb_2022_code"))
icbs_pcsp <- left_join(icbs, icb_coverage_pcsp, by = c("ICB22CD" = "icb_2022_code"))

#HNA coverage map and save
hna_map <- tm_shape(icbs_hna) +
  tm_polygons("percent_group", title = "HNA coverage", palette = "Greens")

tmap_save(hna_map, filename = "hna_map.jpeg", dpi = 300, width = 8, height = 6, units = "in")

#PCSP coverage map and save
pcsp_map <- tm_shape(icbs_pcsp) +
  tm_polygons("percent_group", title = "PCSP coverage", palette = "Greens")

tmap_save(pcsp_map, filename = "pcsp_map.jpeg", dpi = 300, width = 8, height = 6, units = "in")
