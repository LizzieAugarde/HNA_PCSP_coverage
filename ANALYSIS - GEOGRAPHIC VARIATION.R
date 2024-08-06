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
  mutate(percent_group = fct_relevel(percent_group, "Less than 10%")) |>
  arrange(percent_group)

trust_hna_graph <- ggplot(trust_coverage_hna_cats, aes(x = percent_group, y = trust_count)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008D26") + 
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
  mutate(percent_group = fct_relevel(percent_group, "Less than 10%")) |>
  arrange(percent_group)

trust_pcsp_graph <- ggplot(trust_coverage_pcsp_cats, aes(x = percent_group, y = trust_count)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#008D26") + 
  labs(x = "Proportion of patients offered a PCSP", y = "Number of trusts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))  


####HNA coverage map by CA-------
library(sf)
library(tmap)

query <- "select a.tumourid,
                 b.icb_2022_name,
                 b.icb_2022_code
          from analysiselizabethaugarde.hna_pcsps_patient_cohort@casref01 a
          left join av2021.at_geography_england@casref01 b on a.tumourid = b.tumourid"

patients_geog <- dbGetQueryOracle(cas2406, query, rowlimit = NA)

patients_geog <- left_join(patients_geog, patient_level_data, by = c("TUMOURID" = "tumourid"))

icb_coverage_hna <- patients_geog |>
  filter(keepforhna == "INCLUDE") |>
  group_by(ICB_2022_CODE, hna_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ICB_2022_CODE) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(ICB_2022_CODE), hna_status == "Has HNA") |>
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

icb_coverage_pcsp <- patients_geog |>
  filter(keepforpcsp == "INCLUDE") |>
  group_by(ICB_2022_CODE, pcsp_status) |>
  summarise(number_patients = n()) |>
  ungroup() |>
  group_by(ICB_2022_CODE) |>
  mutate(percent = (number_patients/sum(number_patients))*100,
         percent_table = percent((number_patients/sum(number_patients)), accuracy = 0.1)) |>
  filter(!is.na(ICB_2022_CODE), pcsp_status == "Has PCSP") |>
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

#maps
icbs <- read_sf("N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/ICB_JUL_2022_EN_BFC_V2.shp")

icbs_hna <- left_join(icbs, icb_coverage_hna, by = c("ICB22CD" = "ICB_2022_CODE"))
icbs_pcsp <- left_join(icbs, icb_coverage_pcsp, by = c("ICB22CD" = "ICB_2022_CODE"))

hna_map <- tm_shape(icbs_hna) +
  tm_polygons("percent_group", title = "HNA coverage", palette = "Greens")

tmap_save(hna_map, filename = "hna_map.jpeg", dpi = 300, width = 8, height = 6, units = "in")

pcsp_map <- tm_shape(icbs_pcsp) +
  tm_polygons("percent_group", title = "PCSP coverage", palette = "Greens")

tmap_save(pcsp_map, filename = "pcsp_map.jpeg", dpi = 300, width = 8, height = 6, units = "in")
