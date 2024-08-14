################ HNA/PCSP coverage analysis - eHNA comparison ###################

#Analysis for proposal aim XXXX - comparing COSD and eHNA numbers 

#Created July 2024 by Lizzie Augarde 
#Change log:
#25/07/2024 added graph
#01/08/2024 adjusted to use HNAs recorded as occurring in 2021 from COSD
##########################################################################################

#extracting HNA data - occurring in 2021 rather than for those diagnosed in 2021
query <- "select 
  rp.nhsnumber,
  rp.patientid,       
  rp.event_type,
  rp.event_property_1,
  rp.event_property_2,
  rp.event_property_3,
  rp.event_date,
  rp.event_end,
  rp.trust_code
  from analysisncr.at_rapid_pathway@cas2406 rp 
  where rp.event_type = 20 
  and (rp.event_date >= '01-JAN-2021' AND rp.event_date <= '31-DEC-2021')"

hnas_2021_data <- dbGetQueryOracle(cas2406, query, rowlimit = NA)

#pre-processing COSD data
hnas_2021_data <- hnas_2021_data |> 
  clean_names() |>
  unique() |> 
  separate(event_property_1, c("point_of_pathway", "offered_code", "staff_role"), ":") |>
  filter(offered_code == "03")|>
  mutate(month = month.abb[month(event_date)]) |>
  group_by(trust_code, month) |>
  summarise(count = n()) |>
  pivot_wider(names_from = month, values_from = count) |>
  rowwise() |>
  mutate(cosd_total = sum(c_across(where(is.numeric)))) |>
  ungroup() |>
  select(trust_code, cosd_total) |>
  filter(!is.na(cosd_total))

#eHNA data
ehna_data <- read.csv("N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/eHNA_data_2021-limited.csv")

ehna_data <- ehna_data |>
  select(-X) |>
  rename_with(~ gsub("\\.21", "", .)) |>
  rename("diag_trust" = "ODS") |>
  rowwise() |>
  mutate(ehna_total = sum(c_across(where(is.numeric)))) |>
  ungroup() |>
  select(diag_trust, ehna_total) |>
  filter(!is.na(ehna_total), diag_trust != "England")


#number of trusts missing in each dataset 
not_in_ehna <- setdiff(hnas_2021_data$trust_code, ehna_data$diag_trust)
not_in_cosd <- setdiff(ehna_data$diag_trust, hnas_2021_data$trust_code)


#merging and comparing the two datasets 
merged_data <- full_join(hnas_2021_data, ehna_data, 
                         by = c("trust_code" = "diag_trust")) |>
  mutate(difference = cosd_total-ehna_total,
         status = case_when((is.na(ehna_total) & cosd_total > 0) ~ "In COSD, not in eHNA",
                            (is.na(cosd_total) & ehna_total > 0) ~ "In eHNA, not in COSD",
                            (cosd_total > 0 & ehna_total > 0) ~ "In both"),
         ehna_total = ifelse(is.na(ehna_total), 0, ehna_total),
         cosd_total = ifelse(is.na(cosd_total), 0, cosd_total))
        

more_in_cosd <- length(which(merged_data$difference > 0))
more_in_ehna <- length(which(merged_data$difference < 0))


#scatterplot of the number of records in COSD and eHNA by trust
ehna_cosd_comparison_graph <- ggplot(merged_data, 
                                     aes(x = ehna_total, y = cosd_total, 
                                         color = status)) +
  geom_point(stat = "identity", size = 3) + 
  theme_minimal() +
  scale_color_manual(values = c("#008A26", "#02D462", "#03fca5")) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Number of records in Macmillan eHNA", 
       y = "Number of records in COSD", 
       color = "Trust status in each dataset")

  