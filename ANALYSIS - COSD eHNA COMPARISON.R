################ HNA/PCSP coverage analysis - eHNA comparison ###################

#Analysis for proposal aim XXXX - comparing COSD and eHNA numbers 

#Created July 2024 by Lizzie Augarde 
#Change log:
#25/07/2024 added graph
##########################################################################################


#read in and pre-processing
hna_by_trust <- hna_data |>
  mutate(month = month.abb[month(hna_date)]) |>
  filter(hna_offered_code == "03") |>
  group_by(diag_trust, month) |>
  summarise(count = n()) |>
  pivot_wider(names_from = month, values_from = count) |>
  rowwise() |>
  mutate(total_sum = sum(c_across(where(is.numeric)))) |>
  ungroup() |>
  select(diag_trust, total_sum) |>
  filter(!is.na(total_sum))

ehna_data <- read.csv("N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/eHNA_data_2021-limited.csv")

ehna_data <- ehna_data |>
  select(-X) |>
  rename_with(~ gsub("\\.21", "", .)) |>
  rename("diag_trust" = "ODS") |>
  rowwise() |>
  mutate(total_sum = sum(c_across(where(is.numeric)))) |>
  ungroup() |>
  select(diag_trust, total_sum) |>
  filter(!is.na(total_sum), diag_trust != "England")


#trusts missing in each dataset 
not_in_ehna <- setdiff(hna_by_trust$diag_trust, ehna_data$diag_trust) #55 trusts
not_in_cosd <- setdiff(ehna_data$diag_trust, hna_by_trust$diag_trust) #21 trusts


#merging and comparing the two datasets 
merged_data <- inner_join(hna_by_trust, ehna_data, by = "diag_trust")

merged_data <- merged_data |> 
  mutate(difference = total_sum.x-total_sum.y) 

more_in_cosd <- length(which(merged_data$difference > 0))
more_in_ehna <- length(which(merged_data$difference < 0))


#graph of differences between COSD and eHNA numbers of HNAs
differences_graph <- ggplot(merged_data, 
                            aes(x = reorder(diag_trust, difference), y = difference)) +
  geom_bar(stat = "identity") + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(x = "Trust", y = "Difference in number of HNAs\n(COSD minus eHNA)")

  