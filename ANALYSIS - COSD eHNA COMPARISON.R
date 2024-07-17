################ HNA/PCSP coverage analysis - eHNA comparison ###################

#Analysis for proposal aim XXXX - comparing COSD and eHNA numbers 

#Created July 2024 by Lizzie Augarde 
#Change log:
##########################################################################################


#read in and pre-processing
hna_by_trust <- hna_data %>%
  mutate(month = month.abb[month(hna_date)]) %>%
  filter(hna_offered_code == "03") %>%
  group_by(diag_trust, month) %>%
  summarise(count = n()) %>%
  pivot_wider(., names_from = month, values_from = count)

ehna_data <- read.csv("N:/INFO/_LIVE/NCIN/Macmillan_Partnership/HNAs/COSD level 3 analysis/Data/eHNA_data_2021.csv")

ehna_data <- ehna_data %>%
  select(-X) %>%
  rename_with(~ gsub("\\.21", "", .)) %>%
  rename("diag_trust" = "ODS")


#trusts missing in each dataset 
not_in_ehna <- setdiff(hna_by_trust$diag_trust, ehna_data$diag_trust) #55 trusts
not_in_cosd <- setdiff(ehna_data$diag_trust, hna_by_trust$diag_trust) #21 trusts


#merging and comparing the two datasets 
merged_data <- merge(hna_by_trust, ehna_data, by = "diag_trust", suffixes = c("_cosd", "_ehna"))

months_list <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

difference_data <- merged_data
for (month in months_list) {
  difference_data[[paste0(month, "_diff")]] <- merged_data[[paste0(month, "_cosd")]] - merged_data[[paste0(month, "_ehna")]]
}
