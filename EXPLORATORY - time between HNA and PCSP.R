################ HNA/PCSP coverage analysis - EXPLORATORY - time between HNA and PCSP ###################

#Exploring time between first HNA and first PCSP 

#Created August 2024 by Lizzie Augarde 
##########################################################################################


###### TIME BETWEEN HNA AND PCSP ######
time_between_hna_pcsp <- patient_level_data |>
  mutate(time_between = as.numeric(difftime(pcsp_date, hna_date, unit = "days"))) |>
  filter(hna_count != 0 & pcsp_count != 0)

#median time between HNA and PCSP
median_time_between_hna_pcsp <- median(time_between_hna_pcsp$time_between)
mean_time_between_hna_pcsp <- mean(time_between_hna_pcsp$time_between)

median_time_between_hist <- ggplot(time_between_hna_pcsp, aes(x = time_between)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  scale_y_continuous(limits = c(0,80000),
                     breaks = c(20000,40000,60000,80000), 
                     labels = c("20,000", "40,000", "60,000", "80,000")) +
  #scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730), #adding in breaks for time periods of interest
  #labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
  #"6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

median_time_between_hist <- ggplot(time_between_hna_pcsp, aes(x = time_between)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black") +
  #scale_y_continuous(limits = c(0,80000),
  #breaks = c(20000,40000,60000,80000), 
  #labels = c("20,000", "40,000", "60,000", "80,000")) +
  #scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730), #adding in breaks for time periods of interest
  #labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
  #"6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

time_between_summary <- time_between_hna_pcsp |>
  group_by(time_between) |>
  summarise(pat_count = n())
