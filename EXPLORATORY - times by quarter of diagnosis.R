################ HNA/PCSP coverage analysis - EXPLORATORY - times by quarter ###################

#Exploring time between diagnosis and first HNA/PCSP by quarter of diagnosis 
#Follows query from Lesley Smith 19/08/2024

#Created August 2024 by Lizzie Augarde 
##########################################################################################

library(ggplot2)
library(tidyverse) 
library(gridExtra)

patient_level_data <- patient_level_data |>
  mutate(hna_time_diag_event = as.numeric(hna_time_diag_event),
         pcsp_time_diag_event = as.numeric(pcsp_time_diag_event)) |>
  filter(!is.na(diagnosisdatebest)) |>
  mutate(quarter = case_when(month(diagnosisdatebest) %in% c(1,2,3) ~ "Q1",
                             month(diagnosisdatebest) %in% c(4,5,6) ~ "Q2",
                             month(diagnosisdatebest) %in% c(7,8,9) ~ "Q3",
                             month(diagnosisdatebest) %in% c(10,11,12) ~ "Q4"))



hna_hist_q1 <- ggplot(filter(patient_level_data, patient_level_data$quarter == "Q1"), 
                      aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  #scale_y_continuous(limits = c(0,80000),
                     #breaks = c(20000,40000,60000,80000), 
                     #labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

hna_hist_q2 <- ggplot(filter(patient_level_data, patient_level_data$quarter == "Q2"), 
                      aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  #scale_y_continuous(limits = c(0,80000),
  #breaks = c(20000,40000,60000,80000), 
  #labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

hna_hist_q3 <- ggplot(filter(patient_level_data, patient_level_data$quarter == "Q3"), 
                      aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  #scale_y_continuous(limits = c(0,80000),
  #breaks = c(20000,40000,60000,80000), 
  #labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

hna_hist_q4 <- ggplot(filter(patient_level_data, patient_level_data$quarter == "Q4"), 
                      aes(x = hna_time_diag_event)) +
  geom_histogram(binwidth = 10, fill = "#008A26", color = "black", aes(y = cumsum(..count..))) +
  #scale_y_continuous(limits = c(0,80000),
  #breaks = c(20000,40000,60000,80000), 
  #labels = c("20,000", "40,000", "60,000", "80,000")) +
  scale_x_continuous(breaks = c(30,42,56,84,183,365,548,730),
                     labels = c("4 weeks", "6 weeks", "8 weeks", "12 weeks",
                                "6 months", "1 year", "18 months", "2 years")) +
  labs(x = "Time from diagnosis to first HNA", y = "Number of patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

grid.arrange(hna_hist_q1, hna_hist_q2, hna_hist_q3, hna_hist_q4)


median_time <- patient_level_data |>
  group_by(quarter) |>
  summarise(median_time_hna = median(hna_time_diag_event, na.rm=TRUE),
            median_time_pcsp = median(pcsp_time_diag_event, na.rm=TRUE))
