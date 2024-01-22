##### HNAs by month for comparison to Level 2

hna_raw_data <- dbGetQueryOracle(casref01, "SELECT * FROM analysisncr.at_rapid_pathway@cas2311
                             where event_type = 20 and event_date > TO_DATE('31-DEC-20', 'DD-MON-RR')", rowlimit = NA)

hna_data <- hna_raw_data %>%
  clean_names() %>%
  mutate(event_year = year(as.Date(event_date))) %>%
  mutate(event_month = month(as.Date(event_date)))

hnas_by_month <- hna_data %>%
  group_by(event_year, event_month) %>%
  summarise(count = n())
