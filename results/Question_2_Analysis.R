# outliers 
outliers_sub<- biomarker_clean %>%
  mutate(subject = row_number()) %>%
  pivot_longer(-c(group, ados, subject),
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(subject, group) %>%
  summarize(n_outliers = sum(abs(level) > 3, na.rm = TRUE)) %>%
  arrange(desc(n_outliers))



# outliers by group
outliers_sub%>%
  group_by(group) %>%
  summarize(
    total = sum(n_outliers),
    mean = mean(n_outliers),
    subjects = sum(n_outliers > 0),
    n = n()
  )

outliers_sub
