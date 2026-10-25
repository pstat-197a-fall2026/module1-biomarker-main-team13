library(tidyverse)

var_names <- read_csv(
  "data/biomarker-raw.csv",
  col_names = FALSE,
  n_max = 2,
  col_select = -(1:2)
) %>%
  t() %>%
  as_tibble() %>%
  rename(
    name = V1,
    abbreviation = V2
  ) %>%
  na.omit()

biomarker_raw <- read_csv(
  "data/biomarker-raw.csv",
  skip = 2,
  col_select = -2L,
  col_names = c(
    "group",
    "empty",
    pull(var_names, abbreviation),
    "ados"
  ),
  na = c("-", "")
) %>%
  filter(!is.na(group))


set.seed(123)

random_proteins <- biomarker_raw %>%
  select(-any_of(c("group", "empty", "ados"))) %>%
  select(where(is.numeric)) %>%
  names() %>%
  sample(4)

print(random_proteins)

# Histograms
biomarker_raw %>%
  select(all_of(random_proteins)) %>%
  pivot_longer(
    everything(),
    names_to = "protein",
    values_to = "level"
  ) %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 30) +
  facet_wrap(~protein, scales = "free_x") +
  labs(
    title = "Distribution of 4 Randomly Selected Proteins",
    x = "Raw Protein Level",
    y = "Number of Subjects"
  )

#Histograms after log transformation
biomarker_raw %>%
  select(all_of(random_proteins)) %>%
  pivot_longer(
    everything(),
    names_to = "protein",
    values_to = "level"
  ) %>%
  filter(level > 0) %>%
  ggplot(aes(x = log10(level))) +
  geom_histogram(bins = 30) +
  facet_wrap(~protein, scales = "free_x") +
  labs(
    title = "Distribution of 4 Log-Transformed Proteins",
    x = "Log10 Protein Level",
    y = "Number of Subjects"
  )


#Notes: 
# right skewed before, more central now. 
# Most subjects have low proteins and some have high 
# This reduces influence of extremes (compresses extremes but does not remove them)
# This does not mean every protein is normally distributed 