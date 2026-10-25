library(tidyverse)

# Step 1: Get protein names
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

# Step 2: Load raw data (no transformations)
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

# Step 3: Select four proteins
proteins <- biomarker_raw %>%
  select(DERM, RELT, MRC2, IgD)

# Step 4: Make histograms
proteins %>%
  pivot_longer(
    everything(),
    names_to = "protein",
    values_to = "level"
  ) %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 30) +
  facet_wrap(~protein, scales = "free_x") +
  labs(
    title = "Distribution of Raw Protein Levels",
    x = "Raw Protein Level",
    y = "Number of Subjects"
  )
