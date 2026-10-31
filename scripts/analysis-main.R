# Load Packages and Data
library(tidyverse)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')


# Question 1 
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

## Histograms
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

## Histograms after log transformation
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


# Notes: 
# right skewed before, more central now. 
# Most subjects have low proteins and some have high 
# This reduces influence of extremes (compresses extremes but does not remove them)
# This does not mean every protein is normally distributed 


# Question 2
## outliers 
outliers_sub<- biomarker_clean %>%
  mutate(subject = row_number()) %>%
  pivot_longer(-c(group, ados, subject),
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(subject, group) %>%
  summarize(n_outliers = sum(abs(level) > 3, na.rm = TRUE)) %>%
  arrange(desc(n_outliers))



## outliers by group
outliers_sub%>%
  group_by(group) %>%
  summarize(
    total = sum(n_outliers),
    mean = mean(n_outliers),
    subjects = sum(n_outliers > 0),
    n = n()
  )

outliers_sub

# Question 3

## Original results

set.seed(101422)

test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()


rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)


proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)


biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

og_results <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(class = factor(class),
         predicted_class = factor(pred > 0.5)
  ) %>%
  class_metrics(estimate = predicted_class,
                truth = class, pred,
                event_level = 'second'
  )

## Training and Test data design change, 71/29 split

set.seed(197)

subset <- sample(nrow(biomarker_clean), 110)

biomarker_train <- biomarker_clean[subset, ]

biomarker_test <- biomarker_clean[-subset, ]

ttests_out <- biomarker_train %>%
  select(-ados) %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)


predictors <- biomarker_train %>%
  select(-c(group, ados))

response <- biomarker_train %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

rf_out$confusion

proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)


proteins_sstar <- intersect(proteins_s1, proteins_s2)

sstar <- function(df) {df %>%
    select(group, any_of(proteins_sstar)) %>%
    mutate(class = (group == 'ASD')) %>%
    select(-group)}

training_sstar <- sstar(biomarker_train)
testing_sstar <- sstar(biomarker_test)

fit <- glm(class ~ ., 
           data = training_sstar, 
           family = 'binomial')

class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)


train_results <- testing_sstar %>%
  add_predictions(fit, type = 'response') %>%
  mutate(class = factor(class),
         predicted_class = factor(pred > 0.5)
  ) %>%
  class_metrics(estimate = predicted_class,
                truth = class, pred,
                event_level = 'second'
  )

## 25 protein subset selection

set.seed(197)

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 25) %>%
  pull(protein)



predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

rf_out$confusion

proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 25) %>%
  pull(protein)


proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)


biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

num_results <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(class = factor(class),
         predicted_class = factor(pred > 0.5)
  ) %>%
  class_metrics(estimate = predicted_class,
                truth = class, pred,
                event_level = 'second'
  )

print(num_results)

data.frame(metric = og_results$.metric, 
           diff = num_results$.estimate - og_results$.estimate)

## Fuzzy intersection 1

## Assigning protein membership score by min-max normalization

## Rank all proteins by RF importance
rf_membership <- rf_out$importance %>%
  as_tibble(rownames = "protein") %>%
  mutate(
    rf_rank = rank(-MeanDecreaseGini),
    mu_rf = 1 - (rf_rank - 1) / (n() - 1)
  ) %>%
  select(protein, rf_rank, mu_rf)

## Rank all proteins by adjusted t-test p-value
t_membership <- ttests_out %>%
  mutate(
    t_rank = rank(p.adj),
    mu_t = 1 - (t_rank - 1) / (n() - 1)
  ) %>%
  select(protein, t_rank, mu_t)

## Combine memberships using the minimum operator
fuzzy_proteins <- inner_join(
  rf_membership, t_membership,
  by = "protein"
) %>%
  mutate(
    mu_fuzzy = pmin(mu_rf, mu_t)
  ) %>%
  arrange(desc(mu_fuzzy))

## Select the proteins with the highest fuzzy memberships
proteins_fuzzy <- fuzzy_proteins %>%
  slice_max(mu_fuzzy, n = 10, with_ties = FALSE) %>%
  pull(protein)

biomarker_fuzzy <- biomarker_clean %>%
  select(group, any_of(proteins_fuzzy)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

## partition into training and test set

biomarker_split <- biomarker_fuzzy %>%
  initial_split(prop = 0.8)

## fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

## evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

fuzzy_results1 <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(class = factor(class),
         predicted_class = factor(pred > 0.5)
  ) %>%
  class_metrics(estimate = predicted_class,
                truth = class, pred,
                event_level = 'second'
  )

## Fuzzy intersection 2

## takes only the (1-adjusted p-values) and protein names 
## makes into data frame
proteins_s11 <- ttests_out %>% 
  mutate(anti_pvalue = 1-p.adj) %>% 
  select(protein, anti_pvalue)

## keeps every MeanDecreaseGini importance value
proteins_s22 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>% 
  select(protein, MeanDecreaseGini)

## performs a fuzzy intersection 
## takes the average value of (1-p_value) and MeanDecreaseGini for each protein
## takes the top 10 proteins with the highest average values
proteins_sstar2 <- proteins_s11 %>% 
  inner_join(proteins_s22) %>% 
  mutate(avg = (anti_pvalue + MeanDecreaseGini)/2) %>%
  slice_max(avg, n = 10) %>%
  pull(protein)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar2)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)


biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

fuzzy_results2 <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(class = factor(class),
         predicted_class = factor(pred > 0.5)
  ) %>%
  class_metrics(estimate = predicted_class,
                truth = class, pred,
                event_level = 'second'
  )

# Question 3 metric results

## Original

print(og_results)

## Train/Test split

print(train_results)

data.frame(metric = og_results$.metric, 
           diff =  train_results$.estimate - og_results$.estimate)

## 25 Protein subset

print(num_results)

data.frame(metric = og_results$.metric, 
           diff = num_results$.estimate - og_results$.estimate)

## Fuzzy 1

print(fuzzy_results1)

data.frame(metric = og_results$.metric, 
           diff = fuzzy_results1$.estimate - og_results$.estimate)

## Fuzzy 2

print(fuzzy_results2)

data.frame(metric = og_results$.metric, 
           diff = fuzzy_results2$.estimate - og_results$.estimate)

# Question 4

## Setup
set.seed(101422)
baseline_panel <- c('DERM', 'RELT', 'IgD')
biomarker <- biomarker_clean %>%
	select(-ados) %>%
	mutate(class = (group == 'ASD')) %>%
	select(-group)

## Candidate Proteins
p_values <- biomarker %>%
	select(-class) %>%
	map_dbl(~ t.test(.x[biomarker$class], .x[!biomarker$class])$p.value)

candidates <- names(sort(p_values))[1:30]

## Forward Selection
set.seed(101422)
folds <- sample(rep(1:5, length.out = nrow(biomarker)))

cv_auc <- function(panel){
	probs <- numeric(nrow(biomarker))
	for(k in 1:5){
		train <- biomarker[folds != k, c('class', panel)]
		test <- biomarker[folds == k, c('class', panel)]
		fit <- glm(class ~ ., data = train, family = 'binomial')
		probs[folds == k] <- predict(fit, newdata = test, type = 'response')
	}
	roc_auc_vec(truth = factor(biomarker$class), estimate = probs, event_level = 'second')
}

panel <- c()
best_auc <- 0
repeat {
	remaining <- setdiff(candidates, panel)
	scores <- map_dbl(remaining, ~ cv_auc(c(panel, .x)))
	if (max(scores) < best_auc + 0.005) break
	panel <- c(panel, remaining[which.max(scores)])
	best_auc <- max(scores)
	cat('added', remaining[which.max(scores)], '- CV AUC = ', round(best_auc, 3), '\n')
}

new_panel <- panel
panel

## Single Split Benchmark
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)
evaluate_panel <- function(panel, split){
	fit <- glm(class ~ ., data = training(split) %>% select(class, all_of(panel)), 
	family = 'binomial')

	testing(split) %>%
		add_predictions(fit, type = 'response') %>%
		mutate(class = factor(class), predicted_class = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
		class_metrics(estimate = predicted_class, truth = class, pred, event_level = 'second')
}

set.seed(101422)
biomarker_split <- initial_split(biomarker, prop = 0.8)

single_split <- bind_rows(
	evaluate_panel(baseline_panel, biomarker_split) %>% mutate(panel = 'baseline'), 
	evaluate_panel(new_panel, biomarker_split) %>% mutate(panel = 'new')
) %>%
	select(panel, .metric, .estimate) %>%
	pivot_wider(names_from = panel, values_from = .estimate)

## Repeated Split Benchmark
repeated <- map_dfr(1:30, function(seed){
	set.seed(seed)
	split <- initial_split(biomarker, prop = 0.8)
	bind_rows(
		evaluate_panel(baseline_panel, split) %>% mutate(panel = 'baseline'),
		evaluate_panel(new_panel, split) %>% mutate(panel = 'new')
	)
})

repeated_summary <- repeated %>%
	group_by(panel, .metric) %>%
	summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = 'drop') %>%
	pivot_wider(names_from = panel, values_from = c(mean, sd))

print(new_panel)
print(single_split)
print(repeated_summary)