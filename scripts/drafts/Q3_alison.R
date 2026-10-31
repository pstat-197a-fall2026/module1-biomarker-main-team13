# Load Packages and Data
library(tidyverse)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')




#Original Analysis:

#load('data/biomarker-clean.RData')
#load('~/Desktop/Pstat 197A/module1-biomarker-main-team13/data/biomarker-clean.RData')

## MULTIPLE TESTING
####################

# function to compute tests
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

set.seed(101422)
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

set.seed(101422)
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

og_results

# Original Results:
# sensitivity: 0.9375000
# specificity: 0.7333333
# accuracy: 0.8387097
# roc_auc: 0.8958333


# Split Data
set.seed(197)

#80/20 training-test split
biomarker_split <- initial_split(
  biomarker_clean,
  prop = 0.8,
  strata = group_bootstrap
)

biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)


## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_train <- biomarker_train %>%
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

proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)


proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

train_results <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(class = factor(class),
         predicted_class = factor(pred > 0.5)
  ) %>%
  class_metrics(estimate = predicted_class,
                truth = class, pred,
                event_level = 'second'
  )

train_results

