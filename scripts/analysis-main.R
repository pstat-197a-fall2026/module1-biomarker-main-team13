# Load Packages and Data
library(tidyverse)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')













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