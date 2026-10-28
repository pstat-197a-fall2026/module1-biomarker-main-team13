import warnings
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import cross_val_score, train_test_split

df = pyreadr.read_r('data/biomarker-clean.RData')['biomarker_clean']
X = df.drop(columns=['group', 'ados'])
y = (df['group'] == 'ASD').astype(int)


# Forward selection
_, p = stats.ttest_ind(X[y == 1], X[y==0], equal_var=False)
candidates = list(X.columns[np.argsort(p)[:30]])

panel, best_auc = [], 0
while True:
	scores = {c:cross_val_score(LogisticRegression(max_iter=5000), X[panel + [c]], y, cv=5, scoring='roc_auc').mean() 
	          for c in candidates if c not in panel}
	
	c_best = max(scores, key=scores.get)
	if scores[c_best] < best_auc + 0.005:
		break
	panel.append(c_best)
	best_auc = scores[c_best]

print('New panel:', panel)

baseline = ['DERM', 'RELT', 'IgD']
results = []
for seed in range(30):
	X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.2, stratify=y, random_state=seed)

	for name, cols in [('baseline', baseline), ('new panel', panel)]:
		fit = LogisticRegression(max_iter=5000).fit(X_train[cols], y_train)
		prob = fit.predict_proba(X_test[cols])[:, 1]
		results.append({'panel': name,
		               'accuracy': accuracy_score(y_test, prob > 0.5),
					   'roc_auc': roc_auc_score(y_test, prob)})

print(pd.DataFrame(results).groupby('panel').agg(['mean', 'std']).round(3))