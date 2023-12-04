# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/11/13 0:53
import math

from joblib import dump
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve
from sklearn.naive_bayes import GaussianNB
import pandas as pd
import numpy as np
df = pd.read_csv('../data/mutation_matrix2.0.csv')

# Drop rows with any NaN values
df = df.dropna()

# Map string labels to integers
df.iloc[:, -1] = df.iloc[:, -1].map({'human': 0, 'human_origin_avian': 1, 'avian': 2})

df_copy = df.copy()
df_copy.iloc[:, -1] = df_copy.iloc[:, -1].map(({2: 0, 0: 1, 1: 1}))
id_type_dict = dict(zip(df_copy['Strain ID'], df_copy['Sequence Type']))

# 设置索引
df.set_index('Strain ID', inplace = True)

# 转化标签为整数（如果还未转化）
df.iloc[:, -1] = df.iloc[:, -1].astype(int)
X = df.iloc[:, :-1]
y = df.iloc[:, -1]

# Split data into training and validation sets
df0 = df[df['Sequence Type'] == 0]
df1 = df[df['Sequence Type'] == 1]
df2 = df[df['Sequence Type'] == 2]

single_test_size = math.floor(len(df1) / 2)

df0_test = df0.sample(n = single_test_size, replace = False, random_state = 42)
df1_test = df1.sample(n = single_test_size, replace = False, random_state = 42)
df2_test = df2.sample(n = single_test_size * 2, replace = False, random_state = 42)

df_test = pd.concat([df0_test, df1_test, df2_test])
df0_train = df0.drop(df0_test.index).sample(n = int(single_test_size * 2 * 7.5 / 1.5), replace = False,
                                            random_state = 42)
df2_train = df2.drop(df2_test.index).sample(n = int(single_test_size * 2 * 7.5 / 1.5), replace = False,
                                            random_state = 42)
df1_validation = df1.drop(df1_test.index).sample(n = single_test_size, replace = False, random_state = 42)
df0_validation = df0.drop(labels = list(df0_test.index) + list(df0_train.index)).sample(n = single_test_size,
                                                                                        replace = False,
                                                                                        random_state = 42)
df2_validation = df2.drop(labels = list(df2_test.index) + list(df2_train.index)).sample(n = single_test_size * 2,
                                                                                        replace = False,
                                                                                        random_state = 42)
df_train = pd.concat([df0_train, df2_train])
df_validation = pd.concat([df0_validation, df1_validation, df2_validation])


X_train = df_train.iloc[:, 1:-1]
y_train = df_train.iloc[:, -1].replace({2: 0, 0: 1})

X_val = df_validation.iloc[:, 1:-1]
y_val = df_validation.iloc[:, -1].replace({2: 0, 0: 1})

# Train the Random Forest model on the training set
rf = RandomForestClassifier(random_state=42)
rf.fit(X_train, y_train)

# Get the feature importances from the Random Forest model
importances = rf.feature_importances_

# Sort the feature importances in descending order and select the top 12 features
indices = np.argsort(importances)[::-1]
top_features = X_train.columns[indices[:12]]

# Train the Gaussian Naive Bayes model on the top features from the training set
gnb = GaussianNB()
gnb.fit(X_train[top_features], y_train)

# Predict probabilities on the validation set
y_val_proba = gnb.predict_proba(X_val[top_features])[:, 1]
print(top_features)
# Compute ROC curve and select the best threshold
fpr, tpr, thresholds = roc_curve(y_val, y_val_proba)
optimal_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[optimal_idx]

# Save the trained model and the optimal threshold to disk
dump(gnb, '../model/gnb_model.joblib')
dump(optimal_threshold, '../model/optimal_threshold.joblib')
dump(top_features,'../model/top_features.joblib')