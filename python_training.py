#!/usr/bin/python
import os
import sys
import math
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

import uproot
import lightgbm as lgb

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

names = ["Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
         "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
         "Naive Bayes", "QDA","LightGBM"]

classifiers = [
    KNeighborsClassifier(3),
    SVC(kernel="linear", C=0.025),
    SVC(gamma=2, C=1),
    GaussianProcessClassifier(1.0 * RBF(1.0)),
    DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    MLPClassifier(alpha=1, max_iter=1000),
    AdaBoostClassifier(),
    GradientBoostingClassifier(n_estimators=1000, learning_rate=0.01,max_depth=10, random_state=0),
    GaussianNB(),
    QuadraticDiscriminantAnalysis(),
    lgb
]

classifiers =[lgb]

file = uproot.open(os.environ['ML_TRAINING_FILE'])
matchTree = file["matchTree"]

MFT_X      = matchTree.array("MFT_X");
MFT_Y      = matchTree.array("MFT_Y");
MFT_Phi    = matchTree.array("MFT_Phi");
MFT_Tanl   = matchTree.array("MFT_Tanl");
MFT_InvQPt = matchTree.array("MFT_InvQPt");
MFT_Cov00  = matchTree.array("MFT_Cov00");
MFT_Cov01  = matchTree.array("MFT_Cov01");
MFT_Cov11  = matchTree.array("MFT_Cov11");
MFT_Cov02  = matchTree.array("MFT_Cov02");
MFT_Cov12  = matchTree.array("MFT_Cov12");
MFT_Cov22  = matchTree.array("MFT_Cov22");
MFT_Cov03  = matchTree.array("MFT_Cov03");
MFT_Cov13  = matchTree.array("MFT_Cov13");
MFT_Cov23  = matchTree.array("MFT_Cov23");
MFT_Cov33  = matchTree.array("MFT_Cov33");
MFT_Cov04  = matchTree.array("MFT_Cov04");
MFT_Cov14  = matchTree.array("MFT_Cov14");
MFT_Cov24  = matchTree.array("MFT_Cov24");
MFT_Cov34  = matchTree.array("MFT_Cov34");
MFT_Cov44  = matchTree.array("MFT_Cov44");
MCH_X      = matchTree.array("MCH_X");
MCH_Y      = matchTree.array("MCH_Y");
MCH_Phi    = matchTree.array("MCH_Phi");
MCH_Tanl   = matchTree.array("MCH_Tanl");
MCH_InvQPt = matchTree.array("MCH_InvQPt");
MCH_Cov00  = matchTree.array("MCH_Cov00");
MCH_Cov01  = matchTree.array("MCH_Cov01");
MCH_Cov11  = matchTree.array("MCH_Cov11");
MCH_Cov02  = matchTree.array("MCH_Cov02");
MCH_Cov12  = matchTree.array("MCH_Cov12");
MCH_Cov22  = matchTree.array("MCH_Cov22");
MCH_Cov03  = matchTree.array("MCH_Cov03");
MCH_Cov13  = matchTree.array("MCH_Cov13");
MCH_Cov23  = matchTree.array("MCH_Cov23");
MCH_Cov33  = matchTree.array("MCH_Cov33");
MCH_Cov04  = matchTree.array("MCH_Cov04");
MCH_Cov14  = matchTree.array("MCH_Cov14");
MCH_Cov24  = matchTree.array("MCH_Cov24");
MCH_Cov34  = matchTree.array("MCH_Cov34");
MCH_Cov44  = matchTree.array("MCH_Cov44");

MFT_TrackChi2   = matchTree.array("MFT_TrackChi2");
MFT_NClust      = matchTree.array("MFT_NClust");
CorrectMatching = matchTree.array("Truth");

Delta_X      = MCH_X      - MFT_X;
Delta_Y      = MCH_Y      - MFT_Y;
Delta_XY     = np.sqrt((MCH_X-MFT_X)**2 + (MCH_Y-MFT_Y)**2)
Delta_Phi    = MCH_Phi    - MFT_Phi;
Delta_Tanl   = MCH_Tanl   - MFT_Tanl;
Delta_InvQPt = MCH_InvQPt - MFT_InvQPt;
Delta_Cov00  = MCH_Cov00 - MFT_Cov00;
Delta_Cov01  = MCH_Cov01 - MFT_Cov01;
Delta_Cov11  = MCH_Cov11 - MFT_Cov11;
Delta_Cov02  = MCH_Cov02 - MFT_Cov02;
Delta_Cov12  = MCH_Cov12 - MFT_Cov12;
Delta_Cov22  = MCH_Cov22 - MFT_Cov22;
Delta_Cov03  = MCH_Cov03 - MFT_Cov03;
Delta_Cov13  = MCH_Cov13 - MFT_Cov13;
Delta_Cov23  = MCH_Cov23 - MFT_Cov23;
Delta_Cov33  = MCH_Cov33 - MFT_Cov33;
Delta_Cov04  = MCH_Cov04 - MFT_Cov04;
Delta_Cov14  = MCH_Cov14 - MFT_Cov14;
Delta_Cov24  = MCH_Cov24 - MFT_Cov24;
Delta_Cov34  = MCH_Cov34 - MFT_Cov34;
Delta_Cov44  = MCH_Cov44 - MFT_Cov44;

Ratio_X      = MCH_X      - MFT_X;
Ratio_Y      = MCH_Y      - MFT_Y;
Ratio_XY     = np.sqrt((MCH_X-MFT_X)**2 + (MCH_Y-MFT_Y)**2)
Ratio_Phi    = MCH_Phi    - MFT_Phi;
Ratio_Tanl   = MCH_Tanl   - MFT_Tanl;
Ratio_InvQPt = MCH_InvQPt - MFT_InvQPt;
Ratio_Cov00  = MCH_Cov00 - MFT_Cov00;
Ratio_Cov01  = MCH_Cov01 - MFT_Cov01;
Ratio_Cov11  = MCH_Cov11 - MFT_Cov11;
Ratio_Cov02  = MCH_Cov02 - MFT_Cov02;
Ratio_Cov12  = MCH_Cov12 - MFT_Cov12;
Ratio_Cov22  = MCH_Cov22 - MFT_Cov22;
Ratio_Cov03  = MCH_Cov03 - MFT_Cov03;
Ratio_Cov13  = MCH_Cov13 - MFT_Cov13;
Ratio_Cov23  = MCH_Cov23 - MFT_Cov23;
Ratio_Cov33  = MCH_Cov33 - MFT_Cov33;
Ratio_Cov04  = MCH_Cov04 - MFT_Cov04;
Ratio_Cov14  = MCH_Cov14 - MFT_Cov14;
Ratio_Cov24  = MCH_Cov24 - MFT_Cov24;
Ratio_Cov34  = MCH_Cov34 - MFT_Cov34;
Ratio_Cov44  = MCH_Cov44 - MFT_Cov44;

MFT_TrackReducedChi2 = MFT_TrackChi2/MFT_NClust;

training_list = np.stack([MFT_X,MFT_Y,MFT_Phi,MFT_Tanl,MFT_InvQPt,MFT_Cov00,MFT_Cov01,MFT_Cov11,MFT_Cov02,MFT_Cov12,MFT_Cov22,MFT_Cov03,MFT_Cov13,MFT_Cov23,MFT_Cov33,MFT_Cov04,MFT_Cov14,MFT_Cov24,MFT_Cov34,MFT_Cov44,
                          MCH_X,MCH_Y,MCH_Phi,MCH_Tanl,MCH_InvQPt,MCH_Cov00,MCH_Cov01,MCH_Cov11,MCH_Cov02,MCH_Cov12,MCH_Cov22,MCH_Cov03,MCH_Cov13,MCH_Cov23,MCH_Cov33,MCH_Cov04,MCH_Cov14,MCH_Cov24,MCH_Cov34,MCH_Cov44,
                          MFT_TrackChi2,MFT_NClust,MFT_TrackReducedChi2,
                          Delta_X,Delta_Y,Delta_XY,Delta_Phi,Delta_Tanl,Delta_InvQPt,Delta_Cov00,Delta_Cov01,Delta_Cov11,Delta_Cov02,Delta_Cov12,Delta_Cov22,Delta_Cov03,Delta_Cov13,Delta_Cov23,Delta_Cov33,Delta_Cov04,Delta_Cov14,Delta_Cov24,Delta_Cov34,Delta_Cov44,
                          Ratio_X,Ratio_Y,Ratio_XY,Ratio_Phi,Ratio_Tanl,Ratio_InvQPt,Ratio_Cov00,Ratio_Cov01,Ratio_Cov11,Ratio_Cov02,Ratio_Cov12,Ratio_Cov22,Ratio_Cov03,Ratio_Cov13,Ratio_Cov23,Ratio_Cov33,Ratio_Cov04,Ratio_Cov14,Ratio_Cov24,Ratio_Cov34,Ratio_Cov44], 1);

X = training_list
y = CorrectMatching

X = StandardScaler().fit_transform(X)

X_train, X_eval, y_train, y_eval = train_test_split(X, y,test_size=0.2,random_state=0,stratify=y)

evals_result = {}
 
train_data = lgb.Dataset(X_train, label=y_train)
valid_data = lgb.Dataset(X_eval, label=y_eval, reference = train_data)

params = {
    'task':'train',
    'boosting_type':'gbdt',
    'objective': 'binary',
    'learning_rate':0.01,
    'max_depth':10,
    'n_estimators':1000,
    'metric':'auc',
    'verbose': 2,
}

model = lgb.train(
    params,
    train_data,
    valid_sets=[train_data,valid_data],
    valid_names = ['train','eval'],
    num_boost_round=1000,
    early_stopping_rounds=100,
    evals_result=evals_result,
)

train_metric = evals_result['train']['auc']
eval_metric = evals_result['eval']['auc']
plt.plot(train_metric, label='train auc')
plt.plot(eval_metric, label='eval auc')
plt.grid()
plt.legend()
plt.ylim(0, 1.1)

plt.xlabel('rounds')
plt.ylabel('auc')
plt.show()
