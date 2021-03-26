#!/usr/bin/python

# Code source: Gaël Varoquaux
#              Andreas Müller
# Modified for documentation by Jaques Grobler
# License: BSD 3 clause
import os
import sys
import math
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

#from matplotlib.colors import ListedColormap
#from joblib import dump, load

import xgboost as xgb

import ROOT
from ROOT import TH1, TFile

#MFT_X,MFT_Y,MFT_Phi,MFT_Tanl,MFT_InvQPt,MFT_Cov00,MFT_Cov01,MFT_Cov11,MFT_Cov02,MFT_Cov12,MFT_Cov22,MFT_Cov03,MFT_Cov13,MFT_Cov23,MFT_Cov33,MFT_Cov04,MFT_Cov14,MFT_Cov24,MFT_Cov34,MFT_Cov44,MFT_Chi2,MFT_nClu,MCH_X,MCH_Y,MCH_Phi,MCH_Tanl,MCH_InvQPt,MCH_Cov00,MCH_Cov01,MCH_Cov11,MCH_Cov02,MCH_Cov12,MCH_Cov22,MCH_Cov03,MCH_Cov13,MCH_Cov23,MCH_Cov33,MCH_Cov04,MCH_Cov14,MCH_Cov24,MCH_Cov34,MCH_Cov44,Truth

print(os.environ['ML_TRAINING_FILE'])

source_csv      = os.environ['ML_TRAINING_FILE']

df = pd.read_csv(source_csv, header=0)

df['MFT_ReducedTrackChi2'] = df['MFT_Chi2']/df['MFT_nClu']
df['Delta_X']      = df['MCH_X']      -  df['MFT_X']
df['Delta_Y']      = df['MCH_Y']      -  df['MFT_Y']
df['Delta_Phi']    = df['MCH_Phi']    -  df['MFT_Phi']
df['Delta_Tanl']   = df['MCH_Tanl']   -  df['MFT_Tanl']
df['Delta_InvQPt'] = df['MCH_InvQPt'] -  df['MFT_InvQPt']
df['Delta_XY']     = np.sqrt((df['MCH_X']-df['MFT_X'])**2 + (df['MCH_Y']-df['MFT_Y'])**2)

df.loc[df['Delta_Phi'] < -math.pi, 'Delta_Phi'] = df['Delta_Phi'] + 2*math.pi
df.loc[df['Delta_Phi'] >  math.pi, 'Delta_Phi'] = df['Delta_Phi'] - 2*math.pi

target_param        = 'Truth'
training_param_list = ['MFT_ReducedTrackChi2','Delta_X','Delta_Y','Delta_Phi','Delta_Tanl','Delta_InvQPt','Delta_XY']

xgb_params = {
    'objective': 'binary:logistic',
    'learning_rate': 0.01,
    'verbosity' : 2,
    'max_depth' : 6,
    'eval_metric': 'auc',
}
                
X = df.loc[:,training_param_list]
y = df.loc[:,target_param]

X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=0.2,random_state=0,stratify=y)
X_train, X_eval, y_train, y_eval = train_test_split(X_train, y_train,test_size=0.2,random_state=1,stratify=y_train)

xgb_train = xgb.DMatrix(X_train, label=y_train)
xgb_eval  = xgb.DMatrix(X_eval, label=y_eval)
xgb_test  = xgb.DMatrix(X_test, label=y_test)

evals = [(xgb_train, 'train'), (xgb_eval, 'eval')] 

evaluation_results = {}
bst = xgb.train(xgb_params,
                xgb_train,
                num_boost_round=1000,
                early_stopping_rounds=100,
                evals=evals,
                evals_result=evaluation_results,
                verbose_eval=10
                )

pickle.dump(bst, open("xgb_model.pickle", "wb"))

#gbm.plot_importance(gbm)

plt.show()
