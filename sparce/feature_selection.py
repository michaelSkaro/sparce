import os
import re
import sys
import math
import seaborn as sns
import glob
import io
import math
import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import where
import sklearn
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import RFE, SelectFromModel, SelectKBest, chi2
from sklearn.linear_model import LogisticRegression

#from sparce import ars_parse

class feature_selection:

   def __init__(self, X, y, nFeatures,nJobs):
       self.X = X
       self.y = y
       self.nFeatures = nFeatures
       self.nJobs = nJobs
       pass
   
   
   def chiSQ_selector(X, y, nFeatures):
       '''
       :param X: continuous data
       :type X: pandas dataframe
       :param y: class data
       :type y: encoded class data
       
       :return: support
       :rtype: list boolean
       :return: feature_list 
       :rtype: list string
       
       Objective: Use the chi-squared test to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.chi2.html

       '''
       
       feature_list = X.columns
       chiSQ = SelectKBest(chi2, k= nFeatures)
       chiSQ.fit(X, y.ravel())
       chiSQ_support = chiSQ.get_support()

       return chiSQ_support, feature_list
   
   
   def rfR(X, y, nFeatures,nJobs):
       '''
       
       :param X: continuous data
       :type X: pandas dataframe
       :param y: class data
       :type y: encoded class data
       :param nFeatures: number of features to select
       :type nFeatures: int
       :param nJobs: number of jobs to run in parallel
       :type nJobs: int

       
       :return: support
       :rtype: list boolean
       :return: feature_list 
       :rtype: list string
       

       Objective: Use the RandomForestRegressor to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html


       '''


       rfR_features = SelectFromModel(
           RandomForestRegressor(n_estimators=100, n_jobs = nJobs), max_features=nFeatures
       )
       rfR_features.fit(X, y.ravel())
       rfR_support = rfR_features.get_support()

       return rfR_support
   
   
   def recursiveFeatureSelection(X, y, nFeatures,nJobs):
       '''
       :param X: continuous data
       :type X: pandas dataframe
       :param y: class data
       :type y: encoded class data
       :param nFeatures: number of features to select
       :type nFeatures: int
       :param nJobs: number of jobs to run in parallel
       :type nJobs: int
       
       :return: support
       :rtype: list boolean
       :return: feature_list 
       :rtype: list string

       Objective: Use the Recursive Feature Elimination to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html

       '''



       rfe_selector = RFE(
           estimator=LogisticRegression(n_jobs = nJobs),
           n_features_to_select=nFeatures,
           step= math.celi(len(X.columns)/10),
           verbose=0,
       )
       rfe_selector.fit(X, y.ravel())
       rfe_support = rfe_selector.get_support()

       return rfe_support
   
   
   
   def lassoR(X, y, nFeatures,nJobs):
       '''
       :param X: continuous data
       :type X: pandas dataframe
       :param y: class data
       :type y: encoded class data
       :param nFeatures: number of features to select
       :type nFeatures: int
       :param nJobs: number of jobs to run in parallel
       :type nJobs: int
       
       :return: support
       :rtype: list boolean
       :return: feature_list 
       :rtype: list string

       Objective: Use the Lasso Regression to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html
       
       '''



       lR_selector = SelectFromModel(
           LogisticRegression(penalty="l2", n_jobs= nJobs), max_features=nFeatures
       )
       lR_selector.fit(X, y.ravel())
       lR_support = lR_selector.get_support()

       return lR_support
   
   
   
   def rfC(X, y, nFeatures,nJobs):
       '''
       :param X: continuous data
       :type X: pandas dataframe
       :param y: class data
       :type y: encoded class data
       :param nFeatures: number of features to select
       :type nFeatures: int
       :param nJobs: number of jobs to run in parallel
       :type nJobs: int
       
       :return: support
       :rtype: list boolean
       :return: feature_list 
       :rtype: list string

       Objective: Use the RandomForestClassifier to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
       
       '''



       rfC_features = SelectFromModel(
           RandomForestClassifier(n_estimators=100, n_jobs=nJobs), max_features=nFeatures
       )
       rfC_features.fit(X, y.ravel())
       rfC_support = rfC_features.get_support()

       return rfC_support
   
   # cross validation
   
   def cross_validate_feature_selection(feature_list, chiSQ_support, rfe_support, lR_support, rfC_support,rfR_support):
        '''
        :param feature_list: list of features
        :type feature_list: list string
        :param chiSQ_support: support
        :type chiSQ_support: list boolean
        :param rfe_support: support
        :type rfe_support: list boolean
        :param lR_support: support
        :type lR_support: list boolean
        :param rfC_support: support
        :type rfC_support: list boolean
        :param rfR_support: support
        :type rfR_support: list boolean
        
        :return: dataframe quantifying the performance of each feature selection method
        :rtype: dataframe

        Objective: Use the cross validation to select features.
        Reference: https://towardsdatascience.com/feature-selection-techniques-in-machine-learning-with-python-f24e7da3f36e
        
        '''
        df = pd.DataFrame(
            {
                "Feature": feature_list,
                "chi2": chiSQ_support,
                "RFE": rfe_support,
                "Logistics": lR_support,
                "RandomForestClassifier": rfC_support,
                "RandomForstRegression": rfR_support,
            }
        )
        df["Total"] = np.sum(df, axis=1)
        df = df.sort_values(["Total", "Feature"], ascending=False)
        df.index = range(1, len(df) + 1)
        
        return df

   def grade_features(X, y, nFeatures, n_jobs):
       '''
       :param X: continuous data
         :type X: pandas dataframe
       :param y: class data
         :type y: encoded class data
       :param nFeatures: number of features to select
       :type nFeatures: int
       :param nJobs: number of jobs to run in parallel
       :type nJobs: int

       Objective: Call each function in the feature selection class, grade each set of features with each of the algorithms.
       Sum the total boolen supports for each feature and then sort the features by the sum.
       
       Reference: https://towardsdatascience.com/feature-selection-techniques-in-machine-learning-with-python-f24e7da3f36e

       '''
       chiSQ_support, feature_list = feature_selection.chiSQ_selector(X, y, nFeatures=nFeatures)
       rfe_support = feature_selection.recursiveFeatureSelection(X, y, nFeatures=nFeatures, nJobs = nJobs)
       lR_support = feature_selection.lassoR(X, y, nFeatures=nFeatures, nJobs = nJobs)
       rfC_support = feature_selection.rfC(X, y, nFeatures=nFeatures, nJobs = nJobs)
       rfR_support = feature_selection.rfR(X, y, nFeatures=nFeatures, nJobs = nJobs)

       CV = feature_selection.cross_validate_feature_selection(
           feature_list,
           chiSQ-support,
           rfe_support,
           lR_support,
           rfC_support,
           rfR_support,
       )
       CV = CV[1:nFeatures]

       return CV
