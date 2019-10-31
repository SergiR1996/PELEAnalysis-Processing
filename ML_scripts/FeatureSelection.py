# -*- coding: utf-8 -*-

# Imports
from sklearn.svm import SVC, SVR
import warnings # Use to ignore warnings
warnings.filterwarnings("ignore")
import os,sys
import pandas as pd
import numpy as np
import scipy, pickle
from sklearn.feature_selection import RFE, RFECV
from sklearn.externals import joblib
from sklearn.base import BaseEstimator, TransformerMixin
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("tkagg")

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"


class SequentianElimination(BaseEstimator, TransformerMixin):
	"""
	Class that is used to select the best features of a training set
	by using Sequential Feature Selector (SFS).

	PARAMETERS
    ----------
    nfeatures : integer
              the final number of features that want to be selected
    forward_floating : string
              name of the SFS algorithm type recognized by mlxtend
	"""

	def __init__(self, nfeatures, forward_floating):

		self.__kernel = "rbf"
		self.__nfeatures = nfeatures
		self.__forward_floating = forward_floating
		self.__scoring = "accuracy"
		self.__cv = 10


	def fit(self, X,y=None):

		return self

	def transform(self, X, y=None):

		if self.__forward_floating == "SFS":
			forward = True
			floating = False 

		elif self.__forward_floating == "SBS":
			forward = False
			floating = False

		elif self.__forward_floating == "SFFS":
			forward = True
			floating = True

		elif self.__forward_floating == "SBFS":
			forward = False
			floating = True

		sfs = SFS(SVC(kernel = self.__kernel), 
			k_features=self.__nfeatures, 
			forward=forward, 
			floating=floating, 
			verbose=2,
			scoring=self.__scoring,
			cv=self.__cv, n_jobs=4)

		sfs.fit(X,y)
		
		new_columns = list(sfs.k_feature_names_)

		X_final = X.copy()
		for elem in X.columns.values:
			if elem not in new_columns:
				X_final = X_final.drop([elem], axis=1)

		path = os.path.join(os.getcwd(), "ML_results")
		outfile_columns = open(os.path.join(path, "columns_new.pkl"), "wb")
		pickle.dump(X_final.columns.values, outfile_columns)
		df_sfs = pd.DataFrame(X_final, columns=X_final.columns.values)

		print("\nnew {}   length: {}\n".format(new_columns, len(new_columns)))
		
		return df_sfs

class RecursiveFeatureElimination(BaseEstimator, TransformerMixin):
	"""
	Class that is used to select the best features of a training set
	by using Recursive Feature Elimination (RFE).

	PARAMETERS
    ----------
    nfeatures : integer
              the final number of features that want to be selected
	"""

	def __init__(self,nfeatures):

		self.__estimator = SVC(kernel="linear")
		self.__num_features = nfeatures
		self._path = os.path.join(os.getcwd(), "ML_results")


	def fit(self, X,y=None):

		 return self

	def transform(self, X, y=None):

		selector = RFECV(self.__estimator, step=1,verbose=2,min_features_to_select=self.__num_features)
		df_rfe = selector.fit_transform(X,y)
		ranking = selector.ranking_
		new_columns = []
		for rank,feat in zip(ranking, X.columns.values):
			if rank == 1:
				new_columns.append(feat)

		outfile_columns = open(os.path.join(self._path, "columns_new.pkl"), "wb")
		pickle.dump(new_columns, outfile_columns)
		df_rfe = pd.DataFrame(df_rfe,columns=new_columns)

		print("\nnew {}   length: {}\n".format(new_columns, len(new_columns)))
		
		plt.figure()
		plt.xlabel("Number of features selected")
		plt.ylabel("Cross validation score (nb of correct classifications)")
		plt.plot(range(1, len(selector.grid_scores_) + 1), selector.grid_scores_)
		plt.show()

		return df_rfe 
