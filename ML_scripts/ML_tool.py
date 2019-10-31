# -*- coding: utf-8 -*-

# Imports
from sklearn.svm import SVC, SVR
import os,sys
import argparse as ap
import pandas as pd

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV, RandomizedSearchCV,\
										 StratifiedKFold, cross_val_predict, cross_validate, learning_curve, GroupKFold, ShuffleSplit
import numpy as np
import scipy, pickle
from sklearn.feature_selection import RFE, RFECV
from sklearn.externals import joblib
from sklearn.metrics import recall_score, confusion_matrix, classification_report
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("tkagg")
from mlxtend.plotting import plot_decision_regions
from FeatureSelection import *

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"


class Dataset(object):
	"""
	Class that opens a csv file and splits it in train, validation, and 
	test sets.

	PARAMETERS
    ----------
    csv : string
              name of the csv file (dataset)
    test_size : float
              proportion of the total dataset that will be used as validation and test sets
    to_drop : list of strings
              list of column names that want to be removed from the dataset (which could
              include the labels or other useless columns)
    target : string
              name of the column that will act as label
	"""

	def __init__(self, csv, test_size, to_drop, target):

		self.__csv = csv
		self._df = self.__DataFrame()
		self._to_drop = to_drop
		self._target = target
		self.__test_size = test_size
		self._random_state = 5
		self._path = self.__CreateTempDirectory()
		self._train, self._test = self.__SplitDataset()
		self._X_train, self._y_train = self.__SplitTrain()
		self.__SaveTest()

	def __len__(self):

		return self._df.shape[0]

	def __str__(self):

		return "{}".format(self._df)

	def __getitem__(self):

		pass

	def __DataFrame(self):

		return pd.read_csv(self.__csv, sep=",")

	def __SplitDataset(self):

		return train_test_split(self._df, test_size=self.__test_size, random_state=self._random_state)

	def __SplitTrain(self):


		return self._train.drop(self._to_drop, axis=1), self._train[self._target]

	def __SaveTest(self):

		val, test = train_test_split(self._test, train_size=0.40, random_state=self._random_state)
		outfile_val = open(os.path.join(self._path, "val.pkl"), "wb")
		outfile_test = open(os.path.join(self._path, "test.pkl"), "wb")

		if sys.version_info[0]==3:
			pickle.dump(val, outfile_val, protocol=3)
			pickle.dump(test, outfile_test, protocol=3)
		else:
			pickle.dump(val, outfile_val, protocol=2)
			pickle.dump(test, outfile_test, protocol=2)

		outfile_val.close(); outfile_test.close()


	def __CreateTempDirectory(self):

		PATH = os.path.join(os.getcwd(), "ML_results")
		if not os.path.exists(PATH): os.mkdir(PATH)
		else: pass 

		return PATH


	@property
	def RandomState(self):
		return self._random_state

	@RandomState.setter
	def RandomState(self, value):
		self._random_state = value

	@property
	def TestSize(self):
		return self.__test_size

	@TestSize.setter
	def TestSize(self, value):
		self.__test_size = value

	@property
	def TrainSet(self):
		return self._train

	@property
	def TestSet(self):
		return self._test

class TransformDataset(Dataset):
	"""
	Subclass of the Dataset class that uses the methods from the previous class
	and implements the scaler that will be applied to the columns of the train
	set.

	PARAMETERS
    ----------
    scaler : string
              name of the scaler recognized by sci-kit learn
	"""
	
	def __init__(self, csv, to_drop, target, scaler, test_size):

		super(TransformDataset, self).__init__(csv, test_size, to_drop, target)
		self.__scaler_type = scaler
		self._X_transformed = self.__Scale()

	def __Scale(self):

		if self.__scaler_type == "std":
			scaler = StandardScaler().fit(self._X_train)
		elif self.__scaler_type == "minmax":
			scaler = MinMaxScaler().fit(self._X_train)
		elif self.__scaler_type == "robust":
			scaler = RobustScaler().fit(self._X_train)

		outfile = open(os.path.join(self._path, "scaler_{}.pkl".format(self.__scaler_type)), "wb")
		pickle.dump(scaler, outfile)
		outfile.close()

		return scaler.transform(self._X_train)

	@property
	def scaler_type(self):
		return self.__scaler_type
	

class CreateModel(TransformDataset):
	"""
	Subclass of the TransformDataset class that trains the model against the train set
	and also can be validated against both holdout sets.

	PARAMETERS
    ----------
    scaler : string
              name of the scaler recognized by sci-kit learn
	"""

	def __init__(self,to_drop):

		csv, hey, target, nfeatures, scaler, test_size, n_iter, n_jobs, forward_floating, train, rfe, nofeature, test, model_file = self.parseArgs()
		drop = to_drop

		super(CreateModel, self).__init__(csv, to_drop, target, scaler, test_size)
		self.__verbose = 10
		self.__n_jobs=n_jobs
		self.__n_iter = n_iter
		self.__nfeatures = nfeatures
		self.__forward_floating = forward_floating
		self.__train = train
		self.__rfe = rfe
		self.__nofeature = nofeature
		self.__test = test
		self.__model_file = model_file

		if self.__train:

			self._X_transformed = pd.DataFrame(self._X_transformed,columns=self._X_train.columns.values)
			if self.__rfe:
				rfe = RecursiveFeatureElimination(self.__nfeatures)
				rfe.fit(self._X_transformed, self._y_train)
				self.__sf = rfe.transform(self._X_transformed, self._y_train)
			elif self.__nofeature:
				self.__sf = self._X_transformed
			else:
				sfs = SequentianElimination(self.__nfeatures,self.__forward_floating)
				sfs.fit(self._X_transformed, self._y_train)
				self.__sf = sfs.transform(self._X_transformed, self._y_train)

	@property
	def train(self):
		return self.__train

	def parseArgs(self):
		"""
		Parse arguments from command-line
		"""

		parser = ap.ArgumentParser(description='Script used to perform supervised learning with your dataset\
			against a label that can be binary or multiple')
		optional = parser._action_groups.pop()
		parser.add_argument("csv", metavar="FILE",type=str, help="Path of the csv file (dataset)")
		parser.add_argument("target", metavar="Column name",type=str, help="Column name of the label")
		optional.add_argument("-td","--to_drop",metavar="LIST",nargs='*',help="Column names that want to be removed",default=[])
		optional.add_argument("-nf","--nfeatures",metavar="INTEGER",type=int,help="Final number of features for the feature selector algorithm",default=30)
		optional.add_argument("-i","--iterations",metavar="INTEGER",type=int,help="Number of iterations for the cross validation training",default=10000)
		optional.add_argument("-j","--n_jobs",metavar="INTEGER",type=int,help="Number of jobs for the cross validation training",default=10)
		optional.add_argument("-ts","--test_size",metavar="FLOAT",type=float,help="Relative size of the holdout sets comparing with the train set",default=0.3)
		optional.add_argument("-sfs","--SequentialFeatureSelection",metavar="STRING",type=str,help="SFS algorithm type",default="SFS")
		optional.add_argument("-sc","--scaler",metavar="STRING",type=str,help="Scaler type",default="std")
		optional.add_argument("-T", "--train", help="Perform the cross-validation training", action="store_true")
		optional.add_argument("-R", "--rfe", help="Use RecursiveFeatureElimination for feature selection", action="store_true")
		optional.add_argument("-NS", "--nofeature", help="Do not perform feature selection", action="store_true")
		optional.add_argument("-TE", "--test", help="Use the test set with the generated model", action="store_true")		
		optional.add_argument("-mf","--model_file",metavar="FILE",type=str,help="Path of the model pickle file",default="")		
		parser._action_groups.append(optional)
		args = parser.parse_args()

		return args.csv, args.to_drop, args.target, args.nfeatures, args.scaler, args.test_size, args.iterations, args.n_jobs, args.SequentialFeatureSelection, args.train, args.rfe, args.nofeature, args.test, args.model_file

	def __RandomSearch(self):

		param_grid = {"C": scipy.stats.uniform(loc=0, scale=4), "gamma":scipy.stats.expon(scale=.1), 
																"kernel":["linear","rbf","sigmoid"], "class_weight" : ["balanced", None]}
		grid = RandomizedSearchCV(estimator=SVC(), param_distributions=param_grid, 			
													n_jobs=self.__n_jobs,  cv=ShuffleSplit(n_splits=10, test_size = .3, random_state = self._random_state), verbose=self.__verbose, n_iter=self.__n_iter) 
		grid_result = grid.fit(self.__sf, self._y_train)

		return grid_result.best_params_

	def plot(self, estimator, X, y, title="learning curves", ylim=None, cv=None,n_jobs=None, train_sizes=np.linspace(.1, 1.0, 5)):

		plt.figure()
		plt.title(title)
		if ylim is not None:
		    plt.ylim(*ylim)
		plt.xlabel("Training examples")
		plt.ylabel("Score")
		train_sizes, train_scores, test_scores = learning_curve(
		    estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
		train_scores_mean = np.mean(train_scores, axis=1)
		train_scores_std = np.std(train_scores, axis=1)
		test_scores_mean = np.mean(test_scores, axis=1)
		test_scores_std = np.std(test_scores, axis=1)
		plt.grid()

		plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
							train_scores_mean + train_scores_std, alpha=0.1,
							color="r")
		plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
         					test_scores_mean + test_scores_std, alpha=0.1, color="g")
		plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
				label="Training score")
		plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
			label="Cross-validation score")
		plt.legend(loc="best")
		plt.show()

	def CreateSVCModel(self, random_search=True, cross_val=True):

		if random_search:
			best_params = self.__RandomSearch()
			model = SVC(C=best_params["C"], kernel=best_params["kernel"], gamma=best_params["gamma"])
			model.fit(self.__sf, self._y_train)
			joblib.dump(model, os.path.join(self._path, "svm_{}_clf_{}_{}_{}_{}_{}.pkl".format(self._target,best_params["kernel"], best_params["C"], best_params["gamma"], self.__nfeatures , self.__forward_floating)))
			if cross_val:
				scores = cross_val_score(model, self.__sf, self._y_train, cv=ShuffleSplit(n_splits=10, test_size = .3, random_state = self._random_state))
				self.plot(SVC(C=best_params["C"], gamma=best_params["gamma"], kernel=best_params["kernel"], class_weight=best_params["class_weight"]), self.__sf, self._y_train)
				print("scores: {}\n mean_score: {}(+/-{})\n best_params: {}".format(scores, scores.mean(),scores.std(), best_params))
				outfile = open(os.path.join(self._path, "svm_{}_clf_{}_{}_{}_{}_{}.txt".format(self._target,best_params["kernel"], best_params["C"], best_params["gamma"],
					self.__nfeatures,self.__forward_floating)), "w")
				outfile.write("scores: {}\n mean_score: {}(+/-{})\n best_params: {}".format(scores, scores.mean(),scores.std(), best_params)); outfile.close()

				return scores.mean()

	def CheckTest(self):


		if self.__test:
			infile = open(os.path.join(self._path, "test.pkl"), "rb")
		else:
			infile = open(os.path.join(self._path, "val.pkl"), "rb")

		test = pickle.load(infile)
		infile.close()
		X_test, y_test = test.drop(self._to_drop, axis=1), test[self._target]

		scale_file = open(os.path.join(self._path, "scaler_{}.pkl".format(self.scaler_type)), "rb")
		scaler = pickle.load(scale_file)
		scale_file.close()

		X_test_scaled = scaler.transform(X_test)
		X_test_scaled = pd.DataFrame(X_test_scaled, columns=X_test.columns.values)

		if not self.__nofeature:
			columns_file = open(os.path.join(self._path, "columns_new.pkl"), "rb")
			columns = pickle.load(columns_file)
			columns_file.close()
			for elem in X_test_scaled.columns.values:
				if elem not in columns:
					X_test_scaled = X_test_scaled.drop([elem],axis=1)

		model = joblib.load(self.__model_file)
		y_pred = model.predict(X_test_scaled)

		df = pd.DataFrame(columns=["esterase","predicted", "test","exp"])
		df["esterase"] = list(test["esterase"]); df["predicted"] = list(y_pred)
		df["test"] = list(y_test); df["activity"] = list(test["exp"])

		if self.__test:
			df.to_csv("check_test.csv", sep=",")
		else:
			df.to_csv("check_val.csv", sep=",")

		print("Confusion matrix: \n", confusion_matrix(y_test, y_pred))
		print("report: \n", classification_report(y_test, y_pred))

		# Final_df = pd.DataFrame(X_test_scaled, columns=X_test_scaled.columns.values)
		# Values = {}
		# Ranges = {}
		# for i in range(len(Final_df.columns)):
		# 	if i>=2:
		# 		mean = np.mean(Final_df.iloc[:,i])
		# 		Values[i]=abs(mean)
		# 		if mean!=0:
		# 			Ranges[i]=abs(np.std(Final_df.iloc[:,i])/mean)
		# 		else:
		# 			Ranges[i]=1

		# plot_decision_regions(X_test_scaled,np.array(y_test),clf=model,filler_feature_values=Values,filler_feature_ranges=Ranges, legend = 2)
		# plt.show()


if __name__=="__main__":
	#to_drop=["activity","clf","Unnamed: 0","esterase","ID enzyme","compounds"]
	to_drop=["hydrolized","esterase","exp","multi","active_site_type",
	"min_lap_C_C","min_lap_C_O","min_lap_C_N","min_lap_C_H","min_lap_O_O",
	"min_lap_O_N","min_lap_O_H","min_lap_N_N","min_lap_N_H","min_lap_H_H",
	"mean_adj_C_C","mean_adj_C_O","mean_adj_C_N","mean_adj_C_H","mean_adj_O_O",
	"mean_adj_O_N","mean_adj_O_H","mean_adj_N_N","mean_adj_N_H","mean_adj_H_H",
	"min_adj_O_O","min_adj_H_H","min_adj_N_N"]
	model = CreateModel(to_drop)
	if model.train:
		model.CreateSVCModel()
	else:
		model.CheckTest()