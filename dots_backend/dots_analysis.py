# -*- coding: utf-8 -*-
#!/usr/bin/env python

import pandas as pd
import numpy as np
from dots_arrays import Experiment
from sklearn.decomposition import PCA
from itertools import combinations
from scipy.stats import ttest_ind, f_oneway
from statsmodels.sandbox.stats.multicomp import multipletests

## Functions ##

def do_pca(experiment):
	'''Run PCA when given an experiment instance or data frame with expression values.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	Returns:
		A Pandas data frame with results of PCA analysis.

	'''
	
	## The below if/elif checks whether experiment passed to function is an instance of the
	## Experiment class or just a data frame with expression values in.
	if isinstance(experiment, Experiment):
		df = experiment.get_exp_values().T
	elif isinstance(experiment, pd.DataFrame):
		df = experiment.T

	## Run the PCA, get the scores and unzip tuples into separate lists of x and y values.
	pca = PCA(n_components=3)
	pca_fit = pca.fit_transform(df)
	vals = [(x[0], x[1]) for x in pca_fit]
	xvals, yvals = zip(*vals)

	## Convert the data into a dictionary for easy conversion into a Pandas data frame.
	pca_dict = {'xvals': xvals, 'yvals': yvals, 'sampleid': list(df.index), 'group': [x.split('_')[0] for x in list(df.index)]}
	pca_df = pd.DataFrame(pca_dict)

	return pca_df

def do_fold_changes(experiment):
	'''Calculate pairwise fold change and log fold change values.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	Returns:
		A new Pandas data frame with pairwise fold change and log fold change values.

	'''
	groups = experiment.get_groups()
	pairs = list(combinations(groups, 2))
	samples = experiment.get_sampleids()
	df = experiment.df
	new_df = df.ix[:, :5].copy()

	## For each pair, calculate mean values for each group, fold changes and log2 fold changes.
	for pair in pairs:
		name_1, name_2 = pair
		ids_1 = [sample for sample in samples if name_1 == sample.split('_')[0]]
		mean_1 = df[ids_1].mean(axis=1)
		ids_2 = [sample for sample in samples if name_2 == sample.split('_')[0]]
		mean_2= df[ids_2].mean(axis=1)
		new_df['logFC_' + name_1 + '_' + name_2] = mean_1 / mean_2
		new_df['FC_' + name_1 + '_' + name_2] = 2 ** new_df['logFC_' + name_1 + '_' + name_2]
	
	return new_df

def run_stats(experiment):
	'''Run independent T-test or one-way ANOVA dependent on number of groups.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	Returns:
		A new Pandas data frame with p values and adjusted p values.

	'''
	groups = experiment.get_groups()
	samples = experiment.get_sampleids()
	df = experiment.df
	all_vals = []

	## Get values for each group, ready for T-test or ANOVA.
	for group in groups:
		ids = [sample for sample in samples if group == sample.split('_')[0]]
		vals = map(list, df[ids].values)
		all_vals.append(vals)
	
	## Decide whether to use T-test or ANOVA dependent on number of groups.
	if len(groups) == 2:
		p_vals = [ttest_ind(all_vals[0][i], all_vals[1][i])[1] for i in range(len(all_vals[0]))]
	else:
		p_vals = []
		for i in range(len(all_vals[0])):
			row_vals = [all_vals[j][i] for j in range(len(groups))]
			p_val = f_oneway(*row_vals)[1]
			p_vals.append(p_val)
	
	## Adjust the p values and create a new data frame with them in.
	adj_p_vals = list(multipletests(p_vals, method='fdr_bh')[1])
	new_df = df.ix[:, :5].copy()
	new_df['p_val'] = pd.Series(p_vals, index=new_df.index)
	new_df['adj_p_val'] = pd.Series(adj_p_vals, index=new_df.index)
	
	return new_df