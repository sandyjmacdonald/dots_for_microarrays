# -*- coding: utf-8 -*-
#!/usr/bin/env python

import warnings
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as hac
from dots_arrays import Experiment
from sklearn.decomposition import PCA
from itertools import combinations
from scipy.stats import ttest_ind, f_oneway
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.multicomp import MultiComparison
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.cluster import KMeans

## Functions ##

def run_pca(experiment):
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

def get_fold_changes(experiment):
	'''Calculate pairwise fold change and log fold change values.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	Returns:
		A new Pandas data frame with pairwise fold change and log fold change values.

	'''

	groups = experiment.get_groups()
	pairs = map(list, list(combinations(groups, 2)))
	if all([g.isdigit() for g in groups]):
		pairs = sorted(pairs, key=lambda x:x[0])
	samples = experiment.get_sampleids()
	df = experiment.df
	new_df = df.ix[:, :5].copy()

	for group in groups:
		ids = [sample for sample in samples if group == sample.split('_')[0]]
		new_df['mean_' + group] = df[ids].mean(axis=1)

	del df

	## For each pair, calculate mean values for each group, fold changes and log2 fold changes.
	for pair in pairs:
		if all([g.isdigit() for g in pair]):
			pair.sort(key=int, reverse=True)
		else:
			pair.sort()
		name_1, name_2 = pair
		new_df['abs_mean_diff_' + name_1 + '_' + name_2] = abs((2 ** new_df['mean_' + name_1]) - (2 ** new_df['mean_' + name_2]))
		new_df['logFC_' + name_1 + '_' + name_2] = new_df['mean_' + name_1] - new_df['mean_' + name_2]
		new_df['FC_' + name_1 + '_' + name_2] = 2 ** new_df['logFC_' + name_1 + '_' + name_2]

	return new_df

def run_stats(experiment):
	'''Run independent T-test or one-way ANOVA dependent on number of groups.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	Returns:
		A new Pandas data frame with p values, adjusted p values and Tukey HSD
		post-hoc results if there are > 2 groups.

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
	p_val_adj = list(multipletests(p_vals, method='fdr_bh')[1])
	new_df = df.ix[:, :5].copy()
	new_df['p_val'] = pd.Series(p_vals, index=new_df.index)
	new_df['p_val_adj'] = pd.Series(p_val_adj, index=new_df.index)

	## Post-hoc test.

	## Only do the post-hoc test if there are more than 2 groups, duh!
	if len(groups) > 2:
		vals_df = df[samples]
		group_ids = [sample.split('_')[0] for sample in vals_df.columns.values]
		posthoc_results = {}

		## Run the post-hoc test on each row.
		for row in range(len(vals_df)):
			row_vals = vals_df.ix[row]
			mc = MultiComparison(row_vals, group_ids)
			mc_groups = mc.groupsunique
			results = mc.tukeyhsd()
			significant = results.reject
			pairs = zip(*[x.tolist() for x in mc.pairindices])

			## Go through each pair and add results to the posthoc_results dictionary.
			for i in range(len(pairs)):
				pair = list(pairs[i])
				pair.sort()
				pair_name = str(mc_groups[pair[0]]) + '_' + str(mc_groups[pair[1]])
				if pair_name in posthoc_results:
					posthoc_results[pair_name].append(significant[i])
				else:
					posthoc_results[pair_name] = [significant[i]]
		
		## Add the post-hoc results to the data frame.
		for pair_name in posthoc_results:
			new_df['significant_' + pair_name] = posthoc_results[pair_name]

	return new_df

def find_clusters(df, k_range=(3,11), how='hierarchical'):
	'''Find clusters, and if method is k-means run silhouette analysis
	to determine the value of k.

	Args:
		df (data frame): A data frame with normalised expression data.
		k_range (tuple): The range over which to test k.
		how ('hierarchical' or 'kmeans'): Clustering method.

	Returns:
		A list of cluster numbers.

	'''

	## Don't run the silhouette analysis for hierarchical clustering,
	## just calculate the clusters using estimate of k.
	if how == 'hierarchical':
		k = int(np.sqrt((len(df) / 2.0)))
		hc = hac.linkage(df, method='average')
		optimal_clusters = hac.fcluster(hc, t=k, criterion='maxclust')

	## If method is k-means, run silhouette analysis.
	elif how == 'kmeans':
		best_combined_score = 0
		optimal_k = 2
		
		## Try values of k from range and keep track of optimal k according
		## to silhouette score.
		for k in range(*k_range):
			km = KMeans(n_clusters=k, random_state=10)
			clusters = km.fit_predict(df)
			silhouette_avg = silhouette_score(df, clusters)
			sample_silhouette_values = silhouette_samples(df, clusters)
			above_mean = 0
			silhouette_sizes = []

			for i in range(k):
				ith_cluster_silhouette_values = sample_silhouette_values[clusters == i]
				size_cluster_i = ith_cluster_silhouette_values.shape[0]
				silhouette_sizes.append(size_cluster_i)
				if max(ith_cluster_silhouette_values) > silhouette_avg:
					above_mean += 1
			
			## This combined score should pick the best value of k
			above_mean_score = float(above_mean) / k
			std_score = 1.0/np.std(silhouette_sizes) if np.std(silhouette_sizes) > 1.0 else 1.0
			combined_score = (silhouette_avg + above_mean_score + std_score) / 3
			
			## Put the clusters in the new column in the data frame.
			if combined_score > best_combined_score:
				best_combined_score = combined_score
				optimal_k = k
				optimal_clusters = clusters

		optimal_clusters = [cluster + 1 for cluster in optimal_clusters]

	return optimal_clusters

def get_clusters(experiment, how='hierarchical'):
	'''Clusters significantly differentially expressed genes by expression pattern 
	across the samples using hierarchical or k-means clustering and silhouette analysis 
	to pick the value of k (via the find_clusters function).

	Args:
		experiment (Experiment instance): An instance of the Experiment class.
		how ('hierarchical' or 'kmeans'): Clustering method.

	Returns:
		A new Pandas data frame with fold changes, p values and clusters.

	'''

	## Run the stats to filter genes down to significant ones only.
	stats = run_stats(experiment)
	stats = stats[['FeatureNum', 'p_val', 'p_val_adj']].copy()

	## Get the fold changes 
	fcs = get_fold_changes(experiment)	
	keep_cols = [x for x in fcs.columns.values if 'logFC' in x or 'abs_mean_diff' in x]
	fc_cols = [x for x in fcs.columns.values if 'logFC' in x]
	fcs = fcs[['FeatureNum'] + keep_cols].copy()
	
	norm_exp_cols = experiment.get_sampleids()

	abs_mean_diff_cols = [x for x in fcs.columns.values if 'abs_mean_diff' in x]

	## Merge together the stats and fold changes data frames.
	merged_df = pd.merge(experiment.df, stats, on='FeatureNum')
	merged_df = pd.merge(merged_df, fcs, on='FeatureNum')
	
	## Filter the merged data frame to leave only significantly differentially
	## expressed genes (adj. p < 0.05, fold change > +/- 2).
	filtered_df = merged_df[(merged_df['p_val_adj'] < 0.05) & ((abs(merged_df[fc_cols]) > np.log2(1.0)).any(1) == True) & ((merged_df[abs_mean_diff_cols] > 0.5).any(1) == True)].copy()

	## Clean up.
	del merged_df
	del stats
	del fcs

	## A good guesstimate for k.
	k_limit = int(np.sqrt((len(filtered_df) / 2)))

	## Catches numpy warnings about means of empty slices.
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", category=RuntimeWarning)

		## Hierarchical clustering.
		if how == 'hierarchical':
			clusters = find_clusters(filtered_df[norm_exp_cols], how='hierarchical')
			filtered_df['cluster'] = clusters

		## K-means clustering with silhouette analysis to determine value of k.
		elif how == 'kmeans':
			clusters = find_clusters(filtered_df[norm_exp_cols], k_range=(3, k_limit), how='kmeans')
			filtered_df['cluster'] = clusters
	
	## Sort the data frame by cluster and mean expression across samples.
	filtered_df['mean_norm_expression'] = filtered_df[norm_exp_cols].mean(axis=0)
	filtered_df.sort(columns=['cluster', 'mean_norm_expression'], ascending=[True, False], inplace=True)
	filtered_df = filtered_df.reset_index(drop=True)

	return filtered_df

def write_fcs_stats(experiment, outfile='foldchanges_stats.txt'):
	'''Creates a tab-separated table with a full list of fold changes,
	p values, adjusted p values and post hoc results.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.
		outfile (string): The name of the table-separated table to be created.

	'''

	## Run the stats and fold changes and merge them into a single data frame.
	stats = run_stats(experiment)
	posthoc_cols = [colname for colname in stats.columns.values if 'significant' in colname]
	stats = stats[['FeatureNum', 'p_val', 'p_val_adj'] + posthoc_cols]
	fcs = get_fold_changes(experiment)
	fc_cols = [colname for colname in fcs.columns.values if not 'abs_mean_diff_' in colname]
	merged_df = pd.merge(fcs, stats, on='FeatureNum')
	
	## Define the order of the columns in the data frame.
	colnames = list(merged_df.columns.values)
	global col_order
	col_order = ['mean', 'FC', 'logFC', 'abs_mean_diff', 'p_val', 'adj_p_val', 'significant']

	## Function to custom sort the columns.
	def keyfunc(col):
		for c in col_order:
			if col.startswith(c):
				return (col_order.index(c), col.lstrip(c + '_'))
				break

	## Sort the columns.
	sorted_colnames = colnames[:5] + sorted(colnames[5:], key=keyfunc)
	merged_df = merged_df[sorted_colnames]

	## Fix the type of the FeatureNum column and sort it.
	merged_df['FeatureNum'] = merged_df['FeatureNum'].astype(int)
	merged_df.sort(columns='FeatureNum', ascending=True, inplace=True)

	## Write the table.
	merged_df.to_csv(outfile, sep='\t', index=False)

def write_normalised_expression(experiment, outfile='normalised_expression.txt'):
	'''Creates a tab-separated table with all of the normalised expression values.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.
		outfile (string): The name of the table-separated table to be created.

	'''

	## Read in the experiment.
	experiment_df = experiment.df

	## Sort the values columns.
	colnames = list(experiment_df.columns.values)
	sorted_colnames = colnames[:5] + sorted(colnames[5:])
	experiment_df = experiment_df[sorted_colnames]

	## Write the table.
	experiment_df.to_csv(outfile, sep='\t', index=False)