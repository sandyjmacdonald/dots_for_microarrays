# -*- coding: utf-8 -*-
#!/usr/bin/env python

from nose import with_setup
from nose.tools import assert_equals
import pandas as pd
import numpy as np
import glob
from itertools import combinations
from ..dots_backend.dots_arrays import read_experiment
from ..dots_backend.dots_analysis import run_pca, get_fold_changes, run_stats, find_clusters, get_clusters

## Read in all of the arrays and create an Experiment instance to test.
array_filenames = glob.glob('dots_sample_data/*.txt')
sampleids = [fn.split('/')[-1].split('.')[0] for fn in array_filenames]
experiment = read_experiment(array_filenames)

def test_run_pca():
	num_samples = len(experiment.get_sampleids())
	pca = run_pca(experiment)
	assert(isinstance(pca, pd.DataFrame))
	assert_equals(pca.shape, (num_samples, 4))

def test_get_fold_changes():
	fcs = get_fold_changes(experiment)
	num_rows = len(experiment.arrays[0].get_probenames())
	num_cols = 5 + len(experiment.groups) + (3 * len(list(combinations(experiment.groups, 2))))
	assert(isinstance(fcs, pd.DataFrame))
	assert_equals(fcs.shape, (num_rows, num_cols))

def test_run_stats():
	stats = run_stats(experiment)
	num_rows = len(experiment.arrays[0].get_probenames())
	num_cols = 7
	assert(isinstance(stats, pd.DataFrame))
	assert_equals(stats.shape, (num_rows, num_cols))

def test_find_clusters():
	df = experiment.df[experiment.get_sampleids()]
	num_samples = 1000
	rows = np.random.choice(df.index.values, num_samples)
	sampled_df = df.ix[rows]
	hier = find_clusters(sampled_df, how='hierarchical')
	kmeans = find_clusters(sampled_df, k_range=(3,6), how='kmeans')
	assert_equals(len(hier), num_samples)
	assert_equals(len(kmeans), num_samples)

def test_get_clusters():
	cluster_df = get_clusters(experiment)
	print cluster_df.columns.values
	assert(isinstance(cluster_df, pd.DataFrame))
	num_cols = 9 + (2 * len(list(combinations(experiment.groups, 2)))) + len(experiment.get_sampleids())
	assert_equals(len(cluster_df.columns.values), num_cols)