# -*- coding: utf-8 -*-
#!/usr/bin/env python

import pandas as pd
from dots_arrays import Experiment
from sklearn.decomposition import PCA

## Functions ##

def do_pca(experiment):
	"""Run PCA when given an experiment instance or data frame with expression values."""
	
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