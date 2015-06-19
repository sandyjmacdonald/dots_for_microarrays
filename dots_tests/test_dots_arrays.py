# -*- coding: utf-8 -*-
#!/usr/bin/env python

from nose import with_setup
from nose.tools import assert_equals
import pandas as pd
import numpy as np
import glob
from pandas.util.testing import assert_frame_equal
from ..dots_backend.dots_arrays import read_array, Array, read_experiment, Experiment

## Reading in an array and pulling out everything we needs for the tests.
array_fn = 'dots_sample_data/treated_1.txt'
group, replicate = array_fn.split('/')[-1].split('.')[0].split('_')
f = open(array_fn, 'r')
lines = f.readlines()
values = {}
headers = lines[9].split('\t')

## Compile a list of indices for the various items that will be pulled out of the file.
for i, h in enumerate(headers):
	if h == 'FeatureNum':
		feature_num_ind = i
	elif h == 'ProbeName':
		probe_name_ind = i
	elif h == 'GeneName':
		gene_name_ind = i
	elif h == 'SystematicName':
		sys_name_ind = i
	elif h == 'Description':
		desc_name_ind = i
	elif h == 'gProcessedSignal':
		signal_ind = i
	elif h == 'ControlType':
		control_type_ind = i

## Split the lines up and use the indices to pull out the values and add them 
## to the values dictionary.
for l in lines[10:]:
	l = l.split('\t')
	if l[control_type_ind] == '0':
		values[l[feature_num_ind]] = { 'FeatureNum': l[feature_num_ind], \
									   'ProbeName': l[probe_name_ind], \
									   'GeneName': l[gene_name_ind], \
									   'SystematicName': l[sys_name_ind], \
									   'Description': l[desc_name_ind], \
									   'gProcessedSignal': l[signal_ind] }

## Tests for the Array class and associated attributes and methods ##

## Create an array instance.
array = read_array(array_fn, group, replicate)

def test_probenames():
	assert_equals(sorted(array.get_probenames()), 
				  sorted([feat['ProbeName'] for feat in values.values()]))

def test_genenames():
	assert_equals(sorted(array.get_genenames()), 
				  sorted([feat['GeneName'] for feat in values.values()]))

def test_systematicnames():
	assert_equals(sorted(array.get_systematicnames()), 
				  sorted([feat['SystematicName'] for feat in values.values()]))

def test_descriptions():
	assert_equals(sorted(array.get_descriptions()), 
				  sorted([feat['Description'] for feat in values.values()]))

def test_intensities():
	assert_equals(sorted(array.get_intensities()), 
				  sorted([float(feat['gProcessedSignal']) for feat in values.values()]))

def test_normalisation():
	intensities = [float(feat['gProcessedSignal']) if float(feat['gProcessedSignal']) > 1.0 else 1.0 for feat in values.values()]
	logged = np.log2(intensities)
	logged_shifted = logged - np.percentile(logged, 75)
	assert_equals(sorted(array.get_normalised_intensities()), 
				  sorted(logged_shifted))

def test_read_array():
	assert_equals(array.filename, array_fn)
	group, replicate = array_fn.split('/')[-1].split('.')[0].split('_')
	replicate = int(replicate)
	sampleid = group + '_' + str(replicate)
	assert_equals(array.group, group)
	assert_equals(array.replicate, replicate)
	assert_equals(array.sampleid, sampleid)
	assert(isinstance(array, Array))
	assert(isinstance(array.df, pd.DataFrame))

## Tests for the Experiment class and associated attributes and methods ##

## Read in all of the arrays and create a couple of Experiment instances to test.
array_filenames = glob.glob('dots_sample_data/*.txt')
sampleids = [fn.split('/')[-1].split('.')[0] for fn in array_filenames]
experiment = read_experiment(array_filenames)
experiment_2 = read_experiment(array_filenames, baseline=False)

def test_read_experiment():
	assert(isinstance(experiment, Experiment))
	assert_equals(len(experiment.get_sampleids()), len(array_filenames))

def tests_arrays_attribute():
	assert_frame_equal(experiment.arrays[0].df.sort(axis=1), array.df.sort(axis=1), check_names=True)
	assert_equals(sorted(experiment.get_sampleids()), sorted(sampleids))

def test_baseline_to_median_method():
	assert_equals(experiment.baseline_to_median().has_baseline_to_median(), True)

def test_get_exp_values():
	rows, cols = experiment.get_exp_values().shape
	assert_equals(rows, len(array.intensities))
	assert_equals(cols, len(array_filenames))

def test_get_groups():
	groups = list(set(sorted([array_fn.split('/')[-1].split('.')[0].split('_')[0] for array_fn in array_filenames])))
	assert_equals(groups, experiment.get_groups())

def test_get_sample_ids():
	assert_equals(sorted(experiment.get_sampleids()), sorted(sampleids))

def test_remove_sample():
	assert_equals(len(experiment_2.remove_sample(array.sampleid).get_sampleids()), len(sampleids) - 1)