# -*- coding: utf-8 -*-
#!/usr/bin/env python

import pandas as pd
import numpy as np

## Classes ##

class Array:
	'''Holds an array and its associated data'''

	def __init__(self, df, filename, group, replicate):
		'''Set up all of the Array attributes.

		As well as the parameters passed, a number of other attributes are set up
		through class methods.

		Args:
			df (Pandas data frame): Dataframe with array data, preferably created with the read_array() function.
			filename (str): The filename of the raw array data file.
			group (str): The name of the group to which this array belongs.
			replicate (int): The replicate number within the group.

		'''

		self.df = df
		self.df['gProcessedSignal'] = self.df['gProcessedSignal'].astype(float)
		self.filename = filename
		self.probenames = self.get_probenames()
		self.genenames = self.get_genenames()
		self.systematicnames = self.get_systematicnames()
		self.descriptions = self.get_descriptions()
		self.intensities = self.get_intensities()
		self.group = group
		self.replicate = int(replicate)
		self.sampleid = group + '_' + str(replicate)
	
	def get_probenames(self):
		'''Return a list of the probe names from the array.'''

		return self.df['ProbeName'].tolist()
	
	def get_genenames(self):
		'''Return a list of the gene names from the array.'''

		return self.df['GeneName'].tolist()
	
	def get_systematicnames(self):
		'''Return a list of the systematics names from the array.'''

		return self.df['SystematicName'].tolist()
	
	def get_descriptions(self):
		'''Return a list of the descriptions from the array.'''

		return self.df['Description'].tolist()
	
	def get_intensities(self):

		'''Return a list of the intensity values from the array.'''
		return self.df['gProcessedSignal'].tolist()

	def normalise(self):
		'''Normalise the intensity values, i.e. log2 transform and 75th percentile shift.'''

		## Threshold the values to 1 if they are less than 1.
		self.df.ix[self.df['gProcessedSignal'] < 1, 'gProcessedSignal'] = 1

		## Log2 transform the values.
		self.df['gProcessedSignal_log2'] = np.log2(self.df['gProcessedSignal'])

		## 75th percentile shift the values.
		self.df['gProcessedSignal_log2_shifted'] = self.df['gProcessedSignal_log2'] - np.percentile(self.df['gProcessedSignal_log2'], 75)

		return self

	def get_normalised_intensities(self):
		'''Return a list of the normalised intensity values from the array.'''

		self.normalise()

		return self.df['gProcessedSignal_log2_shifted'].tolist()

class Experiment:
	'''Holds an experiment, which is a collection of Array instances.'''

	def __init__(self, arrays):
		'''Set up all of the Experiment attributes.

		Munges together all of the intensity values from the arrays into a single data frame.

		Args:
			arrays (list): A list of Array instances.

		'''

		self.arrays = arrays
		master_array = self.arrays[0].df[['FeatureNum', 'ProbeName', 'GeneName', 'SystematicName', 'Description', 'gProcessedSignal_log2_shifted']].copy()
		
		## Rename the column containing the intesity values with the sample ID.
		master_array.rename(columns={'gProcessedSignal_log2_shifted': self.arrays[0].sampleid}, inplace=True)
		
		## Pull out the group names and set the groups attribute.
		self.groups = list(set([array.group for array in self.arrays]))

		## Loop through the arrays and munge them together into a single data frame.
		for i in range(1, len(self.arrays)):
			array = self.arrays[i].df
			array = array[['FeatureNum', 'gProcessedSignal_log2_shifted']].copy()
			array.rename(columns={'gProcessedSignal_log2_shifted': self.arrays[i].sampleid}, inplace=True)
			master_array = pd.merge(master_array, array, on = 'FeatureNum')

		## The dataframe with all of the data in it.
		self.df = master_array

	def baseline_to_median(self):
		'''For each row in the data frame of an Experiment instance, set the median value to zero.'''

		## Get a list of the sample IDs so that just those columns can be pulled out of the data frame.
		norm_exp_cols = self.get_sampleids()

		## Horrid lambda function to correctly set the basline value for each row to the median.
		self.df[norm_exp_cols] = self.df[norm_exp_cols].apply(lambda x : (x - x.median(axis=0)) if not len(x) % 2 == 0 else (x - sorted(x)[len(x)/2]), axis=1)

		return self

	def has_baseline_to_median(self):
		'''Check whether an Experiment instance has had the baseline set to median. Returns a boolean.'''

		norm_exp_cols = self.get_sampleids()
		norm_exp_vals = self.df[norm_exp_cols]

		## Stacks all of the expression values from all of the columns on top of each other
		## in a single column and count the number of instances of each number.
		stacked = norm_exp_vals.stack().value_counts()

		## If the number of zeroes is the same as the number of rows in the data frame (meaning
		## that there is one zero in each row), then the baseline has been set to median.
		if stacked.index[0] == 0.0 and stacked.iloc[0] == norm_exp_vals.shape[0]:
			return True
		else:
			return False

	def get_exp_values(self):
		'''Return just the expression values from an Experiment instance.'''

		norm_exp_cols = self.get_sampleids()

		return self.df[norm_exp_cols]

	def get_groups(self):
		'''Return a list of the groups from an Experiment instance.'''

		return self.groups

	def get_sampleids(self):
		'''Return a list of the sample IDs from an experiment instance.'''

		return list(self.df.columns.values)[5:]

	def remove_sample(self, sampleid):
		'''Remove a sample from an Experiment instance only if it has not had baseline set to median.'''

		if not self.has_baseline_to_median():
			self.df.drop(sampleid, axis=1, inplace=True)
			return self
		elif self.has_baseline_to_median():
			raise TypeError('Experiments with baseline to median cannot have samples removed')
			return

## Functions ##

def read_array(filename, group, replicate):
	'''Read in a raw array data file and return an Array instance.

	Args:
			filename (str): The filename of the raw array data file.
			group (str): The name of the group to which this array belongs.
			replicate (int): The replicate number within the group.

	'''

	f = open(filename, 'r')
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
			if 'GeneName' in headers and 'Description' in headers: 
				values[l[feature_num_ind]] = { 'FeatureNum': l[feature_num_ind], \
											   'ProbeName': l[probe_name_ind], \
											   'GeneName': l[gene_name_ind], \
											   'SystematicName': l[sys_name_ind], \
											   'Description': l[desc_name_ind], \
											   'gProcessedSignal': l[signal_ind] }
			else:
				values[l[feature_num_ind]] = { 'FeatureNum': l[feature_num_ind], \
											   'ProbeName': l[probe_name_ind], \
											   'GeneName': l[sys_name_ind], \
											   'SystematicName': l[sys_name_ind], \
											   'Description': '', \
											   'gProcessedSignal': l[signal_ind] }
	
	## Use the convenient Pandas function to create a data frame from a dictionary.
	array_df = pd.DataFrame.from_dict(values, orient='index')

	return Array(array_df, filename, group, replicate)

def read_experiment(array_filenames, baseline=True):
	'''Read in a series of array data files and return an Experiment instance.

	Args:
			array_filenames (list): List of filenames of the raw array data files.
			baseline (Boolean): Set the baseline to median?

	'''
	
	arrays = []
	
	## Loop through the array files, get the group and replicate, create a
	## normalised Array instance, and then read all of the arrays into an
	## Experiment instance with the baseline set to the median.
	for a in array_filenames:
		group, replicate = a.split('/')[-1].split('.')[0].split('_')
		replicate = int(replicate)
		norm_array = read_array(a, group, replicate).normalise()
		arrays.append(norm_array)

	experiment = Experiment(arrays)
	if baseline is True:
		experiment = experiment.baseline_to_median()
	
	return experiment