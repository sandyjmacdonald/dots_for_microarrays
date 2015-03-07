# -*- coding: utf-8 -*-
#!/usr/bin/env python

import sys
import uuid
import sqlite3 as lite

## Funtions ##

def make_new_experiment_db():
	'''Create a new SQLite3 database with tables and return its ID.'''
	conn = None
	db_id = None

	## Try to connect and create all of the tables.
	try:
		## Creates a random 16 digit hexadecimal ID. Should always be random.
		db_id = 'db/' + str(uuid.uuid4().hex.upper()[:16]) + '.db'

		## Conect to the database and get a cursor.
		conn = lite.connect(db_id)
		cur = conn.cursor()

		## Create all of the tables.
		cur.execute('CREATE TABLE Samples(SampleUID INTEGER PRIMARY KEY, SampleID TEXT, GroupName TEXT)')
		cur.execute('CREATE TABLE Genes(GeneUID INTEGER PRIMARY KEY, GeneName TEXT, SysName TEXT, GeneDescription TEXT)')
		cur.execute('CREATE TABLE Probes(ProbeID TEXT PRIMARY KEY, GeneUID INTEGER)')
		cur.execute('CREATE TABLE Intensities(IntensityUID INTEGER PRIMARY KEY, ProbeID TEXT, SampleID TEXT, NormIntensity REAL)')
		
		## Commit the tables.
		conn.commit()

	## Catch any errors and roll back database.
	except lite.Error, e:
		if conn:
			conn.rollback()

		print 'Error: %s' % e.args[0]
		sys.exit(1)

	## If there's a connection, close it.
	finally:
		if conn:
			conn.close()

	return db_id

def load_into_db(db_id, table_name, table_data):
	'''Given a database ID, table name and data, read data into the specified database table.

	Args:
		db_id (str): The database identifier, returned by the make_new_experiment_db() function.
		table_name (str): The name of a table in the database into which data will be read.
		table_data (list): A nested list of data to read into the table. Dealt with by the load_experiment_into_db() function.

	Returns:
		True if successful.

	'''
	conn = None
	num_cols = len(table_data[0])

	## Try to connect to the database and get a cursor.
	try:
		## Connect and get a cursor.
		conn = lite.connect(db_id)
		cur = conn.cursor()

		## Build and execute the query.
		query = 'INSERT INTO %s VALUES(%s)' % (table_name, ', '.join(['?'] * num_cols))
		cur.executemany(query, table_data)

		## Commit to the database.
		conn.commit()

	## Catch any errors and roll back database.
	except lite.Error, e:
		if conn:
			conn.rollback()

		print 'Error: %s' % e.args[0]
		sys.exit(1)

	## If there's a connection, close it.
	finally:
		if conn:
			conn.close()

	return True

def load_experiment_into_db(experiment):
	'''Complete function to create an experiment database and read experiment data in.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	Returns:
		Database ID if successful.

	'''
	db_id = make_new_experiment_db()
	exp_vals = experiment.get_exp_values()
	
	## Get the required data.
	ints_table_data = []
	probe_ids = experiment.df['ProbeName'].tolist()
	sample_ids = exp_vals.columns.values
	ints_table_cols = 4
	ints_uid = 1
	
	## Loop through each sample and read into the list of data.
	for sample in sample_ids:
		sample_vals = exp_vals[sample].tolist()
		ints_table_data.extend(zip(range(ints_uid, ints_uid + len(sample_vals)), probe_ids, [sample]*len(sample_vals), [float(v) for v in sample_vals]))
		ints_uid += len(sample_vals)

	## Load into the database table and clean up after.
	load_into_db(db_id, 'Intensities', ints_table_data)
	del ints_table_data

	## Get the required data.
	gene_names = experiment.df['GeneName'].tolist()
	sys_names = experiment.df['SystematicName'].tolist()
	descs = experiment.df['Description'].tolist()
	genes_tuples = zip(probe_ids, gene_names, sys_names, descs)
	genes_tuples = list(set(genes_tuples))
	gene_uids = range(1, len(genes_tuples) + 1)
	probes_table_data = zip([g[0] for g in genes_tuples], gene_uids)

	## Load into the database table and clean up after.
	load_into_db(db_id, 'Probes', probes_table_data)
	del probes_table_data

	## Get the required data.
	genes_tuples = [g[1:] for g in genes_tuples]
	genes_lists = zip(*genes_tuples)
	genes_table_data = zip(gene_uids, *genes_lists)

	## Load into the database table and clean up after.
	load_into_db(db_id, 'Genes', genes_table_data)
	del genes_table_data

	## Get the required data.
	sample_uids = range(1, len(sample_ids) + 1)
	groups = [sample_id.split('_')[0] for sample_id in sample_ids]
	samples_table_data = zip(sample_uids, sample_ids, groups)

	## Load into the database table and clean up after.
	load_into_db(db_id, 'Samples', samples_table_data)
	del samples_table_data

	return db_id