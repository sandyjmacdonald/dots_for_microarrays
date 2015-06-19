# dots for microarrays

dots guides you through uploading your raw microarray data, normalising and 
analysing it, and displays your data in a series of interactive plots.

## Simple analysis of Agilent one-color arrays

dots only supports Agilent one-color arrays at the moment. The workflow is 
similar to the Agilent GeneSpring software in terms of functionality. We plan 
to add support for analysis of Affymetrix arrays in the future.

### Upload

Your raw microarray data are uploaded and archived, and the normalised data are 
stored in a SQLite 3 database. Each experiment is stored in its own database 
and you can easily download all of the normalised data.

### Normalisation

Once your data are uploaded, you can check the QC values and abnormal values are 
flagged up. Data are normalised by log2-transforming and 75th-percentile shifting 
per-sample, then setting the baseline between arrays to the median per-gene.

### Analysis

Pick a group to compare against or compare all-versus-all. A T-test or one-way 
ANOVA is used to identify significantly differentially expressed genes. Visualise 
your data in tables and a range of different plots.

## dots_backend quickstart

The backend of dots can be used to analyse your Agilent array data right in Python,
giving you the flexibility to integrate the functionality of dots with other 
workflows.

Here's how simple it is to read in a series of samples and create a box and whisker
plot showing the distributions of the intensity values in each sample.

```python
# -*- coding: utf-8 -*-
#!/usr/bin/env python

import glob
from dots_plotting import do_boxplot
from dots_arrays import read_experiment

array_filenames = glob.glob('../sample_data/*.txt')
experiment = read_experiment(array_filenames)
experiment = experiment.baseline_to_median()

do_boxplot(experiment)
```

That will create an HTML file called boxplot.html with your box and whisker
plot in it that can be opened in your browser and then saved as an image.

It's equally simple to read that experiment we just created into a SQLite3 
database.

```python
from dots_db import load_experiment_into_db

db_id = load_experiment_into_db(experiment)
```

This reads all of the data into a standalone SQLite3 database for this 
experiment and returns the unqiue identifier (a sixteen digit hex ID) for
the newly created database.

Here's an example of how to query the database and extract the first 10
rows from it.

```python
import sqlite3 as lite

conn = lite.connect(db_id)

def query_db():
	with conn:    
		cur = conn.cursor()    
		cur.execute('SELECT Probes.ProbeID, Genes.GeneName, Intensities.NormIntensity, Samples.GroupName \
					 FROM Probes, Genes, Intensities, Samples \
					 WHERE Probes.GeneUID = Genes.GeneUID \
					 AND Probes.ProbeID = Intensities.ProbeID \
					 AND Intensities.SampleID = Samples.SampleID')
		rows = cur.fetchall()
		return rows[:10]

rows = query_db()
```

More detailed documentation for the dots backend will be added soon.