# Dots for microarrays

Dots is a Python package for working with microarray data. 
Its back-end is a standalone package for reading in, normalisation, statistical 
analysis and plotting of Agilent single-colour microarray data. Its front-end
isn't finished yet (more on that below).

## Installation

Dots has a number of dependencies including NumPy and SciPy and the least painless
way of getting these is to use the 
[https://store.continuum.io/cshop/anaconda/](Anaconda Python distribution) which includes
NumPy and SciPy and a couple of the other required dependencies like Pandas, Scikit-learn
and Bokeh.

Once you've downloaded and installed Anaconda, you can install the Dots package as follows:

```bash
git clone https://github.com/sandyjmacdonald/dots_for_microarrays
cd dots_for_microarrays
sudo python setup.py install
```

Setuptools should take care of the dependencies but, in testing, I've found the NumPy, SciPy
and Scikit-learn installations to be problematic, hence my recommendation of using Anaconda
to relieve those headaches.

Once you have Anaconda, if you'd like to install and use Dots in a fenced-off virtual 
environment that won't interfere with anything else, then you can do so as follows:

```bash
conda create --yes -n dots_env python=2.7
source activate dots_env

git clone https://github.com/sandyjmacdonald/dots_for_microarrays
cd dots_for_microarrays
sudo python setup.py install
```

To test that the installation has worked, you can run Nosetests as follows:

```bash
sudo python setup.py nosetests
```

## What Dots does

1. Reads in a series of Agilent single-color array files.
**It's important that your array files are named correctly, in order for Dots to work out 
to which group and replicate they belong e.g. for treated and untreated groups each with 
three replicates name the files `treated_1.txt, treated_2.txt, treated_3.txt, untreated_1.txt,
untreated_2.txt, untreated_3.txt`.
2. Normalises the data by log2-transforming, 75th percentile-shifting and setting the baseline
to the median for each gene across all samples.
3. Calculates fold changes and log fold changes for all of the pairs of groups.
4. Runs either a T-test or ANOVA (determined automagically by the number of groups) with 
Benjamini-Hochberg p-value adjustment and a Tukey HSD post hoc test to determine signifcant
pairs from the ANOVA.
5. Provides a number of different visualisations of the data: box and whisker plots of the 
normalised data for each sample, a PCA plot of all of the samples, a hierarchically-clustered
(by gene) heatmap for the significantly differentially expressed genes (> +/- 2-fold, p < 0.05),
a plot of k-means clustered groups of genes with similar expression patterns across the samples,
and volcano plots for each pair of samples.
6. Generates tab-separated tables of the normalised data and fold change and statstical analysis.

## What Dots will do in the future

1. Read the array data into an SQLite3 database, signifcantly speeding the whole workflow if 
you re-analyse your array data at a later date.
2. Assess the quality of the arrays.
3. Provide a web front-end to guide you through the workflow.
4. (Possibly) read in Affymetrix array data.
5. Use linear models (similar to Limma) for the stats.
6. Run functional analyses (Gene Ontology, pathways).
7. ~~Make you a cappuccino while you wait.~~

## Quick start

Dots has a handy workflow script that takes as input a folder containing some Agilent array
files (labelled correctly as explained above) and reads in the data, normalises it, and
produces tables of data and all of the various volcano plots, etc.

You can run it, for example, on the sample data included here (in `dots_sample_data`) with:

```bash
python dots_workflow.py dots_sample_data -o sample_data_output
```

The `-o` is an optional argument and, if you don't include it, then they'll be put in a 
folder named `output`.