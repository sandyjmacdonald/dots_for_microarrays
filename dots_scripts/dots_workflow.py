# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os
import glob
import argparse
from dots_backend.dots_plotting import do_boxplot, do_pcaplot, do_volcanoplot, do_heatmap, do_clusters_plot
from dots_backend.dots_arrays import read_experiment
from dots_backend.dots_analysis import get_fold_changes, write_fcs_stats, write_normalised_expression

## Set up argparse arguments
parser = argparse.ArgumentParser()
parser.add_argument('input', help='input folder, .e.g arrays')
parser.add_argument('-o', '--output', help='name of output folder')
args = parser.parse_args()

## Set up output folder.
if args.output:
	outfolder = args.output if args.output.endswith('/') else args.output + '/'
else:
	outfolder = 'output/'

if not os.path.isdir(outfolder):
	os.makedirs(outfolder)

## Read in files and create experiment.
array_filenames = glob.glob(args.input + '*.txt' if args.input.endswith('/') else args.input + '/*.txt')
experiment = read_experiment(array_filenames)
experiment = experiment.baseline_to_median()

## Write tables.
write_fcs_stats(experiment, outfile=outfolder + 'foldchanges_stats.txt')
write_normalised_expression(experiment, outfile=outfolder + 'normalised_expression.txt')

## Do plots.
do_boxplot(experiment, show=False, image=True, html_file=outfolder + 'boxplot.html')
do_pcaplot(experiment, show=False, image=True, html_file=outfolder + 'pcaplot.html')
do_heatmap(experiment, show=False, image=True, html_file=outfolder + 'heatmap.html')
do_clusters_plot(experiment, show=False, image=True, html_file=outfolder + 'clustersplot.html')

## Get fold change columns for volcano plots.
fcs = get_fold_changes(experiment)
fc_cols = [x for x in fcs.columns.values if 'logFC' in x]

## For each pair of groups, create a volcano plot.
for col in fc_cols:
	pair = col[6:]
	groups = tuple(pair.split('_'))
	do_volcanoplot(experiment, groups, show=False, image=True, html_file=outfolder + pair + '_volcanoplot.html')