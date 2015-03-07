# -*- coding: utf-8 -*-
#!/usr/bin/env python

import bokeh.plotting as bp
import bokeh.models as bm
import pandas as pd
from dots_analysis import do_pca
from dots_arrays import Experiment

## Functions ##

def get_n_colours(n):
	'''Create a palette of n contrasting colours.

	Args:
		n (int): The number of colours required.

	Returns:
		A list of hex codes for the n colours required.

	'''
	kelly_colours = [
    '#FFB300', # Vivid Yellow
    '#803E75', # Strong Purple
    '#FF6800', # Vivid Orange
    '#A6BDD7', # Very Light Blue
    '#C10020', # Vivid Red
    '#CEA262', # Grayish Yellow
    '#817066', # Medium Gray

    # The following don't work well for people with defective color vision
    '#007D34', # Vivid Green
    '#F6768E', # Strong Purplish Pink
    '#00538A', # Strong Blue
    '#FF7A5C', # Strong Yellowish Pink
    '#53377A', # Strong Violet
    '#FF8E00', # Vivid Orange Yellow
    '#B32851', # Strong Purplish Red
    '#F4C800', # Vivid Greenish Yellow
    '#7F180D', # Strong Reddish Brown
    '#93AA00', # Vivid Yellowish Green
    '#593315', # Deep Yellowish Brown
    '#F13A13', # Vivid Reddish Orange
    '#232C16', # Dark Olive Green
    ]

	return kelly_colours[:n]

def create_standard_plot(h=600, w=600, title='', x_range=None, tools='previewsave'):
	'''Create a standard plot and set consistent theme. Saves specifiying these every time.

	Args:
		h (int): The height (in pixels) of the plot.
		w (int): The width (in pixels) of the plot.
		title (str): The title of the plot.
		x_range (list): A list of variables for the x axis.
		tools (str): The tools to be displayed at the top of the chart.

	Returns:
		A plot object.

	'''
	if x_range is not None:
		plot = bp.figure(tools=tools, background_fill='#E5E5E5', x_range=x_range, title=title, plot_height=h, plot_width=w)
	else:
		plot = bp.figure(tools=tools, background_fill='#E5E5E5', title=title, plot_height=h, plot_width=w)

	plot.xgrid.grid_line_color = 'white'
	plot.ygrid.grid_line_color = 'white'
	plot.grid.grid_line_width = 1
	plot.xaxis.major_label_text_font_size='10pt'
	plot.yaxis.major_label_text_font_size='10pt'
	plot.xaxis.axis_label_text_font_size='12pt'
	plot.yaxis.axis_label_text_font_size='12pt'

	return plot

def do_boxplot(experiment):
	'''Create a box plot from an experiment instance, one box for each sample.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	'''
	## Checks whether it is an experiment instance or just the expression values.
	if isinstance(experiment, pd.DataFrame):
		groups = experiment.unstack().groupby(level=0)
	elif isinstance(experiment, Experiment):
		groups = experiment.get_exp_values().unstack().groupby(level=0)

	## Gets the quantiles for the boxes.
	q1 = groups.quantile(q=0.25)
	q2 = groups.quantile(q=0.5)
	q3 = groups.quantile(q=0.75)

	## Gets the interquartile range.
	iqr = q3 - q1

	## Gets the upper and lower limits for the whiskers.
	upper = q2 + 1.5*iqr
	lower = q2 - 1.5*iqr

	## Lists to store the values in.
	outx = []
	outy = []

	group_names = []

	## Add all of the data to the lists.
	for group_name, group in groups:
		group_names.append(group_name)
		for value in group[(group > upper[group_name]) | (group < lower[group_name])].get_values():
			outx.append(group_name)
			outy.append(value)

	## Sort the groups into alphabetical order.
	group_names.sort()

	## Create the plot.
	bp.output_file('boxplot.html')
	boxplot = create_standard_plot(h=600, w=900, x_range=group_names)
	boxplot.xgrid.grid_line_color = None
	
	qmin = groups.quantile(q=0.00)
	qmax = groups.quantile(q=1.00)

	## List comps to get the values into the correct format.
	upper_values = [min([x,y]) for (x,y) in zip(qmax.values.tolist(), upper.values.tolist())]
	lower_values = [max([x,y]) for (x,y) in zip(qmin.values.tolist(), lower.values.tolist())]

	## Bokeh doesn't have a built-in box blot type, so we have to hack it together.
	boxplot.segment(group_names, upper_values, group_names, q3.values, line_width=1, line_color='black')
	boxplot.segment(group_names, lower_values, group_names, q1.values, line_width=1, line_color='black')

	boxplot.rect(group_names, (q3.values+q2.values)/2, 0.4, q3.values-q2.values, fill_color='#9FCCFF', line_width=1, line_color='#265B99')
	boxplot.rect(group_names, (q2.values+q1.values)/2, 0.4, q2.values-q1.values, fill_color='#9FCCFF', line_width=1, line_color='#265B99')

	boxplot.rect(group_names, lower_values, 0.2, 0.01, line_color='#265B99')
	boxplot.rect(group_names, upper_values, 0.2, 0.01, line_color='#265B99')

	boxplot.rect(outx, outy, 0.05, 0.01, color='#FF4444', fill_alpha=0.6)

	## Shows the plot.
	bp.show(boxplot)

def do_pcaplot(experiment):
	'''Create a simple scatter plot of PCA scores for the samples.

	Args:
		experiment (Experiment instance): An instance of the Experiment class.

	'''
	## Run the PCA and get the results back in a data frame, set the colours of each value.
	pca_df = do_pca(experiment)
	groups = list(set(pca_df['group']))
	colourmap = dict(zip(groups, get_n_colours(len(groups))))
	pca_df['colour'] = pca_df['group'].map(lambda x: colourmap[x])
	source = bp.ColumnDataSource(data=dict(x=pca_df['xvals'].tolist(), y=pca_df['yvals'].tolist(), group=pca_df['group'].tolist(), sampleid=pca_df['sampleid'].tolist()))

	## Create the plot, label axes and format the points.
	bp.output_file('pcaplot.html')
	pcaplot = create_standard_plot(tools='hover,previewsave')
	pcaplot.xaxis.axis_label = 'Principal Component 1'
	pcaplot.yaxis.axis_label = 'Principal Component 2'
	pcaplot.circle(pca_df['xvals'], pca_df['yvals'], color=pca_df['colour'], fill_alpha=0.6, size=10, source=source)

	## Set up the hover tooltips.
	hover = pcaplot.select(dict(type=bm.HoverTool))
	hover.tooltips = [('x', '@x'), ('y', '@y'), ('Group', '@group'), ('Sample', '@sampleid')]

	## Shows the plot.
	bp.show(pcaplot)