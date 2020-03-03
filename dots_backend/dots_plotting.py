#!/usr/bin/env python3

import os
import subprocess
import random
import io
import math
import pandas as pd
import numpy as np
import bokeh.plotting as bp
import bokeh.models as bm

from bokeh.layouts import gridplot
from bokeh.io import export_png

from PIL import Image
from collections import OrderedDict
from palettable.colorbrewer.diverging import RdYlGn_11

from .dots_analysis import run_pca, run_stats, get_fold_changes, get_clusters
from .dots_arrays import Experiment


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
        plot = bp.figure(tools=tools, background_fill_color='#E5E5E5', x_range=x_range, title=title, plot_height=h, plot_width=w)
    else:
        plot = bp.figure(tools=tools, background_fill_color='#E5E5E5', title=title, plot_height=h, plot_width=w)

    plot.xgrid.grid_line_color = 'white'
    plot.ygrid.grid_line_color = 'white'
    plot.grid.grid_line_width = 1
    plot.xaxis.major_label_text_font_size='10pt'
    plot.yaxis.major_label_text_font_size='10pt'
    plot.xaxis.axis_label_text_font_size='12pt'
    plot.yaxis.axis_label_text_font_size='12pt'

    return plot

def chunks(l, n):
    '''Yield successive n-sized chunks from l.

    Args:
        l (list): The list to chunk.
        n (int): Number of items in each chunk.

    '''
    for i in range(0, len(l), n):
        yield l[i:i+n]

def do_boxplot(experiment, show=False, image=False, html_file='boxplot.html'):
    '''Create a box plot from an experiment instance, one box for each sample.

    Args:
        experiment (Experiment instance): An instance of the Experiment class.
        show (Boolean): Should the plot be shown.
        image (Boolean): Should the plot be saved as an image.
        html_file (string): The name of the html output file.

    '''

    #  Checks whether it is an experiment instance or just the expression values.
    if isinstance(experiment, pd.DataFrame):
        groups = experiment.unstack().groupby(level=0)
    elif isinstance(experiment, Experiment):
        groups = experiment.get_exp_values().unstack().groupby(level=0)

    #  Gets the quantiles for the boxes.
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)

    #  Gets the interquartile range.
    iqr = q3 - q1

    #  Gets the upper and lower limits for the whiskers.
    upper = q2 + 1.5*iqr
    lower = q2 - 1.5*iqr

    #  Lists to store the outliers in.
    outx = []
    outy = []

    group_names = []

    #  Sort the groups into alphabetical order.
    if all([group.split('_')[0].isdigit() for group in group_names]):
            group_names = sorted(group_names, key=lambda x:int(x.split('_')[0]))
    else:
        group_names = sorted(group_names)

    #  Add them to the list of group names.
    for group_name, group in groups:
        group_names.append(group_name)

    #  Create the plot.
    bp.output_file(html_file)
    boxplot = create_standard_plot(h=600, w=900, x_range=group_names)
    boxplot.xgrid.grid_line_color = None
    boxplot.xaxis.major_label_orientation = np.pi/3

    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)

    #  List comps to get the values into the correct format.
    upper_values = [min([x,y]) for (x,y) in zip(qmax.values.tolist(), upper.values.tolist())]
    lower_values = [max([x,y]) for (x,y) in zip(qmin.values.tolist(), lower.values.tolist())]

    #  Bokeh doesn't have a built-in box blot type, so we have to hack it together.
    boxplot.segment(group_names, upper_values, group_names, q3.values, line_width=1, line_color='black')
    boxplot.segment(group_names, lower_values, group_names, q1.values, line_width=1, line_color='black')

    boxplot.rect(group_names, (q3.values+q2.values)/2, 0.4, q3.values-q2.values, fill_color='#9FCCFF', line_width=1, line_color='#265B99')
    boxplot.rect(group_names, (q2.values+q1.values)/2, 0.4, q2.values-q1.values, fill_color='#9FCCFF', line_width=1, line_color='#265B99')

    boxplot.rect(group_names, lower_values, 0.2, 0.01, line_color='#265B99')
    boxplot.rect(group_names, upper_values, 0.2, 0.01, line_color='#265B99')

    #  Hackiness to make sure the plots don't fail with many values.
    #  Randomly samples 250,000 values, which seems to work.

    df_copy = experiment.df.copy()
    exp_shape = experiment.get_exp_values().shape
    num_rows = 250000/exp_shape[1]

    if num_rows < exp_shape[0]:
        experiment.df = experiment.df.ix[random.sample(experiment.df.index, num_rows)]
        groups = experiment.get_exp_values().unstack().groupby(level=0)

    #  Add all of the outliers to the lists.
    for group_name, group in groups:
        for value in group[(group > upper[group_name]) | (group < lower[group_name])]:
            outx.append(group_name)
            outy.append(value)

    boxplot.rect(outx, outy, 0.05, 0.01, color='#FF4444', fill_alpha=0.6)

    #  Shows the plot.
    if show == True:
        bp.show(boxplot)
    else:
        bp.save(obj=boxplot)

    if image == True:
        filename = html_file.split('.')[0] + '.png'
        export_png(boxplot, filename=filename, height=600, width=600, webdriver=None, timeout=5)

    experiment.df = df_copy

    return html_file

def do_pcaplot(experiment, show=False, image=False, html_file='pcaplot.html'):
    '''Create a simple scatter plot of PCA scores for the samples.

    Args:
        experiment (Experiment instance): An instance of the Experiment class.
        show (Boolean): Should the plot be shown.
        image (Boolean): Should the plot be saved as an image.
        html_file (string): The name of the html output file.

    '''

    #  Run the PCA and get the results back in a data frame, set the colours of each value.
    pca_df = run_pca(experiment)
    groups = list(set(pca_df['group']))
    colourmap = dict(list(zip(groups, get_n_colours(len(groups)))))
    pca_df['colour'] = pca_df['group'].map(lambda x: colourmap[x])
    source = bp.ColumnDataSource(data=dict(x=pca_df['xvals'].tolist(), y=pca_df['yvals'].tolist(), colour=pca_df['colour'].tolist(), group=pca_df['group'].tolist(), sampleid=pca_df['sampleid'].tolist()))

    #  Create the plot, label axes and format the points.
    bp.output_file(html_file)
    pcaplot = create_standard_plot(tools='hover,previewsave')
    pcaplot.xaxis.axis_label = 'Principal Component 1'
    pcaplot.yaxis.axis_label = 'Principal Component 2'
    pcaplot.circle(x='x', y='y', color='colour', fill_alpha=0.6, size=10, source=source)

    #  Set up the hover tooltips.
    hover = pcaplot.select(dict(type=bm.HoverTool))
    hover.tooltips = [('x', '@x'), ('y', '@y'), ('Group', '@group'), ('Sample', '@sampleid')]

    #  Shows the plot.
    if show == True:
        bp.show(pcaplot)
    else:
        bp.save(obj=pcaplot)

    if image == True:
        filename = html_file.split('.')[0] + '.png'
        export_png(pcaplot, filename=filename, height=600, width=600, webdriver=None, timeout=5)

    return html_file

def do_volcanoplot(experiment, groups, show=False, image=False, html_file='volcano_plot.html'):
    '''Create a volcano plot for a specified pair of samples.

    Args:
        experiment (Experiment instance): An instance of the Experiment class.
        groups (tuple): The pair of samples to compare, e.g. ('treated', 'untreated').
        show (Boolean): Should the plot be shown.
        image (Boolean): Should the plot be saved as an image.
        html_file (string): The name of the html output file.

    '''

    #  Run the fold changes and stats functions on our experiment.
    fcs = get_fold_changes(experiment)
    stats = run_stats(experiment)

    #  Strip out extraneous data from stats data frame.
    sig_cols = [x for x in stats.columns.values if 'significant' in x]
    stats = stats[['FeatureNum', 'p_val', 'p_val_adj'] + sig_cols].copy()

    #  Merge the fold changes and stats data frames, create new columns for
    #  -log 10 adj. p value and spot colours.
    merged_df = pd.merge(fcs, stats, on = 'FeatureNum')
    fc_cols = [x for x in merged_df.columns.values if 'logFC' in x]
    fc_tuples = [(col[6:].split('_')[0], col[6:].split('_')[1]) for col in fc_cols]
    sorted_groups = tuple(sorted(groups))
    sig_col = 'significant_%s_%s' % sorted_groups
    rev_groups = groups[::-1]
 
    if groups in fc_tuples:
        column = 'logFC_%s_%s' % groups
    elif rev_groups in fc_tuples:
        orig_column = 'logFC_%s_%s' % rev_groups
        column = 'logFC_%s_%s' % groups
        merged_df[column] = merged_df[orig_column] * -1
 
    merged_df['neg_log_10_p_val'] = -1 * np.log10(merged_df['p_val_adj'])
 
    if len(fc_cols) > 1:
        merged_df['colour'] = np.where((abs(merged_df[column]) > 1) & (merged_df['p_val_adj'] < 0.05) & (merged_df[sig_col] == True), 'blue', 'red')
    else:
        merged_df['colour'] = np.where((abs(merged_df[column]) > 1) & (merged_df['p_val_adj'] < 0.05), 'blue', 'red')
 
    source = bp.ColumnDataSource(data=dict(logfc=merged_df[column].tolist(), neglog10pval=merged_df['neg_log_10_p_val'].tolist(), pval=merged_df['p_val_adj'].tolist(), colour=merged_df['colour'].tolist(), gene=merged_df['GeneName'].tolist()))
 
    #  Create the plot, label axes and format the points.
    bp.output_file(html_file)
    volcanoplot = create_standard_plot(tools='hover,pan,wheel_zoom,box_zoom,reset,save')
    volcanoplot.xaxis.axis_label = 'log 2 fold change, ' + column
    volcanoplot.yaxis.axis_label = '- log 10 p value'
    volcanoplot.circle(x='logfc', y='neglog10pval', color='colour', fill_alpha=0.4, line_alpha=0, size=5, source=source)

    #  Set up the hover tooltips.
    hover = volcanoplot.select(dict(type=bm.HoverTool))
    hover.tooltips = [('log 2 FC', '@logfc'), ('adj. p value', '@pval'), ('gene', '@gene')]

    #  Shows the plot.
    if show == True:
        bp.show(volcanoplot)
    else:
        bp.save(obj=volcanoplot)

    if image == True:
        filename = html_file.split('.')[0] + '.png'
        export_png(volcanoplot, filename=filename, height=600, width=600, webdriver=None, timeout=5)

    return html_file

def do_heatmap(experiment, show=False, image=False, html_file='heatmap.html'):
    '''Create a heatmap with normalised expression data with significant
    differential expression.

    Args:
        experiment (Experiment instance): An instance of the Experiment class.
        show (Boolean): Should the plot be shown.
        image (Boolean): Should the plot be saved as an image.
        html_file (string): The name of the html output file.

    '''

    #  Creates a new data frame with nornalised expression data and cluster
    #  numbers from hierarchical clustering.
    cluster_df = get_clusters(experiment, how='hierarchical')

    #  Make a separate data frame with just the normalised expression values.
    norm_exp_cols = experiment.get_sampleids()
    vals_only = cluster_df[norm_exp_cols]

    #  Set up the colour scheme.
    colourmap = RdYlGn_11.hex_colors[::-1]
    limit = np.abs(vals_only.values).max()
    val_to_colour = lambda x: int((x + limit) * (len(colourmap)/(2*limit))) - 1

    #  Empty lists to which we can add the data for the heatmap.
    vals = []
    colour = []
    colname = []
    rowname = []
    cluster = []
    gene = []

    #  Column names and feature IDs for setting the x and y ranges.
    columns = list(vals_only.columns.values)
    feats = cluster_df['FeatureNum'].tolist()

    #  Loop through the columns and values in each column and add them
    #  to the lists for the heatmap.
    for col in columns:
        for i in range(len(vals_only[col])):
            vals.append(vals_only[col][i])
            colour.append(colourmap[val_to_colour(vals_only[col][i])])
            colname.append(col)
            rowname.append(cluster_df['FeatureNum'].iloc[i])
            gene.append(cluster_df['GeneName'].iloc[i])
            cluster.append(cluster_df['cluster'].iloc[i])

    #  For the hover tooltips.
    source = bp.ColumnDataSource(data=dict(colname=colname, rowname=rowname, colour=colour, gene=gene, cluster=cluster, vals=vals))

    height = len(cluster_df) * 2

    #  Create the plot.
    bp.output_file(html_file)
    heatmap = bp.figure(x_range=columns, y_range=list(reversed([str(f) for f in feats])), x_axis_location='above', plot_width=600, plot_height=height, title=None, tools='hover,previewsave')
    heatmap.rect('colname', 'rowname', 1.0, 1.0, source=source, color="colour", line_color=None)
    heatmap.grid.grid_line_color = None
    heatmap.axis.axis_line_color = None
    heatmap.axis.major_tick_line_color = None
    heatmap.xaxis.major_label_text_font_size='10pt'
    heatmap.yaxis.major_label_text_font_size='5pt'
    heatmap.yaxis.major_label_text_color = None
    heatmap.xaxis.axis_label_text_font_size='10pt'
    heatmap.yaxis.axis_label_text_font_size='10pt'
    heatmap.axis.major_label_standoff = 0
    heatmap.xaxis.major_label_orientation = np.pi/3

    #  Set up the hover tooltips.
    hover = heatmap.select(dict(type=bm.HoverTool))
    hover.tooltips = OrderedDict([('Group', '@colname'), ('Gene', '@gene'), ('Cluster', '@cluster'), ('Mean norm. exp.', '@vals')])

    #  Shows the plot.
    if show == True:
        bp.show(heatmap)
    else:
        bp.save(obj=heatmap)

    if image == True:
        filename = html_file.split('.')[0] + '.png'
        export_png(heatmap, filename=filename, height=height, width=600, webdriver=None, timeout=5)

    return html_file

def do_clusters_plot(experiment, show=True, image=False, html_file='clustersplot.html'):
    '''Create a series of line plots from k-means-clustered normalised expression
    data'

    Args:
        experiment (Experiment instance): An instance of the Experiment class.
        show (Boolean): Should the plot be shown.
        image (Boolean): Should the plot be saved as an image.
        html_file (string): The name of the html output file.

    '''

    #  Creates a new data frame with nornalised expression data and cluster
    #  numbers from k-means clustering.
    cluster_df = get_clusters(experiment, how='kmeans')
    samples = experiment.get_sampleids()
    clusters = cluster_df['cluster'].unique()

    #  Set up the plot and create a list to add the plots to.
    bp.output_file(html_file)
    plots = []

    #  For each cluster from the k-means, create a line plot.
    for c in clusters:
        cluster = cluster_df[cluster_df['cluster'] == c]
        vals = cluster[samples]
        xvals = vals.columns.values.tolist()
        yvals = vals.values.tolist()
        clusterplot = create_standard_plot(tools='save', x_range=xvals)
        clusterplot.xaxis.major_label_orientation = np.pi/3
        clusterplot.yaxis.axis_label = 'normalised expression'

        for y in yvals:
            clusterplot.line(xvals, y)

        plots.append(clusterplot)

    #  Make a grid from the plots.
    grid = gridplot(plots, ncols=int(math.sqrt(len(plots))))

    #  Shows the plot.
    if show == True:
        bp.show(grid)
    else:
        bp.save(obj=grid)

    height = 600 * len(clusters)

    if image == True:
        filename = html_file.split('.')[0] + '.png'
        export_png(grid, filename=filename, webdriver=None, timeout=5)

    return html_file