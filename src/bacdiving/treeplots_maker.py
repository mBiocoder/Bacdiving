#!/usr/bin/python3

from .circulartree import ct_plot

import toytree
import toyplot
import toyplot.pdf
import toyplot.svg

import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from ete3 import Tree
from anytree import Node
from anytree.exporter import DotExporter

import matplotlib
matplotlib.style.use('ggplot')
from matplotlib import cm
import logging
from pylab import *
plt.style.use('ggplot')
import pandas as pd


def traverse(df_, a, i, innerl):
    # Traverse tree
    if i+1 < df_.shape[1]:
        a_inner = pd.unique(df_.loc[np.where(df_.iloc[:, i] == a)].iloc[:, i+1])

        desc = []
        for b in a_inner:
            desc.append(traverse(df_, b, i+1, innerl))
        if innerl:
            il = a
        else:
            il = ""
        out = f"({','.join(desc)}){il}"
    else:
        out = a

    return out


def df2newick(df, levels, inner_label=True):
    # Transform dataframe to newick string
    df_tax = df.loc[:, [x for x in levels if x in df.columns]]

    alevel = pd.unique(df_tax.iloc[:, 0])
    strs = []
    for a in alevel:
        strs.append(traverse(df_tax, a, 0, inner_label))

    newick = f"({','.join(strs)});"
    return newick.replace(" ", "_")


def newick_to_linkage(newick: str, label_order: [str] = None) -> (np.ndarray, [str]):
    # Transform newick string to linkage
    tree = Tree(newick, format=8)
    cophenetic_matrix, newick_labels = tree.cophenetic_matrix()
    cophenetic_matrix = pd.DataFrame(cophenetic_matrix, columns=newick_labels, index=newick_labels)

    if label_order is not None:
        # sanity checks
        missing_labels = set(label_order).difference(set(newick_labels))
        superfluous_labels = set(newick_labels).difference(set(label_order))
        assert len(missing_labels) == 0, f'Some labels are not in the newick string: {missing_labels}'
        if len(superfluous_labels) > 0:
            logging.warning(f'Newick string contains unused labels: {superfluous_labels}')

        # reorder the cophenetic_matrix
        cophenetic_matrix = cophenetic_matrix.reindex(index=label_order, columns=label_order)

    # reduce square distance matrix to condensed distance matrices
    pairwise_distances = pdist(cophenetic_matrix)

    # return linkage matrix and labels
    return linkage(pairwise_distances), list(cophenetic_matrix.columns)


def overview_treeplot(resulting_df, pallete = "brg", colormap1 = "bwr", column_name1 = "Culture and growth conditions.culture temp.temperature", column_name2 = "Physiology and metabolism.oxygen tolerance.oxygen tolerance", label_name1 = "Category1", label_name2 = "Category2", colormap2 = "Wistia", fontsize = 14, figsize = [20,10], saveToFile = True, output_dir = "./"):
    """
    Makes overview tree plot showing hierarchical tree structure for all species of input as well as max. 2 BacDive columns of interst.

    :param resulting_df: Resulting dataframe as outputted by bacdive_call().
    :param pallete: Color palette used. Default is "brg".
    :param colormap1: Color map used for first column of interest. Default is "bwr".
    :param column_name1: First column of interest from resulting_df to plot. Default is "Culture and growth conditions.culture temp.temperature".
    :param column_name2: Second column of interest from resulting_df to plot. Default is "Physiology and metabolism.oxygen tolerance.oxygen tolerance".
    :param label_name1: Legend label for first column of interest. Default is "Category1".
    :param label_name2: Legend label for second column of interest. Default is "Category2".
    :param colormap2: Color map for second column of interest. Default is "Wistia".
    :param fontsize: Size of font. Default is 16.
    :param figsize: Size of plot. Default is [20,10].
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: Overview plot
    """

    # convert dataframe to Newick format
    myNewick = df2newick(resulting_df, levels=['Name and taxonomic classification.LPSN.domain',
                                               'Name and taxonomic classification.LPSN.phylum',
                                               'Name and taxonomic classification.LPSN.class',
                                               'Name and taxonomic classification.LPSN.order',
                                               'Name and taxonomic classification.LPSN.family',
                                               'Name and taxonomic classification.LPSN.genus',
                                               'Name and taxonomic classification.LPSN.species'])

    linkage_matrix, labels = newick_to_linkage(newick=myNewick)

    Z2 = sch.dendrogram(linkage_matrix, labels=labels, no_plot=True, leaf_rotation=90)
    labels_Z2_labels = Z2["ivl"]

    cat1_array = []
    cat2_array = []

    for i in labels_Z2_labels:
        if resulting_df.loc[resulting_df["Name and taxonomic classification.species"] == i.replace("_", " ")][column_name1].isnull().values.all() == True:
            cat1_array.append([1., 0., 0., 1.])
        else:
            cat1_array.append([0., 0., 1., 1.])

    for j in labels_Z2_labels:
        if resulting_df.loc[resulting_df["Name and taxonomic classification.species"] == j.replace("_", " ")][column_name2].isnull().values.all() == True:
            cat2_array.append([0.98823529, 0.49803922, 0. , 1. ])
        else:
            cat2_array.append([0.89411765, 1. , 0.47843137, 1. ])

    type_num = 2
    _cmp = cm.get_cmap(colormap1, type_num)  # Setting 2 different colors
    _cmp2 = cm.get_cmap(colormap2, type_num)  # Setting another 2 different colors

    colors_dict = {label_name1: cat1_array,
                   label_name2: cat2_array}

    # Specify the legend of the color labels.
    colors_legends = {label_name1: {"colors": [[0., 0., 1., 1.], [1., 0., 0., 1.]],
                                     "labels": ["yes", "no"]},
                      label_name2: {"colors": [[0.89411765, 1. , 0.47843137, 1. ],[0.98823529, 0.49803922, 0. , 1. ]],
                                     "labels": ["yes", "no"]}}

    if saveToFile == False:
        ct_plot(Z2, pallete=pallete, fontsize=fontsize, figsize=figsize, colorlabels=colors_dict, colorlabels_legend=colors_legends, show=True)
    else:
        ct_plot(Z2, pallete=pallete, fontsize=fontsize, figsize=figsize, colorlabels=colors_dict,
             colorlabels_legend=colors_legends, show=False, output_dir = output_dir)



def circular_treeplot(resulting_df, width = 1400, height = 1400, saveToFile = True, output_format ="pdf", output_dir ="./"):
    """
    Makes tree plot showing hierarchical tree structure for all species of input.

    :param resulting_df:  Resulting dataframe as outputted by bacdive_call().
    :param width: Width of tree plot. Default value is 1400.
    :param height: Height of tree plot. Default value is 1400.
    :param saveToFile: Boolean to save plot or not. Default is True.
    :param output_format: Output file type. Possible file formats include: pdf, svg and html. Default is "pdf".
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: Circular treeplot.
    """

    # convert dataframe to Newick format
    myNewick = df2newick(resulting_df, levels=['Name and taxonomic classification.LPSN.domain',
                                               'Name and taxonomic classification.LPSN.phylum',
                                               'Name and taxonomic classification.LPSN.class',
                                               'Name and taxonomic classification.LPSN.order',
                                               'Name and taxonomic classification.LPSN.family',
                                               'Name and taxonomic classification.LPSN.genus',
                                               'Name and taxonomic classification.LPSN.species'])
    tre1 = toytree.tree(myNewick, tree_format=1)

    if saveToFile == False:  #Attention: if saveToFile == False and the program is run via console the resulting plot may not be visible. Howver, it will work if run on notebook
        tre1.draw(
            layout='c',
            edge_type='c',
            width=width,
            height=height,
        );
    elif saveToFile == True:
        canvas, axes, mark = tre1.draw(
            node_hover=True,
            node_labels=False,
            layout='c',
            edge_type='c',
            node_sizes=[8 if i else 0 for i in tre1.get_node_values()],
            node_style={"stroke": "black"},
            width=width,
            height=height );

        if output_format == "pdf":
            toyplot.pdf.render(canvas, output_dir +"tree-plot.pdf")
        elif output_format == "svg":
            toyplot.svg.render(canvas, output_dir +"tree-plot.svg")
        elif output_format == "html":
            toyplot.html.render(canvas, output_dir +"tree-plot.html")
        else:
            print("This output format is not supported! Please choose from .pdf, .svg or .html.")

