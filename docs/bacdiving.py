#!/usr/bin/python3
import contextlib
import csv
import os

import bacdive
import pandas as pd
import pandas.errors
from pylab import *

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from alive_progress import alive_bar

import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib
from matplotlib.lines import Line2D

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


import matplotlib
matplotlib.style.use('ggplot')
from pylab import *
plt.style.use('ggplot')
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

import worldmap as wm
from bokeh.palettes import Spectral, Turbo256, Bokeh, Category20
import itertools


def bacdive_call(bacdive_id ="", bacdive_password ="", inputs_list = [""], sample_names = [""], print_res_df_ToFile = True, print_access_stats = True, print_flattened_file = False,  columns_of_interest = [""], output_dir ="./"):
    """
    For multiple input files (either all input files or all taxonomy tables) this function reads the input, queries the BacDive database and stores resulting dataframe(s) and access statistics.

    :param bacdive_id: Log in credential: BacDive id. Default is empty string. If entering the log-in credentials multiple times (for each input sample) is not desirable, then it is recommended to input the credentials once as function parameters.
    :param bacdive_password: Log in credential: BacDive password. Default is empty string. If entering the log-in credentials multiple times (for each input sample) is not desirable, then it is recommended to input the credentials once as function parameters.
    :param inputs_list: List which specifies (multiple) strings. Each string has the structure: "<file-path> <file-type> (<content-type>)" and is thus seperated by space(s). Content-type is, however, only required if you have input_via_file; it can have one of the following values: "search_by_id", "search_by_culture_collection", "search_by_taxonomy", "search_by_16S_seq_accession" or "search_by_genome_accession".
    :param sample_names: List of samples names.
    :param print_res_df_ToFile: Print the resulting dataframe with all Bacdive information to file or not. Default is True.
    :param print_access_stats: Print the Bacdive access statistics to file or not. Default is True.
    :param print_flattened_file: Print the flattened Bacdive information for certain columns of interest to file or not. Default is False.
    :param columns_of_interest: Specify in this list which columns from BacdiveInformation.tsv you want to include in the flattened file. Default is empty list of strings.
    :param output_dir: Path to where resulting dataframe should be saved. Default is current directory ("./").
    :return: List containing the resulting dataframe(s)  with all strain-level BacDive information for all those multiple inputs.
    """

    print('If you have not registered for Bacdive web services yet, please do so using the following link before running this package: ' + 'https://sso.dsmz.de/auth/realms/DSMZ/protocol/openid-connect/auth?response_type=code&redirect_uri=https%3A%2F%2Fapi.bacdive.dsmz.de%2Flogin&client_id=api.bacdive&nonce=fc4f465de78388722385d9bd3f82de1f&state=cab5a9d249f9502bd01f7715cc5d99bd&scope=openid')

    resulting_list_with_all_res_dfs = []
    for element in inputs_list:
        sample_name = sample_names[inputs_list.index(element)]
        element_list = element.split(" ")
        file = element_list[0]
        filetype_first = element_list[1]

        if len(element_list) == 3:
            filetype_second = element_list[2]

        if filetype_first == "taxtable_input":
            resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, taxtable_input=1, taxtable_file_path=file, sample_name=sample_name, output_dir=output_dir,
                                     print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
            resulting_list_with_all_res_dfs.append(resulting_df)
            if print_flattened_file == True:
                flattened_full_file_maker(resulting_df, columns_of_interest= columns_of_interest, sample_name=sample_name, taxtable_file_path=file, output_dir=output_dir)

        elif filetype_first == "input_via_file":
            if filetype_second == "search_by_id":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_id= True, sample_name=sample_name, output_dir=output_dir,
                                         print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype_second == "search_by_culture_collection":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_culture_collection=True, sample_name=sample_name,
                                         output_dir=output_dir, print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype_second == "search_by_taxonomy":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_taxonomy=True,
                                         output_dir=output_dir, sample_name=sample_name,
                                         print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype_second == "search_by_16S_seq_accession":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_16S_seq_accession=True,
                                         output_dir=output_dir, sample_name=sample_name,
                                         print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype_second == "search_by_genome_accession":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_genome_accession=True,
                                         output_dir=output_dir, sample_name=sample_name, print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)
        else:
            print("Please check if the filetype, sample names and search parameters are set correctly!")
            exit(1)

    return resulting_list_with_all_res_dfs

############################# treeplots_maker ################################

def overview_treeplot(resulting_df, pallete = "brg", colormap1 = "bwr", column_name1 = "Culture and growth conditions.culture temp.temperature", column_name2 = "Physiology and metabolism.oxygen tolerance.oxygen tolerance", label_name1 = "Category1", label_name2 = "Category2", colormap2 = "Wistia", fontsize = 14, figsize = [20,10], saveToFile = True, output_dir = "./"):
    """
    Makes overview tree plot showing hierarchical tree structure for all species of input as well as maximum 2 BacDive columns of interest.

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



def pieplot_maker(resulting_df, plot_column, title = " ", ylabel_name = " ", saveToFile = False, output_dir = "./", figsize = [6.4, 4.8]):
    """
    Makes pie plot for any categorical column of interest from resulting dataframe

    :param resulting_df:  Resulting dataframe as outputted by bacdive_call().
    :param plot_column: (Categorical) Column of interest from resulting_df. Default is empty string.
    :param title:  Title for this plot. Default is empty string.
    :param ylabel_name: y-axis label name. Default is emtpy string.
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :param figsize: Size of the resulting plot. Default is [6.4, 4.8].
    :return: Pie plot
    """
    try:
        number_Nan = sum(pd.isnull(resulting_df[plot_column]))
        title_string = "excluding " + str(number_Nan) + " Nan-values"

        fig, ax = plt.subplots(figsize=figsize)
        resulting_df.groupby(plot_column).size().plot(kind='pie', autopct='%1.0f%%')
        ax.set_ylabel(ylabel_name, size=15)
        ax.set_title(title)
        figtext(.5, .85, title_string, fontsize=10, ha='center')
        if saveToFile == False:
            plt.show()
        elif saveToFile == True:
            fig.savefig(output_dir + '%s piechart.pdf' % title, bbox_inches="tight")
    except KeyError:
        print("No information available for " + plot_column)


def barplot_maker(resulting_df, plot_column = " ", title = " ", ylabel_name = " ", xlabel_name = " ", color= "green", species_list = [], saveToFile=False, output_dir = "./", figsize = [15,10]):
    """
    Makes bar plot for any continuous column of interest from resulting dataframe

    :param resulting_df: Resulting dataframe as outputted by bacdive_call().
    :param plot_column: (Categorical) Column of interest from resulting_df. Default is empty string.
    :param title: Title for this plot. Default is empty string.
    :param ylabel_name: y-axis label name. Default is emtpy string.
    :param xlabel_name: x-axis label name. Default is emtpy string.
    :param color: Color of bars. Default is "green".
    :param species_list: List of species of interest. Default is empty list.
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :param figsize:  Size of the resulting plot. Default is [15, 10].
    :return: Bar plot
    """
    resulting_df['contains_species'] = resulting_df["Name and taxonomic classification.species"].isin(species_list)
    resulting_df = resulting_df.loc[resulting_df['contains_species'] == True]
    resulting_df2 = resulting_df[resulting_df[plot_column].notna()]

    number_Nan = sum(pd.isnull(resulting_df[plot_column]))
    title_string = "excluding " + str(number_Nan) + " Nan-values"

    fig, ax = plt.subplots(figsize = figsize)
    data_values = float_seperator(resulting_df2[plot_column].str.replace(r'µm', '').dropna())

    bars = resulting_df2["Name and taxonomic classification.species"]
    y_pos = np.arange(len(data_values))
    plt.bar(y_pos, data_values, color=color)
    ax.set_ylabel(ylabel_name, size=14)
    ax.set_xlabel(xlabel_name, size=14)
    ax.set_title(title + "\n" + title_string)
    plt.xticks(y_pos, bars, rotation=90)
    if saveToFile == False:
        plt.show()
    elif saveToFile == True:
        fig.savefig(output_dir + '%s barplot.pdf' % ylabel_name, bbox_inches="tight")

def freqplot_maker(resulting_df, plot_column = " ", title = " ", ylabel_name = " ", saveToFile=False, output_dir = "./", figsize = [15, 10]):
    """
    Makes freq plot for any categorical column of interest from resulting dataframe

    :param resulting_df: Resulting dataframe as outputted by bacdive_call().
    :param plot_column: (Categorical) Column of interest from resulting_df. Default is empty string.
    :param title: Title for this plot. Default is empty string.
    :param ylabel_name: y-axis label name. Default is emtpy string.
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :param figsize: Size of the resulting plot. Default is [15, 10].
    :return: Frequency plot
    """
    number_Nan = sum(pd.isnull(resulting_df[plot_column]))
    title_string = "excluding " + str(number_Nan) + " Nan-values"

    fig, ax = plt.subplots(figsize= figsize)
    resulting_df[plot_column].dropna().value_counts().plot(ax=ax, kind='barh')
    ax.set_title(title+ "\n" + title_string)
    ax.set_ylabel(ylabel_name, size=12)
    ax.set_xlabel("frequency", size=12)
    if saveToFile == False:
        plt.show()
    if saveToFile == True:
        fig.savefig(output_dir + '%s freqplot.pdf' % ylabel_name, bbox_inches="tight")


def access_list_df_objects(resulting_df, plot_column = " ", plot_category = " ", temp = 0, pH = 0, halophily = 0, species_list = []):
    """
    Access all categories of interest only for the pH, temperature and halophily columns from the resulting dataframe

    :param resulting_df:  Resulting dataframe as outputted by bacdive_call().
    :param plot_column: Column of interest from resulting_df; should be either pH, temperature or halophily. Default is empty string.
    :param plot_category: Category of interest from column of interest from resulting_df. Default is empty string.
    :param temp: Either one of temp, pH or halophily can be accessed. If temp = 1, temp will be accessed. Default is temp = 0.
    :param pH: Either one of temp, pH or halophily can be accessed. If pH = 1, pH will be accessed. Default is pH = 0.
    :param halophily: Either one of temp, pH or halophily can be accessed. If halophily = 1, halophily will be accessed. Default is halophily = 0.
    :param species_list: List of species of interest. Default is empty list.
    :return: Dictionary: <species> : <values>
    """
    try:
        getValues = lambda key, inputData: [subVal[key] for subVal in inputData if key in subVal]

        resulting_df['contains_species'] = resulting_df["Name and taxonomic classification.species"].isin(species_list)
        resulting_df = resulting_df.loc[resulting_df['contains_species'] == True]
        number_Nan = sum(pd.isnull(resulting_df[plot_column]))
        print("excluding " + str(number_Nan) + " Nan-values" + " from " + plot_column)
        resulting_df2 = resulting_df[resulting_df[plot_column].notna()]
        resulting_dict = {}

        for index, row in resulting_df2.iterrows():
            if temp == 1 and pH == 0 and halophily == 0:
                temperature_list = getValues(plot_category, resulting_df2[plot_column][index])
                temp_data = convertToInt(temperature_list)
                species_name = resulting_df2["Name and taxonomic classification.species"][index]
                resulting_dict.update({species_name: temp_data})
            elif pH == 1 and temp == 0 and halophily == 0:
                pH_list = getValues(plot_category, resulting_df2[plot_column][index])
                pH_data = convertToInt(pH_list)
                species_name = resulting_df2["Name and taxonomic classification.species"][index]
                resulting_dict.update({species_name: pH_data})
            elif halophily == 1 and temp == 0 and pH == 0:
                halophily_list = getValues(plot_category, resulting_df2[plot_column][index])
                res = symbol_remover(halophily_list, " %")
                halophily_data = convertToInt(res)
                species_name = resulting_df2["Name and taxonomic classification.species"][index]
                resulting_dict.update({species_name: halophily_data})
            else:
                print("Please set the parameters correctly.")
                exit(1)
        return resulting_dict
    except KeyError:
        print("No information available for " + plot_column)


def fatty_acid_profile(resulting_df, species = "", title = "Fatty acid profile plot", figsize= [10, 10], barwidth = 0.05, fontsize = 6, saveToFile = True, output_dir = "./"):
    """
    Makes fatty acid profile plot for any one fatty acid of interest of interest from resulting dataframe

    :param resulting_df: Resulting dataframe as outputted by bacdive_call().
    :param species: Species of interest (e.g. Bacteroides vulgatus). Default is empty string.
    :param title: Title for this plot. Default is "Fatty acid profile plot".
    :param figsize: Size of the resulting plot. Default is [10, 10].
    :param barwidth: Width of the bars. Default is 0.05 cm.
    :param fontsize: Size of the font. Default is 6.
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: Fatty acid profile plot
    """
    resulting_dict = {}

    resulting_df = resulting_df.loc[resulting_df["Name and taxonomic classification.species"] == species]
    list1 = resulting_df["Physiology and metabolism.fatty acid profile.fatty acids"].dropna().to_numpy()

    for i in list1:
        for key in i:
            fatty_acid_name = key["fatty acid"]
            percentage = key["percentage"]
            ECL_value = key["ECL"]
            resulting_dict[fatty_acid_name] = [percentage, ECL_value]

    fig, ax = plt.subplots(figsize= figsize)
    labels = []
    x_axis = []
    y_axis = []
    for k in list(resulting_dict.items()):
        labels.append(k[0])
        x_axis.append(k[1][0])
        y_axis.append(k[1][1])

    ax.bar(x_axis, y_axis, width = barwidth)
    ax.set_xlabel("ECL (equivalent chain length)")
    ax.set_ylabel("Percentage of fatty acids")
    ax.set_title(title)

    rects = ax.patches
    for rect, label in zip(rects, labels):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2, height + 0.4, label, ha="center", va="bottom", fontdict = None, size= fontsize)

    if saveToFile == False:
        plt.show()
    if saveToFile == True:
        fig.savefig(output_dir + '%s_fatty_acid_profile.pdf' % species, bbox_inches="tight")


def get_resulting_df_values(resulting_df, plot_column = " ", plot_category = " ", species_list = []):
    """
    Access all categories of interest only for a column of interest from the resulting dataframe

    :param resulting_df: Resulting dataframe as outputted by bacdive_call().
    :param plot_column: Column of interest from resulting_df. Default is empty string.
    :param plot_category: Category of interest from column of interest from resulting_df. Default is empty string.
    :param species_list: List of species of interest. Default is empty list.
    :return: Dictionary: <species> : <values>
    """
    getValues = lambda key, inputData: [subVal[key] for subVal in inputData if key in subVal]
    resulting_df['contains_species'] = resulting_df["Name and taxonomic classification.species"].isin(species_list)
    resulting_df = resulting_df.loc[resulting_df['contains_species'] == True]
    resulting_df2 = resulting_df[resulting_df[plot_column].notna()] #remove nan rows for column of interest

    number_Nan = sum(pd.isnull(resulting_df[plot_column]))
    title_string = "excluding " + str(number_Nan) + " Nan-values"
    print(title_string)

    resulting_dict = {}

    for index, row in resulting_df2.iterrows():
        list = getValues(plot_category, resulting_df2[plot_column][index])
        species_name = resulting_df2["Name and taxonomic classification.species"][index]
        resulting_dict.update({species_name: list})
    return resulting_dict

def boxplot_maker(resulting_dict, title = " ", xlabel_name = " ", ylabel_name = " ", saveToFile = False, output_dir = "./", figsize = [15, 10]):
    """
    Makes box plot given a dictionary with values of interest.

    :param resulting_dict: Dictionary input with values (e.g. temperature or pH).
    :param title: Title for this plot. Default is emtpy string.
    :param xlabel_name: x-axis label name. Default is emtpy string.
    :param ylabel_name: y-axis label name. Default is emtpy string.
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :param figsize: Size of the resulting plot. Default is [15, 10].
    :return: Box plot
    """
    fig, ax = plt.subplots(figsize = figsize)
    fig.set_tight_layout(True)
    try:
        labels, data = [*zip(*resulting_dict.items())]  # 'transpose' items to parallel key, value lists

        plt.boxplot(data)
        ax.set_title(title)
        ax.set_xlabel(xlabel_name)
        ax.set_ylabel(ylabel_name)
        plt.xticks(range(1, len(labels) + 1), labels, rotation=90)
        if saveToFile == False:
            plt.show()
        elif saveToFile == True:
            fig.savefig(output_dir + 'All_species %s.pdf' % ylabel_name, bbox_inches="tight")
    except ValueError:
        print("No values are found to plot for " + title)
    except AttributeError:
        print("No information available for " + title)


def stacked_barplot_relative_abundance(resulting_list_with_all_res_dfs, top_x = 15, sample_names = [], plot_column ="", title =" ", title_label = "", saveToFile = True, output_dir ="./", figsize = [15, 10]):
    """
    Makes stacked bar plot for any taxonomy level from resulting dataframe.

    :param resulting_list_with_all_res_dfs: Resulting dataframe as outputted by bacdive_access_for_multiple_inputs().
    :param top_x: Limit for how many different color categories should be seen in the plot. Default is 15, meaning every other catagory will be fall under the category "other".
    :param sample_names: List of names for each sample in the resulting_list_with_all_res_dfs. Default is empty list.
    :param plot_column: Taxonomy level of interest (e.g. Name and taxonomic classification.genus). Default is empty string.
    :param title: Title for this plot. Default is emtpy string.
    :param title_label: Title for legend (e.g. Genus). Default is empty string.
    :param saveToFile: Boolean to save plot as a .pdf file or not. Default is True.
    :param output_dir: Path to where plot should be saved if saveToFile is set to True. Default is current directory ("./").
    :param figsize: Size of the resulting plot. Default is [15, 10].
    :return: Stacked bar plot
    """
    resulting_table = {}
    counter = 0
    for i in resulting_list_with_all_res_dfs:
        counter = counter + 1
        value_for_dict = i[plot_column].dropna().value_counts(normalize=True).to_dict()
        resulting_table[counter] = value_for_dict
    resulting_table = pd.DataFrame.from_dict(resulting_table).transpose().fillna(0)
    resulting_table["samples"] = sample_names
    a = resulting_table[resulting_table.columns[:-1]].values
    a.sort(axis=1)  # no ascending argument
    a = a[:, ::-1]  # reverse
    resulting_table = pd.DataFrame(a, resulting_table[resulting_table.columns[:-1]].index, resulting_table[resulting_table.columns[:-1]].columns)
    resulting_table = resulting_table.iloc[:, :top_x]
    resulting_table["other"] = 1 - resulting_table.sum(axis=1)
    resulting_table["samples"] = sample_names
    all_genera = resulting_table[resulting_table.columns[:-1]]
    all_genera_count = len(all_genera.columns)

    color_iter = itertools.cycle(Category20[20])
    colors = [next(color_iter) for genus in all_genera]

    ax = resulting_table.plot(x='samples', kind='bar', stacked=True, title=title, figsize= figsize, color = colors)
    ax.set_ylabel('Relative Abundance in %')
    plt.legend(title=title_label, bbox_to_anchor=(1.0, 1), loc='upper left')
    if saveToFile == False:
        plt.show()
    elif saveToFile == True:
        fig = ax.get_figure()
        fig.savefig(output_dir + 'Relative_abundance_stacked_barplot.pdf', bbox_inches="tight")


def worldmap_maker(resulting_df):
    """
    Makes world map displaying all countries where species from the input originate from.

    :param resulting_df: Resulting dataframe as outputted by bacdive_call()
    :return: World map plot
    """
    number_Nan = sum(pd.isnull(resulting_df["Isolation, sampling and environmental information.isolation.country"]))
    title_string = "excluding " + str(number_Nan) + " Nan-values"
    print(title_string)
    country_list = resulting_df["Isolation, sampling and environmental information.isolation.country"].dropna()
    county_names = country_list.replace("Pacific Ocean", "").replace("Atlantic Ocean", "").replace("Indian Ocean","").replace("Arctic Ocean", "")
    #print(county_names)
    wm.plot(county_names, map_name='world', cmap='Set1')

