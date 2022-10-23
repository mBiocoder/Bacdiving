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



def bacdive_call(bacdive_id = "", bacdive_password = "", input_via_file = 0, input_file_path =" ", search_by_id = False, search_by_culture_collection = False, search_by_taxonomy = False, search_by_16S_seq_accession = False, search_by_genome_accession = False, taxtable_input = 0, taxtable_file_path =" ", print_res_df_ToFile = 1, print_access_stats = 0, output_dir ="./"):
    """
    Reads input, queries BacDive database and stores resulting dataframe and access statistics.

    :param bacdive_id: Log-in credential: BacDive id. Default is empty string. Either bacdive_id and bacdive_pw can be set as function parameters or if left as empty strings, users will be prompted to input credentials.
    :param bacdive_password: Log-in credential: BacDive password. Default is empty string. Either bacdive_id and bacdive_pw can be set as function parameters or if left as empty strings, users will be prompted to input credentials.
    :param input_via_file: If input is a file (.csv, .tsv, or .txt) with one entry per row, then set input_via_file = 1. Default is input_via_file = 0.
    :param input_file_path: If input_via_file = 1, specify the input file path. Default is empty string.
    :param search_by_id: If input_via_file = 1 and the content of the input file are BacDive ids, set search_by_id = True. Default is False.
    :param search_by_culture_collection: If input_via_file = 1 and the content of the input file are culture collection ids, set search_by_culture_collection = True. Default is False.
    :param search_by_taxonomy: If input_via_file = 1 and the content of the input file are species, set search_by_taxonomy = True. Default is False.
    :param search_by_16S_seq_accession: If input_via_file = 1 and the content of the input file are 16S sequence accession ids (e.g. SILVA ids), set search_by_16S_seq_accession = True. Default is False.
    :param search_by_genome_accession: If input_via_file = 1 and the content of the input file are genome accession ids, set search_by_genome_accession = True. Default is False.
    :param taxtable_input: If input is a taxonomy table (e.g. extracted from phyloseq-object), then set taxtable_input = 1. Default is taxtable_input = 0.
    :param taxtable_file_path: If taxtable_input = 1, specify the taxtable file path. Default is empty string.
    :param print_res_df_ToFile: If the resulting dataframe should be printed to file (as BacdiveInformation.tsv), then set print_res_df_ToFile = 1. Default is print_res_df_ToFile = 1.
    :param print_access_stats: If the resulting BacDive access statistics should be printed to file, then set print_access_stats = 1. Default is print_access_stats = 0.
    :param output_dir: Path to where resulting dataframe should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: Resulting dataframe with all strain-level BacDIve information.
    """

    print('If you have not registered for Bacdive web services yet, please do so using the following link before running this package: ' + 'https://sso.dsmz.de/auth/realms/DSMZ/protocol/openid-connect/auth?response_type=code&redirect_uri=https%3A%2F%2Fapi.bacdive.dsmz.de%2Flogin&client_id=api.bacdive&nonce=fc4f465de78388722385d9bd3f82de1f&state=cab5a9d249f9502bd01f7715cc5d99bd&scope=openid')

    try:
        if bacdive_id == "" and bacdive_password == "":
            # Get credentials
            bacdive_id = input("Please enter your BacDive id here: ")
            bacdive_password = input("Please enter your BacDive password here: ")

        if input_via_file == 1 and taxtable_input == 0:
            dfs = []
            client = bacdive.BacdiveClient(bacdive_id, bacdive_password)  # Access Bacdive

            total_number_silva_ids_input = 0
            access_failed = 0
            not_found = []

            # Reading in all SILVA ids and running query in Bacdive
            with open(input_file_path) as file:
                txt_file = csv.reader(file, delimiter="\t")
                with alive_bar(force_tty=True) as bar:
                    for line in txt_file:
                        query_id = "".join(line)
                        total_number_silva_ids_input += 1

                        bar()
                        # Search Bacdive in various ways:
                        if search_by_16S_seq_accession == True and search_by_id == False and search_by_culture_collection == False and search_by_taxonomy == False and search_by_genome_accession == False:
                            query = {"16s": query_id}
                        elif search_by_taxonomy == True and search_by_genome_accession == False and search_by_16S_seq_accession == False and search_by_id == False and search_by_culture_collection == False:
                            query = {"taxonomy": query_id}
                        elif search_by_genome_accession == True and search_by_16S_seq_accession == False and search_by_id == False and search_by_culture_collection == False and search_by_taxonomy == False:
                            query = {"genome": query_id}
                        elif search_by_id == True and search_by_culture_collection == False and search_by_taxonomy == False and search_by_genome_accession == False and search_by_16S_seq_accession == False:
                            query = {"id": query_id}
                        elif search_by_culture_collection == True and search_by_taxonomy == False and search_by_genome_accession == False and search_by_16S_seq_accession == False and search_by_id == False:
                            query = {"culturecolno": query_id}
                        else:
                            print("Please make sure your parameters for your input file are correct.")
                            exit(1)
                        # run query
                        try:
                            with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                                result = client.search(**query)
                                for strain in client.retrieve():
                                    # print(strain)
                                    tmp = pd.json_normalize(strain)
                                    df = pd.DataFrame.from_dict(tmp)
                                    dfs.append(df)  # append the data frame to the list
                        except KeyError:
                            access_failed += 1
                            not_found.append(query_id)

            # Create resulting dataframe
            resulting_df = pd.concat(dfs, ignore_index=True)  # concatenate all the data frames in the list

            # Bacdive access statistics
            if print_access_stats == 1:
                access_stats_for_silva_input(total_number_silva_ids_input, access_failed, not_found,
                                             out_file_path=output_dir + "access_statistics_SILVA_input.txt")

            if print_res_df_ToFile == 1:
                resulting_df_to_file_printer(resulting_df, output_dir=output_dir)

            return resulting_df

        elif taxtable_input == 1 and input_via_file == 0:

            dfs = []
            client = bacdive.BacdiveClient(bacdive_id, bacdive_password)

            access_failed = 0
            not_found = []

            final_df = pd.read_table(taxtable_file_path, index_col=0)

            # Select rows which do not have NaN value in column 'Species'
            selected_rows = final_df[~final_df['Species'].isnull()]
            # Merge values from Genus and Species column for Bacdive query
            df_new = selected_rows.Genus.str.cat(selected_rows.Species, sep=' ')
            # Write resulting species found from input file to an output file
            df_new.to_csv(output_dir + 'Species_names_from_taxtable_file.csv', index=False, header=None)
            # read output file and return this dataframe
            try:
                downstream_df = pd.read_csv(output_dir + "Species_names_from_taxtable_file.csv", header=None)
            except pandas.errors.EmptyDataError:
                print("There are no species information available for this dataset!")
                sys.exit(0)

            total_number_species_input = downstream_df[downstream_df.columns[0]].count()

            with alive_bar(force_tty=True, total=int(total_number_species_input)) as bar:
                for species in downstream_df[0]:
                    bar()
                    # run query
                    try:
                        with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                            result = client.search(taxonomy=species)
                        for strain in client.retrieve():
                            tmp = pd.json_normalize(strain)
                            df = pd.DataFrame.from_dict(tmp)
                            dfs.append(df)  # append the data frame to the list
                    except KeyError:
                        access_failed += 1
                        not_found.append(species)

                resulting_df = pd.concat(dfs, ignore_index=True)  # concatenate all the data frames in the list

                # Bacdive access statistics
                if print_access_stats == 1:
                    access_stats_for_taxtable_input(final_df, selected_rows, total_number_species_input, access_failed,
                                                    not_found,
                                                    out_file_path=output_dir + "access_statistics_taxtable_input.txt")

                if print_res_df_ToFile == 1:
                    resulting_df_to_file_printer(resulting_df, output_dir=output_dir)

                return resulting_df
        else:
            print(
                'If you have not registered for Bacdive web services yet, please do so using the following link before running this package: ' + 'https://sso.dsmz.de/auth/realms/DSMZ/protocol/openid-connect/auth?response_type=code&redirect_uri=https%3A%2F%2Fapi.bacdive.dsmz.de%2Flogin&client_id=api.bacdive&nonce=fc4f465de78388722385d9bd3f82de1f&state=cab5a9d249f9502bd01f7715cc5d99bd&scope=openid')
            print('Only try calling this function after gaining access.')
            exit(1)
    except AttributeError:
        print("Your login credentials were wrong. Please try again!")
        sys.exit(1)


def bacdive_access_for_multiple_inputs(bacdive_id = "", bacdive_password = "", input_lists = {}, output_dir = "./"):
    """
    For multiple input files (either all input files or all taxonomy tables) a single resulting dataframe will be created and stored.

    :param bacdive_id: Log in credential: BacDive id. Default is empty string. If entering the log-in credentials multiple times (for each input sample) is not desirable, then it is recommended to input the credtials once as function parameters.
    :param bacdive_password: Log in credential: BacDive id. Default is empty string. If entering the log-in credentials multiple times (for each input sample) is not desirable, then it is recommended to input the credtials once as function parameters.
    :param input_lists: Dictionary which specifies {<file-path> : <file-type>}. In case of input_via_file as the inpu type, then the exact content types has to be specified. The possible content types are the same as with bacdive_call(). For example: {"input_data/SILVA_ids.txt" : ["input_via_file", "search_by_16S_seq_accession"], "./input_data/taxtable_from_phyloseq/nagel_taxtab.tsv" : ["taxtable_input"]}
    :param output_dir: Path to where resulting dataframe should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: Resulting dataframe with all strain-level BacDive information for all those multiple inputs.
    """
    resulting_list_with_all_res_dfs = []
    for file, filetype in input_lists.items():
        if filetype[0] == "taxtable_input":
            resulting_df = bacdive_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, taxtable_input=1, taxtable_file_path=file, output_dir=output_dir,
                                        print_res_df_ToFile=1, print_access_stats=1)
            resulting_list_with_all_res_dfs.append(resulting_df)
        elif filetype[0] == "input_via_file":
            if filetype[1] == "search_by_id":
                resulting_df = bacdive_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_id= True, output_dir=output_dir,
                                            print_res_df_ToFile=1, print_access_stats=1)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_culture_collection":
                resulting_df = bacdive_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_culture_collection=True,
                                            output_dir=output_dir,
                                            print_res_df_ToFile=1, print_access_stats=1)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_taxonomy":
                resulting_df = bacdive_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_taxonomy=True,
                                            output_dir=output_dir,
                                            print_res_df_ToFile=1, print_access_stats=1)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_16S_seq_accession":
                resulting_df = bacdive_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_16S_seq_accession=True,
                                            output_dir=output_dir,
                                            print_res_df_ToFile=1, print_access_stats=1)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_genome_accession":
                resulting_df = bacdive_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_genome_accession=True,
                                            output_dir=output_dir,
                                            print_res_df_ToFile=1, print_access_stats=1)
                resulting_list_with_all_res_dfs.append(resulting_df)

        else:
            print("Please check if the filetype and search parameters are set correctly!")
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

