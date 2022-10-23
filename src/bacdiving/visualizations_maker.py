#!/usr/bin/python3

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



def pieplot_maker(resulting_df, plot_column, title = " ", ylabel_name = " ", saveToFile = False, output_dir = "./", figsize = [6.4, 4.8]):
    """
    Makes pie plot for any categorical column of interest (e.g. Morphology.cell morphology.motility) from resulting dataframe

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
    Makes bar plot for any continuous column of interest (e.g. Morphology.cell morphology.cell length) from resulting dataframe

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
    Makes world map displaying all countries where species from the input originate from
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



def symbol_remover(inputList, symbol):
    # Given a list as input, the output is the same list but now the symbol of interest is removed for each list element
    result = []
    for part in inputList:
        split_string = part.split(symbol, 1)
        substring = split_string[0]
        result.append(substring)
    return result


def float_seperator(x):
    # Given a list x as input with float value ranges, the output is the same list but now all value ranges are replaced with the mean value instead
    result = []
    for part in x:
        if '-' in part:
            a, b = part.split('-')
            a, b = float(a), float(b)
            mean = (a + b) / 2
            result.append(mean)
        else:
            a = float(part)
            result.append(a)
    return result


def convertToInt(x):
    #:param x: List with value ranges like '30-40'
    #:return: List with value ranges transformed into seperate int values
    result = []
    for part in x:
        if '-' in part:
            a, b = part.split('-')
            a, b = int(float(a)), int(float(b))
            result.extend(range(a, b + 1))
        else:
            a = int(float(part))
            result.append(a)
    return result


