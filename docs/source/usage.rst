Usage
=====

About Bacdiving
---------------

Bacdiving is a Python package which can access and retrieve information from the world’s largest database for standardized bacterial phenotypic information: `BacDive <https://bacdive.dsmz.de/>`_.
Additionally, Bacdiving provides access statistics and options to visualize this information.

The following figure gives an overview of Bacdiving and how this package could be used:

.. figure:: Workflow.png
    :width: 700px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center


As depicted in this workflow, Bacdiving can deal with two types of inputs: either a taxonomy
table (e.g. resulting from a phyloseq-object) or a file input.

Let us begin with the file input type (in p x 1 format). The
file can contain p rows of either all BacDive-IDs, culture collection numbers, taxonomy
information (either as full name or as list with genus name, species with optional epithet,
and optional subspecies) or sequence accession numbers (either 16S sequences, SILVA-IDs
or genome sequences). For each one of the p rows BacDive is then queried and all strain-level information is stored in a single dataframe.
This resulting dataframe is of the format p x r with r being the number
of BacDive columns for which we have information. 
This dataframe, along with the corresponding access statistics, can then be written to file. All other core functions in BacDiving rely on the resulting dataframe.

The second possible input type is a taxonomy table (in p x 7 format) which has the following 7 taxonomic ranks: kingdom, phylum, class, order, family, genus, species.
Bacdiving first filters out all rows for which the species is unknown. This results in a
”new” taxonomy table (in m x 7 format). Each one of these m species will then be checked
if it can be found on BacDive or not. If information is available for a given species, then BacDive data for all of its known strains will be
appended into a single dataframe. In the end, this resulting dataframe (m x r) will contain all strain-level information for all species of the taxonomy table. Thus, the resulting dataframe and the 
corresponding BacDive access statistics results can be written to file. Additionally, given a taxonomy table as input, one may wish to see Bacdive information across all strains of a given species in order to perform 
various other downstream statistical tasks (e.g. matrix completion, regression, etc.). Therefore, given specific BacDive features of interest, a flattened file (in p x r' format) can be outputted which contains the number
of strains per species found on BacDive as well as the majority values across all strains for given BacDive features of interest.

Depending on the research question, either BacDiving’s visualization options can be used
or custom visualizations can be made using the resulting dataframe. There is really no limit on
how you can extend and make use of the resulting dataframe.

For instance, this resulting dataframe along with metadata (which is often stored in
phyloseq-objects) could be used in tools like `NetCoMi <https://github.com/stefpeschel/NetCoMi>`_ to construct various types of networks. 
The nodes of these networks could be colored with specific phylogenetic information from BacDive as stored in the resulting dataframe file which in turn may explain why
a given network looks the way it does. In other words, coloring the nodes in a network
based on phylogentic information may explain the underlying correlation between various features and conditions in a dataset.

Accessing BacDive
-----------------

As soon as you have registered on BacDive, you can use your credentials to run Bacdiving’s most central function :py:func:`bacdiving.bacdive_call()`:

.. autofunction:: bacdiving.bacdive_call()

The first thing :py:func:`bacdiving.bacdive_call()` does is, it will prompt you to input your login credentials prior to querying BacDive, if you did not input your credentials via the function parameters ``"bacdive_id"`` and ``"bacdive_password"``. 

After that, it generates the resulting dataframe(s) (BacdiveInformation.tsv) with all strain-level information and 
it can output the BacDive access statistics (if the parameter is set) as a .txt-file which gives information on the percentage of input
species found on BacDive and also lists all species which could not be found on BacDive. Additional files (like Species_names_from_taxtable_file.csv or Flattened_Bacdive_data.tsv) may as well be outputted if your input was a taxonomy table. Note that the file Species_names_from_taxtable_file.csv lists all species from the taxonomy table, even prior to querying BacDive.

For accessing specific data entries in your resulting dataframe you can either run :py:func:`bacdiving.get_resulting_df_values()` or :py:func:`bacdiving.access_list_df_objects()`. 

.. autofunction:: bacdiving.get_resulting_df_values()

.. autofunction:: bacdiving.access_list_df_objects()

However, :py:func:`bacdiving.access_list_df_objects()` is only designed to be used if you are interested in retrieving information for either pH, temperature or halophily (e.g. prior to making a box plot), whereas :py:func:`bacdiving.get_resulting_df_values` is more generic.

Visualizations
--------------

Bacdiving supports 8 different visualization types:

1. Circular hierarchical taxonomic tree plot (also referred to as overview tree plot since it gives information on which species have what kind of BacDive information):

.. autofunction:: bacdiving.overview_treeplot()

A similar circular hierarchical tree plot but without BacDive information can be created as well:

.. autofunction:: bacdiving.circular_treeplot()

2. Stacked bar plot to show relative abundance (of e.g. different genera) per sample:

.. autofunction:: bacdiving.stacked_barplot_relative_abundance()

3. Pie chart to plot information like oxygen tolerance:

.. autofunction:: bacdiving.pieplot_maker()

4. World map to show all countries (not water bodies!) of origin for a given set of species:

.. autofunction:: bacdiving.worldmap_maker()

5. Fatty acid profile plot for a fatty acid of interest:

.. autofunction:: bacdiving.fatty_acid_profile()

6. Frequency plot (of e.g. most frequent sampling type):

.. autofunction:: bacdiving.freqplot_maker()

7. Box plot to compare e.g. optimal temperature ranges for various species

.. autofunction:: bacdiving.boxplot_maker()

8. Bar plot to compare e.g. cell length of different species

.. autofunction:: bacdiving.barplot_maker()

