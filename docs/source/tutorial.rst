Tutorial
========

We start by importing Bacdiving:

.. code-block:: python

	from bacdiving import bacdive_caller as bc
   	from bacdiving import treeplots_maker as tm
   	from bacdiving import visualizations_maker as vm

Now assume we have the following taxonomy table taxtab.tsv (which we prior extracted from a phyloseq-object):

.. csv-table:: taxtab.tsv
   :header: " ", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
   :widths: 15, 20, 20, 20, 20, 20, 20, 20

   "ASV1", "Bacteria",	"Bacteroidota",	"Bacteroidia",	"Bacteroidales", "Bacteroidaceae", "Bacteroides", "vulgatus"
   "ASV2", "Bacteria", 	"Firmicutes",	"Clostridia",	"Lachnospirales", "Lachnospiraceae",	"Blautia", "NA"
   "ASV3", "Bacteria",	"Bacteroidota",	"Bacteroidia",	"Bacteroidales", "Bacteroidaceae", "Bacteroides", "NA"
   "ASV4", "Bacteria",	"Bacteroidota",	"Bacteroidia",	"Bacteroidales", "Bacteroidaceae", "Bacteroides", "uniformis"
   "ASV5", "Bacteria",	"Firmicutes",	"Clostridia",	"Oscillospirales", "Ruminococcaceae",	"Faecalibacterium", "NA"

.. note:: This taxonomy table stems from Nagel et al. (2016) but you can use any taxonomy table in this .tsv format to follow along with this tutorial. For demonstration purposes this table only shows the first 5 rows of the taxonomy table. 

.. warning:: Note how species-level information is simply not always known for all ASVs in a taxonomy table. On average, you can expect 80% - 95% of all species from a given taxonomy table to be documented in BacDive. The exact percentage is mentioned in the access statistics output file which is generated after running :py:func:`bacdiving.bacdive_call()`.

To get the resulting dataframe with all strain-level BacDive information for all the species in this taxonomy table you can run:

.. code-block:: python

	#Taxonomy table input (e.g. as extracted from phyloseq-object)
	resulting_df = bc.bacdive_call(taxtable_input = 1, taxtable_file_path='./input_data/taxtable_from_phyloseq/taxtab.tsv', output_dir="./", print_res_df_ToFile= 1, print_access_stats=1)

Assuming you do not have a taxonomy table, but have a simple file as input instead which looks something like this:

.. csv-table:: SILVA_ids.txt
   :widths: 15

   "AB681649"
   "AB121974"
   "EU847536"
   "D30768"
   "L35516"
   "AB681086"
   "AB052706"
   "AF144407"
   "AF363064"
   "AJ430586"  

.. note:: For demonstration purposes this SILVA_ids.txt file only contains 10 SILVA-ids. Note that other input types for querying BacDive are possible as well, e.g. taxonomy (as in a list of all species names of interest), Bacdive-ids, culture collection ids or genome accession ids. However, it is important that a single file can not contain multiple input types. 

Given your input file, you can run the following, depending on your input type:

.. code-block:: python

	#SILVA id query
	resulting_df = bc.bacdive_call(input_via_file=1, input_file_path="./input_data/SILVA_ids.txt",search_by_16S_seq_accession=True, print_res_df_ToFile=1, print_access_stats=1)

	#Taxonomy query
	resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="./input_data/Taxonomy_ids.txt", search_by_taxonomy = True, output_dir="./", print_res_df_ToFile= 1, print_access_stats=1)

	#Bacdive id query
	resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="./input_data/Bacdive_ids.txt", search_by_id = True, output_dir="./", print_res_df_ToFile= 1)

	#Culture collection query
	resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="Culture_col_ids.txt", search_by_culture_collection = True, taxtable_input = 0, output_dir="./", print_res_df_ToFile= 1)

	#Genome accession query
	resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="Genome_accession_ids.txt", search_by_genome_accession = True, output_dir="./", print_res_df_ToFile= 1)


If you have multiple inputs of possibly different input types, you can run the following command:

.. code-block:: python

	# Run for multiple inputs (of possibly different input types)
	resulting_list_with_all_res_dfs = bc.bacdive_access_for_multiple_inputs(input_lists={"./input_data/SILVA_ids.txt" : ["input_via_file", "search_by_16S_seq_accession"], "./input_data/taxtable_from_phyloseq/taxtab.tsv" : ["taxtable_input"]})


Now that you have the resulting dataframe at hand, you are almost ready to start visualizing the data. 

.. warning:: If you try to plot information for a column which is not present in the resulting dataframe or if your parameters are set incorrectly or do not match the resulting dataframe, you may get an error. So, make sure to get to know your resulting dataframe (and especially its columns) beforehand.

Let's first take a look at the resulting dataframe:

.. code-block:: python

	#Dataframe exploration
   	print(resulting_df.head()) # prints head of resulting dataframe
  	print(len(resulting_df.index)) #print number of resulting_df rows
   	print(resulting_df.keys())  #print resulting_df column names
   	print(resulting_df.iloc[0:5, 2:4]) #print all specific column information via column index
   	print(resulting_df["Physiology and metabolism.oxygen tolerance.oxygen tolerance"].unique()) #print unique values in a given column
   	print(resulting_df.loc[resulting_df["Name and taxonomic classification.species" ] == "Bacteroides uniformis"])  # print all strains and all columns for Bacteroides uniformis from resulting_df
   	print(resulting_df.loc[resulting_df["Name and taxonomic classification.species"] == "Helicobacter pylori"]["Culture and growth conditions.culture temp"])  # print all strains for column "Culture and growth conditions.culture temp" for Helicobacter pylori from resulting_df
   	print(resulting_df.loc[(resulting_df["Name and taxonomic classification.species"] == "Helicobacter pylori") & (resulting_df["Isolation, sampling and environmental information.isolation.country"] == "Germany")]) #Subset resulting_df to certain parameters
   	print(len(resulting_df.loc[resulting_df["Name and taxonomic classification.species"] == 'Zhihengliuella alba'].index)) #find out how many strains are present for a given species
   	print(get_resulting_df_values(resulting_df, "Culture and growth conditions.culture pH", "pH", species_list=["Helicobacter pylori", "Bacteroides clarus", "Actinomyces odontolyticus", "Bacteroides salyersiae", "Zhihengliuella alba"])) #Given a list of species of interest, access elements in resulting_df which are nested

Great, now we know the basic structure of our resulting dataframe and what kind of BacDive information we have, so it is finally time to start plotting! 

.. note:: There are many possibilites for which columns from the resulting dataframe can be plotted for each plotting function. This tutorial shall only demonstrate a few examples.


In order to first gain some overview of our data, let us start with Bacdiving's overview treeplot. Assume we want to know for which species BacDive information on temperature and oxygen tolerance is known or not. We can do this by running the following command:

.. code-block:: python

	#Overview treeplot
	tm.overview_treeplot(resulting_df, label_name1="Temperature", label_name2="Oxygen tolerance", saveToFile=True, output_dir="./")

This results in the following plot:

.. figure:: nagel_overview_treeplot.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center

If you do not want the BacDive information to be shown and just prefer the hierarchical taxonomy tree plot, then run:

.. code-block:: python

	#Circular treeplot
	tm.circular_treeplot(resulting_df, output_dir="./")

This results in the following plot:

.. figure:: nagel_treeplot.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center

Let's say we are interested in generating a fatty acid profile plot for "Achromobacter denitrificans":

.. code-block:: python

	#Fatty acid profile plot
	vm.fatty_acid_profile(resulting_df, species = "Achromobacter denitrificans",  figsize=[20, 15], saveToFile=True, output_dir="./")

This results in the following plot:

.. figure:: nagel_fatty_acid.png
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center


We can also make pie plots to look at the motility of our species:

.. code-block:: python

	#Pie plot
	vm.pieplot_maker(resulting_df,"Morphology.cell morphology.motility", title="Motility for all species", saveToFile = True, output_dir="./")

This results in the following plot:

.. figure:: nagel_motility.png
    :width: 300px
    :align: center
    :height: 300px
    :alt: alternate text
    :figclass: align-center



If we are interesting in knowing from which countries the species in our dataset originate from we can create a world map:

.. code-block:: python

	#World map
	vm.worldmap_maker(resulting_df)

This results in the following plot:

.. figure:: nagel_worldmap.png
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center


Going from this world map if we want to know which country is the most frequent, we can run:

.. code-block:: python

	#Frequency plot
	vm.freqplot_maker(resulting_df, "Isolation, sampling and environmental information.isolation.country", title="Countries of origin", ylabel_name = "All countries", saveToFile=True, output_dir="./")


This results in the following plot:

.. figure:: nagel_countries.png
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center



Next, we want to make a bar plot to visualize the differences in cell width across various species:

.. code-block:: python

	#Species list for ALL species in resulting_df, not for a subset
	species_list = resulting_df["Name and taxonomic classification.species"].tolist()

	#Barplot
        vm.barplot_maker(resulting_df, "Morphology.cell morphology.cell width", "Cell width", "Width in µm", figsize=[20,10], species_list=species_list, saveToFile=True, output_dir="./")

This results in the following plot:

.. figure:: nagel_cell_width.png
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center


Let's make a boxplot which shows the optimal pH range for all the species in our dataset:

.. code-block:: python

	#Species list for ALL species in resulting_df, not for a subset
	species_list = resulting_df["Name and taxonomic classification.species"].tolist()

	#Boxplot
   	value_dict = vm.access_list_df_objects(resulting_df, "Culture and growth conditions.culture pH", "pH", pH= 1, species_list=species_list)
   	vm.boxplot_maker(value_dict, title= "Optimal pH for species", xlabel_name= "species", ylabel_name="Opt. Culture pH",figsize=[20, 10], saveToFile=True, output_dir="./")

	
This results in the following plot:

.. figure:: nagel_optimal_pH.png
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center


Lastly, since we also ran :py:func:`bacdiving.bacdive_access_for_multiple_inputs()` before, we can compare the relative abundances of e.g. the genera for our SILVA-ids.txt and our taxonomy table input in a stacked bar plot:

.. code-block:: python

	# Run for multiple inputs
	resulting_list_with_all_res_dfs = bc.bacdive_access_for_multiple_inputs({"input_data/SILVA_ids.txt" : ["input_via_file", "search_by_16S_seq_accession"], "input_data/taxtable_from_phyloseq/nagel_taxtab.tsv" : ["taxtable_input"]}, output_dir="./")
 
This results in the following plot:

.. figure:: Abundance_plot.png
    :width: 550px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center

In effect, this plot shows us the genera composition for all those species (for which BacDive information is available) in the resulting dataframe.

This concludes this tutorial for Bacdiving but feel free to use the resulting dataframe to either generate your own custom visualizations or to use it as an input for other tools!

