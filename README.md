# Bacdiving

Bacdiving accesses and retrieves information from the world's largest database for standardized bacterial phenotypic information: BacDive.
Additionally, Bacdiving provides several options to visualize this information.  

Before using Bacdiving please register (for free) on [BacDive](https://api.bacdive.dsmz.de/).
Using your BacDive credentials you can dive into Bacdiving. 

In general, Bacdiving can deal with two types of input data: a taxonomy table (e.g. as extracted from a phyloseq-object) or an input file (.csv, .txt, .tsv) with one query-type per row.
Possible BacDive query types include: BacDive id, taxonomy (as in species name), 16S sequencing accession id (e.g. SILVA id), culture collection accession id or genome sequence accession id.
However, the input file should be consistant with only contain one (!) query type for all of its rows.

Here is a minimal example on how to use Bacdiving, please refer to the full documentationfor more details:

```
from bacdiving import bacdive_caller as bc
from bacdiving import treeplots_maker as tm
from bacdiving import visualizations_maker as vm

### Retrieve and access information stored on BacDive ###
#SILVA id query
resulting_df = bc.bacdive_call(input_via_file=1, input_file_path="./input_data/SILVA_ids.txt",search_by_16S_seq_accession=True, print_res_df_ToFile=1, print_access_stats=1)

#Taxonomy query
resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="./input_data/Taxonomy_ids.txt", search_by_taxonomy = True, output_dir="./output_data/SILVA/", print_res_df_ToFile= 1, print_access_stats=1)

#Bacdive id query
resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="./input_data/Bacdive_ids.txt", search_by_id = True, taxtable_input = 0, taxtable_file_path='./input_data/taxtable_from_phyloseq/mars_taxtab.tsv', output_dir="./", print_res_df_ToFile= 1)

#Culture collection query
resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="Culture_col_ids.txt", search_by_culture_collection = True, taxtable_input = 0, output_dir="./", print_res_df_ToFile= 1)

#Genome accession query
resulting_df = bc.bacdive_call(input_via_file= 1, input_file_path="Genome_accession_ids.txt", search_by_genome_accession = True, output_dir="./", print_res_df_ToFile= 1)

#Taxonomy table input (e.g. as extracted from phyloseq-object)
resulting_df = bc.bacdive_call(taxtable_input = 1, taxtable_file_path='./input_data/taxtable_from_phyloseq/taxtab.tsv', output_dir="./", print_res_df_ToFile= 1, print_access_stats=1)

# Run for multiple inputs (of possibly different input types)
resulting_list_with_all_res_dfs = bc.bacdive_access_for_multiple_inputs(input_lists={"./input_data/SILVA_ids.txt" : ["input_via_file", "search_by_16S_seq_accession"], "./input_data/taxtable_from_phyloseq/taxtab.tsv" : ["taxtable_input"]})
```

```
### Some possible visualizations ###

#Relative abundance plot
vm.stacked_barplot_relative_abundance(resulting_list_with_all_res_dfs, sample_names=["Silva_input", "Taxtab_input"], plot_column="Name and taxonomic classification.genus", title="Relative abundance", saveToFile = True, output_dir="./")

#Tree plots
tm.overview_treeplot(resulting_df, label_name1="Temperature", label_name2="Oxygen tolerance", saveToFile=True, output_dir="./")
tm.circular_treeplot(resulting_df, output_dir="./")

#Fatty acid profile plot
vm.fatty_acid_profile(resulting_df, species = "Achromobacter denitrificans",  figsize=[20, 15], saveToFile=True, output_dir="./")

#Pie plot
vm.pieplot_maker(resulting_df,"Morphology.cell morphology.motility", title="Motility for all species", saveToFile = True, output_dir="./")

#World map
vm.worldmap_maker(resulting_df)

#Frequency plot
vm.freqplot_maker(resulting_df, "Isolation, sampling and environmental information.isolation.country", title="Countries of origin", ylabel_name = "All countries", saveToFile=True, output_dir="./")

#Species list for ALL species in resulting_df, not for a subset
species_list = resulting_df["Name and taxonomic classification.species"].tolist()

#Barplot
vm.barplot_maker(resulting_df, "Sequence information.GC content.GC-content", "GC-content", "GC-content", figsize=[40,10],  species_list=species_list, saveToFile=True, output_dir="./")

#Boxplot
value_dict = vm.access_list_df_objects(resulting_df, "Culture and growth conditions.culture temp", "temperature", temp= 1, species_list=species_list)
vm.boxplot_maker(value_dict, title= "Optimal temperature for species", xlabel_name= "species", figsize=[20, 10], ylabel_name="Opt. Culture Temp. $C^{o}$", saveToFile=True, output_dir="./")
```
