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


def base_call(bacdive_id ="", bacdive_password ="", input_via_file = 0, input_file_path =" ", search_by_id = False, search_by_culture_collection = False, search_by_taxonomy = False, search_by_16S_seq_accession = False, search_by_genome_accession = False, taxtable_input = 0, taxtable_file_path =" ", sample_name = "", print_res_df_ToFile = True, print_access_stats = True, output_dir ="./"):
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
    :param sample_name: Name of the dataset
    :param print_res_df_ToFile: If the resulting dataframe should be printed to file (as BacdiveInformation.tsv), then set print_res_df_ToFile = 1. Default is print_res_df_ToFile = 1.
    :param print_access_stats: If the resulting BacDive access statistics should be printed to file, then set print_access_stats = 1. Default is print_access_stats = 0.
    :param output_dir: Path to where resulting dataframe should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: Resulting dataframe with all strain-level BacDive information.
    """

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
            if print_access_stats == True:
                access_stats_for_file_input(total_number_silva_ids_input, access_failed, not_found,
                                            out_file_path=output_dir, sample_name = sample_name)

            if print_res_df_ToFile == True:
                resulting_df_to_file_printer(resulting_df, output_dir=output_dir, sample_name = sample_name)

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
            df_new.to_csv(output_dir + '%s_Species_names_from_taxtable_file.csv' %sample_name, index=False, header=None)
            # read output file and return this dataframe
            try:
                downstream_df = pd.read_csv(output_dir + "%s_Species_names_from_taxtable_file.csv" %sample_name, header=None)
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
                if print_access_stats == True:
                    access_stats_for_taxtable_input(final_df, selected_rows, total_number_species_input, access_failed,
                                                    not_found, out_file_path=output_dir, sample_name = sample_name)
                if print_res_df_ToFile == True:
                    resulting_df_to_file_printer(resulting_df, output_dir=output_dir, sample_name = sample_name)

                return resulting_df
        else:
            print(
                'If you have not registered for Bacdive web services yet, please do so using the following link before running this package: ' + 'https://sso.dsmz.de/auth/realms/DSMZ/protocol/openid-connect/auth?response_type=code&redirect_uri=https%3A%2F%2Fapi.bacdive.dsmz.de%2Flogin&client_id=api.bacdive&nonce=fc4f465de78388722385d9bd3f82de1f&state=cab5a9d249f9502bd01f7715cc5d99bd&scope=openid')
            print('Only try calling this function after gaining access.')
            exit(1)
    except AttributeError:
        print("Your login credentials were wrong. Please try again!")
        sys.exit(1)


def flattened_full_file_maker(resulting_df, columns_of_interest = [""], sample_name = "", taxtable_file_path = "", output_dir = "./"):
    # Write columns to file of the following structure: taxonomic ranks -> # of strains per species found on Bacdive -> column of interest flattened
    out_file_path = output_dir + '%s_Flattened_Bacdive_data.tsv' %sample_name
    file = open(out_file_path, "w")
    file.write(" " + "\t" + "Kingdom" + "\t" + "Phylum" + "\t" + "Class" + "\t" + "Order" + "\t" + "Family" + "\t" + "Genus" + "\t" + "Species" + "\t" + "Number of strains" + "\t" + "\t".join(str(interested_col) for interested_col in columns_of_interest) + "\n")

    #Iterate throguh taxtable and if species unknown or not in Bacdive, fill everything with NA, else stored respective flattened values for the columns of interest
    final_df = pd.read_table(taxtable_file_path, index_col=0)
    for ind in final_df.index:
        kingdom = str(final_df['Kingdom'][ind])
        phylum = str(final_df['Phylum'][ind])
        class_rank = str(final_df['Class'][ind])
        order = str(final_df['Order'][ind])
        family = str(final_df['Family'][ind])
        genus = str(final_df['Genus'][ind])
        species = str(final_df['Species'][ind])

        full_species =  genus + " " + species

        # If full species name is contained in resulting dataframe, then get strain count and flattened values
        if (resulting_df["Name and taxonomic classification.species"] == full_species).any() == True:
            strain_number = len(resulting_df.loc[resulting_df["Name and taxonomic classification.species"] ==  full_species].index)
            list_results_interested_columns = []
            # Get majority of non-na values and nan if information for all strains is nan
            for interested_col in columns_of_interest:
                first_cat_col = (resulting_df.loc[resulting_df["Name and taxonomic classification.species"] == full_species][interested_col]).dropna().tolist()
                if len(first_cat_col) != 0:
                    value = max(set(first_cat_col), key=first_cat_col.count)
                    list_results_interested_columns.append(value)
                else: # if no information is there for any strain of a given spieces then add nan
                    value = "nan"
                    list_results_interested_columns.append(value)
            # Write to file
            file.write(str(ind) + "\t" + kingdom + "\t" + phylum + "\t" + class_rank + "\t" + order + "\t" + family + "\t" + genus + "\t" + species + "\t" + str(strain_number) + "\t" + "\t".join(str(list_results_interested_columns[item]) for item in range(len(list_results_interested_columns))) + "\n")
        # If not contained in bacdive or species is unknown, then fill up all following columns with NA's
        else:
            file.write(str(ind) + "\t" + kingdom + "\t" + phylum + "\t" + class_rank + "\t" + order + "\t" + family + "\t" + genus + "\t" + species + "\t" + "nan" + "\t" + "\t".join("nan" for item in range(len(columns_of_interest))) + "\n")
    file.close()


def access_stats_for_file_input(total_number_inputs, access_failed, not_found, sample_name ="", out_file_path ="./"):
    # Prints access statistics for input file types
    access_worked = total_number_inputs - access_failed
    percentage_found = access_worked / total_number_inputs * 100
    out_file_path =  out_file_path + '%s_Access_stats_file_input.tsv' %sample_name
    file = open(out_file_path, "w")
    file.write("-- Access statistics --" + "\n")
    file.write("-> Total number of items are: " + str(total_number_inputs) + ", out of which " + str(round(percentage_found, 2)) + "% were found on Bacdive. Therefore, " + str(access_failed) + " items were not found in Bacdive." + "\n")
    file.write("The following items were not found on Bacdive: " + "\n")
    for i in not_found:
        file.write(i + "\n")
    file.close()


def access_stats_for_taxtable_input(final_df, selected_rows, total_number_species_input, access_failed, not_found, sample_name = "", out_file_path ="./"):
    # Prints access statistics for input taxonomy tables
    access_worked = total_number_species_input - access_failed
    percentage_found = access_worked / total_number_species_input * 100
    out_file_path = out_file_path + '%s_Access_stats_taxtable_input.tsv' %sample_name
    file = open(out_file_path, "w")
    file.write("-- Access statistics --" + "\n")
    file.write("-> Your input taxtable file contains " + str(len(final_df)) + " rows. After removing NaN species-values from the dataframe we are left with " + str(len(selected_rows)) + " rows." + "\n")
    file.write("-> Total number of rows are: " + str(total_number_species_input) + ", out of which " + str(round(percentage_found, 2)) + "% were found on Bacdive. Therefore, " + str(access_failed) + " species were not found in Bacdive." + "\n")
    file.write("\n")
    file.write("The following species were not found on Bacdive: " + "\n")
    for i in not_found:
        file.write(i + "\n")
    file.close()


def resulting_df_to_file_printer(resulting_df,sample_name = "", output_dir = "./"):
    # Writes resulting dataframe with all BacDive information to file
    resulting_df.to_csv(output_dir + '%s_BacdiveInformation.tsv' %sample_name, sep='\t', encoding='utf-8', index=True)


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
