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


def base_call(bacdive_id ="", bacdive_password ="", input_via_file = 0, input_file_path =" ", search_by_id = False, search_by_culture_collection = False, search_by_taxonomy = False, search_by_16S_seq_accession = False, search_by_genome_accession = False, taxtable_input = 0, taxtable_file_path =" ", sample_name = "", print_res_df_ToFile = 1, print_access_stats = 0, output_dir ="./"):
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
            if print_access_stats == 1:
                access_stats_for_silva_input(total_number_silva_ids_input, access_failed, not_found,
                                             out_file_path=output_dir, sample_name = sample_name)

            if print_res_df_ToFile == 1:
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
            df_new.to_csv(output_dir + '%s Species_names_from_taxtable_file.csv' %sample_name, index=False, header=None)
            # read output file and return this dataframe
            try:
                downstream_df = pd.read_csv(output_dir + "%s Species_names_from_taxtable_file.csv" %sample_name, header=None)
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
                                                    not_found, out_file_path=output_dir, sample_name = sample_name)
                if print_res_df_ToFile == 1:
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


def access_stats_for_silva_input(total_number_silva_ids_input, access_failed, not_found, sample_name = "", out_file_path = "./"):
    # Prints access statistics for input file types
    access_worked = total_number_silva_ids_input - access_failed
    percentage_found = access_worked / total_number_silva_ids_input * 100
    out_file_path =  out_file_path + '%s access_stats_SILVA_input.tsv' %sample_name
    file = open(out_file_path, "w")
    file.write("-- Access statistics --" + "\n")
    file.write("-> Total number of SILVA ids are: " + str(total_number_silva_ids_input) + ", out of which " + str(round(percentage_found, 2)) + "% were found on Bacdive. Therefore, " + str(access_failed) + " SILVA ids were not found in Bacdive." + "\n")
    file.write("The following SILVA ids were not found on Bacdive: " + "\n")
    for i in not_found:
        file.write(i + "\n")
    file.close()


def access_stats_for_taxtable_input(final_df, selected_rows, total_number_species_input, access_failed, not_found, sample_name = "", out_file_path ="./"):
    # Prints access statistics for input taxonomy tables
    access_worked = total_number_species_input - access_failed
    percentage_found = access_worked / total_number_species_input * 100
    out_file_path = out_file_path + '%s access_stats_taxtable_input.tsv' %sample_name
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
    resulting_df.to_csv(output_dir + '%s BacdiveInformation.tsv' %sample_name, sep='\t', encoding='utf-8', index=True)


def bacdive_call(bacdive_id ="", bacdive_password ="", input_lists = {}, sample_names = [], print_res_df_ToFile = 1, print_access_stats = 1, output_dir ="./"):
    """
    For multiple input files (either all input files or all taxonomy tables) this function reads the input, queries the BacDive database and stores resulting dataframe(s) and access statistics.

    :param bacdive_id: Log in credential: BacDive id. Default is empty string. If entering the log-in credentials multiple times (for each input sample) is not desirable, then it is recommended to input the credtials once as function parameters.
    :param bacdive_password: Log in credential: BacDive id. Default is empty string. If entering the log-in credentials multiple times (for each input sample) is not desirable, then it is recommended to input the credtials once as function parameters.
    :param input_lists: Dictionary which specifies {<file-path> : <file-type>}. In case of input_via_file as the input type, then the exact content types has to be specified. The dictionary keys are thus the respective file paths; the dictionary values are lists with two elements. The first list element specifies the file type: "taxtable_input" or "input_via_file"; the second list element is only neccessary if the file type is "input_via_file" in order to further specify the input type of this file: "search_by_id", "search_by_culture_collection", "search_by_taxonomy", "search_by_16S_seq_accession" or "search_by_genome_accession". For example: {"input_data/SILVA_ids.txt" : ["input_via_file", "search_by_16S_seq_accession"], "./input_data/taxtable_from_phyloseq/nagel_taxtab.tsv" : ["taxtable_input"]}
    :param sample_names: List of samples names
    :param output_dir: Path to where resulting dataframe should be saved if saveToFile is set to True. Default is current directory ("./").
    :return: List containing the resulting dataframe(s)  with all strain-level BacDive information for all those multiple inputs.
    """
    print('If you have not registered for Bacdive web services yet, please do so using the following link before running this package: ' + 'https://sso.dsmz.de/auth/realms/DSMZ/protocol/openid-connect/auth?response_type=code&redirect_uri=https%3A%2F%2Fapi.bacdive.dsmz.de%2Flogin&client_id=api.bacdive&nonce=fc4f465de78388722385d9bd3f82de1f&state=cab5a9d249f9502bd01f7715cc5d99bd&scope=openid')

    resulting_list_with_all_res_dfs = []
    for file, filetype in input_lists.items():
        index = list(input_lists).index(file)
        sample_name = sample_names[index]

        if filetype[0] == "taxtable_input":
            resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, taxtable_input=1, taxtable_file_path=file, sample_name=sample_name, output_dir=output_dir,
                                     print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
            resulting_list_with_all_res_dfs.append(resulting_df)
        elif filetype[0] == "input_via_file":
            if filetype[1] == "search_by_id":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_id= True, sample_name=sample_name, output_dir=output_dir,
                                         print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_culture_collection":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_culture_collection=True, sample_name=sample_name,
                                         output_dir=output_dir, print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_taxonomy":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_taxonomy=True,
                                         output_dir=output_dir, sample_name=sample_name,
                                         print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_16S_seq_accession":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_16S_seq_accession=True,
                                         output_dir=output_dir, sample_name=sample_name,
                                         print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)

            elif filetype[1] == "search_by_genome_accession":
                resulting_df = base_call(bacdive_id=bacdive_id, bacdive_password=bacdive_password, input_via_file=1, input_file_path=file, search_by_genome_accession=True,
                                         output_dir=output_dir, sample_name=sample_name, print_res_df_ToFile=print_res_df_ToFile, print_access_stats=print_access_stats)
                resulting_list_with_all_res_dfs.append(resulting_df)
        else:
            print("Please check if the filetype, sample names and search parameters are set correctly!")
            exit(1)

    return resulting_list_with_all_res_dfs
