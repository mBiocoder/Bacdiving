import unittest
import pandas as pd

try:
    from  BacDiving.bacdive_caller import bacdive_access_for_multiple_inputs
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_bacdive_access_for_multiple_inputs(self):
        #Note: This function only takes an input list of identical input types (either all input files or all taxonomy tables)

        ######################### Input file ############################
        #Run bacdive_call for a small input .txt file containing 2 SILVA ids and a bigger .txt file containing 10 SILVA ids
        resulting_list_with_all_res_dfs = bacdive_access_for_multiple_inputs(
            ["C:/Users/mahim/PycharmProjects/Bacdiving/input_data/SILVA_ids_2.txt",
             "C:/Users/mahim/PycharmProjects/Bacdiving/input_data/SILVA_ids.txt"])

        list_with_all_input_names = ["SILVA_ids_2", "SILVA_ids"]

        #Iterate through list and check for each dataframe inside the list if the dimensions and exact column names are identical with the ground truth data
        for elem in resulting_list_with_all_res_dfs:
            for item in list_with_all_input_names:
                elem.index = list(elem.index)
                # read in the resulting file
                read_in_df = pd.read_table("C:/Users/mahim/PycharmProjects/Bacdiving/output_data/%s/BacdiveInformation.tsv" %item, index_col=0)

                # Check if number of rows and columns of the resulting dataframe is same
                self.assertEqual(len(elem.columns), len(read_in_df.columns))
                self.assertEqual(len(elem.columns), len(read_in_df.columns))

                # Check if column names are identical
                comp_value = True if len(read_in_df.columns) == len(set(elem.columns).intersection(set(read_in_df.columns))) else False
                self.assertEqual(True, comp_value)

        ######################### Taxonomy table ############################
        # Run bacdive_call for two taxonomy tables
        resulting_list_with_all_res_dfs = bacdive_access_for_multiple_inputs(
            ["C:/Users/mahim/PycharmProjects/Bacdiving/input_data/taxtable_from_phyloseq/nagel_taxtab.tsv",
             "C:/Users/mahim/PycharmProjects/Bacdiving/input_data/taxtable_from_phyloseq/mars_taxtab.tsv"])

        list_with_all_input_names = ["nagel", "mars"]

        # Iterate through list and check for each dataframe inside the list if the dimensions and exact column names are identical with the ground truth data
        for elem in resulting_list_with_all_res_dfs:
            for item in list_with_all_input_names:
                elem.index = list(elem.index)
                # read in the resulting file
                read_in_df = pd.read_table(
                    "C:/Users/mahim/PycharmProjects/Bacdiving/output_data/%s/BacdiveInformation.tsv" % item,
                    index_col=0)

                # Check if number of rows and columns of the resulting dataframe is same
                self.assertEqual(len(elem.columns), len(read_in_df.columns))
                self.assertEqual(len(elem.columns), len(read_in_df.columns))

                # Check if column names are identical
                comp_value = True if len(read_in_df.columns) == len(
                    set(elem.columns).intersection(set(read_in_df.columns))) else False
                self.assertEqual(True, comp_value)


if __name__ == '__main__':
    unittest.main()

