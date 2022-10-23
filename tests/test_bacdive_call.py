import unittest
import pandas as pd

try:
    from  BacDiving.bacdive_caller import bacdive_call
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_bacdive_call(self):
        ######################### Input file ############################
        #Run bacdive_call for a small input .txt file containing 2 SILVA ids
        resulting_df = bacdive_call(input_via_file= 1, input_file_path="C:/Users/mahim/PycharmProjects/Bacdiving/input_data/SILVA_ids_2.txt", search_by_16S_seq_accession = True, output_dir="output_data/SILVA_output/", print_res_df_ToFile= 0, print_access_stats=0)
        resulting_df.index = list(resulting_df.index)

        #read in the resulting file and check if it is same as ground truth
        read_in_df = pd.read_table("C:/Users/mahim/PycharmProjects/Bacdiving/output_data/SILVA/BacdiveInformation.tsv", index_col= 0)

        #Check if number of rows and columns of the resulting dataframe is same
        self.assertEqual(len(resulting_df.columns), len(read_in_df.columns))
        self.assertEqual(len(resulting_df.columns), len(read_in_df.columns))

        #Check if column names are identical
        comp_value = True if len(read_in_df.columns) == len(set(resulting_df.columns).intersection(set(read_in_df.columns))) else False
        self.assertEqual(True, comp_value)

        ######################### Taxonomy table ############################
        resulting_df = bacdive_call(input_via_file=1,
                                    input_file_path="C:/Users/mahim/PycharmProjects/Bacdiving/input_data/taxtable_from_phyloseq/nagel_taxtab.tsv",
                                    search_by_16S_seq_accession=True, output_dir="output_data/nagel/",
                                    print_res_df_ToFile=0, print_access_stats=0)
        resulting_df.index = list(resulting_df.index)

        # read in the resulting file and check if it is same as ground truth
        read_in_df = pd.read_table("C:/Users/mahim/PycharmProjects/Bacdiving/output_data/nagel/BacdiveInformation.tsv",
                                   index_col=0)

        # Check if number of rows and columns of the resulting dataframe is same
        self.assertEqual(len(resulting_df.columns), len(read_in_df.columns))
        self.assertEqual(len(resulting_df.columns), len(read_in_df.columns))

        # Check if column names are identical
        comp_value = True if len(read_in_df.columns) == len(
            set(resulting_df.columns).intersection(set(read_in_df.columns))) else False
        self.assertEqual(True, comp_value)


if __name__ == '__main__':
    unittest.main()

