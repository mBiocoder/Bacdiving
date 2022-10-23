import unittest

try:
    from  BacDiving.bacdive_caller import bacdive_call
    from BacDiving.visualizations_maker import access_list_df_objects
except ImportError:
    raise ImportWarning('Bacdiving not installed.')

class TestUtil(unittest.TestCase):
    def test_access_list_df_objects(self):
        #Run bacdive_call() to retrieve the resulting df
        resulting_df = bacdive_call(input_via_file=1,
                                    input_file_path="C:/Users/mahim/PycharmProjects/Bacdiving/input_data/SILVA_ids.txt",
                                    search_by_16S_seq_accession=True, output_dir="output_data/SILVA/",
                                    print_res_df_ToFile=0, print_access_stats=0)

        species_list = resulting_df["Name and taxonomic classification.species"].tolist()

        # Note: This function deals solely with temperature, pH and halophily data
        value_dict = access_list_df_objects(resulting_df, "Culture and growth conditions.culture temp", "temperature",
                                            temp=1, species_list=species_list)
        #Ground truth dict
        res_dict = {'Zeaxanthinibacter enoshimensis': [30, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 28, 29, 30, 37, 37], 'Zobellia amurskyensis': [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 24, 25, 22, 23, 24, 25], 'Zhihengliuella alba': [28, 28, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 32, 28, 37], 'Acetobacter aceti': [30, 26, 30, 30], 'Acetobacter orientalis': [30, 30, 30], 'Alishewanella fetalis': [35, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 37, 37, 37], 'Helicobacter pylori': [37, 37, 37]}

        #Check if the resulting dictionaries are identical
        self.assertDictEqual(value_dict, res_dict)


if __name__ == '__main__':
    unittest.main()





