import unittest

try:
    from BacDiving.bacdive_caller import bacdive_call
    from BacDiving.visualizations_maker import get_resulting_df_values
except ImportError:
    raise ImportWarning('Bacdiving not installed.')


class TestUtil(unittest.TestCase):
    def test_get_resulting_df_values(self):
        # Run bacdive_call() to retrieve the resulting df
        resulting_df = bacdive_call(input_via_file=1,
                                    input_file_path="C:/Users/mahim/PycharmProjects/Bacdiving/input_data/SILVA_ids.txt",
                                    search_by_16S_seq_accession=True, output_dir="output_data/SILVA/",
                                    print_res_df_ToFile=0, print_access_stats=0)

        species_list = resulting_df["Name and taxonomic classification.species"].tolist()
        value_dict = get_resulting_df_values(resulting_df, plot_column="Culture and growth conditions.culture pH", plot_category="pH", species_list=species_list)
        # Ground truth dict
        res_dict = {'Zeaxanthinibacter enoshimensis': ['5.5-11', '07-08'], 'Zhihengliuella alba': ['5.0-9.0', '7']}

        # Check if the resulting dictionaries are identical
        self.assertDictEqual(value_dict, res_dict)


if __name__ == '__main__':
    unittest.main()





