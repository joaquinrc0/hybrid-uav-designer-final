import unittest
from scripts.utils import Utils

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.utils = Utils(file_number=2)

    def test_save_lists_from_input(self):
        list_of_lists = [["Alpha", "Cd"], ["Alpha", "CL"], [], ["Alpha", "Cm"]]
        input_strings = ["Alpha", "CL"]
        result = self.utils.save_lists_from_input(list_of_lists, input_strings)
        self.assertEqual(result, [["Alpha", "CL"], []]) 

    def test_extract_columns(self):
        data = [["1.0", "2.0"], ["3.0", "4.0"], ["5.0", "6.0"]]
        result = self.utils.extract_columns(data)
        self.assertEqual(result, [[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]])

    def test_interpolate_fun(self):
        data_dict = {"Alpha_list": [0, 1, 2], "Alpha_CL_list": [0.0, 0.5, 1.0]}
        result = self.utils.interpolate_fun(data_dict, "Alpha_list", "Alpha_CL_list", 1.5)
        self.assertAlmostEqual(result, 0.75)

    def test_round_to_significant_figures(self):
        result = self.utils.round_to_significant_figures(123.456, 3)
        self.assertEqual(result, 123)

if __name__ == "__main__":
    unittest.main()