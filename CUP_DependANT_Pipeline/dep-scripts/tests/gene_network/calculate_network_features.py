import unittest
import pandas as pd
from math import tanh, log

from tests.mock_classes.mock_paths import TestStaticPaths
from app.calculate_network_features import CalculateNetworkFeatures
from app.tools.dataframes import get_file, json_acquisition


class MyTestCase(unittest.TestCase):

    def test_convert_expression_genes(self):
        # C1orf100 should map to SPMIP3
        # COG8 is not recognised so should be removed along with corresponding value
        data = [3.211012193485512, 6.998759650910368, 10.80445373565528, 7.812562400731847, 2.6870606883398924,
                4.208673319629471, 5.004950759951045, 1.454175893185802, 11.495105510007784, 0.925999418556223,
                0.1243281350022017, 0.0]
        index = ["HDAC5", "MRPS33", "RPS19", "MAPRE1", "POLH", "COG8", "XIAP", "ANK1", "FTL", "KCNQ2", "GSTA2",
                 "C1orf100"]
        expressions = pd.Series(data=data, index=index)
        hugo_file = TestStaticPaths.hugo_conversion
        hugo = get_file(hugo_file, "Hugo conversion dict", json_acquisition)
        relevant_expressions = CalculateNetworkFeatures.convert_expression_genes(expressions, hugo)
        expected_data = [3.211012193485512, 6.998759650910368, 10.80445373565528, 7.812562400731847, 2.6870606883398924,
                         5.004950759951045, 1.454175893185802, 11.495105510007784, 0.925999418556223,
                         0.1243281350022017, 0.0]
        expected_index = ["HDAC5", "MRPS33", "RPS19", "MAPRE1", "POLH", "XIAP", "ANK1", "FTL", "KCNQ2", "GSTA2",
                          "SPMIP3"]
        self.assertEqual(expected_index, relevant_expressions.index.tolist())
        self.assertEqual(expected_data, relevant_expressions.tolist())

    def test_lofs_go_to_zero(self):
        data = [1, 1, 1, 1, 1]
        index = ["HDAC5", "MRPS33", "RPS19", "MAPRE1", "POLH"]
        lofs = {"MRPS33", "POLH"}
        expressions = pd.Series(data=data, index=index)
        cnf = CalculateNetworkFeatures(expressions, lofs, set(), TestStaticPaths)
        adjusted = cnf.get_scores()
        expected_index = ["HDAC5", "MRPS33", "RPS19", "MAPRE1", "POLH"]
        expected_data = [0.5, 0, 0.5, 0.5, 0]
        self.assertEqual(expected_index, adjusted.index.tolist())
        self.assertEqual(expected_data, adjusted.tolist())

    def test_gofs(self):
        data = [0, 0.1, 1]
        index = ["HDAC5", "MRPS33", "RPS19"]
        gofs = {"HDAC5", "MRPS33", "RPS19"}
        expressions = pd.Series(data=data, index=index)
        cnf = CalculateNetworkFeatures(expressions, set(), gofs, TestStaticPaths)
        adjusted = cnf.get_scores()
        x = 0.5 - 0.5 * (tanh(log(10)))
        expected_index = ["HDAC5", "MRPS33", "RPS19"]
        expected_data = [0, 0.5, x]
        self.assertEqual(expected_index, adjusted.index.tolist())
        self.assertEqual(expected_data, adjusted.tolist())

    def test_everything(self):

        def score(x):
            if x == 0:
                return 0
            else:
                return 0.5 - 0.5 * (tanh(log(x)))

        data = [0, 0.1, 1, 2, 3, 4, 5, 6, 7, 8]
        index = ["HDAC5", "MRPS33", "RPS19", "MAPRE1", "POLH", "COG8", "XIAP", "ANK1", "FTL", "C1orf100"]
        gofs = {"HDAC5", "POLH", "ANK1"}
        lofs = {"MRPS33", "XIAP"}
        expressions = pd.Series(data=data, index=index)
        cnf = CalculateNetworkFeatures(expressions, lofs, gofs, TestStaticPaths)
        adjusted = cnf.get_scores()
        expected_index = ["HDAC5", "MRPS33", "RPS19", "MAPRE1", "POLH", "XIAP", "ANK1", "FTL", "SPMIP3"]
        expected_data = [0, 0, score(1), score(2), score(30), 0, score(60), score(7), score(8)]
        self.assertEqual(expected_index, adjusted.index.tolist())
        self.assertEqual(expected_data, adjusted.tolist())



