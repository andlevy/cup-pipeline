import unittest
import pandas as pd

from app.lofs_and_gofs.alpha_missense import AlphaMissense
from tests.mock_classes.mock_paths import TestStaticPaths, TestPaths


class Database:

    def __init__(self, file):
        self.database = pd.read_csv(file, sep="\t")

    def query(self, chrom, start, end):
        records = self.database.loc[(self.database["1"] == chrom) & (self.database["2"] == (start+1))]
        return records.values


class MyTestCase(unittest.TestCase):
    paths = TestPaths()

    def test_19_1226638_pathogenic_isoform(self):
        alpha_missense = AlphaMissense("chr19",
                                       1226638,
                                       "G",
                                       "T",
                                       "ENST00000586243.5")
        AlphaMissense.alpha_missense_isoforms = Database(TestStaticPaths.alpha_missense_isoforms)
        AlphaMissense.alpha_missense = Database(TestStaticPaths.alpha_missense)

        self.assertEqual((True, 'K430N'), alpha_missense.double_query())

    def test_15_66490592_missing_isoform(self):
        alpha_missense = AlphaMissense("chr15",
                                       66490592,
                                       "C",
                                       "T",
                                       "ENST00000307102.10")
        AlphaMissense.alpha_missense_isoforms = Database(TestStaticPaths.alpha_missense_isoforms)
        AlphaMissense.alpha_missense = Database(TestStaticPaths.alpha_missense)

        self.assertEqual(alpha_missense.double_query(), (False, ''))

    def test_7_27165453_pathogenic_canonical(self):
        alpha_missense = AlphaMissense("chr7",
                                       27165453,
                                       "G",
                                       "T",
                                       "ENST00000343483.7")
        AlphaMissense.alpha_missense_isoforms = Database(TestStaticPaths.alpha_missense_isoforms)
        AlphaMissense.alpha_missense = Database(TestStaticPaths.alpha_missense)

        self.assertEqual((True, 'A2D'), alpha_missense.double_query())

    def test_12_6606363_benign(self):
        alpha_missense = AlphaMissense("chr12",
                                       6606363,
                                       "C",
                                       "A",
                                       "ENST00000544040.7")
        AlphaMissense.alpha_missense_isoforms = Database(TestStaticPaths.alpha_missense_isoforms)
        AlphaMissense.alpha_missense = Database(TestStaticPaths.alpha_missense)

        self.assertEqual((False, ''), alpha_missense.double_query())

    def test_missing_data_benign(self):
        alpha_missense = AlphaMissense("chr20",
                                       664,
                                       "C",
                                       "A",
                                       "ENST00000541350.3")
        AlphaMissense.alpha_missense_isoforms = Database(TestStaticPaths.alpha_missense_isoforms)
        AlphaMissense.alpha_missense = Database(TestStaticPaths.alpha_missense)

        self.assertEqual((False, ''), alpha_missense.double_query())


if __name__ == '__main__':
    unittest.main()
