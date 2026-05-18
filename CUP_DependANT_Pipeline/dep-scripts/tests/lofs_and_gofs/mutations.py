import unittest
import pandas as pd

from app.lofs_and_gofs.alpha_missense import AlphaMissense
from app.lofs_and_gofs.mutations import Genes, Variant
from app.tools.columns import VepResultDf, VCF
from app.tools.file_reading import read_vcf_as_tsv
from tests.lofs_and_gofs.alpha_missense import Database
from tests.mock_classes.mock_paths import TestStaticPaths, TestPaths


class MyTestCase(unittest.TestCase):

    @staticmethod
    def setup():
        AlphaMissense.alpha_missense_isoforms = Database(TestStaticPaths.alpha_missense_isoforms)
        AlphaMissense.alpha_missense = Database(TestStaticPaths.alpha_missense)
        paths = TestPaths()
        mutations = read_vcf_as_tsv(paths.mutations_file)
        original_vcf = mutations[[VCF.CHROM, VCF.POS, VCF.REF, VCF.ALT]]
        return Genes(TestStaticPaths, original_vcf)

    def test_17_7588365_benign(self):
        df = pd.DataFrame(data=[["17",
                                 "7588365",
                                 {"stop_retained_variant", "3_prime_UTR_variant"},
                                 "ENST00000250055.3",
                                 "SOX15"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "SOX15")

        self.assertEqual(Variant.NONE, genes.classify(gene))

    def test_17_7572865_lof_consequence(self):
        df = pd.DataFrame(data=[["17",
                                 "7572865",
                                 {"splice_donor_variant"},
                                 "ENST00000293831.13",
                                 "EIF4A1"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "EIF4A1")

        self.assertEqual(Variant.LOF, genes.classify(gene))

    def test_17_7572845_missense_benign(self):
        df = pd.DataFrame(data=[["17",
                                 "7572845",
                                 {"missense_variant"},
                                 "ENST00000293831.13",
                                 "EIF4A1"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "EIF4A1")

        self.assertEqual(Variant.NONE, genes.classify(gene))

    def test_EIF4A1_canon_lof_gof_expect_lof(self):
        df = pd.DataFrame(data=[["17",
                                 "7578434",
                                 {"missense_variant"},
                                 "ENST00000293831.13",
                                 "EIF4A1"],
                                ["17",
                                 "7572865",
                                 {"splice_donor_variant"},
                                 "ENST00000293831.13",
                                 "EIF4A1"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "EIF4A1")

        self.assertEqual(Variant.LOF, genes.classify(gene))

    def test_19_1226638_missense_pathogenic_tsg_kinase_lof(self):
        df = pd.DataFrame(data=[["19",
                                 "1226638",
                                 {"missense_variant"},
                                 "ENST00000326873.11",
                                 "STK11"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "STK11")

        self.assertEqual(Variant.LOF, genes.classify(gene))

    def test_15_66490592_missense_pathogenic_onco_kinase_gof(self):
        df = pd.DataFrame(data=[["15",
                                 "66490592",
                                 {"missense_variant"},
                                 "ENST00000307102.9",
                                 "MAP2K1"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "MAP2K1")

        self.assertEqual(Variant.GOF, genes.classify(gene))

    def test_7_27165451_missense_pathogenic_tsg_onco_lof(self):
        df = pd.DataFrame(data=[["7",
                                 "27165451",
                                 {"missense_variant"},
                                 "ENST00000343483.7",
                                 "HOXA9"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "HOXA9")

        self.assertEqual(Variant.NONE, genes.classify(gene))

    def test_20_63537278_missense_pathogenic_tsg_onco_kinase_lof(self):
        df = pd.DataFrame(data=[["20",
                                 "63537278",
                                 {"missense_variant"},
                                 "ENST00000542869.3",
                                 "PTK6"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "PTK6")

        self.assertEqual(Variant.LOF, genes.classify(gene))

    def test_3_12434141_missense_pathogenic_tsg_lof(self):
        df = pd.DataFrame(data=[["3",
                                 "12434141",
                                 {"missense_variant"},
                                 "ENST00000652431.1",
                                 "PPARG"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "PPARG")

        self.assertEqual(Variant.LOF, genes.classify(gene))

    def test_7_143407672_missense_pathogenic_kinase_gof(self):
        df = pd.DataFrame(data=[["7",
                                 "143407672",
                                 {"missense_variant"},
                                 "ENST00000275815.4",
                                 "EPHA1"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "EPHA1")

        self.assertEqual(Variant.GOF, genes.classify(gene))

    def test_12_6606363_missense_pathogenic_onco_gof(self):
        df = pd.DataFrame(data=[["12",
                                 "6606364",
                                 {"missense_variant"},
                                 "ENST00000544484.6",
                                 "CHD4"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "CHD4")

        self.assertEqual(Variant.GOF, genes.classify(gene))

    def test_7_5529654_missense_pathogenic_no_group_lof(self):
        df = pd.DataFrame(data=[["7",
                                 "5529654",
                                 {"missense_variant"},
                                 "ENST00000646664.1",
                                 "ACTB"]],
                          columns=[VCF.CHROM,
                                   VCF.POS,
                                   VepResultDf.CONSEQUENCE,
                                   VepResultDf.HGVS,
                                   VepResultDf.SYMBOL])

        genes = MyTestCase.setup()
        gene = genes.get_gene(df, "ACTB")

        self.assertEqual(Variant.LOF, genes.classify(gene))


if __name__ == '__main__':
    unittest.main()
