import tabix

from app.tools.columns import AlphaMissenseTSV
from app.tools.config_parser import config
from app.tools.logger import create_logger

logger = create_logger()


class AlphaMissense:
    logger.info("Opening alpha missense using tabix")
    alpha_missense = tabix.open(config.alpha_missense_file)
    alpha_missense_isoforms = tabix.open(config.alpha_missense_isoforms)
    nucleotide_map = {"T": "A", "A": "T", "G": "C", "C": "G"}

    def __init__(self, chrom, pos, ref, alt, versioned_transcript):
        self.chrom = AlphaMissense.chrom_format(chrom)
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.versioned_transcript = versioned_transcript

    def double_query(self):
        if self.versioned_transcript:
            logger.debug("Looking for canonical missense mutation")
            status, aa_substitution = self.single_query()
            if status:
                return status, aa_substitution
        logger.debug("Looking for isoform missense mutation")
        status, aa_substitution = self.single_query(canonical=False)
        return status, aa_substitution

    def single_query(self, canonical=True):
        logger.debug(f"AlphaMissense.single_query for {self.chrom}, {self.pos}")
        database = AlphaMissense.alpha_missense if canonical else AlphaMissense.alpha_missense_isoforms
        aa_substitution_col = 7 if canonical else 6
        transcript_col = 6 if canonical else 5
        records = database.query(self.chrom,
                                 self.pos - 1,
                                 self.pos)
        oth_ref, oth_alt = self.other_strand()
        for record in records:
            if record[-1] == AlphaMissenseTSV.CL_PATH and record[transcript_col] == self.versioned_transcript:
                if (record[2], record[3]) == (self.ref, self.alt):
                    logger.debug("Matching record found on same strand")
                    return True, record[aa_substitution_col]
                if (record[2], record[3]) == (oth_ref, oth_alt):
                    logger.debug("Matching record found on opposite strand")
                    return True, record[aa_substitution_col]
        return False, ""

    def other_strand(self):
        return AlphaMissense.nucleotide_map[self.ref], AlphaMissense.nucleotide_map[self.alt]

    @staticmethod
    def chrom_format(chrom):
        chrom = str(chrom)
        if len(chrom) < 4:
            return "chr" + chrom
        if chrom[:3] == "chr":
            return chrom
        raise ValueError(f"Invalid chromosome format: {chrom}")
