class VCF:
    CHROM = "#CHROM"
    POS = "POS"
    REF = "REF"
    ALT = "ALT"
    INFO = "INFO"
    ID = "ID"
    QUALITY = "QUAL"
    FILTER = "FILTER"


class ChromoGene:
    START = "Gene start (bp)"
    END = "Gene end (bp)"
    CHROM = "Chromosome/scaffold name"
    STRAND = "Strand"
    NAME = "Gene name"


class Chembl:
    DRUG = "Drug"
    MAX_PHASE = "max_phase"
    GENE = "Gene"
    UNIPROT = "UniProt"
    IC50 = "IC50"
    DOC = "supporting_doc_reference"


class GenesToDrugs:
    DRUG = "Drug"
    GENE = "Gene"
    CONFIDENCE = "confidence_score"


class GeneDf:
    CHROM = 'Chromosome'
    POS = 'Position'
    NAME = 'Gene name'
    STRAND = 'Strand'
    REF = 'Ref'
    ALT = 'Alt'
    IS_ONCO = "IsOncogene"


class VepResultDf:
    UPLOAD_VAR = "#Uploaded_variation"  # 1_396703_A/G
    LOCATION = "Location"  # 1:396703
    ALLELE = "Allele"  # G
    GENE = "Gene"  # ENSG00000290385
    FEATURE = "Feature"  # ENST00000455207
    FEAT_TYPE = "Feature_type"  # Transcript
    CONSEQUENCE = "Consequence"  # intron_variant,non_coding_transcript_variant
    C_DNA = "cDNA_position"
    CDS = "CDS_position"
    PROT_POS = "Protein_position"
    AMINO = "Amino_acids"
    CODONS = "Codons"
    EXISTING_VAR = "Existing_variation"
    EXTRA = "Extra"  # IMPACT=MODIFIER;STRAND=-1

    # acquired by breaking down the Extra column
    HGVSC = "HGVSc"
    SYMBOL = "SYMBOL"
    # acquired by splitting the HGVSc column
    HGVS = "HGVS"

    lof_consequences = {"splice_acceptor_variant",
                        "splice_donor_variant",
                        "stop_gained",
                        "frameshift_variant",
                        "stop_lost",
                        "start_lost"}
    missense_consequence = "missense_variant"


class AlphaMissenseTSV:
    CHROM = "#CHROM"  # in format chr1
    POS = "POS"
    REF = "REF"
    ALT = "ALT"
    GENOME = "genome"  # hg38
    UNIPROT = "uniprot_id"
    TRANSCRIPT = "transcript_id"  # ensembl transcript id, e.g. ENST0000102
    PROTEIN = "protein_variant"
    PATHOGENIC = "am_pathogenicity"
    CLASS = "am_class"
    # results options for am_class
    CL_PATH = "likely_pathogenic"
    CL_BENIGN = "likely_benign"
    CL_AMB = "ambiguous"


class Biomart:
    ENSEMBL_ID = "Gene stable ID"
    ENSEMBL_ID_VERS = "Gene stable ID version"
    ENSEMBL_TRANSCRIPT = "Transcript stable ID"
    TRANSCRIPT_VERS = "Transcript stable ID version"
    HUGO = "Gene name"
    SOURCE = "Source of gene name"
    VERSION = "Version(gene)"
    GENE_SOURCE = "Source(gene)"
    STRAND = "Strand"
