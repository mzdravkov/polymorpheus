import os
import shlex
import re
import csv
import gzip
import vcf
import pandas as pd
import numpy as np

from subprocess import PIPE
from subprocess import Popen

from enum import Enum

from config import CONFIG

import utils


ACCEPTED_REFERENCE_GENOMES = ('GRCh38', 'GRCh37', 'hg19', 'hg38')


class VCF_COLUMNS(Enum):
    CHROM  = 0
    POS    = 1
    ID     = 2
    REF    = 3
    ALT    = 4
    QUAL   = 5
    FILTER = 6
    INFO   = 7
    FORMAT = 8


class ANN_COLUMNS(Enum):
    # In case of multiple ALT fields, this helps to identify which ALT we are referring to. 
    ALT = 0

    # Annotated using Sequence Ontology terms. Multiple effects can be concatenated using '&'.
    # TODO: link effect to sequence ontology browser and show descriptions.
    EFFECT = 1

    # A simple estimation of putative impact / deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}
    IMPACT = 2

    # Common gene name (HGNC). Optional: use closest gene when the variant is "intergenic".
    GENE = 3

    # Gene ID
    GENE_ID = 4

    # Which type of feature is in the next field (e.g. transcript, motif, miRNA, etc.).
    # It is preferred to use Sequence Ontology (SO) terms, but 'custom' (user defined) are allowed.
    FEATURE_TYPE = 5

    # Depending on the annotation, this may be:
    # Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc.
    # Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID).
    FEATURE_ID = 6

    # The bare minimum is at least a description on whether the transcript is {"Coding", "Noncoding"}.
    # Whenever possible, use ENSEMBL biotypes.
    # https://uswest.ensembl.org/info/genome/genebuild/biotypes.html
    TRANSCRIPT_BIOTYPE = 7

    # Exon or Intron rank / total number of exons or introns.
    RANK_TO_TOTAL = 8

    # Variant using HGVS notation (DNA level)
    HGVS_DNA = 9

    # If variant is coding, this field describes the variant using HGVS notation (Protein level).
    # Since transcript ID is already mentioned in 'feature ID', it may be omitted here.
    HGVS_PROTEIN = 10

    # Position in cDNA and trancript's cDNA length (one based).
    CDNA_POS_TO_CDNA_LEN = 11

    # Position and number of coding bases (one based includes START and STOP codons).
    CDS_POS_TO_CDS_LEN = 12

    # Position and number of AA (one based, including START, but not STOP).
    PROT_POS_TO_PROT_LEN = 13

    # All items in this field are options, so the field could be empty.
    # Up/Downstream: Distance to first / last codon
    # Intergenic: Distance to closest gene
    # Distance to closest Intron boundary in exon (+/- up/downstream). If same, use positive number.
    # Distance to closest exon boundary in Intron (+/- up/downstream)
    # Distance to first base in MOTIF
    # Distance to first base in miRNA
    # Distance to exon-intron boundary in splice_site or splice _region
    # ChipSeq peak: Distance to summit (or peak center)
    # Histone mark / Histone state: Distance to summit (or peak center)
    DISTANCE_TO_FEATURE = 14

    # Error, warning or info code
    NOTE = 15


class VCFParsingException(Exception):
    pass


def __construct_pipe(*args):
    if len(args) == 0:
        return

    p1 = Popen(args[0], stdout=PIPE)
    if len(args) == 1:
        return p1

    p2 = None
    for arg in args[1:]:
        p2 = Popen(arg, stdin=p1.stdout, stdout=PIPE)
        p1.stdout.close()
        p1 = p2

    return p2
    

def __get_annotation_cmd(file, ref_genome):
    snpeff_path = os.path.join(CONFIG['snpEff_path'], 'snpEff.jar')
    if not os.path.isabs(snpeff_path):
         snpeff_path = os.path.join(os.getcwd(), snpeff_path)
    cmd = "java -Xmx25g -jar {} ann -noStats {} {}".format(snpeff_path, ref_genome, file)
    return shlex.split(cmd)


def __get_filter_by_genes_cmd(gene_names_file):
    snpsift = os.path.join(CONFIG['snpEff_path'], 'SnpSift.jar')
    cmd = "java -jar {} filter -s {} \"ANN[0].GENE in SET[0]\"".format(snpsift, gene_names_file)
    return shlex.split(cmd)


def __get_gene_HGNC(vcf_line):
    if match := re.search(r"ANN\=.*?\|.*?\|.*?\|(.*?)\|", vcf_line):
        return match.group(1)


def create_annotated_vcf_files_for_genes(file, ref_genome, gene_names_file):
    """
    Takes a VCF file, reference genome name and file containing one gene HGNC per line and
    parses the VCF file to annotate it and split it into a set of annotated per-gene VCF files.
    """
    commands = [
        __get_annotation_cmd(file, ref_genome),
        __get_filter_by_genes_cmd(gene_names_file)
    ]

    proc = __construct_pipe(*commands)

    files = {}

    header_lines = []

    while line := proc.stdout.readline():
        line = line.decode('utf-8')

        if len(line) == 0:
            continue

        if line[0] == '#':
            header_lines += line

        # TODO: handle multiple ANN annotations with different genes
        gene = __get_gene_HGNC(line)

        if not gene:
            continue

        if gene in files:
            files[gene].write(line)
        else:
            dest_dir = utils.get_data_dir(file)
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            files[gene] = open(os.path.join(dest_dir, gene + '.vcf'), 'w')
            # TODO: add customized header that describes our filtering
            for header_line in header_lines:
                files[gene].write(header_line)
            files[gene].write(line)

    for file in files.values():
        file.close()

    return files


def __remove_key(hashmap, key):
    del hashmap[key]
    return hashmap


def __get_vcf_row_dict(vcf_record):
    vcf_dict = vars(vcf_record)
    vcf_dict['var_type'] = vcf_record.var_type
    vcf_dict['var_subtype'] = vcf_record.var_subtype
    return vcf_dict


def parse_vcf(file):
    """
    Parses the given VCF file and returns two dataframes:
    - The first one - variants dataframe - includes all colums from the VCF
    (with the INFO column as map), but without the INFO.ANN field.
    - The second one - annotations dataframe - contains the INFO.ANN annotations,
    where each annotation is a separate column.
    Note that the annotations dataframe contains indices to the row number in the
    variants dataframe.
    """
    reader = vcf.Reader(open(file))
    df = pd.DataFrame([__get_vcf_row_dict(r) for r in reader])

    # Get snpEff annotations as dataframe.
    # Each row in the vcf may have multiple annotations
    # so each annotation keeps the index to the row (like a foreign key) 
    annotation_series = df.INFO.map(lambda info: info['ANN']).rename('ANN')
    denorm_annotations_df = pd.DataFrame(annotation_series).explode('ANN')
    denorm_annotations_df.reset_index(inplace=True)

    variation_annotation_ids = []
    i = 0
    prev = None
    for variation in denorm_annotations_df['index']:
        i += 1
        if prev not in (None, variation):
            i = 1
        prev = variation
        variation_annotation_ids.append(i)

    denorm_annotations_df.insert(1, 'gene_variation', variation_annotation_ids)
    annotations = denorm_annotations_df.ANN.str.split('|', expand=True)
    annotations_df = denorm_annotations_df.drop('ANN', axis=1).merge(annotations, left_index=True, right_index=True)

    # delete the raw ANN field from info
    df['INFO'] = df['INFO'].apply(lambda info: __remove_key(info, 'ANN'))
    # delete sample indexes as we don't need it
    df = df.drop('_sample_indexes', axis=1)

    # drop samples
    # TODO: we probably will need this in order to get some useful statistics
    # but it is hard to parse, so we can ignore it for now
    df = df.drop('samples', axis=1)

    col_names = [col.name.lower() for col in ANN_COLUMNS]
    col_names.insert(0, 'gene_variation')
    col_names.insert(1, 'variation_annotation')
    annotations_df.columns = col_names

    df.columns = [col + '_pos' if col in ('start', 'end') else col.lower() for col in df.columns]

    return df, annotations_df


def __get_matching_references(line, refs):
    matches = lambda ref: re.search(ref, line, flags=re.IGNORECASE)
    return filter(matches, refs)


def get_header_lines(file):
    header_lines = []

    input_file = None
    if file.endswith('.gz'):
        input_file = gzip.open(file, 'rt')
    else:
        input_file = open(file, 'r')

    with input_file as vcf:
        for line in vcf:
            if not line.startswith('##'):
                break
            header_lines.append(line)

    return header_lines


def validate_vcf_version(header_lines):
    fileformat_line = next((l for l in header_lines if l.startswith('##fileformat')), None)

    if not fileformat_line:
        raise VCFParsingException("Cannot find fileformat header.")

    version_match = re.match(r"##fileformat=VCFv(?P<version>[0-9.]+)", fileformat_line)
    if version_match:
        try:
            version = float(version_match.group('version'))
        except:
            raise VCFParsingException("Cannot parse ##fileformat header")

        if version < 4:
            raise VCFParsingException("VCF files with version older than 4 are not supported.")
    else:
        raise VCFParsingException("Cannot parse ##fileformat header")


def validate_and_get_genome_reference(header_lines):
    reference_line = next((l for l in header_lines if l.startswith('##reference')), None)

    if not reference_line:
        raise VCFParsingException("Cannot find reference header.")

    ref_matches = __get_matching_references(reference_line, ACCEPTED_REFERENCE_GENOMES)
    try:
        return next(ref_matches)
    except StopIteration:
        raise VCFParsingException("Unsupported genome reference version. Should be one of {}".format(ACCEPTED_REFERENCE_GENOMES))


def get_filtered_vcf_name(file_name):
    filtered = file_name.removesuffix('.vcf')
    return filtered + '_filtered.vcf'


def create_filtered_vcf_file(vcf_file):
    filtered_vcf_name = get_filtered_vcf_name(vcf_file)
    filtered_vcf = open(filtered_vcf_name, 'w')

    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#') or re.search(r"HIGH|MODERATE|LOW", line):
                filtered_vcf.write(line)

    return filtered_vcf_name


# create_annotated_vcf_files_for_genes("/home/me/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf", "GRCh38.105", "/home/me/Downloads/snpEff/test_gene_names.csv")
# print(parse_vcf('data/intermediary/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf/IL9R.vcf'))