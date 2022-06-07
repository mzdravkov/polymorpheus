import os
import shlex
import re
import csv

from subprocess import PIPE
from subprocess import Popen

from enum import Enum

from config import CONFIG


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
    snpeff = os.path.join(CONFIG['snpEff_path'], 'snpEff.jar')
    cmd = "java -Xmx8g -jar {} {} {}".format(snpeff, ref_genome, file)
    return shlex.split(cmd)


def __get_filter_by_genes_cmd(gene_names_file):
    snpsift = os.path.join(CONFIG['snpEff_path'], 'SnpSift.jar')
    cmd = "java -jar {} filter -s {} \"ANN[0].GENE in SET[0]\"".format(snpsift, gene_names_file)
    return shlex.split(cmd)


def __get_gene_HGNC(vcf_line):
    if match := re.search(r"ANN\=.*?\|.*?\|.*?\|(.*?)\|", vcf_line):
        return match.group(1)


def create_annotated_vcf_files_for_genes(file, ref_genome, gene_names_file):
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
            files[gene] = open("data/intermediary/{}.vcf".format(gene), 'w')
            # TODO: add customized header that describes our filtering
            for header_line in header_lines:
                files[gene].write(header_line)
            files[gene].write(line)


def parse_vcf(file):
    data = []
    with open(file) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if len(row) == 0 or row[0][0] == '#':
                continue
            info = row[VCF_COLUMNS.INFO.value]
            annotations = re.findall(r"ANN\=.*?(?:;|$)", info)
            print(row[:9])
            for annotation in annotations:
                annotation_data = annotation[4:].split('|')
                print(annotation_data)





# create_annotated_vcf_files_for_genes("/home/me/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf", "GRCh38.105", "/home/me/Downloads/snpEff/test_gene_names.csv")
parse_vcf('data/intermediary/IL9R')