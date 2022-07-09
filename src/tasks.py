import os
import pysam
from datetime import datetime
from vcf_processing import create_annotated_vcf_files_for_genes
from vcf_processing import parse_vcf
from vcf_processing import get_header_lines
from vcf_processing import validate_vcf_version
from vcf_processing import validate_and_get_genome_reference
from vcf_processing import create_filtered_vcf_file

from db import get_file_by_sha, save_file, save_gene_data, update_file_status
from utils import sha256sum, get_data_dir


def __get_snpeff_genome_reference(genome_reference):
    if genome_reference.startswith('GRCh38'):
        return 'GRCh38.105'
    elif genome_reference.startswith('GRCh37'):
        return 'GRCh37.75'
    return genome_reference


def parse(vcf_file, genes_file):
    vcf_sha = sha256sum(vcf_file)

    existing_row = get_file_by_sha(vcf_sha)

    if existing_row and existing_row['status'] == 'processed':
        data_dir = get_data_dir(vcf_file)
        print('Already parsed. Reading annotated per-gene VCFs from ' + data_dir)
    else:
        header = get_header_lines(vcf_file)
        validate_vcf_version(header)
        genome_reference = validate_and_get_genome_reference(header)
        snpeff_ref = __get_snpeff_genome_reference(genome_reference)

        if not existing_row:
            print('Saving file hash')
            save_file(os.path.basename(vcf_file), vcf_sha, vcf_file, datetime.now())

        print('Annotating the VCF')
        gene_to_vcf = create_annotated_vcf_files_for_genes(vcf_file, snpeff_ref, genes_file)
        print('Parsing the data and saving it to the database')
        for gene in gene_to_vcf:
            gene_vcf = gene_to_vcf[gene].name

            filtered_vcf = create_filtered_vcf_file(gene_vcf)

            variants, annotations = parse_vcf(gene_vcf)
            save_gene_data(vcf_sha, gene, variants, annotations)

            # create a tabix index for the vcf and filtered vcf
            # and compress them with gzip
            pysam.tabix_index(gene_vcf, preset='vcf', force=True)
            pysam.tabix_index(filtered_vcf, preset='vcf', force=True)
        update_file_status(vcf_sha, 'processed')
    print('Processed ' + vcf_file)