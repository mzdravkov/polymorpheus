import os
from datetime import datetime
from vcf_processing import create_annotated_vcf_files_for_genes
from vcf_processing import parse_vcf

from db import get_file_by_sha, save_file, save_gene_data, update_file_status
from utils import sha256sum, get_data_dir


def parse(vcf_file, genes_file):
    vcf_sha = sha256sum(vcf_file)

    existing_row = get_file_by_sha(vcf_sha)

    if existing_row and existing_row['status'] == 'processed':
        data_dir = get_data_dir(vcf_file)
        print('Already parsed. Reading annotated per-gene VCFs from ' + data_dir)
    else:
        if not existing_row:
            print('Saving file hash')
            save_file(os.path.basename(vcf_file), vcf_sha, vcf_file, datetime.now())

        print('Annotating the VCF')
        # TODO: support different reference genomes
        gene_to_vcf = create_annotated_vcf_files_for_genes(vcf_file, "GRCh38.105", genes_file)
        print('Parsing the data and saving it to the database')
        for gene in gene_to_vcf:
            gene_vcf = gene_to_vcf[gene].name
            variants, annotations = parse_vcf(gene_vcf)
            
            save_gene_data(vcf_sha, gene, variants, annotations)
        update_file_status(vcf_sha, 'processed')