from vcf_processing import create_annotated_vcf_files_for_genes
from vcf_processing import parse_vcf

from db import get_file_by_sha, save_file, save_gene_data
from utils import sha256sum, get_data_dir


def parse(vcf_file, genes_file):
    vcf_sha = sha256sum(vcf_file)

    already_parsed = get_file_by_sha(vcf_sha)

    if already_parsed:
        data_dir = get_data_dir(vcf_file)
        print('Already parsed. Reading annotated per-gene VCFs from ' + data_dir)
    else:
        print('Saving file hash')
        save_file(vcf_file, vcf_sha)

        print('Annotating the VCF')
        # TODO: support different reference genomes
        gene_to_vcf = create_annotated_vcf_files_for_genes(vcf_file, "GRCh38.105", genes_file)
        print('Parsing the data and saving it to the database')
        for gene in gene_to_vcf:
            gene_vcf = gene_to_vcf[gene].name
            variants, annotations = parse_vcf(gene_vcf)
            
            save_gene_data(vcf_sha, gene, variants, annotations)