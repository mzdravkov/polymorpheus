from vcf_processing import create_annotated_vcf_files_for_genes
from vcf_processing import parse_vcf

def parse(vcf_file, genes_file):
    # TODO: support different reference genomes
    gene_to_vcf = create_annotated_vcf_files_for_genes(vcf_file, "GRCh38.105", genes_file)
    for gene in gene_to_vcf:
        vcf = gene_to_vcf[gene].name
        print(parse_vcf(vcf))