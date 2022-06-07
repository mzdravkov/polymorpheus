java -Xmx8g -jar snpEff.jar GRCh38.105 ~/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf > test_chrX.vcf

cat test_chrX.vcf | java -jar SnpSift.jar filter -s test_gene_names.csv "ANN[0].GENE in SET[0]"
