# Release 3.43.0
# 2021/1/19


# X_gen.fasta is genomic DNA seq
# X_nuc.fasta is CDS DNA seq



# whole gene DNA seq
https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/A_gen.fasta
https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/B_gen.fasta
https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/C_gen.fasta

cat A_gen.fasta B_gen.fasta C_gen.fasta >hla_abc_gen.fasta

# CDS
cat A_nuc.fasta B_nuc.fasta C_nuc.fasta >hla_abc_cds.fasta

