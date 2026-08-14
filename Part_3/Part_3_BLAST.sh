#!/bin/bash

# LIFE748 Assessment 2 - Part 3
# Identification of GN3 homologues of selected E. coli K-12 transcription factors Cra and FeaR


mkdir -p Part_3/BLAST
mkdir -p Part_3/Sequences

# Check reference protein sequences
ls -lh Part_3/Reference/Cra_K12.faa
ls -lh Part_3/Reference/FeaR_K12.faa

# Record BLAST version
blastp -version

# Build BLAST database from GN3 Bakta protein annotation
makeblastdb \
  -in Part_1/Annotation/Bakta/GN3_Bakta.faa \
  -dbtype prot \
  -out Part_3/BLAST/GN3_Bakta


# Search for GN3 homologue of Cra

blastp \
  -query Part_3/Reference/Cra_K12.faa \
  -db Part_3/BLAST/GN3_Bakta \
  -out Part_3/BLAST/Cra_vs_GN3.tsv \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -max_target_seqs 10


#Search for GN3 homologue of FeaR

blastp \
  -query Part_3/Reference/FeaR_K12.faa \
  -db Part_3/BLAST/GN3_Bakta \
  -out Part_3/BLAST/FeaR_vs_GN3.tsv \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -max_target_seqs 10


# Display top matches
echo "===== Cra top matches ====="
column -t Part_3/BLAST/Cra_vs_GN3.tsv | head

echo
echo "===== FeaR top matches ====="
column -t Part_3/BLAST/FeaR_vs_GN3.tsv | head


# Search for GN3 homologue of TreR

blastp \
  -query Part_3/Reference/TreR_K12.faa \
  -db Part_3/BLAST/GN3_Bakta \
  -out Part_3/BLAST/TreR_vs_GN3.tsv \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -max_target_seqs 10

echo
echo "===== TreR top matches ====="
column -t Part_3/BLAST/TreR_vs_GN3.tsv | head



# EXtract GN3 Cra and TreR protein sequences

awk '
/^>/ {
    keep = ($0 ~ /^>PLOCGI_01064([[:space:]]|$)/)
}
keep
' Part_1/Annotation/Bakta/GN3_Bakta.faa \
> Part_3/Sequences/GN3_Cra.faa


awk '
/^>/ {
    keep = ($0 ~ /^>PLOCGI_01344([[:space:]]|$)/)
}
keep
' Part_1/Annotation/Bakta/GN3_Bakta.faa \
> Part_3/Sequences/GN3_TreR.faa


echo
echo "===== Extracted GN3 sequences ====="

grep "^>" Part_3/Sequences/GN3_Cra.faa
grep "^>" Part_3/Sequences/GN3_TreR.faa

echo
echo "Cra amino-acid length:"
grep -v "^>" Part_3/Sequences/GN3_Cra.faa | tr -d '\n' | wc -c

echo "TreR amino-acid length:"
grep -v "^>" Part_3/Sequences/GN3_TreR.faa | tr -d '\n' | wc -c
