#!/bin/bash

# LIFE748 Assessment 2 - Part 1
# Genome assembly, assembly assessment and genome annotation

# ============================================================
# 1. INPUT DATA
# ============================================================

# PacBio HiFi reads:
# Data/GN3_hifix30.fastq

# Count sequencing reads
echo "Number of reads:"
echo $(( $(wc -l < Data/GN3_hifix30.fastq) / 4 ))


# ============================================================
# 2. FLYE ASSEMBLY
# ============================================================

mkdir -p Assembly/Flye Results

/usr/bin/time -v flye \
    --pacbio-hifi Data/GN3_hifix30.fastq \
    --out-dir Assembly/Flye \
    --threads 4 \
    2> Results/flye_runtime.txt


# ============================================================
# 3. FLYE ASSEMBLY ASSESSMENT
# ============================================================

mkdir -p Results/QUAST_Flye

quast.py \
    Assembly/Flye/assembly.fasta \
    -o Results/QUAST_Flye \
    --threads 4


# ============================================================
# 4. SPADES ASSEMBLY
# ============================================================

mkdir -p Assembly/SPAdes

# A restricted k-mer series was used to accommodate the
# computational resources available.

/usr/bin/time -v spades.py \
    --isolate \
    -s Data/GN3_hifix30.fastq \
    -o Assembly/SPAdes \
    -t 4 \
    -m 5 \
    -k 21,33,55,77 \
    2> Results/spades_runtime.txt


# ============================================================
# 5. ASSEMBLY COMPARISON
# ============================================================

mkdir -p Results/QUAST_Comparison

quast.py \
    Assembly/Flye/assembly.fasta \
    Assembly/SPAdes/scaffolds.fasta \
    -o Results/QUAST_Comparison \
    --threads 4 \
    --labels Flye,SPAdes


# ============================================================
# 6. PROKKA ANNOTATION
# ============================================================

mkdir -p Annotation/Prokka

/usr/bin/time -v prokka \
    --outdir Annotation/Prokka \
    --prefix GN3_Prokka \
    --cpus 4 \
    --force \
    Assembly/Flye/assembly.fasta \
    2> Results/prokka_runtime.txt


# ============================================================
# 7. BAKTA ANNOTATION
# ============================================================

# Bakta v1.12.0 was run in a separate Conda environment.
# Bakta light database v6.0 was used.

mkdir -p Annotation/Bakta

/usr/bin/time -v bakta \
    --db ~/life748_part1/BaktaDB/db-light \
    --output Annotation/Bakta \
    --prefix GN3_Bakta \
    --threads 4 \
    --force \
    Assembly/Flye/assembly.fasta \
    2> Results/bakta_runtime.txt
