#!/bin/bash

# Variables
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
REF_SDF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf
TRUTH_VCF=HG00438.GRCh38_no_alt.deepvariant.vcf.gz


# Inputs
SAMPLE=$1
QUERY_VCF=$2
REGIONS_BED=$3
STRATIFICATION_TSV=$4
OUTPUT_DIR=$5
THREADS=$6


echo "VCF renaming"
zcat "$QUERY_VCF" | sed 's/^grch38#//g' | bgzip -c > "$QUERY_VCF".renamed.vcf.gz && tabix "$QUERY_VCF".renamed.vcf.gz


echo "VCF processing"
bash vcf_preprocess.sh "$QUERY_VCF".renamed.vcf.gz "$SAMPLE" 50
rm "$QUERY_VCF".renamed.vcf.gz

FNAME=$(basename "$QUERY_VCF".renamed.vcf.gz)
PREFIX="${FNAME%.vcf.gz}"
NORMALIZED_VCF=${PREFIX}.max50.chr1-22.vcf.gz

echo "VCF evaluation"
mkdir -p "$OUTPUT_DIR"

OUTPUT_PREFIX="$OUTPUT_DIR"/"$SAMPLE"

docker run -v "${PWD}"/:/data paramost/hap.py /opt/hap.py/bin/hap.py \
     --threads "$THREADS"                             \
     -r data/$REF                                     \
     -o data/"$OUTPUT_PREFIX"                         \
     -f data/"$REGIONS_BED"                           \
     --stratification data/"$STRATIFICATION_TSV"      \
     --no-leftshift                                   \
     --no-decompose                                   \
     --engine vcfeval                                 \
     --engine-vcfeval-template data/$REF_SDF          \
     data/$TRUTH_VCF                                  \
     data/"$NORMALIZED_VCF"