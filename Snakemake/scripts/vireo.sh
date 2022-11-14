#!/bin/bash

REF=/data/share/htp/CZI_kidney/mouse_experiments/experiment1_data/cellSNP/reference_vcf/reference_3strains_chr1to19.vcf.gz
CELLSNP_PATH=$1
VIREO_PATH=$2
VIREO_VCF=$3

/opt/bin/bcftools view $REF -R $CELLSNP_PATH/cellSNP.base.vcf > $VIREO_VCF
gzip -k $CELLSNP_PATH/cellSNP.base.vcf
vireo -c $CELLSNP_PATH -d $VIREO_VCF -o $VIREO_PATH -t GT
