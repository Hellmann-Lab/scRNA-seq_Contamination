#!/bin/bash
workdir=/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/Snakemake_backup
bam=$1/possorted_genome_bam.bam
REGION_VCF=/data/share/htp/CZI_kidney/mouse_experiments/experiment1_data/cellSNP/reference_vcf/reference_3strains_chr1to19.vcf.gz
cellBCs=$2
outdir=$3

sbatch -J cellsnp_sc --mem=20G --dependency afterok:143401 --cpus-per-task=19 --workdir=$workdir --wrap="srun cellsnp-lite -s $bam -b $cellBCs -O $outdir/ -p 19 -R $REGION_VCF --chrom=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 --minMAF 0.1 --minCOUNT 20 --cellTAG CB --UMItag UB"