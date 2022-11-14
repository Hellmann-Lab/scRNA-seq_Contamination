#!/bin/bash
#SBATCH --cpus-per-task=16

CELLRANGER=$1
CB_OUT=$2
EXP_CELLS=$3
TOTAL_DROP=$4
FPR=$5

cd  /opt/anaconda3/bin/
source activate
conda activate CellBender

cd /data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/Snakemake_backup

time cellbender remove-background \
	     --input $CELLRANGER/raw_feature_bc_matrix.h5 \
	     --output $CB_OUT \
	     --expected-cells $EXP_CELLS \
	     --total-droplets-included $TOTAL_DROP \
	     --fpr $FPR \
	     --epochs 150

conda deactivate