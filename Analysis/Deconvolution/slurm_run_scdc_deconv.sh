#!/bin/bash

# settings
base="/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution_new/"
slurms="$base/slurm";


declare -a scripts=("CAST_multiref" "CAST_singleref")

for s in "${scripts[@]}"
do

cat $base/scripts/bash/slurm_deconv_"$s".txt | while read ID sce eset single goodID contID emptyID experiment assay pseudobulk outpath outname;

do
# MAKING THE HEADER
echo '#!/bin/bash' > "$ID"_"$s".sh
echo '#SBATCH -n 1' >> "$ID"_"$s".sh
echo '#SBATCH --workdir='$slurms'' >> "$ID"_"$s".sh
echo '#SBATCH --error='$slurms/"$ID"_"$s"'.%J.err' >>"$ID"_"$s".sh
echo '#SBATCH --output='$slurms/"$ID"_"$s"'.%J.out' >> "$ID"_"$s".sh
echo '#SBATCH --cpus-per-task=2' >> "$ID"_"$s".sh
echo '#SBATCH --mem=10000' >> "$ID"_"$s".sh
# echo '#SBATCH --nodelist=gorilla1' >> $ID.sh
# echo '#SBATCH --exclude=gorilla1' >> $ID.sh
echo '#SBATCH -p normal' >> "$ID"_"$s".sh

# THE ACTUAL COMMANDS
echo "echo $outname $ID" >> "$ID"_"$s".sh
echo "Rscript $base/scripts/Rscript/run_deconv_scdc.R --sce=$sce --eset=$eset --single=$single --goodID=$goodID --contID=$contID --emptyID=$emptyID --experiment=$experiment --assay=$assay --pseudobulk=$pseudobulk --outpath=$outpath --outname=$outname" >> "$ID"_"$s".sh

sbatch "$ID"_"$s".sh

done
done

