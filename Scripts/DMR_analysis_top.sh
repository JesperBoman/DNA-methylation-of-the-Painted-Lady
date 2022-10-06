#!/bin/bash -l

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# Top-level script for finding DMRs. Produces all possible pairwise comparisons.
# ==========================================================================================================
# Jesper Boman                      6 oct 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Dependencies
#R - with bsseq package installed (https://www.bioconductor.org/packages/release/bioc/html/bsseq.html)



dir="" #Add working directory
declare -a arr=("HDAL" "HDLI" "LDAL")
declare -a arr2=("HDAL" "HDLI" "LDAL")



p=0
for Group1 in "${arr[@]}"
do

#Remove the focal treatment from this array
unset 'arr2[$p]'
p=$(($p+1))


for Group2 in "${arr2[@]}"
do

cd $dir



argsG1=$(ls bsseq_formatted_data | grep -E "^$Group1" | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
argsG2=$(ls bsseq_formatted_data | grep -E "^$Group2" | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')


if [ $Group1 == "LDAL" ] ; then
argsG1$(ls bsseq_formatted_data | grep -E "^$Group1" | grep -v "b.bed" | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
fi

if [ $Group2 == "LDAL" ] ; then
argsG2=$(ls bsseq_formatted_data | grep -E "^$Group2" | grep -v "b.bed" | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
fi

Rscript DMR_analysis.R $argsG1 $argsG2



mv BS_data.rda BS_data_${Group1}v${Group1}.rda
mv BS_data.fit.rda BS_data_${Group1}v${Group2}.fit.rda

mv BS_data_${Group1}v${Group2}.rda bsseq_read_rda
mv BS_data_${Group1}v${Group2}.fit.rda bsseq_fit_rda

mkdir $dir/results/${Group1}v${Group2}

cd $dir/results/${Group1}v${Group2}

Rscript ../../DMR_analysis_findDMRs.R $dir/bsseq_fit_rda/BS_data_${Group1}v${Group2}.fit.rda $dir/bsseq_read_rda/BS_data_${Group1}v${Group2}.rda


done
done
