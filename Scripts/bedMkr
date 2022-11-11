#!/bin/bash -l

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to transform Bismark-outputted Stranded CpG reports to gzipped .bed-files
# ==========================================================================================================
# Jesper Boman                      9 nov 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #


dir="/crex/proj/uppstore2017185/b2014034_nobackup/Jesper/VanessaDNAm/stranded_CpG_reports"
ls "$dir" > sample_list

while IFS= read -r samp
do
	shortname=$(echo $samp | sed -E 's/TJ-2719-(.*)_S.*gz/\1/g' )
	zcat $dir/${samp} | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}'  > ${shortname}.bed
	gzip ${shortname}.bed
	mv ${shortname}.bed.gz ../bed_wgbs_data
done < "sample_list"
