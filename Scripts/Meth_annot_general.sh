#!/bin/bash -l
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to calculate methylation level per feature in .bed-files
# ==========================================================================================================
# Jesper Boman                      9 nov 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #



module load bioinfo-tools BEDTools/2.29.2 R_packages/4.1.1

dir="/VanessaDNAm/bed_wgbs_data"

mkdir various_meth_bedfiles

while IFS= read -r annot_double
do

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annotation=$(cut -f2 <(echo $annot_double) -d ",")

inputDir="various_meth_bedfiles"

ls -l $dir| awk 'NR!=1{print $9}' | sed 's/.bed.gz//g'  > sample_list

awk '{print $1 "\t" $2 "\t" $3 "\t" $1 "_" $2 "_" $3}' $annotation| sort -k1,1 -k2,2n  > tmp_annot
annotation=tmp_annot


p=0

while IFS= read -r sample
do
	


if [ ! -f  ${inputDir}/temp_${shortAnnot}_${sample}.bed.gz ]; then
        bedtools intersect -a <(zcat $dir/${sample}.bed.gz | awk '{if($7=="CG"){print $0}}') \
        -b $annotation -wa -wb > $inputDir/temp_${shortAnnot}_${sample}.bed &
fi


#Parallelization code: make 5 bedfiles at a time
p=$(($p+1))
remainder=$(( p  % 5 ))
if [[ $remainder -eq 0 ]] ; then
wait
fi


done < "sample_list"
wait


##CLEAN UP FROM PREVIOUS RUN - HASHTAG TO SKIP##
rm meth_${shortAnnot}_data
################################################



while IFS= read -r sample
do

ind=$(basename -s ".bed.gz" $sample)



#Gzip the files for future use
gzip $inputDir/temp_${shortAnnot}_${sample}.bed




#DNA methylation per annotation feature - 2/12 - 2020
#Output: Chromosome Start End  Total_CpG Post_filter_CpG Methylated_reads_per_feature Unmethylated_reads_per_feature  Sum_of_proportions_of_methylated_over_total_reads_per_CpG Mean_per_dinuc (i.e. Methylation level) Individual Tissue




#annot=$9 "\t" $10 "\t" $11 - for methylation per single feature per sample
#annot="Feature" - for average across all features of a sample


awk -v ind="$ind" 'BEGIN {OFS="\t"} {if(!(annot in pos)){skip[annot]=="FALSE"}; annot=$9 "\t" $10 "\t" $11 "\t" $12; 
if(pos[annot]==$2 && $7=="CG" && prev=="CG" && NR!=1 && skip[annot]=="FALSE"){ 
        if(($5+$6+prev_total_reads[annot])>5 && ($5+$6+prev_total_reads[annot])<200){skip[annot]="TRUE";
        a[annot]+=1;  b[annot]+=($5+prev_meth_reads[annot]); c[annot]+=($6+prev_unmeth_reads[annot]); rat[annot]+=(($5+prev_meth_reads[annot])/($5+$6+prev_total_reads[annot])) }; 
}
    else{skip[annot]="FALSE"}; 
prev=$7; if(skip[annot]=="FALSE"){pos[annot]=$3} ; 
if($7=="CG"){
        totCpG[annot]+=0.5; prev_meth_reads[annot]=$5; prev_unmeth_reads[annot]=$6; prev_total_reads[annot]=$5+$6}
}
    
    END{for (wind in totCpG){if (a[wind]==""){print wind, totCpG[wind], "NA", "NA", "NA", "NA", "NA", ind}
    else{print wind, totCpG[wind], a[wind], b[wind], c[wind], rat[wind], rat[wind]/a[wind], ind}}}'  \
<(zcat ${inputDir}/temp_${shortAnnot}_${sample}.bed.gz)  >> meth_${shortAnnot}_data


done < "sample_list"


rm tmp_annot
done < "annotation.list.csv"
