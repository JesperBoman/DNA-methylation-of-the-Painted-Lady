#!/bin/bash -l

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==============================================================================================================================
# A bootstrapping method to calculate the empirical p-value of base pair overlap between DMRs and a bed file of e.g. genes
# ==============================================================================================================================
# Jesper Boman                      6 oct 2022
# ==============================================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


## DMR overlap analysis ##

#mkdir annotation_overlap

ml bioinfo-tools BEDTools/2.29.2
faidx="fasta index file"
dir="./results"


while IFS= read -r annot_double
do

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annot=$(cut -f2 <(echo $annot_double) -d ",")

#Bedtools merge here is a convenient heuristic since annotation features can overlap but are bound to a specific gene such as is the case for UTRs
#Genes can also overlap and it is probably fine to consider them as one long gene region instead of counting them as a double overlap.
awk '{print $1 "\t" $2 "\t" $3}' $annot | sort -k1,1 -k2,2n | bedtools merge > tmp_annot_sp
annot=tmp_annot_sp

annot_len=$(awk '{sum+=($3-$2)} END{print sum}' $annot)
genome_len=$(awk '{sum+=$2} END{print sum}' $faidx)



RANDOM=$(date +%N | cut -b4-9)


HDALvHDLI_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dir/HDALvHDLI/dmrs.bed)
HDALvLDAL_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dir/HDALvLDAL/dmrs.bed)

HDALvHDLI_overlap=$(bedtools intersect -header -wo -a $dir/HDALvHDLI/dmrs.bed -b $annot | cut -f20 | awk '{sum+=$1}END{print sum}')
HDALvLDAL_overlap=$(bedtools intersect -header -wo -a $dir/HDALvLDAL/dmrs.bed -b $annot | cut -f20 | awk '{sum+=$1}END{print sum}')

resamples=1000

for i in $(seq 1 $resamples)
do


bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/HDALvHDLI/dmrs.bed) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> HDALvHDLI_resample_overlap &

bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/HDALvLDAL/dmrs.bed) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> HDALvLDAL_resample_overlap &
	wait
done

HDALvHDLI_ln=$(cat HDALvHDLI_resample_overlap <(echo $HDALvHDLI_overlap) | sort -n | grep -n -w "$HDALvHDLI_overlap" | cut -f1 -d ":" | head -n1 )
HDALvLDAL_ln=$(cat HDALvLDAL_resample_overlap <(echo $HDALvLDAL_overlap) | sort -n | grep -n -w "$HDALvLDAL_overlap" | cut -f1 -d ":" | head -n1 )


HDALvHDLI_prop=$(awk -v HDALvHDLI_ln=$HDALvHDLI_ln -v resamples=$resamples 'BEGIN{r=resamples+1-HDALvHDLI_ln; print r/(resamples)}')
HDALvLDAL_prop=$(awk -v HDALvLDAL_ln=$HDALvLDAL_ln -v resamples=$resamples 'BEGIN{r=resamples+1-HDALvLDAL_ln; print r/(resamples)}')


# DMR total length, DMR annotation overlap, annotation as a fraction of the genome, odds ratio of DMRs overlapping annotation, empirical proportion of samples with greater overlap than observed. Can be transformed to two-tailed p-values through the following formula: if prop <= 0.5 (i.e. lower tail) prop*2, if prob > 0.5 then take (1-prop)*2
awk -v annot_len=$annot_len -v genome_len=$genome_len -v HDALvHDLI_DMR_len=$HDALvHDLI_DMR_len -v HDALvHDLI_overlap=$HDALvHDLI_overlap -v HDALvHDLI_prop=$HDALvHDLI_prop -v shortAnnot=$shortAnnot 'BEGIN{print HDALvHDLI_DMR_len "\t" HDALvHDLI_overlap  "\t" annot_len/genome_len "\t" (HDALvHDLI_overlap/HDALvHDLI_DMR_len)/(annot_len/genome_len) "\t" HDALvHDLI_prop "\t" "HDALvHDLI" "\t" shortAnnot}' >> annotation_overlap/stats_HDALvHDLI
awk -v annot_len=$annot_len -v genome_len=$genome_len -v HDALvLDAL_DMR_len=$HDALvLDAL_DMR_len -v HDALvLDAL_overlap=$HDALvLDAL_overlap -v HDALvLDAL_prop=$HDALvLDAL_prop -v shortAnnot=$shortAnnot 'BEGIN{print HDALvLDAL_DMR_len "\t" HDALvLDAL_overlap  "\t" annot_len/genome_len "\t" (HDALvLDAL_overlap/HDALvLDAL_DMR_len)/(annot_len/genome_len) "\t" HDALvLDAL_prop "\t" "HDALvLDAL" "\t" shortAnnot}' >> annotation_overlap/stats_HDALvLDAL


rm HDALvHDLI_resample_overlap
rm HDALvLDAL_resample_overlap


rm tmp_annot_sp

done < "annotation.list.csv" #Annotation list in comma-separated values format
