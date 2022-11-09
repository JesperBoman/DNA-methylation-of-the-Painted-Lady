#!/bin/bash -l

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to annotate DMRs with gene ontology terms from an annotation
# ==========================================================================================================
# Jesper Boman                      9 nov 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #



module load bioinfo-tools BEDTools/2.29.2

dir="VanessaDNAm/DMRs/results"
mRNA_list="VanessaDNAm/annotation/mRNA.gff"
prom_list="VanessaDNAm/annotation/promoters.list"
#HDALvHDLI
#HDALvLDAL


comp="HDALvHDLI"

while IFS= read -r annot_double
do

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annot=$(cut -f2 <(echo $annot_double) -d ",")

#Bedtools merge here is a convenient heuristic since annotation features can overlap but are bound to a specific gene such as is the case for UTRs
#Genes can also overlap and it is probably fine to consider them as one long gene region instead of counting them as a double overlap.
awk '{print $1 "\t" $2 "\t" $3}' $annot | sort -k1,1 -k2,2n | bedtools merge > tmp_annot_sp
annot=tmp_annot_sp


if [[ "$shortAnnot" !=  "Promoters" ]]
then

bedtools intersect -a $dir/$comp/dmrs.bed -b $annot -wa > tmp_dmrs.bed
bedtools intersect -a tmp_dmrs.bed -b $mRNA_list -wb > tmp_dmrs_mrna.bed

cut -f25 tmp_dmrs_mrna.bed |  sed -E 's/ID.*;Parent=(.*);Name.*;/\1/g' > tmpIDs
cut -f16 tmp_dmrs_mrna.bed > tmpMethDirection

awk 'NR == FNR{geneIDs[NR]=$0 } NR != FNR{print geneIDs[FNR] "\t" $0}' tmpIDs tmpMethDirection> $comp.$shortAnnot.geneIDs.list


else

#Promoter route

bedtools intersect -a $dir/$comp/dmrs.bed -b $annot -wb | cut -f16- > tmp_dmrs.bed
awk 'NR==FNR{a[$1 "\t" $2-1 "\t" $3]=$5 } NR!=FNR{if( ($2 "\t" $3 "\t" $4) in a ){print a[$2 "\t" $3 "\t" $4] "\t" $1}}' $prom_list tmp_dmrs.bed > $comp.$shortAnnot.geneIDs.list

fi

grep "hyper" $comp.$shortAnnot.geneIDs.list > GO_analysis/$comp.$shortAnnot.geneIDs.hyper.list
grep "hypo" $comp.$shortAnnot.geneIDs.list > GO_analysis/$comp.$shortAnnot.geneIDs.hypo.list

mv $comp.$shortAnnot.geneIDs.list GO_analysis


done < "annotation.list.geneFeatures.csv"

rm tmp_dmrs.bed
rm tmp_dmrs_mrna.bed
rm tmpIDs
rm tmpMethDirection
rm tmp_annot_sp
