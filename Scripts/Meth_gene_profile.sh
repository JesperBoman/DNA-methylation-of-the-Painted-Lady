#!/bin/bash -l

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to make a dataset of gene profile from an annotation and associate external data, such as methylation level
# ==========================================================================================================
# Jesper Boman                      6 oct 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Create a .bed-file with genes and exons
awk -f create_genes_and_exons_bedfile.awk <(  sort -k1,1 -k4,4n -k5,5n ../annotation/vcard.gtf) | sort -k1,1 -k2,2n -k3,3n > genes_and_exons.bed


#Filter for genes that have both 5' and 3' UTRs 


awk '{if($3 == "gene" && $7 == "+"){split($10, a, "\""); plus_gene[a[2]]= $1 "\t" $4-1 "\t" $5 "\t" $7 "\t" a[2]};
  if($3 == "gene" && $7 == "-"){split($10, a, "\""); minus_gene[a[2]]= $1 "\t" $4-1 "\t" $5 "\t" $7 "\t" a[2]; minus_gene_end[$1,$4]};
  if($3 == "five_prime_utr" && $7 == "+"){split($10, a, "\""); plus_5utr[a[2]]};
  if($3 == "five_prime_utr" && $7 == "-"){split($10, a, "\""); minus_5utr[a[2]]};
  if($3 == "three_prime_utr" && $7 == "+"){split($10, a, "\""); plus_3utr[a[2]]};
  if($3 == "three_prime_utr" && $7 == "-"){split($10, a, "\""); minus_3utr[a[2]]};
  }
  END{for (gene in plus_gene){if( (gene in plus_5utr) && (gene in plus_3utr) ){print plus_gene[gene]}};
  for (gene in minus_gene){if( (gene in minus_5utr) && (gene in minus_3utr) ){print minus_gene[gene]}}}' ../annotation/vcard.gtf > genes_with_both_five_and_three_prime_UTRs_IDs.bed 
  


  
awk 'FNR==NR{a[$5]} FNR!=NR{if($4 in a)print $0 }' genes_with_both_five_and_three_prime_UTRs_IDs.bed  genes_and_exons.bed >  genes_and_exons_5n3p.bed
  
  
#Gene profile - split up intervals in 100 bp segment 

#Chop up into 100 bp pieces
ml bioinfo-tools BEDOPS/2.4.39

#Obtain 5kb flanks (obs: does not check if it is beyond the upper limit of the chromosome)

awk '{if($5 == "gene"){if(($2-5000)<0){start=0}else{start=$2-5000}; print $1 "\t" start "\t" $3+5000 "\t" $4 "\t" $5 "\t" $6 "\t" $7}}' genes_and_exons_5n3p.bed > genes_5n3p_5kb_flanks.bed



while IFS=$'\t' read -r gene_entry
do
	gene=$(echo "$gene_entry" | cut -f4)
	strand=$(echo "$gene_entry" | cut -f7)
	bedops --chop 100 <(echo "$gene_entry") > temp_gene_split
	awk -v gene="$gene" -v strand="$strand" '{print $0 "\t" gene "\t" strand "\t" "seg_" NR}' temp_gene_split >> genes_5n3p_5kb_flanks_100bp_seg.bed
done < "genes_5n3p_5kb_flanks.bed" 


#Check
cat genes_5n3p_5kb_flanks_100bp_seg.bed | cut -f4 | sort | uniq | wc -l
#3685



###INTERSECTION WITH WGBS DATA


#create a sample_list

module load bioinfo-tools BEDTools/2.29.2

bed_wgbs_data_dir=""

while IFS= read -r sample
do
        ind=$(echo $sample | cut -d- -f1)
        tissue=$(echo $sample | cut -d. -f1 | cut -d- -f2)

#Intersection of gene features
        bedtools intersect -a <(zcat $bed_wgbs_data_dir/$sample) \
        -b genes_5n3p_5kb_flanks_100bp_seg.bed -wa -wb > temp_prof_sample.bed
        

#Output: Chromosome Start End Gene_name Strand Total_CpG Post_filter_CpG Methylated_reads_per_feature Unmethylated_reads_per_feature  Sum_of_proportions_of_methylated_over_total_reads_per_CpG Mean_per_dinuc (i.e. Methylation level) Individual Tissue

awk -v ind="$ind" -v tissue="$tissue" 'BEGIN {OFS="\t"} {if(!(annot in pos)){skip[annot]=="FALSE"}; annot=$9 "\t" $10 "\t" $11 "\t" $12; 
if(pos[annot]==$2 && $7=="CG" && prev=="CG" && NR!=1 && skip[annot]=="FALSE"){ 
        if(($5+$6+prev_total_reads[annot])>5 && ($5+$6+prev_total_reads[annot])<200){skip[annot]="TRUE";
        a[annot]+=1;  b[annot]+=($5+prev_meth_reads[annot]); c[annot]+=($6+prev_unmeth_reads[annot]); rat[annot]+=(($5+prev_meth_reads[annot])/($5+$6+prev_total_reads[annot])) }; 
        ;  strand[annot]=$13; seg_num[annot]=$14}
    else{skip[annot]="FALSE"}; 
prev=$7; if(skip[annot]=="FALSE"){pos[annot]=$3} ; 
if($7=="CG"){
        totCpG[annot]+=0.5; prev_meth_reads[annot]=$5; prev_unmeth_reads[annot]=$6; prev_total_reads[annot]=$5+$6; strand[annot]=$13; seg_num[annot]=$14}
}
    
    END{for (wind in totCpG){if (a[wind]==""){print wind, seg_num[wind], strand[wind], totCpG[wind], "NA", "NA", "NA", "NA", "NA", ind, tissue}
    else{print wind, seg_num[wind], strand[wind], totCpG[wind], a[wind], b[wind], c[wind], rat[wind], rat[wind]/a[wind], ind, tissue}}}' temp_prof_sample.bed >> meth_gene_5kb_flanks_100bp_seg_data
#Cleaning up
        rm temp_prof_sample.bed
done < "sample_list"








#PROFILE COORDINATES

#Script for positions along the gene and ordering of upstream and downstream segments


awk '{if($5 == "gene")print}'  genes_and_exons_5n3p.bed >  genes_5n3p.bed


#Output: Chromosome Start End Gene_name Segment_# Strand Total_CpG Post_filter_CpG Methylated_reads_per_feature Unmethylated_reads_per_feature  Sum_of_proportions_of_methylated_over_total_reads_per_CpG Mean_per_dinuc (i.e. Methylation level) Individual Tissue Segment_type Segment_number (Coordinate from 0 to 1 if a genic element)


awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)}; FNR==NR{gl[$4]=$3-$2} FNR!=NR{split($5,seg_num,"_"); if(seg_num[2]<=50){if($6 == "+"){print $0 "\t" "Upstream" "\t" seg_num[2] "\t" "NA"} else{print $0 "\t" "Downstream" "\t" 51-seg_num[2] "\t" "NA"}}
         else if(50 < seg_num[2] && (seg_num[2] < (ceil(gl[$4]/100)+51))){if($6 == "+"){print $0 "\t" "Genic" "\t" seg_num[2]-50 "\t" ((seg_num[2]-50.5)*100)/gl[$4]} else{print $0 "\t" "Genic" "\t" seg_num[2]-50 "\t" 1-(((seg_num[2]-50.5)*100)/gl[$4])} }
         else if((ceil(gl[$4]/100)+50) < seg_num[2]){if($6 == "+"){print $0 "\t" "Downstream" "\t" seg_num[2]-(ceil(gl[$4]/100)+50) "\t" "NA" } else{print $0 "\t" "Upstream" "\t" 50-(seg_num[2]-(ceil(gl[$4]/100)+51)) "\t" "NA" }}}' genes_5n3p.bed meth_gene_5kb_flanks_100bp_seg_data > meth_gene_5kb_flanks_100bp_seg_data_coords




### INTERSECTION WITH EXPR (in this case: TPM) data

expr_data="salmon.merged.gene_tpm.tsv"


awk ' FNR==NR{if(NR!=1){HDALX1[$1]=$3; HDALX18[$1]=$4; HDALX3[$1]=$5; HDLIX18[$1]=$6; HDLIX3[$1]=$7; LDALX1[$1]=$8; 
	LDALX18[$1]=$9; LDALX3[$1]=$10;}} FNR!=NR{ 
	if($13 == "HDAL-X1"){print $0 "\t" HDALX1[$4]}
	if($13 == "HDAL-X18"){print $0 "\t" HDALX18[$4]}
	if($13 == "HDAL-X3"){print $0 "\t" HDALX3[$4]}
	if($13 == "HDLI-X18"){print $0 "\t" HDLIX18[$4]}
	if($13 == "HDLI-X3"){print $0 "\t" HDLIX3[$4]}
	if($13 == "LDAL-X1"){print $0 "\t" LDALX1[$4]}
	if($13 == "LDAL-X18"){print $0 "\t" LDALX18[$4]}
	if($13 == "LDAL-X3"){print $0 "\t" LDALX3[$4]}
	if($13 == "HDLI-X1"){print $0 "\t" "NA"}}' $expr_data  meth_gene_5kb_flanks_100bp_seg_data_coords > meth_gene_5kb_flanks_100bp_seg_data_coords_expr







# Exon/Intron annotation #


module load bioinfo-tools BEDTools/2.29.2

awk '{if($5 == "exon")print $0}' genes_and_exons_5n3p.bed > exons_5n3p.bed

bedtools subtract -a genes_5n3p.bed -b exons_5n3p.bed > introns_raw_5n3p.bed


sort -k1,1 -k2,2n -k3,3n introns_raw_5n3p.bed > introns_raw_5n3p_sort.bed

awk 'NR==FNR{gene_strand[$4]=$7; intron_num[$4]++} NR!=FNR{if ($7 == "+"){ a[$4]+=1; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "intron" "\t" "intron_"a[$4] "\t" $7}
	else{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "intron" "\t" "intron_"intron_num[$4]-a[$4] "\t" $7; a[$4]+=1}}' introns_raw_5n3p_sort.bed introns_raw_5n3p_sort.bed > intron_intronNumbers_5n3p.bed



segment="genes_5n3p_5kb_flanks_100bp_seg.bed"

bedtools intersect -wo -f 0.9 -a $segment -b intron_intronNumbers_5n3p.bed  > genes_5n3p_5kb_flanks_100bp_seg_intronOverlap.bed
bedtools intersect -wo -f 0.9 -a $segment -b exons_5n3p.bed > genes_5n3p_5kb_flanks_100bp_seg_exonOverlap.bed


meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_expr"

awk 'FNR==NR{if($11 == "intron"){intSeg[$1 "\t" $2 "\t" $3]= $12} else{exSeg[$1 "\t" $2 "\t" $3]= $12}}
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in intSeg){print $0 "\t" intSeg[seg]} else if(seg in exSeg){print $0 "\t" exSeg[seg]} else{print $0 "\t" "NA"}}' <(cat genes_5n3p_5kb_flanks_100bp_seg_intronOverlap.bed genes_5n3p_5kb_flanks_100bp_seg_exonOverlap.bed)  $meth_table > ${meth_table}_ExIn





# CpG density and GC content annotation #

module load bioinfo-tools BEDTools/2.29.2

genome="/crex/proj/uppstore2017185/b2014034_nobackup/Jesper/VanessaDNAm/reference/GCA_905220365.1_ilVanCard2.1_genomic_chroms.fasta"

bedtools nuc -bed genes_5n3p_5kb_flanks_100bp_seg.bed -fi $genome -pattern CG | awk 'NR>1{totbases=$9+$10+$11+$12; if(totbases != 0 && ($10 < 95) && ($11 < 95) ){print $1 "\t" $2 "\t" $3 "\t"  ($10+$11)/totbases "\t" ($10*$11)/(totbases) "\t" $16} else{print $1 "\t" $2 "\t" $3 "\t" "NA" "\t" "NA" "\t" "NA"}}' > genes_5n3p_5kb_flanks_100bp_seg_GC_content_CpGOE.bed

meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_expr"

awk 'FNR==NR{a[$1 "\t" $2 "\t" $3]=$4 "\t" $5 "\t" $6}
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in a){print $0 "\t" a[seg]} else{print $0 "\t" "NA" "\t" "NA" "\t" "NA"}}' genes_5n3p_5kb_flanks_100bp_seg_GC_content_CpGOE.bed ${meth_table}_ExIn > ${meth_table}_ExIn_CpG


