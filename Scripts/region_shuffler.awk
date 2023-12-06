#!/usr/bin/awk -f

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to randomly place regions in the genome (provided in a .bed - file) of fixed lengths
# ==========================================================================================================
# Jesper Boman                      6 dec 2023
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Usage: awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx file.bed
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #

# -v genome_len=$genome_len  #Total length of the genome
# -v RANDOM = $RANDOM #To set seed 

## NOTES ##
#Chromosome positions are sampled with replacement - this might conflict with your null model

#One modification which could be reasonable is to not sample chromosome positions with replacement
#But at the same time if two randomly placed regions overlap eachother 
#and an e.g. annotation feature .bed file then both will be counted in the intersect so it might not be a problem



#First file faidx of genome assembly fasta

BEGIN{srand(RANDOM); prev=0}

FNR==NR{
chrom_cumprob[$1]=prev+($2/genome_len);
chrom_len[$1]=$2;
sequence[n++] = $1;
prev=prev+($2/genome_len)
}

#Second file should be a bed file with a header, change it if doesn't have a header

FNR!=NR && FNR > 1{

DMR_len=$3-$2;

#Pick a chromosome randomly
randchrom_num=rand();

for (i = 0; i < n; i++){
if(randchrom_num < chrom_cumprob[sequence[i]]){
chrpick=sequence[i];
rand_end=int(rand()*chrom_len[chrpick]);
rand_start=rand_end-DMR_len;

#This is an approximation that works well if you are shuffling regions much smaller than the average chromosome size.
#What could happen here is that the end overshoots, since we do not double check with annotation. Bedtools doesn't seem to mind.
if(rand_start<0){
rand_start=rand_end;
rand_end=rand_end+DMR_len}

break 
}



}
print chrpick "\t" rand_start "\t" rand_end


}
