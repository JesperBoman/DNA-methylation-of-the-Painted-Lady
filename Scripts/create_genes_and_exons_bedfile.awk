#!/usr/bin/awk -f
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# awk -f create_genes_and_exons_bedfile.awk <(  sort -k1,1 -k4,4n -k5,5n ../annotation/vcard.gtf) | sort -k1,1 -k2,2n -k3,3n > genes_and_exons.bed
# ==========================================================================================================
# Jesper Boman                      9 Nov 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #


NR>1{
if($3 == "gene"){
	split($10, a, "\""); 
	gene[a[2]] = $1 "\t" $4-1 "\t" $5 "\t" a[2] "\t" "gene" "\t" "NA" "\t" $7 
} 
if($3 == "exon"){split($10, b, "\"");  
	exon_count[b[2]]++; 
	if($7 == "+"){plus_exon[b[2], exon_count[ b[2]]] = $1 "\t"$4-1 "\t" $5 "\t" b[2] "\t" "exon" "\t" "exon_" exon_count[b[2]] "\t" $7 }
	 else{minus_exon[b[2] "|" exon_count[b[2]]] = $1 "\t" $4-1 "\t" $5 "\t" b[2] "\t" "exon" "\t" "exon_" }}  
}

END{ 
for(i in minus_exon){
	split(i, c, "|"); print minus_exon[i] exon_count[c[1]] - c[2] + 1 "\t" "-"  } 
for(i in plus_exon){
	print plus_exon[i]}
for(i in gene){print gene[i]}
}
