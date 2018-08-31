#!/bin/sh
AF_cutoff=0.0001
input="Data/simulation_db/gnomad_PASS_sort_subpop_snpEff_ns_All_CADD_all"
output=${input}.screen_$AF_cutoff
 grep -v -P "gene_fusion|\s+upstream|\s+downstream|\s+3_prime_UTR_variant|\s+5_prime_UTR_variant|\s+intron|\s+intergenic|\s+non_coding|\s+sequence_feature|\s+synonymous_variant"  $input | awk '{split($1,a,":"); split($12,z,":");split($13,x,":");split($14,y,":");split($15,k,":");split($16,b,":");split($17,l,":");split($22,n,":");split($23,m,":");split($24,o,":");split($25,c,":"); split($26,d,":");split($27,e,":");split($28,f,":");split($29,g,":"); split($37,h,":");  h1= 1-(10^(-(h[2]/10))); if(l[2]<=$AF_cutoff){print a[1]"\t"a[2]"\t.\t"a[3]"\t"a[4]"\tmax:"k[2]":"b[2]":"l[2]"\tAdj:"x[2]":"z[2]":"y[2]"\tAFR:NA:NA:"n[2]"\tAMR:NA:NA:"m[2]"\tASJ:NA:NA:"o[2]"\tEAS:NA:NA:"c[2]"\tFIN:NA:NA:"d[2]"\tNFE:NA:NA:"e[2]"\tSAS:NA:NA:"f[2]"\tOTH:NA:NA:"g[2]"\t"$36"\t"$37"\tCADD_transf:"h1"\t"$6"\t"$2}}'   > $output
