#!/usr/bin/perl -w
use Cwd;
use strict;
use Getopt::Long;
use List::Util qw( min max );

use vars qw( $opt_help $opt_input $opt_output $opt_cutoff $opt_num  $opt_eth $opt_ratio);

if (!GetOptions(
                "help|h",
                "input|i:s",
                "output|o:s",
                "cutoff|c:s",
		"num|n:s",
		"eth|e:s",
		"ratio|r:s",
            )) {
    &usage(1);
}


if ($opt_help) {
    &usage(0);
}

if (($opt_input eq "")||($opt_output eq "")||($opt_num eq "")) {&usage(0); exit;}
if($opt_cutoff eq ""){
$opt_cutoff = 0.005;
}
if($opt_ratio eq ""){
$opt_ratio=0.5;
}
if(($opt_eth eq "") || (!(($opt_eth eq "max")||($opt_eth eq "AFR")||($opt_eth eq "AMR")||($opt_eth eq "Adj")||($opt_eth eq "EAS")||($opt_eth eq "FIN")||($opt_eth eq "NFE")||($opt_eth eq "OTH")||($opt_eth eq "SAS")||($opt_eth eq "raw")||($opt_eth eq "min")))){
  print "the ethnicity option is invalidate\n";
  &usage(0);

}
#for (my $i=0; $i< $opt_num; $i++){ 
my $female_num = int($opt_num*$opt_ratio);
print "total_num:$opt_num\tfemale_num:$female_num\n";
my %hash = ();

open(INPUT,$opt_input);
my %hash_sum;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
#my $chr = $info[0];
#my $pos = $info[1];
#my $ref = $info[2];
my $id = "$info[0]"."_$info[1]";
my $id1 = "$info[3]"."_$info[4]";	
my $allele_sum=0;
#for(my $j=4; $j<=$#info-4; $j++){
for(my $j=4; $j<=$#info; $j++){

        #my @info1 = split(/\;/,$info[$j]);
	#	for( my $k=1; $k<=$#info1; $k++){
	#	chr1    13372   G       C       max:2:5052:0.000395882818685669042      AFR:0:770:0     AMR:0:134:0     Adj:2:8432:0.000237191650853889943      EAS:0:254:0     FIN:0:16:0      NFE:0:2116:0    OTH:0:90:0      SAS:2:5052:0.000395882818685669042
			my @info2=split(/\:/,$info[$j]);
			if($info2[0] eq "$opt_eth"){
			#	if(($info2[3] <= $opt_cutoff)&&($info2[3]>0)){
			        
				if((defined $info2[3]) && ($info2[3] !~ /^\.$/) && ($info2[3]>0)){
#			   	$hash{$id}->{$info1[0]}->{indi_freq}= $info2[7];
#				$hash{$id}->{$info1[0]}->{het_hem} = $info2[8]+$info2[10];
#	       	 		$hash{$id}->{$info1[0]}->{hom} = $info2[9];
				$hash{$id}->{$id1}->{allele_freq} = $info2[3];
				#$hash_sum{$id} += $info2[3];
#				my @raw = split(/:/,$info[15]);
				my @raw = split(/:/,$info[-5]);

				$hash{$id}->{$id1}->{raw_score} =$raw[1];
#				my @phred = split(/:/,$info[16]);
 				my @phred = split(/:/,$info[-4]);

				$hash{$id}->{$id1}->{phred_score} =$phred[1];
#				my @transf = split(/:/,$info[17]);

				my @transf = split(/:/,$info[-3]);
				$hash{$id}->{$id1}->{transf_score} = $transf[1];
#				$hash{$id}->{$id1}->{variant} = $info[18];
#				$hash{$id}->{$id1}->{gene} = $info[19];

				$hash{$id}->{$id1}->{variant} = $info[-2];
				$hash{$id}->{$id1}->{gene} = $info[-1];


				}
			last;
			}
		#}
}
#if($allele_sum >0){
#	my $temp_sum =0;
#}
}

#for my $key1 (keys %hash){
#	for my $key (keys %{$hash{$key1}}){
#	my $temp_sum = $hash{$key1}->{$key}->{allele_freq};
#	$hash{$key1}->{$key}->{nuc_freq}=$temp_sum/$hash_sum{$key1};
#	}
#}
my $count=0;
#`mkdir $opt_output`;
for (my $n =0 ; $n< $opt_num ; $n++){
#my $output = "$opt_output/Exac"."_$n";
my $output = "$opt_output"."_Exac"."_$n";

open(OUTPUT,">$output");
print OUTPUT "#CHROM\tPOS\tREF\tALT\tGENE\tSCORE\n";
for my $key (keys %hash ){
	my $l=0;
	my @arr;
	my @id = split(/\_/,$key);
	
	my $alt="NA";
	my $sum=0;
	for my $key1 (keys %{$hash{$key}} ) {
		$sum += $hash{$key}->{$key1}->{allele_freq};
		$arr[$l]->{nuc} = $key1;
		$arr[$l]->{nuc_freq} = $sum; 
# 		$arr[$l]->{raw_score} = $hash{$key}->{$key1}->{raw_score};
#		$arr[$l]->{phred_score} = $hash{$key}->{$key1}->{phred_score};
#		$arr[$l]->{trans_score} = $hash{$key}->{$key1}->{transf_score};
#		$arr[$l]->{variant} = $hash{$key}->{$key1}->{variant};
#		$arr[$l]->{gene} = $hash{$key}->{$key1}->{gene};
		$l++;
	}	

	$count++;

##
#the first allele
##
	#$id[0] =~ s/chr//g;

	if(($id[0]=~/Y/)&&($count<=$female_num)){
	;
	}else{
	my $ran_nuc = rand();
	if($ran_nuc < $arr[0]->{nuc_freq}){
	   $alt = $arr[0]->{nuc};
	}
	for(my $m=1; $m<=$#arr; $m++){
	if(($ran_nuc < $arr[$m]->{nuc_freq})&& ($ran_nuc >= $arr[$m-1]->{nuc_freq})){	
		$alt = $arr[$m]->{nuc};
		last;
	}
	}
	if($alt ne "NA"){
	if($hash{$key}->{$alt}->{allele_freq} <= $opt_cutoff){
	my @id1 = split(/\_/,$alt);
	print OUTPUT "$id[0]\t$id[1]\t$id1[0]\t$id1[1]\t$hash{$key}->{$alt}->{gene}\t$hash{$key}->{$alt}->{transf_score}\t$hash{$key}->{$alt}->{raw_score}\t$hash{$key}->{$alt}->{phred_score}\n";
#	print OUTPUT "$id[0]\t$id[1]\t$id1[0]\t$id1[1]\t$hash{$key}->{$alt}->{gene}\t$hash{$key}->{$alt}->{raw_score}\t$hash{$key}->{$alt}->{phred_score}\t$hash{$key}->{$alt}->{transf_score}\n";

	}
	}
	}
####
#the second allele
#####
	if((($count>$female_num)&&($id[0] =~ /X/))||($id[0] =~ /Y/)){
	;
	}else{
	$alt = "NA";
	my $ran_nuc1 = rand();
	if($ran_nuc1 < $arr[0]->{nuc_freq}){
	   $alt = $arr[0]->{nuc};
	}
	for(my $m=1; $m<=$#arr; $m++){
	if(($ran_nuc1 < $arr[$m]->{nuc_freq})&& ($ran_nuc1 >= $arr[$m-1]->{nuc_freq})){	
		$alt = $arr[$m]->{nuc};
		last;
	}
	}
	if($alt ne "NA"){
	if($hash{$key}->{$alt}->{allele_freq} <= $opt_cutoff){
	my @id1 = split(/\_/,$alt);
	print OUTPUT "$id[0]\t$id[1]\t$id1[0]\t$id1[1]\t$hash{$key}->{$alt}->{gene}\t$hash{$key}->{$alt}->{transf_score}\t$hash{$key}->{$alt}->{raw_score}\t$hash{$key}->{$alt}->{phred_score}\n";

#	print OUTPUT "$id[0]\t$id[1]\t$id1[0]\t$id1[1]\t$hash{$key}->{$alt}->{gene}\t$hash{$key}->{$alt}->{raw_score}\t$hash{$key}->{$alt}->{phred_score}\t$hash{$key}->{$alt}->{transf_score}\n";


#	print OUTPUT "$id[0]\t$id[1]\t.\t$id1[0]\t$id1[1]\t$hash{$key}->{$alt}->{raw_score}\t$hash{$key}->{$alt}->{phred_score}\n";
	}
	}
	}

}
close(OUTPUT);
}




sub usage {


                print "
       vcf simulation tool
       Authors: Jun Wang;  Rui Chen Lab
       Date: 8-30-2018

       Usage: -i: input file, vcf4 format
              -o: output file, vcf4 format 
              -c: allele frequency cutoff, default is 0.005
              -n: integer number of simulation
	      -e: ethnicity, need to one of the parameters below
			     max: the population with maxium allele frequency
			     AFR: African/African American
			     AMR: American 
			     Adj: Adjusted all 
			     EAS: East Asian
			     FIN: Finnish
			     NFE: Non-Finnish European
			     OTH: Other
			     SAS: South Asian
			     raw: raw count	
	      -r: female vs. male ratio, the number of simulated females will be
	          the integer number of (total_num * female_male_ratio)		     
	\n";
                exit;

}

