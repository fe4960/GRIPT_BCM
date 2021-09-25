#!/usr/bin/perl -w
use Cwd;
use strict;
use Getopt::Long;
use List::Util qw( min max );

use vars qw( $opt_help $opt_remove_snp $opt_snpeff $opt_input $opt_output $opt_add_CADD $opt_miss_CADD $opt_tabix_output $opt_filter_gnomad $opt_freq_cutoff $opt_read_flag $opt_CADD_flag $opt_snpeff_jar $opt_snpeff_conf $opt_readtype $opt_CADD_snp_db $opt_CADD_indel_db $opt_gnomad_db $opt_R_dir $opt_tabix);

if (!GetOptions(
                "help|h",
                "remove_snp|rml:s",
                "snpeff|sf:s",
		"input|i:s",
		"output|o:s",
                "add_CADD|ad_cad:s",
                "miss_CADD|mis_cad:s",
                "tabix_output|tabix_out:s",
                "filter_gnomad|f_gn:s",
                "freq_cutoff|fc:s",
                "filter_read|f_r:s",
                "CADD_flag|cf:s",
                "snpeff_jar|sfj:s",
                "snpeff_conf|sfc:s",
                "readtype|rt:s",
                "CADD_snp_db|cd_snp:s",
                "CADD_indel_db|cd_ind:s",
                "gnomad_db|gn:s",
                "R_dir|r:s",
                "tabix|tab:s"
            )) {
    &usage(1);
}


if ($opt_help) {
    &usage(0);
}

sub usage {


                print "
       GRIPT preprocessing script
       Authors: Jun Wang;  Rui Chen Lab
       Date: 8-30-2018

       Usage: -i:         input folder of vcf files, vcf4 format
              -o:         output folder of files for GRIPT analysis
              -ad_cad:    a CADD score file generated from CADD website by uploading the variant vcf file to CADD website
              -rml:       a list of variants to be removed
              -mis_cad:   a output list of variants without CADD score. default = GRIPT_CADD_missing
              -tabix_out: a temp file for tabix output. default = output folder namne + _tabix_tmp
              -f_gn:      a binary value indicating whether applying gnomad allele freq filtering. default = 1
              -fc:        a variant freq filtering cutoff. default = 0.005
              -f_r:       a binary value indicating whether applying read coverage filtering for variant calling in vcf file. default =1
              -cf:        a binary value indicating whether annotating CADD score as the variant score. default = 1
              -sf:        a binary value indicating whether using snpEff to annotate input variants. default = 1
              -sfj:       the directory of snpeff_jar. default = snpEff/snpEff.jar
              -sfc:       the directory of snpeff_conf. default = snpEff/snpEff.config
              -rt:        the variant caller of the input vcf, two option: GATK or Atlas, default = GATK
              -cd_snp:    the directory of snp CADD score gz file downloaded from CADD website. default = CADD_score/whole_genome_SNVs.tsv.gz
              -cd_ind:    the directory of indel CADD score gz file downloaded from CADD website. default = CADD_score/InDels.tsv.gz
              -gn:        the directory of variant frequency from gnomAD database. default = gnomad/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz
              -r:         the path of R. default = R
              -tab:       the path of tabix. default = tabix
	\n";
                exit;

}


if (( ! $opt_input )||( ! $opt_output)) {&usage(0); exit;}

if (! $opt_read_flag) {$opt_read_flag=1;}

if(! $opt_freq_cutoff){ $opt_freq_cutoff = 0.005; }

if(! $opt_CADD_flag ){$opt_CADD_flag=1;}

if(! $opt_filter_gnomad){$opt_filter_gnomad=1;}

if(! $opt_CADD_snp_db ){$opt_CADD_snp_db = "CADD_score/whole_genome_SNVs.tsv.gz";}
if(! $opt_CADD_indel_db){$opt_CADD_snp_db = "CADD_score/InDels.tsv.gz";}
if(! $opt_gnomad_db){$opt_gnomad_db = "gnomad/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz";}
if(! $opt_readtype){$opt_readtype = "GATK";}
if(! $opt_snpeff_conf ){$opt_snpeff_conf = "snpEff/snpEff.config";}
if(! $opt_snpeff_jar){$opt_snpeff_jar = "snpEff/snpEff.jar";}
if(! $opt_tabix ){$opt_tabix = "tabix";}
if(! $opt_R_dir){$opt_R_dir = "R";}
if(! $opt_snpeff ){$opt_snpeff=1;}
if(! $opt_tabix_output ){$opt_tabix_output = "$opt_output"."_tabix_tmp";}
if(! $opt_miss_CADD ){$opt_miss_CADD = "GRIPT_CADD_missing";}
##################
#
#
###########
#$opt_tabix = "tabix"; #$ARGV[1];
#$opt_snpeff_jar = "software/snpEff/snpEff.jar"; 
#$opt_snpeff_conf = "software/snpEff/snpEff.config";
#$opt_readtype = "GATK";
#$opt_CADD_snp_db = "/storage/novaseq/backup_jw/CADD_score/whole_genome_SNVs.tsv.gz";
#$opt_CADD_indel_db = "/storage/novaseq/backup_jw/CADD_score/InDels.tsv.gz";
#$opt_gnomad_db = "gnomad/data/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz";
#$opt_R_dir = "R";
#$opt_read_flag=1;
#$opt_freq_cutoff = 0.005;
#$opt_CADD_flag=1;
#$opt_filter_gnomad=1;
#$opt_add_CADD ="";
#$opt_remove_snp = "";
#$opt_tabix_output = "$opt_output"."_tabix_tmp";
my $fisher_tmp = "R_temp".time(); 
#my $filter_db = "";
my $outputlog = "GRIPT_".time().".log";
open(OUTPUTlog,">$outputlog");
############################################
open(OUTPUTerr,">$opt_miss_CADD");
my $cur_time=`date`;
print OUTPUTlog "start at: $cur_time\n";
print OUTPUTlog "input dir: $opt_input\n";
print OUTPUTlog "filter with gnomad: $opt_gnomad_db\tfreq cutoff: $opt_freq_cutoff\n";
print OUTPUTlog "filter with read coverage: $opt_read_flag\n";
print OUTPUTlog "annotate CADD score: $opt_CADD_flag\n";
print OUTPUTlog "CADD_snp_db: $opt_CADD_snp_db\n";
print OUTPUTlog "CADD_indel_db: $opt_CADD_indel_db\n";

my %select_ID;
my %remove;
my %linkage;
my %gnomad_freq;
my %CADD_score;
my %add_CADD;
my %hash_CADD;
############read the add CADD score#######
if($opt_add_CADD){
open(INPUTaddc,$opt_add_CADD);
#CHROM  POS     REF     ALT     RawScore        PHRED
#1       14946   G       A       0.101439        3.623

while(my $line = <INPUTaddc>){
if($line =~ /^#/){
next;
}
chomp $line;
my @info = split(/\s+/,$line);
my $var = $info[0]."_".$info[1]."_".$info[2]."_".$info[3];
#$add_CADD{$var}->{pred} =$info[-1];
#$add_CADD{$var}->{raw} = $info[-2];
my $converted_score = convert_CADD($info[-1]);
$add_CADD{$var}->{total} = "$converted_score\t$info[-2]\t$info[-1]\n";
}
print OUTPUTlog "read added CADD score:$opt_add_CADD\n";
}




##########read the removed snp/indel into %remove###########
if($opt_remove_snp){
open(INPUTrm,"$opt_remove_snp");

while(my $liner = <INPUTrm>){
chomp $liner;
my @infor = split(/\s+/,$liner);
$infor[0] =~ s/chr//g;
my $ref = uc($infor[2]);
my $alt = uc($infor[3]);
my $ID = "$infor[0]"."_"."$infor[1]"."_"."$ref"."_"."$alt";
$remove{$ID}=1;
}
print OUTPUTlog "read removed variant list:$opt_remove_snp\n";
}
###########check whether temp_dir exists. If yes, backup the older temp_dir folder######
my $temp_dir = "$opt_output"."_temp";
if(-d $temp_dir){
my $backup_temp_dir = $temp_dir."_bak.".time();

print OUTPUTlog "$temp_dir exists, is backuped as $backup_temp_dir\n";

`mv $temp_dir $backup_temp_dir`;
}
`mkdir $temp_dir`;
print OUTPUTlog "make temp dir:$temp_dir\n";

###########check whether output_dir exists. If yes, backup the older output_dir folder######
if(-d $opt_output){
my $backup_out_dir = $opt_output."_bak.".time();

print OUTPUTlog "$opt_output exists, is backuped as $backup_out_dir\n";

`mv $opt_output $backup_out_dir`;
}

`mkdir $opt_output`;
print OUTPUTlog "make output dir:$opt_output\n";


###############
opendir(DIR,$opt_input);
my  @dir = grep {-f} map {"$opt_input/$_"} readdir(DIR); ####### read each file in the input_dir #############
print OUTPUTlog "start to process each vcf\n";
for my $in_file (@dir){  ######### process each vcf file in the input_dir############
my $temp_file;
my $out_file;
if($in_file =~ /$opt_input\/(\S+).vcf/){
my $file_name = "$1";   ######### get the file title, which will be used in the following##########

$temp_file = "$temp_dir/$file_name"; 

open(INPUT,"$in_file");
open(OUTPUTt,">$temp_file");

line1: while(my $line = <INPUT>){

if($line =~ /^#/){ ##########skip the annotation header lines in vcf file ############
next;
}
chomp $line;
$line =~ s/chr//g;
$line = uc($line);
my @info = split(/\s+/,$line);
my $variant = "$info[0]"."_"."$info[1]"."_"."$info[3]"."_"."$info[4]";
############if the filter against gnomad flag is on, then call the get_genome_freq sub_rotine############
if( $opt_filter_gnomad ==1){ 
my $max_freq;
if((exists $gnomad_freq{$variant})&&($gnomad_freq{$variant} =~ /\d+/)){
$max_freq = $gnomad_freq{$variant};
}else{
my $region = "$info[0]:$info[1]-$info[1]";  
$max_freq = get_gnomad_freq($opt_tabix,$opt_gnomad_db,$variant,$opt_tabix_output);
$gnomad_freq{$variant} = $max_freq;
}
if(($max_freq ne '.')&&($max_freq > $opt_freq_cutoff)){
goto line1;
}
}

##################if the read coverage filter is on, then filter the variants according to the read coverage######### 
if($opt_read_flag ==1){
my $vr=0;
my $rr=0;
my @read = split(/\:/,$info[9]);

if($opt_readtype eq "Atlas"){  # the read format is: 1/1:34:0:34:.
$vr = $read[1];
$rr = $read[2];
}elsif($opt_readtype eq "GATK"){ # the read format is: 1/1:0,43:43:99:1107,128,0
my @read1 = split(/\,/,$read[1]);
$vr = $read1[1];
$rr = $read1[0];
}

if((($vr+$rr)>=5)&&((($vr == 3) && (($vr/($vr+$rr))>=0.4))||(($vr>=4)&&(($vr/($vr+$rr))>=0.2))||(($vr>=10)&&(($vr/($vr+$rr))>=0.15)))){ #5-24-2017

if(!(exists $remove{$variant})){
if($read[0] =~ /1\/1/){

if(($vr/($vr+$rr))>=0.85){ #5-24-2017

$select_ID{$file_name}->{$variant}=$info[9];
$linkage{$variant}->{$file_name}=1;
}else{
$info[9] =~ s/1\/1/0\/1/g;
$select_ID{$file_name}->{$variant} = $info[9];
$linkage{$variant}->{$file_name}=1;
}
}else{
if(($vr/($vr+$rr))>=0.85){ #5-24-2017
$info[9] =~ s/\d+\/\d+/1\/1/g;
}
$select_ID{$file_name}->{$variant}=$info[9];
$linkage{$variant}->{$file_name}=1;
}
}
}else{
goto line1;
}
}
#################if add CADD score is on,

if($opt_CADD_flag == 1){
if(( defined $hash_CADD{$variant}->{total}) &&($hash_CADD{$variant}->{total} =~ /\d+/)){
my $new_line = "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$hash_CADD{$variant}->{total}"; 
print OUTPUTt "$new_line\n";
	if($info[9] =~ /1\/1/){
	print OUTPUTt "$new_line\n";
	}
}elsif(( defined $add_CADD{$variant}->{total}  ) && ($add_CADD{$variant}->{total} =~ /\d+/)){
my $new_line = "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$add_CADD{$variant}->{total}"; 
print OUTPUTt "$new_line\n";
	if($info[9] =~ /1\/1/){
	print OUTPUTt "$new_line\n";
	}
$hash_CADD{$variant}->{total} = $add_CADD{$variant}->{total};
}else{

my $CADD_score = get_CADD_score($opt_tabix,$variant,$opt_tabix_output);

if($CADD_score =~ /\d+/){
my $new_line = "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$CADD_score"; 
print OUTPUTt "$new_line\n";
	if($info[9] =~ /1\/1/){
	print OUTPUTt "$new_line\n";
	}
$hash_CADD{$variant}->{total} = $CADD_score;
}else{
print OUTPUTerr "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\n";
}
}
}

}


#`/software/jdk1.8.0_45/bin/java -Xmx4g -jar /data/fgi9/jwang/software/snpEff/snpEff.jar -c /data/fgi9/jwang/software/snpEff/snpEff.config -v hg19 $temp_file > $temp_file.snpEff`;
`java -Xmx4g -jar $opt_snpeff_jar -c $opt_snpeff_conf -v hg19 $temp_file > $temp_file.snpEff`;
if($opt_snpeff==1){
#`cat $temp_file.snpEff|tr ';' '\t'|tr '|' '\t'|grep -v "#"|cut -f 1,2,4,5,6,7,8,10,12| awk '{if((\$8 !~/ncRNA/) &&(\$8 !~ /UTR/)&&(\$8 !~ /^intergenic/)){print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$9"\t"\$5"\t"\$6"\t"\$7"\t"\$8}}' > $temp_file.score`;
`cat $temp_file.snpEff|tr ';' '\t'|tr '|' '\t'|grep -v "#"|cut -f 1,2,4,5,6,7,8,10,12| awk '{if(\$8 !~/LOW/){print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$9"\t"\$5"\t"\$6"\t"\$7"\t"\$8}}' > $temp_file.score`;

}else{
`cat $temp_file.snpEff|tr ';' '\t'|tr '|' '\t'|grep -v "#"|cut -f 1,2,4,5,6,7,8,10,12| awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$9"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' > $temp_file.score`;
}



}
}

opendir(DIR,$opt_input);

my  @dir1 = grep {-f} map {"$opt_input/$_"} readdir(DIR);
for my $in_file1 (@dir1){

my $out_file;
if($in_file1 =~ /$opt_input\/(\S+).vcf/){
my $file_name = "$1";
my $temp_file = "$temp_dir/$file_name";
my %remove_link=();

open(INPUTts,"$temp_file.score");
my %gene_id;
while(my $linets = <INPUTts>){
chomp $linets;
my @infos = split(/\s+/,$linets);
my $id = "$infos[0]"."_"."$infos[1]"."_"."$infos[2]"."_"."$infos[3]";
$gene_id{$id}->{id} = $infos[4]; ## the gene id of each variant
$gene_id{$id}->{score} = $infos[7];
}
my $final_output = "$opt_output/$file_name".".score";
my $count=0;
my %file_id=();
my %summary=();
for my $id1 (sort keys %{$select_ID{$file_name}}){ #for each individual, loop through all the variants
if((!(defined $hash_CADD{$id1}->{total})) || ($hash_CADD{$id1}->{total} !~ /\d+/)){
next;
}
   for my $file_name1 (sort keys %{$linkage{$id1}}){ ##for each variant, loop through all the individual
    if(! (exists $summary{$file_name1} )){ #index the inidividual with unique digit 
        $file_id{$count} = $file_name1; 
	$count++;
    }
   $summary{$file_name1}->{$id1}=$gene_id{$id1}->{id}; # for each individual, give the key "id1" the value of genename,
   }
}
my %score=();
	my $filename = $file_name;
       	  for(my $j=0; $j<$count; $j++){ ## for each individual####
		if($file_id{$j} ne $filename){ ########if the indi is not the current indi 
      		my @ids = keys %{$summary{$filename}}; #### get all the variant of the current inidi
		my %link_allele=();
		my %overlap=();
		       for my $ids (@ids){ ####for each variant of the current indi
			  if(defined $summary{$file_id{$j}}->{$ids}){ ###### if the looped inid has the same variant
			     my $locl = $summary{$filename}->{$ids}; #### record the gene of that same variant
			     $link_allele{$locl}++;  #### count how many time  the same gene have hits 
			     $overlap{$locl}->{$ids}=$gene_id{$ids}->{score};	####record the same gene with the overlapped variant the CADD score######
			  }
			}
		for my $subkey (sort keys %link_allele){ ##### for each gene which overlaps with the current indi and looped indi
		if($link_allele{$subkey} >=2){  ##### if the gene have >=2 overlap variants between the current and looped indi

		  my $uniq=0;
		  for my $id_key (sort {$overlap{$subkey}->{$b}<=> $overlap{$subkey}->{$a}} keys %{$overlap{$subkey}}){ ###only keep the variant with the highest CADD score, and remove the rest variants.
		  if($uniq >=1){
		  $remove_link{$filename}->{$id_key}=1;  
		  }
		  $uniq++;
		  }
		}
		}
            }
	}
my @array=();
my $ct=0;
for my $id2 (sort {$a cmp $b} keys %{$select_ID{$file_name}}){ ###loop each variant for cur inidividual 
if((!(defined $hash_CADD{$id2}->{total})) || ($hash_CADD{$id2}->{total} !~ /\d+/)){
next;
}
	    $array[$ct]->{id} = $gene_id{$id2}->{id}; ####give each variant the uniq digital id, and give the genename  
            $array[$ct]->{score} = $gene_id{$id2}->{score}; ####give the CADD score
     	    $array[$ct]->{pos}= $id2; ######give the variant info
	    my @read2 = split(/:/,$select_ID{$file_name}->{$id2}); #####parse the genotype info
	  #  $array[$ct]->{var} = $read[1]; ## if it atlas
	  #  $array[$ct]->{ref} = $read[2];

	    if($opt_readtype eq "Atlas"){  # the read format is: 1/1:34:0:34:.
            $array[$ct]->{var} = $read2[1];
            $array[$ct]->{ref} = $read2[2];
            }elsif($opt_readtype eq "GATK"){ # the read format is: 1/1:0,43:43:99:1107,128,0
            my @read3 = split(/\,/,$read2[1]);
            $array[$ct]->{var} = $read3[1];
            $array[$ct]->{ref} = $read3[0];
            }

	    $ct++;
}
	for(my $i=0; $i<=$#array; $i++){ #### loop each variant of cur individual
	   my @cur_pos = split(/_/,$array[$i]->{pos}); 
	   my $cur_var = $array[$i]->{var};
           my $cur_ref = $array[$i]->{ref};
           my $cur_score = $array[$i]->{score};
           my $cur_id = $array[$i]->{id};
              for(my $j=$i+1; $j<$#array; $j++){ ####loop each variant which is not the cur variant
		my @next_pos = split(/_/,$array[$j]->{pos});
		my $next_var = $array[$j]->{var};
		my $next_ref = $array[$j]->{ref};
		my $next_score = $array[$j]->{score};
                my $next_id = $array[$j]->{id};
                #if(!$next_id){
                #print $array[$j]->{pos};
                #exit;
                #}
		if(($cur_pos[0] eq $next_pos[0])&&(abs($next_pos[1] -$cur_pos[1]) <=100)&&($next_id eq $cur_id)){ ### if the distance between two variants <100 and in the same gene
		my $fisher_p = fisher_test($cur_var,$cur_ref, $next_var,$next_ref,$opt_R_dir,$fisher_tmp); #########conduct a fisher_test
		if($fisher_p >=0.4){ ####### if no difference in the read ratio  
		if($next_score >= $cur_score){ #### compare the two jacent variants, and remove the one with lower CADD score
  		 $remove_link{$file_name}->{$array[$i]->{pos}}=1;
		}else{
		  $remove_link{$file_name}->{$array[$j]->{pos}}=1;
		}
		 }	
		}
              }
         }
open(OUTPUTfinal,">$final_output");
print OUTPUTfinal "#CHROM\tPOS\tREF\tALT\tGENE\tSCORE\tCADD_RAW\tCADD_PHRED\tID\n";
open(INPUTts,"$temp_file.score");

while(my $linets = <INPUTts>){
chomp $linets;
my @infos = split(/\s+/,$linets);
my $id = "$infos[0]"."_"."$infos[1]"."_"."$infos[2]"."_"."$infos[3]";

my $region = "$infos[0]:$infos[1]-$infos[1]";


if((exists $select_ID{$file_name}->{$id})&&($infos[7]>=0)&&(! (exists $remove_link{$file_name}->{$id}))){####if the variant pass previous criteria and was not removed, then print out it in the final output

print OUTPUTfinal "$linets\t$file_name\n";
}
}
}
}
my $end_time=`date`;
print OUTPUTlog "end time: $end_time\n";

####################this subrotine is used to do fisher_test#############
sub fisher_test{
my $values = "$_[0],$_[1],$_[2],$_[3]";
my $R_dir = "$_[4]";
my $R_tmp_in = "$_[5]".".R";
my $R_tmp_out = "$_[5]".".out";
open OUTFILE, ">$R_tmp_in";

print OUTFILE

"TeaTasting<-matrix(c($values),nrow=2,dimnames = list(Guess = c('Milk','Tea'),Truth=c('Milk','Tea')))\n",
"p<-fisher.test(TeaTasting)\$p\n",
"write(p,('$R_tmp_out'))";
close OUTFILE;
`$R_dir < $R_tmp_in --no-save --slave`;
open INFILE,"< $R_tmp_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_tmp_in $R_tmp_out`;
close INFILE;
return($results[0]);
}


##############this subrotine is used to filter the allele freq based on MAXPOP of gnomad WES db###########

sub get_gnomad_freq{ 
my ($tabix_dir_sub,$filter_db_sub,$var_sub,$tabix_output_sub) = (@_);
my @var_info = split(/_/,$var_sub);
my $region_sub = "$var_info[0]:$var_info[1]-$var_info[1]";
# ~/gnomad/data/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz
`$tabix_dir_sub $filter_db_sub $region_sub > $tabix_output_sub`;
my $max_freq_sub=0;
open(INPUTtbx,"$tabix_output_sub");
my $flag=0;
while(my $linetbx = <INPUTtbx>){
chomp $linetbx;
my @infotbx = split(/\s+/,$linetbx);
if(("$infotbx[0]" eq "$var_info[0]")&&($infotbx[1] == $var_info[1])&&(uc($infotbx[2]) eq uc($var_info[2]))&&(uc($infotbx[3]) eq uc($var_info[3]))){

my @max_freq = split(/:/,$infotbx[9]);
$max_freq_sub = $max_freq[1];
last;
}
}
return($max_freq_sub);
}


###################this subrotine is used to annotate CADD score ###############
sub get_CADD_score{ 
my ($tabix_dir_sub,$var_sub,$tabix_output_sub) = (@_);
my @var_info = split(/_/,$var_sub);
my $CADD_db_sub;
if((length($var_info[2])==1)&&(length($var_info[3]))){
$CADD_db_sub = $opt_CADD_snp_db; 
}else{
$CADD_db_sub = $opt_CADD_indel_db; 
}
my $region_sub = "$var_info[0]:$var_info[1]-$var_info[1]";
`$tabix_dir_sub $CADD_db_sub $region_sub > $tabix_output_sub`;
my $CADD_sub="NA";
open(INPUTtbx,"$tabix_output_sub");
my $flag=0;
while(my $linetbx = <INPUTtbx>){
chomp $linetbx;
my @infotbx = split(/\s+/,$linetbx);
if(("$infotbx[0]" eq "$var_info[0]")&&($infotbx[1] == $var_info[1])&&(uc($infotbx[2]) eq uc($var_info[2]))&&(uc($infotbx[3]) eq uc($var_info[3]))){
#1       36244   G       C       -0.428735       0.329
my $convert = convert_CADD($infotbx[-1]);
$CADD_sub = "$convert\t$infotbx[-2]\t$infotbx[-1]";
last;
}
}
return($CADD_sub);
}

sub convert_CADD {
my $CADD_pred = shift(@_);
my $CADD_convert = 1-(10**(-($CADD_pred/10)));
return($CADD_convert);
}
