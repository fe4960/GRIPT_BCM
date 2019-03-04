GRIPT: A novel case-control analysis method for Mendelian disease gene discovery

Wang, Jun; Zhao, Li; Wang, Xia; Chen, Yong; Xu, Mingchu; Soens, Zachry T.; Ge, Zhongqi; Wang, Peter Ronghan; Wang, Fei; Chen, Rui

Baylor College of Medicine

The GRIPT package is distributed under the GNU General Public License v3.0 

####################
1. Download the GRIPT_BCM from github

2.Download Data.tar.gz and Examples.tar.gz from https://zenodo.org/record/2583149#.XH1sucBKiJA

tar -xvf Data.tar.gz

tar -xvf Examples.tar.gz 

3. move the unzipped "Data" and "Examples" folders (from 2.) to "GRIPT_BCM" folder

mv Data  GRIPT_BCM
mv Examples GRIPT_BCM

########################################
GRIPT
########################################

-to see the parameters/options of the running GRIPT script

perl 2run_GRIPT.pl

-example command to run GRIPT

perl 2run_GRIPT.pl -case_in $PWD/Examples/case/AF0.005_N600/RPE65/ -case_out $PWD/Examples/case/AF0.005_N600/RPE65_out -control_in $PWD/Examples/control/AF0.005_N5000/ -control_out $PWD/Examples/control/AF0.005_N5000_out

-example of GRIPT input folder:

case folder: 

Examples/case/AF0.005_N600/RPE65/

Examples/case/AF0.0001_N600/TINF2/

control folder:

Examples/control/AF0.005_N5000/
               
Examples/control/AF0.0001_N5000/ 

-example of GRIPT output folder

Examples/case/AF0.005_N600/RPE65_out

Examples/control/AF0.005_N5000_out

-example of the prioritized gene list by GRIPT

Examples/case/AF0.005_N600/RPE65_out/recessive_model.txt


##############################
preprocessing step
##############################

-tools and data needed to run the preprocessing script:

1. download the CADD score from the CADD websit: https://cadd.gs.washington.edu/download
including "All possible SNVs of GRCh37/hg19" and "48M InDels to initiate a local setup"  

2. install tabix

3. install R

4. install snpEff

-to see the parameters/options of the preprocessing script

perl 1preprocess_vcf.pl 

-example command to run the preprocessing script

perl 1preprocess_vcf.pl -i Examples/case/AF0.02_N2_vcf  -o Examples/case/AF0.02_N2_score -sfj snpEff/snpEff.jar -sfc snpEff/snpEff.config -rt Atlas -cd_snp CADD_score/whole_genome_SNVs.tsv.gz  -cd_ind CADD_score/InDels.tsv.gz  -gn Data/gnomad/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz


#####################################
WES simulator
####################################

-simulator script

perl 3simulation.pl

-example of shell command to run simulator

sh 4run_simulation_example.sh

-the input data files of the simulator script:

1) the protein change or splicing affect variants with max population allele frequency = 0.00001 in ExAC database with CADD score and annotated by snpEff

Data/simulation_db/ExAC.r0.3.1_sort_nf.snpEff.score.screen.max0.00001

2) the protein change or splicing affect variants with max population allele frequency = 0.0001  in ExAC database with CADD score and annotated by snpEff

Data/simulation_db/ExAC.r0.3.1_sort_nf.snpEff.score.screen.max0.0001

3)the protein change or splicing affect variants with max population allele frequency = 0.005 in ExAC database with CADD score and annotated by snpEff

Data/simulation_db/ExAC.r0.3.1_sort_nf.snpEff.score.screen.max0.005

4)the protein change or splicing affect variants with max population allele frequency = 0.0001 in gnomAD WES database with CADD score and annotated by snpEff

Data/simulation_db/gnomad_PASS_sort_subpop_snpEff_ns_All_CADD_all.screen_0.0001

5)all the variants in gnomAD WES database with CADD score and annotated by snpEff

Data/simulation_db/gnomad_PASS_sort_subpop_snpEff_ns_All_CADD_all

-example shell command to generate the protein change or splicing affect variants with a given max population allele frequency in gnomAD WES database with CADD score and annotated by snpEff

sh 5make_simulation_db_example.sh
