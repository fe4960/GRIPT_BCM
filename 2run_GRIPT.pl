#!/usr/bin/perl -w
use strict;
use Cwd;
use strict;
use Getopt::Long;
use List::Util qw( min max );

use vars qw( $opt_help $opt_case_input $opt_case_output $opt_control_input $opt_control_output  $opt_model  $opt_var_cutoff $opt_gript_dir);

if (!GetOptions(
                "help|h",
                "case_input|case_in:s",
                "case_output|case_out:s",
		"control_input|control_in:s",
		"control_output|control_out:s",
                "model|mod:s",
                "var_cutoff|vc:s",
                "gript_dir|gript"
            )) {
    &usage(1);
}


if ($opt_help) {
    &usage(0);
}

sub usage {


                print "
       GRIPT running script
       Authors: Jun Wang;  Rui Chen Lab
       Date: 8-30-2018

       Usage: -case_in:         the absolute directory of case input folder 
              -case_out:        the absolute directory of case output folder 
              -control_in:      the absolute directory of control input folder
              -control_out:     the absolute directory of control output folder 
              -gript:           the absolute directory of GRIPT folder
              -vc:              variant score cutoff. default  = 1.
              -mod:             inheritance model. options: \"recessive_model\" and \"dominant_model\". default = recessive_model
	\n";
                exit;

}


if (($opt_case_input eq "")||($opt_case_output eq "") || ($opt_control_input eq "") || ($opt_control_output eq "")) {&usage(0); exit;}

if ( ! $opt_var_cutoff ){ $opt_var_cutoff=0;}

if (! $opt_model ){$opt_model = "recessive_model";}

if (! $opt_gript_dir) {$opt_gript_dir = "GRIPT";}


system("cd $opt_gript_dir; java Score_cutoff -casein $opt_case_input -controlin $opt_control_input -inheritance $opt_model -caseout $opt_case_output -controlout $opt_control_output -cutoff $opt_var_cutoff");

my $case_in1 = "$opt_case_output/geneScore";
my $control_in1 = "$opt_control_output/geneScore";

system("cd $opt_gript_dir; java Statistic_Fisher -casein $case_in1 -controlin $control_in1 -inheritance $opt_model -caseout $opt_case_output -controlout $opt_control_output");

my $case_in2 = "$case_in1"."Matrix";
my $control_in2 = "$control_in1"."Matrix";

 system("cd $opt_gript_dir; java Identify_Fisher_new -casein $case_in2 -inheritance $opt_model -out $opt_case_output");

my $case_in3 = "$opt_case_output/variantScore";

#`rm -r $opt_control_output`;
#`tar -cvf $case_in1.tar $case_in1`;
#`tar -cvf $case_in2.tar $case_in2`;
#`tar -cvf $case_in3.tar $case_in3`;
#`gzip $case_in1.tar`;
#`gzip $case_in2.tar`;
#`gzip $case_in3.tar`;
#`rm -r $case_in1`;
#`rm -r $case_in2`;
#`rm -r $case_in3`;

