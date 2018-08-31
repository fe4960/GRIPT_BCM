#!/bin/sh

exec 1> 4run_simulation_example.log
exec 2>> 4run_simulation_example.log

num=10  
ethn=NFE
sex_ratio=0.5
AF_cutoff=0.2
sim_db="Data/simulation_db/gnomad_PASS_sort_subpop_snpEff_ns_All_CADD_all.screen_0.0001"
output="Examples/simulation10_NFE"
mkdir $output

date_start=$(date +%s)
echo -e "startTime:::::::\c" ; date

perl 3simulation.pl -i $sim_db -o $output/gnomad -c $AF_cutoff -n $num -e $ethn -r $sex_ratio

echo -e "endTime:::::::\c" ; date
date_end=$(date +%s)
echo "total time: $((date_end-date_start)) s"


