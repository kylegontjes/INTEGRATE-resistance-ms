#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/data/amrfinder

# Penn KPC 
awk '(NR==1)' PCMP_H7_amrfinder.out > Penn_KPC_amrfinder_results.txt

for isolate in $(ls PCMP_H* )
do 
        echo $isolate
        awk '(NR>1)' $isolate  >> Penn_KPC_amrfinder_results.txt
done
 
# INTEGRATE
awk '(NR==1)' INT_CRE_46_amrfinder.out > INTEGRATE_amrfinder_results.txt

for isolate in $(ls INT_CRE_* )
do 
        echo $isolate
        awk '(NR>1)' $isolate  >> INTEGRATE_amrfinder_results.txt
done