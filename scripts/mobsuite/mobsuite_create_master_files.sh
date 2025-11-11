#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/data/mobsuite

# Penn KPC
## Mobtyper
awk '(NR==1)' ./PCMP_H9/mobtyper_results.txt > Penn_KPC_mobtyper_results.txt

for isolate in $(ls -d PCMP_H*/ | cut -d/ -f1)
do 
        echo $isolate
        awk '(NR>1)' ./$isolate/mobtyper_results.txt  >> Penn_KPC_mobtyper_results.txt
done

## Contig report
awk '(NR==1)' ./PCMP_H9/contig_report.txt > Penn_KPC_MOB_suite_contig_report.txt

for isolate in $(ls -d PCMP_H*/ | cut -d/ -f1)
do
        echo $isolate
        awk '(NR>1)' ./$isolate/contig_report.txt >> Penn_KPC_MOB_suite_contig_report.txt
done

# INTEGRATE
## Mobtyper
awk '(NR==1)' ./INT_CRE_143/mobtyper_results.txt > INTEGRATE_mobtyper_results.txt

for isolate in $(ls -d INT_CRE*/ | cut -d/ -f1)
do 
        echo $isolate
        awk '(NR>1)' ./$isolate/mobtyper_results.txt  >> INTEGRATE_mobtyper_results.txt
done

## Contig report
awk '(NR==1)' ./INT_CRE_143/contig_report.txt > INTEGRATE_MOB_suite_contig_report.txt

for isolate in $(ls -d INT_CRE*/ | cut -d/ -f1)
do
        echo $isolate
        awk '(NR>1)' ./$isolate/contig_report.txt >> INTEGRATE_MOB_suite_contig_report.txt
done