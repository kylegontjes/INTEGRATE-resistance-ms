#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/data/public_genomes/kleborate/results

# Penn KPC
## Mobtyper
awk '(NR==1)' ./partition_120/klebsiella_pneumo_complex_output.txt > NCBI_pathogens_Kpn_kleborate.txt

for isolate in $(ls -d partition*/ | cut -d/ -f1 | sort)
do 
        echo $isolate
        awk '(NR>1)' ./$isolate/klebsiella_pneumo_complex_output.txt  >> NCBI_pathogens_Kpn_kleborate.txt
done