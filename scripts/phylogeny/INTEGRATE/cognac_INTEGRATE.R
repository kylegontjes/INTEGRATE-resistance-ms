
# Cognac Run
# Start Date: 12/05/25
# Population: INTEGRATE isolates
# Activate
library(tidyverse)
library(cognac)
Sys.info()
sessionInfo()

setwd('/nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/')

outDir = '/nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/data/phylogeny/INTEGRATE/'

genomic_id = readRDS('./scripts/phylogeny/INTEGRATE/INTEGRATE_isolate_list.RDS')  

fasta_files = paste0('/nfs/turbo/umms-esnitkin/Project_INTEGRATE/Sequence_data/assembly/illumina/prokka/',genomic_id,'/',genomic_id,'.fna')
gff_files = paste0('/nfs/turbo/umms-esnitkin/Project_INTEGRATE/Sequence_data/assembly/illumina/prokka/',genomic_id,'/',genomic_id,'.gff')

# Run Conac requiring at least 500 genes are included in the alignment 
cognac(
  fastaFiles = fasta_files ,
  featureFiles = gff_files,
  outDir        = outDir,
  minGeneNum    = 500,
  maxMissGenes  = 0.05, 
  njTree     = TRUE,
  mapNtToAa = TRUE,
  keepTempFiles = TRUE
)
 
