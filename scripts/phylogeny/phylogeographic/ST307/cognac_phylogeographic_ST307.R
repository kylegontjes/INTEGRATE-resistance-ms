
# Cognac Run
# Start Date: 12/16/25
# Population: ST307 isolates
# Activate
library(tidyverse)
library(cognac)
Sys.info()
sessionInfo()

setwd('/nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/')

outDir = '/nfs/turbo/umms-esnitkin/Project_INTEGRATE/Analysis/Resistance_dynamics/INTEGRATE-resistance-ms/data/phylogeny/phylogeographic/ST307/'

genomic_id = readRDS('./scripts/phylogeny/phylogeographic/ST307/ST307_isolate_list.RDS')  
fasta_files = readRDS('./scripts/phylogeny/phylogeographic/ST307/ST307_genomes_path.RDS')
gff_files =  readRDS('./scripts/phylogeny/phylogeographic/ST307/ST307_annotations_path.RDS')

# Run Conac requiring at least 500 genes are included in the alignment 
cognac(
  fastaFiles = fasta_files ,
  featureFiles = gff_files,
  outDir        = outDir,
  minGeneNum    = 500,
  maxMissGenes  = 0.05, 
  njTree     = TRUE,
  mapNtToAa = TRUE,
  keepTempFiles = TRUE,
  outGroup = 'INT_CRE_1243'
)
 
