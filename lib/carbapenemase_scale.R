# Genotype scale
dataset_for_creating_cb_scale <- readRDS("./data/dataset/CA_LTACH_first_clinical_per_quarter_CRKP_genomes_merged_studies_metadata.RDS") 
cab_dataset <- readRDS("./data/carbapenemase_content/carbapenemase_content_merged_studies.RDS")
dataset_for_creating_cb_scale <- left_join(dataset_for_creating_cb_scale,cab_dataset)

format_blaKPC_md <- function(gene) {
  gene <- trimws(gene,which='both')
  if(grepl("^bla", gene)){
    allele <- gsub("bla", "", gene)
    paste0("<i>bla</i><sub>", allele, "</sub>") %>% trimws(.,which='both')
  }  else {
    if(grepl("No carbapenemase",gene) == T ){
      paste0("No carbapenemase") 
    } else {
      paste0("") 
    }
  } 
}
carbapenemase_markdown <- str_split(dataset_for_creating_cb_scale$carbapenemase,pattern = "&",simplify=T)
carbapenemase_markdown[,1] <- sapply(carbapenemase_markdown[,1],format_blaKPC_md)
carbapenemase_markdown[,2] <- sapply(carbapenemase_markdown[,2],format_blaKPC_md)
colnames(carbapenemase_markdown) <- c("c1","c2")
dataset_for_creating_cb_scale$carbapenemase_md <- apply(carbapenemase_markdown,1,FUN=function(x){
  if(x['c2']==''){
    paste0(x['c1'])
  } else {
    paste0(x['c1']," & ",x['c2'])
  }
})
genotypes <- c(table(dataset_for_creating_cb_scale$carbapenemase_md %>% subset(.!="No carbapenemase")) %>% sort(decreasing = T) %>% names,"No carbapenemase")
genotype_colors <- setNames(c(hues::iwanthue(length(genotypes)-1,hmin = 5,cmin = 5,lmin = 5),"gray"), genotypes)

genotype_scale <- scale_fill_manual(breaks = genotypes,values = genotype_colors,name = "Carbapenemase", guide = guide_legend(ncol=4, title.position = "top", label.position = "right",order=2)) 

genotype_scale_4_row <- scale_fill_manual(breaks = genotypes,values = genotype_colors,name = "Carbapenemase", guide = guide_legend(nrow=4, title.position = "top", label.position = "right",order=2)) 
