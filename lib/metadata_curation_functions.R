# Metadata curation functions

## Recode antibiotic exposure history data
### Recode history
binarize_abx_history <- function(variable,cutoff,df){
  history <- ifelse(df %>% select(variable) %>% unlist  >=cutoff,1,0) %>% `names<-`(NULL)
  return(history)
}

# Factorize minimum-inhibitory concentration (MIC) data
factorize_MIC <- function(dataset,resistance){
  values <- as.vector(unique(dataset[,resistance]))
  combination <- ifelse(length(grep("/",values))>0,"yes","no")
  if(combination == "no"){
    level <- c(
      #Endpoint 
      values[grep("≤|<|<=",values)], 
      #if non endpoint included 
      #non-endpoint sorted
      as.character(sort(as.numeric(values[!grepl("≤|<|>|>=|≥", values)]))),
      #Endpoint
      values[grep(">|>=|≥",values)])
  }
  if(combination == "yes"){
    values_no_in <- gsub("/.*","",values)
    in_value <- gsub(".*\\/","",values)
    level <- c(
      #Endpoint 
      values_no_in[grep("≤",values_no_in)], 
      #if non endpoint included 
      #non-endpoint sorted
      as.character(sort(as.numeric(values_no_in[!grepl("≤|>", values_no_in)]))),
      #Endpoint
      values_no_in[grep(">",values_no_in)])
    level <- paste0(level,"/",in_value)
  }    
  factorized <- factor(as.vector(dataset[,resistance]),levels = level)
  return(factorized)
}

# Binarize MIC data
binarize_MIC <- function(dataset,resistance){
  binarized <- ifelse(dataset[,resistance]=="Susceptible","Susceptible","Non-Susceptible")
  return(binarized)
}