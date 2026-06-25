# Common functions
library(colorspace)

## Create tableone output into a dataframe for export as a kable output
convert_tableone_into_df <- function(dataset,
                                     vars,
                                     strata=NULL, 
                                     factorVars=NULL,
                                     exact_vars=NULL,
                                     nonnormal_vars=NULL,
                                     includeNA = FALSE,
                                     test = TRUE,
                                     testApprox = chisq.test,
                                     argsApprox = list(correct = TRUE),
                                     testExact = fisher.test,
                                     argsExact = list(workspace = 2 * 10^5),
                                     testNormal = oneway.test,
                                     argsNormal = list(var.equal = TRUE),
                                     testNonNormal = kruskal.test,
                                     argsNonNormal = list(NULL),
                                     outcome_names,
                                     smd = TRUE,
                                     addOverall = FALSE){
  if(is.null(strata)==T){
    overall <- tableone::CreateTableOne(vars = vars, 
                                        data = dataset,
                                        factorVars = factorVars,
                                        includeNA = includeNA, 
                                        testApprox = testApprox,
                                        argsApprox = argsApprox,
                                        testExact = testExact,
                                        argsExact = argsExact,
                                        testNormal = testNormal,
                                        argsNormal = argsNormal,
                                        testNonNormal = testNonNormal,
                                        argsNonNormal = argsNonNormal,
                                        smd = TRUE,
                                        addOverall = addOverall)
    
    bound_table <-  capture.output(x <- print(overall, quote = FALSE, noSpaces = TRUE ))
    names <- x %>% as.matrix()  %>% rownames
    values <- x %>% as.data.frame %>% `rownames<-`(NULL)
    table <- cbind.data.frame(names,x %>% as.data.frame()) %>% `rownames<-`(NULL)
    colnames(table) <- c("Variable",paste0("Overall (n=",table[1,2],")"))
    rownames(table) <- NULL
    table <- table[-1,]
    return(table)
  }
  if(is.null(strata)==F){
    overall <- tableone::CreateTableOne(vars = vars,
                                        strata = strata,
                                        data = dataset,
                                        factorVars = factorVars,
                                        includeNA = includeNA,
                                        test = test,
                                        testApprox = testApprox,
                                        argsApprox = argsApprox,
                                        testExact = testExact,
                                        argsExact = argsExact,
                                        testNormal = testNormal,
                                        argsNormal = argsNormal,
                                        testNonNormal = testNonNormal,
                                        argsNonNormal = argsNonNormal,
                                        smd = TRUE,
                                        addOverall = addOverall) 
    
    bound_table <-  capture.output(x <- print(overall, quote = FALSE, noSpaces = TRUE,exact=exact_vars,nonnormal = nonnormal_vars))
    table <- x[,]  %>% as.data.frame()
    names <- x %>% as.matrix()  %>% rownames
    table <- cbind.data.frame(names,table %>% as.data.frame()) %>% `rownames<-`(NULL) 
    if(addOverall ==T) {
      colnames(table) <- c("Variable",paste0(c("Overall",outcome_names)," (n=",table[1,2:c(ncol(table)-2)],")"),"P-value",'Test')
    } else {
      colnames(table) <- c("Variable",paste0(c(outcome_names)," (n=",table[1,2:c(ncol(table)-2)],")"),"P-value",'Test')
    }
    
    rownames(table) <- NULL
    table <- table[-1,]
    return(table)
  }
}


# Get presence absence data
get_presence_absence_matrix <- function(variable,df){
  presences <- str_split(unlist(c(df %>% dplyr::select(paste0(variable)))),";") %>% unlist(.)  %>% unique(.) %>% .[!. %in% "-"] %>% sort
  vector <- as.vector(tidyr::unite(df %>% dplyr::select(paste0(variable)),determinants,sep =";")) %>% unlist() %>% as.vector
  results <- matrix(nrow = nrow(df),ncol=length(presences))  %>% `colnames<-`(presences)
  for(i in presences){
    results[,paste0(i)] <- as.vector(ifelse(unlist(lapply(str_split(as.matrix(df[,paste0(variable)]),";"), function(x) i %in% x))==TRUE,1,0))
  }
  rownames(results) <- rownames(df)
  return(results)
}
# Mobsuite plasmid analyses
## Cluster
get_mobtyper_cluster_matrix <- function(df,mobtyper){
  plasmid_cluster <- mobtyper$primary_cluster_id %>% unique
  lapply(plasmid_cluster,FUN =function(x,df,mobtyper){
    isolates_with_plasmid <-  subset(mobtyper,primary_cluster_id == x) %>% select(isolate_no) %>% unlist
    row <- assign(paste0(x),ifelse(df$isolate_no %in% isolates_with_plasmid,1,0))
    return(row)},df,mobtyper) %>% do.call(bind_cols,.) %>% as.data.frame %>% `colnames<-`(plasmid_cluster) %>% `rownames<-`(df$isolate_no)
}

## Nearest neighbor
get_mobtyper_mash_nn_matrix <- function(df,mobtyper){
  nn_plasmid <- mobtyper$mash_nearest_neighbor %>% unique
  lapply(nn_plasmid,FUN =function(x,df,mobtyper){
    isolates_with_plasmid <-  subset(mobtyper,mash_nearest_neighbor == x) %>% select(isolate_no) %>% unlist
    row <- assign(paste0(x),ifelse(df$isolate_no %in% isolates_with_plasmid,1,0))
    return(row)},df,mobtyper) %>% do.call(bind_cols,.) %>% as.data.frame %>% `colnames<-`(nn_plasmid) %>% `rownames<-`(df$isolate_no)
}

# Genotype recode
recode_genotypes <- function(genotype_string){
  recode(genotype_string,
       "ompK35_rare_truncation"='Rare ompK35 truncation',
       "ompK36_rare_truncation"='Rare ompK36 truncation', 
       "OmpK36_c25t" = "ompK36 C25T",
       "Rare_ompK35_insertion" = "Rare ompK35 insertion",
       'Rare_ompK36_insertion' = 'Rare ompK36 insertion',
       "ompK36_T136TDT" = 'ompK36 T136TDT',
       "ompK36_D135DGD" = "ompK36 D135DGD",
       "ompK36_D135DD" ="ompK36 D135DD",
       "blaSHV_C-112A" = 'blaSHV C112A') %>% gsub("Rare_","Rare ",.)
}




# Format carbapenemase markdown

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

format_carbapenemase <- function(carbapenemase_str){

# Split
carbapenemase_markdown <- str_split(carbapenemase_str,pattern = "&",simplify=T)
carbapenemase_markdown[,1] <- sapply(carbapenemase_markdown[,1],format_blaKPC_md)
carbapenemase_markdown[,2] <- sapply(carbapenemase_markdown[,2],format_blaKPC_md)
colnames(carbapenemase_markdown) <- c("c1","c2")

# Apply throughout string
results <- apply(carbapenemase_markdown,1,FUN=function(x){
  if(x['c2']==''){
    paste0(x['c1'])
  } else {
    paste0(x['c1']," & ",x['c2'])
  }
})

return(results)

}

# Format all
format_blaKPC_genotypes_md <- function(gene) {
  gene <- trimws(gene,which='both')
  if(grepl("^bla", gene)){
    allele <- gsub("bla", "", gene)
    final <- paste0("<i>bla</i><sub>", allele, "</sub>") %>% trimws(.,which='both')
  } 
  if(grepl("Rare bla",gene)==T){
    allele <- gsub("Rare bla", "", gene)
    final <- paste0("Rare <i>bla</i><sub>",allele,"</sub>")   %>% trimws(.,which='both') 
  } 
  if(grepl("ompK",gene) == T){
    final <- gsub("ompK35","<i>ompK35</i>",gene) %>% gsub("ompK36","<i>ompK36</i>",.)
  } 
  
  if(grepl("OmpK",gene)==T){
    final <- gene
  }
  final 
}
 
build_color_palette <- function(df,variable,other_groups_named_color,iwanthues_specs,color_num = total_number,starting_method,distance){
  
  # Get ranking
  df_split <-  df %>%
    separate_rows(!!sym(variable), sep = " & ")
  
  count_of_var <- df_split %>% 
    count(study,.data[[variable]]) %>% 
    complete(study,.data[[variable]],fill = list(n=0)) %>% 
    group_by(study) %>% 
    mutate(perc = n / sum(n) * 100) 
  
  new_order <- count_of_var %>% subset(study =='2021-23') %>% arrange(-n) %>% subset(n>0) %>% .[[variable]] %>% subset(!. %in% names(other_groups_named_color))
  zero_new_order_old <- zero_new_order_old <- count_of_var %>% subset(study =='2014-15') %>% subset(!get(variable) %in% c(new_order,names(other_groups_named_color))) %>% arrange(-n) %>% .[[variable]]
 
  total_colors_to_pick <- length(c(new_order,zero_new_order_old))
  if(color_num < total_colors_to_pick){
    stop("Not enough requested numbers")
  }
  
  
  # Get colors
  if(color_num %in% c(NULL,"min")){
    color_num = total_colors_to_pick
  }
  
  large_pal <- do.call(
    hues::iwanthue,
    c(list(n = color_num), iwanthues_specs))
  
  # Get distance matrix 
  if(distance == "normal"){
    
    lab <- colorspace::coords(as(colorspace::hex2RGB(large_pal), "LAB"))
    
    dmat <- as.matrix(dist(lab))
  }
  if(distance == ("colorblind")){
    cb <- colorspace::simulate_cvd(large_pal,
                                   interpolate_cvd_transform(tritanomaly_cvd, severity = 0.6))
    
    cb_lab <-  colorspace::coords(as(colorspace::hex2RGB(cb), "LAB"))
    
    dmat <- as.matrix(dist(cb_lab))
  }
  
  # Get max distance value
  pick_distant <- function(candidates,selected,dmat){
    if(length(selected) == 0 ){
      return(candidates[1])
    }
    
    # For each candidate, compute its distance to the closest already-selected color
    min_dist <- apply(dmat[candidates,selected,drop=F],1,min,na.rm=T)
    
    candidates[which.max(min_dist)]
    
  }
  
  # Select 2021-23 color palete
  used <- rep(F,nrow(dmat))
  selected_202123 <- integer(length(new_order))
  
  ## Starting value
  diag(dmat) <- NA 
  
  # Most distinct   (purple ->  green -> red ->  cyan -> yellow)
  if(starting_method == "max_mean"){
    seed <- which.max(rowMeans(dmat, na.rm = TRUE))
  } 
  
  # Most central (olive-brown → purple → cyan → orange → green → pink)
  if(starting_method == "min_mean") {
    seed <- which.min(rowMeans(dmat, na.rm = TRUE))
  } 
  
  # Most isolated (teal → magenta → yellow → blue → orange → green)
  if(starting_method == "farthest_neighbor") {
    seed <- which.max(apply(dmat, 1, min, na.rm = TRUE))
  }
  # Broadly distinct (blue → orange → green → purple → cyan → red)
  if(starting_method == "max_median"){
    seed <- which.max(apply(dmat, 1, median,na.rm = TRUE))
  }
  # Broadly central (gray-purple → green → orange → blue → pink → cyan)
  if(starting_method == "min_median"){
    seed <- which.min(apply(dmat, 1, median,na.rm = TRUE))
  }
  if(!starting_method %in% c('max_mean','min_mean','farthest_neighbor','max_median','min_median')){
    stop(paste0("Must be in: 'max_mean','min_mean','farthest_neighbor','max_median','min_median'"))
  }
  
  selected_202123[1] <- seed 
  used[seed] <- TRUE
  
  ## Get distinct
  for (i in 2:length(new_order)) {
    candidates <- which(!used)
    selected_202123[i] <- pick_distant(candidates, selected = selected_202123[seq_len(i - 1)], dmat)
    used[selected_202123[i]] <- TRUE
  }
  
  ## Add names
  selected_202123_values <- large_pal[selected_202123]
  names(selected_202123_values) <- new_order
  
  # Select 2014-15 color palette
  selected_201415 <- integer(length(zero_new_order_old))
  selected_201415_values <- setNames(character(0), character(0))  
  
  if(length(selected_201415)>0){
    
  for (i in 1:length(zero_new_order_old)) {
    candidates <- which(!used)
    selected_201415[i] <- pick_distant(candidates, selected_201415[seq_len(i - 1)], dmat)
    used[selected_201415[i]] <- TRUE
  } 
  
  ## Add names
  selected_201415_values <- large_pal[selected_201415]
  names(selected_201415_values) <- zero_new_order_old 
  
  }
  
  ## Create named palette
  palette_final <- c(selected_202123_values, selected_201415_values,other_groups_named_color) 
  
  return(palette_final)
  
}

# Build a global color palette
build_color_palette_global <- function(df,variables,other_groups,iwanthues_specs,color_num = total_number,starting_method,distance,
                                       group_locked = F,group_sizes=NULL,group_method,string_split=NULL){
    # Get ranking
  if(is.null(string_split) ==F){
    total_colors_to_pick <- length(df[,c(variables)] %>% unlist %>% subset(!. %in% other_groups) %>% unique %>% str_split(.,pattern = string_split) %>% unlist %>% unique) 
  } else {
    total_colors_to_pick <- length(df[,c(variables)] %>% unlist %>% subset(!. %in% other_groups) %>% unique)
  }
  
  if(color_num < total_colors_to_pick){
    stop(paste0("Not enough requested numbers. Consider a minimum of: ",total_colors_to_pick))
  }
  
  
  # Get colors
  if(color_num %in% c(NULL,"min")){
    color_num = total_colors_to_pick
  }
  
  large_pal <- do.call(
    hues::iwanthue,
    c(list(n = color_num), iwanthues_specs))
  
  # Get distance matrix 
  if(distance == "normal"){
    
    lab <- colorspace::coords(as(colorspace::hex2RGB(large_pal), "LAB"))
    
    dmat <- as.matrix(dist(lab))
  }
  if(distance == ("colorblind")){
    cb <- colorspace::simulate_cvd(large_pal,
                                   interpolate_cvd_transform(tritanomaly_cvd, severity = 0.6))
    
    cb_lab <-  colorspace::coords(as(colorspace::hex2RGB(cb), "LAB"))
    
    dmat <- as.matrix(dist(cb_lab))
  }
  
  # Get max distance value
  pick_distant <- function(candidates,selected,dmat){
    if(length(selected) == 0 ){
      return(candidates[1])
    }
    
    # For each candidate, compute its distance to the closest already-selected color
    min_dist <- apply(dmat[candidates,selected,drop=F],1,min,na.rm=T)
    
    candidates[which.max(min_dist)]
    
  }
  
  # Select color palete
  used <- rep(F,nrow(dmat))
  selected <- integer(total_colors_to_pick)
  
  ## Starting value
  diag(dmat) <- NA 
  
  if(!group_locked){
  seed <- get_seed(dmat,starting_method) 
  selected[1] <- seed 
  used[seed] <- TRUE
  
  ## Get distinct
  for (i in 2:total_colors_to_pick) {
    candidates <- which(!used)
    selected[i] <- pick_distant(candidates, selected = selected[seq_len(i - 1)], dmat)
    used[selected[i]] <- TRUE
  }
  
  ## Create named palette
  palette_final <- large_pal[selected]
  
  }
  
  if(group_locked==T){
    
    if (is.null(group_method)) {
      stop("group_method must be provided when group_locked = TRUE")
    }
    
    if(is.null(group_sizes)){
      stop("group_sizes must be provided when group_locked = TRUE")
    }
    
    used <- rep(FALSE, nrow(dmat))
    selected <- list()
    
    # shared greedy context across all groups
    global_selected <- integer(0)
    
    if(group_method == 'sequential'){
    for(g in seq_along(group_sizes)){
      
      selected[[g]] <- integer(group_sizes[g])
      
      # ---- group seed (best remaining global candidate)
      pool <- which(!used)
      
      remaining_dmat <- dmat[pool, pool, drop = FALSE]
      
      seed_g <- get_seed(remaining_dmat,starting_method) 
      
      selected[[g]][1] <- seed_g
      used[seed_g] <- TRUE
      global_selected <- c(global_selected, seed_g)
      
      # ---- greedy fill within group, but global constraint
      for(i in 2:group_sizes[g]){
        
        pool <- which(!used)
        
        selected[[g]][i] <- pick_distant(
          candidates = pool,
          selected = global_selected,
          dmat = dmat
        )
        
        used[selected[[g]][i]] <- TRUE
        global_selected <- c(global_selected, selected[[g]][i])
      }
    }
    
    palette_final <- large_pal[unlist(selected)]
    }
    if(group_method == "connected"){
      
      if (is.null(group_sizes)) {
        stop("group_sizes must be provided when group_locked = TRUE")
      }
      
      selected_groups <- vector("list", length(group_sizes))
      
      global_selected <- integer(0)
      
      for (g in seq_along(group_sizes)) {
        
        g_size <- group_sizes[g]
        selected_groups[[g]] <- integer(g_size)

        available <- setdiff(seq_len(nrow(dmat)), global_selected)
        
        # ---- group seed from FULL remaining space (not a shrinking matrix)
        seed_candidates <- available
        
        seed_sub_dmat <- dmat[seed_candidates, seed_candidates, drop = FALSE]
        seed_local <- get_seed(seed_sub_dmat, starting_method)
        
        seed_g <- seed_candidates[seed_local]
        
        selected_groups[[g]][1] <- seed_g
        global_selected <- c(global_selected, seed_g)
        
        # ---- greedy fill WITH global constraint, but stable pool
        for (i in 2:g_size) {
          
          available <- setdiff(seq_len(nrow(dmat)), global_selected)
          
          selected_groups[[g]][i] <- pick_distant(
            candidates = available,
            selected = global_selected,
            dmat = dmat
          )
          
          global_selected <- c(global_selected, selected_groups[[g]][i])
        }
      } 
      
      palette_final <-  large_pal[unlist(selected_groups)]
    }
  }
  
  return(palette_final)
  
}

get_seed <- function(dmat,starting_method){
  # Most distinct   (purple ->  green -> red ->  cyan -> yellow)
  if(starting_method == "max_mean"){
    seed <- which.max(rowMeans(dmat, na.rm = TRUE))
  } 
  
  # Most central (olive-brown → purple → cyan → orange → green → pink)
  if(starting_method == "min_mean") {
    seed <- which.min(rowMeans(dmat, na.rm = TRUE))
  } 
  
  # Most isolated (teal → magenta → yellow → blue → orange → green)
  if(starting_method == "farthest_neighbor") {
    seed <- which.max(apply(dmat, 1, min, na.rm = TRUE))
  }
  # Broadly distinct (blue → orange → green → purple → cyan → red)
  if(starting_method == "max_median"){
    seed <- which.max(apply(dmat, 1, median,na.rm = TRUE))
  }
  # Broadly central (gray-purple → green → orange → blue → pink → cyan)
  if(starting_method == "min_median"){
    seed <- which.min(apply(dmat, 1, median,na.rm = TRUE))
  }
  if(!starting_method %in% c('max_mean','min_mean','farthest_neighbor','max_median','min_median')){
    stop(paste0("Must be in: 'max_mean','min_mean','farthest_neighbor','max_median','min_median'"))
  }
  return(seed)
}

# Get ranking
get_ranking <- function(df,variable,other_groups){
# Get ranking
non_split <- df[[variable]] %>% subset(grepl("&",.)) %>% unique %>% sort
df_split <-  df %>%
  separate_rows(!!sym(variable), sep = " & ")

count_of_var <- df_split %>% 
  count(study,.data[[variable]]) %>% 
  complete(study,.data[[variable]],fill = list(n=0)) %>% 
  group_by(study) %>% 
  mutate(perc = n / sum(n) * 100) 

new_order <- count_of_var %>% subset(study =='2021-23') %>% arrange(-n) %>% subset(n>0) %>% .[[variable]] %>% subset(!. %in% other_groups)
zero_new_order_old <- count_of_var %>%  dplyr::filter(study == "2014-15") %>%  dplyr::filter(!.data[[variable]] %in% c(new_order, other_groups)) %>%  dplyr::arrange(desc(n)) %>%   dplyr::pull(.data[[variable]])
order <- unique(c(new_order, zero_new_order_old, non_split, other_groups))
return(order) 
}

# Resistance related functions
recode_resistance_variables <- function(variables){
  recode(variables,
         "CZA_dich_num"='Ceftazidime-avibactam',
         "CST_dich_num"='Colistin', 
         'MERO_dich_num'='Meropenem',
         'IMI_dich_num'='Imipenem',
         "blbli_dich_num"='KPC-inhibiting BL/BLIs',
         "MVB_dich_num"='Meropenem-vaborbactam',
         "IR_dich_num"='Imipenem-relebactam',
         "PLZ_dich_num"='Plazomicin',
         "TMP_SMX_dich_num"='Trimethoprim-Sulfamethoxazole',
         "OMC_dich_num" = "Omadacycline",
         "DLX_dich_num" = "Delafloxacin",
         "ERV_dich_num" = "Eravacycline",
         "FOS_dich_num" = "Fosfomycin",
         "CT_dich_num" = "Ceftolozane-tazobactam",
         "FDC_dich_num"="Cefiderocol")
}

# Functions for comparing resistance
get_frequency_stats <- function(variable,df,name){
  freq <- df[[variable]] %>% sum(.,na.rm=T)
  tested <- df[[variable]] %>% subset(is.na(.)==F) %>% length
  not_tested <- nrow(df) - tested
  prop <- round(freq / tested * 100,2) 
  data.frame(variable = variable,study = name,freq,tested,prop,not_tested)
}

get_frequency_stats_compare <- function(variable,df,comparitor){
  renamed_variable <- recode_resistance_variables(variable)
  data_table <- table(df[[variable]],df[[comparitor]],useNA = 'no') 
  test <- if(sum(unlist(data_table) <5)>0){
    fisher.test(data_table,simulate.p.value = T)
  } else {
    chisq.test(data_table)
  }
  group1 <- colnames(data_table)[1]
  group2 <- colnames(data_table)[2]
  
  group1_total <- sum(data_table[,1])
  group2_total <- sum(data_table[,2])
  
  group1_count <- data_table[,1] %>% subset(names(.)==1)
  group2_count <- data_table[,2] %>% subset(names(.)==1)
  
  group1_perc <- round(c(group1_count/group1_total)*100,1)
  group2_perc <- round(c(group2_count/group2_total)*100,1)
  
  sig_test <- ifelse(class(test)=="htest" & test$method == "Fisher's Exact Test for Count Data","Fisher's Exact Test","Chi-squared")
  p_value <- test$p.value
  
  data.frame(name = renamed_variable, variable=variable,group1=group1,group1_total=group1_total,group1_count=group1_count,group1_perc=group1_perc,
             group2=group2,group2_total=group2_total,group2_count=group2_count,group2_perc=group2_perc,
             test=sig_test,p_value=p_value)
}