#################################################################
# Prepare four sets of feed efficiency and RNA data
#################################################################
source("data_for_asso_fe_rna.R")

#################################################################
# 3. Run negative binomial regression analysis
#################################################################

#
# 3.1 Negative binomial regression models
#

# negative binomial regression model
calc_gene_stat <- function(dat.y, mod, sequencing_depth) {
  my_dispersions <- estimateDisp(dat.y, design = mod)
  my_fits <- glmFit(dat.y, design = mod, 
                    dispersion = my_dispersions$tagwise.dispersion, 
                    offset = log(sequencing_depth$N))
  my_tests <- glmLRT(my_fits, coef = "zFER")
  my_stats <- my_tests$table %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(gene = rowname, 
           p.value = PValue) %>% 
    tbl_df() %>% 
    select(gene, logFC, p.value) %>% 
    left_join(ensembl_ncbi_final, by = "gene") %>% 
    arrange(gene)
  
  my_stats$FDR <- p.adjust(my_stats$p.value, method = "BH")
  
  return(my_stats)
}

# negative binomial regression for 1M and 3M group
build_mod_1 <- function(dat.x) {
  fer <- dat.x$fer
  zFER <- ( fer - mean(fer) ) / sd(fer)
  
  Batch <- factor(dat.x$Reagent)
  mod <- model.matrix(~Batch + zFER)
  
  return(mod)
}

run_nb_1 <- function(dat.x, dat.y, sequencing_depth) {
  mod <- build_mod_1(dat.x)
  my_stats <- calc_gene_stat(dat.y, mod, sequencing_depth)
  
  return(my_stats)
}

# negative binomial regression for 2M group
build_mod_2 <- function(dat.x) {
  fer <- dat.x$fer
  zFER <- ( fer - mean(fer) ) / sd(fer)
  
  Batch <- factor(dat.x$Batch)
  mod <- model.matrix(~Batch + zFER)
  
  return(mod)
}

run_nb_2 <- function(dat.x, dat.y, sequencing_depth) {
  mod <- build_mod_2(dat.x)
  
  my_stats <- calc_gene_stat(dat.y, mod, sequencing_depth)
  
  return(my_stats)
}

# negative binomial regression for 6MA and 6MC together
build_mod_3 <- function(dat.x) {
  fer <- dat.x$fer
  zFER <- ( fer - mean(fer) ) / sd(fer)
  
  Group <- factor(dat.x$group, levels = c("6MA", "6MC"))
  Batch <- factor(dat.x$Reagent)
  mod <- model.matrix(~Group + Batch + zFER)
  
  return(mod)
}

run_nb_3 <- function(dat.x, dat.y, sequencing_depth) {
  mod <- build_mod_3(dat.x)
  
  my_stats <- calc_gene_stat(dat.y, mod, sequencing_depth)
  
  return(my_stats)
}

#
# 3.2 Run negative binomial regression
#

# 1M group
xiang_cluster_1 <- new_cluster(n = 22) %>% 
  cluster_library("edgeR") %>% 
  cluster_library("limma") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("purrr") %>% 
  cluster_library("tibble") %>% 
  cluster_copy("build_mod_1") %>% 
  cluster_copy("run_nb_1") %>% 
  cluster_copy("dat.rna.1M") %>% 
  cluster_copy("ensembl_ncbi_final") %>% 
  cluster_copy("sequencing_depth_1M") %>% 
  cluster_copy("calc_gene_stat")

fer_1M_nested <- fer_1M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_1) %>% 
  mutate(result = map(data, ~run_nb_1(dat.x = .x, dat.y = dat.rna.1M, sequencing_depth = sequencing_depth_1M))) %>% 
  collect() %>% 
  as_tibble()

fer_1M_nested$group <- "1 month"

print("1 month group analysis is finished")

# 2M group
xiang_cluster_2 <- new_cluster(n = 22) %>% 
  cluster_library("edgeR") %>% 
  cluster_library("limma") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("purrr") %>% 
  cluster_library("tibble") %>%
  cluster_copy("build_mod_2") %>% 
  cluster_copy("run_nb_2") %>%
  cluster_copy("dat.rna.2M") %>% 
  cluster_copy("ensembl_ncbi_final") %>% 
  cluster_copy("sequencing_depth_2M") %>% 
  cluster_copy("calc_gene_stat") 

fer_2M_nested <- fer_2M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_2) %>% 
  mutate(result = map(data, ~run_nb_2(dat.x = .x, dat.y = dat.rna.2M, sequencing_depth = sequencing_depth_2M))) %>% 
  collect() %>% 
  as_tibble()

fer_2M_nested$group <- "2 months"

print("2 month group analysis is finished")

# 3M group
xiang_cluster_3 <- new_cluster(n = 22) %>% 
  cluster_library("edgeR") %>% 
  cluster_library("limma") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("purrr") %>% 
  cluster_library("tibble") %>% 
  cluster_copy("build_mod_1") %>% 
  cluster_copy("run_nb_1") %>% 
  cluster_copy("dat.rna.3M") %>% 
  cluster_copy("ensembl_ncbi_final") %>% 
  cluster_copy("sequencing_depth_3M") %>% 
  cluster_copy("calc_gene_stat")

fer_3M_nested <- fer_3M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_3) %>% 
  mutate(result = map(data, ~run_nb_1(dat.x = .x, dat.y = dat.rna.3M, sequencing_depth = sequencing_depth_3M))) %>% 
  collect() %>% 
  as_tibble()

fer_3M_nested$group <- "3 months"

print("3 month group analysis is finished")

# 6M group
xiang_cluster_4 <- new_cluster(n = 22) %>% 
  cluster_library("edgeR") %>% 
  cluster_library("limma") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("purrr") %>% 
  cluster_library("tibble") %>% 
  cluster_copy("build_mod_3") %>% 
  cluster_copy("run_nb_3") %>% 
  cluster_copy("dat.rna.6M") %>% 
  cluster_copy("ensembl_ncbi_final") %>% 
  cluster_copy("sequencing_depth_6M") %>% 
  cluster_copy("calc_gene_stat")

fer_6M_nested <- fer_6M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_4) %>% 
  mutate(result = map(data, ~run_nb_3(dat.x = .x, dat.y = dat.rna.6M, sequencing_depth = sequencing_depth_6M))) %>% 
  collect() %>% 
  as_tibble()

fer_6M_nested$group <- "6 months"

print("6 month group analysis is finished")

# combine all the gene-level statistics
fer.nested <- fer_1M_nested %>% 
  bind_rows(fer_2M_nested) %>% 
  bind_rows(fer_3M_nested) %>% 
  bind_rows(fer_6M_nested)

write_rds(fer.nested, "association_fer_transcriptomics_gene.rds")


