library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(tidybayes)
library(tibble)
library(purrr)
library(multidplyr)
library(broom)
library(stringr)
library(forcats)
library(ggrepel)

###########################################
# 1. Import data
###########################################

# load gene expression data
gene.annotation <- read_xlsx("1605_Paalvast_Lexogen.expression.genelevel.v82.htseq.xlsx", range = "A1:B22898")
gene.annotation <- gene.annotation %>% 
  rename(gene = probe, 
         symbol = `gene name`)

rna <- read_tsv("1605_Paalvast_Lexogen.expression.genelevel.v82.htseq.txt.table")

rna.batch1 <- read_xlsx("batch.xlsx")
rna.batch2 <- read_tsv("keys.txt")

rna.batch2 <- rna.batch2 %>% 
  filter(!is.na(`Liver Weight`)) %>% 
  rename(mouse = Samples) %>% 
  select(mouse, Batch)

rna.batch <- rna.batch1 %>% 
  rename(R1 = `batch 1'`, 
         R2 = `batch 2'`) %>% 
  gather(Reagent, mouse) %>% 
  left_join(rna.batch2, by = "mouse")

# load mouse birth and sacrifice data
paalvast_mice <- read_rds("Paalvast_lifespan.rds")
paalvast_mice <- paalvast_mice %>% 
  mutate(day = interval(`Start HFD`, Sacrifice) / ddays(1), 
         Y_sacrifice = year(Sacrifice), 
         M_sacrifice = month(Sacrifice), 
         age_day = interval(`Date of Birth`, Sacrifice) / ddays(1)) %>% 
  arrange(mouse)

mouse.df <- tibble(mouse_id = as.character(seq(from = 1, to = 105)), 
                   mouse = paalvast_mice$mouse, 
                   group = paalvast_mice$group, 
                   Y_sacrifice = as.character(paalvast_mice$Y_sacrifice), 
                   M_sacrifice = as.character(paalvast_mice$M_sacrifice), 
                   day = paalvast_mice$day, 
                   age_day = paalvast_mice$age_day)

# total counts for each RNA samples
sequencing_depth <- rna %>% 
  gather(mouse, y, -probe, -`gene name`) %>% 
  filter(mouse %in% mouse.df$mouse) %>% 
  group_by(mouse) %>% 
  summarise(N = sum(y)) %>% 
  arrange(mouse)

# merge mouse data with gene expression data
mouse_rna.df <- mouse.df %>% 
  semi_join(sequencing_depth, by = "mouse") %>% 
  left_join(rna.batch, by = "mouse") %>% 
  mutate(M_sacrifice = as.character(M_sacrifice), 
         Batch = as.character(Batch))

#########################################################################
# 2. Preprocessing
#########################################################################

# load feed efficiency data
#mingled.mcmc <- read_stan_csv("/home/xiang/stan/cmdstan-2.19.1/MINGLeD2/output.csv")
mingled.mcmc <- read_stan_csv("/home/xiang/cmdstan-2.19.1/mingled/output.csv")

tidy(mingled.mcmc, pars = c("fer_max", "bw_max"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

tidy(mingled.mcmc, pars = c("mu_fi", "sigma_fi"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

post <- tidy_draws(mingled.mcmc)

fer <- post %>% 
  select(.iteration, fer_sacrifice.1:fer_sacrifice.105) %>% 
  gather(what, fer, -.iteration) %>% 
  separate(what, into = c("what", "mouse_id"), sep = "\\.") %>% 
  select(-what) %>% 
  inner_join(mouse_rna.df, by = "mouse_id") %>% 
  mutate(group = factor(group, levels = c("1M", "2M", "3M", "6MA", "6MC"))) %>% 
  arrange(group)

# Because RNA isolation was not well randomized
# the analysis will be performed within each cohort
removed_samples <- tibble(
  group = c("1M", "3M", "6MA", "6MC"), 
  Reagent = c("R1", "R1", "R1", "R2"), 
  Batch = c("2", "1", "1", "2")
)

mouse_rna.df %>% 
  count(group, Reagent, Batch) %>% 
  ggplot(aes(Batch, Reagent)) + 
  geom_point(aes(size = n)) + 
  theme_bw() + 
  facet_wrap(~group)

mouse_rna.df %>% 
  anti_join(removed_samples, by = c("group", "Reagent", "Batch")) %>% 
  count(group, Reagent, Batch) %>% 
  ggplot(aes(Batch, Reagent)) + 
  geom_point(aes(size = n)) + 
  theme_bw() + 
  facet_wrap(~group)
    
# In 1M group, 
# if a RNA sample was isolated on day 1, reagent 1 was used
# if a RNA sample was isolated on day 3, reagent 2 was used 
mouse_rna_1M <- mouse_rna.df %>% 
  filter(group == "1M", Batch != 2)

fer_1M <- fer %>% 
  semi_join(mouse_rna_1M, by = "mouse")

# 2M group RNA samples were isolated by reagent 2 only
# the RNA samples were isolated on day 2 or day 3
mouse_rna_2M <- mouse_rna.df %>% 
  filter(group == "2M", Batch != 1)

fer_2M <- fer %>% 
  semi_join(mouse_rna_2M, by = "mouse")

# In 3M, 6MA and 6MC group,
# if a RNA sample was isolated on day 2, reagent 1 was used
# if a RNA sample was isolated on day 3, reagent 2 was used
mouse_rna_3M <- mouse_rna.df %>% 
  filter(group == "3M", Batch != 1)

fer_3M <- fer %>% 
  semi_join(mouse_rna_3M, by = "mouse")

# 6MA and 6MC groups will be combined
# as 6M group
mouse_rna_6MA <- mouse_rna.df %>% 
  filter(group == "6MA", Batch != 1)

mouse_rna_6MC <- mouse_rna.df %>% 
  filter(group == "6MC", (Reagent == "R1" & Batch == "2") | (Reagent == "R2" & Batch == "3"))

mouse_rna_6M <- mouse_rna_6MA %>% 
  bind_rows(mouse_rna_6MC)

fer_6M <- fer %>% 
  semi_join(mouse_rna_6M, by = "mouse")

# combine feed efficiency data of different groups
# this data set is used for Figure 2
dat_fer <- fer_1M %>% 
  bind_rows(fer_2M) %>% 
  bind_rows(fer_3M) %>% 
  bind_rows(fer_6M) %>% 
  mutate(group = factor(group, levels = c("1M", "2M", "3M", "6MA", "6MC"), 
                        labels = c("1 month", "2 months", '3 months', "6 months", "6 months")))

# Figure 2: Predicted feed efficiency
ggplot(dat_fer, aes(fct_inorder(mouse), fer)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Mouse", y = "Feed Efficiency (g body weight gain/g food intake)") + 
  theme_bw() + 
  theme(legend.position = c(0.9,0.9), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 90))

ggsave("fe_sacrifice.png", height = 200, width = 220, units = "mm", dpi = 300)

# genes with too many zero counts will be excluded
nonzero_counts <- rna %>% 
  gather(mouse, y, -probe, -`gene name`) %>% 
  filter(mouse %in% dat_fer$mouse) %>% 
  group_by(probe) %>% 
  summarize(n_nonzero = sum(y != 0))

keep_genes <- nonzero_counts %>% 
  filter(n_nonzero >= 40)

# gene expression data filtering
# genes with too many zero counts
dat.rna <- rna %>% 
  semi_join(keep_genes, by = "probe") %>% 
  rename(gene = `gene name`) %>% 
  column_to_rownames("probe") %>% 
  select(-gene) %>% 
  as.matrix()

ensembl_ncbi <- read_xlsx("1605_Paalvast_Lexogen.expression.genelevel.v82.htseq.xlsx", sheet = "keys")
ensembl_ncbi <- ensembl_ncbi %>% 
  select(`Ensembl Gene ID`, `Associated Gene Name`, `EntrezGene ID`) %>% 
  rename(ensembl_id = `Ensembl Gene ID`, 
         gene_symbol = `Associated Gene Name`, 
         entrez_id = `EntrezGene ID`) %>% 
  filter(ensembl_id %in% row.names(dat.rna)) %>% 
  arrange(ensembl_id)

#
# To perform GSEA, Entrez ids are required
#
for (i in 1:nrow(dat.rna)) {
  ensembl_id <- row.names(dat.rna)[i]
  # If an Ensembl id has no Entrez id, the Ensembl id will be kept as the row name 
  if (sum(ensembl_id %in% ensembl_ncbi$ensembl_id) > 0 & sum(is.na(ensembl_ncbi$entrez_id[which(ensembl_ncbi$ensembl_id == ensembl_id)]) == 0)) {
    # If an Ensembl id has multiple Entrez ids, the first Entrez id will be used
    if (length(ensembl_ncbi$entrez_id[which(ensembl_ncbi$ensembl_id == ensembl_id)]) > 1) {
      row.names(dat.rna)[i] <- ensembl_ncbi$entrez_id[which(ensembl_ncbi$ensembl_id == ensembl_id)][1]
    } else {
      row.names(dat.rna)[i] <- ensembl_ncbi$entrez_id[which(ensembl_ncbi$ensembl_id == ensembl_id)] 
    }
  }
  ensembl_ncbi$gene[i] <- row.names(dat.rna)[i]
}

#
# Subset the gene expression data corresponding to different groups
#

subset_rna_dat <- function(rna_dat, fer_dat) {
  rna_dat_subgroup <- rna_dat[,colnames(rna_dat) %in% fer_dat$mouse]
  rna_dat_subgroup <- rna_dat_subgroup[,sort(colnames(rna_dat_subgroup))]
  return(rna_dat_subgroup)
} 

dat.rna.1M <- subset_rna_dat(dat.rna, fer_1M)
dat.rna.2M <- subset_rna_dat(dat.rna, fer_2M)
dat.rna.3M <- subset_rna_dat(dat.rna, fer_3M)
dat.rna.6M <- subset_rna_dat(dat.rna, fer_6M)

#
# Nest the feed efficiency data
#

nest_fer_dat <- function(dat.fer) {
  dat.fer.nested <- dat.fer %>% 
    arrange(mouse) %>% 
    group_by(.iteration) %>% 
    nest()
  
  return(dat.fer.nested)
}

fer_1M_nested <- nest_fer_dat(fer_1M)
fer_2M_nested <- nest_fer_dat(fer_2M)
fer_3M_nested <- nest_fer_dat(fer_3M)
fer_6M_nested <- nest_fer_dat(fer_6M)

sequencing_depth_1M <- sequencing_depth %>% semi_join(mouse_rna_1M, by = "mouse")
sequencing_depth_2M <- sequencing_depth %>% semi_join(mouse_rna_2M, by = "mouse")
sequencing_depth_3M <- sequencing_depth %>% semi_join(mouse_rna_3M, by = "mouse")
sequencing_depth_6M <- sequencing_depth %>% semi_join(mouse_rna_6M, by = "mouse")

# Check if the sample ids are matched
is_sampleid_matched <- function(sequencing_depth_dat, fer_dat, rna_dat) {
  if ((identical(sequencing_depth_dat$mouse, colnames(rna_dat))) & (identical(fer_dat$data[[1]]$mouse, sequencing_depth_dat$mouse))) {
    print("Sample IDs are matched")
  } else {
    print("Sample IDs are not matched!")
  }
}

is_sampleid_matched(sequencing_depth_1M, fer_1M_nested, dat.rna.1M)
is_sampleid_matched(sequencing_depth_2M, fer_2M_nested, dat.rna.2M)
is_sampleid_matched(sequencing_depth_3M, fer_3M_nested, dat.rna.3M)
is_sampleid_matched(sequencing_depth_6M, fer_6M_nested, dat.rna.6M)

#################################################################
# 3. Run negative binomial regression analysis
#################################################################

#
# 3.1 Negative binomial regression models
#

# negative binomial regression model
calc_gene_stat <- function(dat.y, mod, sequencing_depth) {
  my_dispersions <- estimateDisp(dat.y, design = mod)
  my_fits <- glmFit(dat.y, design = mod, dispersion = my_dispersions$tagwise.dispersion, offset = log(sequencing_depth$N))
  my_tests <- glmLRT(my_fits, coef = "zFER")
  my_stats <- my_tests$table %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(gene = rowname, 
           p.value = PValue) %>% 
    tbl_df() %>% 
    select(gene, logFC, p.value) %>% 
    left_join(ensembl_ncbi, by = "gene") %>% 
    arrange(gene)
  
  my_stats$FDR <- p.adjust(my_stats$p.value, method = "BH")
  
  return(my_stats)
}

calc_path_stat <- function(dat.y, mod, sequencing_depth) {
  my_dispersions <- estimateDisp(dat.y, design = mod)
  my_fits <- glmFit(dat.y, design = mod, dispersion = my_dispersions$tagwise.dispersion, offset = log(sequencing_depth$N))
  my_tests <- glmLRT(my_fits, coef = "zFER")
  # Pathway analysis (pathways with at most 30 genes)
  my_gsea <- kegga(my_tests, species = 'Mm')
  my_gsea_stats <- topKEGG(my_gsea, num = Inf, truncate.path = 30)
  
  return(my_gsea_stats)
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

run_path_1 <- function(dat.x, dat.y, sequencing_depth) {
  mod <- build_mod_1(dat.x)
  
  my_stats <- calc_path_stat(dat.y, mod, sequencing_depth)
  
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

run_path_2 <- function(dat.x, dat.y, sequencing_depth) {
  mod <- build_mod_2(dat.x)
  
  my_stats <- calc_path_stat(dat.y, mod, sequencing_depth)
  
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

run_path_3 <- function(dat.x, dat.y, sequencing_depth) {
  mod <- build_mod_3(dat.x)
  
  my_stats <- calc_path_stat(dat.y, mod, sequencing_depth)
  
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
  cluster_copy("run_path_1") %>% 
  cluster_copy("dat.rna.1M") %>% 
  cluster_copy("ensembl_ncbi") %>% 
  cluster_copy("sequencing_depth_1M") %>% 
  cluster_copy("calc_gene_stat") %>% 
  cluster_copy("calc_path_stat")

fer_1M_nested <- fer_1M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_1) %>% 
  mutate(result_gene = map(data, ~run_nb_1(dat.x = .x, dat.y = dat.rna.1M, sequencing_depth = sequencing_depth_1M)),
         result_pathway = map(data, ~run_path_1(dat.x = .x, dat.y = dat.rna.1M, sequencing_depth = sequencing_depth_1M))) %>% 
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
  cluster_copy("run_path_2") %>% 
  cluster_copy("dat.rna.2M") %>% 
  cluster_copy("ensembl_ncbi") %>% 
  cluster_copy("sequencing_depth_2M") %>% 
  cluster_copy("calc_gene_stat") %>% 
  cluster_copy("calc_path_stat") 

fer_2M_nested <- fer_2M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_2) %>% 
  mutate(result_gene = map(data, ~run_nb_2(dat.x = .x, dat.y = dat.rna.2M, sequencing_depth = sequencing_depth_2M)), 
         result_pathway = map(data, ~run_path_2(dat.x = .x, dat.y = dat.rna.2M, sequencing_depth = sequencing_depth_2M))) %>% 
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
  cluster_copy("run_path_1") %>% 
  cluster_copy("dat.rna.3M") %>% 
  cluster_copy("ensembl_ncbi") %>% 
  cluster_copy("sequencing_depth_3M") %>% 
  cluster_copy("calc_gene_stat") %>% 
  cluster_copy("calc_path_stat")

fer_3M_nested <- fer_3M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_3) %>% 
  mutate(result_gene = map(data, ~run_nb_1(dat.x = .x, dat.y = dat.rna.3M, sequencing_depth = sequencing_depth_3M)), 
         result_pathway = map(data, ~run_path_1(dat.x = .x, dat.y = dat.rna.3M, sequencing_depth = sequencing_depth_3M))) %>% 
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
  cluster_copy("run_path_3") %>% 
  cluster_copy("dat.rna.6M") %>% 
  cluster_copy("ensembl_ncbi") %>% 
  cluster_copy("sequencing_depth_6M") %>% 
  cluster_copy("calc_gene_stat") %>% 
  cluster_copy("calc_path_stat")

fer_6M_nested <- fer_6M_nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster_4) %>% 
  mutate(result_gene = map(data, ~run_nb_3(dat.x = .x, dat.y = dat.rna.6M, sequencing_depth = sequencing_depth_6M)), 
         result_pathway = map(data, ~run_path_3(dat.x = .x, dat.y = dat.rna.6M, sequencing_depth = sequencing_depth_6M))) %>% 
  collect() %>% 
  as_tibble()

fer_6M_nested$group <- "6 months"

print("6 month group analysis is finished")

# combine all the gene-level statistics
fer.nested <- fer_1M_nested %>% 
  bind_rows(fer_2M_nested) %>% 
  bind_rows(fer_3M_nested) %>% 
  bind_rows(fer_6M_nested)

write_rds(fer.nested, "association_fer_transcriptomics.rds")
#fer.nested <- read_rds("association_fer_transcriptomics.rds")

################################################################################
# 4. Visualizing associations between feed efficiency and liver gene expression
################################################################################

#fer2rna <- fer.nested %>% 
#  select(-data) %>% 
#  unnest(result)

#fer2rna.sig <- fer2rna %>% 
#  group_by(symbol, group) %>% 
#  summarise(`Mean logFC` = mean(logFC), 
#            `Low logFC` = quantile(logFC, probs = 0.025), 
#            `High logFC` = quantile(logFC, probs = 0.975), 
#            `Mean pval` = mean(p.value), 
#            `Low pval` = quantile(p.value, probs = 0.025), 
#            `High pval` = quantile(p.value, probs = 0.975), 
#            `Mean FDR` = mean(FDR), 
#            `Low FDR` = quantile(FDR, probs = 0.025), 
#            `High FDR` = quantile(FDR, probs = 0.975), 
#            `How many times FDR is below 0.05` = sum(FDR < 0.05))

#write_csv(fer2rna.sig, "Association_feed_efficiency_liver_gene_expression.csv")

#fer2rna.sig2 <- fer2rna.sig %>% 
#  filter(`How many times FDR is below 0.05` >= 950 & abs(`Mean logFC`) >= 0.6)

# volcano plot
#ggplot(fer2rna.sig, aes(`Mean logFC`, -log10(`Mean pval`))) + 
#  geom_point(alpha = 0.1) + 
#  geom_point(data = fer2rna.sig2, color = "red") + 
#  geom_text_repel(data = fer2rna.sig2, aes(label = symbol)) + 
#  labs(x = "Mean log(Fold change)", y = "-log10(Mean P value)") + 
#  facet_wrap(~group) + 
#  theme_bw() + 
#  theme(axis.title = element_text(size = 12), 
#        strip.text = element_text(size = 12))

ggsave("asso_fer_liver_gene.pdf", height = 200, width = 200, units = "mm", dpi = 300)


