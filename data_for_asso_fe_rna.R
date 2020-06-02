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
         entrez_id = `EntrezGene ID`)

#
# To perform GSEA, Entrez ids are required
# The RNA data used Ensembl id, I have to convert into Entrez IDs
#

ori_gene_ids <- tibble(ensembl_id = row.names(dat.rna))
gene_id_mapping <- ori_gene_ids %>% 
  left_join(ensembl_ncbi, by = "ensembl_id")

# There are 133 ensembl ids having multiple entrez ids

gene_id_mapping %>% 
  count(ensembl_id) %>% 
  filter(n > 1) %>% 
  arrange(desc(n))

# After manual check, I found that the Entrez id should be the one 
# corresponding to the Ensembl id
# Therefore, I have to remove the redundant rows in the ensembl_ncbi table

ensembl_ncbi_new <- gene_id_mapping %>% 
  group_by(ensembl_id, gene_symbol) %>% 
  summarize(entrez_id = min(entrez_id))

identical(ensembl_ncbi_new$ensembl_id, row.names(dat.rna))

row.names(dat.rna) <- ensembl_ncbi_new$entrez_id
# if a Ensembl id has no Entrez id, then the Ensembl id will be used as the final
# gene id
gene_without_entrez <- which(is.na(row.names(dat.rna)))
row.names(dat.rna)[gene_without_entrez] <- ensembl_ncbi_new$ensembl_id[gene_without_entrez]

ensembl_ncbi_final <- ensembl_ncbi_new %>% 
  ungroup() %>% 
  mutate(gene = row.names(dat.rna))

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
