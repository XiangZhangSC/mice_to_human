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

# negative binomial regression model
run_edgeR <- function(dat.x, dat.y) {
  fer <- dat.x$fer
  
  zFER <- ( fer - mean(fer) ) / sd(fer)
  
  Batch <- factor(dat.x$Reagent)
  Group <- factor(dat.x$group, levels = c("1 month", "3 month", "6 month"))
  #zAge <- (dat.x$age_day - mean(dat.x$age_day)) / sd(dat.x$age_day)
  #mod1 <- model.matrix(~Batch + zAge + zFER)
  mod1 <- model.matrix(~Batch + zFER)
  my_dispersions <- estimateDisp(dat.y, design = mod1)
  my_fits <- glmFit(dat.y, design = mod1, dispersion = my_dispersions$tagwise.dispersion, offset = log(sequencing_depth$N))
  my_tests <- glmLRT(my_fits, coef = "zFER")
  my_stats <- my_tests$table %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(gene = rowname, 
           p.value = PValue, 
           estimate = logFC) %>% 
    tbl_df() %>% 
    select(gene, estimate, p.value) %>% 
    left_join(gene.annotation, by = "gene") %>% 
    arrange(gene)
  
  my_stats$FDR <- p.adjust(my_stats$p.value, method = "BH")
  
  return(my_stats)
}

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

# load mouse data
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
                   age_day = paalvast_mice$age_day) %>% 
  filter(group != "2M") # 2 month feeding mice were processed only Reagent 2

#########################################################################
# 2. Preprocessing
#########################################################################

# load feed efficiency data
#mingled.mcmc <- read_stan_csv("/home/xiang/stan/cmdstan-2.19.1/MINGLeD2/output.csv")
mingled.mcmc <- read_stan_csv("/home/xiang/cmdstan-2.19.1/mingled/output.csv")

tidy(mingled.mcmc, pars = c("fer_max", "bw_max"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

tidy(mingled.mcmc, pars = c("mu_fi", "sigma_fi"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

post <- tidy_draws(mingled.mcmc)

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

fer <- post %>% 
  select(.iteration, fer_sacrifice.1:fer_sacrifice.105) %>% 
  gather(what, fer, -.iteration) %>% 
  separate(what, into = c("what", "mouse_id"), sep = "\\.") %>% 
  select(-what) %>% 
  inner_join(mouse_rna.df, by = "mouse_id") %>% 
  mutate(group = factor(group, levels = c("1M", "2M", "3M", "6MA", "6MC"), 
                        labels = c("1 month", "2 months", "3 months", "6 months", "6 months"))) %>% 
  arrange(group)

ggplot(fer, aes(fct_inorder(mouse), fer)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Mouse", y = "Feed Efficiency (g BW gain/g Food)") + 
  theme_bw() + 
  theme(legend.position = c(0.9,0.9), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 90))

ggsave("fe_sacrifice.pdf", height = 200, width = 220, units = "mm", dpi = 300)

# genes with too many zero counts will be excluded
nonzero_counts <- rna %>% 
  gather(mouse, y, -probe, -`gene name`) %>% 
  filter(mouse %in% fer$mouse) %>% 
  group_by(probe) %>% 
  summarize(n_nonzero = sum(y != 0))

keep_genes <- nonzero_counts %>% 
  filter(n_nonzero >= 30)

dat.rna <- rna %>% 
  semi_join(keep_genes, by = "probe") %>% 
  rename(gene = `gene name`) %>% 
  column_to_rownames("probe") %>% 
  select(-gene) %>% 
  as.matrix()

dat.rna <- dat.rna[,colnames(dat.rna) %in% fer$mouse]
dat.rna <- dat.rna[,sort(unique(fer$mouse))]

#################################################################
# 3. Run negative binomial regression analysis
#################################################################
fer.nested <- fer %>% 
  arrange(mouse) %>% 
  group_by(.iteration) %>% 
  nest()

identical(sequencing_depth$mouse, colnames(dat.rna))
identical(fer.nested$data[[1]]$mouse, sequencing_depth$mouse)

xiang_cluster <- new_cluster(n = 22) %>% 
  cluster_library("edgeR") %>% 
  cluster_library("limma") %>% 
  cluster_library("biobroom") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("purrr") %>% 
  cluster_library("tibble") %>% 
  cluster_copy("run_edgeR") %>% 
  cluster_copy("dat.rna") %>% 
  cluster_copy("gene.annotation") %>% 
  cluster_copy("sequencing_depth")

fer.nested <- fer.nested %>% 
  group_by(.iteration) %>% 
  partition(cluster = xiang_cluster) %>% 
  mutate(result = map(data, ~run_edgeR(dat.x = .x, dat.y = dat.rna))) %>% 
  collect() %>% 
  as_tibble()

write_rds(fer.nested, "association_fer_transcriptomics.rds")
#fer.nested <- read_rds("association_fer_transcriptomics.rds")

fer2rna <- fer.nested %>% 
  select(-data) %>% 
  unnest(result)

fer2rna.sig <- fer2rna %>% 
  group_by(symbol) %>% 
  summarise(`Mean beta` = mean(estimate), 
            `Low beta` = quantile(estimate, probs = 0.025), 
            `High beta` = quantile(estimate, probs = 0.975), 
            `Mean pval` = mean(p.value), 
            `Low pval` = quantile(p.value, probs = 0.025), 
            `High pval` = quantile(p.value, probs = 0.975), 
            `Mean FDR` = mean(FDR), 
            `Low FDR` = quantile(FDR, probs = 0.025), 
            `High FDR` = quantile(FDR, probs = 0.975), 
            `How many times FDR is below 0.05` = sum(FDR < 0.05))

write_csv(fer2rna.sig, "Association_feed_efficiency_liver_gene_expression.csv")

fer2rna.sig2 <- fer2rna.sig %>% 
  filter(`How many times FDR is below 0.05` == 1000)

# volcano plot
ggplot(fer2rna.sig, aes(`Mean beta`, -log10(`Mean pval`))) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = fer2rna.sig2, color = "red") + 
  geom_text_repel(data = fer2rna.sig2, aes(label = symbol)) + 
  labs(x = "Mean log(Fold change)", y = "-log10(Mean P value)") + 
  facet_warp(~group, scale = "free_y", ncol = 3)
  theme_bw() + 
  theme(axis.title = element_text(size = 12), 
        strip.text = element_text(size = 12))

ggsave("asso_fer_liver_gene.pdf", height = 100, width = 300, units = "mm", dpi = 300)


