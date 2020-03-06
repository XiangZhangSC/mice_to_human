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


run_edgeR <- function(dat.x, dat.y) {
  x <- ( dat.x$fer - mean(dat.x$fer) ) / sd(dat.x$fer)
  Batch <- factor(dat.x$Reagent)
  zAge <- (dat.x$age_day - mean(dat.x$age_day)) / sd(dat.x$age_day)
  mod1 <- model.matrix(~Batch + zAge + x)
  my_dispersions <- estimateDisp(dat.y, design = mod1)
  my_fits <- glmFit(dat.y, design = mod1, dispersion = my_dispersions$tagwise.dispersion, offset = log(sequencing_depth$N))
  my_tests <- glmLRT(my_fits, coef = "x")
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

gene.annotation <- read_xlsx("1605_Paalvast_Lexogen.expression.genelevel.v82.htseq.xlsx", range = "A1:B22898")
gene.annotation <- gene.annotation %>% 
  rename(gene = probe, 
         symbol = `gene name`)

rna <- read_tsv("1605_Paalvast_Lexogen.expression.genelevel.v82.htseq.txt.table")

sequencing_depth <- rna %>% 
  gather(mouse, y, -probe, -`gene name`) %>% 
  group_by(mouse) %>% 
  summarise(N = sum(y)) %>% 
  arrange(mouse)

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

mouse_rna.df <- mouse.df %>% 
  semi_join(sequencing_depth, by = "mouse") %>% 
  left_join(rna.batch, by = "mouse") %>% 
  mutate(M_sacrifice = as.character(M_sacrifice), 
         Batch = as.character(Batch))


# Preprocessing

nonzero_counts <- rna %>% 
  gather(mouse, y, -probe, -`gene name`) %>% 
  group_by(probe) %>% 
  summarize(n_nonzero = sum(y != 0))

keep_genes <- nonzero_counts %>% 
  filter(n_nonzero >= 43)


mingled.mcmc <- read_stan_csv("~/cmdstan-2.19.1/mingled/output.csv")

tidy(mingled.mcmc, pars = c("fer_max", "bw_max"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

tidy(mingled.mcmc, pars = c("mu_fi", "sigma_fi"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

tidy(mingled.mcmc, pars = c("mu_bw_max", "sigma_bw_max"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

tidy(mingled.mcmc, pars = c("fer_max_a", "fer_max_b"), estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

post <- tidy_draws(mingled.mcmc)

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

#ggsave("fe_sacrifice.pdf", height = 200, width = 220, units = "mm", dpi = 300)

dat.rna <- rna %>% 
  semi_join(keep_genes, by = "probe") %>% 
  rename(gene = `gene name`) %>% 
  column_to_rownames("probe") %>% 
  select(-gene) %>% 
  as.matrix()

dat.rna <- dat.rna[,sort(unique(fer$mouse))]

fer.nested <- fer %>% 
  group_by(.iteration) %>% 
  nest()

identical(sequencing_depth$mouse, colnames(dat.rna))

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

#write_rds(fer.nested, "association_fer_transcriptomics.rds")
fer.nested <- read_rds("association_fer_transcriptomics.rds")

fer2rna <- fer.nested %>% 
  select(-data) %>% 
  unnest(result)

fer2rna.sig <- fer2rna %>% 
  group_by(symbol) %>% 
  summarise(`Median beta` = median(estimate), 
            `How many times FDR is below 0.05` = sum(FDR < 0.05)) %>% 
  arrange(desc(`How many times FDR is below 0.05`))

fer2rna.sig

write_rds(fer2rna.sig, "result_genes_associated_FE.rds")

fer2rna %>% 
  filter(symbol %in% "Txnip") %>% 
  ggplot(aes(estimate)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  labs(title = "Txnip", 
       x = "Mean difference in expression due to 1 standard deviation increase in feed efficiency", 
       y = NULL) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 15), 
        axis.text.x = element_text(size = 12), 
        plot.title = element_text(size = 30))

phenotype.paalvast <- read_rds("Paalvast2017.rds")
liver.weight <- phenotype.paalvast %>% 
  filter(what == "liver_g") %>% 
  semi_join(mouse_rna.df, by = "mouse") %>% 
  select(mouse, y) %>% 
  rename(liver_g = y) %>% 
  arrange(mouse)

dat.x <- fer.nested$data[[1]]
x <- (dat.x$fer - 0.1) / 0.01
Batch <- factor(dat.x$Reagent)
mod1 <- model.matrix(~Batch + x)

Y <- voom(dat.rna, design = mod1, normalize.method = "none")
myfit <- lmFit(Y, design = mod1, weights = Y$weights)
myfit <- eBayes(myfit)
myres <- tidy(myfit) %>% 
  filter(term == "x") %>% 
  select(gene, estimate, statistic, p.value) %>% 
  inner_join(gene.annotation, by = "gene")

myres$FDR <- p.adjust(myres$p.value, method = "BH")

#mod0 <- model.matrix(~Batch)
#svseq <- svaseq(dat.rna,mod1,mod0)
#modSv <- cbind(mod1,svseq$sv)
mod1 <- model.matrix(~Reagent + x)
my_dispersions <- estimateDisp(dat.rna, design = mod1)
my_fits <- glmFit(dat.rna, design = mod1, dispersion = my_dispersions$tagwise.dispersion, offset = log(colSums(dat.rna)))

my_betas <- my_fits$coefficients %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  tbl_df() %>% 
  rename(gene = rowname) %>% 
  select(gene, x)

my_tests <- glmLRT(my_fits, coef = "x")
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

my_stats$fdr <- p.adjust(my_stats$p.value, method = "BH")

my_stats %>% 
  filter(symbol == "Gm22748")

fi.post <- tidy(mingled.mcmc, pars = "fi", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")
colnames(fi.post) <- str_c("fi_", colnames(fi.post))

fe.post <- tidy(mingled.mcmc, pars = "fer_max", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")
colnames(fe.post) <- str_c("fer_", colnames(fe.post))

fife <- fi.post %>% 
  bind_cols(fe.post)

fife_rna <- fife %>% 
  mutate(mouse_id = str_sub(fi_term, 4, -2)) %>% 
  semi_join(mouse_rna.df, by = "mouse_id")

ggplot(fife, aes(x = fer_estimate, y = fi_estimate)) + 
  geom_point(shape = 1, size = 3) + 
  geom_point(data = fife_rna , size = 3, color = "blue") + 
  geom_errorbarh(aes(xmin = fer_conf.low, xmax = fer_conf.high), alpha = 0.3) + 
  geom_errorbar(aes(ymin = fi_conf.low, ymax = fi_conf.high), alpha = 0.3) + 
  labs(x = "Posterior mean of maximum feed efficiency (g BW gain/g Food)", 
       y = "Posterior mean of food intake (g/day)") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

tidy(mingled.mcmc, pars = c("fer_max", "bw_max"))
