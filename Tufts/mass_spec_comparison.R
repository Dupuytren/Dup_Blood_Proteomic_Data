library(tidyverse)
library(readxl)
library(broom)
library(gee)
library("gridExtra")

work_dir <- "C:/Users/agurinovich/Box/Sebastiani_BERD/charles_eaton/"

somascan <- read_csv(paste0(work_dir, "analysis/analysis2_overall.csv"))
somascan <- somascan %>%
  filter(p.value.adj < 0.25) %>%
  pull(CohortDup)
somascan <- str_split_i(somascan, pattern = "_", 2)
str(somascan)

mass_spec <- read_csv(paste0(work_dir, "original_data/MacCoss Mass Spec Proteomic Data/Dupytren-protein_group-dec2021.csv"))
mass_spec <- mass_spec %>%
  pull(Protein)
mass_spec <- str_split_i(mass_spec, pattern = "\\|", 2)
str(mass_spec)

proteins <- somascan[which(somascan %in% mass_spec)]
# Out of 51 Proteins that had a 25% FDR significant association with the disease: 6 are in the Mass Spectrometry data: "P02787" "Q15063" "P01042" "Q16853" "Q8NBP7" "P01008"

mass_spec <- read_csv(paste0(work_dir, "original_data/MacCoss Mass Spec Proteomic Data/Dupytren-protein_group-dec2021.csv"))
mass_spec$Protein_id <- mass_spec %>%
  pull(Protein) %>%
  str_split_i(pattern = "\\|", 2)
mass_spec <- mass_spec %>%
  select(Protein_id, everything(), -Peptide, -Protein) %>%
  filter(Protein_id %in% proteins)
mass_spec <- data.frame(mass_spec)
row.names(mass_spec) <- mass_spec$Protein_id
mass_spec <- mass_spec[,-1]
mass_spec <- data.frame(t(mass_spec))
mass_spec$Subject_id <- row.names(mass_spec)
mass_spec <- as_tibble(mass_spec)
mass_spec <- mass_spec %>%
  select(Subject_id, everything())
mass_spec$Subject_id <- str_replace(str_remove(mass_spec$Subject_id, pattern = "_.*"), pattern = "\\.", replacement = "-")
mass_spec

#across samples: replace 0s with NAs
mass_spec[mass_spec == 0] <- NA

# histograms of all proteins

p1 <- ggplot(mass_spec, aes(x=P02787))+
  geom_histogram()
p2 <- ggplot(mass_spec, aes(x=Q15063))+
  geom_histogram()
p3 <- ggplot(mass_spec, aes(x=P01042))+
  geom_histogram()
p4 <- ggplot(mass_spec, aes(x=Q16853))+
  geom_histogram()
p5 <- ggplot(mass_spec, aes(x=Q8NBP7))+
  geom_histogram()
p6 <- ggplot(mass_spec, aes(x=P01008))+
  geom_histogram()

grid.arrange(p1,p2,p3,p4,p5,p6)

#mass_spec[, 2:ncol(mass_spec)] <- log(mass_spec[, 2:ncol(mass_spec)])

pheno.data <- read_excel(paste0(work_dir, "original_data/Rosetta_Stone_080923.xlsx"))
pheno.data <- pheno.data %>%
  select(`MacCoss ID Number`, `Subject ID`, Draw, Cohort, `Age at draw`, Gender) %>%
  rename(Subject_id = `MacCoss ID Number`, Redraw_id = `Subject ID`, Age = `Age at draw`)
pheno.data

pheno.data <- pheno.data %>% 
  inner_join(mass_spec, by = "Subject_id")
pheno.data

# sort the data!
pheno.data <- pheno.data %>%
  arrange(Redraw_id)
table(pheno.data$Cohort)

pheno.data$Sample_id <- pheno.data$Subject_id

pheno.data$Subject_id <- round(rank(pheno.data$Redraw_id))

pheno.data <- pheno.data %>%
  select(Sample_id, Redraw_id, Subject_id, everything())

model_list <- list()

#i <- 1
for (i in 1:length(proteins)) {
  p <- proteins[i]
  model1 <- gee(formula = paste0(p, " ~ Cohort + Age + Gender"), data = pheno.data, id = Subject_id)
  summary(model1)
  model2 <- tidy(model1)
  names(model2) <- c("term", "Estimate", "Naive_S.E.", "Naive_z", "Robust_S.E.", "Robust_z")
  model2 <-  filter(model2, term == "CohortDup")
  names(model2)[1] <- "CohortDup"
  model2[1,1] <- p
  model2$P_value <- (1-pnorm(abs(model2$Robust_z)))*2
  model_list[[i]] <- model2
}

model1 <- bind_rows(model_list)
rm(model_list)

model1 <- model1 %>%
  arrange(P_value) %>%
  mutate(p.value.adj = p.adjust(P_value, "BH")) 


colnames(model1)[2:ncol(model1)] <- paste0(names(model1[2:ncol(model1)]), "_MassSpec")
model1

model2 <- read_csv(paste0(work_dir, "analysis/analysis2_overall.csv"))
model2 <- model2 %>%
  filter(p.value.adj < 0.25)
colnames(model2)[2:ncol(model2)] <- paste0(names(model2[2:ncol(model2)]), "_Somascan")
model2$CohortDup <- str_split_i(model2$CohortDup, pattern = "_", 2)
model2

model1 <- model1 %>%
  inner_join(model2, by = "CohortDup")

#log_b(a) = log_d(a)/log_d(b)
model1$Estimate_MassSpec <- log(2)*model1$Estimate_MassSpec 
model1

protein_ids <- read_csv(paste0(work_dir, "analysis/Table_51_soma_proteins.csv"))
protein_ids <- protein_ids %>%
  select(Aptamer_id, Uniprot, Protein_name)
protein_ids

model1 <- model1 %>%
  inner_join(protein_ids, by = c("CohortDup" = "Uniprot")) %>%
  rename(Uniprot = CohortDup) %>%
  select(Aptamer_id, Uniprot, Protein_name, everything())
model1

write_csv(model1, paste0(work_dir, "analysis/MassSpec_overlap.csv"))

model1 <- read_csv(paste0(work_dir, "analysis/MassSpec_overlap.csv"))
model1

jpeg(paste0(work_dir, "analysis/MassSpec_overlap.jpg"))
ggplot(model1, aes(x = Estimate_MassSpec, y = Estimate_Somascan)) +
  geom_point() +
  xlab("Ln(beta_MassSpec)") +
  ylab("Ln(beta_Somascan)") +
  theme(text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )
dev.off()

