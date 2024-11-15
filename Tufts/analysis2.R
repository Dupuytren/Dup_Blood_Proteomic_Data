# Load packages
library(tidyverse)
library(readxl)
library(tidyr)
library(broom)
library(gee)

work_dir <- "C:/Users/agurinovich/Box/Sebastiani_BERD/charles_eaton/"

prot.data <- read_csv(paste0(work_dir,  "analysis/SS-2216784_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.QCed.human.proteins.csv"))

# To log-transform the data:
prot.data[2:ncol(prot.data)] <- log(prot.data[2:ncol(prot.data)])
prot.data

prot.data.df <- data.frame(prot.data)
row.names(prot.data.df) <- prot.data$Assay
prot.data.df <- prot.data.df[,-1]
data.intensity <- data.frame(t(prot.data.df))
dim(data.intensity)
data.intensity$`Somalogic ID Number` <- row.names(data.intensity)
data.intensity <- as_tibble(data.intensity) %>%
  select(`Somalogic ID Number`, everything()) %>%
  rename(Sample_id = `Somalogic ID Number`)
data.intensity

pheno.data <- read_excel(paste0(work_dir, "original_data/Rosetta_Stone_080923.xlsx"))
pheno.data

pheno.data <- pheno.data %>%
  select(`Somalogic ID Number`, `Subject ID`, Draw, Cohort, `Age at draw`, Gender) %>%
  rename(Sample_id = `Somalogic ID Number`, Redraw_id = `Subject ID`, Age = `Age at draw`)
pheno.data

# to remove samples that are not in the QCed protein data
samples_prot <-  names(prot.data)[2:ncol(prot.data)]
pheno.data <- pheno.data %>%
  filter(Sample_id %in% samples_prot)

pheno.data$Sample_id <- paste0("X", pheno.data$Sample_id)
# substitute space and parentheses on dots:
pheno.data$Sample_id <- str_replace(str_replace(str_replace(pheno.data$Sample_id, " ", "."), "\\(", "."), "\\)", ".")
pheno.data

pheno.data <- pheno.data %>% 
  inner_join(data.intensity, by = "Sample_id")
pheno.data

proteins <- names(data.intensity)[2:ncol(data.intensity)]

####### Redraw analysis

pheno.data_redraw <- pheno.data %>%
  filter(Draw == "Redraw")
table(pheno.data_redraw$Cohort)
pheno.data_redraw

model_redraw_list <- list()

#i <- 1
for (i in 1:length(proteins)) {
  p <- proteins[i]
  model1 <- lm(formula = paste0(p, " ~ Cohort + Age + Gender"), data = pheno.data_redraw)
  model1 <- tidy(model1)
  model1 <-  filter(model1, term == "CohortDup")
  names(model1)[1] <- "CohortDup"
  model1[1,1] <- p
  model_redraw_list[[i]] <- model1
}

model1 <- bind_rows(model_redraw_list)
rm(model_redraw_list)

model1 <- model1 %>%
  arrange(p.value) %>%
  mutate(p.value.adj = p.adjust(p.value, "BH")) 
model1

write_csv(model1, paste0(work_dir, "analysis/analysis2_redraw.csv"))

####### Initial analysis

pheno.data_initial <- pheno.data %>%
  filter(Draw == "Initial")
pheno.data_initial
table(pheno.data_initial$Cohort)

model_initial_list <- list()

#i <- 1
for (i in 1:length(proteins)) {
  p <- proteins[i]
  model1 <- lm(formula = paste0(p, " ~ Cohort + Age + Gender"), data = pheno.data_initial)
  model1 <- tidy(model1)
  model1 <-  filter(model1, term == "CohortDup")
  names(model1)[1] <- "CohortDup"
  model1[1,1] <- p
  model_initial_list[[i]] <- model1
}

model1 <- bind_rows(model_initial_list)
rm(model_initial_list)

model1 <- model1 %>%
  arrange(p.value) %>%
  mutate(p.value.adj = p.adjust(p.value, "BH")) 
model1

write_csv(model1, paste0(work_dir, "analysis/analysis2_initial.csv"))

####### Redraw + Initial analysis

# sort the data!
pheno.data <- pheno.data %>%
  arrange(Redraw_id)
table(pheno.data$Cohort)

pheno.data$Subject_id <- round(rank(pheno.data$Redraw_id))

pheno.data <- pheno.data %>%
  select(Sample_id, Redraw_id, Subject_id, everything())

model_list <- list()

# pheno.data <- read_csv(paste0(work_dir, "analysis/file_for_Ana.csv"))

#i <- 851
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
model1

write_csv(model1, paste0(work_dir, "analysis/analysis2_overall.csv"))

##########

model1 <- read_csv(paste0(work_dir, "analysis/analysis2_overall.csv"))
model1 <- model1 %>%
  separate_wider_delim(cols = CohortDup, delim = "_", names = c("Aptamer_id", "Uniprot", "Protein_name"))

model1 <- model1 %>%
  filter(p.value.adj < 0.25)

model1

write_csv(model1, paste0(work_dir, "analysis/Table_51_soma_proteins.csv"))


##########


library(readxl)

model1 <- read_csv(paste0(work_dir, "analysis/analysis2_overall.csv"))
model1 <- model1 %>%
  separate_wider_delim(cols = CohortDup, delim = "_", names = c("Aptamer_id", "Uniprot", "Protein_name"))
model1



prot_list <- read_excel(paste0(work_dir, "analysis/Top-down Somalogic SeqIDs.xlsx"))
prot_list <- prot_list %>%
  rename(Aptamer_id = `445 Unique SeqIDs`, Uniprot2 = `327 Unique Uniprots associated with these 445 Unique SeqIDs`)
prot_list

model1$Aptamer_id <- str_replace(str_replace(model1$Aptamer_id, "seq.", ""), "\\.", "-")
model1  

prot_list %>%
  inner_join(model1, by = "Aptamer_id")
#422 out of 444 matched (22 didn't match)

prot_list %>%
  anti_join(model1, by = "Aptamer_id") %>%
  print(n=22)

prot_list2 <- prot_list %>%
  inner_join(model1, by = "Aptamer_id") %>%
  select(-Uniprot2)
prot_list2 <- prot_list2 %>%
  arrange(P_value) %>%
  mutate(p.value.adj = p.adjust(P_value, "BH")) 
prot_list2

write_csv(prot_list2, paste0(work_dir, "analysis/Updated_protein_list_adjusted_pvalues.csv"))
