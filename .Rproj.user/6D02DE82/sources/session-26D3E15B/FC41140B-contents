library(readr)
library(dplyr)
library(tidyverse)
library(stringr)
library(strex)
library(tidyverse)
library(edgeR)
library(dplyr)
library(purrr)
library(data.table)
library(rlang)
library(ggplotify)
library(ggplot2)
library(wesanderson)
library(ggrepel)
library(gprofiler2)
library(kableExtra)
#install.packages("BiocManager")
library(Biostrings)
library(ggpubr)


################################################################################ 
change_esid12_anno<-function(anno_in){
  anno_in<-anno_in%>%dplyr::select(-c(ESid))
  anno_in$ESid<-anno_in$ESid2
  anno_in<-anno_in%>%dplyr::select(-c(ESid2))
  anno_in
}


################################################################################ 
gencode<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/All_Table/gencode.v38.primary_assembly.tx2gene.tsv.gz")
gencode$GENEID<-gencode$geneid
gencode$GENENAME<-gencode$genename
gencode<-gencode%>%dplyr::select(c(GENEID,GENENAME))

############### ANNOTATION
M1_anno_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_All.tsv")
M1_anno_all<-change_esid12_anno(M1_anno_all)
M1_anno_all_cohort<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_All.tsv")
M1_anno_all_cohort<-change_esid12_anno(M1_anno_all_cohort)
M1_anno_all_cohort$Cohort<-c("Discovery")
############### DEG
M1_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Gene_Differential/LPS_DEG_Toptable.tsv")
M1_DEG_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Gene_Differential/IFNb_DEG_Toptable.tsv")
############### DES
M1_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Editing_Differential/LPS_1_Toptable.tsv")
M1_DES_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Editing_Differential/IFNb_1_Toptable.tsv")


############### ANNOTATION
M2_anno_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Jacusa/Processed_Annotation_All.tsv")
M2_anno_all<-change_esid12_anno(M2_anno_all)
M2_anno_all_cohort<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Jacusa/Processed_Annotation_All.tsv")
M2_anno_all_cohort<-change_esid12_anno(M2_anno_all_cohort)
M2_anno_all_cohort$Cohort<-c("Alternation")
############### DEG
M2_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Toptable/Gene_Differential/LPS_DEG_Toptable.tsv")
M2_DEG_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Toptable/Gene_Differential/IFNg_DEG_Toptable.tsv")
############### DES
M2_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Toptable/Editing_Differential/LPS_1_Toptable.tsv")
M2_DES_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Toptable/Editing_Differential/IFNg_1_Toptable.tsv")


############### ANNOTATION
Micro_anno_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Jacusa/Processed_Annotation_All.tsv")
Micro_anno_all<-change_esid12_anno(Micro_anno_all)
Micro_anno_all_cohort<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Jacusa/Processed_Annotation_All.tsv")
Micro_anno_all_cohort<-change_esid12_anno(Micro_anno_all_cohort)
Micro_anno_all_cohort$Cohort<-c("Replication")
############### DEG
Micro_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Toptable/Gene_Differential/LPS_DEG_Toptable.tsv")
Micro_DEG_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Toptable/Gene_Differential/IFNg_DEG_Toptable.tsv")
############### DES
Micro_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Toptable/Editing_Differential/LPS_1_Toptable.tsv")
Micro_DES_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Toptable/Editing_Differential/IFNg_1_Toptable.tsv")

# save RData
################################################################################ 
################################################################################
################################################################################
all_esid<-rbind(M1_anno_all, M2_anno_all, Micro_anno_all)
all_esid<-all_esid%>%filter(!duplicated(across(-ESid)))
################################################################################

all_esid<- all_esid %>%
  mutate(Discovery_DES_LPS = ifelse(ESid %in% M1_DES_LPS$ESid, 
                                M1_DES_LPS$DE_Index[match(ESid, M1_DES_LPS$ESid)], 0),
         Discovery_DES_IFN = ifelse(ESid %in% M1_DES_IFNb$ESid, 
                                 M1_DES_IFNb$DE_Index[match(ESid, M1_DES_IFNb$ESid)], 0),
         Replication_DES_LPS = ifelse(ESid %in% M2_DES_LPS$ESid, 
                                  M2_DES_LPS$DE_Index[match(ESid, M2_DES_LPS$ESid)], 0),
         Replication_DES_IFN = ifelse(ESid %in% M2_DES_IFNg$ESid, 
                                   M2_DES_IFNg$DE_Index[match(ESid, M2_DES_IFNg$ESid)], 0),
         Alternation_DES_LPS = ifelse(ESid %in% Micro_DES_LPS$ESid, 
                                  Micro_DES_LPS$DE_Index[match(ESid, Micro_DES_LPS$ESid)], 0),
         Alternation_DES_IFN = ifelse(ESid %in% Micro_DES_IFNg$ESid, 
                                   Micro_DES_IFNg$DE_Index[match(ESid, Micro_DES_IFNg$ESid)], 0))

all_esid<-all_esid%>%dplyr::select(c(ESid, Editing_Index, Gene, 
                              Discovery_DES_LPS,Discovery_DES_IFN,
                              Replication_DES_LPS,Replication_DES_IFN,
                              Alternation_DES_LPS,Alternation_DES_IFN,
                              everything()))

# de_dataframes <- list()
# conditions <- c("M1_DES_LPS","M2_DES_LPS", "Micro_DES_LPS", 
#                 "M1_DES_IFNb", "M2_DES_IFNg", "Micro_DES_IFNg")
# 
# # Loop through conditions to create dataframes
# for (condition in conditions) {
#   df <- get(condition) %>%
#     filter(grepl("A:G_DE|C:T_DE", DE_Index)) %>%
#     select(ESid, logFC) %>%
#     rename(!!paste0("Discovery_", condition, "_logFC") := !!sym("logFC"))
#   
#   # Store dataframe in the list
#   de_dataframes[[condition]] <- df
#}

make_de_des<-function(des_in, colname_in){
  des_de<-des_in%>%
    filter(grepl("A:G_DE|C:T_DE", DE_Index)) %>%
    select(ESid, logFC)
  des_de
}


M1_DES_LPS_de <-make_de_des(M1_DES_LPS)%>%
  dplyr::rename(Discovery_DES_LPS_logFC = logFC)

M2_DES_LPS_de <-make_de_des(M2_DES_LPS)%>%
  dplyr::rename(Replication_DES_LPS_logFC = logFC)

Micro_DES_LPS_de <-make_de_des(Micro_DES_LPS)%>%
  dplyr::rename(Alternation_DES_LPS_logFC = logFC)

M1_DES_IFN_de <-make_de_des(M1_DES_IFNb)%>%
  dplyr::rename(Discovery_DES_IFN_logFC = logFC)

M2_DES_IFN_de <-make_de_des(M2_DES_IFNg)%>%
  dplyr::rename(Replication_DES_IFN_logFC = logFC)

Micro_DES_IFN_de <-make_de_des(Micro_DES_IFNg)%>%
  dplyr::rename(Alternation_DES_IFN_logFC = logFC)

# Merge dataframes into all_esid_merge
all_esid_merge <- 
  left_join(all_esid, M1_DES_LPS_de, by = "ESid") %>%
  left_join(M2_DES_LPS_de, by = "ESid") %>%
  left_join(Micro_DES_LPS_de, by = "ESid") %>%
  left_join(M1_DES_IFN_de, by = "ESid") %>%
  left_join(M2_DES_IFN_de, by = "ESid") %>%
  left_join(Micro_DES_IFN_de, by = "ESid")


all_esid_merge<-all_esid_merge%>%select(c(
  ESid, Gene, Editing_Index, Discovery_DES_LPS, Discovery_DES_LPS_logFC,
  Discovery_DES_IFN, Discovery_DES_IFN_logFC, Replication_DES_LPS,
  Replication_DES_LPS_logFC, Replication_DES_IFN, Replication_DES_IFN_logFC,
  Alternation_DES_LPS, Alternation_DES_LPS_logFC, Alternation_DES_IFN,
  Alternation_DES_IFN_logFC, Location, Mutation, known_a_i, rep_type
))

################################################################################ 
################################################################################
################################################################################

update_and_select <- function(df, pattern, suffix) {
  df %>%
    mutate(!!paste0(pattern, suffix) := DE_Direction) %>%
    select(GENEID, !!paste0(pattern, suffix), GENENAME)}

M1_DEG_LPS_de <- update_and_select(M1_DEG_LPS,"Discovery", "_DEG_LPS")
M1_DEG_IFNb_de <- update_and_select(M1_DEG_IFNb, "Discovery","_DEG_IFNb")
M2_DEG_LPS_de <- update_and_select(M2_DEG_LPS,"Replication","_DEG_LPS")
M2_DEG_IFNg_de <- update_and_select(M2_DEG_IFNg, "Replication","_DEG_IFNg")
Micro_DEG_LPS_de <- update_and_select(Micro_DEG_LPS,"Alternation", "_DEG_LPS")
Micro_DEG_IFNg_de <- update_and_select(Micro_DEG_IFNg, "Alternation","_DEG_IFNg")

data_frames <- list(M1_DEG_LPS_de, M1_DEG_IFNb_de,M2_DEG_LPS_de, 
                    M2_DEG_IFNg_de, Micro_DEG_LPS_de, Micro_DEG_IFNg_de)

################################################################################
all_esid_gene <- Reduce(function(df1, df2) {
  merge(df1, df2, by.x = "Gene", by.y = "GENENAME", all.x = TRUE)
}, data_frames, init = all_esid_merge)

all_esid_gene <- all_esid_gene %>%
  select(-starts_with("GENEID"))%>%select(-starts_with("GENENAME"))   

all_esid_gene<- merge(all_esid_gene, gencode, by.x = "Gene", by.y = "GENENAME", all.x = TRUE)
all_esid_gene<-all_esid_gene%>%distinct(ESid, .keep_all=TRUE)
all_esid_gene<-all_esid_gene%>%dplyr::select(c(ESid, Gene, Editing_Index,everything()))

all_esid_gene<-all_esid_gene%>%dplyr::filter(grepl("A:G|C:T|T:C|G:A",Editing_Index))

################################################################################
################################################################################
# # Create M1_DES_LPS_de dataframe
# make_de_deg<-function(deg_in, colname_in){
#   deg_de<-deg_in%>%
#     filter(grepl("UP|DOWN", DE_Direction)) 
#   deg_de$Gene <- deg_de$GENENAME
#   deg_de<- deg_de%>% select(Gene, logFC) 
#   deg_de
# }
# 
# M1_DEG_LPS_de <-make_de_deg(M1_DEG_LPS)%>%
#   dplyr::rename(Discovery_DEG_LPS_logFC = logFC)
# 
# M2_DEG_LPS_de <-make_de_deg(M1_DEG_LPS)%>%
#   dplyr::rename(Replication_DEG_LPS_logFC = logFC)
# 
# Micro_DEG_LPS_de <-make_de_deg(Micro_DEG_LPS)%>%
#   dplyr::rename(Alternation_DEG_LPS_logFC = logFC)
# 
# M1_DEG_IFN_de <-make_de_deg(M1_DEG_IFNb)%>%
#   dplyr::rename(Discovery_DEG_IFN_logFC = logFC)
# 
# M2_DEG_IFN_de <-make_de_deg(M2_DEG_IFNg)%>%
#   dplyr::rename(Replication_DEG_IFN_logFC = logFC)
# 
# Micro_DEG_IFN_de <-make_de_deg(Micro_DEG_IFNg)%>%
#   dplyr::rename(Alternation_DEG_IFN_logFC = logFC)


# Merge dataframes into all_esid_merge
# all_esid_merge <- 
#   left_join(all_esid_merge, M1_DEG_LPS_de, by = "Gene") %>%
#   left_join(M2_DEG_LPS_de, by = "Gene") %>%
#   left_join(Micro_DEG_LPS_de, by = "Gene") %>%
#   left_join(M1_DEG_IFN_de, by = "Gene") %>%
#   left_join(M2_DEG_IFN_de, by = "Gene") %>%
#   left_join(Micro_DEG_IFN_de, by = "Gene")
################################################################################
################################################################################

all_esid_cohort<-rbind(M1_anno_all_cohort, M2_anno_all_cohort, Micro_anno_all_cohort)
all_esid_cohort<-all_esid_cohort%>%filter(!duplicated(across(-ESid)))

all_esid_cohort_count <- all_esid_cohort %>%
  group_by(ESid) %>%
  mutate(
    Editing_Count = n(),
    Editing_Cohort = paste(unique(Cohort), collapse = ", ")
  ) %>%
  ungroup()%>%
  distinct(ESid, .keep_all=TRUE)

all_esid_cohort_count<-all_esid_cohort_count%>%
  dplyr::select(c("ESid","Editing_Count","Editing_Cohort"))

################################################################################
all_esid_gene<-merge(all_esid_gene, all_esid_cohort_count, by='ESid')
all_esid_gene<-all_esid_gene%>%dplyr::select(
  ESid, Gene, Editing_Index, Editing_Count, Editing_Cohort,
  everything())

write_tsv(all_esid_gene, "All_Table/All_DES_DEG/BOTH_DES_DEG_All.tsv")
################################################################################
lps_all_esid_gene<-all_esid_gene%>%dplyr::select(
  ESid, Gene, Editing_Index, Editing_Count, Editing_Cohort,
  Discovery_DES_LPS, Discovery_DES_LPS_logFC, 
  Replication_DES_LPS, Replication_DES_LPS_logFC, 
  Alternation_DES_LPS, Alternation_DES_LPS_logFC,
  Discovery_DEG_LPS, Replication_DEG_LPS, Alternation_DEG_LPS,
  everything())%>%
  select(-matches("IFN|IFNb|IFNg"))

write_tsv(lps_all_esid_gene, "All_Table/All_DES_DEG/LPS_DES_DEG_All.tsv")

ifn_all_esid_gene<-all_esid_gene%>%dplyr::select(
  ESid, Gene, Editing_Index, Editing_Count, Editing_Cohort,
  Discovery_DES_IFN, Discovery_DES_IFN_logFC,
  Replication_DES_IFN, Replication_DES_IFN_logFC,
  Alternation_DES_IFN, Alternation_DES_IFN_logFC,
  Discovery_DEG_IFNb, Replication_DEG_IFNg, Alternation_DEG_IFNg,
  everything())%>%
  select(-matches("LPS"))

write_tsv(ifn_all_esid_gene, "All_Table/All_DES_DEG/IFN_DES_DEG_All.tsv")

################################################################################ 
################################################################################
################################################################################
deg_des_in<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/All_Table/All_DES_DEG/LPS_DES_DEG_All.tsv")
deg_des_in_dis<-deg_des_in%>%dplyr::filter(grepl("Discovery",Editing_Cohort))
deg_des_in_dis<-deg_des_in_dis%>%dplyr::filter(!grepl("Non_DE",Discovery_DES_LPS))
deg_des_in_dis<-deg_des_in_dis%>%dplyr::filter(!grepl("FALSE",Discovery_DEG_LPS))
deg_des_in_dis$Editing_Index<-gsub("T:C","A:G",deg_des_in_dis$Editing_Index)
deg_des_in_dis$Editing_Index<-gsub("G:A","C:T",deg_des_in_dis$Editing_Index)
# > length(unique(deg_des_in_dis$ESid))
# [1] 263
# > length(unique(deg_des_in_dis$Gene))
# [1] 107
deg_des_in_dis <- deg_des_in_dis %>%
  group_by(Gene) %>%
  filter(abs(Discovery_DES_LPS_logFC) == max(abs(Discovery_DES_LPS_logFC))) %>%
  ungroup() 
deg_des_in_dis<-deg_des_in_dis%>%dplyr::select(c(ESid, Gene, Editing_Index, Discovery_DES_LPS_logFC, Discovery_DEG_LPS, Mutation, Location))
dis_deg<-read_tsv("Discovery/Toptable/Gene_Differential/LPS_DEG_Toptable.tsv")
dis_deg<-dis_deg%>%dplyr::select(c(GENENAME, logFC))
dis_deg$Gene<-dis_deg$GENENAME
deg_des_in_dis<-merge(dis_deg, deg_des_in_dis, by="Gene")
deg_des_in_dis$DES_logFC<-deg_des_in_dis$Discovery_DES_LPS_logFC
write_tsv(deg_des_in_dis, "/Users/hyominseo/Desktop/RAJ_RNA_Editing/All_Table/All_DES_DEG/Discovery_LPS_DES_DEG_All.tsv")

deg_des_in<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/All_Table/All_DES_DEG/IFN_DES_DEG_All.tsv")
deg_des_in_dis<-deg_des_in%>%dplyr::filter(grepl("Discovery",Editing_Cohort))
deg_des_in_dis<-deg_des_in_dis%>%dplyr::filter(!grepl("Non_DE",Discovery_DES_IFN))
deg_des_in_dis<-deg_des_in_dis%>%dplyr::filter(!grepl("FALSE",Discovery_DEG_IFNb))
deg_des_in_dis$Editing_Index<-gsub("T:C","A:G",deg_des_in_dis$Editing_Index)
deg_des_in_dis$Editing_Index<-gsub("G:A","C:T",deg_des_in_dis$Editing_Index)
# > length(unique(deg_des_in_dis$ESid))
# [1] 1377
# > length(unique(deg_des_in_dis$Gene))
# [1] 432
deg_des_in_dis <- deg_des_in_dis %>%
  group_by(Gene) %>%
  filter(abs(Discovery_DES_IFN_logFC) == max(abs(Discovery_DES_IFN_logFC))) %>%
  ungroup() 
deg_des_in_dis<-deg_des_in_dis%>%dplyr::select(c(ESid, Gene, Editing_Index, Discovery_DES_IFN_logFC, Discovery_DEG_IFNb, Mutation, Location))
dis_deg<-read_tsv("Discovery/Toptable/Gene_Differential/IFNb_DEG_Toptable.tsv")
dis_deg<-dis_deg%>%dplyr::select(c(GENENAME, logFC))
dis_deg$Gene<-dis_deg$GENENAME
deg_des_in_dis<-merge(dis_deg, deg_des_in_dis, by="Gene")
deg_des_in_dis$DES_logFC<-deg_des_in_dis$Discovery_DES_IFN_logFC

write_tsv(deg_des_in_dis, "/Users/hyominseo/Desktop/RAJ_RNA_Editing/All_Table/All_DES_DEG/Discovery_IFN_DES_DEG_All.tsv")

################################################################################
################################################################################
all_esid<-rbind(M1_anno_all, M2_anno_all, Micro_anno_all)
all_esid<-all_esid%>%filter(!duplicated(across(-ESid)))

### Get all RAW ANNOTAITON AND FILTER BY
dis_anno<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/all_sites_pileup_annotation.tsv")
dis_anno$cohort<-c("Discovery")

rep_anno<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Jacusa/all_sites_pileup_annotation.tsv")
rep_anno$cohort<-c("Replication")

alt_anno<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Alternation/Jacusa//all_sites_pileup_annotation.tsv")
alt_anno$cohort<-c("Alternation")

all_anno<-rbind(dis_anno, rep_anno, alt_anno)
all_anno<-change_esid12_anno(all_anno)

all_esid_gene<-read_tsv("All_Table/All_DES_DEG/BOTH_DES_DEG_All.tsv")
all_anno_filt<-all_anno%>%dplyr::filter(ESid %in% all_esid_gene$ESid)

all_anno_filt <- all_anno_filt %>%
  group_by(ESid) %>%
  mutate(
    DuplicatedCount = n(),
    CohortValues = paste(unique(cohort), collapse = ", ")
  ) %>%
  ungroup()

# Just Cohort counts
cohort_count<-all_anno_filt%>%dplyr::select(c("ESid","DuplicatedCount", "CohortValues"))

all_esid_gene<-merge(all_esid_gene, cohort_count, by='ESid')
all_esid_gene<-all_esid_gene%>%distinct(ESid, .keep_all=TRUE)

################################################################################
################################################################################

library(strex)
exonic_recoding <- all_anno_filt%>%mutate(aa = str_after_first(AAChange.refGene, "p."),
                                          AAsub = gsub(",.*", "", aa),
                                          original_aa = str_extract(AAsub, "^[A-Z]+"),
                                          new_aa = str_sub(AAsub,-1),
                                          original_aa_class = case_when(original_aa %in% c("A","G","I","L","P","V","M") ~ "hydrophilic",
                                                                        original_aa %in% c("F","W","Y") ~ "hydrophobic",
                                                                        original_aa %in% c("R","H","K") ~ "basic",
                                                                        original_aa %in% c("D","E") ~ "acidic",
                                                                        original_aa %in% c("S","T","N","Q","C") ~ "polar uncharged"),
                                          new_aa_class = case_when(new_aa %in% c("A","G","I","L","P","V","M") ~ "hydrophilic",
                                                                   new_aa %in% c("F","W","Y") ~ "hydrophobic",
                                                                   new_aa %in% c("R","H","K") ~ "basic",
                                                                   new_aa %in% c("D","E") ~ "acidic",
                                                                   new_aa %in% c("S","T","N","Q","C") ~ "polar uncharged"),
                                          aa_change = ifelse(!original_aa_class == new_aa_class, "sig", "benign"),
                                          novelty = ifelse(REDIportal_info == ".", "novel", "known")) %>%
  dplyr::select(ESid, Gene.refGene, AAsub, original_aa_class, new_aa_class, aa_change, novelty,cohort)

exonic_recoding<-na.omit(exonic_recoding)

all_esid_gene_exo<- merge(exonic_recoding, all_esid_gene, by = "ESid", all.x = TRUE, all.y = TRUE)

all_esid_gene_exo<-all_esid_gene_exo%>%distinct(ESid, .keep_all=TRUE)
write_tsv(all_esid_gene_exo, "All_Table/all_esid_gene_exonic.tsv")

################################################################################
################################################################################
################################################################################
################################################################################
all_esid_gene_exo_re_table<-all_esid_gene_exo_re%>%
  dplyr::select(c("ESid","AAsub","original_aa_class","new_aa_class","aa_change",
                  "novelty","Gene","CohortValues",'DuplicatedCount'))
all_esid_gene_exo_re_table<- kable(all_esid_gene_exo_re_table, format = "html")  %>%
  kable_classic(full_width = T, html_font = "Cambria")

all_esid_gene_exo_re_table_html <- kable_styling(all_esid_gene_exo_re_table,full_width = FALSE)

all_esid_gene_exo_re_table<-kbl(all_esid_gene_exo_re_table_html, escape=FALSE)

# Save the HTML table to a file
writeLines(all_esid_gene_exo_re_table, "all_esid_gene_exo_re_table.html")
#################################################################



