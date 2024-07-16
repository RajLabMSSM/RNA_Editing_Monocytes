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
library(strex)

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
################################################################################
# ESid Identification Cohort Counts 
################################################################################
M1_anno_all_cohort<-M1_anno_all
M1_anno_all_cohort$Cohort<-"Discovery"
M2_anno_all_cohort<-M2_anno_all
M2_anno_all_cohort$Cohort<-"Replication"
Micro_anno_all_cohort<-Micro_anno_all
Micro_anno_all_cohort$Cohort<-"iMicroglia"

all_esid<-rbind(M1_anno_all_cohort, M2_anno_all_cohort, Micro_anno_all_cohort)
all_esid_cohort_count <- all_esid %>%
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
# Raw annotation binding for Recoding 
################################################################################

M1_raw_anno<-read_tsv("Discovery/Jacusa/all_sites_pileup_annotation.tsv")
M2_raw_anno<-read_tsv("Replication/Jacusa/all_sites_pileup_annotation.tsv")
mic_raw_anno<-read_tsv("Alternation/Jacusa/all_sites_pileup_annotation.tsv")

all_raw_anno<- bind_rows(M1_raw_anno, M2_raw_anno, mic_raw_anno) %>%
  distinct(ESid, .keep_all = TRUE)
all_raw_anno<-change_esid12_anno(all_raw_anno)

################################################################################
# Adding DEG logFC 
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
all_esid_gene$Editing_Index<-gsub("G:A","C:T", all_esid_gene$Editing_Index)
all_esid_gene$Editing_Index<-gsub("T:C","A:G", all_esid_gene$Editing_Index)
all_esid_gene<-all_esid_gene%>%dplyr::select(ESid, contains("DEG"))
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

all_esid$Editing_Index<-gsub("G:A","C:T", all_esid$Editing_Index)
all_esid$Editing_Index<-gsub("T:C","A:G", all_esid$Editing_Index)
all_esid<-all_esid%>%dplyr::filter(grepl("C:T|A:G", Editing_Index))

################################################################################
make_de_des<-function(des_in, colname_in){
  des_de<-des_in%>%
    filter(grepl("A:G_DE|C:T_DE", DE_Index)) %>%
    select(ESid, logFC)
  des_de
}
################################################################################

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

################################################################################
all_esid_merge<-all_esid_merge%>%select(c(
  ESid, Gene, Editing_Index, Discovery_DES_LPS, Discovery_DES_LPS_logFC,
  Discovery_DES_IFN, Discovery_DES_IFN_logFC, Replication_DES_LPS,
  Replication_DES_LPS_logFC, Replication_DES_IFN, Replication_DES_IFN_logFC,
  Alternation_DES_LPS, Alternation_DES_LPS_logFC, Alternation_DES_IFN,
  Alternation_DES_IFN_logFC, Location, Mutation, known_a_i, rep_type
))

all_esid_merge<-merge(all_esid_merge,all_esid_cohort_count, by ="ESid")
all_esid_merge<-merge(all_esid_merge,all_esid_gene, by ="ESid")
################################################################################
write_tsv(all_esid_merge, "All_Table/All_esid.tsv")
################################################################################

all_esid_exonic<-all_esid_merge%>%dplyr::filter(Location == "exonic")
all_raw_anno_exonic<-all_raw_anno%>%dplyr::filter(ESid %in% all_esid_exonic$ESid)
all_raw_anno_exonic<-all_raw_anno_exonic%>%dplyr::select(c(ESid, AAChange.refGene,
                                                           REDIportal_info,Gene.refGene))

all_esid_exonic<-merge(all_esid_exonic, all_raw_anno_exonic, by ="ESid")

exonic_recoding <- all_esid_exonic%>%mutate(aa = str_after_first(AAChange.refGene, "p."),
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
  dplyr::select(ESid, Gene.refGene, AAsub, original_aa_class, new_aa_class, aa_change, novelty,
                everything())#,cohort)

de_esid<- apply(exonic_recoding, 1, function(row) {
  any(grepl("A:G_DE|C:T_DE", row))
})

# 173 ESid that is DE at least in one category,
# associated with 107 genes
recoding_de <- exonic_recoding[de_esid, ]
# how many times was it DE 
# DE 2 times, 133
# DE 4 times, 29
# DE 6 times, 11
recoding_de$DE_Cohort_Number <- apply(recoding_de, 1, function(row) {
  sum(grepl("_DE", row))
})
write_tsv(recoding_de, "All_Table/DE_All_Recoding.tsv")

recoding_de_sig<-recoding_de%>%dplyr::filter(grepl("sig",aa_change))
#write_tsv(recoding_de_sig, "All_Table/DE_Sig_Recoding.tsv")
################################################################################
# DE SIG Recoding all cohort table
recoding_de_sig<-recoding_de_sig%>%dplyr::select(c(ESid, Gene, Editing_Index, 
                                                        aa_change, novelty, Location,
                                                        Mutation, Editing_Count))
write_tsv(recoding_de_sig, "All_Table/DE_Sig_Recoding.tsv")
################################################################################
################################################################################
# Sitmulation-wise DEG DES
all_esid_merge<-read_tsv("All_Table/All_esid.tsv")
write_tsv(all_esid_merge, "All_Table/All_DES_DEG/BOTH_DES_DEG_All.tsv")
################################################################################
lps_all_esid_gene<-all_esid_merge%>%dplyr::select(
  ESid, Gene, Editing_Index, Editing_Count, Editing_Cohort,
  Discovery_DES_LPS, Discovery_DES_LPS_logFC, 
  Replication_DES_LPS, Replication_DES_LPS_logFC, 
  Alternation_DES_LPS, Alternation_DES_LPS_logFC,
  Discovery_DEG_LPS, Replication_DEG_LPS, Alternation_DEG_LPS,
  everything())%>%
  select(-matches("IFN|IFNb|IFNg"))

write_tsv(lps_all_esid_gene, "All_Table/All_DES_DEG/LPS_DES_DEG_All.tsv")

ifn_all_esid_gene<-all_esid_merge%>%dplyr::select(
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
# pathway GP 
################################################################################

multi_gp = gost(
  list("all-regulated" = re_gene$GENENAME),
  multi_query = FALSE,
  evcodes = TRUE, 
  user_threshold = 0.05, # P val 
  sources = c("GO:BP", "KEGG") 
)


gp_mod <- multi_gp$result[,c("query", "source", "term_id","term_name",
                             "p_value","query_size","intersection_size",
                             "term_size","effective_domain_size","intersection")]
# Only getting GO:BP and KEGG source
gp_mod <- gp_mod%>%dplyr::filter(grepl("GO:BP|KEGG",gp_mod$source))
# Making Gene Ratio and BGRatio
gp_mod$GeneRatio <- paste0(gp_mod$intersection_size, "/",gp_mod$query_size)
gp_mod$BgRatio <- paste0(gp_mod$term_size, "/",gp_mod$effective_domain_size)

names(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust",
                   "query_size", "Count", "term_size","effective_domain_size",
                   "geneID", "GeneRatio","BgRatio")

gp_mod$geneID <- gsub(",", "/", gp_mod$geneID) 

message("Writing GP Enrich result Matrix")
#write_tsv(gp_mod, "Result/gp_enrich_result_trial.tsv")

gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod)
gp_mod_enrich  <- new("enrichResult", result = gp_mod)

gp_mod_enrich<-as.data.frame(gp_mod_enrich@result)
write_tsv(gp_mod_enrich, "All_Table/exonic_sig_recoding_gene_gp_enrich.tsv")
################################################################################ 
########################## GP ENRICH PLOT SETTING ############################## 
################################################################################
text_size=9
legend_text_size=8
legend_size=5
number_category=15
gp_barplot<-barplot(gp_mod_enrich, showCategory = number_category, font.size = text_size) + 
  ggplot2::ylab("Intersection size")+
  labs(x="Count",y ="Intersection Size", title="Sig Recoding Genes")+theme_bw()+
  theme(axis.text.y = element_text(color = "black", size = text_size),
        axis.title = element_text(size = 8),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=legend_text_size),
        legend.text = element_text(size=legend_text_size),
        plot.title = element_text(size=text_size,face = "bold",hjust=0.5))+ 
  ggeasy::easy_rotate_y_labels(angle = 30, side = c("right"))+ 
  guides(color = guide_legend(override.aes = list(size=legend_size)))
gp_barplot
message("Making Plots")
ggsave(plot=gp_barplot,filename="/Users/hyominseo/Desktop/RAJ_RNA_Editing/Figure/Temp/Sig_Recoding_Gene_90_Pathway.jpg",
       width = 8, height = 10,dpi = 600)

