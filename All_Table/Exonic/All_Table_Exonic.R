################################# DISCOVERY ####################################
# start from RAW ANNOTAION this needs to be re-written

################################################################################
############### METADATA
M1_meta<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Metadata/165_Monocytes_Metadata.tsv")

############### TPM
M1_tpm<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Metadata/Gene_Expression/Processed_Normalized_TPM.tsv")%>%column_to_rownames("GENEID")
names(M1_tpm) <- gsub(x = names(M1_tpm), pattern = "\\.", replacement = "-")

############### DEG
M1_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Gene_Differential/LPS_DEG_Toptable.tsv")
M1_DEG_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Gene_Differential/IFNb_DEG_Toptable.tsv")

############### ANNOTATION
M1_anno_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Local_Pipeline_Result/165_Navarro/Processed_Annotation_All.tsv")
M1_anno_all<-change_esid12_anno(M1_anno_all)

M1_anno_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Local_Pipeline_Result/165_Navarro/Processed_Annotation_Type.tsv")
M1_anno_type<-change_esid12_anno(M1_anno_type)

############### EDITING
M1_editing_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Local_Pipeline_Result/165_Navarro/Processed_Editing_All.tsv")
M1_editing_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Local_Pipeline_Result/165_Navarro/Processed_Editing_Type.tsv")

M1_AG<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Local_Pipeline_Result/165_Navarro/Processed_AG_Editing.tsv")
M1_CT<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/Local_Pipeline_Result/165_Navarro/Processed_CT_Editing.tsv")

############### DES
M1_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Editing_Differential/LPS_1_Toptable.tsv")
lps_DE<-M1_DES_LPS%>%dplyr::filter(!grepl("Non_DE", DE_Index))
lps_DE_anno<-M1_anno_type%>%dplyr::filter(ESid %in% lps_DE$ESid)
lps_DE_anno$Location<-sub("_", "\n", lps_DE_anno$Location) 
lps_DE_anno$Mutation<-sub(" ", "\n", lps_DE_anno$Mutation) 

M1_DES_Model2_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Editing_Differential/LPS_2_Toptable.tsv")

M1_DES_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Editing_Differential/IFNb_1_Toptable.tsv")
ifnb_DE<-M1_DES_IFNb%>%dplyr::filter(!grepl("Non_DE", DE_Index))
ifnb_DE_anno<-M1_anno_type%>%dplyr::filter(ESid %in% ifnb_DE$ESid)
ifnb_DE_anno$Location<-sub("_", "\n", ifnb_DE_anno$Location) 
ifnb_DE_anno$Mutation<-sub(" ", "\n", ifnb_DE_anno$Mutation) 


M1_DES_Model2_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Editing_Differential/IFNb_2_Toptable.tsv")

############### DES-DEG
M1_anno_deg_des<-M1_anno_all%>%dplyr::select(c(ESid,Location, Mutation))

make_deg_des(M1_DEG_LPS,M1_DES_LPS,M1_anno_deg_des, "/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/LPS_DEG_DES.tsv")
make_deg_des(M1_DEG_IFNb,M1_DES_IFNb,M1_anno_deg_des, "/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/IFNb_DEG_DES.tsv")

M1_lps_deg_des<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/LPS_DEG_DES.tsv")
M1_ifnb_deg_des<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/IFNb_DEG_DES.tsv")



################################################################################ 
################################ REPLICATION ###################################
################################################################################
############### METADATA
M2_meta<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Metadata/105_Monocytes_Metadata.tsv")

############### TPM
M2_tpm<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Metadata/Gene_Expression/Processed_Normalized_TPM.tsv")%>%column_to_rownames("GENEID")
names(M2_tpm) <- gsub(x = names(M2_tpm), pattern = "\\.", replacement = "-")

############### DEG
M2_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/Gene_Differential/LPS_DEG_Toptable.tsv")
M2_DEG_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/Gene_Differential/IFNg_DEG_Toptable.tsv")

############### ANNO
M2_anno_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Local_Pipeline_Result/105_Navarro/Processed_Annotation_All.tsv")
M2_anno_all<-change_esid12_anno(M2_anno_all)

M2_anno_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Local_Pipeline_Result/105_Navarro/Processed_Annotation_Type.tsv")
M2_anno_type<-change_esid12_anno(M2_anno_type)

############### EDITING
M2_editing_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Local_Pipeline_Result/105_Navarro/Processed_Editing_All.tsv")
M2_editing_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Local_Pipeline_Result/105_Navarro/Processed_Editing_Type.tsv")

M2_AG<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Local_Pipeline_Result/105_Navarro/Processed_AG_Editing.tsv")
M2_CT<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/Local_Pipeline_Result/105_Navarro/Processed_CT_Editing.tsv")

############### DES
M2_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/Editing_Differential/LPS_1_Toptable.tsv")
M2_DES_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/Editing_Differential/IFNg_1_Toptable.tsv")

############### DES-DEG
M2_anno_deg_des<-M2_anno_all%>%dplyr::select(c(ESid,Location, Mutation))

make_deg_des(M2_DEG_LPS,M2_DES_LPS,M2_anno_deg_des, "/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/LPS_DEG_DES.tsv")
make_deg_des(M2_DEG_IFNg,M2_DES_IFNg,M2_anno_deg_des, "/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/IFNg_DEG_DES.tsv")

M2_lps_deg_des<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/LPS_DEG_DES.tsv")
M2_ifng_deg_des<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/105_Stim_Monocytes/TopTable/IFNg_DEG_DES.tsv")




################################################################################ 
################################ ALTERNATION ###################################
################################################################################
############### METADATA
micro_meta<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/Metadata/efthymiou_meta.tsv")

############### TPM
############### DEG
Micro_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/TopTable/Gene_Differential/LPS_DEG_Toptable.tsv")
Micro_DEG_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/TopTable/Gene_Differential/IFNg_DEG_Toptable.tsv")

############### ANNOTATION
Micro_anno_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/Local_Pipeline_Result/Efthymiou/Processed_Annotation_All.tsv")
Micro_anno_all<-change_esid12_anno(Micro_anno_all)

Micro_anno_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/Local_Pipeline_Result/Efthymiou/Processed_Annotation_Type.tsv")
Micro_anno_type<-change_esid12_anno(Micro_anno_type)

############### EDITING
Micro_editing_all<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/Local_Pipeline_Result/Efthymiou/Processed_Editing_All.tsv")

Micro_AG<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/Local_Pipeline_Result/Efthymiou/Processed_AG_Editing.tsv")
Micro_CT<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/Local_Pipeline_Result/Efthymiou/Processed_CT_Editing.tsv")

############### DES
Micro_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/TopTable/Editing_Differential/LPS_1_Toptable.tsv")
Micro_DES_IFNg<-read_tsv("/Users/hyominseo/Desktop/RAJ_Microglia/Efthymiou/TopTable/Editing_Differential/IFNg_1_Toptable.tsv")

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


exonic_recoding<-exonic_recoding%>%dplyr::select(-c(cohort))

all_esid_gene_exo<- merge(exonic_recoding, all_esid_gene, by = "ESid", all.x = TRUE, all.y = TRUE)

all_esid_gene_exo<-all_esid_gene_exo%>%distinct(ESid, .keep_all=TRUE)
write_tsv(all_esid_gene_exo, "All_Table/all_esid_gene_exonic.tsv")

# duplicated_rows <- all_esid_gene_exo %>%
#   filter(duplicated(ESid) | duplicated(ESid, fromLast = TRUE))

recoding<-read_tsv("All_Table/all_esid_gene_exonic.tsv")
recoding<-recoding%>%dplyr::filter(grepl("sig",aa_change))
re_gene<-recoding%>%distinct(Gene.refGene.keep_all=TRUE)
write_tsv(re_gene, "All_Table/exonic_sig_recoding_gene.tsv")
################################################################################
################################################################################
#################################################################################
#################################################################################

bar_plot <- ggplot(all_esid_gene_exo, aes(x = CohortValues, y=DuplicatedCount)) +
  geom_col(aes(color=DuplicatedCount, fill=DuplicatedCount))+#, position = position_dodge(0.8), width = 0.7)+
  #scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
  #scale_fill_manual(values = c("#0073C2FF", "#EFC000FF"))+
  geom_text(aes(label = DuplicatedCount), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3, color = "black") +
  labs(title = "Sites Found in at least # Cohorts",
       x = "Cohort Occurrences",
       y = "Site Count") +
  theme_bw()+theme(legend.title = element_text(size=text_size),
                   legend.text = element_text(size=text_size),
                   axis.text.x = element_text(size = axis_size,angle = 20, hjust = 1),
                   axis.text.y = element_text(size = axis_size),
                   axis.title = element_text(size = axis_size),
                   plot.title = element_text(size=title_size, hjust=0.5),
                   legend.margin=margin(0,0,0,0),
                   legend.box.margin=margin(-5,-5,-5,-5))+
  theme(legend.position = none)


combined_color<-c("sig"="#ff8565", "Unidentified"= "#d6d6d6","benign" = "#40B5AD")
all_esid_gene_exo_re<-all_esid_gene_exo%>%dplyr::filter(grepl("sig|benign",aa_change))

anno_plot<- ggplot(all_esid_gene_exo_re, aes(x = "", fill = aa_change)) +
  geom_bar(stat = "count", width = 0.8) +
  scale_fill_manual(values=combined_color)+
  geom_text(
    aes(label = after_stat(count)),stat = "count",
    position = position_stack(vjust = 0.5),
    size = 3,color = "black") +
  labs(title = "445 Exonic Sites Category",x = "",y = "Count") +
  theme_bw()+theme(legend.title = element_text(size=text_size),
                   legend.text = element_text(size=text_size),
                   axis.text.x = element_text(size = axis_size),
                   axis.text.y = element_text(size = axis_size),
                   axis.title = element_text(size = axis_size),
                   plot.title = element_text(size=title_size, hjust=0.5),
                   legend.margin=margin(0,0,0,0),
                   legend.box.margin=margin(-5,-5,-5,-5))+
  theme(legend.position = "right")

bar_plots<-ggarrange(bar_plot, anno_plot, 
                     labels=c("A","B"),ncol=2,
                     font.label = list(color = "black", size = text_size,
                                       margin=c(label_margin,label_margin)))

ggsave(plot=bar_plots,filename="/Users/hyominseo/Desktop/Editing_Figures/Temp_Figure/All_sites_bars.jpg",
       width = 10, height = 5,dpi = 600)

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

