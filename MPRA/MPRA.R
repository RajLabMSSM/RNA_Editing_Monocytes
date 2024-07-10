# Discovery 
# I can probably just delete this script
#```{r stim_1, echo= TRUE, message = FALSE, warning=FALSE, eval=TRUE}
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")
out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_raj_discovery<-out_raj%>%dplyr::filter(grepl("165", cohorts))

dis_lps<-out_raj_discovery%>%dplyr::filter(grepl("LPS",cohorts))
dis_lps$significant<-gsub("Y","LPS_Y",dis_lps$significant)
dis_ifn<-out_raj_discovery%>%dplyr::filter(grepl("IFN",cohorts))
dis_ifn$significant<-gsub("Y","IFNb_Y",dis_ifn$significant)

dis_lps_com<-dis_lps%>%dplyr::select(c(ESid, logFC, logP, significant))
dis_lps_com_de<-dis_lps_com%>%dplyr::filter(grepl("LPS_Y",significant))
dis_ifn_com<-dis_ifn%>%dplyr::select(c(ESid, logFC, logP, significant))
dis_ifn_com_de<-dis_ifn_com%>%dplyr::filter(grepl("IFNb_Y",significant))

both<-as.data.frame(intersect(dis_lps_com_de$ESid, dis_ifn_com_de$ESid))
colnames(both)<-c("ESid")
#colnames(dis_ifn_com)[4]<-"significant_ifn"

dis_com<-rbind(dis_lps_com, dis_ifn_com)
dis_com<-dis_com%>%mutate(significant = case_when(
  ESid %in% both$ESid ~ "LPS_IFNb",
  TRUE ~ significant
))

mpra_dis_com<-vol_plot_mpra(dis_com, "MPRA: Discovery LPS/ IFNb Sites")
################################################################################ 
M1_anno_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_Type.tsv")
M1_anno_type<-change_esid12_anno(M1_anno_type)
M1_anno_type<-rename_location(M1_anno_type)
M1_anno_type$GENENAME<-M1_anno_type$Gene
M1_anno_type<-M1_anno_type%>%dplyr::select(c(GENENAME,Editing_Index, 
                                             Location, Mutation, known_a_i, 
                                             rep_type, ESid))

mpra_in<-read_tsv("MPRA/Input/All_Monocytes_DESites_3UTR_Model1.tsv")
# 3067
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")
col_names<-colnames(mpra_out)
col_names_new <- ifelse(col_names %in% c("ESid","Count"), 
                        col_names, paste0(col_names, "_MPRA"))
colnames(mpra_out) <- col_names_new
################################################################################ 

M1_lps<-merge(M1_DEG_LPS, M1_anno_type, by ="GENENAME", .keep_all=TRUE)
M1_lps<-merge(M1_lps,mpra_out, by = "ESid", .keep_all=TRUE)
M1_lps <- M1_lps[!is.na(M1_lps$logFC_MPRA), ]
M1_lps$significant_MPRA<-gsub("Y","LPS_Y",M1_lps$significant_MPRA)

M1_ifnb<-merge(M1_DEG_IFNb, M1_anno_type, by ="GENENAME", .keep_all=TRUE)
M1_ifnb<-merge(M1_ifnb,mpra_out, by = "ESid", .keep_all=TRUE)
M1_ifnb <- M1_ifnb[!is.na(M1_ifnb$logFC_MPRA), ]
M1_ifnb$significant_MPRA<-gsub("Y","IFNb_Y",M1_ifnb$significant_MPRA)

################################################################################ 
deg_lps<-deg_comp_plot_mpra(M1_lps, "Discovery_LPS: DEG")
deg_ifnb<-deg_comp_plot_mpra(M1_ifnb, "Discovery_IFNb: DEG")



mpra_deg<-ggarrange(deg_lps,deg_ifnb,nrow=2)
second_row<-ggarrange(mpra_deg, mpra_dis_com,mpra_deg,labels=c("D","E","F"),ncol=3,
                      font.label = list(color = "black", size = text_size,
                                        margin=c(label_margin,label_margin)))

second_row
################################################################################ 
################################################################################ 
mpra_color <- c("Basal" = "#c5c5c5", "LPS_Y" = "#FFA800", 
                "Discover_LPS" ="#FFA800", 
                "IFNb_Y" = "#3A867F","N"="#D3D3D3",
                "Discover_IFNb" = "#3A867F",
                "LPS_MPRA"="#FFA800","IFNb_MPRA" ="#3A867F",
                "LPS_IFNb" = "#F78520",
                "Y" ="#F78520"  , "N" ="#D3D3D3" )


vol_plot_mpra<-function(table_in, plot_title){
  vol_plot<- ggplot(table_in, aes(x=logFC, y=logP,col = miRNA,shape=significant))+ 
    geom_point(size =point_size+1, alpha=0.8)+
    scale_shape_manual(values = c("N" = 16, "LPS_Y" = 17, "IFNb_Y" = 16, "LPS_IFNb" = 18))+
    scale_color_manual(values =mpra_color)+
    labs(x = "logFC", y = "-log10(P.Value)",title=plot_title)+
    geom_vline(xintercept=0,  col ="#818589", linetype="longdash")+
    theme_classic()+
    theme(legend.title = element_text(size=text_size),
          legend.text = element_text(size=text_size),
          axis.text.x = element_text(size = axis_size,color="black"),
          axis.text.y = element_text(size = axis_size,color="black"),
          axis.title = element_text(size = axis_size),
          plot.title = element_text(size=title_size, hjust=0.5),
          plot.subtitle = element_text(size=title_size, hjust=0.5),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-3,-3,-3,-3))+
    theme(legend.position = "top")+
    guides(color = guide_legend(override.aes = list(size=legend_size)))
  vol_plot
}

################################################################################ 

new_mpra_table<-read_tsv("MPRA/Result/monocyte_miRNA_seed_sites_info.txt")
new_mpra<-read_tsv("MPRA/Result/monocyte_MPRA_results_Winston.txt")

mpra_out<-new_mpra
out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_raj_discovery<-out_raj%>%dplyr::filter(grepl("165", cohorts))

dis_lps<-out_raj_discovery%>%dplyr::filter(grepl("LPS",cohorts))
dis_lps$significant<-gsub("Y","LPS_Y",dis_lps$significant)
dis_ifn<-out_raj_discovery%>%dplyr::filter(grepl("IFN",cohorts))
dis_ifn$significant<-gsub("Y","IFNb_Y",dis_ifn$significant)

dis_lps_com<-dis_lps%>%dplyr::select(c(ESid, logFC, logP, significant,miRNA))
dis_lps_com_de<-dis_lps_com%>%dplyr::filter(grepl("LPS_Y",significant))
dis_ifn_com<-dis_ifn%>%dplyr::select(c(ESid, logFC, logP, significant,miRNA))
dis_ifn_com_de<-dis_ifn_com%>%dplyr::filter(grepl("IFNb_Y",significant))

both<-as.data.frame(intersect(dis_lps_com_de$ESid, dis_ifn_com_de$ESid))
colnames(both)<-c("ESid")
#colnames(dis_ifn_com)[4]<-"significant_ifn"

dis_com<-rbind(dis_lps_com, dis_ifn_com)
dis_com<-dis_com%>%mutate(significant = case_when(
  ESid %in% both$ESid ~ "LPS_IFNb",
  TRUE ~ significant
))

mpra_dis_com<-vol_plot_mpra(dis_com, "MPRA: Discovery LPS/ IFNb Sites")
ggsave(plot=mpra_dis_com,filename="MPRA/Result/MPRA_miRNA_Sig_vol.jpg",width = 8, height = 8,dpi = 600)

sig_mirna<-dis_com%>%dplyr::filter(grepl("Y",miRNA))


text_size=10
legend_text_size=10
legend_size=8
number_category=15
library(clusterProfiler)

generate_path<-function( df_in,file_name, title_in, plot_name){
  #df_in<-mg_pro_gene_5
  
  df<-df_in
  multi_gp = gost(
    list("all-regulated" = df$GENENAME),
    multi_query = FALSE,
    evcodes = TRUE, 
    user_threshold = 0.05, # P val 
    sources = c("GO:BP", "KEGG") 
  )
  gp_mod <- multi_gp$result[,c("query", "source", "term_id","term_name",
                               "p_value","query_size","intersection_size",
                               "term_size","effective_domain_size","intersection")]
  
  gp_mod <- gp_mod%>%dplyr::filter(grepl("GO:BP|KEGG",gp_mod$source))
  
  gp_mod$GeneRatio <- paste0(gp_mod$intersection_size, "/",gp_mod$query_size)
  gp_mod$BgRatio <- paste0(gp_mod$term_size, "/",gp_mod$effective_domain_size)
  
  names(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust",
                     "query_size", "Count", "term_size","effective_domain_size",
                     "geneID", "GeneRatio","BgRatio")
  
  gp_mod$geneID <- gsub(",", "/", gp_mod$geneID) 
  
  message("Writing GP Enrich result Matrix")
  
  gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod)
  gp_mod_enrich  <- new("enrichResult", result = gp_mod)
  
  gp_mod_enrich_df<-as.data.frame(gp_mod_enrich@result)
  write_tsv(gp_mod_enrich_df, file_name)
  
  gp_barplot<-barplot(gp_mod_enrich, showCategory = number_category, 
                      font.size = text_size) + 
    ggplot2::ylab("Intersection size")+
    labs(x="Count",y ="Intersection Size", title=title_in)+theme_bw()+
    theme(axis.text.y = element_text(color = "black", size = text_size),
          axis.title = element_text(size = text_size),
          legend.key.height= unit(0.4, 'cm'),
          legend.key.width= unit(0.2, 'cm'),
          legend.title = element_text(size=legend_text_size),
          legend.text = element_text(size=legend_text_size),
          plot.title = element_text(size=text_size,face = "bold",hjust=0.5))+ 
    ggeasy::easy_rotate_y_labels(angle = 30, side = c("right"))+ 
    guides(color = guide_legend(override.aes = list(size=legend_size)))
  #gp_barplot
  message("Making Plots")
  ggsave(plot=gp_barplot,filename=plot_name,width = 8, height = 10,dpi = 600)
}

#df_in,file_name, title_in, plot_name
generate_path(gene_sig_mirna, "MPRA/Pathway/Discovery_miRNA_Sig_Pathway_enrich.tsv",
              "89 Genes for 187 Sig miRNA MPRA Sites", "MPRA/Pathway/Discovery_miRNA_Sig_Pathway_enrich.jpg")

sig_mirna_path<-read_tsv("MPRA/Pathway/Discovery_miRNA_Sig_Pathway_enrich.tsv")

################################################################################ 
M1_anno_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_Type.tsv")
M1_anno_type<-change_esid12_anno(M1_anno_type)
M1_anno_type<-rename_location(M1_anno_type)
M1_anno_type$GENENAME<-M1_anno_type$Gene
M1_anno_type<-M1_anno_type%>%dplyr::select(c(GENENAME,Editing_Index, 
                                             Location, Mutation, known_a_i, 
                                             rep_type, ESid))
################################################################################ 
anno_sig_mirna<-M1_anno_type%>%dplyr::filter(ESid %in% sig_mirna$ESid)
gene_sig_mirna<-anno_sig_mirna%>%distinct(GENENAME)

mpra_in<-read_tsv("MPRA/Input/All_Monocytes_DESites_3UTR_Model1.tsv")
# 3067
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")
col_names<-colnames(mpra_out)
col_names_new <- ifelse(col_names %in% c("ESid","Count"), 
                        col_names, paste0(col_names, "_MPRA"))
colnames(mpra_out) <- col_names_new
################################################################################ 

M1_lps<-merge(M1_DEG_LPS, M1_anno_type, by ="GENENAME", .keep_all=TRUE)
M1_lps<-merge(M1_lps,mpra_out, by = "ESid", .keep_all=TRUE)
M1_lps <- M1_lps[!is.na(M1_lps$logFC_MPRA), ]
M1_lps$significant_MPRA<-gsub("Y","LPS_Y",M1_lps$significant_MPRA)

M1_ifnb<-merge(M1_DEG_IFNb, M1_anno_type, by ="GENENAME", .keep_all=TRUE)
M1_ifnb<-merge(M1_ifnb,mpra_out, by = "ESid", .keep_all=TRUE)
M1_ifnb <- M1_ifnb[!is.na(M1_ifnb$logFC_MPRA), ]
M1_ifnb$significant_MPRA<-gsub("Y","IFNb_Y",M1_ifnb$significant_MPRA)

################################################################################ 
deg_lps<-deg_comp_plot_mpra(M1_lps, "Discovery_LPS: DEG")
deg_ifnb<-deg_comp_plot_mpra(M1_ifnb, "Discovery_IFNb: DEG")



mpra_deg<-ggarrange(deg_lps,deg_ifnb,nrow=2)
second_row<-ggarrange(mpra_deg, mpra_dis_com,mpra_deg,labels=c("D","E","F"),ncol=3,
                      font.label = list(color = "black", size = text_size,
                                        margin=c(label_margin,label_margin)))

second_row






M1_DEG_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Gene_Differential/LPS_DEG_Toptable.tsv")
M1_DEG_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Gene_Differential/IFNb_DEG_Toptable.tsv")

col_names<-colnames(M1_DEG_LPS)
col_names_new <- ifelse(col_names %in% c("ESid", "GENEID", "GENENAME"), 
                        col_names, paste0(col_names, "_DEG"))
colnames(M1_DEG_LPS) <- col_names_new

col_names<-colnames(M1_DEG_IFNb)
col_names_new <- ifelse(col_names %in% c("ESid", "GENEID", "GENENAME"), 
                        col_names, paste0(col_names, "_DEG"))
colnames(M1_DEG_IFNb) <- col_names_new

# M1_DES_LPS<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Editing_Differential/LPS_1_Toptable.tsv")
# M1_DES_IFNb<-read_tsv("/Users/hyominseo/Desktop/RAJ_Monocytes/165_Stim_Monocytes/TopTable/Editing_Differential/IFNb_1_Toptable.tsv")
# 
# col_names<-colnames(M1_DES_LPS)
# col_names_new <- ifelse(col_names %in% c("ESid", "Gene"), 
#                         col_names, paste0(col_names, "_DES"))
# colnames(M1_DES_LPS) <- col_names_new
# 
# col_names<-colnames(M1_DES_IFNb)
# col_names_new <- ifelse(col_names %in% c("ESid", "Gene"), 
#                         col_names, paste0(col_names, "_DES"))
# colnames(M1_DES_IFNb) <- col_names_new

############### ANNOTATION
M1_anno_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_Type.tsv")
M1_anno_type<-change_esid12_anno(M1_anno_type)
M1_anno_type<-rename_location(M1_anno_type)
M1_anno_type$GENENAME<-M1_anno_type$Gene
M1_anno_type<-M1_anno_type%>%dplyr::select(c(GENENAME,Editing_Index, 
                                             Location, Mutation, known_a_i, 
                                             rep_type, ESid))


# M1_editing_type<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Jacusa/Processed_Editing_Type.tsv")
# 
# mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")
# col_names<-colnames(mpra_out)
# col_names_new <- ifelse(col_names %in% c("ESid","Count"), 
#                         col_names, paste0(col_names, "_MPRA"))
# colnames(mpra_out) <- col_names_new


mpra_in<-read_tsv("MPRA/Input/All_Monocytes_DESites_3UTR_Model1.tsv")
# 3067
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")
col_names<-colnames(mpra_out)
col_names_new <- ifelse(col_names %in% c("ESid","Count"), 
                        col_names, paste0(col_names, "_MPRA"))
colnames(mpra_out) <- col_names_new

out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_raj_dis<-out_raj%>%dplyr::filter(grepl("*165",cohorts))
out_raj_dis<-na.omit(out_raj_dis)


M1_lps<-merge(M1_DEG_LPS, M1_anno_type, by ="GENENAME", .keep_all=TRUE)
#M1_lps<-merge(M1_lps, M1_DES_LPS, by = "ESid", .keep_all=TRUE)
M1_lps<-merge(M1_lps,mpra_out, by = "ESid", .keep_all=TRUE)
M1_lps <- M1_lps[!is.na(M1_lps$logFC_MPRA), ]

M1_ifnb<-merge(M1_DEG_IFNb, M1_anno_type, by ="GENENAME", .keep_all=TRUE)
#M1_ifnb<-merge(M1_ifnb, M1_DES_IFNb, by = "ESid", .keep_all=TRUE)
M1_ifnb<-merge(M1_ifnb,mpra_out, by = "ESid", .keep_all=TRUE)
M1_ifnb <- M1_ifnb[!is.na(M1_ifnb$logFC_MPRA), ]


deg_comp_plot <- function(comp_in, plot_title) {
  log_comp_plot <- ggplot(comp_in, aes(x = logFC_DEG, y = logFC_MPRA, 
                                       col = significant_MPRA)) +
    geom_point(size = point_size, alpha = 0.7) +
    #scale_color_manual(values = combined_color)+
    #stat_cor(label.x.npc="left",label.y.npc ="top") +
    geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +  
    geom_vline(xintercept = c(0), col = "black", linetype="longdash") +
    geom_hline(yintercept = c(0), col = "black", linetype="longdash") +
    labs(x = "DEG logFC:Discovery", y = "MPRA logFC:Discovery",title=plot_title) +
    theme_classic() +
    theme(legend.title = element_text(size=text_size),
          legend.text = element_text(size=text_size),
          axis.text.x = element_text(size = axis_size-0.5,color="black"),
          axis.text.y = element_text(size = axis_size-0.5,color="black"),
          axis.title = element_text(size = axis_size),
          plot.title = element_text(size=title_size, hjust=0.5),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-3,-3,-3,-3))+
    theme(legend.position = "top")+
    guides(color = guide_legend(override.aes = list(size=legend_size)))
  log_comp_plot
}


deg_lps<-deg_comp_plot(M1_lps, "Discovery_LPS: DEG")
deg_ifnb<-deg_comp_plot(M1_ifnb, "Discovery_IFNb: DEG")

des_comp_plot <- function(comp_in, plot_title) {
  log_comp_plot <- ggplot(comp_in, aes(x = logFC_DES, y = logFC_MPRA, 
                                       col = significant_MPRA)) +
    geom_point(size = point_size, alpha = 0.7) +
    #scale_color_manual(values = combined_color)+
    #stat_cor(label.x.npc="left",label.y.npc ="top") +
    geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +  
    geom_vline(xintercept = c(0), col = "black", linetype="longdash") +
    geom_hline(yintercept = c(0), col = "black", linetype="longdash") +
    labs(x = "DES logFC:Discovery", y = "MPRA logFC:Discovery",title=plot_title) +
    theme_classic() +
    theme(legend.title = element_text(size=text_size),
          legend.text = element_text(size=text_size),
          axis.text.x = element_text(size = axis_size-0.5,color="black"),
          axis.text.y = element_text(size = axis_size-0.5,color="black"),
          axis.title = element_text(size = axis_size),
          plot.title = element_text(size=title_size, hjust=0.5),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-3,-3,-3,-3))+
    theme(legend.position = "top")+
    guides(color = guide_legend(override.aes = list(size=legend_size)))
  log_comp_plot
}


des_lps<-des_comp_plot(M1_lps, "Discovery_LPS: DES")
des_ifnb<-des_comp_plot(M1_ifnb, "Discovery_IFNb: DES")

M1_deg_mpra_plots<-ggarrange(deg_lps, deg_ifnb, ncol=2)

ggsave(plot=M1_mpra_plots,filename="Figure/MPRA/MPRA_figures.jpg",width = 12, height = 12,
       dpi = 600)

# 2150
mpra_in<-read_tsv("MPRA/Input/All_Monocytes_DESites_3UTR_Model1.tsv")
# 3067
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")

out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_raj_dis<-out_raj%>%dplyr::filter(grepl("*165",cohorts))
out_raj_dis<-na.omit(out_raj_dis)


table(mpra_out$category)
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")
out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_n_cont<-mpra_out%>%dplyr::filter(grepl("negative_control",category))
out_p_cont<-mpra_out%>%dplyr::filter(grepl("positive_control",category))

# 2006 RAJ sites, that has LogFC 
# FDR >5 (non sig) 1473
# Positive Sig sites : 369
# Negative Sig sites : 288
mpra_vol_plot<-vol_plot_mpra(out_raj, "MPRA: All Dis/Rep/IFNb and LPS, 2007 ESid","UP: 369, Down: 288, Nonsig: 1473")
ggsave(plot=mpra_vol_plot,filename="Figure/MPRA/MPRA_volplot.jpg",width = 8, 
       height = 8,dpi = 600)

# ifnb_DE_anno  just do it again 
M1_anno_all$GENENAME<-M1_anno_all$Gene
M1_anno_all_select<-M1_anno_all%>%dplyr::select(c(GENENAME,Editing_Index, 
                                                  Location, Mutation, known_a_i, 
                                                  rep_type, ESid))
M1_DEG_IFNb_anno<-merge(M1_DEG_IFNb, M1_anno_all_select, by ='GENENAME')
M1_DEG_IFNb_anno<-M1_DEG_IFNb_anno%>%column_to_rownames("ESid")

colnames(M1_DEG_IFNb_anno) <- paste0(colnames(M1_DEG_IFNb_anno), "_DEG")

colnames(M1_DEG_IFNb_anno) <- paste0(colnames(M1_DEG_IFNb_anno), "_DEG")
# Just getting Non-Control sites in mpra_out
M1_DEG_IFNb_anno_mpra<-M1_DEG_IFNb_anno%>%dplyr::filter(ESid %in% mpra_out$ESid)
M1_DEG_DES_IFNb_anno_mpra<-merge(M1_DEG_IFNb_anno_mpra, M1_DES_IFNb, by ="ESid")
colnames(des_in) <- paste0(colnames(des_in), "_DES")




mpra_in<-read_tsv("MPRA/Input/All_Monocytes_DESites_3UTR_Model1.tsv")
# 3067
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")

out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_raj_dis<-out_raj%>%dplyr::filter(grepl("*165",cohorts))
out_raj_dis<-na.omit(out_raj_dis)

out_raj_dis<-out_raj_dis%>%dplyr::select(c(ESid, cohorts, logFC, pval, gene, significant))
mpra_ifn<-out_raj_dis%>%dplyr::filter(grepl("*IFN*", cohorts))
write_tsv(mpra_ifn, "MPRA/MPRA_MFE_Discovery_IFNb_Genes.tsv")



################################################################################
# trying MFE
################################################################################
prep_stim_ed<-function(ed_in, gene_freq, max_flank){
  
  ed<-ed_in
  
  ed$Location <- gsub("^chr[0-9XY]+:(\\d+):.*", "\\1", ed$ESid)
  Gene_counts <- table(ed$gene)
  ed$Gene_Frequency <- Gene_counts[ed$gene]
  ed$Gene_Frequency<-as.numeric(ed$Gene_Frequency)
  ed$Location<-as.numeric(ed$Location)
  ed$Gene<-ed$gene
  ed<- ed %>%
    mutate(Editing_Index = sub(".*:(.*:.*)$", "\\1", ESid))
  ed<-ed%>%dplyr::select(c(Gene, Gene_Frequency, Location, ESid, Editing_Index))#, everything()))
  #ed<-ed%>%dplyr::filter(!grepl("Name=*", Gene))
  
  # maybe just to more than 2
  ed<-ed%>%dplyr::filter(ed$Gene_Frequency > 2)
  
  location_sum <- aggregate(Location ~ Gene, data = ed, 
                            FUN = function(x) c(min = min(x), max = max(x)))
  
  # Print
  for (i in 1:nrow(location_sum)) {
    flank_size <- location_sum$Location[i, "max"] - location_sum$Location[i, "min"]
    if (flank_size < max_flank) {
      cat("For category", location_sum$Gene[i], "in col 'Gene', the minimum location is", 
          location_sum$Location[i, "min"], "and the maximum location is",
          location_sum$Location[i, "max"], ". This is a flank of size", flank_size, "\n")
    }
  }
  
  # Subset lps_ed to keep only the corresponding genes less then 2000 flank
  selected_genes <- location_sum$Gene[location_sum$Location[, "max"] - location_sum$Location[, "min"] < 2000]
  ed <- ed %>% dplyr::filter(Gene %in% selected_genes)
  
  ed_bed<-ed%>%dplyr::select(c(Gene, ESid, Location, Editing_Index, Gene_Frequency))
  ed_bed
}


mpra_dis_ifn<-prep_stim_ed(mpra_ifn,2,2000)


