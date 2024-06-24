# Discovery 
# I can probably just delete this script

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


