BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE"))
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(gprofiler2)
library(ggplotify)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggpubr)
recoding<-read_tsv("All_Table/all_esid_gene_exonic.tsv")
recoding<-recoding%>%dplyr::filter(grepl("sig",aa_change))
re_gene<-recoding%>%distinct(Gene.refGene,.keep_all=TRUE)
# here before distinct, get the sites, seperate by cohort, merge with sites in the stimulation-mean matrix. 
# do box plots
write_tsv(re_gene, "All_Table/exonic_sig_recoding_gene.tsv")



re_gene<-read_tsv("All_Table/exonic_sig_recoding_gene.tsv")
re_gene$GENENAME<-re_gene$Gene.refGene


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


################################################################################
################################################################################
M1_AG<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_AG_Editing.tsv")
M1_CT<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_CT_Editing.tsv")
M2_AG<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Jacusa/Processed_AG_Editing.tsv")
M2_CT<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Replication/Jacusa/Processed_CT_Editing.tsv")


recoding<-read_tsv("All_Table/all_esid_gene_exonic.tsv")
recoding<-recoding%>%dplyr::filter(grepl("sig",aa_change))
recoding<-recoding%>%distinct(ESid, .keep_all=TRUE)
################################################################################
columns_to_keep <- c("Sample", "Stimulation", "Mean")
################################################################################
comparison_dis <- list( c("Basal", "IFNb"), c("IFNb", "LPS"), c("Basal", "LPS"))
comparison_rep <- list( c("B", "IFNg"), c("IFNg", "LPS"), c("B", "LPS"))

fun_mean <- function(x){
  return(data.frame(y=round(mean(x),3),label=mean(round(x,3),na.rm=T)))}

box_plot<-function(df_in, x_in,y_in, fill_in,shape_in,my_comparisons,plot_title){
  box_plot<-ggplot(df_in,aes(x = x_in , y= y_in,fill = x_in))+
    geom_jitter(aes(colour = x_in,alpha = 0.7, width = 0.25, size = point_size,shape=shape_in)) +
    scale_fill_manual(values = stim_color)+
    geom_boxplot(fill = NA, outlier.shape=NA) +
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",label.x.npc = "center",
                       label.y.npc = "top") +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=20, face="bold",size=point_size, show.legend=TRUE)+
    stat_summary(aes(label=round(..y..,3)),fun.data = fun_mean, geom="text", 
                 size=inplot_size, face="bold",vjust=2)+
    labs(x = "", y = "Editing Rate", title = plot_title)+
    theme_classic()+theme(legend.title = element_text(size=text_size,color="black"),
                          legend.text = element_text(size=text_size,color="black"),
                          axis.text.x = element_text(size = axis_size,color="black"),
                          axis.text.y = element_text(size = axis_size,color="black"),
                          axis.title = element_text(size = axis_size,color="black"),
                          plot.title = element_text(size=title_size, hjust=0.5,color="black"),
                          legend.margin=margin(0,0,0,0),
                          legend.box.margin=margin(-3,-3,-3,-3))+
    theme(legend.position = "none")
  box_plot
}

resites_rate<-function( cohort_in, editing_in, comparison_in,title_in){
  editing_in<-M1_AG
  cohort_in<-"Discovery"
  
  editing_info<-editing_in%>%dplyr::select(c(Sample, Stimulation))
  recoding_select<-recoding%>%dplyr::filter(grepl( cohort_in, cohort))
  samples<-intersect(recoding_select$ESid, colnames(editing_in))
  editing_select<-editing_in%>%select(any_of(samples))
  editing_select<-merge(editing_select, editing_info, by=0)  
  
  editing_select<-editing_select%>%select(-c(Row.names))
  editing_select$Mean <- rowMeans(editing_select[, !(names(editing_select) %in% c("Stimulation", "Sample"))])
  editing_select<-editing_select%>%select(c(Sample, Stimulation,Mean, everything()))
  editing_select_mean<-aggregate(editing_select[, 4:24], list(editing_select$Stimulation), mean)
  editing_select_mean<-editing_select_mean%>%column_to_rownames("Group.1")
  editing_select_mean<-as.data.frame(t(editing_select_mean))
  editing_select_mean<-editing_select_mean%>%rownames_to_column("ESid")
  
  editing_select_mean<-editing_select_mean%>%
    pivot_longer(!ESid, names_to = "Stimulation", values_to="Mean")
  
  filt<-editing_select_mean%>%dplyr::filter()
  # editing_select<-merge(recoding_select, editing_select, by=0)
  # editing_select<-editing_select%>%select(c(Sample, Stimulation,Mean, everything()))
  # editing_select<-editing_select%>%select(-c(Row.names))
  # editing_check<-editing_select%>%dplyr::filter(!grepl("^Non_DE$",Discovery_DES_IFN))
  # editing_check<-editing_check%>%dplyr::filter(!grepl("Basal",Stimulation))
  # 
  # editing_ch
  
  #print(samples)
  plot<-box_plot(editing_select_mean, editing_select_mean$Stimulation, editing_select_mean$Mean, 
                 editing_select_mean$Stimulation,editing_select_mean$ESid,comparison_dis, "title_in")
  plot                                     
}

dis_ag<-resites_rate("Discovery",M1_AG,comparison_dis,"30 sites in M1 AG")
dis_ct<-resites_rate("Discovery",M1_CT,comparison_dis,"33 sites in M1 CT")


rep_ag<-resites_rate("Replication",M2_AG,comparison_rep,"14 sites in M2 AG")
rep_ct<-resites_rate("Discovery",M2_CT,comparison_rep,"6 sites in M2 CT")

re_rate_plots<-ggarrange(dis_ag, dis_ct, rep_ag, rep_ct, nrow=2, ncol=2)
ggsave(plot=re_rate_plots,filename="/Users/hyominseo/Desktop/RAJ_RNA_Editing/Figure/Temp/Sig_Recoding_Gene_Editing_Rate.jpg",
       width = 10, height = 10,dpi = 600)

