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

args <-commandArgs(trailingOnly = TRUE)
input_file<-args[1]

diff_in<-read_tsv(input_file)
diff_df<-diff_in%>%dplyr::select(c(GENEID, GENENAME, DE_Direction))

################################################################################
################################## GP ENRICH ################################### 
################################################################################
# Getting up GENENAME and down GENENAME seperately 
up<-diff_df%>%dplyr::filter(grepl("UP", DE_Direction))
down<-diff_df%>%dplyr::filter(grepl("DOWN", DE_Direction))
  
message("Getting  Multi-GP")
multi_gp = gost(list("up-regulated" = up$GENENAME, 
                     "down-regulated" = down$GENENAME),
                    multi_query = FALSE, evcodes = TRUE)

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
write_tsv(gp_mod, "Result/gp_enrich_result_trial.tsv")

gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod)
gp_mod_enrich  <- new("enrichResult", result = gp_mod)

################################################################################ 
########################## GP ENRICH PLOT SETTING ############################## 
################################################################################
text_size=9
legend_text_size=8
legend_size=5
number_category=15

message("Making Plots")
################################################################################
############################## GP ENRICH PLOT ##################################
################################################################################
gp_mod_enrich_up<-gp_mod_enrich%>%dplyr::filter(grepl("up-regulated",Cluster))
gp_barplot_up<-barplot(gp_mod_enrich_up, showCategory = number_category, font.size = text_size) + 
  ggplot2::ylab("Intersection size")+
  labs(x="Count",y ="Intersection Size", title="UP-regulated")+theme_bw()+
  theme(axis.text.y = element_text(color = "black", size = text_size),
        axis.title = element_text(size = 8),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=legend_text_size),
        legend.text = element_text(size=legend_text_size),
        plot.title = element_text(size=text_size,face = "bold",hjust=0.5))+ 
  ggeasy::easy_rotate_y_labels(angle = 30, side = c("right"))+ 
  guides(color = guide_legend(override.aes = list(size=legend_size)))

gp_mod_enrich_down<-gp_mod_enrich%>%dplyr::filter(grepl("down-regulated",Cluster))
gp_barplot_down<-barplot(gp_mod_enrich, showCategory = number_category,  font.size = text_size) + 
  ggplot2::ylab("Intersection size")+
  labs(x="Count",y ="Intersection Size", title="DOWN-regulated")+theme_bw()+
  theme(axis.text.y = element_text(color = "black", size = text_size),
        axis.title = element_text(size = 8),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=legend_text_size),
        legend.text = element_text(size=legend_text_size),
        plot.title = element_text(size=text_size,  face = "bold",hjust=0.5))+ 
  ggeasy::easy_rotate_y_labels(angle = 30, side = c("right"))+ 
  guides(color = guide_legend(override.aes = list(size=legend_size)))

gp_dotplot<-enrichplot::dotplot(gp_mod_enrich,x = "GeneRatio",
                                showCategory = number_category, font.size=text_size)+theme_bw()+
  labs( title="GeneRatio")+
  theme(axis.text.y = element_text(color = "black", size = text_size),
        axis.title = element_text(size = 8),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=legend_text_size),
        legend.text = element_text(size=legend_text_size),
        plot.title = element_text(size=text_size,  face = "bold",hjust=0.5))+
  ggeasy::easy_rotate_y_labels(angle = 30, side = c("right"))+ 
  guides(color = guide_legend(override.aes = list(size=legend_size)))

gp_plots<-ggarrange(gp_dotplot, gp_barplot_up,gp_barplot_down, ncol=3)
gp_plots<-annotate_figure(gp_plots, top=text_grob("Gene Enrichment", face = "bold", size =text_size))

ggsave(plot=gp_plots,filename="Result/Figure/GP_Enrich_DOT_BAR_Trial.jpg",
       width = 12, height = 8,dpi = 600)
message("Finish running GP_Enrich.R")

################################################################################
############################## GP ENRICH PLOT ##################################
################################################################################

