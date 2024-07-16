library(dplyr)
library(readr)
library(tidyr)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract command line arguments
input_folder <- args[1]
output_folder <- args[2]
threshold <- as.numeric(args[3])

# Print the values of command line arguments (for demonstration)
cat("Input Folder:", input_folder, "\n")
cat("Output Folder:", output_folder, "\n")
cat("Threshold:", threshold, "\n")

##############################################################################
##############################################################################

dtu <- read_tsv(file.path(input_folder,"all_sites_pileup_dtu.tsv.gz"))
dtu <- dplyr::mutate(dtu, isoform_id = paste0(ESid, ":", allele) ) %>%
  mutate(gene_id = ESid) %>%
  dplyr::select(-ESid, -allele) %>%
  dplyr::select(gene_id, isoform_id, everything())
dtu[is.na(dtu)] <- 0
write_tsv(dtu,file.path(output_folder,"Processed_DTU.tsv"))

### Editing 
editing<-read_tsv(file.path(input_folder,"all_sites_pileup_editing.tsv.gz"))
editing[c('chr','num','Editing_Index')]<-str_split_fixed(editing$ESid, ":",3)
editing<-subset(editing, select = -c(chr,num))%>%
  dplyr::select(ESid,Editing_Index, everything())
editing[is.na(editing)]<-0
editing<-editing%>%column_to_rownames("ESid")
editing$Mean<-rowMeans(editing[,-1])
editing<-editing%>%dplyr::filter(Mean > editing_mean_threahold)
editing<-editing%>%rownames_to_column("ESid")%>%
  dplyr::select(ESid,Editing_Index,Mean, everything())

message("Writing Processed_Editing_All")
write_tsv(editing,file.path(output_folder,"Processed_Editing_All.tsv"))


### Annotation by Editing sites 
anno <- read_tsv(file.path(input_folder,"all_sites_pileup_annotation.tsv.gz"))
anno <- anno %>%
  mutate(Editing_Index = paste0(Ref, ":", Alt))%>%
  dplyr::select(-contains("AF") )%>%
  dplyr::select(ESid, ESid2, Editing_Index, Gene = Gene.refGene, 
                Location = Func.refGene, 
                Mutation = ExonicFunc.refGene, REDIportal_info,rmsk,)%>%
  mutate(known_a_i = REDIportal_info != ".") %>%
  mutate(rep_type = case_when(
    grepl("Alu", rmsk) ~ "Alu",grepl("\\=L1", rmsk) ~ "LINE",
    grepl(")n", rmsk) ~ "Simple repeat",rmsk == "." ~ "None",
    TRUE ~ "Other"
  ))%>%
  mutate(Function= case_when(grepl("ncRNA", Location) ~ "ncRNA",
                             TRUE ~ Location
  ))
anno<-anno%>%dplyr::filter(ESid2 %in% editing$ESid)
anno$Mutation[anno$Mutation == '.'] <- 'unknown'
message("Writing Processed_Annotation_All (ESid2 is the Real ESid)")
write_tsv(anno,file.path(output_folder,"Processed_Annotation_All.tsv"))

### Annotation by Index 
anno_type <- anno %>% dplyr::filter(grepl('A:G|T:C|C:T|G:A',Editing_Index))
anno_type<- anno_type%>%mutate(Editing_Index =  case_when(
  Editing_Index == 'A:G'~ 'A:G',Editing_Index == 'T:C' ~'A:G',
  Editing_Index == 'C:T'~ 'C:T', Editing_Index == 'G:A' ~ 'C:T'))

message("Writing Processed_Annotation_Type (A2G and C2T + conversion) (ESid2 is the Real ESid)")
write_tsv(anno_type,file.path(output_folder,"Processed_Annotation_Type.tsv"))

### Editing by Index
editing<-read_tsv("Local_Pipeline_Result/Patani_Editing_All.tsv")
editing_type <-editing %>% dplyr::filter(grepl('A:G|T:C|C:T|G:A',Editing_Index))%>%
  mutate(Editing_Index = case_when(
    Editing_Index == 'A:G'~ 'A:G',Editing_Index == 'T:C' ~'A:G',
    Editing_Index == 'C:T'~ 'C:T', Editing_Index == 'G:A' ~ 'C:T'))

message("Writing Processed_Editing_Type (A2G and C2T + conversion)")
write_tsv(editing_type,file.path(output_folder,"Processed_Editing_Type.tsv"))

# ## Index-Condition
# editing_ag_condition<-editing_type%>%dplyr::filter(grepl("A:G", Editing_Index))
# ##
# make_editing_condition<-function(editing_type_in, tsv_title){
#   editing_condition<-editing_type_in%>% dplyr::select(-c("Editing_Index","Mean","ESid"))
#   editing_condition<-as.data.frame(t(editing_condition))
#   editing_condition<-editing_condition%>%rownames_to_column("Sample")
#   meta_treatment<-meta%>%dplyr::select(c(Sample, sample_treatment))
#   editing_condition<-merge(editing_condition, meta_treatment, by ="Sample")
#   editing_condition<-editing_condition%>%dplyr::select(-c("Sample"))
#   write_tsv(editing_condition,tsv_title)
# }
# 
# make_editing_condition(editing_ag_condition,"Local_Pipeline_Result/AG_Editing_Condition.tsv" )