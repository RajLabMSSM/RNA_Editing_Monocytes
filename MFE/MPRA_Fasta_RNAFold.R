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
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)


###################################################################################
###################################################################################
gene_strand<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/MFE/gencode.v38.primary_assembly.gene_strand.tsv")
gene_strand$Gene<-gene_strand$gene_name
gene_strand<-gene_strand%>%dplyr::select(c(Gene, strand))
###################################################################################
###################################################################################
# STEP 1: PREP Editing differential toptable
###################################################################################

# Taking Editing Differential toptable and getting
# Genes that are edited more than twice within 2000 pair seq
change_esid12_anno<-function(anno_in){
  anno_in<-anno_in%>%dplyr::select(-c(ESid))
  anno_in$ESid<-anno_in$ESid2
  anno_in<-anno_in%>%dplyr::select(-c(ESid2))
  anno_in
}

prep_stim_ed<-function(ed_in){
  ed<-ed_in
  ed$Location <- gsub("^chr[0-9XY]+:(\\d+):.*", "\\1", ed$ESid)
  Gene_counts <- table(ed$Gene)
  ed$Gene_Frequency <- Gene_counts[ed$Gene]
  ed$Gene_Frequency<-as.numeric(ed$Gene_Frequency)
  ed$Location<-as.numeric(ed$Location)
  ed<-ed%>%dplyr::select(c(Gene, Gene_Frequency, Location, ESid, Editing_Index, everything()))
  ed<-ed%>%dplyr::filter(!grepl("Name=*", Gene))
 
  ed_bed<-ed%>%dplyr::select(c(Gene, ESid, Location, Editing_Index, Gene_Frequency))
  ed_bed
}

# Bed files contains all the sites for genes
# If a gene is edited 5 times, all 5 rows of editing site information are in this
###################################################################################
mpra_out<-read_tsv("MPRA/Result/monocyte_MPRA_results.txt")

out_raj<-mpra_out%>%dplyr::filter(!grepl("*control",category))
out_raj_dis<-out_raj%>%dplyr::filter(grepl("*165",cohorts))
out_raj_dis<-na.omit(out_raj_dis)

out_raj_dis<-out_raj_dis%>%dplyr::select(c(ESid, cohorts, logFC, pval, gene, significant))
mpra_ifn<-out_raj_dis%>%dplyr::filter(grepl("*IFN*", cohorts))
mpra_ifn<-mpra_ifn%>%dplyr::filter(grepl("Y", significant))
mpra_ifn$Gene<-mpra_ifn$gene
mpra_ifn <-mpra_ifn %>%
  mutate(Editing_Index = str_extract(ESid, "(?<=:)[^:]+:[^:]+$"))
mpra_dis_ifn<-mpra_ifn%>%dplyr::select(-c(gene))
write_tsv(mpra_dis_ifn, "MPRA/MPRA_MFE_Discovery_IFNb_Genes.tsv")
###################################################################################
mpra_dis_ifn<-read_tsv("MPRA/MPRA_MFE_Discovery_IFNb_Genes.tsv")
mpra_bed<-prep_stim_ed(mpra_dis_ifn)
write_tsv(mpra_bed, "MFE/Bed_files/MPRA_Disc_IFNb_Seperate_Ed_Sites_Location.tsv")
###################################################################################
mpra_lps<-out_raj_dis%>%dplyr::filter(grepl("*LPS*", cohorts))
mpra_lps<-mpra_lps%>%dplyr::filter(grepl("Y", significant))
mpra_lps$Gene<-mpra_lps$gene
mpra_lps <-mpra_lps %>%
  mutate(Editing_Index = str_extract(ESid, "(?<=:)[^:]+:[^:]+$"))
mpra_dis_lps<-mpra_lps%>%dplyr::select(-c(gene))
write_tsv(mpra_dis_lps, "MPRA/MPRA_MFE_Discovery_LPS_Genes.tsv")
###################################################################################
mpra_dis_lps<-read_tsv("MPRA/MPRA_MFE_Discovery_LPS_Genes.tsv")
mpra_bed<-prep_stim_ed(mpra_dis_lps)
write_tsv(mpra_bed, "MFE/Bed_files/MPRA_Disc_LPS_Seperate_Ed_Sites_Location.tsv")
###################################################################################
###################################################################################
# STEP 2: MAKE BED files from previous PREP editing toptable
###################################################################################
# Making bed files from previous Editing metrix to be put into getFasta in minerva
# Considering gene stranding problem 

make_bed<-function(bed_in,rough_bed_name,final_bed_name){
  # Gives the fisrt and last editing sites of the gene
  #bed_in<-"MFE/Bed_files/MPRA_Disc_IFNb_Seperate_IFNb_Ed_Sites_Location.tsv"
  bed_df<-read_tsv(bed_in)
  # the location of the editing is at 151th
  bed_df$Start<-(bed_df$Location - 151)
  bed_df$End<-(bed_df$Location + 150)
  bed_df$flank_size<- bed_df$End - bed_df$Start
  
  bed_df$Chr <- sub("^(chr[0-9XY]+):.*", "\\1", bed_df$ESid)
  
  
  bed_df<-bed_df%>%dplyr::select(c(Gene, ESid,Chr,Start, End,
                                    flank_size, 
                                    Gene_Frequency,Editing_Index))
  bed_df$Name<-bed_df$ESid
  bed_df$Base_1<-c('.')
  bed_df$Base_2<-c('.')
  
  bed_strand<-merge(bed_df, gene_strand, by ="Gene")
  
  # considering stranding with gene info 
  bed_strand$Edited_strand<-case_when(
    bed_strand$Editing_Index == "A:G" & bed_strand$strand == "-" ~ "-",
    bed_strand$Editing_Index == "C:T" & bed_strand$strand == "-" ~ "-",
    bed_strand$Editing_Index == "A:G" & bed_strand$strand == "+" ~ "+",
    bed_strand$Editing_Index == "C:T" & bed_strand$strand == "+" ~ "+",
    
    bed_strand$Editing_Index == "G:A" & bed_strand$strand == "-" ~ "+",
    bed_strand$Editing_Index == "T:C" & bed_strand$strand == "-" ~ "+",
    bed_strand$Editing_Index == "G:A" & bed_strand$strand == "+" ~ "-",
    bed_strand$Editing_Index == "T:C" & bed_strand$strand == "+" ~ "-")
  
  rough_bed<-bed_strand%>%dplyr::select(c(Gene,Chr, Start, End, Name, Edited_strand))
  write_tsv(rough_bed,rough_bed_name)
  
  final_bed<-bed_strand%>%dplyr::select(c(Chr, Start, End, Base_1, Base_2, Edited_strand))
  write.table(final_bed,final_bed_name, 
              sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
}

###################################################################################
make_bed("MFE/Bed_files/MPRA_Disc_IFNb_Seperate_Ed_Sites_Location.tsv", 
         "MFE/Bed_files/All_MPRA_IFNb_Ed_Sites_rough.tsv",
         "MFE/Bed_files/All_MPRA_IFNb_Ed_Sites.bed")

rough_mpra_ifnb_bed<-read_tsv("MFE/Bed_files/All_MPRA_IFNb_Ed_Sites_rough.tsv")

make_bed("MFE/Bed_files/MPRA_Disc_LPS_Seperate_Ed_Sites_Location.tsv", 
         "MFE/Bed_files/All_MPRA_LPS_Ed_Sites_rough.tsv",
         "MFE/Bed_files/All_MPRA_LPS_Ed_Sites.bed")

rough_mpra_lps_bed<-read_tsv("MFE/Bed_files/All_MPRA_LPS_Ed_Sites_rough.tsv")
###################################################################################
###################################################################################

# must be upload the _.bed files to the minerva and run 
# bedtools getfasta -fi /sc/arion/projects/ad-omics/data/references/GRCh38_references/GRCh38.primary_assembly.genome.fa
# -bed All_Esid_DES_DEG_Discovery_LPS.bed 
# -fo All_Esid_DES_DEG_Discovery_LPS.fa -s

###################################################################################
###################################################################################
# STEP 3: FOR RAW EDITING FILE, just go to step # RNAFold in Minerva 
###################################################################################
# STEP 3: FOR MULTI EDITED EDITING FILE, following below to make multi edited fasta files
# STEP 3: Making MULTI EDITED fasts sequence 
# download generated fa files back into the directory from minerva
###################################################################################
###################################################################################


# Reading in generated Fasta fa files with support data (rough bed) files
make_fasta_comp<-function( fasta_in,rough_bed_in,seperate_location_info,final_df_name){

  fasta_lines <- readLines(fasta_in)
  sequences <- list()
  current_seq_name <- NULL
  current_seq <- ""

  for (line in fasta_lines) {
    # line starts with '>'
    if (substr(line, 1, 1) == ">") {
      # when new sequence header is encountered, store the previous sequence
      if (!is.null(current_seq_name)) {
        sequences[[current_seq_name]] <- current_seq
      }
      # Extract the sequence - following row
      current_seq_name <- substr(line, 2, nchar(line))
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }
  sequences[[current_seq_name]] <- current_seq
  
  fasta_df <- data.frame(seq = unlist(sequences), row.names = names(sequences))
  
  # Downstream organization 
  fasta_df<-fasta_df%>%rownames_to_column("Chr")
  fasta_df[c('Chr','Location')]<-str_split_fixed(fasta_df$Chr, ":",2)
  fasta_df[c('Start',"End")]<-str_split_fixed(fasta_df$Location, "-",2)
  
  # Making Strand col from EndSeq(-)
  split_position <- str_split_fixed(fasta_df$End, "\\(", n = 2)
  fasta_df$End <- split_position[,1]
  
  fasta_df$Strand <- paste0("(", split_position[,2])
  fasta_df$Fasta_seq<-fasta_df$seq
  fasta_df<-fasta_df%>%dplyr::select(c(Chr, Start, End, Strand, Fasta_seq))
  ###########FLIPPING the fasta seq if the gene is - strand!!!!!1 
  fasta_df$Fasta_seq <- ifelse(fasta_df$Strand == '(-)', 
                               sapply(strsplit(fasta_df$Fasta_seq, ""), 
                                      function(x) paste0(rev(x), collapse = "")),
                               fasta_df$Fasta_seq)
  
  # Getting Gene name information from rough_bed 
  rough_bed<-read_tsv(rough_bed_in)
 
  rough_bed<-rough_bed%>%dplyr::select(c(Start, End,Gene, Name))
  merged_fasta<- merge(fasta_df, rough_bed, by = c("Start","End"), all = TRUE)
  merged_fasta$Start<-as.numeric(merged_fasta$Start)
  merged_fasta$End<-as.numeric(merged_fasta$End)
  
  merged_fasta<-merged_fasta%>%dplyr::select(c("Gene","Name","Chr",'Start',"End",
                                               "Strand","Fasta_seq"))
  

  bed_in<-read_tsv(seperate_location_info)

  bed_in$Editing_Index<-NULL
  bed_in[c("Chr","Location","Raw_Index","Edited_Index")]<-str_split_fixed(bed_in$ESid,":",4)
  bed<-bed_in%>%dplyr::select(c( ESid,Location, Raw_Index, Edited_Index,Gene_Frequency))
  bed$Location<-as.numeric(bed$Location)
  bed<-bed%>%distinct(ESid, .keep_all = TRUE)
  bed$Name<-bed$ESid
  
  merged_fasta_gene<-merge(merged_fasta, bed, by ="Name", .keep_all=TRUE)
  # for MPRA
  #merged_fasta_gene<-merged_fasta_gene%>%distinct(Location, .keep_all=TRUE)
  merged_fasta_gene$Edited_seq_num<-merged_fasta_gene$Location-merged_fasta_gene$Start
  
  merged_fasta_gene<-merged_fasta_gene%>%dplyr::select(c(Gene, Chr,Start,End,
                                                         ESid,Location, Raw_Index,Edited_Index,
                                                         Edited_seq_num, Gene_Frequency,
                                                         Strand, Name,Fasta_seq))

  write.table(merged_fasta_gene, final_df_name, 
              sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
  }

###################################################################################
###################################################################################

make_fasta_comp("MFE/Raw_fasta/All_MPRA_IFNb_Ed_Sites.fa",
                "MFE/Bed_files/All_MPRA_IFNb_Ed_Sites_rough.tsv",
                "MFE/Bed_files/MPRA_Disc_IFNb_Seperate_Ed_Sites_Location.tsv",
                "MFE/Raw_fasta/All_MPRA_IFNb_Ed_Sites_multi.tsv")
mpra_multi<-read_tsv("MFE/Raw_fasta/All_MPRA_IFNb_Ed_Sites_multi.tsv")

make_fasta_comp("MFE/Raw_fasta/All_MPRA_LPS_Ed_Sites.fa",
                "MFE/Bed_files/All_MPRA_LPS_Ed_Sites_rough.tsv",
                "MFE/Bed_files/MPRA_Disc_LPS_Seperate_Ed_Sites_Location.tsv",
                "MFE/Raw_fasta/All_MPRA_LPS_Ed_Sites_multi.tsv")
mpra_multi<-read_tsv("MFE/Raw_fasta/All_MPRA_LPS_Ed_Sites_multi.tsv")
###################################################################################
###################################################################################
# STEP 4: Checking Raw Index = Original Index in Generated Fasta 
###################################################################################
##################################################################################
check_sequence_match_again <- function(df) {
  for (gene in unique(df$Gene)) {
    group_df <- subset(df, Gene == gene)
    
    for (row in 1:nrow(group_df)) {
      editing_seq_num <- group_df$Edited_seq_num[row]
      raw_editing_seq <- group_df$Raw_Index[row]
      if (substr(group_df$Fasta_seq[row], editing_seq_num, editing_seq_num) == raw_editing_seq) {
        print(paste("Raw index matches raw fasta seq at", editing_seq_num, "th element for gene", gene))
      } else {
        print(paste("Mismatch found at", editing_seq_num, "for gene", gene))
      }
    }
  }
}


check_sequence_match <- function(df) {
  passed_rows <- data.frame()  # Initialize an empty dataframe to store passed rows
  for (gene in unique(df$Gene)) {
    group_df <- subset(df, Gene == gene)
    
    for (row in 1:nrow(group_df)) {
      editing_seq_num <- group_df$Edited_seq_num[row]
      raw_editing_seq <- group_df$Raw_Index[row]
      
      if (substr(group_df$Fasta_seq[row], editing_seq_num, editing_seq_num) == raw_editing_seq) {
        passed_rows <- rbind(passed_rows, group_df[row, ])  # Add the passed row to the dataframe
      }
    }
  }
  
  return(passed_rows)
}

check_sequence_match_again(mpra_multi)
# Everything passes the checking after flipping seq read for - strand genes 
mpra_match<-check_sequence_match(mpra_multi)
write_tsv(mpra_match,"MFE/Edited_fasta/Edited_All_MPRA_LPS_Ed_Sites.tsv")
check_sequence_match_again(mpra_match)

##################################################################################
# SUMMARY 
# Threshold 
# Gene must be edited more then Twice within 2000 base seq on the read
# the Raw_Index (pre-editing base) must match with Original fasta seq's base 
# there are like 10 gene~ mismatch.. for each LPS and IFNb
# LPS : There are 97 unique Gene matched with 97 unique sites (input fasta seq and passed all threshold)
# IFNb : There are 243 unique Gene matched with 243 unique sites (input fasta seq and passed all threshold)
###################################################################################
##################################################################################
# STEP 5 : on Python, make MULTI EDITED FASTA Seq metrix 
###################################################################################
# MFE_Multi_Editing.ipynb
# Using Matched fa files as input 
##################################################################################
##################################################################################
# STEP 6 : Make FASTA.fa format from MULTI EDITING python output
###################################################################################
# Making Fasta.fa format to be uploaded to Minerva \

make_multi_edit_fasta<-function(multi_edit_in, support_name, file_path_name){

  #df<-read.table("MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites.tsv", sep = ",", header = TRUE)
  df<-read.table(multi_edit_in, sep = ",", header = TRUE)
  df[, 1] <-NULL
  df$Strand<-c("+")
  df$Name<-paste0(">",df$Chr,":",df$Start,"-",df$End,df$Strand)
  
  support<- write_tsv(df,support_name)
  
  fasta<-df%>%select(c(Name,Edited_seq))#%>%
    #distinct(Name,.keep_all=TRUE)
  t_fasta <- t(fasta)
  
  # Write fasta format 
  file_path <- file_path_name
  file_conn <- file(file_path, "w")
  for (i in 1:length(t_fasta)) {
    writeLines(as.character(t_fasta[i]), file_conn)
    }
  close(file_conn)
  }

make_multi_edit_fasta("MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites.tsv",
                      "MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites_support.tsv",
                      "MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites.fa")

make_multi_edit_fasta("MFE/Edited_fasta/Edited_All_MPRA_LPS_Ed_Sites.tsv",
                      "MFE/Edited_fasta/Edited_All_MPRA_LPS_Ed_Sites_support.tsv",
                      "MFE/Edited_fasta/Edited_All_MPRA_LPS_Ed_Sites.fa")
################################################################################
################################################################################

##################################################################################
##################################################################################
# STEP 7 : MINERVA RNAFold 
# for RAW Fasta, skip step 3,4,5,6 and do this step 
# for Edited Fasts, must go through step 3,4,5,6 and upload generated multi edited fasts to minerva
##################################################################################
# module load viennarna/2.4.12
# RNAfold --infile=test_fasta.fa --outfile=test_fasta_out.txt 
# --jobs=20 --partfunc --noPS --noDP
###################################################################################
##################################################################################
# STEP 7 : RNAFold output Organization 
##################################################################################
# Organize RNAfold output txt to be readable , and add support data

make_rnafold_organized<-function(txt_in, support_in,file_name){
  #txt_in<-"MFE/RNAfold/Edited_All_MPRA_IFNb_Ed_Sites_RNAfold.txt"
  rnafold_out <- readLines(txt_in)
  print(head(rnafold_out,8))
  
  filtered_lines <- rnafold_out[seq(1, length(rnafold_out), by = 6)]
  
  rnafold_df <- data.frame(
    Name = filtered_lines,
    Seq = rnafold_out[seq(2, length(rnafold_out), by = 6)],
    Str = rnafold_out[seq(3, length(rnafold_out), by = 6)]
  )
  
  rnafold_df$Name<-gsub(">","",rnafold_df$Name)
  rnafold_df$Name <- gsub("\\(\\+\\)|\\(\\-\\)", "", rnafold_df$Name)
  
  rnafold_df[c("Chr", "Location")]<-str_split_fixed(rnafold_df$Name,":",2)
  rnafold_df[c("Start", "Location")]<-str_split_fixed(rnafold_df$Location,"-",2)
  rnafold_df$End <- as.numeric(gsub("\\D", "", rnafold_df$Location))
  rnafold_df$Stranding <- gsub("\\d", "", rnafold_df$Location)
  # GPNMB
  #rnafold_df$Stranding <-c("+")
  
  extract_num <- function(str) {
    last_open_index <- max(gregexpr("\\(", str)[[1]])
    last_close_index <- max(gregexpr("\\)", str)[[1]])
    if (last_open_index > 0 && last_close_index > 0 && last_close_index > last_open_index) {
      return(substr(str, last_open_index + 1, last_close_index - 1))
    } else {
      return(NA)
    }
  }
  
  rnafold_df$MFE <- sapply(rnafold_df$Str, extract_num)
  rnafold_df$Start<-as.numeric(rnafold_df$Start)
  rnafold_df$Flank_size<-rnafold_df$End-rnafold_df$Start
  rnafold_df$MFE<-as.numeric(rnafold_df$MFE)
  rnafold_df$Name<-gsub("\\+$","",rnafold_df$Name)
  
  #support_in<-"MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites_support.tsv"
  support<-read_tsv(support_in)
  # to get the name 
  #support<-support%>%distinct(Name, .keep_all=TRUE)
  support<-support%>%select(c(Gene, Gene_Frequency,Name))
  support$Name<-gsub(">","",support$Name)
  support$Name<-gsub("\\+$","",support$Name)
  #support<-support%>%
  #  mutate(Name = substr(Name,1,nchar(Name)-1))
  #support<-support%>%select(c(Gene, Gene_Frequency,Name))
  
  rnafold_support<-merge(rnafold_df, support, by="Name")%>%
    distinct(Name, .keep_all=TRUE)
  
  rnafold_support<-rnafold_support%>%dplyr::select(c(Gene,Gene_Frequency, Flank_size, 
                                           MFE,Name,Chr,Start,End,Seq,Str))
  write_tsv(rnafold_support,file_name)
}

###################################################################################
make_rnafold_organized("MFE/RNAfold/All_MPRA_IFNb_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites_support.tsv",
                       "MFE/RNAfold/All_MPRA_IFNb_Ed_Sites_RNAfold_Organized.tsv")

raw_mpra<-read_tsv("MFE/RNAfold/All_MPRA_IFNb_Ed_Sites_RNAfold_Organized.tsv")
################################################################################

make_rnafold_organized("MFE/RNAfold/Edited_All_MPRA_IFNb_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Edited_All_MPRA_IFNb_Ed_Sites_support.tsv",
                       "MFE/RNAfold/Edited_All_MPRA_IFNb_Ed_Sites_RNAfold_Organized.tsv")

edited_mpra<-read_tsv("MFE/RNAfold/Edited_All_MPRA_IFNb_Ed_Sites_RNAfold_Organized.tsv")
################################################################################
make_rnafold_organized("MFE/RNAfold/All_MPRA_LPS_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Edited_All_MPRA_LPS_Ed_Sites_support.tsv",
                       "MFE/RNAfold/All_MPRA_LPS_Ed_Sites_RNAfold_Organized.tsv")

raw_mpra<-read_tsv("MFE/RNAfold/All_MPRA_LPS_Ed_Sites_RNAfold_Organized.tsv")
################################################################################

make_rnafold_organized("MFE/RNAfold/Edited_All_MPRA_LPS_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Edited_All_MPRA_LPS_Ed_Sites_support.tsv",
                       "MFE/RNAfold/Edited_All_MPRA_LPS_Ed_Sites_RNAfold_Organized.tsv")

edited_mpra<-read_tsv("MFE/RNAfold/Edited_All_MPRA_LPS_Ed_Sites_RNAfold_Organized.tsv")
################################################################################

# Checking 
# 6th row 109th element gin Edited_seq
edited_nth_element <- substr(edited_limk_ifnb$Edited_seq[6], 109, 109)
raw_nth_element <- substr(edited_limk_ifnb$Raw_seq[6], 109, 109)
print(nth_element)
print(raw_nth_element)
################################################################################
# STEP 8 : Make Delta MFE metrics (Raw MFE - Edited MFE for all matched gene and site)
################################################################################
make_delta_mfe<-function(raw_in, edited_in,cohort_in){
  #raw_in<-raw_lps
  #edited_in<-edited_lps
  raw<- raw_in[order(raw_in$Chr), ]
  edited <- edited_in[order(edited_in$Chr), ]
  
  raw<-raw%>%dplyr::filter(Name %in% edited$Name)%>%
    distinct(Name, .keep_all=TRUE)
  edited<-edited%>%
    distinct(Name, .keep_all=TRUE)
  
  common<-raw%>%dplyr::select(c(Gene, Gene_Frequency,Flank_size,Name,
                                Chr,Start,End))

  colnames(raw) <- paste(colnames(raw), "Raw", sep = "_")
  colnames(edited) <- paste(colnames(edited), "Edited", sep = "_")
  
  all<-cbind(common,raw, edited)
  
  # Making Delta MFE Raw-Edited
  all$MFE_Delta<-  (all$MFE_Raw - all$MFE_Edited)
  all$Cohort<-c(cohort_in)
  all<-all%>%dplyr::select(c(Gene, Gene_Frequency,Flank_size,Name,
                             Chr,Start,End,MFE_Delta,MFE_Raw,MFE_Edited,
                             Cohort,Seq_Raw,Seq_Edited))
  all
}


delta_mpra<-make_delta_mfe(raw_mpra, edited_mpra,"IFNb_MPRA")
write_tsv(delta_mpra, "MFE/Delta_MFE/MPRA_IFNb_Delta_MFE.tsv")

delta_mpra<-make_delta_mfe(raw_mpra, edited_mpra,"LPS_MPRA")
write_tsv(delta_mpra, "MFE/Delta_MFE/MPRA_LPS_Delta_MFE.tsv")
################################################################################
################################################################################