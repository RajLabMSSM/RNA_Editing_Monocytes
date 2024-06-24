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
###################################################################################
make_des_ed<-function( ed_in, anno_in){
  des<-read_tsv(ed_in)
  des<-des%>%dplyr::filter(!grepl("Non_DE",DE_Index))
  
  anno<-read_tsv(anno_in)
  anno<-change_esid12_anno(anno)
  anno<-anno%>%dplyr::filter(ESid %in% des$ESid)
  anno<-anno%>%dplyr::select(-c(Gene, Editing_Index))
  
  anno_des<-merge(des, anno, by ='ESid')
  anno_des
}
###################################################################################
lps_des<-make_des_ed("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Editing_Differential/LPS_1_Toptable.tsv",
                     "/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_All.tsv")

ifnb_des<-make_des_ed("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Editing_Differential/IFNb_1_Toptable.tsv",
                     "/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Jacusa/Processed_Annotation_All.tsv")

deg_check<-read_tsv("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/Toptable/Gene_Differential/IFNb_DEG_Toptable.tsv")

prep_stim_ed<-function(ed_in, gene_freq, max_flank){

  ed<-ed_in
  ed$Location <- gsub("^chr[0-9XY]+:(\\d+):.*", "\\1", ed$ESid)
  Gene_counts <- table(ed$Gene)
  ed$Gene_Frequency <- Gene_counts[ed$Gene]
  ed$Gene_Frequency<-as.numeric(ed$Gene_Frequency)
  ed$Location<-as.numeric(ed$Location)
  ed<-ed%>%dplyr::select(c(Gene, Gene_Frequency, Location, ESid, Editing_Index, everything()))
  ed<-ed%>%dplyr::filter(!grepl("Name=*", Gene))
  
  # maybe just to more than 2
  ed<-ed%>%dplyr::filter(ed$Gene_Frequency > gene_freq)
  
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

# Bed files contains all the sites for genes
# If a gene is edited 5 times, all 5 rows of editing site information are in this
lps_bed<-prep_stim_ed(lps_des,
                      2,
                      2000)
write_tsv(lps_bed, "MFE/Bed_files/Seperate_LPS_Ed_Sites_Location.tsv")

ifnb_bed<-prep_stim_ed(ifnb_des,
                      2,
                      2000)
write_tsv(ifnb_bed, "MFE/Bed_files/Seperate_IFNb_Ed_Sites_Location.tsv")

###################################################################################
###################################################################################
###################################################################################
###################################################################################
# STEP 2: MAKE BED files from previous PREP editing toptable
###################################################################################
# Making bed files from previous Editing metrix to be put into getFasta in minerva
# Considering gene stranding problem 

make_bed<-function(bed_in,rough_bed_name,final_bed_name){
  # Gives the fisrt and last editing sites of the gene
  bed_in<-"MFE/Bed_files/MPRA_Disc_IFNb_Seperate_IFNb_Ed_Sites_Location.tsv"
  bed_df<-read_tsv(bed_in)
  subset_df <- bed_df %>%
    group_by(Gene) %>%
    filter(Location == max(Location) | Location == min(Location))
  
  subset_df$End<-subset_df$Location
  subset_df$Location<- NULL
  
  # in one row per gene, start (min location), end (max location)
  summary_df <- bed_df %>%
    group_by(Gene) %>%
    summarize(Start = min(Location), End = max(Location)) %>%
    ungroup()
  
  summary_df$flank_size<-summary_df$End - summary_df$Start
  
  result_df <- merge(subset_df, summary_df, by = c("Gene", "End"))%>%
    distinct(Gene, .keep_all=TRUE)
  
  # check on this result_df for more info
  result_df$Chr <- sub("^(chr[0-9XY]+):.*", "\\1", result_df$ESid)
  result_df<-result_df%>%dplyr::select(c(Gene, ESid,Chr,Start, End,
                                         flank_size,
                                         Gene_Frequency,Editing_Index))
  bed<-result_df
  bed$Name<-result_df$ESid
  bed$Chr<-result_df$Chr
  bed$Start<- (result_df$Start - 101)
  bed$End<- (result_df$End + 100)
  bed$Base_1<-c('.')
  bed$Base_2<-c('.')
  
  bed<-bed%>%dplyr::select(c(Chr, Start, End, Name, Base_1, Base_2,Gene,Editing_Index))
  
  bed_strand<-merge(bed, gene_strand, by ="Gene")
  
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
make_bed("MFE/Bed_files/Seperate_LPS_Ed_Sites_Location.tsv", 
         "MFE/Bed_files/All_LPS_Ed_Sites_rough.tsv",
         "MFE/Bed_files/All_LPS_Ed_Sites.bed")

make_bed("MFE/Bed_files/Seperate_IFNb_Ed_Sites_Location.tsv", 
         "MFE/Bed_files/All_IFNb_Ed_Sites_rough.tsv",

# Making support files for the bed files 
# 106 sites
rough_lps_bed<-read_tsv("MFE/Bed_files/All_LPS_Ed_Sites_rough.tsv")
# 263 sites
rough_ifnb_bed<-read_tsv("MFE/Bed_files/All_IFNb_Ed_Sites_rough.tsv")


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
  bed<-bed_in%>%dplyr::select(c(Gene, ESid,Location, Raw_Index, Edited_Index,Gene_Frequency))
  bed$Location<-as.numeric(bed$Location)
  
  merged_fasta_gene<-merge(merged_fasta, bed, by ="Gene", .keep_all=TRUE)
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

make_fasta_comp("MFE/Raw_fasta/All_LPS_Ed_Sites.fa",
                "MFE/Bed_files/All_LPS_Ed_Sites_rough.tsv",
                "MFE/Bed_files/Seperate_LPS_Ed_Sites_Location.tsv",
                "MFE/Raw_fasta/All_LPS_Ed_Sites_multi.tsv")
lps_multi<-read_tsv("MFE/Raw_fasta/All_LPS_Ed_Sites_multi.tsv")

make_fasta_comp("MFE/Raw_fasta/All_IFNb_Ed_Sites.fa",
                "MFE/Bed_files/All_IFNb_Ed_Sites_rough.tsv",
                "MFE/Bed_files/Seperate_IFNb_Ed_Sites_Location.tsv",
                "MFE/Raw_fasta/All_IFNb_Ed_Sites_multi.tsv")
ifnb_multi<-read_tsv("MFE/Raw_fasta/All_IFNb_Ed_Sites_multi.tsv")

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

check_sequence_match_again(lps_multi)
# Everything passes the checking after flipping seq read for - strand genes 
lps_match<-check_sequence_match(lps_multi)
write_tsv(lps_match,"MFE/Edited_fasta/Edited_All_LPS_Ed_Sites.tsv")
check_sequence_match_again(lps_match)

check_sequence_match_again(ifnb_multi)
# Everything passes the checking after flipping seq read for - strand genes 
ifnb_match<-check_sequence_match(ifnb_multi)
write_tsv(ifnb_match,"MFE/Edited_fasta/Edited_All_IFNb_Ed_Sites.tsv")
check_sequence_match_again(ifnb_match)

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

  #df<-read.table("MFE/Edited_fasta/Multi_Edited_All_LPS_Ed_Sites.tsv", sep = ",", header = TRUE)
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

make_multi_edit_fasta("MFE/Edited_fasta/Multi_Edited_All_LPS_Ed_Sites.tsv",
                      "MFE/Edited_fasta/Multi_Edited_LPS_support.tsv",
                      "MFE/Edited_fasta/Multi_Edited_All_LPS_Ed_Sites.fa")

make_multi_edit_fasta("MFE/Edited_fasta/Multi_Edited_All_IFNb_Ed_Sites.tsv",
                      "MFE/Edited_fasta/Multi_Edited_IFNb_support.tsv",
                      "MFE/Edited_fasta/Multi_Edited_All_IFNb_Ed_Sites.fa")

################################################################################
################################################################################
# LIMK2 has positive strand
#df$Strand<-c("+")
make_multi_edit_fasta("MFE/Edited_fasta/LIMK2_Edited_All_IFNb_Ed_Sites.tsv",
                      "MFE/Edited_fasta/LIMK2_Edited_IFNb_support.tsv",
                      "MFE/Edited_fasta/LIMK2_Edited_All_IFNb_Ed_Sites.fa")

# limk_check<-read_tsv("MFE/Edited_fasta/LIMK2_Edited_All_IFNb_Ed_Sites.tsv")
# limk_check_support<-read_tsv("MFE/Edited_fasta/LIMK2_Edited_IFNb_support.tsv")
# limk_check_fa<-read_tsv("MFE/Edited_fasta/LIMK2_Edited_All_IFNb_Ed_Sites.fa")

################################################################################
#df$Strand<-c("+")
make_multi_edit_fasta("MFE/Edited_fasta/GPR141_Edited_All_IFNb_Ed_Sites.tsv",
                      "MFE/Edited_fasta/GPR141_Edited_IFNb_support.tsv",
                      "MFE/Edited_fasta/GPR141_Edited_All_IFNb_Ed_Sites.fa")

check<-read_tsv("MFE/Edited_fasta/GPR141_Edited_IFNb_support.tsv")
##################################################################################
##################################################################################
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
  #txt_in<-"MFE/RNAfold/All_MPRA_IFNb_Ed_Sites_RNAfold.txt"
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

  support<-read_tsv(support_in)
  # to get the name 
  support<-support%>%distinct(Name, .keep_all=TRUE)
  support<-support%>%select(c(Gene, Gene_Frequency,Name))
  support$Name<-gsub(">","",support$Name)
  
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
make_rnafold_organized("MFE/RNAfold/All_LPS_Ed_Sites_RNAfold.txt", 
                       "MFE/Edited_fasta/Multi_Edited_LPS_support.tsv",
                       "MFE/RNAfold/All_LPS_Ed_Sites_RNAfold_Organized.tsv")

raw_lps<-read_tsv("MFE/RNAfold/All_LPS_Ed_Sites_RNAfold_Organized.tsv")


make_rnafold_organized("MFE/RNAfold/All_IFNb_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Multi_Edited_IFNb_support.tsv",
                       "MFE/RNAfold/All_IFNb_Ed_Sites_RNAfold_Organized.tsv")

raw_ifnb<-read_tsv("MFE/RNAfold/All_IFNb_Ed_Sites_RNAfold_Organized.tsv")
################################################################################
make_rnafold_organized("MFE/RNAfold/Multi_Edited_All_LPS_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Multi_Edited_LPS_support.tsv",
                       "MFE/RNAfold/Multi_Edited_All_LPS_Ed_Sites_RNAfold_Organized.tsv")

edited_lps<-read_tsv("MFE/RNAfold/Multi_Edited_All_LPS_Ed_Sites_RNAfold_Organized.tsv")

make_rnafold_organized("MFE/RNAfold/Multi_Edited_All_IFNb_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/Multi_Edited_IFNb_support.tsv",
                       "MFE/RNAfold/Multi_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv")

edited_ifnb<-read_tsv("MFE/RNAfold/Multi_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv")
################################################################################
################################################################################

gene_make_rnafold_organized<-function(txt_in, support_in,file_name){
  #txt_in<-"MFE/RNAfold/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_RNAfold.txt"
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

  #support<-read.table("MFE/Edited_fasta/GPR141_Edited_IFNb_support.tsv",sep='\t',header=TRUE)
  support<-read.table(support_in,sep='\t',header=TRUE)
  support<-rownames_to_column(support, "row_names")
  support<-support%>%dplyr::select(c('row_names','Edited_seq_num','Edited_Index',
                                    'Raw_seq','Edited_seq'))
                                     
  rnafold_df<-rownames_to_column(rnafold_df, "row_names")
  
  combined_df<-merge(rnafold_df, support, by ='row_names')
  combined_df<-combined_df%>%dplyr::select(c(row_names,Name,Edited_seq_num,
                                             Edited_Index, MFE, Chr, Location,
                                             Start, End, Stranding,
                                             Flank_size, Raw_seq, Edited_seq,Str))
  write.table(combined_df, file = "MFE/RNAfold/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  }


make_rnafold_organized("MFE/RNAfold/Gene_specific/LIMK2_Edited_All_IFNb_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/LIMK2_Edited_IFNb_support.tsv",
                       "MFE/RNAfold/Gene_specific/LIMK2_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv")

limk_edited_ifnb<-read_tsv("MFE/RNAfold/Gene_specific/LIMK2_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv")


make_rnafold_organized("MFE/RNAfold/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_RNAfold.txt",
                       "MFE/Edited_fasta/GPR141_Edited_IFNb_support.tsv",
                       "MFE/RNAfold/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv")

gpr_edited_ifnb<-read_tsv("MFE/RNAfold/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv")

################################################################################
################################################################################


edited_limk_ifnb<-read.table("MFE/RNAfold/Gene_specific/LIMK2_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv",
                            sep="\t",header=TRUE)
edited_gpr_ifnb<-read.table("MFE/RNAfold/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_RNAfold_Organized.tsv",
                            sep="\t",header=TRUE)
# Checking 
# 6th row 109th element gin Edited_seq
edited_nth_element <- substr(edited_limk_ifnb$Edited_seq[6], 109, 109)
raw_nth_element <- substr(edited_limk_ifnb$Raw_seq[6], 109, 109)
print(nth_element)
print(raw_nth_element)

gene_make_delta_mfe<-function(raw_in,edited_in, gene_in,file_name ){
  raw_gene<-raw_in%>%dplyr::filter(grepl(gene_in, Gene))
  edited_in$MFE_Edited<-as.numeric(edited_in$MFE)
  edited_in$MFE_Raw<-c(as.numeric(raw_gene$MFE))
  edited_gene<-edited_in%>%dplyr::select(-c(MFE))
  
  edited_gene$Delta_MFE<-(edited_gene$MFE_Raw - edited_gene$MFE_Edited)
  edited_gene$Edited_Frequency<- str_count(edited_gene$Edited_Index, "[CGA]")
  
  edited_gene<-edited_gene%>%dplyr::select(c(Name, Chr, Location, Start, End, Stranding,
                                             Flank_size, Edited_seq_num, Edited_seq, 
                                             Edited_Index,MFE_Edited, MFE_Raw, Delta_MFE, 
                                             Edited_Frequency, Edited_seq, Raw_seq, Str))
  
  write.table(edited_gene, file = file_name , 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
}

gene_make_delta_mfe(raw_ifnb, edited_limk_ifnb, "LIMK2",
                    "MFE/Delta_MFE/Gene_specific/LIMK2_Edited_All_IFNb_Ed_Sites_Delta_MFE.tsv")
gene_make_delta_mfe(raw_ifnb, edited_gpr_ifnb, "GPR141",
                    "MFE/Delta_MFE/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_Delta_MFE.tsv")


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

delta_lps<-make_delta_mfe(raw_lps, edited_lps,"LPS")
write_tsv(delta_lps, "MFE/Delta_MFE/LPS_Delta_MFE.tsv")

delta_ifnb<-make_delta_mfe(raw_ifnb, edited_ifnb,"IFNb")
write_tsv(delta_ifnb, "MFE/Delta_MFE/IFNb_Delta_MFE.tsv")

delta_lps<-read_tsv("MFE/Delta_MFE/LPS_Delta_MFE.tsv")
delta_ifnb<-read_tsv("MFE/Delta_MFE/IFNb_Delta_MFE.tsv")

################################################################################
# STEP 9 : All plots 
################################################################################
inplot_size = 4
point_size= 2
axis_size = 12
ticks_size =12
title_size = 12
text_size = 12
label_text_size =15
table_text_size =5
anno_legend_size= 3
legend_size = 3
label_margin = 10

cohort_color <- c("A:G" = "#F12C17","C:T" = "#215BEC",
                  "LPS" = "#fea706","IFNb" = "#3a857e","IFNb_MPRA" = "#3a857e",
                  "GPNMB_IFNb" = "#215BEC",
                  "Outlier"="grey",
                  "MFE_Raw" = "#F12C17","MFE_Edited" = "#215BEC")
################################################################################
qc_plot <- function(df_in, Cohort, y_in,plot_title) {
  qc_plot <- ggplot(data = df_in, aes(x = Cohort, y = y_in)) +
    scale_color_manual(values = cohort_color) +
    geom_jitter(aes(color=Cohort),alpha = 0.7, width = 0.25, size = point_size) +
    geom_boxplot(fill = NA,outlier.shape=NA) +
    #stat_compare_means(label = "p.signif", label.x.npc = "center", label.y.npc = "center") +
    labs(x = "", y = plot_title)+
    #stat_summary(aes(label = round(..y.., 3)),
    #  fun.data = fun_mean,geom = "text",
    #  size = inplot_size,vjust = 2) 
    theme_classic() + theme(
      legend.title = element_text(size = text_size,color="black"),
      legend.text = element_text(size = text_size,color="black"),
      axis.text.y = element_text(size = axis_size,color="black"),
      axis.text.x = element_text(size = axis_size,color='black'),
      axis.title = element_text(size = axis_size,color="black"),
      plot.title = element_text(size = title_size, hjust = 0.5,color="black"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-3,-3,-3,-3))+
    guides(color = guide_legend(override.aes = list(size = anno_legend_size)))+
    theme(legend.position = "top")
  qc_plot
}


line_plot<-function(comp_in,x_in,y_in,
                    comp_select_in, comp_select_in_x, comp_select_in_y,x_title){
  discovery_plot<-
    ggplot(comp_in, aes(x=x_in, y=y_in,col=Cohort))+
    geom_point(size= point_size+0.5, alpha=0.7)+
    scale_color_manual(values = cohort_color)+
    geom_smooth(method="lm", level=0.98, se = FALSE, show.legend = FALSE)+
    stat_cor(method="pearson")+
    labs(x = x_title, y = "Delta MFE")+
    theme_classic() +
    #geom_vline(xintercept = (0), linetype="dashed", color = "grey", size=0.3)+
    #geom_hline(yintercept = (0), linetype="dashed", color = "grey", size=0.3)+
    theme(legend.title = element_text(size=text_size),
          legend.text = element_text(size=text_size),
          axis.text.x = element_text(size = axis_size),
          axis.text.y = element_text(size = axis_size),
          axis.title = element_text(size = axis_size),
          plot.title = element_text(size=title_size, hjust=0.5),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-3,-3,-3,-3))+
    theme(legend.position = "top")+
    guides(color = guide_legend(override.aes = list(size=legend_size)))+
    geom_point(data =comp_select_in, 
               aes(x = comp_select_in_x, y = comp_select_in_y), color="black",size = point_size+0.5)+
    geom_label_repel(data = comp_select_in,  
                     aes(x =comp_select_in_x, y =comp_select_in_y, label = Gene),
                     box.padding = 0, point.padding = 0, force = 50,
                     segment.size = 0.2,segment.color = "black",
                     angle = 180,point.size = point_size+1,color="black",
                     nudge_x = 0.3,hjust =0,size = point_size)
  discovery_plot
}

################################################################################
############ the MORE Negative MFE is the more stable it is ####################
################################################################################
delta_lps_mfe<-delta_lps%>%dplyr::select(c(Name,Cohort, MFE_Delta)) 
delta_ifnb_mfe<-delta_ifnb%>%dplyr::select(c(Name,Cohort, MFE_Delta)) 
delta_all<-rbind(delta_lps_mfe, delta_ifnb_mfe)

delta_all <-delta_all%>%
  pivot_longer(where(is.numeric),names_to="MFE_category", values_to="Delta_MFE")

################################################################################

mfe_all<-qc_plot(delta_all,delta_all$Cohort, delta_all$Delta_MFE,
                 "Raw-Edited MFE")

mfe_all

################################################################################
delta_lps_select<-delta_lps%>%
  dplyr::filter(grepl("^LIMK2$",Gene))

delta_ifnb_select<-delta_ifnb%>%
  dplyr::filter(grepl("^GPR141$",Gene))

mfe_gene_lps<-line_plot(delta_lps,delta_lps$Gene_Frequency, delta_lps$MFE_Delta,
                        delta_lps_select, 
                        delta_lps_select$Gene_Frequency, delta_lps_select$MFE_Delta,
                 "Gene Editing Frequency")

mfe_gene_ifnb<-line_plot(delta_ifnb,delta_ifnb$Gene_Frequency, delta_ifnb$MFE_Delta,
                         delta_ifnb_select, 
                         delta_ifnb_select$Gene_Frequency, delta_ifnb_select$MFE_Delta,
                    "Gene Editing Frequency")
mfe_gene<-ggarrange(mfe_gene_lps, mfe_gene_ifnb, ncol=2)

mfe_gene
################################################################################
#delta_lps$Flank_Frequency<-(delta_lps$Gene_Frequency)/(delta_lps$Flank_size)
#delta_ifnb$Flank_Frequency<-(delta_ifnb$Gene_Frequency)/(delta_ifnb$Flank_size)

delta_lps$Flank_Frequency<-(delta_lps$Flank_size)/(delta_lps$Gene_Frequency)
delta_ifnb$Flank_Frequency<-(delta_ifnb$Flank_size)/(delta_ifnb$Gene_Frequency)

delta_lps_select<-delta_lps%>%
  dplyr::filter(grepl("^CCR2$|^GPNMB$|^CD80$|^LIMK2$",Gene))

delta_ifnb_select<-delta_ifnb%>%
  dplyr::filter(grepl("^CASP10$|^GPNMB$|^TRIM56$|^GPR141$",Gene))


mfe_flank_lps<-line_plot(delta_lps,
                         delta_lps$Flank_Frequency, delta_lps$MFE_Delta,
                         delta_lps_select, 
                         delta_lps_select$Flank_Frequency, delta_lps_select$MFE_Delta,
                        "Flank size/ Gene Frequency")+
  theme(legend.position = "none")

mfe_flank_ifnb<-line_plot(delta_ifnb,
                          delta_ifnb$Flank_Frequency, delta_ifnb$MFE_Delta,
                          delta_ifnb_select, 
                          delta_ifnb_select$Flank_Frequency, delta_ifnb_select$MFE_Delta,
                          "Flank size/ Gene Frequency")+
  theme(legend.position = "none")
mfe_flank<-ggarrange(mfe_flank_lps, mfe_flank_ifnb, ncol=2)

mfe_flank


################################################################################
# plots with DEG information 
################################################################################

mfe_deg<-function(deg_in, delta_mfe_in){
  deg<-read_tsv(deg_in)
  delta_mfe<-read_tsv(delta_mfe_in)

  deg<-deg%>%dplyr::filter(GENENAME %in% delta_mfe$Gene)
  deg$Gene<-deg$GENENAME
  deg$DEG_Direction<-deg$DE_Direction
  deg<-deg%>%dplyr::select(c(Gene,logFC,DEG_Direction))
  
  delta_mfe_deg<-merge(delta_mfe, deg, by="Gene")
  print(length(unique(delta_mfe_deg$Gene)))
  delta_mfe_deg
}

delta_mfe_lps_deg<-mfe_deg("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/TopTable/Gene_Differential/LPS_DEG_Toptable.tsv",
                           "MFE/Delta_MFE/LPS_Delta_MFE.tsv")
delta_mfe_lps_deg_select<-delta_mfe_lps_deg%>%
  dplyr::filter(grepl("^CCR2$|^GPNMB$|^CD80$|^LIMK2$",Gene))


delta_mfe_ifnb_deg<-mfe_deg("/Users/hyominseo/Desktop/RAJ_RNA_Editing/Discovery/TopTable/Gene_Differential/IFNb_DEG_Toptable.tsv",
                           "MFE/Delta_MFE/IFNb_Delta_MFE.tsv")
delta_mfe_ifnb_deg_select<-delta_mfe_ifnb_deg%>%
  dplyr::filter(grepl("^CASP10$|^GPNMB$|^TRIM56$|^GRP141$",Gene))

# LPS : CCR2, GPNMB, CD80 
# IFNb : GPNMB, CASP10, TRIM56

# -log(Delta MFE )
################################################################################
mfe_logfc_lps <- line_plot(delta_mfe_lps_deg, delta_mfe_lps_deg$logFC, 
                           delta_mfe_lps_deg$MFE_Delta, 
                           delta_mfe_lps_deg_select, delta_mfe_lps_deg_select$logFC, 
                           delta_mfe_lps_deg_select$MFE_Delta,"DEG Log Fold Change") + 
  theme(legend.position = "none")


mfe_logfc_ifnb<-line_plot(delta_mfe_ifnb_deg,delta_mfe_ifnb_deg$logFC, 
                         delta_mfe_ifnb_deg$MFE_Delta,
                         delta_mfe_ifnb_deg_select, delta_mfe_ifnb_deg_select$logFC, 
                         delta_mfe_ifnb_deg_select$MFE_Delta,"DEG Log Fold Change") + 
  theme(legend.position = "none")



mfe_logfc<-ggarrange(mfe_logfc_lps, mfe_logfc_ifnb, ncol=2)

mfe_logfc

################################################################################
#gpn_plot<-ggarrange(delta_ifnb_gpn_mfe_frq, check_plot, ncol=2)
new_all_mfe_plot<-ggarrange(mfe_all, mfe_gene, mfe_flank, mfe_logfc,
                            #delta_ifnb_gpn_mfe_frq, check_plot, 
                        ncol=2,nrow=2)
ggsave(plot=new_all_mfe_plot,filename="Figure/MFE/MFE_figures.jpg",width = 12, height = 12,
       dpi = 600)
################################################################################
################################################################################
################################################################################
# A. gene selection for TIEBRUSH and seperate-editing event bed files. 

make_gene_bed<-function(sep_ed_in, gene_name,  filename){
  # have to change stranding in manually 
  
  ed_in<-read_tsv(sep_ed_in)
  ed_gene<-ed_in%>%dplyr::filter(grepl(gene_name,Gene))
  ed_gene$base_1<-c('.')
  ed_gene$base_2<-c('.')
  ed_gene$stranding<-c('-')
  ed_gene$Start <- (ed_gene$Location -1)
  
  ed_gene<-ed_gene%>%select(Chr, Start, Location, base_1, base_2, stranding)
  write_tsv(ed_gene, filename, col_names=FALSE)
}

## GBP4 hg38 chr1:89,181,144-89,198,942
## Testing GBP4
make_gene_bed('/Users/hyominseo/Desktop/RAJ_RNA_Editing/MFE/Edited_fasta/Edited_All_IFNb_Ed_Sites.tsv',
              "GBP4",# - stranding
              "MFE/Bed_files/Gene_Specific/GBP4_IFNb_Ed_Sites.bed")

# Call the function for each gene
genes <- c( "LIMK2", "PECAM1", "GPR141", "HLA-DOA", "CASP10")
for (gene in genes) {
  make_gene_bed(
    "/Users/hyominseo/Desktop/RAJ_RNA_Editing/MFE/Edited_fasta/Edited_All_IFNb_Ed_Sites.tsv",
    gene,
    paste0("MFE/Bed_files/Gene_Specific/", gene, "_IFNb_Ed_Sites.bed")
  )
}

##LIMK2 chr22:31,248,487-31,277,638
##PECAM1 chr17:64,319,415-64,390,860
##GPR141 chr7:37683766-37743835
##HLA-DOA chr6:33004182-3300991
##CASP10  chr2:201182872-201200728



# 2. make subset of select gene in all bam files in minerva 
#  for file in *.bam; do
#     filename=$(basename -- "$file")
#     filename_no_ext="${filename%.*}"
#     output_file="GPNMB_bams/${filename_no_ext}_GPNMB.bam"
#     samtools view -bh "$file" chr7:23215922-23318402 > "$output_file"
#     samtools index "$output_file"
# done

# 3. load that subset bam file in subset folder symulink to tiebrush directory by stimulation 
# find "$PWD" -name "*-Basal_LIMK2.bam" -exec sh -c 'echo "$1"; ln -s "$1" 
# /sc/arion/projects/ad-omics/flora/RNA-pipelines/tiebrush/LIMK2/all_control/' _ {} \;

# 4. ran tiebrush by editing the config file 
# conda activate snakemake
# snakemake -s Snakefile --configfile config.yaml -c2 

# 5. load the seperate editing event bed file and control.lps.ifm merged bam on IGV 



# B. select gene for MFE combindation analysis on python and RNA fold 
# 1. go to python script


group_color <- c("1 event" = "#6c97be", "2 events" = "#5486b4",
                 "3 events" = "#3b75a9","4 events" = "#23649e",
                 "5 events" = "#0b5394", "All events" ="#063158" )

subgroup_plot <- function(df_in, group_in, y_in, plot_title) {
  qc_plot <- ggplot(data = df_in, aes(x = group_in, y = y_in, color = group_in)) +
    geom_point(size = point_size + 0.5, alpha = 0.7) +
    geom_smooth(method = "lm", level = 0.98, se = FALSE, show.legend = FALSE) +
    geom_boxplot(fill=NA)+
    stat_cor(method = "pearson") +
    labs(x = "", y = plot_title) +
    scale_color_manual(values = group_color) +  # Assuming cohort_color is defined somewhere
    theme_classic() + theme(
      legend.title = element_text(size = text_size, color = "black"),
      legend.text = element_text(size = text_size, color = "black"),
      axis.text.y = element_text(size = axis_size, color = "black"),
      axis.text.x = element_text(size = axis_size, color = "black"),
      axis.title = element_text(size = axis_size, color = "black"),
      plot.title = element_text(size = title_size, hjust = 0.5, color = "black"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-3, -3, -3, -3)) +
    guides(color = guide_legend(override.aes = list(size = anno_legend_size))) +
    theme(legend.position = "none")
  
  return(qc_plot)
}

# LIMK2
delta_limk_mfe<-read.table("MFE/Delta_MFE/Gene_specific/LIMK2_Edited_All_IFNb_Ed_Sites_Delta_MFE.tsv",
                           sep='\t',header=TRUE)

delta_limk_mfe$Cohort<-c("LIMK2_IFNb")
     
delta_limk_group<-delta_limk_mfe%>%dplyr::select(c(Edited_seq_num,Delta_MFE, Edited_Frequency))
delta_limk_group$Subgroup <- case_when(
  delta_limk_group$Edited_Frequency == 1 ~ "1 event",
  delta_limk_group$Edited_Frequency == 2 ~  "2 events",
  delta_limk_group$Edited_Frequency == 3 ~  "3 events",
  delta_limk_group$Edited_Frequency == 4 ~  "4 events",
  delta_limk_group$Edited_Frequency == 5 ~  "All events"
)

delta_limk_group$Subgroup <- factor(delta_limk_group$Subgroup, 
                                    levels = c("1 event", "2 events", "3 events", 
                                               "4 events","All events"))
delta_limk_plot<- subgroup_plot(delta_limk_group,delta_limk_group$Subgroup,
                               delta_limk_group$Delta_MFE," ")

# GPR141
delta_gpr_mfe<-read.table("MFE/Delta_MFE/Gene_specific/GPR141_Edited_All_IFNb_Ed_Sites_Delta_MFE.tsv",
                           sep='\t',header=TRUE)

delta_gpr_mfe$Cohort<-c("GPR141_IFNb")

delta_gpr_group<-delta_gpr_mfe%>%dplyr::select(c(Edited_seq_num,Delta_MFE, Edited_Frequency))
delta_gpr_group$Subgroup <- case_when(
  delta_gpr_group$Edited_Frequency == 1 ~ "1 event",
  delta_gpr_group$Edited_Frequency == 2 ~  "2 events",
  delta_gpr_group$Edited_Frequency == 3 ~  "3 events",
  delta_gpr_group$Edited_Frequency == 4 ~  "4 events",
  delta_gpr_group$Edited_Frequency == 5 ~  "5 events",
  delta_gpr_group$Edited_Frequency == 6 ~  "All events"
)

delta_gpr_group$Subgroup <- factor(delta_gpr_group$Subgroup, 
                                    levels = c("1 event", "2 events", "3 events", 
                                               "4 events","5 events", "All events"))

delta_gpr_plot<- subgroup_plot(delta_gpr_group,delta_gpr_group$Subgroup,
                               delta_gpr_group$Delta_MFE," ")

delta_gene_plot<-ggarrange(delta_limk_plot, delta_gpr_plot, nrow=2)#

new_all_mfe_plot<-ggarrange(mfe_all, mfe_gene, mfe_flank, mfe_logfc,
                            delta_limk_plot, delta_gpr_plot, ncol=2,nrow=3,
                            labels=c("A","B","C","D","E","F"),
                            font.label = list(color = "black", size = label_text_size,
                                              margin=c(label_margin,label_margin)))

new_all_mfe_plot<-ggarrange(mfe_all, mfe_all,mfe_gene,
                            mfe_gene, mfe_gene,mfe_gene,
                            mfe_gene, delta_gene_plot, ncol=2,nrow=3,
                            labels=c("A","B","C","D","E","F"),
                            font.label = list(color = "black", size = label_text_size,
                                              margin=c(label_margin,label_margin)))

mfe_plot_1<-ggarrange(mfe_all, mfe_all,mfe_gene,
                       mfe_gene, mfe_gene,mfe_gene, ncol=3,nrow=2,
                            labels=c("A","B","C","D","D","D"),
                            font.label = list(color = "black", size = label_text_size,
                                              margin=c(label_margin,label_margin)))

mfe_plot_2<-ggarrange(delta_limk_plot,delta_limk_plot,
                      delta_gpr_plot,delta_gpr_plot, ncol=2,nrow=2,
                      labels=c("E","F","G","H"),
                      font.label = list(color = "black", size = label_text_size,
                                        margin=c(label_margin,label_margin)))

mfe_plot_all<-ggarrange(mfe_plot_1, mfe_plot_2, nrow=2)

ggsave(plot=mfe_plot_all,filename="Figure/MFE/Final_MFE_figures.jpg",width = 12, height = 12,
       dpi = 600)


#limk_line_plot<-line_plot(delta_limk_mfe,delta_limk_mfe$Edited_Frequency, delta_limk_mfe$Delta_MFE,"LIMk2")

################################################################################
################################################################################
################################################################################
################################################################################