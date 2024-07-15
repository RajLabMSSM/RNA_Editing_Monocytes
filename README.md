* ADD MPRA ANALYSIS authors 
### Cytosine-to-Uracil RNA editing is upregulated by pro-inflammatory stimulation of myeloid cells

* Hyomin Seo1-5, Winston Cuddleston1-5, Madison Parks1-5 Elisa Navarro1-5, Towfique Raj1-5, Jack Humphrey1-5

1. Department of Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai, New York, NY, USA
2. Nash Family Department of Neuroscience & Friedman Brain Institute, Icahn School of Medicine at Mount Sinai, New York, NY, USA
3. Ronald M. Loeb Center for Alzheimer’s Disease, Icahn School of Medicine at Mount Sinai, New York, NY, USA
4. Icahn Institute for Data Science and Genomic Technology, Icahn School of Medicine at Mount Sinai, New York, NY, USA
5. Estelle and Daniel Maggin Department of Neurology, Icahn School of Medicine at Mount Sinai, New York, NY, USA

-----
## Raj Lab, Hyomin Seo. 2024 July

1) [Sample Cohorts](#sample-cohorts)
   - Discovery, Replication, and iMicroglia (Alternative) cohorts.
2) [Rmarkdowns_PCA/DES/DEG](#Rmarkdowns-(PCA/DES/DEG)) 
   - Editing_Paper_Figure.Rmd
   - Editing_Paper_Functions.Rmd
   - Discovery/Discovery_Analysis & _Data_Prep.Rmd
   - Replication/Replication_Analysis & _Data_Prep.Rmd
   - Alternative/Alternative_Analysis & _Data_Prep.Rmd
3) [Jacusa outputs](#Jacusa-outputs)

4) [MFE: Minimum Free Energy](#MFE-analysis)
   - Fasta_RNAFold.R & MPRA_Fasta_RNAFold.R
   - Python/MFE_Multi_Editing.ipynb, _Editing_MPRA.ipynb, GENENAME.ipynb
   - Bed_files/, Raw_fasta/, Edited_fasta/, RNAfold/, Delta_MFE/

5) [MPRA: Massively Parallel Reporter Assay](#MPRA-analysis)
   - MPRA.R
   - MPRA_MFE_Discovery_IFNb & LPS_Genes.tsv
   - MPRA/Input, MPRA/Result
   - MPRA/Pathway
  
6) [Tiebrush](#tiebrush)
   - Gene-specific files for IGV viewing
  
6) [Figures](#figures)
   - Instead of having knitted markdown, we created sets of organized figures
   - Main_Figure/ , Supp_Figure/, MFE/
_____
## Sample Cohorts 
### Discovery Cohort
* 55 Individual donors

| Stimulation  | Basal   | IFNbeta |  LPS    |   SUM    |
|     :--:     |  :---:  |  :---:  |  :---:  |   :---:  | 
| Sample       | 55      | 55      |  55     |   165    | 


### Replication Cohort
* 35 Individual donors

| Stimulation  | Basal   | IFNgamma |  LPS    |  SUM    |
|     :--:     |  :---:  |  :---:   |  :---:  |  :---:  | 
| Sample       | 35      | 35       |  35     |  105    | 


### iMicroglia Cohort (Alternative) 
* 3 Individual donors ( I am not even sure) 

| Stimulation  | Basal   | IFNgamma |  LPS    |  SUM    |
|     :--:     |  :---:  |  :---:   |  :---:  |  :---:  | 
| Sample       | 35      | 35       |  35     |  105    | 
_____
## Rmarkdowns (PCA/DES/DEG)
* All the makrdowns are in the form of Rmd but **not configured to be an organized-knitted html**, because the focus was to produce publication-worth figures.
1) Editing_Paper_Figure.Rmd & Editing_Paper_Functions.Rmd
   - Markdown where all main figures were produced/ all functions used to produce the figures in seperate rmd. 
2) [Discovery/Discovery_Analysis.Rmd & _Data_Prep.Rmd](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/Discovery)
   - Markdown for **PCA analysis, DES and DEG** for discovery cohort/ sets of scripts to organize the raw data for discovery cohorts
3) [Replication/Replication_Analysis.Rmd &_Data_Prep.Rmd](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/Replication)
   - Same as 2), for Replication cohort 
3) [Alternative/Alternative_Analysis.Rmd & _Data_Prep.Rmd](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/Alternation)
   - Same as 2), for Alternative (iMicroglia) cohort 
___
## Jacusa outputs

___
## [MFE Analysis](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/MFE)
* This can be organized into one snakemake pipeline - for future 
1) Fasta_RNAFold.R & MPRA_Fasta_RNAFold.R
   - MFE analysis Rscript for IFNb and LPS genes (2000bp)
   - MFE analysis Rscript for MPRA specific genes (300bp)
     
2) Python/MFE_Multi_Editing.ipynb, _Editing_MPRA.ipynb, GENENAME.ipynb
   - Python script for making combination MFE analysis for IFNb and LPS genes
   - for MPRA specific genes
   - for gene-specific combination
     
3) Bed_files/, Raw_fasta/, Edited_fasta/, RNAfold/, Delta_MFE/
   - folders with all sub-data created from R and Python scripts

___
## [MPRA Analysis](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/MPRA)
1) MPRA.R
   - Script for MPRA analysis for overall significant-volcano plot, miRNA-significant-volcano plot, pathway analysis for miRNA_significant genes (nothing significant), getting list of  IFNb/LPS significant MPRA genes for MFE analysis. 

2) MPRA_MFE_Discovery_IFNb & LPS_Genes.tsv
   - list of IFNb/LPS significant MPRA genes for MFE analysis. 
3) MPRA/Input
   - The list of all 3'UTR sites we sent out for MPRA validation.
5) MPRA/Result
   - Result of MPRA validation
7) MPRA/Pathway
   - Result for pathway analysis on miRNA-significant genes, nothing significant was found 

___
## [Tiebrush](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/Tiebrush)
1) Gene specific folders
   - For specific genes, bam, bai ,bigwig, bed files for Basal samples, IFN-stimulated samples, and LPS-stimulated samples, to be loaded in IGV to visually see Editing rates. 
  
___  
## [Figures](https://github.com/RajLabMSSM/RNA_Editing_Monocytes/tree/main/Figure)
1) Main_Figure
   - Figure 1~5 (.png), created in Editing_Paper_Figure.Rmd + Adobe illustrator
     
2) Supp_Figure
   - Supp figure (PCA figures) for discovery, replication, alternation (iMicroglia), created seperatley in cohort-wise analysis.Rmd
   - MPRA supplement figure, created in Editing_Paper_Figure.Rmd + Adobe illustrator

3) MFE
   - Two MFE figures (new_Final_MFE_Figure, new_Supp_MFE_Figure.png) 
___  
## Abastract 

RNA editing is a posttranscriptional modification that produces RNA molecules with different sequence information from what is genomically encoded. The consequences of RNA editing include not only the production of alternative protein products but also changes in gene expression, alterations in RNA stability, and modifications to splicing patterns. The two most frequent types of RNA editing are Cytosine to Uracil (C-to-U, read as Thymine), facilitated by APOBEC enzyme genes, and Adenosine to Inosine (A-to-G, read as Guanosine), facilitated by the ADAR enzyme family of genes. It has been demonstrated that ADAR-mediated RNA editing critically regulates the innate immune interferon (IFN) pathways, and genetic variants associated with several immunoinflammatory diseases dysregulate this vital function of RNA editing. 

Neurodegenerative diseases such as Alzheimer’s Disease (AD) and Parkinson’s Disease (PD) are triggered by DNA variants that affect the genes expressed in myeloid cells, a type of white blood cell that plays a critical role in the immune system. Genome-wide association study (GWAS) has identified a number of genetic variants linked to AD and PD disease risk. Further studies showed that these identified variants were explicitly expressed in myeloid cells of the innate immune system, including monocytes in the periphery and microglia in the central nervous system. 

We hypothesized that RNA editing plays a critical role in myeloid cells and that GWAS-risk variants for AD and PD regulate RNA editing activity in these cells. To gain a deeper understanding of the relationship between RNA editing in myeloid cells and neurodegenerative disease, we conducted a study investigating how RNA editing in monocyte samples contributes to the genetic risk of these diseases. We performed differential analysis on editing and gene expression levels in two monocyte sample cohorts and one induced microglia sample cohort that were subjected to external stimuli or diagnosed with neurodegenerative diseases. 
