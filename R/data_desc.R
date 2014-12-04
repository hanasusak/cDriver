#' Data that are included in cDriver package
#'
#' @name lawrence.df
#' @title Lawrence's paper information
#' @description This data frame containts information used in package which is provided in Lawrence paper.
#' There is 18859 rows and each row represents gene. Rownames are Hugo symbols for genes. They were corrected to most recent used Hugo gene name.
#' Data frame has 4 columns.
#' @section Variables:
#' \itemize{
#'  \item N_nonsilent number of all possible nonsilent substitutions
#'  \item N_silent number of all possible silent substitutions
#'  \item N_noncoding number of all possible noncoding substitutions
#'  \item noncoding_mutation_rate measured noncoding mutation rate - somatic background mutation rate
#'  }
#' @docType data
#' @usage lawrence.df
#' @format data grame
#' @source Suplementary matirial from Lawrence paper: \url{http://www.ncbi.nlm.nih.gov/pubmed/23770567}.
# @author Hana Susak, 1/4/2014
#' @keywords data
NULL


#' @name length.genes
#' @title Genes length by NCBI, length which was covered by sequencing,  
#' @description This data frame containts information about genes length. 
#' There is 16485 rows and each row represents gene. This set of genes were analysed as part of CLL somatic driver genes at CRG and ICGC. 
#' Rownames are Hugo symbols for genes. They were corrected to most recent used Hugo gene name.
#' Data frame has 5 columns.
#' @section Variables:
#' \itemize{
#'  \item Hugo_Symbol same as rowname. Hugo symbol as unique indentifier for genes.
#'  \item Coverd_len length of the gene which was coverd by sequncing CLL project.
#'  \item Length length of genes reported in NCBI
#'  \item procentege which was covered of total gene (Coverd_len / Length)
#'  \item p_by_len  probability that gene get one mutation, if each base pair had same chance; Coverd_len/sum(length.genes$Coverd_len)
#'  }
#' @docType data
#' @usage length.genes
#' @format data grame
# @source Suplementary matirial from Lawrence paper: \url{http://www.ncbi.nlm.nih.gov/pubmed/23770567}.
# @author Hana Susak, 1/4/2014
#' @keywords data
NULL


#' @name all.genes.lengths
#' @title Coding genes length by Ensamble, length is calculated as sum of all coding exons (merged before so each position only once ocunted). 
#' @description This data frame containts information about genes length. 
#' There is 19202 rows and each row represents gene. This set of genes is acquired from HGNC website and only approved and protein coding genes are included.. 
#' Rownames are Hugo symbols for genes. 
#' Data frame has 2 columns.
#' @section Variables:
#' \itemize{
#'  \item Hugo_Symbol same as rowname. Hugo symbol as unique indentifier for genes.
#'  \item Length length of genes reported in NCBI
#'  }
#' @docType data
#' @usage all.genes.lengths
#' @format data grame
#' @source http://www.genenames.org/ (Most recent HUGO approved names ) and ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens (info for genes, exons, lenghts)
# @author Hana Susak, 1/4/2014
#' @keywords data
NULL

#' @name sample.genes.mutect
#' @title Genes length by NCBI, length which was covered by sequencing,  
#' @description This data frame containts information about SNVs and InDels that are reported as part of CLL project. 
#' There is 7860 rows and each row represents annotated SNV/InDel reported in one of 337 patients with CLL. 
#' These are filtered SNVs and InDels for segmental duplication, exonix or splicing, intersection of used kits region, etc.
#' They follow MAF like format with additional columns like ploidy, CCF_CNV, etc.
#' Data frame has 23 columns.
#' @section Variables:
#' \itemize{
#'  \item Chromosome 
#'  \item Start_Position 
#'  \item Reference_Allele
#'  \item Tumor_Seq_Allele2 
#'  \item Hugo_Symbol
#'  \item Tumor_Sample_Barcode
#'  \item VAF
#'  \item Function.Refseq
#'  \item ExonicFunction.Refseq 
#'  \item ...
#'  \item ploidy
#'  \item CCF_CNV
#'  \item Damege_score
#'  \item gender
#'  \item Variant_Classification
#'  \item Variant_Type
#'  }
#' @docType data
#' @usage sample.genes.mutect
#' @format data grame
# @source Suplementary matirial from Lawrence paper: \url{http://www.ncbi.nlm.nih.gov/pubmed/23770567}.
# @author Hana Susak, 1/4/2014
#' @keywords data
NULL