######################################################################################
# CRG 
# Hana SUSAK
######################################################################################

# Function to calculate ploidy
# @param VAF - Variant allele frequency observed in reads; 
# @param ploidy - ploidy in position of reported variant (optional, default = 2 ). In other words, is this variant together with CNV;
# @param ccf_cnv - Cancer Cell Fraction of this ploidy. For germline CNVs its 1, and for somatic CNVs it can take values in interval (0,1] (optional, default = 1);
# @param purity - purity of cancer tissue and it is value in interval (0,1] but is expected to be high, much closer to 1 then 0.  (optional, default = 1)
ccfPloidy <- function (vaf, ploidy = 2, ccf_cnv = 1, purity = 1) {
    if (sum(is.na(ploidy))){
        ploidy[is.na(ploidy)] <- 2
    }
    if (sum(is.na(ccf_cnv))){
        ccf_cnv[is.na(ccf_cnv)] <- 1
    }  
    if (sum(is.na(purity))){
        purity[is.na(purity)] <- 1
    }  
    ccf <- ((2 + (ploidy-2)*ccf_cnv)*vaf)/purity    
    return(ccf)
}


# function to correct CCF above 1
# asumptions considered to correct CCF:
#   1) 1 < ccf <= 1.2 and ploidy = 2 =>  rough estimation of baf which should be 0.5 therefore CCF should be 1
#   2) 1.2 < ccf and ploidy = 2   =>  missing deletion
#   3) ploidy != 2 and ccf > 1 => CNV and SNV are found in fraction of same cells, so estimation is overestimated as above 1, and should be 1.
# In case there is no ploidy column, it will be assumed as 2
# @param sample.mutations - Data Frame with columns: 'VAF', 'ploidy', and 'CCF'
ccfCorrection <- function(sample.mutations){  
    if  (!'purity' %in% colnames(sample.mutations)){
        
        # correct BAF between 0.5 and 0.6 and diploid
        if ( 'ploidy' %in% colnames(sample.mutations) ){
            condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & sample.mutations$ploidy == 2   
        } else {
            condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition, ]$CCF <- 1        
        }
        
        # correct BAF between 0.6 and 1 and diploid
        if ( 'ploidy' %in% colnames(sample.mutations) ){    
            condition <- sample.mutations$vaf > 0.6  &  (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))
        } else {
            condition <- sample.mutations$vaf > 0.6 
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1)
        }
        
        # correct ploidy != 2 and ccf >1
        if ( 'ploidy' %in% colnames(sample.mutations) ){   
            condition <- sample.mutations$CCF > 1  & (sample.mutations$ploidy != 2   | is.na(sample.mutations$ploidy ))
        } else {
            condition <- sample.mutations$CCF > 1      
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition, ]$CCF <- 1
        }
        
    } else {
        if (sum(is.na(sample.mutations$purity))){
            sample.mutations[is.na(sample.mutations$purity),'purity'] <- 1
        } 
        
        # correct BAF between 0.5 and 0.6 and diploid
        if ( 'ploidy' %in% colnames(sample.mutations) ){
            condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))  
        } else {
            condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition, ]$CCF <-  min( (sample.mutations[condition, ]$vaf*2  / sample.mutations[condition, ]$purity ) , 1)    
        }
        
        # correct BAF between 0.6 and 1 and diploid
        if ( 'ploidy' %in% colnames(sample.mutations) ){    
            condition <- sample.mutations$CCF > 1.2  &  (sample.mutations$ploidy == 2  | is.na(sample.mutations$ploidy ))
        } else {
            condition <- sample.mutations$CCF > 1.2 
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1, purity=sample.mutations[condition ,]$purity)
        }
        
        # correct ploidy != 2 and ccf >1
        if ( 'ploidy' %in% colnames(sample.mutations) ){   
            condition <- sample.mutations$CCF > 1  #& sample.mutations$ploidy != 2   
        } else {
            condition <- sample.mutations$CCF > 1      
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition, ]$CCF <- 1
        }
        
    }
    
    
    sample.mutations
}




# function to correct purrity 
# asumptions considered to correct purity:
#   1) 95% of snps are in interval 0-1.2
#   2) 3 or more snps are above 1.2
# In case SNVs after correction for purity (and CNV if provided) are violating 2 mentioned condicions, purity is not used.
purityCorrection <- function(sample.mutations){  
    if  (!'purity' %in% colnames(sample.mutations)){
       stop('There need to be purity column for correction by purity')      
    } else {
        ## check conditions for each patient, less then 5% and less then 3 SNVs above 1.2 CCF estmated
        
        
        if (sum(is.na(sample.mutations$purity))){
            sample.mutations[is.na(sample.mutations$purity),'purity'] <- 1
        } 
        
        # correct BAF between 0.5 and 0.6 and diploid
        if ( 'ploidy' %in% colnames(sample.mutations) ){
            condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))  
        } else {
            condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition, ]$CCF <-  min( (sample.mutations[condition, ]$vaf*2  / sample.mutations[condition, ]$purity ) , 1)    
        }
        
        # correct BAF between 0.6 and 1 and diploid
        if ( 'ploidy' %in% colnames(sample.mutations) ){    
            condition <- sample.mutations$CCF > 1.2  &  (sample.mutations$ploidy == 2  | is.na(sample.mutations$ploidy ))
        } else {
            condition <- sample.mutations$CCF > 1.2 
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1, purity=sample.mutations[condition ,]$purity)
        }
        
        # correct ploidy != 2 and ccf >1
        if ( 'ploidy' %in% colnames(sample.mutations) ){   
            condition <- sample.mutations$CCF > 1  #& sample.mutations$ploidy != 2   
        } else {
            condition <- sample.mutations$CCF > 1      
        }
        if (sum(condition,  na.rm = T)) {
            condition[is.na(condition)] <- FALSE
            sample.mutations[condition, ]$CCF <- 1
        }
        
    }
    
    
    sample.mutations
}




#' Calculation of Cancer Cell Fraction (CCF) for SNVs from allele frequency (VAF).
#' @description
#'   \code{CCF} function calculates  CCF for each variant based on its 
#'  allele frequency, CNV/ploidy context, cancer cell fraction of reporeted CNVS within variant position and purity of tumor tissue.
#' @param sample.mutations Data Frame which should follow MAF format. Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item VAF variant allele frequncey for reported SNV
#'      \item ploidy (optional, default = 2) ploidy within reoported SNV. 
#'      For example if SNV is reporeted in Y chromosome and with no CNV in this position, ploidy should be 1.
#'      If gender is not known, than recomandation is to to exclude all SNVs with X chromosome.
#'      \item CCF_CNV (optional, default = 1) cancer cell fraction of somatic SNV in region with reported SNV. 
#'      \item purity (optional, default = 1) purity for sample in which  SNV is reported.
#' } 
#' If not provided they need to be specifed as paramiters of the CCF function.
#' @param VAF (optional) integer/numeric value indicating column in \code{sample.mutations} representing variant allele frequncey for reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param ploidy (optional) integer/numeric value indicating column in \code{sample.mutations} representing ploidy context of reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column, or default value of 2 is taken)
#' @param CCF_CNV (optional) integer/numeric value indicating column in \code{sample.mutations} representing CCF of CNV which is reportedin region of reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column, or default value of 1 is taken)
#' @param purity (optional) integer/numeric value indicating column in \code{sample.mutations} representing purity of tumor tissue for sample with reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column, or default value of 1 is taken)
#' @param correct (optional, default = TRUE) Correction to perform on SNVs for which CCF is calculated as larger then 1. 
#'      This is justifed with rough estimation of VAF values, missing CNVs and 
#'      violation of mutal exclusivit assumption (two mutatations in same gene/patient are in different cancer frations ). 
#'      It is recomanted to keep this parameter to TRUE value, othervise unrealistic CCF (> 1) values can be returned for some SNVs.
#' @return a data frame with one additional column, giving CCF vlaues for each SNV in intial \code{sample.mutations} data frame.
#' @keywords CCF
#' @examples
#' # Simulate some VAF, ploidy and CCF_CNV values
#' df <- data.frame(VAF=runif(100, min=0.05, max=0.75), 
#'                  ploidy=sample(c(1:4), 100, replace=TRUE, prob=c(0.4,0.9,0.5,0.1)), 
#'                  CCF_CNV=runif(100, min=0.1,max=1))
#' df[df$ploidy == 2, 'CCF_CNV'] <- 1
#' # call CCF function
#' df2 <- CCF(df)
#' head(df2)
#' @export
CCF <- function(sample.mutations, VAF = NULL, ploidy = NULL, CCF_CNV = NULL, purity = NULL, correct=TRUE){
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    } 
    
    if (!is.null(VAF)){
        sample.mutations <- assign.columns(sample.mutations, VAF, "VAF")
    }
    if (!is.null(ploidy)){
        sample.mutations <- assign.columns(sample.mutations, ploidy, "ploidy")
    }
    if (!is.null(CCF_CNV)){
        sample.mutations <- assign.columns(sample.mutations, CCF_CNV, "CCF_CNV")
    }
    if (!is.null(purity)){
        sample.mutations <- assign.columns(sample.mutations, purity, "purity")
    }
    
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    # check if BAF column is there
    if ( 'vaf' %in% colnames(sample.mutations) ){
        if  (!is.numeric(sample.mutations$vaf)){
            stop("VAF column is not numeric!")
        }            
    } else {
        stop("There is no mandatory VAF column!")
    }
    
    if ( 'ploidy' %in% colnames(sample.mutations) ){
        if  (!is.numeric(sample.mutations$ploidy)){
            stop("Ploidy column is not numeric!")
        }   
        if ( 'ccf_cnv' %in% colnames(sample.mutations) ){
            if  (!is.numeric(sample.mutations$ccf_cnv)){
                stop("CCF_CNV column is not numeric!")
            }
            if ('purity' %in% colnames(sample.mutations) ) {
                # calculate CCF as ploidy is 2 
                sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv, purity=sample.mutations$purity)
            } else {
                # calculate CCF! there is baf, ploidy and ccf of cnv
                sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv)
            }
        } else {
            if ('purity' %in% colnames(sample.mutations) ) {
                # calculate CCF as ploidy is 2 
                sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy,  purity=sample.mutations$purity)
            } else {
                # calculate CCF! there is baf, ploidy and ccf of cnv
                sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy)  
            }
        }           
    } else {
        if ('purity' %in% colnames(sample.mutations) ) {
            # calculate CCF as ploidy is 2 
            sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, purity=sample.mutations$purity)
        } else {
            # calculate CCF as ploidy is 2 
            sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf)
        }

    }

    if (correct){
        sample.mutations <- ccfCorrection(sample.mutations)
    }
    
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    sample.mutations
    
}


