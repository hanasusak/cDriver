######################################################################################
# CRG 
# Hana SUSAK
#------------------------------------------------------------------------------------
# Bayes models
######################################################################################
#' @import ggplot2
#' @import data.table
## @import reshape2
#' @import Rmpfr
#' @importFrom stats sd aggregate as.formula prcomp xtabs
#' @importFrom grDevices dev.off pdf
#' @importFrom utils methods write.table
## @importFrom base apply
 
# Function to create matrix genes versus patients. 
# Values of matrix can be counts of mutations, or sum of some other variable.
# Also it user can provide upper limit for sum of this value for gene-patienet pair.
# @param sample.mutations data frame to be converted gene-patient matrix
# @param genes vector of unique genes of interest
# @param valueCol Column for which cuming gonna be done (like CCF, Demage), should be numeric or character (exact name of column).
#                   If NULL, it will count occurences.
# @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs.
#       Default is NULL value (in this case \code{sample.mutations} should already have this column)
# @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs. 
#      Default is NULL value (in this case \code{sample.mutations} should already have this column)
# @param sample.gene.lim a numeric value specifying upper limit when summing is perfomed for gene-sample pair
# @param mode a charechter value indicationg how to solve when in one gene-sample pair there are multiple mutations. Options are SUM, MAX and ADVANCE
# @param epsilon a numeric value. If mode is ADVANCE, epsilone value will be threshold for CCF difference to decide if they are in same or different clone. 
pat.vs.genes <- function (sample.mutations, genes = NULL, valueCol = NULL, Hugo_Symbol=NULL, Tumor_Sample_Barcode=NULL, sample.gene.lim = NULL, mode='MAX', epsilon= 0.05 ){   
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    mode <- toupper(mode)
    if ( !mode %in% c('SUM', 'MAX', 'ADVANCE')) {
        stop("mode mast be or SUM or MAX or ADVACE! ", call. = FALSE)
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mutations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(sample.gene.lim) & (length(sample.gene.lim) > 1 | !is.numeric(sample.gene.lim))) {
        stop("sample.gene.lim must be a single numeric variable", call. = FALSE)
    }
    
    if (!is.null(valueCol) & (length(valueCol) > 1 )) {
        stop("valueCol must be a single variable, charachter or numeric", call. = FALSE)
    } 
    
    if( !is.null(valueCol) ){
        if (is.character(valueCol) & (! valueCol %in% colnames(sample.mutations))) {
            stop("valueCol must one of column names in sample.mutations", call. = FALSE)
        } 
        
        if (is.numeric(valueCol) & (valueCol > ncol(sample.mutations))) {
            stop("valueCol must be numeric value indicating valid column", call. = FALSE)
        }   
    }
    
    #if (is.null(valueCol) ){
    #    message('As valueCol is not specified only counting of mutations for sample-gene pairs will be done.')
    #}
      
    if (is.null(valueCol) & !is.null(sample.gene.lim) ){
        message('With limit value for sample-gene pair, there need to be specified targeted column (valueCol). Limit value will be ignored.')
    }
    
    #if ( mode!='ADVANCE' & !is.null(epsilon) ){
    #    message('Epsilon value is only used with mode ADVANCE. Epsilon will be ignored.')
    #}
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    
    if (is.null(genes)){
        genes <- as.character(unique(sample.mutations$hugo_symbol))
    }
    
    advance.ccf.calc <- function(x, epsilon){
        x <- sort(x, decreasing=T)
        sum.ccf <- x[1]
        if (length(x) >1){
            for (i in 2:length(x)) {
                if ((x[i-1] - x[i]) > epsilon){
                    sum.ccf <- sum.ccf + x[i]
                }
            }  
        }
        return(sum.ccf)
    }
    
    # solve NA values for valueCol, set to 1
    sample.mutations2 <- sample.mutations
    if (!is.null(valueCol)){
        ind <- is.na(sample.mutations2[,valueCol])
        if (sum(ind)) {
            warning("There are ", sum(ind), " NA values in column ", valueCol, " which are replaced by 1! \n")
            sample.mutations2[ind, valueCol] <- 1
        }
    }

    #test <-  acast(sample.mutations2, Hugo_Symbol ~ Tumor_Sample_Barcode, sum )
    
    
    ### creat data table, to be faster ### 
    sample.gene.dt <- data.table(sample.mutations2)
    keycols <- c("hugo_symbol","tumor_sample_barcode")
    setkeyv(sample.gene.dt,keycols)
    
    if(!is.factor(sample.gene.dt$hugo_symbol)){
        sample.gene.dt$hugo_symbol <- factor(sample.gene.dt$Hugo_Symbol, levels = (genes))
    }
    if(!is.factor(sample.gene.dt$tumor_sample_barcode)){
        sample.gene.dt$tumor_sample_barcode <- factor(sample.gene.dt$tumor_sample_barcode, levels= sort(as.vector(unique(sample.gene.dt$tumor_sample_barcode))) )
    }
  

    
    if (mode == 'SUM' & (!is.null(valueCol) )){
    # if we sum values of column valueCol, for same sample-gene pairs  
        if (is.null(sample.gene.lim)){
            sample.gene.lim <- Inf
        }
        if (is.character(valueCol)){
            f <- paste(valueCol, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            sample.gene.dt <- aggregate(as.formula(f) , data=sample.gene.dt, function(x) min(sum(x),sample.gene.lim))         
        } else if (is.numeric(valueCol)){
            name.col <- colnames(sample.mutations2)[as.integer(valueCol)]                 
            f <- paste(name.col, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            sample.gene.dt <- aggregate(as.formula(f) , data=sample.gene.dt, function(x) min(sum(x),sample.gene.lim))                                            
        } else {
            stop('No valid sum column (should be CCF or Damage score - like Condel)')
        }
        
    } else if (mode == 'MAX' & (!is.null(valueCol) )) {
    # if we take max value of column valueCol, for same sample-gene pairs    
        if (is.null(sample.gene.lim)){
            sample.gene.lim <- Inf
        }
        if (is.character(valueCol)){
            f <- paste(valueCol, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            sample.gene.dt <- aggregate(as.formula(f) , data=sample.gene.dt, function(x) min(max(x, na.rm=T),sample.gene.lim))         
        } else if (is.numeric(valueCol)){
            name.col <- colnames(sample.mutations2)[as.integer(valueCol)]                 
            f <- paste(name.col, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            sample.gene.dt <- aggregate(as.formula(f) , data=sample.gene.dt, function(x) min(max(x, na.rm=T),sample.gene.lim))                                            
        } else {
            stop('No valid sum column (should be CCF or Damage score - like Condel)')
        }
    } else if  (mode == 'ADVANCE' & (!is.null(valueCol) ) ){
    # if we take max or sum value of column valueCol, for same sample-gene pairs, depending if in same or different clone   
        if (is.null(sample.gene.lim)){
            sample.gene.lim <- Inf
        }
        #sample.gene.dt <- sample.gene.dt[order(sample.gene.dt$Hugo_Symbol, sample.gene.dt$Tumor_Sample_Barcode, -sample.gene.dt$CCF),]
        
        if (is.character(valueCol)){
            f <- paste(valueCol, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            sample.gene.dt <- aggregate(as.formula(f) , data=sample.gene.dt,  function(x) min(advance.ccf.calc(x, epsilon),sample.gene.lim))         
        } else if (is.numeric(valueCol)){
            name.col <- colnames(sample.mutations2)[as.integer(valueCol)]                 
            f <- paste(name.col, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            sample.gene.dt <- aggregate(as.formula(f) , data=sample.gene.dt, function(x) min(advance.ccf.calc(x, epsilon),sample.gene.lim))                                            
        } else {
            stop('No valid sum column (should be CCF or Damage score - like Condel)')
        }
    } 

    
    # creating matrix
    if (is.null(valueCol)){
        ## just counts
        mat<- xtabs(~ hugo_symbol + tumor_sample_barcode, data=sample.gene.dt)
    } else {        
        if (is.character(valueCol)){
            f <- paste(valueCol, " ~ hugo_symbol + tumor_sample_barcode",  collapse="")
            mat <- xtabs(as.formula(f), data=sample.gene.dt)     
        } else if (is.numeric(valueCol)){
            name.col <- colnames(sample.mutations2)[as.integer(valueCol)] 
            f <- paste(name.col, " ~  hugo_symbol + tumor_sample_barcode",  collapse="")
            mat <- xtabs(as.formula(f), data=sample.gene.dt)        
        } else {
            stop('No valid sum column (should be CCF or Damage score - like Condel)')
        }
    }
    
    
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    class(mat) <- 'matrix'
    return(mat)
    
}


#' Bayesina risk inference model.
#' @description
#'   \code{bayes.risk} function performs by runing Bayesian risk inference model when priors are set by user and liklihood is calculated from given data of SNVs/InDels.
#' @param sample.mutations data frame with SNVs and InDels in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item CCF numeric column produce by \code{CCF} function.
#'      \item Damage_score numeric column with values between 0 and 1, where 1 means very damaging SNV/IndDel and 0 not damaging SNV/InDel
#' } 
#' @param bcgr.prob a numeric vector, same lenght as genes (should be same orderd also) which gives probability of gene having somatic mutation in healfy population.
#'          There are functions for obtaining this vector: \code{bcgr}, \code{bcgr.lawrence} and \code{bcgr.combine}.
#' @param genes vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV/Indel. in given cohort). Default NULL.
#' @param prior.sick a numeric value representing incidence of tumor in population. Set by default to 0.0045 .
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs/Indels.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Damage_score (optional) integer/numeric value indicating column in \code{sample.mutations} which contain damage score for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param mode a charechter value indicationg how to solve when in one gene-sample pair there are multiple mutations. Options are SUM, MAX and ADVANCE
#' @param epsilon a numeric value. If mode is ADVANCE, epsilone value will be threshold for CCF difference to decide if they are in same or different clone. 
#' @return a data frame with ranked genes by posteriory probability of gene beeing risk factor for developing tumor. 
#' Additional columns with usefull info are contained in data frame.
#' @seealso \code{\link{bcgr}}, \code{\link{bcgr.lawrence}}  and \code{\link{bcgr.combine}} for obtaining bcgr.prob variable.
#' @keywords risk
#' @examples 
#' \donttest{
#' # first calculate CCF
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # then somatic background probability
#' bcgr.prob <- bcgr.combine(sample.genes.mutect)
#' # bayes risk model
#' risk.genes <- bayes.risk(sample.genes.mutect,  bcgr.prob, prior.sick = 0.00007) 
#' head(risk.genes)  
#' }
#' @export
bayes.risk <- function(sample.mutations,  bcgr.prob, genes=NULL, prior.sick = 0.0045, 
    Variant_Classification=NULL, Hugo_Symbol=NULL, Tumor_Sample_Barcode=NULL, CCF=NULL, Damage_score=NULL , mode='MAX', epsilon=0.05 ) {
    
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    mode <- toupper(mode)
    if ( !mode %in% c('SUM', 'MAX', 'ADVANCE')) {
        stop("mode mast be or SUM or MAX or ADVACE! ", call. = FALSE)
    }
    
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mutations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(CCF)){
        sample.mutations <- assign.columns(sample.mutations, CCF, "CCF")
    }
    
    if (!is.null(Damage_score)){
        sample.mutations <- assign.columns(sample.mutations, Damage_score, "Damage_score")
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    
    if (is.null(genes)){
        genes <- as.character(unique(sample.mutations$hugo_symbol))
    }
    
    
    if(!is.factor(sample.mutations$hugo_symbol)){
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=genes)
    }
    if(!is.factor(sample.mutations$tumor_sample_barcode)){
        sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=unique(sample.mutations$tumor_sample_barcode))
    }
    prior.healthy <- 1 - prior.sick
    
    # take only exonic
    sample.mutations <-  exonic.only(sample.mutations) 
    
    # matrices genes vs. samples
    mat.sample.gene.ns <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification != 'Silent',], genes)
    mat.sample.gene.s <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification == 'Silent',], genes )
    if ('damage_score' %in% colnames(sample.mutations)){
        mat.sample.gene.ns.Damage <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification != 'Silent',], genes, valueCol='damage_score', sample.gene.lim=1, mode=mode, epsilon)
    } else {
        mat.sample.gene.ns.Damage <- matrix(1, nrow=nrow(mat.sample.gene.ns), ncol=ncol(mat.sample.gene.ns))
    }
    mat.sample.gene.ns.ccf <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification != 'Silent',], genes, valueCol='ccf', sample.gene.lim=1, mode=mode, epsilon)
    mat.sample.gene.s.ccf <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification == 'Silent',], genes, valueCol='ccf', sample.gene.lim=1, mode=mode, epsilon)
    
    #condtitionals
    mat.ns.somatic <- (mat.sample.gene.ns.Damage*mat.sample.gene.ns.ccf)
    
    gene.mutated.if.sick <- rowSums(mat.ns.somatic)/length(unique(sample.mutations$tumor_sample_barcode))
    
    gene.mutated.if.healthy <-  bcgr.prob
    
    bayes.risk <- (prior.sick * gene.mutated.if.sick[genes])/(prior.sick*gene.mutated.if.sick[genes] + prior.healthy*gene.mutated.if.healthy[genes])
    
    #additional columns
    observed.mut.ns <- base::apply(as.matrix(mat.sample.gene.ns),c(1),function(x) sum(x >0))
    total.mut.ns <- rowSums(mat.sample.gene.ns)    
    observed.mut.s <- base::apply(as.matrix(mat.sample.gene.s),c(1),function(x) sum(x >0))
    observed.mut.ns.ccf <- rowSums(mat.sample.gene.ns.ccf)
    observed.mut.s.ccf <- rowSums(mat.sample.gene.s.ccf)
    n <- rowSums(mat.ns.somatic)
    n.indels <- table(sample.mutations[sample.mutations$variant_type  %in% c('DEL', 'INS'),'hugo_symbol'], exclude=F)
    
    top <- as.character(names(sort(bayes.risk, decreasing=T)))
    sample.ns.mut <- paste('sample_ns_mut(',length(levels(sample.mutations$tumor_sample_barcode)),')',sep="")
    sample.s.mut <- paste('sample_s_mut(',length(levels(sample.mutations$tumor_sample_barcode)),')',sep="")
    
    final.df.bayes1 <- data.frame('Hugo_Symbol' = top,
                                  'postProb' = bayes.risk[top],
                                  'CCF_x_damage' = n[top],
                                  'fold' = bayes.risk[top]/prior.sick,
                                  'sample_ns_mut' = observed.mut.ns[top],
                                  'total_ns_mut'= total.mut.ns[top],
                                  'nonsilent_CCF_sum' = observed.mut.ns.ccf[top],
                                  'sample_s_mut'= observed.mut.s[top],
                                  'silent_CCF_sum' = observed.mut.s.ccf[top],                                  
                                  'indels'= as.numeric(n.indels[top]),
                                  'background_mut_prob'=bcgr.prob[top],
                                  'ccf_mean'= base::apply(as.matrix(mat.sample.gene.ns.ccf[top,]),1,function(x) mean(x[x!=0])),
                                  'ccf_median'= base::apply(as.matrix(mat.sample.gene.ns.ccf[top,]),1,function(x) median(x[x!=0])),
                                  'ccf_sd'= base::apply(as.matrix(mat.sample.gene.ns.ccf[top,]),1,function(x) sd(x[x!=0])),
                                  'rank'=1:length(top) 
                                  )
    
    colnames(final.df.bayes1)[5] <- sample.ns.mut
    colnames(final.df.bayes1)[8] <- sample.s.mut
    
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
        
    return(final.df.bayes1)    
}


#' Bayesina risk inference model.
#' @description
#'   \code{bayes.driver} function performs by runing Bayesian driver inference model where each observed mutation in gene is taken as proof for being true driver.
#' @param sample.mutations data frame with SNVs and InDels in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item CCF numeric column produce by \code{CCF} function.
#'      \item Damage_score numeric column with values between 0 and 1, where 1 means very damaging SNV/IndDel and 0 not damaging SNV/InDel
#' } 
#' @param bcgr.prob a numeric vector, same lenght as genes (should be same orderd also) which gives probability of gene having somatic mutation in healfy population.
#'          There are functions for obtaining this vector: \code{bcgr}, \code{bcgr.lawrence} and \code{bcgr.combine}.
#' @param genes a vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV/Indel. in given cohort). Default NULL.
#' @param prior.driver a numeric value representing prior probability that random gene is dirver. 
#'          Default is set to \code{length(driver.genes)}/20000, as it assumed there is ~20000 protein goding genes.
#' @param gene.mut.driver a numeric value or named vector representing likelihood that gene is mutated if it is knowen to be driver. 
#'          Gene does not need to be mutated if it is driver, as cancers are heterogenious. Default is set to NULL and driver.genes are considered as drivers.
#' @param driver.genes a character vector of genes which are considered as drivers for this cancer. If NULL then used set is \code{driver.genes.concensus}.
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs/Indels.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Damage_score (optional) integer/numeric value indicating column in \code{sample.mutations} which contain damage score for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param mode a charechter value indicationg how to solve when in one gene-sample pair there are multiple mutations. Options are SUM, MAX and ADVANCE
#' @param epsilon a numeric value. If mode is ADVANCE, epsilone value will be threshold for CCF difference to decide if they are in same or different clone. 
#' @return a data frame with ranked genes by posteriory probability of gene beeing true driver. 
#' Additional columns with usefull info are contained in data frame.
#' @seealso \code{\link{bcgr}}, \code{\link{bcgr.lawrence}}  and \code{\link{bcgr.combine}} for obtaining bcgr.prob variable.
#' @keywords risk
#' @examples 
#' \donttest{
#' # first calculate CCF
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # then somatic background probability
#' bcgr.prob <- bcgr.combine(sample.genes.mutect)
#' # bayes risk model
#' driver.genes <- bayes.driver(sample.genes.mutect, driver = 0.001, gene.mut.driver=1/50) 
#' head(driver.genes)  
#' }
#' @export
bayes.driver <- function(sample.mutations,  bcgr.prob, genes=NULL, prior.driver = NULL, gene.mut.driver=NULL, driver.genes=NULL,
    Variant_Classification=NULL, Hugo_Symbol=NULL, Tumor_Sample_Barcode=NULL, CCF=NULL, Damage_score=NULL , mode='MAX', epsilon=0.05 ) {
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    mode <- toupper(mode)
    if ( !mode %in% c('SUM', 'MAX', 'ADVANCE')) {
        stop("mode mast be or SUM or MAX or ADVACE! ", call. = FALSE)
    }
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mut9ations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(CCF)){
        sample.mutations <- assign.columns(sample.mutations, CCF, "CCF")
    }
    
    if (!is.null(Damage_score)){
        sample.mutations <- assign.columns(sample.mutations, Damage_score, "Damage_score")
    }
      
    
    if(is.null(driver.genes) & is.null(gene.mut.driver)){
        driver.genes <- cDriver::driver.genes.concensus
    } 
    
    if (is.null(prior.driver)){
        if(is.null(gene.mut.driver)){
            prior.driver <- length(driver.genes)/20000     
        } else {
            prior.driver <- length(gene.mut.driver)/20000   
        }
        
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))

    
    
    if (is.null(genes)){
        genes <- as.character(unique(sample.mutations$hugo_symbol))
    }
    
    if(!is.factor(sample.mutations$hugo_symbol)){
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=genes)
    }
    if(!is.factor(sample.mutations$tumor_sample_barcode)){
        sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=unique(sample.mutations$tumor_sample_barcode))
    }   
    
    prior.passenger <- 1 - prior.driver
    
    # take only exonic
    sample.mutations <-  exonic.only(sample.mutations) 
    
    #condtitionals
    if(is.null(gene.mut.driver)){
        num.mut <- nrow(sample.mutations[sample.mutations$hugo_symbol %in% driver.genes,])
        num.pat <- length(unique(as.character(sample.mutations$tumor_sample_barcode)))
        num.canc.genes <- length(unique(driver.genes))
        gene.mut.if.driver <- c()
        gene.mut.if.driver[genes]  <- num.mut/(num.pat*num.canc.genes)
    } else {
        if(is.null(names(gene.mut.driver))){
            gene.mut.if.driver <- c()
            gene.mut.if.driver[genes] <- gene.mut.driver
        } else {
            gene.mut.if.driver[genes] <- gene.mut.driver[genes]   
            gene.mut.if.driver[is.na(gene.mut.if.driver)] <- 0
            names(gene.mut.if.driver) <- genes
        }
        
    }
    
    # make sure that probability that gene is mutated if a driver is bigger or equal then if passenger 
    #gene.mut.if.driver <- mapply(max, gene.mut.if.driver[genes], bcgr.prob[genes])

    gene.mut.if.passenger <- bcgr.prob
    gene.notmut.if.passenger <- 1 - gene.mut.if.passenger

    gene.notmut.if.driver <- 1 - gene.mut.if.driver
    
    # matrices genes vs. samples
    mat.sample.gene.ns <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification != 'Silent',], genes)
    mat.sample.gene.s <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification == 'Silent',], genes)
    if ('damage_score' %in% colnames(sample.mutations) ){
        mat.sample.gene.ns.Damage <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification != 'Silent',], genes, valueCol='damage_score', sample.gene.lim=1, mode=mode, epsilon)
    } else {
        mat.sample.gene.ns.Damage <- matrix(1, nrow=nrow(mat.sample.gene.ns), ncol=ncol(mat.sample.gene.ns))
    }
    mat.sample.gene.ns.ccf <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification != 'Silent',], genes, valueCol='ccf', sample.gene.lim=1, mode=mode, epsilon)
    mat.sample.gene.s.ccf <- pat.vs.genes(sample.mutations[sample.mutations$variant_classification == 'Silent',], genes, valueCol='ccf', sample.gene.lim=1, mode=mode, epsilon)

    # ccf corrected counts
    n <-  rowSums(mat.sample.gene.ns.Damage*mat.sample.gene.ns.ccf)
    m <-  ncol(mat.sample.gene.ns.Damage) - n
    
    
    if (FALSE){
        logxpy <- function(lx,ly) max(lx,ly) + log1p(exp(-abs(lx-ly)))    
 
        # A/A+B
        lA <- log(prior.driver) + n[genes]*log(gene.mut.if.driver[genes] ) + m[genes]*log(gene.notmut.if.driver[genes])
        #A <- log(prior.driver) + n*log(gene.mut.if.driver)  + m*log(gene.notmut.if.driver)
        #print(lA['TP53'])
        lB <- log(prior.passenger) + n[genes]*log(gene.mut.if.passenger[genes]) + m[genes]*log(gene.notmut.if.passenger[genes])
        #B <- log(prior.passenger) + n*log(gene.mut.if.passenger['TP53']) + m*log(gene.notmut.if.passenger['TP53'])
        #print(lB['TP53'])
        
        gene.driver.given.patients <-  exp(lA - 
                                               (lB + log1p( (prior.driver/prior.passenger)*
                                                                ((gene.mut.if.driver[genes]/gene.mut.if.passenger[genes])^n[genes])*
                                                                ((gene.notmut.if.driver[genes]/gene.notmut.if.passenger[genes])^m[genes])  )))
        #print(gene.driver.given.patients['TP53'])
    }
        
    # 
    ind <- (gene.mut.if.driver[genes]^n[genes] == 0 ) | (gene.notmut.if.driver[genes]^m[genes] == 0) |  (gene.mut.if.passenger[genes]^n[genes] == 0) |  (gene.notmut.if.passenger[genes]^m[genes] == 0)
    if (sum(ind)){
        genes1 <- genes[!ind] 
        genes2 <- genes[ind]  
    
        gene.driver.given.patients <- (prior.driver * gene.mut.if.driver[genes1]^n[genes1]*gene.notmut.if.driver[genes1]^m[genes1])/
            ( (prior.driver* gene.mut.if.driver[genes1]^n[genes1]*gene.notmut.if.driver[genes1]^m[genes1])+
                  (prior.passenger* gene.mut.if.passenger[genes1]^n[genes1]*gene.notmut.if.passenger[genes1]^m[genes1]) )  
        
    
        # only to be able to run as Rscript
        requireNamespace(methods)
        
        gene.driver.given.patients2 <- (prior.driver * mpfr(gene.mut.if.driver[genes2], precBits = 128)^n[genes2] * mpfr(gene.notmut.if.driver[genes2], precBits = 128)^m[genes2])/
            ( (prior.driver* mpfr(gene.mut.if.driver[genes2], precBits = 128)^n[genes2]*mpfr(gene.notmut.if.driver[genes2], precBits = 128)^m[genes2])+
                  (prior.passenger* mpfr(gene.mut.if.passenger[genes2], precBits = 128)^n[genes2]*mpfr(gene.notmut.if.passenger[genes2], precBits = 128)^m[genes2]) ) 
        
        gene.driver.given.patients2 <- as.numeric(gene.driver.given.patients2)
        names(gene.driver.given.patients2) <- genes2
        gene.driver.given.patients <- c(gene.driver.given.patients, gene.driver.given.patients2)
        gene.driver.given.patients <- gene.driver.given.patients[genes]
    
    } else {
        gene.driver.given.patients <- (prior.driver * gene.mut.if.driver[genes]^n[genes]*gene.notmut.if.driver[genes]^m[genes])/
            ( (prior.driver* gene.mut.if.driver[genes]^n[genes]*gene.notmut.if.driver[genes]^m[genes])+
                  (prior.passenger* gene.mut.if.passenger[genes]^n[genes]*gene.notmut.if.passenger[genes]^m[genes]) )  
    }

   
    
    #additional columns
    observed.mut.ns <- base::apply(mat.sample.gene.ns,1,function(x) sum(x >0))
    total.mut.ns <- rowSums(mat.sample.gene.ns)    
    observed.mut.s <- base::apply(mat.sample.gene.s,1,function(x) sum(x >0))
    observed.mut.ns.ccf <- rowSums(mat.sample.gene.ns.ccf)
    observed.mut.s.ccf <- rowSums(mat.sample.gene.s.ccf)
    n.indels <- table(sample.mutations[sample.mutations$variant_type %in% c('DEL', 'INS'),'hugo_symbol'], exclude=F)
    
    df.temp <- as.data.frame(cbind(gene.driver.given.patients[genes],n[genes]))
    df.temp <- df.temp[ order(round(df.temp$V1, digits=4), df.temp$V2, decreasing=T), ]
    
    
    #top <- as.character(names(sort(gene.driver.given.patients, decreasing=T)))
    top <- as.character(rownames(df.temp))
    sample.ns.mut <- paste('sample_ns_mut(',length(levels(sample.mutations$tumor_sample_barcode)),')',sep="")
    sample.s.mut <- paste('sample_s_mut(',length(levels(sample.mutations$tumor_sample_barcode)),')',sep="")
    
    final.df.bayes2 <- data.frame('Hugo_Symbol' = top,
                                  'postProb' = round(gene.driver.given.patients[top], digits=4),
                                  'CCF_x_damage' = n[top],
                                  'fold' = gene.driver.given.patients[top]/prior.driver,
                                  'sample_ns_mut' = observed.mut.ns[top],
                                  'total_ns_mut'= total.mut.ns[top],
                                  'nonsilent_CCF_sum' = observed.mut.ns.ccf[top],
                                  'silent'= observed.mut.s[top],
                                  'silent_CCF_sum' = observed.mut.s.ccf[top],                                  
                                  'indels'= as.numeric(n.indels[top]),
                                  'background_mut_prob'=bcgr.prob[top],
                                  'ccf_mean'= base::apply(mat.sample.gene.ns.ccf[top,],1,function(x) mean(x[x!=0])),
                                  'ccf_median'= base::apply(mat.sample.gene.ns.ccf[top,],1,function(x) median(x[x!=0])),
                                  'ccf_sd'= base::apply(mat.sample.gene.ns.ccf[top,],1,function(x) sd(x[x!=0])),
                                  'gene_mut_if_driver'=gene.mut.if.driver[top],
                                  'rank'=1:length(top)
    )

    
    colnames(final.df.bayes2)[5] <- sample.ns.mut
    colnames(final.df.bayes2)[8] <- sample.s.mut
    
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    # correct to NA postProb for genes without any nonsilent mutation
    #final.df.bayes2[final.df.bayes2$total_ns_mut == 0, 'postProb'] <- NA
    #final.df.bayes2 <- final.df.bayes2[order(final.df.bayes2$postProb, decreasing=T, na.last=T), ]
    
    return(final.df.bayes2)
}


#' Combining mulitple methods for ranking genes.
#' @description
#'   \code{combine.ranking} function combines mulitple rankings to one final. All the rankings takes the same importans for final decision.
#' @param rankings a list containing at least two elements, first one data frame output of bayes.risk and second  data frame output of bayes.driver functions.  
#'          All next elements of the list, should be named numeric vectors, where nemas are Hugo_Symbol and value represent any score where higer values are scoring better for gene to be driver.
#' @param genes a vector of genes which were sequenced/analysed. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV/Indel. in given cohort).
#' @param min.mut a numeric value, threshold for filtering genes based on minimum number of nonsilent mutations.
#' @return a data frame with ranked genes by combining different methods ranking. 
#' Additional columns with usefull info are contained in resulting data frame also.
#' @seealso \code{\link{bayes.risk}}, \code{\link{bayes.driver}}   for obtaining bayes.risk.rank and bayes.driver.rank variables.
#' @keywords combine
# @examples 
#' @export
combine.ranking <- function(rankings=list(bayes.risk.rank, bayes.driver.rank), genes = NULL, min.mut = 2 ){
    bayes.risk.rank <- rankings[[1]]
    bayes.driver.rank <- rankings[[2]]
   
    
    if (is.null(genes)){
        if(!'Hugo_Symbol' %in% colnames(bayes.driver.rank)){
            stop(paste("If genes not provided, then there must be column Hugo_Symbol specified in first dataframe in list rankings "), call. = FALSE)
        }
        genes <- as.character(unique(bayes.driver.rank$Hugo_Symbol))
    }
    
    rank.df <- data.frame('Hugo_Symbol' = genes,
                          'bayes_risk' = (bayes.risk.rank[genes,]$postProb),
                          'bayes_driver' = (bayes.driver.rank[genes,]$postProb),
                          'risk_rank' = (bayes.risk.rank[genes,]$rank),
                          'driver_rank' = (bayes.driver.rank[genes,]$rank),
                          'sample_ns_mut' =  bayes.driver.rank[genes,5],
                          'total_ns_mut'= bayes.driver.rank[genes,]$total_ns_mut,
                          'ccf_sum' =  bayes.driver.rank[genes,]$nonsilent_CCF_sum ,
                          'CCFxDamage'= bayes.driver.rank[genes,]$CCF_x_damage,
                          'indels_count'= bayes.driver.rank[genes,]$indels ,                            
                          'mut_rate' = bayes.driver.rank[genes,]$background_mut_prob
    )
    
    colnames(rank.df)[6] <- colnames(bayes.driver.rank)[5]
    
    if (length(rankings) > 2) {
        n <- length(rankings) 
        for (i in 3:n){
            rank.df[,paste('score', i, sep='', collapse='')] <- rankings[[i]][genes]
        }
         
    } 
    
    suppressWarnings(
    rank.df <- rank.df[rank.df$sample_ns_mut >= min.mut & !is.na(rank.df$sample_ns_mut),  ]
    )
    
    rank.df$rank_bayes_risk <- rank( ( rank.df$risk_rank), ties.method="min")
    rank.df$rank_bayes_driver <- rank( ( rank.df$driver_rank), ties.method="min")
    
    if (length(rankings) > 2) {
        n <- length(rankings) 
        for (i in 3:n){
          rank.df[,paste('rank', i, sep='', collapse='')] <- rank( (- rank.df[,paste('score', i, sep='', collapse='')] ), ties.method="min")              
        }       
    } 
    
    #rank.df <- rank.df[with(rank.df, order(rank.bayes.driver, decreasing=F)),]
    if (length(rankings) <= 2) {
        prc.rank <- prcomp(~rank_bayes_risk+rank_bayes_driver,rank.df, center=F)
    } else {
        n <- length(rankings) 
        formula<- paste('~rank_bayes_risk+rank_bayes_driver', paste('rank',3:n, sep='', collapse='+'), sep='+', collapse='')
        prc.rank <- prcomp(as.formula(formula),rank.df, center=F)
        
    }
    
    ind <- prc.rank$rotation[1,1] < 0 
    genes.order <-  names(sort(prc.rank$x[,1], decreasing=ind))
    
    rank.df <- rank.df[genes.order, ]
    rank.df$comb_rank <-  rank( abs(sort((prc.rank$x[,1]), decreasing=ind)), ties.method="min")
    
 
    
    return(rank.df)
} 