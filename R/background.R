######################################################################################
# CRG 
# Hana SUSAK
#------------------------------------------------------------------------------------
# background model functions
######################################################################################



# Function to assign columns in case no header is provided or naming is not standard
# (better practice would be to already have named columns as specified)
# @param sample.mutations - data frame to which column need to be named
# @param column - integer/numeric value for which we are assigning name
# @param columnName - charachter value, new name of the column
assign.columns <- function(sample.mutations, column, columnName){

    if (is.numeric(column)){
        column <- as.integer(column)  
    }
    if (!is.integer(column)){
        stop(paste(column," must be integer (or numeric which will be cast to integer)"), call. = FALSE)
    }
    if (length(column) > 1) {
        stop("column must be a single integer/numeric variable", call. = FALSE)
    }
    colnames(sample.mutations)[column] <- as.character(columnName)
    return(sample.mutations)
}

# Function to take only exonic mutations
# @param sample.mutations - data frame to which filtration is done
# @param Variant_Classification Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV . 
#      Default is NULL value (in this case \code{sample.mutations} should already have this column)
exonic.only <- function(sample.mutations, Variant_Classification=NULL){
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "variant_classification")
    }
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    # take only exonic
    num.del <- nrow(sample.mutations[!sample.mutations$variant_classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',
                                                                                     'Missense_Mutation', 'Nonsense_Mutation', 'Silent', 
                                                                                     'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation'), ])
    if (num.del ){
        warning(paste(num.del, "mutations excluded from analysis as they are classified as one of: 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR , Intron, RNA, Targeted_Region"))
        
        sample.mutations <- sample.mutations[sample.mutations$variant_classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',      
                                                                                            'Missense_Mutation', 'Nonsense_Mutation', 'Silent', 
                                                                                            'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation'), ]
        
    }
    
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    return(sample.mutations)
}



#' Background mutation rate calculated based on silent mutations.
#' @description
#'   \code{bcgr} function calculates  background probability that gene is mutated based on frequency od silent mutations.
#' @param sample.mutations data frame in MAF like format with nonsilent and silent mutations.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item CCF numeric column produce by \code{CCF} function.
#' } 
#' @param genes vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV in given cohort).
#' Default is NULL value and then list of unique genes is takend from \code{sample.mutations}.
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column). 
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @details With assumption of neutral selections, function estimates expected number of nonsilent mutations from observed number of silent mutations.
#'  Na (number of all possible nonsilent substitutions) and Ns (number of all possible silent substitutions) are taken from Lawrence paper.
#'  They are provdied in this package in file lawrence.RData.
#'  When expected number of nonsilent mutations for each gene is known, probability to get nonsilent mutation in each gene is calculated.
#'  This is based on 
#' @return a named numeric vector of probabilites that gene has nonsilent mutation (not caused by cancer).
#' @keywords background
#' @examples 
#' \donttest{
#' # We first need CCF column
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' somatic.background <- bcgr(sample.genes.mutect, length.genes$Hugo_Symbol)
#' head(somatic.background)
#' }
#' @export
bcgr <- function(sample.mutations, genes=NULL, Variant_Classification=NULL, Hugo_Symbol=NULL, Tumor_Sample_Barcode=NULL, CCF=NULL){
   
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
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
    
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    if(!'hugo_symbol' %in% colnames(sample.mutations)){
        stop(paste("There must be column Hugo_Symbol specified in sample.mutations data frame "), call. = FALSE)
    }
    
    if(!'tumor_sample_barcode' %in% colnames(sample.mutations)){
        stop(paste("There must be column Tumor_Sample_Barcode specified in sample.mutations data frame "), call. = FALSE)
    }
    
    if(!'variant_classification' %in% colnames(sample.mutations)){
        stop(paste("There must be column Variant_Classification specified in sample.mutations data frame "), call. = FALSE)
    }
    
    if (is.null(genes)){
        message('In the function bcgr you did not provide list of genes which were sequenced, union of all protein coding genes (first column of all.genes.lengths) and your genes is taken')
        #genes <- unique(sample.mutations$hugo_symbol)
        genes <-  as.character(union(cDriver::all.genes.lengths$Hugo_Symbol, unique(sample.mutations$hugo_symbol)))
        
    }
    
    
    if(!is.factor(sample.mutations$hugo_symbol)){
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=genes)
    }
    if(!is.factor(sample.mutations$tumor_sample_barcode)){
        sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=unique(sample.mutations$tumor_sample_barcode))
    }
    
    sample.mutations <-  exonic.only(sample.mutations)
  
    # estimete number of nonsilent mutations which are high clonal (before cancer occure)
    nonsilent.df <- sample.mutations[sample.mutations$variant_classification != 'Silent' & !is.na(sample.mutations$variant_classification), ]
    silent.df <- sample.mutations[sample.mutations$variant_classification == 'Silent' & !is.na(sample.mutations$variant_classification), ]
    
    # mutations with CCF above  0.8 are considered clonal | taken as maximum likelihood estimate of lambda (poisson distribution)
    #n.healthy <- round( sum( table(nonsilent.df[nonsilent.df$ccf >= 0.85,'tumor_sample_barcode']) )/length(unique(sample.mutations$tumor_sample_barcode)))
    n.healthy <- floor(as.numeric(median( table(nonsilent.df[nonsilent.df$ccf >= 0.85,'tumor_sample_barcode']) )))
    n.healthy <- max(1, n.healthy)
    #silent
    mat.sample.gene.s.ccf <- pat.vs.genes(silent.df, genes, valueCol='ccf', sample.gene.lim=1)
    observed.mut.s.ccf <- rowSums(mat.sample.gene.s.ccf)      
    
    # correction values
    prim.corr.s.zero <- cDriver::lawrence.df[genes,]$N_silent/cDriver::lawrence.df[genes,]$N_nonsilent
    names(prim.corr.s.zero) <- genes
  
    # indicator which genes are not provided in lawrence table
    ind <- names(prim.corr.s.zero[is.na(prim.corr.s.zero)])
        
    avg.corr <- mean(prim.corr.s.zero, na.rm=T) # for missing values, where there is no N_silent and N_nonsilent
    if (sum(is.na(prim.corr.s.zero)) > 0){
        prim.corr.s.zero[(is.na(prim.corr.s.zero))] <- avg.corr      
    }
    #ind <- (!names(observed.mut.s.ccf) %in% (names(prim.corr.s.zero)))
    #prim.corr.s.zero[names(observed.mut.s.ccf)[ind]] <- avg.corr
    
    # calculate corrected (no zeros) somatic silent CCF sums
    observed.mut.s.ccf.cor <- c()
    observed.mut.s.ccf.cor <- sapply(genes, function(x) max(observed.mut.s.ccf[x], prim.corr.s.zero[x], na.rm=T))
    names(observed.mut.s.ccf.cor) <- genes
    
    # estimation of n_a per gene from KaKs 
    mean.NaNs <- mean((cDriver::lawrence.df[genes,'N_nonsilent']/cDriver::lawrence.df[genes,'N_silent']), na.rm=T)
    
    n_a <- observed.mut.s.ccf.cor[genes]*(cDriver::lawrence.df[genes,'N_nonsilent']/cDriver::lawrence.df[genes,'N_silent'])
    names(n_a) <- genes
    if (sum(is.na(n_a)) > 0){
        warning(paste('For ',sum(is.na(n_a)),' genes there was no Na and Ns values (lawrence.df). Mean ratio of Na/Ns is taken for these genes'))
        n_a[is.na(n_a)] <- observed.mut.s.ccf.cor[names(n_a[is.na(n_a)])]*mean.NaNs
    }
    
    n_a.prob <- n_a/sum(n_a)
    
    gene.mutated.if.healthy <-  sapply(n_a.prob, function(x) sum(dbinom(1:n.healthy, n.healthy, x )))

    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    gene.mutated.if.healthy  
    
}

#' Background mutation rate based on Lawrence's background mutation rate estimated on cohort of 12 different cancers.
#' @description
#'   \code{bcgr.lawrence} function calculates  background probability that gene is mutated based on background somatic mutation rate 
#'   provided in Lawrenc paper.
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item CCF numeric column produce by \code{CCF} function.
#' } 
#' @param genes vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV in given cohort).
#' Default is NULL value and then list of unique genes is takend from \code{sample.mutations}.
#' @param lengthGenes numeric vector of lengths (suquenced) for \code{genes}. Should be given in same order as variable genes.
#' Default is NULL value and then Length of genes is taken from data set \code{length.genes} (form package \code{cDriver}) as column Length. 
#' If gene is not found in this data frame, then median value is taken from listed genes in this data frame.
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @details Function uses \code{sample.mutations} just to estimate number of nonsilent mutations per patient.
#' @return a named numeric vector of probabilites that gene has nonsilent mutation (not caused by cancer).
#' @keywords Lawrence
#' @examples 
#' bcgr.L <- bcgr.lawrence(sample.genes.mutect, length.genes$Hugo_Symbol, length.genes$Coverd_len)
#' head(bcgr.L)
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23770567}.
#' @export
bcgr.lawrence <- function(sample.mutations, genes=NULL, lengthGenes=NULL, Variant_Classification=NULL, Hugo_Symbol=NULL, Tumor_Sample_Barcode=NULL, CCF=NULL){
    #load('data/lawrence.RData')    
    #ow <- options("warn")
    
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
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
    
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    
    if (is.null(lengthGenes)){
        message('You did not provide lenght of sequenced genes, so sum of lenghts of all exons for each gene is taken (all.genes.lengths).')
        
        lengthGenes <- cDriver::all.genes.lengths$Length
        names(lengthGenes) <- cDriver::all.genes.lengths$Hugo_Symbol
        
    } else if(  is.null(names(lengthGenes))) {
        if (!is.null(genes)){
            names(lengthGenes) <- genes
        } else {
            stop('When providing length of genes, you need also to provide their names (or as named vector or in genes paramether). ')
        }
            
    }
      
    if (is.null(genes)){
        message('In the function bcgr.lawrence you did not provide list of genes which were sequenced, union of all protein coding genes (first column of all.genes.lengths) and your genes is taken')
        
        genes <- as.character(union(cDriver::all.genes.lengths$Hugo_Symbol, unique(sample.mutations$hugo_symbol)))
    }
    
    if(!is.factor(sample.mutations$hugo_symbol)){
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=genes)
    }
    if(!is.factor(sample.mutations$tumor_sample_barcode)){
        sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=unique(sample.mutations$tumor_sample_barcode))
    }
    
    # take only exonic
    suppressWarnings(    sample.mutations <-  exonic.only(sample.mutations) )

  
    ## ovde moras napraviti neki warning, vidi koliko gena ima duzinu. 
    ## ako nema daj prosjecnu ali i warning
    
    # From lawrenc
    mut.rate <- cDriver::lawrence.df[genes,]$noncoding_mutation_rate
    names(mut.rate) <- genes
    if(sum(is.na(mut.rate))){
        proc.miss.mut.rate <- sum(is.na(mut.rate))/length(genes)
        warning("In your list of genes there are genes with unkown mutation rate (lawrence estimates). Median value is taken fot those. 
                From ",length(genes), " genes there is no mutation rate provided for ", sum(is.na(mut.rate))," \n")
    }
    mut.rate[is.na(mut.rate)] <- median(mut.rate, na.rm=T)
    
    bcgr <- c()
    bcgr.prob <- c()
    
    if(sum(is.na(lengthGenes[genes]))){
        # to be sure of order and to add lengths which are in genes but not ine lengthGenes
        lengthGenes <-  lengthGenes[genes]
        names(lengthGenes) <- genes
        miss.len <- sum(is.na(lengthGenes))
        lengthGenes[is.na(lengthGenes)] <- median(lengthGenes, na.rm=T)
        if (miss.len){
            proc.mis <- miss.len/length(genes)
            #options(warn = w); cat("\n warn =", w, "\n")
            warning("In your list of genes there are genes with unkown lenght. Medain value is taken for those. 
                From ",length(genes), " genes there is no length provided for ",miss.len," \n")
        }
    } else {
        lengthGenes <-  lengthGenes[genes]
        names(lengthGenes) <- genes
        
    }
    
    Ns.Na <- cDriver::lawrence.df[genes,'N_silent']/cDriver::lawrence.df[genes,'N_nonsilent']
    if (sum(is.na(Ns.Na)) > 0){
        warning(paste('For ',sum(is.na(Ns.Na)),' genes there was no Na and Ns values (lawrence.df). Mean ratio of Ns/Na is taken for these genes'))
        Ns.Na[is.na(Ns.Na)] <- mean(Ns.Na, na.rm=T)
    }
   
    
    bcgr[genes] <- (mut.rate[genes]*lengthGenes[genes])/(1+Ns.Na)
    
    bcgr.prob[genes] <- bcgr[genes] / sum(bcgr[genes])
    
  
    
    
    # estimete number of nonsilent mutations which are high clonal (before cancer occure)
    nonsilent.df <- sample.mutations[sample.mutations$variant_classification != 'Silent' & !is.na(sample.mutations$variant_classification), ]
    # mutations with CCF above  0.8 are considered clonal | taken as maximum likelihood estimate of lambda (poisson distribution)
    #n.healthy <- round( sum( table(nonsilent.df[nonsilent.df$ccf >= 0.85,'tumor_sample_barcode']) )/length(unique(sample.mutations$tumor_sample_barcode)))
    n.healthy <- floor(as.numeric(median( table(nonsilent.df[nonsilent.df$ccf >= 0.8,'tumor_sample_barcode']) )))
    
    gene.mutated.if.healthy <-  sapply(bcgr.prob, function(x) sum(dbinom(1:n.healthy, n.healthy, x )))
    #warnings()
    #options(ow)  # reset
    
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    
    gene.mutated.if.healthy  
}

#' Combining two somatic background mutation probability, outputs of functions bcgr.lawrence and bcgr
#' @description
#'   \code{bcgr.combine} function calculates both somatic background mutation probabilities (using bcgr.lawrence and bcgr functions) and take average value for each gene
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item CCF numeric column produce by \code{CCF} function.
#' } 
#' @param genes vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV in given cohort).
#' Default is NULL value and then list of unique genes is takend from \code{sample.mutations}.
#' @param lengthGenes numeric vector of lengths (suquenced) for \code{genes}. Should be given in same order as variable genes.
#' #' Default is NULL value and then Length of genes is taken from data set \code{length.genes} (form package \code{cDriver}) as column Length. 
#' If gene is not found in this data frame, then median value is taken from listed genes in this data frame.
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#'      Column with this name should not already exist in \code{sample.mutations}.
#' @return a named numeric vector of probabilites that gene has nonsilent mutation (not caused by cancer).
#' @keywords Lawrence
#' @examples
#' # First we need to calculate CCF
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # Calculate somatic nonsilent background mutation probability
#'  background <- bcgr.combine(sample.genes.mutect, length.genes$Hugo_Symbol, length.genes$Coverd_len)
#'  head(background)
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23770567}.
#' @export
bcgr.combine <- function(sample.mutations, genes=NULL, lengthGenes=NULL, Variant_Classification=NULL, Hugo_Symbol=NULL, Tumor_Sample_Barcode=NULL, CCF=NULL){
    na.p.lawrence <- bcgr.lawrence(sample.mutations, genes=genes, lengthGenes=lengthGenes, Variant_Classification=Variant_Classification, Hugo_Symbol=Hugo_Symbol, Tumor_Sample_Barcode=Tumor_Sample_Barcode, CCF=CCF)
    na.p.silent <- bcgr(sample.mutations, genes=genes, Variant_Classification=Variant_Classification, Hugo_Symbol=Hugo_Symbol, Tumor_Sample_Barcode=Tumor_Sample_Barcode, CCF=CCF)
    if (is.null(genes)){
        # make it not sensitive to lower/upper case in column names
        original.col.names <- colnames(sample.mutations)
        num.col <- ncol(sample.mutations)
        colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
        
        #genes <- unique(sample.mutations$hugo_symbol)
        genes <- as.character(union(cDriver::all.genes.lengths$Hugo_Symbol, unique(sample.mutations$hugo_symbol)))
        
        #return original names
        colnames(sample.mutations)[1:num.col] <- original.col.names
        
    }
    #together <- as.data.frame(cbind(na.p.lawrence[genes], na.p.silent[genes]))
    #names(together) <- c('noncoding_lawrence','nonsilent_kaks')
    #gene.mutated.if.healthy <- apply(as.matrix(together),1,mean)
    gene.mutated.if.healthy <- (na.p.lawrence[genes] + na.p.silent[genes])/2
    names(gene.mutated.if.healthy) <- genes
    return(gene.mutated.if.healthy)
}

