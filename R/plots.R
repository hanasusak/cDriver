######################################################################################
# CRG 
# Hana SUSAK
#------------------------------------------------------------------------------------
# plotc for cDriver
######################################################################################




#' Plot mutations per sample distribution.
#' @description
#'   \code{plot.samplesMut} Function to plot samples mutation counts.  
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item Variant_Type columns pecifed by MAF format; is mutation SNV or InDel
#' } 
#' @param indels a boolean value indicating should indels be counted. By default it is True.
#' @param silent a boolean value indicating should silent mutations be counted. By default it is True.
#' @param fill a boolean value indicating should plot only represent proportion between 3 types of mutations, not counts. By default it is False.
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs/Indels.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Type (optional) integer/numeric value indicating column in \code{sample.mutations} which contains indormation if mutations is SNV or InDel .
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @return ggplot2 object 
#' @examples 
#' \donttest{
#' # plot sample's mutations , all of them
#' plotSamplesMut(sample.genes.mutect)
#' # plot proportion of silent and nonsilent, without indels
#' plotSamplesMut(sample.genes.mutect, indels=FALSE, fill=TRUE)
#' }
#' @export
plotSamplesMut <- function(sample.mutations, indels=TRUE, silent=TRUE, fill=FALSE,  Tumor_Sample_Barcode=NULL, Hugo_Symbol=NULL, Variant_Classification=NULL, Variant_Type=NULL ) {
    
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    
    if (!is.null(Variant_Type)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Type, "Variant_Type")
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mutations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    
    if(!is.factor(sample.mutations$Tumor_Sample_Barcode)){
        sample.mutations$Tumor_Sample_Barcode <- factor(sample.mutations$Tumor_Sample_Barcode, levels=unique(sample.mutations$Tumor_Sample_Barcode))
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    # take only exonic
    suppressWarnings( sample.mutations <-  exonic.only(sample.mutations))
    
   
    
    if (indels){
        #indel
        sample.gene.indels <- sample.mutations[sample.mutations$variant_type %in% c('INS', 'DEL' ),]
        sample.gene.indels$typeMut <- 'Indel'
        ord.ind <- names(sort(table(sample.gene.indels$tumor_sample_barcode),decreasing=T))
    }  else {
        sample.gene.indels <- c()
    }   
    if (silent){
        #silent
        sample.gene.silent <- sample.mutations[sample.mutations$variant_type %in% c('SNP', 'DNP','TNP' ) & sample.mutations$variant_classification == "Silent",]
        sample.gene.silent$typeMut <-  'Silent'
        if (indels){
            sample.gene.silent$tumor_sample_barcode <- factor(sample.gene.silent$tumor_sample_barcode, levels=ord.ind)
        }
        ord.sil <- names(sort(table(sample.gene.silent$tumor_sample_barcode),decreasing=F))
    } else {
        sample.gene.silent <- c()
    }
    
    #nonsilent
    sample.gene.non.silent <- sample.mutations[sample.mutations$variant_type %in% c('SNP', 'DNP','TNP' ) & 
                                                sample.mutations$variant_classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 
                                                                                            'Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site', 
                                                                                            'Translation_Start_Site', 'Nonstop_Mutation'),]
    sample.gene.non.silent$typeMut <- 'NonSilent'
    
    #all
    sample.gene.all <- rbind(sample.gene.non.silent, sample.gene.silent, sample.gene.indels)
    if (silent){
        sample.gene.all$tumor_sample_barcode <- factor(sample.gene.all$tumor_sample_barcode, levels=ord.sil)      
    } else if (indels) {
        sample.gene.all$tumor_sample_barcode <- factor(sample.gene.all$tumor_sample_barcode, levels=ord.ind)      
    }
    
    order <- names(sort(table(sample.gene.all$tumor_sample_barcode),decreasing=T))
    sample.gene.all$tumor_sample_barcode <- factor(sample.gene.all$tumor_sample_barcode, levels=order)
    
    max.y <- max(table(sample.gene.all$tumor_sample_barcode))
    
    if (fill) {
        p1 <- ggplot(data=sample.gene.all, aes_string('tumor_sample_barcode', fill='typeMut')) + geom_bar(position='fill') +
            scale_y_continuous( name='Proportions', expand = c(0, 0))  
    } else {
        p1 <- ggplot(data=sample.gene.all, aes_string('tumor_sample_barcode', fill='typeMut')) + geom_bar(position='stack') + 
            scale_y_discrete(name='Count', breaks=round(seq(0, max.y, length.out=10)), expand = c(0, 0))  
    }
    p1 <- p1 +  scale_x_discrete(expand = c(0, 0), drop=FALSE) +
        theme_bw(base_size = 10) + 
        scale_fill_manual(values=c( "orange4",'orange1',"steelblue3")) + 
        xlab('Sample ID' ) +
        theme( axis.text.x = element_text(size = 10 *0.5, angle = 90, hjust = 0.8, vjust=0.5, colour = "grey50")) + 
        theme(axis.title.x=element_text(angle=0, vjust=+0.5, size = 10 )) +
        guides(fill = guide_legend(title = "Mutation type"))
        #scale_fill_discrete(guide = guide_legend(title = "Mutation type"))
        #theme( axis.text.x =element_blank()) + 
        #theme(axis.title.x=element_blank())   
  
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    
    p1
}
    



#' Plot mutations pattern frequncies.
#' @description
#'   \code{plotMutChange} Function to plot frequency of mutation changes.  
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item Variant_Type columns pecifed by MAF format; is mutation SNV or InDel
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Reference_Allele column specifed by MAF format, represents referent allele
#'      \item Tumor_Seq_Allele2 column specifed by MAF format, represents alternative allele
#' } 
#' @param silent a boolean value indicating should silent mutations be counted. By default it is True.
#' @param fill a boolean value indicating should plot only represent proportion between 3 types of mutations, not counts. By default it is True
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Type (optional) integer/numeric value indicating column in \code{sample.mutations} which contains indormation if mutations is SNV or InDel .
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Reference_Allele (optional) integer/numeric value indicating column in \code{sample.mutations} which contain reference alleles. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Tumor_Seq_Allele2 (optional) integer/numeric value indicating column in \code{sample.mutations} which contain alternative alleles. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @return ggplot2 object 
#' @examples 
#' \donttest{
#' # plot SNVs changes patterns
#' plotMutChange(sample.genes.mutect)
#' # plot SNVs changes patterns only for non silent mutaions
#' plotMutChange(sample.genes.mutect, silent=FALSE, fill=FALSE)
#' }
#' @export
plotMutChange <- function(sample.mutations, silent=TRUE, fill=TRUE,  Tumor_Sample_Barcode=NULL, Variant_Type=NULL, Variant_Classification=NULL, Reference_Allele=NULL, Tumor_Seq_Allele2=NULL ) {
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(Variant_Type)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Type, "Variant_Type")
    }
    
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    
    if (!is.null(Reference_Allele)){
        sample.mutations <- assign.columns(sample.mutations, Reference_Allele, "Reference_Allele")
    }
    
    if (!is.null(Tumor_Seq_Allele2)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Seq_Allele2, "Tumor_Seq_Allele2")
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    
    # take only exonic
    suppressWarnings( sample.mutations <-  exonic.only(sample.mutations) )
    
  
    if(!is.factor(sample.mutations$tumor_sample_barcode)){
        sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=unique(sample.mutations$tumor_sample_barcode))
    }
    
    if(!silent){
        sample.mutations <- sample.mutations[sample.mutations$variant_classification != 'Silent',]
    }
    
    
    
    sample.mutations <- sample.mutations[sample.mutations$variant_type %in% c('SNP', 'DNP','TNP' ),]
    sample.mutations$Change <- apply(sample.mutations, 1, function(x)paste(x['reference_allele'], x['tumor_seq_allele2'], sep=';') )
    
    sample.mutations[sample.mutations$Change == 'G;T','Change'] <- 'C;A'
    sample.mutations[sample.mutations$Change == 'G;C','Change'] <- 'C;G'
    sample.mutations[sample.mutations$Change == 'G;A','Change'] <- 'C;T'
    sample.mutations[sample.mutations$Change == 'A;G','Change'] <- 'T;C'
    sample.mutations[sample.mutations$Change == 'A;T','Change'] <- 'T;A'
    sample.mutations[sample.mutations$Change == 'A;C','Change'] <- 'T;G'
    sample.mutations$Change <- unlist(lapply(sample.mutations$Change, function(x) gsub(';','->',x)))
    sample.mutations$Change <- factor((sample.mutations$Change), levels=c( 'C->T', 'T->C', 'C->A', 'C->G', 'T->A','T->G'))
    
    sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=names(sort(table(sample.mutations$tumor_sample_barcode), decreasing=T) ))
    
    
    p2 <- ggplot(sample.mutations, aes_string('tumor_sample_barcode', fill='Change')) + geom_bar(position=position_fill()) + theme_bw(base_size = 10) +
        scale_y_continuous(name='Proportion',expand=c(0,0)) + 
        scale_x_discrete(name='Sample ID',drop=FALSE) + scale_fill_brewer(palette="Set2") +
        theme( axis.text.x = element_text(size = 10 *0.5, angle = 90, hjust = 0.8, vjust=0.5, colour = "grey50")) + 
        theme(axis.title.x=element_text(angle=0, vjust=+0.5, size = 10 )) 
    
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
     
    p2    

}





#' Plot boxplot of mutations for top genes.
#' @description
#'   \code{boxplotCCF.mutations} Function to plot CCF's of all mutations for top genes.  
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item CCF numeric column produce by \code{CCF} function.
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Variant_Type columns pecifed by MAF format; is mutation SNV or InDel
#' } 
#' @param result.df a data frame with Hugo_Symbol column. Genes in this column should be ordered by importance.
#' @param topGenes a numeric/integer value; How meny top genes from results will be ploted.
#' @param silent a boolean value indicating should silent mutations be counted. By default it is True.
#' @param indels a boolean value indicating should indels be counted. By default it is True.
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs/Indels.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Type (optional) integer/numeric value indicating column in \code{sample.mutations} which contains indormation if mutations is SNV or InDel .
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param color a numeric or character value indicating column to color the samples.
#' @param shape a numeric or character value indicating column to use as shape for the samples.
#' @return ggplot2 object 
#' @examples 
#' \donttest{
#' #get CCF column
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # get background
#' bcgr.prob <- bcgr.combine(sample.genes.mutect, length.genes$Hugo_Symbol, length.genes$Coverd_len)
#' # get ranking
#' df1 <- bayes.risk(sample.genes.mutect, bcgr.prob, prior.sick = 0.00007) 
#' df2 <- bayes.driver(sample.genes.mutect, bcgr.prob,  prior.driver = 0.001) 
#' df.final <- combine.ranking(list(df1, df2),  min.mut = 2 )
#' # plot boxplot of CCF for top genes
#' boxplotCCF(sample.mutations=sample.genes.mutect, result.df=df.final)
#' }
#' @export
boxplotCCF.mutations <- function(sample.mutations, result.df, topGenes=20, silent=FALSE, indels=TRUE, 
                       Tumor_Sample_Barcode=NULL, Hugo_Symbol=NULL,  CCF=NULL, Variant_Classification=NULL,  Variant_Type = NULL, color=NULL, shape=NULL ) {
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mutations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    if (!is.null(CCF)){
        sample.mutations <- assign.columns(sample.mutations, CCF, "CCF")
    }
    
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    

    
    if (!is.null(Variant_Type)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Type, "Variant_Type")
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    
    if (!is.null(color) & is.character(color)){
        color <- tolower(color)
    }
    if (!is.null(shape) & is.character(shape)){
        shape <- tolower(shape)
    }
    # take only exonic
    suppressWarnings( sample.mutations <-  exonic.only(sample.mutations) )
    
    
    if(!'Hugo_Symbol' %in% colnames(result.df)){
        stop('There must be Hugo_Symbol column in result.df data frame!')
    }
    
    if(!is.factor(sample.mutations$hugo_symbol)){
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=unique(sample.mutations$hugo_symbol))
    }
    
    plotGenes <- as.character(result.df$Hugo_Symbol[1:min(topGenes,nrow(result.df))])
    sample.mutations <- sample.mutations[sample.mutations$hugo_symbol %in% plotGenes, ]
   
   if(!silent){
        sample.mutations <- sample.mutations[sample.mutations$variant_classification != 'Silent',]
    }
    if(!indels){
        sample.mutations <- sample.mutations[! sample.mutations$variant_type %in% c('DEL', 'INS' ),]
    }
    
   order.gene <- as.table(by(sample.mutations[,'ccf'],sample.mutations[,'hugo_symbol'],median))
   sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol , levels=names(sort(order.gene, decreasing=T)))
    
   if (is.null(color) & is.null(shape)){
        p <- ggplot(sample.mutations, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
            xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') + 
            geom_jitter(size=2.5, alpha=0.65, position = position_jitter(width = .3, height=0))
    } 
   #just color 
   if (!is.null(color) & is.null(shape)){
      if((is.numeric(color) | is.integer(color)) & color <= ncol(sample.mutations) ){
        col.var <- colnames(sample.mutations)[color]  
      } else if (is.character(color) & color %in% colnames(sample.mutations)){
          col.var <- color 
      } else {
          stop('color need to be or numeric indicator of column in sample.mutations data frame or charater with column name.')
      }
       p <- ggplot(sample.mutations, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
           xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') + 
           geom_jitter(aes_string(color=col.var),size=2.5, alpha=0.65, position = position_jitter(width = .3, height=0)) 
       
    }    
    # just shape
    if (!is.null(shape) & is.null(color)){
        if((is.numeric(shape) | is.integer(shape)) & shape <= ncol(sample.mutations) ){
            sh.var <- colnames(sample.mutations)[shape]  
        } else if (is.character(shape) & shape %in% colnames(sample.mutations)){
            sh.var <- shape 
        } else {
            stop('shape need to be or numeric indicator of column in sample.mutations data frame or charater with column name.')
        }
        n.sh <- length(unique(sample.mutations[,sh.var]))
        if (n.sh <= 11 ) {
            sh.seq <- 15:(n.sh+14)
        } else {
            sh.seq <- c(15:25,1:(n.sh-11))
        }
        p <- ggplot(sample.mutations, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
            xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') + 
            geom_jitter(aes_string(shape=sh.var),size=2.5, alpha=0.65, position = position_jitter(width = .3, height=0)) +
            scale_shape_manual(values=sample(sh.seq, replace=F) ) 
            
        
    }
   if (!is.null(shape) & !is.null(color)){
       if((is.numeric(shape) | is.integer(shape)) & shape <= ncol(sample.mutations) ){
           sh.var <- colnames(sample.mutations)[shape]  
       } else if (is.character(shape) & shape %in% colnames(sample.mutations)){
           sh.var <- shape 
       } else {
           stop('shape need to be or numeric indicator of column in sample.mutations data frame or charater with column name.')
       }
       n.sh <- length(unique(sample.mutations[,sh.var]))
       if (n.sh <= 11 ) {
           sh.seq <- 15:(n.sh+14)
       } else {
           sh.seq <- c(15:25,1:(n.sh-11))
       }
       if((is.numeric(color) | is.integer(color)) & color <= ncol(sample.mutations) ){
           col.var <- colnames(sample.mutations)[color]  
       } else if (is.character(color) & color %in% colnames(sample.mutations)){
           col.var <- color 
       } else {
           stop('color need to be or numeric indicator of column in sample.mutations data frame or charater with column name.')
       }
       p <- ggplot(sample.mutations, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
           xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') + 
           geom_jitter(aes_string(color=(col.var), shape=(sh.var)),size=2.5, alpha=0.65, position = position_jitter(width = .15, height=0)) +#, drop=FALSE
           scale_shape_manual(values=sample(sh.seq, replace=F)) 
   }
   
   
   #return original names
   colnames(sample.mutations)[1:num.col] <- original.col.names
   
   
   return(p)
}
    



#' Plot boxplot of patients mutations CCF's for top genes.
#' @description
#'   \code{boxplotCCF.patients} Function to plot CCF's of  aggregated patient-gene mutatuons, for top genes.  
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item CCF numeric column produce by \code{CCF} function.
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Variant_Type columns pecifed by MAF format; is mutation SNV or InDel
#' } 
#' @param result.df a data frame with Hugo_Symbol column. Genes in this column should be ordered by importance.
#' @param topGenes a numeric/integer value; How meny top genes from results will be ploted.
#' @param silent a boolean value indicating should silent mutations be counted. By default it is True.
#' @param indels a boolean value indicating should indels be counted. By default it is True.
#' @param mode a charechter value indicationg how to solve when in one gene-sample pair there are multiple mutations. Options are SUM, MAX and ADVANCE
#' @param epsilon a numeric value. If mode is ADVANCE, epsilone value will be threshold for CCF difference to decide if they are in same or different clone. 
#' @param sample.gene.lim a numeric value specifying upper limit when summing is perfomed for gene-sample pair
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs/Indels.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param CCF (optional) integer/numeric value indicating column in \code{sample.mutations} which have cancer cell fraction information for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Type (optional) integer/numeric value indicating column in \code{sample.mutations} which contains indormation if mutations is SNV or InDel .
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param color (optional) charachter value indicating column in \code{sample.mutations} to used as color dots.
#'      Default is NULL value and grey dots gonna be ploted. This variable should be unique by patients.
#' @param shape (optional) charachter value indicating column in \code{sample.mutations} to used as shape for dots .
#'      Default is NULL value and grey dots gonna be ploted. This variable should be unique by patients.
#' @return ggplot2 object 
#' @examples 
#' \donttest{
#' #get CCF column
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # get background
#' bcgr.prob <- bcgr.combine(sample.genes.mutect)
#' # get ranking
#' df1 <- bayes.risk(sample.genes.mutect, bcgr.prob,  prior.sick = 0.00007) 
#' df2 <- bayes.driver(sample.genes.mutect, bcgr.prob,  prior.driver = 0.001) 
#' df.final <- combine.ranking(list(df1, df2),  min.mut = 2 )
#' # plot boxplot of CCF for top genes
#' boxplotCCF(sample.mutations=sample.genes.mutect, result.df=df.final)
#' }
#' @export
boxplotCCF.patients <- function(sample.mutations, result.df, topGenes=20, silent=FALSE, indels=TRUE, mode='MAX', epsilon=0.05, sample.gene.lim=1,
                                 Tumor_Sample_Barcode=NULL, Hugo_Symbol=NULL,  CCF=NULL, Variant_Classification=NULL,  Variant_Type = NULL, color=NULL, shape=NULL) {
    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    mode <- toupper(mode)
    if ( !mode %in% c('SUM', 'MAX', 'ADVANCE')) {
        stop("mode mast be or SUM or MAX or ADVACE! ", call. = FALSE)
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mutations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    if (!is.null(CCF)){
        sample.mutations <- assign.columns(sample.mutations, CCF, "CCF")
    }
    
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    
    
    if (!is.null(Variant_Type)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Type, "Variant_Type")
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    # take only exonic
    suppressWarnings( sample.mutations <-  exonic.only(sample.mutations ))
    
    if ( (!is.null(color)) & is.character(color)  ) { 
        if (tolower(color) %in% colnames(sample.mutations)){
            color <- tolower(color)  
        } else {
            stop('color variable must be one of the column names in sample.mutations data frame')
        }       
    } else if (!is.null(color)){
        stop('color variable must be character type')
    }  
    
    if ( (!is.null(shape)) & is.character(shape)   ) { 
        if (tolower(shape) %in% colnames(sample.mutations)){
            shape <- tolower(shape)  
        } else {
            stop('shape variable must be one of the column names in sample.mutations data frame')
        }      } else if (!is.null(shape)){
        stop('shape variable must be character type')
    } 
    
    
    if(!'Hugo_Symbol' %in% colnames(result.df)){
        stop('There must be Hugo_Symbol column in result.df data frame!')
    }
    
    if(!is.factor(sample.mutations$hugo_symbol)){
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=unique(sample.mutations$hugo_symbol))
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
    
    plotGenes <- as.character(result.df$Hugo_Symbol[1:min(topGenes,nrow(result.df))])
    sample.mutations <- sample.mutations[sample.mutations$hugo_symbol %in% plotGenes, ]
    
    if(!silent){
        sample.mutations <- sample.mutations[sample.mutations$variant_classification != 'Silent',]
    }
    if(!indels){
        sample.mutations <- sample.mutations[! sample.mutations$variant_type %in% c('DEL', 'INS' ),]
    }
    
 
    # AGGREGATE SNPs which are in same sample-gene pair
    
    if (mode == 'SUM' ){
        # if we sum values of column CCF, for same sample-gene pairs  
        if (is.null(sample.gene.lim)){
            sample.gene.lim <- Inf
        }
        f <- paste( "ccf ~  hugo_symbol + tumor_sample_barcode",  collapse="")
        sample.mutations2 <- aggregate(as.formula(f) , data=sample.mutations, function(x) min(sum(x),sample.gene.lim))         
        if (!is.null(color)){
            df.temp <- sample.mutations[,c(color,'tumor_sample_barcode')]
            df.temp <- df.temp[!duplicated(df.temp), ]
            first.n <- nrow(sample.mutations2)
            sample.mutations2 <- merge(x=sample.mutations2, y=df.temp, by='tumor_sample_barcode',all.x=T)
            second.n <- nrow(sample.mutations2)
            if(first.n!=second.n){
                warning('You choose to color by paramether which is not unique by patients!')
            }
        }
        if (!is.null(shape)){
            df.temp <- sample.mutations[,c(shape,'tumor_sample_barcode')]
            df.temp <- df.temp[!duplicated(df.temp), ]
            first.n <- nrow(sample.mutations2)
            sample.mutations2 <- merge(x=sample.mutations2, y=df.temp, by='tumor_sample_barcode',all.x=T)
            second.n <- nrow(sample.mutations2)
            if(first.n!=second.n){
                warning('You choose to shape by paramether which is not unique by patients!')
            }
        }
    } else if (mode == 'MAX' ) {
        # if we take max value of column CCF, for same sample-gene pairs    
        if (is.null(sample.gene.lim)){
            sample.gene.lim <- Inf
        }
        f <- paste("ccf ~  hugo_symbol + tumor_sample_barcode",  collapse="")
        sample.mutations2 <- aggregate(as.formula(f) , data=sample.mutations, function(x) min(max(x, na.rm=T),sample.gene.lim))         
        if (!is.null(color)){
            df.temp <- sample.mutations[,c(color,'tumor_sample_barcode')]
            df.temp <- df.temp[!duplicated(df.temp), ]
            first.n <- nrow(sample.mutations2)
            sample.mutations2 <- merge(x=sample.mutations2, y=df.temp, by='tumor_sample_barcode',all.x=T)
            second.n <- nrow(sample.mutations2)
            if(first.n!=second.n){
                warning('You choose to color by paramether which is not unique by patients!')
            }
        }        
        if (!is.null(shape)){
            df.temp <- sample.mutations[,c(shape,'tumor_sample_barcode')]
            df.temp <- df.temp[!duplicated(df.temp), ]
            first.n <- nrow(sample.mutations2)
            sample.mutations2 <- merge(x=sample.mutations2, y=df.temp, by='tumor_sample_barcode',all.x=T)
            second.n <- nrow(sample.mutations2)
            if(first.n!=second.n){
                warning('You choose to shape by paramether which is not unique by patients!')
            }
        }
    } else if  (mode == 'ADVANCE'  ){
        # if we take max or sum value of column CCF, for same sample-gene pairs, depending if in same or different clone   
        if (is.null(sample.gene.lim)){
            sample.gene.lim <- Inf
        }
        #sample.gene.dt <- sample.gene.dt[order(sample.gene.dt$Hugo_Symbol, sample.gene.dt$Tumor_Sample_Barcode, -sample.gene.dt$CCF),]
        f <- paste( "ccf ~  hugo_symbol + tumor_sample_barcode",  collapse="")
        sample.mutations2 <- aggregate(as.formula(f) , data=sample.mutations,  function(x) min(advance.ccf.calc(x, epsilon),sample.gene.lim))          
        if (!is.null(color)){
            df.temp <- sample.mutations[,c(color,'tumor_sample_barcode')]
            df.temp <- df.temp[!duplicated(df.temp), ]
            first.n <- nrow(sample.mutations2)
            sample.mutations2 <- merge(x=sample.mutations2, y=df.temp, by='tumor_sample_barcode',all.x=T)
            second.n <- nrow(sample.mutations2)
            if(first.n!=second.n){
                warning('You choose to color by paramether which is not unique by patients!')
            }
        }
        if (!is.null(shape)){
            df.temp <- sample.mutations[,c(shape,'tumor_sample_barcode')]
            df.temp <- df.temp[!duplicated(df.temp), ]
            first.n <- nrow(sample.mutations2)
            sample.mutations2 <- merge(x=sample.mutations2, y=df.temp, by='tumor_sample_barcode',all.x=T)
            second.n <- nrow(sample.mutations2)
            if(first.n!=second.n){
                warning('You choose to shape by paramether which is not unique by patients!')
            }
        }
    } 
    
    order.gene <- as.table(by(sample.mutations2[,'ccf'],sample.mutations2[,'hugo_symbol'],median))
    sample.mutations2$hugo_symbol <- factor(sample.mutations2$hugo_symbol , levels=names(sort(order.gene, decreasing=T)))
    
    if ( (!is.null(color)) &  (!is.null(shape))){
        p <- ggplot(sample.mutations2, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
            xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') +
            geom_jitter(aes_string(color=(color), shape=(shape)),size=3, alpha=0.75, position = position_jitter(width = .3, height=0))
    }  else if ( !is.null(color) ){
        p <- ggplot(sample.mutations2, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
            xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') +
            geom_jitter(aes_string(color=(color)),size=3, alpha=0.75, position = position_jitter(width = .3, height=0))
    }  else if ( !is.null(shape) ){
        p <- ggplot(sample.mutations2, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
            xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') +
            geom_jitter(aes_string(shape=(shape)),size=3, alpha=0.75, position = position_jitter(width = .3, height=0))
    } else {
        p <- ggplot(sample.mutations2, aes_string('hugo_symbol','ccf', label='tumor_sample_barcode')) + geom_boxplot(outlier.size=0) + theme_bw()  +
            xlab('Top  Genes') + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust=1)) + ylab('Cancer Cell Fraction') +
            geom_jitter(size=3, alpha=0.75, position = position_jitter(width = .3, height=0))
    }
    
   
    
    #return original names
    colnames(sample.mutations)[1:num.col] <- original.col.names
    
    return(p)
}




#' Plot heatmap for genes-sample pairs .
#' @description
#'   \code{plotStaircase} Function to plot heatmap for genes-sample pairs orderd by genes which explain most of the patients .  
#' @param sample.mutations data frame in MAF like format.  
#' Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item Tumor_Sample_Barcode column specifed by MAF format, reporting for each SNV in wich patient was found. 
#'      \item Hugo_Symbol column specifed by MAF format, which reports gene for each SNV.
#'      \item Variant_Classification column specifed by MAF format, used to distinguish between silent and nonsilent SNVs
#'      \item Variant_Type columns pecifed by MAF format; is mutation SNV or InDel
#' } 
#' @param result.df a data frame with Hugo_Symbol column. Genes in this column should be ordered by importance.
#' @param topGenes a numeric/integer value; How meny top genes from results will be ploted.
#' @param silent a boolean value indicating should silent mutations be counted. By default it is True.
#' @param indels a boolean value indicating should indels be counted. By default it is True.
#' @param allSamples a boolean value indicating if all samples should be plotted.
#' @param Tumor_Sample_Barcode (optional) integer/numeric value indicating column in \code{sample.mutations} which have sample ids for SNVs/Indels. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Hugo_Symbol (optional) integer/numeric value indicating column in \code{sample.mutations} having gene names for reported SNVs/Indels.
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Classification (optional) integer/numeric value indicating column in \code{sample.mutations} which contain classification for SNV (Silent or not). 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param Variant_Type (optional) integer/numeric value indicating column in \code{sample.mutations} which contains indormation if mutations is SNV or InDel .
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param color a numeric or character value indicating column to color the samples.
#' @param groupMax a numeric value which will be taken as maximum value when grouping color column for samples and genes pairs.
#' @param order a character value which can be or \code{'frequency'} or \code{'CCF_sum'}.
#' @return ggplot2 object 
#' @examples 
#' \donttest{
#' #get CCF column
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # get background
#' bcgr.prob <- bcgr.combine(sample.genes.mutect)
#' # get ranking
#' df1 <- bayes.risk(sample.genes.mutect, bcgr.prob, prior.sick = 0.00007) 
#' df2 <- bayes.driver(sample.genes.mutect, bcgr.prob, prior.driver = 0.001) 
#' df.final <- combine.ranking(list(df1, df2), min.mut = 2 )
#' # plot boxplot of CCF for top genes
#' plotStaircase(sample.mutations=sample.genes.mutect, result.df=df.final)
#' }
#' @export
plotStaircase <- function(sample.mutations, result.df, topGenes=40, silent=FALSE, indels=TRUE, allSamples= FALSE,
                       Tumor_Sample_Barcode=NULL, Hugo_Symbol=NULL, Variant_Classification=NULL,  Variant_Type = NULL, 
                       color='CCF' , groupMax=1, order='frequency') {

    if (is.atomic(sample.mutations)) {
        sample.mutations <- data.frame(x = sample.mutations)
    }
    
    if (!is.null(Tumor_Sample_Barcode)){
        sample.mutations <- assign.columns(sample.mutations, Tumor_Sample_Barcode, "Tumor_Sample_Barcode")
    }
    
    if (!is.null(Hugo_Symbol)){
        sample.mutations <- assign.columns(sample.mutations, Hugo_Symbol, "Hugo_Symbol")
    }
    
    
    if (!is.null(Variant_Classification)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Classification, "Variant_Classification")
    }
    

    if (!is.null(Variant_Type)){
        sample.mutations <- assign.columns(sample.mutations, Variant_Type, "Variant_Type")
    }
    
    ##########
    # make it not sensitive to lower/upper case in column names
    original.col.names <- colnames(sample.mutations)
    num.col <- ncol(sample.mutations)
    colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
    
    # take only exonic
    suppressWarnings( sample.mutations <-  exonic.only(sample.mutations) )
    
    if (!is.null(color) & is.character(color)){
        color <- tolower(color)
    }
    
    if(!'Hugo_Symbol' %in% colnames(result.df)){
        stop('There must be Hugo_Symbol column in result.df data frame!')
    }
    
    if(!is.factor(sample.mutations$tumor_sample_barcode)){
        sample.mutations$tumor_sample_barcode <- factor(sample.mutations$tumor_sample_barcode, levels=unique(sample.mutations$tumor_sample_barcode))
    }
    
    plotGenes <- as.character(result.df$Hugo_Symbol[1:min(topGenes,nrow(result.df))])
    sample.mutations <- sample.mutations[sample.mutations$hugo_symbol %in% plotGenes, ]
    
    if(!silent){
        sample.mutations <- sample.mutations[sample.mutations$variant_classification != 'Silent',]
    }
    if(!indels){
        sample.mutations <- sample.mutations[! sample.mutations$variant_type %in% c('DEL', 'INS' ),]
    }
    
    if (order == 'CCF_sum'){
        ccf.sums <- aggregate(ccf ~ hugo_symbol, data=sample.mutations, sum)
        ccf.sums <- ccf.sums[order(ccf.sums$ccf, decreasing=F),]
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=ccf.sums$hugo_symbol )
        
    } else {
        sample.mutations$hugo_symbol <- factor(sample.mutations$hugo_symbol, levels=names( sort(table(sample.mutations$hugo_symbol), decreasing=F )) )
    }
    
    if((is.numeric(color) | is.integer(color)) & color <= ncol(sample.mutations) ){
        col.var <- colnames(sample.mutations)[color]  
    } else if (is.character(color) & color %in% colnames(sample.mutations)){
        col.var <- color 
    } else {
        stop('color need to be or numeric indicator of column in sample.mutations data frame or charater with column name (usually CCF column).')
    }
    
    formula<- paste(col.var,' ~ tumor_sample_barcode + hugo_symbol',  sep='', collapse='')
    
    sample.mutations <- aggregate(as.formula(formula), data = sample.mutations, function (x) min(sum(x),groupMax))   
    
    sample.mutations <- sample.mutations[order( sample.mutations$hugo_symbol, sample.mutations[,col.var], decreasing=T),]
    
    if (allSamples) {
        order.samples <- unique(c(as.character(sample.mutations$tumor_sample_barcode), levels(sample.mutations$tumor_sample_barcode)))
        sample.mutations$tumor_sample_barcode <- factor( sample.mutations$tumor_sample_barcode, levels= order.samples)
    } else {
        sample.mutations$tumor_sample_barcode <- factor( sample.mutations$tumor_sample_barcode, levels=unique(as.character(sample.mutations$tumor_sample_barcode)))
    } 
    
    p <- ggplot(sample.mutations, aes_string(x='tumor_sample_barcode', y='hugo_symbol', fill=col.var))
    
    p <- p  + geom_tile() + scale_fill_gradient(low="white", high="red") +
        theme_bw() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4)) +
        scale_x_discrete(name='Samples',drop=FALSE) + ylab('Genes')
    
  
    
    return(p)
    
}
    
    