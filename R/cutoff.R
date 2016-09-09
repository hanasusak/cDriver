######################################################################################
# CRG 
# Hana SUSAK
#------------------------------------------------------------------------------------
# cut-off suggestion functions
######################################################################################

#' @title Calculation of suggested cut off for bayesian risk model
#' @description
#'   \code{cutoff.risk} function runs Bayesian risk inference model n times, but with randomly generated gene names 
#'   (probablity of gene beeing mutated is taken from background model)
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
#' @param n a integer number indicating how many random genes mutations (by background probablity) tests will be done. Default is 100.
#' @param fdr expected false discover rate. Value can be between 0 and 1, while closer to 0 less false discoveries will be allowed. 
#'      Default value is 0.1 (10\% of ranked genes before suggested cut off are expected to be false postives).
#' @param simulation.quantile represent numeric value between 0 and 1 that will take for each ranking that qunantile from n permutations. 
#'      Default value is 0.5 (median).
#' @param genes vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV/Indel. in given cohort). Default NULL.
#' @param prior.sick a numeric value representing incidence of tumor in population. Set by default to 0.0045
#' @param plot.save a boolean variable to indicate if plot should be saved
#' @param permutationResults.save a boolean variable to indicate if n permutations results should be saved
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
#' @param mode a charechter value indicationg how to solve when in one gene sample pair there are multiple mutations. Options are SUM, MAX and ADVANCE
#' @param epsilon a numeric value. If mode is ADVANCE, epsilone value will be threshold for CCF difference to decide if they are in same or different clone. 
#' @return a integer value, where suggested cut off for ranking is.
#' @seealso \code{\link{CCF}}, \code{\link{bcgr}}, \code{\link{bcgr.lawrence}}, \code{\link{bcgr.combine}}   and \code{\link{bayes.risk}} 
#' @examples 
#' \donttest{
#' # first calculate CCF
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # then somatic background probability
#' bcgr.prob <- bcgr.combine(sample.genes.mutect)
#' # bayes risk model suggested cut off
#' suggested.cut.off <- cutoff.risk(sample.genes.mutect,  bcgr.prob, prior.sick = 0.00007) 
#' print(suggested.cut.off)  
#' }
#' @export
cutoff.risk <- function(sample.mutations,  bcgr.prob, n=100, fdr=0.1, simulation.quantile=0.5,  genes=NULL, 
                        prior.sick = 0.0045, plot.save=FALSE, permutationResults.save=FALSE, 
                        Variant_Classification=NULL, Hugo_Symbol=NULL,
                        Tumor_Sample_Barcode=NULL, CCF=NULL,Damage_score=NULL, mode='MAX', epsilon=0.05){

    suppressWarnings(res1 <- bayes.risk(sample.mutations,  bcgr.prob=bcgr.prob, genes=genes,  prior.sick=prior.sick,
                             Variant_Classification=Variant_Classification,Hugo_Symbol=Hugo_Symbol,
                              Tumor_Sample_Barcode=Variant_Classification, CCF=CCF, Damage_score=Damage_score,
                             mode=mode, epsilon=epsilon ) )
    cols.pos <- match(c('Hugo_Symbol','postProb'),colnames(res1))
    true.model <- res1[,cols.pos]
    true.model <- true.model[order(true.model$postProb, decreasing=T),]
    true.model$rank <- 1:nrow(true.model)
    true.model$type <- 'TRUE'
    true.model$perm_number <- 'cDriver_run'
    rm(res1)
    
    cancer.maf.rand <- sample.mutations
    n.nonsilent <- nrow(cancer.maf.rand[cancer.maf.rand$Variant_Classification !=  'Silent', ])
    list_b1 <- list()
    for(i in 1:n){
        message(paste('Running permutation',i))
        cancer.maf.rand$Hugo_Symbol <- factor(cancer.maf.rand$Hugo_Symbol, levels=names(bcgr.prob))
        cancer.maf.rand[cancer.maf.rand$Variant_Classification != 'Silent', 'Hugo_Symbol'] <- sample(names(bcgr.prob),replace=T, prob=bcgr.prob, size=n.nonsilent)
        suppressWarnings(res1.rand <- bayes.risk(cancer.maf.rand,  bcgr.prob=bcgr.prob, genes=genes,  prior.sick=prior.sick,
                                Variant_Classification=Variant_Classification,Hugo_Symbol=Hugo_Symbol,
                                Tumor_Sample_Barcode=Variant_Classification, CCF=CCF, Damage_score=Damage_score,
                                mode=mode, epsilon=epsilon )   )             
        df.postProb.rand <- res1.rand[,cols.pos]
        colnames(df.postProb.rand) <- c('Hugo_Symbol',paste0('bayes_risk_',i))
        list_b1[[i]] <- df.postProb.rand
        rm(res1.rand)             
    }
    rm(i)
    rm(n.nonsilent)
    df.rand.bayes_risk <- Reduce(function(...) merge(...,by='Hugo_Symbol', all=T), list_b1)
    if(permutationResults.save){
        write.table(df.rand.bayes_risk, file=paste0('Bayes_risk_',n,'_permutations_postProb.txt'),
                    sep='\t', row.names = F, col.names = T, quote = F)
    }
    
    df.random.melt <- melt(df.rand.bayes_risk)
    colnames(df.random.melt) <- c('Hugo_Symbol','perm_number','postProb')
    df.random.melt <- df.random.melt[order(df.random.melt$perm_number, df.random.melt$postProb,  decreasing=T,na.last=T),]
    df.random.melt$rank <- rep(1:nrow(df.rand.bayes_risk),n)
    df.random.melt$type <- 'RANDOM'
    df.random.melt <- df.random.melt[,c('Hugo_Symbol','postProb','rank','type','perm_number')]
    
    #  random postProb median per ranking, from 100 lines
    df.random.final<- data.frame(rank=1:nrow(df.rand.bayes_risk))
    df.random.final$postProb <- NA
    for (i in 1:nrow(df.rand.bayes_risk)){
        df.random.final[df.random.final$rank == i, 'postProb'] <- quantile(df.random.melt[df.random.melt$rank == i, 'postProb'],simulation.quantile, na.rm=T)
    }
    df.random.final$type <- paste('Random ',simulation.quantile*100,'% Quantile', sep='')
    df.random.final$Hugo_Symbol <- NA
    df.random.final$perm_number <- 'median_perms'
    df.random.final<- df.random.final[,c('Hugo_Symbol','postProb','rank','type','perm_number')]

       
    df.rank.combine <- merge(x=true.model, y=df.random.final, by='rank', all.t=T)
    df.rank.combine$fdr.value <- apply(df.rank.combine,1, function (x) sum(df.rank.combine$postProb.y >= as.numeric(x["postProb.x"] ),na.rm = T )/ sum(df.rank.combine$postProb.x >=  as.numeric(x["postProb.x"]),na.rm = T ) )
    cutt.off <- max(which(df.rank.combine$fdr.value < fdr))
    print(paste0('Suggested cut-off for expected ',fdr*100,'% FDR based on ',n,' random simulations, is at gene rank: ',cutt.off))
    
    post.prob.cut <- df.rank.combine[cutt.off,'postProb.x']
    
    df.all <- rbind.data.frame(true.model,df.random.melt,df.random.final)
    if(plot.save){
        p1 <- ggplot(df.all, aes_string('rank', 'postProb', group='perm_number', colour = 'type', alpha='type', size='type')) +
            geom_line() + geom_hline(yintercept=post.prob.cut,color='red') +
            scale_size_manual(values = c(1,1,1.5)) +   
            scale_alpha_manual(values = c(0.5,1,1)) +
            scale_color_manual(values = c("darkgray", "black","blue")) + 
            annotate("text", x = max(cutt.off*1.5,200), y = post.prob.cut*1.05, label = paste0(cutt.off, ' ranking'), vjust=0, hjust=0, size=6) +
            ggtitle("Bayesian Hazard (Risk) model with permutations") + theme_bw()
        
        name <- paste0('Bayes-Risk_Cut-off-suggestion_',n,'-simulations_',fdr*100,'-fdr.pdf')
        pdf(name, width = 10, height = 8, useDingbats = FALSE )
        suppressWarnings(print(p1))
        dev.off()
        
        p2 <- ggplot(df.all, aes_string('rank', 'postProb', group='perm_number', colour = 'type', alpha='type', size='type')) +
            geom_segment(x = cutt.off, y = 0, xend = cutt.off, yend = post.prob.cut, colour = "red", size=1, show_guide=FALSE) +
            geom_line() + geom_hline(yintercept=post.prob.cut,color='red') +
            scale_size_manual(values = c(1,1,1.5)) +   
            scale_alpha_manual(values = c(0.5,1,1)) +
            scale_color_manual(values = c("darkgray", "black","blue")) + 
            annotate("text", x = cutt.off*1.05, y = post.prob.cut*1.05, label = paste0(cutt.off, ' ranking'), vjust=0, hjust=0, size=6) +
            xlim(0,cutt.off*3)  + ggtitle("Bayesian Hazard (Risk) model with permutations") + theme_bw()
        name <-paste0('ZOOM_Bayes-Risk_Cut-off-suggestion_',n,'-simulations_',fdr*100,'-fdr.pdf')
        pdf(name, width = 10, height = 8, useDingbats = FALSE )
            suppressWarnings(print(p2))
        dev.off()
    }
    return(as.integer(cutt.off))
}


#' @title Calculation of suggested cut off for bayesian risk model
#' @description
#'   \code{cutoff.driver} function runs Bayesian driver inference model n times, but with randomly generated gene names 
#'   (probablity of gene beeing mutated is taken from background model)
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
#' @param n a integer number indicating how many random genes mutations (by background probablity) tests will be done. Default is 100.
#' @param fdr expected false discover rate. Value can be between 0 and 1, while closer to 0 less false discoveries will be allowed. 
#'      Default value is 0.1 (10\% of ranked genes before suggested cut off are expected to be false postives).
#' @param simulation.quantile represent numeric value between 0 and 1 that will take for each ranking that qunantile from n permutations. 
#'      Default value is 0.5 (median).
#' @param genes vector of genes which were sequenced. 
#' They should be unique values of Hugo_Symbol column (with possibility of more additional genes which did not have any SNV/Indel. in given cohort). Default NULL.
#' @param prior.driver a numeric value representing prior probability that random gene is dirver. 
#'          Default is set to \code{length(driver.genes)}/20000, as it assumed there is ~20000 protein goding genes.
#' @param gene.mut.driver a numeric value or named vector representing likelihood that gene is mutated if it is knowen to be driver. 
#'          Gene does not need to be mutated if it is driver, as cancers are heterogenious. Default is set to NULL and driver.genes are considered as drivers.
#' @param driver.genes a character vector of genes which are considered as drivers for this cancer. If NULL then used set is \code{driver.genes.concensus}.
#' @param plot.save a boolean variable to indicate if plot should be saved
#' @param permutationResults.save a boolean variable to indicate if n permutations results should be saved
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
#' @param mode a charechter value indicationg how to solve when in one gene sample pair there are multiple mutations. Options are SUM, MAX and ADVANCE
#' @param epsilon a numeric value. If mode is ADVANCE, epsilone value will be threshold for CCF difference to decide if they are in same or different clone. 
#' @return a integer value, where suggested cut off for ranking is.
#' @seealso \code{\link{CCF}}, \code{\link{bcgr}}, \code{\link{bcgr.lawrence}}, \code{\link{bcgr.combine}}   and \code{\link{bayes.driver}} 
#' @examples 
#' \donttest{
#' # first calculate CCF
#' sample.genes.mutect <- CCF(sample.genes.mutect)
#' # then somatic background probability
#' bcgr.prob <- bcgr.combine(sample.genes.mutect)
#' # bayes risk model suggested cut off
#' suggested.cut.off <- cutoff.driver(sample.genes.mutect,  bcgr.prob) 
#' print(suggested.cut.off)  
#' }
#' @export
cutoff.driver <- function(sample.mutations,  bcgr.prob, n=100, fdr=0.1, simulation.quantile=0.5,  genes=NULL, 
                          prior.driver = NULL, gene.mut.driver=NULL, driver.genes=NULL, plot.save=FALSE, 
                          permutationResults.save=FALSE, Variant_Classification=NULL, Hugo_Symbol=NULL,
                        Tumor_Sample_Barcode=NULL, CCF=NULL, Damage_score=NULL, mode='MAX', epsilon=0.05){
    
    suppressWarnings( res1 <- bayes.driver(sample.mutations,  bcgr.prob=bcgr.prob, genes=genes,  prior.driver = prior.driver,
                        gene.mut.driver=gene.mut.driver, driver.genes=driver.genes,
                       Variant_Classification=Variant_Classification,Hugo_Symbol=Hugo_Symbol,
                       Tumor_Sample_Barcode=Variant_Classification, CCF=CCF, Damage_score=Damage_score,
                       mode=mode, epsilon=epsilon ) )
    cols.pos <- match(c('Hugo_Symbol','postProb'),colnames(res1))
    true.model <- res1[,cols.pos]
    true.model <- true.model[order(true.model$postProb, decreasing=T),]
    true.model$rank <- 1:nrow(true.model)
    true.model$type <- 'TRUE'
    true.model$perm_number <- 'cDriver_run'
    rm(res1)
    
    cancer.maf.rand <- sample.mutations
    n.nonsilent <- nrow(cancer.maf.rand[cancer.maf.rand$Variant_Classification !=  'Silent', ])
    list_b1 <- list()
    for(i in 1:n){
        message(paste('Running permutation',i))
        cancer.maf.rand$Hugo_Symbol <- factor(cancer.maf.rand$Hugo_Symbol, levels=names(bcgr.prob))
        cancer.maf.rand[cancer.maf.rand$Variant_Classification != 'Silent', 'Hugo_Symbol'] <- sample(names(bcgr.prob),replace=T, prob=bcgr.prob, size=n.nonsilent)
        suppressWarnings( res1.rand <- bayes.driver(cancer.maf.rand,  bcgr.prob=bcgr.prob, genes=genes, prior.driver = prior.driver, 
                                gene.mut.driver=gene.mut.driver, driver.genes=driver.genes ,
                                Variant_Classification=Variant_Classification,Hugo_Symbol=Hugo_Symbol,
                                Tumor_Sample_Barcode=Variant_Classification, CCF=CCF, Damage_score=Damage_score,
                                mode=mode, epsilon=epsilon ) )                
        df.postProb.rand <- res1.rand[,cols.pos]
        colnames(df.postProb.rand) <- c('Hugo_Symbol',paste0('bayes_driver_',i))
        list_b1[[i]] <- df.postProb.rand
        rm(res1.rand)             
    }
    rm(i)
    rm(n.nonsilent)
    df.rand.bayes_driver <- Reduce(function(...) merge(...,by='Hugo_Symbol', all=T), list_b1)
    if(permutationResults.save){
        write.table(df.rand.bayes_driver, file=paste0('Bayes_driver_',n,'_permutations_postProb.txt'),
                    sep='\t', row.names = F, col.names = T, quote = F)
    }
    
    df.random.melt <- melt(df.rand.bayes_driver)
    colnames(df.random.melt) <- c('Hugo_Symbol','perm_number','postProb')
    df.random.melt <- df.random.melt[order(df.random.melt$perm_number, df.random.melt$postProb,  decreasing=T,na.last=T),]
    df.random.melt$rank <- rep(1:nrow(df.rand.bayes_driver),n)
    df.random.melt$type <- 'RANDOM'
    df.random.melt <- df.random.melt[,c('Hugo_Symbol','postProb','rank','type','perm_number')]
    
    #  random postProb median per ranking, from 100 lines
    df.random.final<- data.frame(rank=1:nrow(df.rand.bayes_driver))
    df.random.final$postProb <- NA
    for (i in 1:nrow(df.rand.bayes_driver)){
        df.random.final[df.random.final$rank == i, 'postProb'] <- quantile(df.random.melt[df.random.melt$rank == i, 'postProb'],simulation.quantile, na.rm=T)
    }
    df.random.final$type <- paste('Random ',simulation.quantile*100,'% Quantile', sep='')
    df.random.final$Hugo_Symbol <- NA
    df.random.final$perm_number <- 'median_perms'
    df.random.final<- df.random.final[,c('Hugo_Symbol','postProb','rank','type','perm_number')]
    
    
    df.rank.combine <- merge(x=true.model, y=df.random.final, by='rank', all.t=T)
    df.rank.combine$fdr.value <- apply(df.rank.combine,1, function (x) sum(df.rank.combine$postProb.y >= as.numeric(x["postProb.x"] ),na.rm = T )/ sum(df.rank.combine$postProb.x >=  as.numeric(x["postProb.x"] ),na.rm = T ) )
    cutt.off <- max(which(df.rank.combine$fdr.value < fdr))
    print(paste0('Suggested cut-off for expected ',fdr*100,'% FDR based on ',n,' random simulations, is at gene rank: ',cutt.off))
    
    post.prob.cut <- df.rank.combine[cutt.off,'postProb.x']
    df.all <- rbind.data.frame(true.model,df.random.melt,df.random.final)
   
    if(plot.save){
        p1 <- ggplot(df.all, aes_string('rank', 'postProb', group='perm_number', colour = 'type', alpha='type', size='type')) +
            geom_line() + geom_hline(yintercept=post.prob.cut,color='red') +
            scale_size_manual(values = c(1,1,1.5)) +   
            scale_alpha_manual(values = c(0.5,1,1)) +
            scale_color_manual(values = c("darkgray", "black","blue")) + 
            annotate("text",  x = max(cutt.off*1.5,200), y = post.prob.cut*1.05, label = paste0(cutt.off, ' ranking'), vjust=0, hjust=0, size=6) +
            ggtitle("Bayesian Driver model with permutations")  + theme_bw()
        name <- paste0('Bayes-Driver_Cut-off-suggestion_',n,'-simulations_',fdr*100,'-fdr.pdf')
        pdf(name, width = 10, height = 8, useDingbats = FALSE )
            suppressWarnings(print(p1))
        dev.off()
        p2 <-ggplot(df.all, aes_string('rank', 'postProb', group='perm_number', colour = 'type', alpha='type', size='type')) +
            geom_segment(x = cutt.off, y = 0, xend = cutt.off, yend = post.prob.cut, colour = "red", size=1, show_guide=FALSE) +
            geom_line() + geom_hline(yintercept=post.prob.cut,color='red') +
            scale_size_manual(values = c(1,1,1.5)) +   
            scale_alpha_manual(values = c(0.5,1,1)) +
            scale_color_manual(values = c("darkgray", "black","blue")) + 
            annotate("text", x = cutt.off*1.05, y = post.prob.cut*1.05, label = paste0(cutt.off, ' ranking'), vjust=0, hjust=0, size=6) +
            xlim(0,cutt.off*3)  + ggtitle("Bayesian Driver model with permutations")  + theme_bw()
        name <- paste0('ZOOM_Bayes-Driver_Cut-off-suggestion_',n,'-simulations_',fdr*100,'-fdr.pdf')
        pdf(name, width = 10, height = 8, useDingbats = FALSE )
        suppressWarnings(print(p2))
        dev.off()
    }
    return(as.integer(cutt.off))
}


