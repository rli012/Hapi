##' @title Filter out hetSNPs with potential genotyping errors
##' @description Filter out hetSNPs with potential genotyping errors
##' @param gmt a dataframe of genotype data of gamete cells
##' @param hmm a list containing probabilities of a HMM. Default is \code{NULL}
##' @return a dataframe of genotype data of gamete cells
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' 
##' gmt <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##'     gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##'     stringsAsFactors = FALSE)
##'     
##' idx <- sort(sample(seq_len(500), 10, replace = FALSE))
##' gmt[idx,1] <- 1
##' 
##' gmtDa <- hapiFilterError(gmt = gmt)
hapiFilterError <- function(gmt, hmm=NULL) {
    gmt <- data.frame(gmt)
    total <- nrow(gmt)
    if (is.null(hmm)) {
        hmm <- initHMM(States=c('F','M'), Symbols=c('f','m'), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    nSNP <- 10000000
    while (nrow(gmt) < nSNP) {
        nSNP <- nrow(gmt)
        for (i in 1:ncol(gmt)) {
            genoError <- lapply(gmt[,-i], function(x) 
                filterErrorFun(gmt[,i],x,hmm=hmm))
            filter <- sort(unique(unlist(genoError)))
            
            #message (paste(length(filter), 'hetSNPs 
            #with potential genotyping errors are filtered out !',
            #sep=' '))
            
            if (length(filter) == 0) {
                gmt <- gmt
            } else {
                gmt <- gmt[-filter,]
            }
        }
    }
    
    final <- nrow(gmt)
    
    message (paste(total-final, 
        'hetSNPs with potential genotyping errors are filtered out !',
        sep=' '))
    
    return (gmt)
}



##' @title Selection of hetSNPs to form a framework
##' @description Selection of hetSNPs to form a framework
##' @param gmt a dataframe of genotype data of gamete cells
##' @param n a numeric value of the minumum number of gametes with 
##' observed genotypes at a locus
##' @return a dataframe of genotype data of gamete cells
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' 
##' gmt <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' 
##' idx <- sort(sample(seq_len(500), 10, replace = FALSE))
##' 
##' gmt[idx,1] <- NA
##' gmt[idx,2] <- NA
##' gmt[idx,3] <- NA
##' 
##' gmtFrame <- hapiFrameSelection(gmt = gmt, n = 3)
hapiFrameSelection <- function(gmt, n=3) {
    idx <- which(apply(gmt, 1, function(y) sum(!is.na(y)))>=n)
    message (paste('Number of hetSNPs in the framework: ', 
        length(idx), sep=''))
    gmt <- gmt[idx,]
    return (gmt)
}







####
filterErrorFun <- function(gmt1, gmt2, hmm=NULL) {
    if (is.null(hmm)) {
        hmm <- initHMM(States=c('F','M'), Symbols=c('f','m'), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    idComp <- gmt1 == gmt2
    idKnownPos <- which(!is.na(idComp))
    
    if (length(idKnownPos)<=1) {
        genoError <- fastCorrectIdentityFun(c(1,1), position=NULL, hmm=hmm)
        return (as.numeric(idKnownPos[genoError]))
    }
    
    idComp <- as.numeric(idComp[!is.na(idComp)])
    genoError <- fastCorrectIdentityFun(idComp, position=NULL, hmm=hmm)
    
    #if (length(genoError)==0) {
    #  return (NULL)
    #}
    return (as.numeric(idKnownPos[genoError]))
}


###############


#####
fastCorrectIdentityFun <- function(genoIdentity, position=NULL, hmm=NULL) {
    if (is.null(hmm)) {
        hmm <- initHMM(States=c('F','M'), Symbols=c('f','m'), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    genoSymbol <- ifelse(genoIdentity==0,hmm$Symbols[1],hmm$Symbols[2])
    
    correctGeno <- viterbi(hmm, genoSymbol)
    correctGeno <- ifelse(correctGeno==hmm$States[1], 0, 1)
    
    genoError <- correctGeno-genoIdentity
    genoError <- which(genoError != 0)
    
    return (genoError)
}


