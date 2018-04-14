##' @title Automatic inference of haplotypes
##' @description Automatic inference of haplotypes
##' @param gmt a dataframe of genotype data of gamete cells
##' @param code a character indicating the code style of genotype data. 
##' One of \code{'atcg'} and \code{'01'}. Default is \code{'atcg'}
##' @return a dataframe of inferred consensus haplotypes
##' @importFrom HMM initHMM
##' @importFrom HMM viterbi
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- sample(c('A','T'),500, replace=TRUE)
##' alt <- sample(c('C','G'),500, replace=TRUE)
##' 
##' gmt <- data.frame(chr=rep(1,500), pos=seq_len(500),
##'     ref=ref, alt=alt, gmt1=ref, gmt2=alt, gmt3=ref,
##'     gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##'     stringsAsFactors = FALSE)
##'     
##' hapOutput <- hapiAutoPhase(gmt=gmt, code='atcg')
hapiAutoPhase <- function(gmt, code='atcg') {
    rownames(gmt) <- gmt$pos
    hetDa <- gmt[,1:4]
    
    if (code=='atcg') {
        ref <- hetDa$ref
        alt <- hetDa$alt
        gmtDa <- gmt[,-(1:4)]
        gmtDa <- base2num(gmt=gmtDa, ref=ref, alt=alt)
    } else {
        gmtDa <- gmt[,-(1:4)]
    }
    
    hmm = initHMM(States=c("S","d"), Symbols=c("S","d"), 
        transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
        emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
        startProbs = c(0.5,0.5))
    
    gmtDa <- hapiFilterError(gmt=gmtDa, hmm = hmm)
    
    gmtFrame <- hapiFrameSelection(gmtDa, 3) ###
    
    imputedFrame <- hapiImupte(gmtFrame, nSPT=2, allowNA=0)
    
    ### NEW framework construction method ###
    draftHap <- hapiPhase(imputedFrame) ###
    cvCluster <- hapiCVCluster(draftHap, cvlink=2)
    
    filter <- c()
    for (i in 1:nrow(cvCluster)) {
        filter <- c(filter, which (rownames(draftHap) >= cvCluster$left[i] 
            & rownames(draftHap) <= cvCluster$right[i]))
    }
    
    if (length(filter) > 0) {
        imputedFrame <- imputedFrame[-filter, ]
        draftHap <- hapiPhase(imputedFrame)
    } 
    
    finalDraft <- hapiBlockMPR(draftHap, gmtFrame, cvlink = 2)
    
    consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)
    consensusHap <- hapiAssembleEnd(gmt = gmtDa, draftHap = finalDraft, 
        consensusHap = consensusHap, k = 300)
    
    snp <- which(rownames(hetDa) %in% rownames(consensusHap))
    
    if (code=='atcg') {
        ref <- hetDa$ref[snp]
        alt <- hetDa$alt[snp]
        
        consensusHap <- num2base(hap=consensusHap, ref=ref, alt=alt)
        hapOutput <- data.frame(gmt[snp,], consensusHap)
    } else {
        hapOutput <- data.frame(gmt[snp,], consensusHap)
    }
    
    return(hapOutput)
}