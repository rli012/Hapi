
##' @title Convert genotype coded in A/T/C/G to 0/1
##' @description Convert base (A/T/C/G) coded genotype to numeric (0/1) coded 
##' @param gmt a dataframe of genotype data of gamete cells
##' @param ref a character represents reference allele
##' @param alt a character represents alternative allele
##' @return a dataframe containing converted genotype
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
##' gmtDa <- base2num(gmt=gmt[5:9], ref=ref, alt=alt)
base2num <- function(gmt, ref, alt) {
    newGMT <- apply(gmt, 2, function(x) x != ref)
    newGMT[newGMT==TRUE] <- 1
    newGMT[newGMT==FALSE] <- 0
    return (data.frame(newGMT))
}


##' @title Convert genotype coded in 0/1 to A/T/C/G
##' @description Convert numeric (0/1) coded genotype to base (A/T/C/G) coded 
##' @param hap a dataframe of consensus haplotypes
##' @param ref a character represents reference allele
##' @param alt a character represents alternative allele
##' @return a dataframe containing converted haplotypes
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- sample(c('A','T'),500, replace=TRUE)
##' alt <- sample(c('C','G'),500, replace=TRUE)
##' 
##' consensusHap <- data.frame(hap1=rep(0,500),hap2=rep(1,500),
##'     total=rep(5,500),rate=rep(1,500),
##'     confidence=rep('F',500),
##'     stringsAsFactors = FALSE)
##' rownames(consensusHap) <- seq_len(500)
##' 
##' hap <- num2base(hap=consensusHap, ref=ref, alt=alt)
num2base <- function(hap, ref, alt) {
    
    hap$hap1[hap$hap1==0] <- ref[hap$hap1==0]
    hap$hap1[hap$hap1=='1'] <- alt[hap$hap1=='1']
    hap$hap1[hap$hap1=='7'] <- NA
    
    hap$hap2[hap$hap2==0] <- ref[hap$hap2==0]
    hap$hap2[hap$hap2=='1'] <- alt[hap$hap2=='1']
    hap$hap2[hap$hap2=='7'] <- NA
    
    return (hap)
}



### Flip one allele to the althernative allele ###
flipFun <- function(v){
    v2 <- ifelse(v==7,7,ifelse(v==0, 1, 0))
    return (v2)
}

