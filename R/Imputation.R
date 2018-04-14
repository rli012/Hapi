##' @title Imputation of missing genotypes in the framework
##' @description Imputation of missing genotypes in the framework
##' @param gmt a dataframe of genotype data of gamete cells in the framework
##' @param nSPT a numeric value of the minumum number of supports 
##' for an imputation
##' @param allowNA a numeric value of the maximum number of gametes with 
##' NA at a locus
##' @return a dataframe of imputed genotypes in the framework
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' 
##' gmtFrame <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' 
##' idx1 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' idx2 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' idx3 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' 
##' gmtFrame[idx1,1] <- NA
##' gmtFrame[idx2,2] <- NA
##' gmtFrame[idx3,3] <- NA
##' imputedFrame <- hapiImupte(gmtFrame, nSPT=2, allowNA=0)
hapiImupte <- function(gmt, nSPT=2, allowNA=0) {
    
    total <- nrow(gmt)
    
    minNA2 <- 10000000
    newNA2 <- sum(apply(gmt, 1, function(y) sum(is.na(y)))>=1)
    
    while (newNA2 < minNA2) {
        minNA2 <- newNA2
        
        gmt <- imputationFun1(gmt, nSPT=2)
        gmt <- imputationFun2(gmt, nSPT=2)
        
        newNA2 <- sum(apply(gmt, 1, function(y) sum(is.na(y)))>=1)
    }
    
    #gmt <- imputationFun3(gmt, nSPT=2)
    
    ### filter out failed imputation
    
    filter <- which(apply(gmt, 1, function(y) sum(is.na(y)))>allowNA)
    
    message (paste('Number of hetSNPs after imputation: ', 
        total-length(filter), sep=''))
    
    if (length(filter)>0) {
        gmt <- gmt[-filter,]
        return (gmt)
        
    } else {
        return(gmt)
    }
}




### Determine the allele of imputation (only one) ###
imptDetermineFun <- function(v,nSPT=2) {
    if (sum(!is.na(v))==0) {
        return (7)
    }
    
    v <- v[which(!is.na(v))]
    
    if (length(v) >= nSPT & length(unique(v)) == 1) {
        return (unique(v))
    } else {
        return (7)
    }
}



####################################################################
####################### IMPUTATION

naAdjacentFun <- function(gmt) {
    nonNAPos <- which(!is.na(gmt))
    nonNAPosDiff <- diff(nonNAPos)
    nonNAPosLeft <- nonNAPos[which(nonNAPosDiff != 1)]
    nonNAPosRight <- nonNAPos[which(nonNAPosDiff != 1)+1]
    nonNAMatrix <- cbind(nonNAPosLeft, nonNAPosRight)
    
    return (nonNAMatrix)
}


# two adjacent SNPs both are identical or opposite to those in other pollens #
# supported by at least 2 pollens ###
imputationFun1 <- function(pollens, nSPT=2) {
    
    snpName <- rownames(pollens)
    polName <- colnames(pollens)
    
    polNum <- ncol(pollens)
    snpNumTotal <- nrow(pollens)
    
    impute.final <- matrix(rep(NA,snpNumTotal*polNum), snpNumTotal, polNum)
    dim(impute.final)
    
    minNA <- 10000000
    newNA <- sum(apply(pollens, 1, function(y) sum(is.na(y)))>=1)
    
    pols <- 1:polNum
    while (newNA < minNA) {
        minNA <- newNA
        
        for (k in pols) {
            pol.k <- pollens[,k]
            nonNAMatrix <- naAdjacentFun(pol.k)
            
            impute.final[,k] <- pol.k
            
            if (nrow(nonNAMatrix) == 0) {
                next
                
            } else {
                for (i in 1:nrow(nonNAMatrix)) {
                    l <- nonNAMatrix[i,1]
                    r <- nonNAMatrix[i,2]
                    
                    polPart <- pollens[l:r,]
                    ll <- 1
                    rr <- nrow(polPart)
                    
                    snpNum <- rr-ll+1
                    
                    polPart.k <- polPart[,k]
                    impute.k <- matrix(rep(NA, snpNum*polNum), snpNum, polNum)
                    impute.k[,k] <- polPart[,k]
                    
                    for (j in pols[pols!=k]) {
                        polPart.j <- polPart[,j]
                        if (is.na(polPart.j[ll]) | is.na(polPart.j[rr])) {
                            next
                        } else if (polPart.k[ll] == polPart.j[ll] & 
                            polPart.k[rr] == polPart.j[rr]) {
                            impute.k[ll:rr,j] <- polPart.j[ll:rr]
                        } else if (polPart.k[ll] == flipFun(polPart.j[ll]) & 
                            polPart.k[rr] == flipFun(polPart.j[rr])) {
                            impute.k[ll:rr,j] <- flipFun(polPart.j[ll:rr])
                        } else {
                            next
                        }
                    }
                    #impute.final[l:r,k] <- 
                    #apply(impute.k, 1, function(x) mean(x, na.rm=T))
                    
                    ### only consider the missing data
                    impute.k <- impute.k[-c(1,nrow(impute.k)),] 
                    
                    
                    impute.k <- matrix(impute.k, ncol=polNum)
                    
                    if (is.matrix(impute.k)) {
                        impute.final[(l+1):(r-1),k] <- 
                            apply(impute.k, 1, function(x) 
                                imptDetermineFun(x,nSPT=nSPT))
                        impute.final[(l+1):(r-1),k][impute.final[(l+1):(r-1),k]==7] <- NA
                    } else { ### can be deleted
                        impute.final[(l+1):(r-1),k] <- 
                            imptDetermineFun(impute.k)
                        impute.final[(l+1):(r-1),k][impute.final[(l+1):(r-1),k]==7] <- NA
                    }
                    
                    
                }
            }
        }
        pollens <- impute.final
        newNA <- sum(apply(pollens, 1, function(y) sum(is.na(y)))>=1)
    }
    
    rownames(pollens) <- snpName
    colnames(pollens) <- polName
    return (pollens)
}


##########################################################
multipleNaAdjacentFun <- function(gmt) {
    nonNAPos <- which(apply(gmt,1,function(x) sum(is.na(x))) == 0)
    nonNAPosDiff <- diff(nonNAPos)
    nonNAPosLeft <- nonNAPos[which(nonNAPosDiff != 1)]
    nonNAPosRight <- nonNAPos[which(nonNAPosDiff != 1)+1]
    nonNAMatrix <- cbind(nonNAPosLeft, nonNAPosRight)
    
    return (nonNAMatrix)
}


### solve missing genotypes in a missing region where ###
### pairwise imputation doesn't work ###
### supported by at least 2 pollens ###

imputationFun2 <- function(pollens, nSPT=2) {
    
    snpName <- rownames(pollens)
    polName <- colnames(pollens)
    
    polNum <- ncol(pollens)
    pols <- 1:polNum
    
    nonNAMatrix <- multipleNaAdjacentFun(pollens)
    
    if (nrow(nonNAMatrix) == 0) {
        return (pollens)
        
    } else {
        for (i in 1:nrow(nonNAMatrix)) {
            l <- nonNAMatrix[i,1]
            r <- nonNAMatrix[i,2]
            
            polPart <- pollens[l:r,]
            ll <- 1
            rr <- nrow(polPart)
            
            if (sum(is.na(polPart[ll,])) + sum(is.na(polPart[rr,])) != 0) {
                next
            } else {
                snpNum <- rr-ll+1
                impute.final <- matrix(rep(NA,snpNum*polNum), snpNum, polNum)
                for (k in pols) {
                    polPart.k <- polPart[,k]
                    impute.k <- matrix(rep(NA, snpNum*polNum), snpNum, polNum)
                    impute.k[,k] <- polPart[,k]
                    for (j in pols[pols!=k]) {
                        polPart.j <- polPart[,j]
                        if (polPart.k[ll] == polPart.j[ll] & 
                            polPart.k[rr] == polPart.j[rr]) {
                            impute.k[ll:rr,j] <- polPart.j[ll:rr]
                        } else if (polPart.k[ll] == flipFun(polPart.j[ll]) 
                            & polPart.k[rr] == flipFun(polPart.j[rr])) {
                            impute.k[ll:rr,j] <- flipFun(polPart.j[ll:rr])
                        } else {
                            next
                        }
                    }
                    #impute.final[,k] <- apply(impute.k, 1, 
                    #function(x) mean(x, na.rm=T))
                    
                    impute.final[,k] <- apply(impute.k, 1, function(x) 
                        imptDetermineFun(x,nSPT=nSPT))
                    impute.final[,k][impute.final[,k]==7] <- NA
                    nonNA.k <- which(!is.na(polPart.k))
                    impute.final[nonNA.k,k] <- polPart.k[nonNA.k]
                }
                pollens[l:r,] <- impute.final
            }
        }
        rownames(pollens) <- snpName
        colnames(pollens) <- polName
        return (pollens)
    }
}


### impute missing genotypes in the crossover regions ###
### supported by at least 2 pollens (or more????) ###
imputationFun3 <- function(pollens, nSPT=2) {
    
    snpName <- rownames(pollens)
    polName <- colnames(pollens)
    
    polNum <- ncol(pollens)
    pols <- 1:polNum
    
    nonNAMatrix <- multipleNaAdjacentFun(pollens)
    
    if (nrow(nonNAMatrix) == 0) {
        return (pollens)
        
    } else {
        for (i in 1:nrow(nonNAMatrix)) {
            l <- nonNAMatrix[i,1]
            r <- nonNAMatrix[i,2]
            
            polPart <- pollens[l:r,]
            ll <- 1
            rr <- nrow(polPart)
            
            if (sum(is.na(polPart[ll,])) + sum(is.na(polPart[rr,])) != 0) {
                next
            } else {
                snpNum <- rr-ll+1
                impute.final <- matrix(rep(NA,snpNum*polNum), snpNum, polNum)
                for (k in pols) {
                    polPart.k <- polPart[,k]
                    impute.k <- matrix(rep(NA, snpNum*polNum), snpNum, polNum)
                    impute.k[,k] <- polPart[,k]
                    for (j in pols[pols!=k]) {
                        polPart.j <- polPart[,j]
                        
                        ### may not be needed ?????????
                        if (polPart.k[ll] == polPart.j[ll] & 
                            polPart.k[rr] == polPart.j[rr]) {
                            impute.k[ll:rr,j] <- polPart.j[ll:rr]
                        } else if (polPart.k[ll] == flipFun(polPart.j[ll]) & 
                            polPart.k[rr] == flipFun(polPart.j[rr])) {
                            impute.k[ll:rr,j] <- flipFun(polPart.j[ll:rr])
                            ###  
                            
                        } else if (polPart.k[ll] == polPart.j[ll] & 
                            polPart.k[rr] == flipFun(polPart.j[rr])) {
                            mid <- as.integer((ll+rr)/2)
                            impute.k[ll:mid,j] <- polPart.j[ll:mid]
                            impute.k[(mid+1):rr,j] <- 
                                flipFun(polPart.j[(mid+1):rr])
                        } else if (polPart.k[ll] == flipFun(polPart.j[ll]) & 
                            polPart.k[rr] == polPart.j[rr]) {
                            mid <- as.integer((ll+rr)/2)
                            impute.k[ll:mid,j] <- flipFun(polPart.j[ll:mid])
                            impute.k[(mid+1):rr,j] <- polPart.j[(mid+1):rr]
                        } else {
                            next
                        }
                    }
                    #impute.final[,k] <- apply(impute.k, 1, 
                    #function(x) mean(x, na.rm=T))
                    impute.final[,k] <- apply(impute.k, 1, function(x) 
                        imptDetermineFun(x, nSPT=nSPT))
                    impute.final[,k][impute.final[,k]==7] <- NA
                    nonNA.k <- which(!is.na(polPart.k))
                    impute.final[nonNA.k,k] <- polPart.k[nonNA.k]
                }
                pollens[l:r,] <- impute.final
            }
        }
        rownames(pollens) <- snpName
        colnames(pollens) <- polName
        return (pollens)
    }
}

