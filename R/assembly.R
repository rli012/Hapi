
############################################################

######################### Assembly #########################

##' @title Consensus haplotype assembly
##' @description Assemble the consensus high-resolution haplotypes
##' @param gmt a dataframe of genotype data of gamete cells
##' @param draftHap a dataframe with draft haplotype information
##' @param keepLowConsistency logical, if low-consistent gamete cells 
##' should be kept
##' @param consistencyThresh a numeric value of the threshold determining 
##' low-consistent gamete cells compared with the draft haplotype. 
##' Default is 0.85
##' @return a dataframe containing phased haplotypes
##' @export
##' @author Ruidong Li
##' @examples 
##' finalDraft <- rep(0,500)
##' names(finalDraft) <- seq_len(500)
##' 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' 
##' gmtDa <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' 
##' idx1 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' idx2 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' idx3 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' 
##' gmtDa[idx1,1] <- NA
##' gmtDa[idx2,2] <- NA
##' gmtDa[idx3,3] <- NA
##' 
##' consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)
hapiAssemble <- function(gmt, draftHap, keepLowConsistency=TRUE, 
    consistencyThresh=0.85) {

    ovlp <- Reduce(intersect, list(rownames(gmt),names(draftHap)))
    hap <- rep(NA, nrow(gmt))
    names(hap) <- rownames(gmt)
    hap[ovlp] <- draftHap[ovlp]
    
    #nSNP <- nrow(gmt)
    #nSample <- ncol(gmt)
    
    #coord <- rownames(gmt)
    #sam <- colnames(gmt)
    
    haps <- apply(gmt,2,function(x) hapInferFun(hap, x, consistencyThresh))
    
    filter <- c()
    for (i in 1:ncol(haps)) {
        iden <- haps[,i] == hap
        ratio.t <- sum(iden==TRUE, na.rm=TRUE)/sum(!is.na(iden))
        
        if (ratio.t < consistencyThresh) {
            filter <- c(filter, i)
            message (paste('Warning: bad consistence between draft haplotype 
                and', colnames(haps)[i], ': ', ratio.t, '\n', sep=''))
        } else {
            message (paste('Good consistence between draft haplotype and', 
                colnames(haps)[i], ': ', ratio.t, '\n', sep=''))
        }
    }
    
    if (keepLowConsistency==FALSE & length(filter)>0) {
        haps <- haps[,-filter]
    }
    
    newHap <- apply(haps, 1, function(v) alleleDetermineFun(v))
    newHap <- data.frame(t(newHap), stringsAsFactors=FALSE)
    colnames(newHap) <- c('hap1','hap2','total','rate')
    newHap$confidence <- ifelse(newHap$total>1 & newHap$rate>=0.6, 'H', 'L')
    
    fixedPos <- which(!is.na(hap))
    
    newHap$hap1[fixedPos] <- hap[fixedPos]
    newHap$hap2[fixedPos] <- flipFun(hap[fixedPos])
    newHap$confidence[fixedPos] <- 'F'
    
    return (newHap)
}


##' @title Assembly of haplotypes in regions at the end of a chromosome
##' @description Assembly of haplotypes in regions at the end of a chromosome
##' @param draftHap a dataframe with draft haplotype information
##' @param gmt a dataframe of genotype data of gamete cells
##' @param consensusHap a dataframe of the consensus haplotype information
##' @param k a numeric value for the number of hetSNPs that will be combined
##' with markers beyond the framework for assembly. Default is 300
##' @return a dataframe containing phased haplotypes
##' @export
##' @author Ruidong Li
##' @examples 
##' finalDraft <- rep(0,500)
##' names(finalDraft) <- seq_len(500)
##' 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' 
##' gmtDa <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' 
##' idx1 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' idx2 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' idx3 <- sort(sample(seq_len(500), 30, replace = FALSE))
##' 
##' gmtDa[idx1,1] <- NA
##' gmtDa[idx2,2] <- NA
##' gmtDa[idx3,3] <- NA
##' 
##' consensusHap <- data.frame(hap1=rep(0,500),hap2=rep(1,500),
##' total=rep(5,500),rate=rep(1,500),
##' confidence=rep('F',500),
##' stringsAsFactors = FALSE)
##' rownames(consensusHap) <- seq_len(500)
##' 
##' consensusHap <- hapiAssembleEnd(gmt = gmtDa, draftHap = finalDraft, 
##' consensusHap = consensusHap, k = 300)
hapiAssembleEnd <- function(gmt, draftHap, consensusHap, k=300) {
    
    headPos <- which(rownames(gmt) == names(draftHap)[1])
    tailPos <- which(rownames(gmt) == names(draftHap)[length(draftHap)])
    
    
    if (headPos==1) {
        consensusHap <- consensusHap
    } else {
        headBlock <- gmt[1:(headPos+k-1),]
        headHap <- headTailHapAssembly(headBlock)
        
        sum11 <- sum(headHap[headPos:(headPos+k-1),]$hap1 == 
            consensusHap[headPos:(headPos+k-1),]$hap1)
        sum22 <- sum(headHap[headPos:(headPos+k-1),]$hap1 != 
            consensusHap[headPos:(headPos+k-1),]$hap1)
        
        if (sum11>=sum22) {
            ratio <- sum11/k
        } else {
            ratio <- sum22/k
        }
        
        if (ratio < 0.9) {
            consensusHap <- consensusHap
        } else {
            if (sum11>=sum22) {
                consensusHap[1:(headPos-1),] <- headHap[1:(headPos-1),]
            } else {
                headHap$hap1 <- flipFun(headHap$hap1)
                headHap$hap2 <- flipFun(headHap$hap2)
                consensusHap[1:(headPos-1),] <- headHap[1:(headPos-1),]
            }
        }
    }
    
    if (tailPos==length(draftHap)) {
        consensusHap <- consensusHap
    } else {
        tailBlock <- gmt[(tailPos-k+1):nrow(gmt),]
        
        tailHap <- headTailHapAssembly(tailBlock)
        
        sum11 <- sum(tailHap[1:k,]$hap1 == 
            consensusHap[(tailPos-k+1):tailPos,]$hap1)
        sum22 <- sum(tailHap[1:k,]$hap1 != 
            consensusHap[(tailPos-k+1):tailPos,]$hap1)
        
        
        if (sum11>=sum22) {
            ratio <- sum11/k
        } else {
            ratio <- sum22/k
        }
        
        if (ratio < 0.9) {
            consensusHap <- consensusHap
        } else {
            if (sum11>=sum22) {
                consensusHap[(tailPos+1):nrow(gmt),] <- 
                    tailHap[(k+1):nrow(tailHap),]
            } else {
                tailHap$hap1 <- flipFun(tailHap$hap1)
                tailHap$hap2 <- flipFun(tailHap$hap2)
                consensusHap[(tailPos+1):nrow(gmt),] <- 
                    tailHap[(k+1):nrow(tailHap),]
            }
        }

    }
    
    return (consensusHap)
}








hapInferFun <- function(hap, gmt, consistencyThresh=0.85) {
    
    window <- defineWindowFun(hap, gmt)
    
    newHap <- rep(NA,length(hap))
    
    for (i in 1:length(window$top)) {
        
        start = window$top[i]
        end = window$bottom[i]
        
        naPos <- which(is.na(hap))
        
        if(i%%2!=0){
            newHap[start:end] <- gmt[start:end]
            
        } else {
            newHap[start:end] <- flipFun(gmt[start:end])
            
        }
    }
    
    names(newHap) <- names(hap)
    
    iden <- newHap == hap
    
    if (sum(is.na(iden))==length(newHap)) {
        return (newHap)
        
    } else {
        ratio.t <- sum(iden==TRUE, na.rm=TRUE)/sum(!is.na(iden))
        ratio.f <- sum(iden==FALSE, na.rm=TRUE)/sum(!is.na(iden))
        
        if (ratio.t >= ratio.f) {
            return (newHap)
        } else if (ratio.f > ratio.t) {
            return (flipFun(newHap))
        }
    }
}



#### return coordinates
switchPosFun.L2 <- function(gmt1,gmt2){
    idComp <- gmt1 == gmt2
    idKnownPos <- which(!is.na(idComp))
    
    idKnown <- idComp[idKnownPos]
    idKnown
    
    v1 <- idKnown[-1]
    v2 <- idKnown[-length(idKnown)]
    vd <- v1-v2
    cvPos <- which(vd != 0)
    cvCoord <- idKnownPos[cvPos]
    
    return (as.numeric(cvCoord))
}

######## the original switchPosFun
switchPosFun.R2 <- function(gmt1,gmt2){
    idComp <- gmt1 == gmt2
    idKnownPos <- which(!is.na(idComp))
    
    idKnown <- idComp[idKnownPos]
    idKnown
    
    v1 <- idKnown[-1]
    v2 <- idKnown[-length(idKnown)]
    vd <- v1-v2
    cvPos <- which(vd != 0)+1
    cvCoord <- idKnownPos[cvPos]
    
    return (as.numeric(cvCoord))
}



defineWindowFun <- function(gmt,hap) {
    top <- c(1,switchPosFun.R2(gmt,hap))
    bottom <- c(switchPosFun.L2(gmt,hap), length(hap))
    
    window <- list(top=top, bottom=bottom)
    
    return (window)
    
}



alleleDetermineFun <- function(v) {
    if (sum(!is.na(v))==0) {
        majorBase <- 7
        total <- 0
        majorRate <- 0
        
    } else {
        v <- v[which(!is.na(v))]
        total <- length(v)
        
        if (sum(v) > total/2) {
            majorBase <- 1
            majorRate <- sum(v)/total
        } else if (sum(v) < total/2) {
            majorBase <- 0
            majorRate <- 1-sum(v)/total
        } else {
            majorBase <- 7
            majorRate <- 0.5
        }
    }
    
    v <- v[which(!is.na(v))]
    total <- length(v)
    
    if (sum(v) > total/2) {
        majorBase <- 1
        majorRate <- sum(v)/total
    } else if (sum(v) < total/2) {
        majorBase <- 0
        majorRate <- 1-sum(v)/total
    } else {
        majorBase <- 7
        majorRate <- 0.5
    }
    
    altBase <- flipFun(majorBase)
    majorRate <- round(majorRate, 2)
    return (c(majorBase, altBase, total, majorRate))
}




####################################################################

hmmCVCountFun <- function(gmt1, gmt2, hmm=NULL) {
    if (is.null(hmm)) {
        hmm = initHMM(States=c("F","M"), Symbols=c("f","m"), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    idComp <- gmt1 == gmt2
    idComp <- as.numeric(idComp[!is.na(idComp)])
    
    if (length(idComp)<=1) {
        return (0)
    }
    
    cv <- realHMMCVFun(idComp, hmm=hmm)
    return (nrow(cv))
}


#########################################

cvCountFun <- function(gmt1, gmt2) {
    
    idComp <- gmt1 == gmt2
    idComp <- as.numeric(idComp[!is.na(idComp)])
    
    if (length(idComp)<=1) {
        return (0)
    }
    
    v1 <- idComp[-1]
    v2 <- idComp[-length(idComp)]
    vd <- v1-v2
    
    cv <- sum(vd != 0)
    return (cv)
}



headTailHapAssembly <- function(block, consistencyThresh=0.95) {
    nCV <- c()
    
    for (i in 1:ncol(block)) {
        ref <- block[,i]
        cv <- apply(block, 2, function(gmt) hmmCVCountFun(ref, gmt))
        nCV <- rbind(nCV, cv)
    }
    
    rownames(nCV) <- colnames(nCV)
    
    cv <- apply(nCV, 2, function(v) sum(v==min(nCV)))
    ref <- names(cv)[cv==max(cv)][1]
    
    refs <- nCV[,ref]==min(nCV)
    message('Number of reference pollens: ', sum(refs), '\n', sep='')
    
    
    hapFrame <- block[,refs]
    hapFrame
    
    
    if (sum(refs)==1) {
        refHap <- hapFrame
    } else {
        comp <- apply(hapFrame, 2, function(gmt) 
            pairwiseRefFun(hapFrame[,ref], gmt))
        
        #keep <- comp[1,] >= consistencyThresh
        
        #hapFrame <- hapFrame[,keep]
        #comp <- comp[,keep]
        
        for (i in 1:ncol(hapFrame)) {
            if (comp[2,i] == -1) {
                hapFrame[,i] <- flipFun(hapFrame[,i])
            }
        }
        
        refHap <- apply(hapFrame, 1, mean, na.rm=TRUE)
        refHap[refHap >0 & refHap<1] <- NA
        
        filter <- apply(hapFrame,1, function(v) sum(!is.na(v))<=1)
        refHap[filter] <- NA
    }
    
    cvList <- apply(block, 2, function(v) hmmCVFun(refHap,v))
    names(cvList) <- colnames(block)
    
    position <- as.integer(rownames(block))
    
    for (i in 1:length(cvList)) {
        cvList[[i]] <- apply(cvList[[i]], 1, function(v) 
            as.integer((position[v[1]]+position[v[2]])/2))
    }
    
    haps <- lapply(1:ncol(block), 
        function(i) pairwiseHapInferFun2(gmt=block[,i], 
            cv= cvList[[i]],
            pos=position))
    
    
    haps <- do.call(cbind, haps)
    haps <- data.frame(haps, stringsAsFactors=FALSE)
    
    rownames(haps) <- rownames(block)
    colnames(haps) <- colnames(block)
    
    comp <- apply(haps, 2, function(gmt) pairwiseRefFun(refHap, gmt))
    
    keep <- comp[1,] >= consistencyThresh
    
    #if (sum(keep) <= 0.5*ncol(block)) {
    #  return (FALSE)
    #}
    
    haps <- haps[,keep]
    comp <- comp[,keep]
    
    for (i in 1:ncol(haps)) {
        if (comp[2,i] == -1) {
            haps[,i] <- flipFun(haps[,i])
        }
    }
    
    hap <- apply(haps, 1, function(v) alleleDetermineFun(v))
    hap <- data.frame(t(hap), stringsAsFactors=FALSE)
    colnames(hap) <- c('hap1','hap2','total','rate')
    
    hap$confidence <- ifelse(hap$total>2 & hap$rate>=0.6, 'H', 'L')
    return (hap)
    
}


#####
realHMMCVFun <- function(genoIdentity, hmm=NULL) {
    if (is.null(hmm)) {
        hmm = initHMM(States=c("F","M"), 
            Symbols=c("f","m"), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    if (length(genoIdentity)<=1) {
        vd <- c(0,0)
        cvPos.L <- which(vd != 0)
        cvPos.R <- cvPos.L+1
        return (data.frame(cvPos.L, cvPos.R, stringsAsFactors = FALSE))
    }
    
    genoSymbol <- ifelse(genoIdentity==0,hmm$Symbols[1],hmm$Symbols[2])
    
    correctGeno <- viterbi(hmm, genoSymbol)
    correctGeno <- ifelse(correctGeno==hmm$States[1], 0, 1)
    
    v1 <- correctGeno[-1]
    v2 <- correctGeno[-length(correctGeno)]
    vd <- v1-v2
    
    cvPos.L <- which(vd != 0)
    #cvPos.L <- trueCVFun(cvPos.L, nSNP = nSNP)
    
    cvPos.R <- cvPos.L+1
    
    return (data.frame(cvPos.L, cvPos.R, stringsAsFactors = FALSE))
    
}


hmmCVFun <- function(gmt1, gmt2, hmm=NULL) {
    if (is.null(hmm)) {
        hmm = initHMM(States=c("F","M"), 
            Symbols=c("f","m"), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    idComp <- gmt1 == gmt2
    idKnownPos <- which(!is.na(idComp))
    
    idComp <- as.numeric(idComp[!is.na(idComp)])
    #genoError <- fastCorrectIdentityFun(idComp, position=NULL, hmm=hmm)
    cv <- realHMMCVFun(idComp, hmm=hmm)
    #if (length(genoError)==0) {
    #  return (NULL)
    #}
    
    cvCoord.L <- idKnownPos[cv[,1]]
    cvCoord.R <- idKnownPos[cv[,2]]
    
    return (data.frame(cvCoord.L, cvCoord.R, stringsAsFactors = FALSE))
}


pairwiseRefFun <- function(ref, hap) {
    
    comp <- ref == hap
    
    same <- sum(comp, na.rm=TRUE)
    total <- sum(! is.na(comp))
    
    if (total==0) {
        return (c(0.5, 0))
    }
    
    ratio <- same/total
    
    if (ratio > 0.5) {
        relationship <- 1
    } else {
        relationship <- -1
        ratio <- 1 - ratio
    }
    
    return (c(ratio, relationship))
    
}



########### 
pairwiseHapInferFun2 <- function(gmt,cv, pos) {
    
    if (length(cv)==0) {
        hap1 <- gmt
    } else {
        pos <- as.integer(pos)
        cv <- as.integer(cv)
        
        hap.s <- gmt
        hap.c <- flipFun(gmt)
        
        lth <- length(gmt)
        
        nCV <- length(cv)
        
        if (nCV == 0) {
            hap1 <- hap.s
            #hap2 <- hap.c
            
        } else {
            
            hap1 <- NULL
            ind <- 1
            start <- 0
            for(i in 1:(nCV+1)){
                if(i<=nCV){
                    hapPos <- which(pos >= start & pos < cv[i])
                    
                    if(ind%%2!=0){
                        hap1 <- c(hap1, hap.s[hapPos])
                        ind <- ind + 1
                        start <- cv[i]
                    }else{
                        hap1 <- c(hap1, hap.c[hapPos])
                        ind <- ind + 1
                        start <- cv[i]
                    }
                }else{
                    hapPos <- which(pos >= start & pos <= pos[length(pos)])
                    if(ind%%2!=0){
                        hap1 <- c(hap1, hap.s[hapPos])
                    }else{
                        hap1 <- c(hap1, hap.c[hapPos])
                    }
                }
            }
            #hap2 <- flipFun(hap1)
        }
    }
    
    #return (data.frame(hap1,hap2))
    return (hap1)
}
