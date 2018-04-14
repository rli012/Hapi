
##' @title Phase draft haplotypes by majority voting
##' @description Phase draft haplotypes by majority voting
##' @param gmt a dataframe of imputed genotype data of gamete cells
##' @return a dataframe of inferred draft haplotypes
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' imputedFrame <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref, 
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' draftHap <- hapiPhase(gmt=imputedFrame)
hapiPhase <- function(gmt) {
    n <- ncol(gmt)
    
    lastLink = list(link=0,haplink=n,tlink=n)
    newLink = list(link=0,haplink=n,tlink=n)
    hap_1 = 0
    haplink = n
    tlink = n
    ratio <- 0
    cvlink <- 0
    
    for (i in 2:nrow(gmt)) {
        newLink <- walkingFun(lastLink$link, gmt[i-1,], gmt[i,])
        hap_1 <- c(hap_1, newLink$link)
        haplink <- c(haplink, newLink$haplink)
        tlink <- c(tlink, newLink$tlink)
        ratio <- c(ratio, newLink$ratio)
        cvlink <- c(cvlink, newLink$cvlink)
        
        if (is.na(newLink$link)) {
            lastLink <- lastLink
        } else {
            lastLink <- newLink
        }
    }
    
    hap_1[which(is.na(hap_1))] <- 7
    #hap_2 <- flipFun(hap_1)
    
    #hap <- data.frame(hap_1, hap_2)
    #rownames(hap) <- rownames(gmt)
    
    #names(hap_1) <- rownames(gmt)
    #names(haplink) <- rownames(gmt)
    #names(tlink) <- rownames(gmt)
    
    #return (list(hap_1=hap_1,haplink=haplink,tlink=tlink))
    
    hap <- data.frame(hap=hap_1, haplink, cvlink, tlink, ratio)
    rownames(hap) <- rownames(gmt)
    return (hap)
    
}




###########################################
############# make framework ##############
walkingFun <- function(lastLink, site1, site2){
    
    sumSite <- site1 + site2
    sum00 <- sum(sumSite==0 | sumSite==2,na.rm=TRUE)
    sum01 <- sum(sumSite==1,na.rm=TRUE)
    
    
    if (lastLink == 0) {
        if (sum00 > sum01) {
            link <- 0
            haplink <- sum00
            cvlink <- sum01
            tlink <- sum00+sum01
            ratio <- sum01/sum00
            
            
        } else if (sum00 < sum01) {
            link <- 1
            haplink <- sum01
            cvlink <- sum00
            tlink <- sum00+sum01
            ratio <- sum00/sum01
            
        } else if (sum00 == sum01) {
            link <- NA
            haplink <- sum00
            cvlink <- sum01
            tlink <- sum00+sum01
            ratio <- 1
        }
        
    } else if (lastLink == 1) {
        if (sum00 > sum01) {
            link <- 1
            haplink <- sum00
            cvlink <- sum01
            tlink <- sum00+sum01
            ratio <- sum01/sum00
            
        } else if (sum00 < sum01) {
            link <- 0
            haplink <- sum01
            cvlink <- sum00
            tlink <- sum00+sum01
            ratio <- sum00/sum01
            
        } else if (sum00 == sum01) {
            link <- NA
            haplink <- sum00
            cvlink <- sum01
            tlink <- sum00+sum01
            ratio <- 1
        }
    }
    return(list(link=link,haplink=haplink,cvlink=cvlink,
        tlink=tlink, ratio=ratio))
}





##' @title Filter out hetSNPs in potential complex regions
##' @description Filter out hetSNPs in potential complex regions
##' @param draftHap a dataframe with draft haplotype information
##' @param minDistance a numeric value of the distance between two 
##' genomic positions with cv-links. Default is \code{1000000}
##' @param cvlink a numeric value of number of cvlinks. Default is \code{2}
##' @return a dataframe of regions to be filtered out
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- rep(0,500)
##' alt <- rep(1,500)
##' 
##' imputedFrame <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' 
##' draftHap <- hapiPhase(imputedFrame)
##' cvCluster <- hapiCVCluster(draftHap = draftHap, cvlink=2)
hapiCVCluster <- function(draftHap, minDistance=1000000, cvlink=2) {
    
    cvPos <- as.numeric(rownames(draftHap[draftHap$cvlink>=cvlink,]))
    #posAll <- as.numeric(rownames(draftHap))
    
    left <- c()
    right <- c()
    
    if (length(cvPos) < 2) {
        cvCluster <- data.frame(left=cvPos,right=cvPos)
    } else {
        i <- 1
        left <- c(left, cvPos[i])
        
        while (i < length(cvPos)) {
            cvDiff1 <- cvPos[i+1]-cvPos[i]
            
            if (cvDiff1 < minDistance) {
                i <- i + 1
            } else {
                right <- c(right, cvPos[i])
                left <- c(left, cvPos[i+1])
                i <- i + 1
            }
            
            if (i == length(cvPos)) {
                right <- c(right, cvPos[i])
            }
        }
        cvCluster <- data.frame(left,right)
    }
    
    return (cvCluster)
}


##' @title Maximum Parsimony of Recombination (MPR) for 
##' proofreading of draft haplotypes
##' @description Maximum Parsimony of Recombination (MPR) for 
##' proofreading of draft haplotypes
##' @param draftHap a dataframe with draft haplotype information
##' @param gmtFrame a dataframe of raw genotype data in the framework
##' @param cvlink a numeric value of number of cvlinks. Default is \code{2}
##' @param smallBlock a numeric value determining the size of small blocks 
##' that should be excluded from the draft haplotypes 
##' @return a dataframe of draft haplotypes after proofreading
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
##' 
##' imputedFrame <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##' gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##' stringsAsFactors = FALSE)
##' 
##' draftHap <- hapiPhase(imputedFrame)
##' 
##' finalDraft <- hapiBlockMPR(draftHap, gmtFrame, cvlink=2, smallBlock=100)
hapiBlockMPR <- function(draftHap, gmtFrame, cvlink=2, smallBlock=100) {
    
    blockPoint <- which(draftHap$cvlink >= cvlink)
    
    if (length(blockPoint)==0) {
        message('No region is required for proofreading !')
        hap <- draftHap$hap
        names(hap) <- rownames(draftHap)
        return (hap)
    } else {
        start <- 1
        hapBlock <- list()
        for (i in 1:length(blockPoint)) {
            end <- blockPoint[i]-1
            hapBlock[[i]] <- draftHap[start:end,]
            
            start <- blockPoint[i]
            
            if (i == length(blockPoint)) {
                end <- nrow(draftHap)
                hapBlock[[i+1]] <- draftHap[start:end,]
            }
        }
        
        message(paste('Size of blocks: ', 
            paste(unlist(lapply(hapBlock, nrow)), collapse=','), 
            '\n', sep=''))
        
        
        filter <- which(lapply(hapBlock, nrow) < smallBlock)
        message(paste(length(filter), ' blocks are removed !\n',sep=''))
        
        hapBlock[filter] <- NULL
        
        gmt <- gmtFrame[rownames(draftHap),]
        
        currentHap <- hapBlock[[1]]
        for (i in 2:(length(hapBlock))) {
            currentHap <- MPRFun(gmt, currentHap, hapBlock[[i]])
        }
        
        hap <- currentHap$hap
        names(hap) <- rownames(currentHap)
        
        return (hap)
    }
}


MPRFun <- function(gmt, hapBlock1, hapBlock2, nSNP=100) {
    
    filter <- which(hapBlock1$hap==7)
    if (length(filter) > 0) {
        hapBlock1 <- hapBlock1[-filter,]
    }
    
    
    filter <- which(hapBlock2$hap==7)
    if (length(filter) > 0) {
        hapBlock2 <- hapBlock2[-filter,]
    }
    
    
    if (nrow(hapBlock1)>nSNP) {
        sites <- rownames(hapBlock1)[(nrow(hapBlock1)-nSNP+1):nrow(hapBlock1)]
        raw1 <- gmt[sites,]
        hap1 <- hapBlock1[(nrow(hapBlock1)-nSNP+1):nrow(hapBlock1),]
        
    } else {
        raw1 <- gmt[rownames(hapBlock1),]
        hap1 <- hapBlock1
    }
    
    if (nrow(hapBlock2)>nSNP) {
        sites <- rownames(hapBlock2)[1:nSNP]
        raw2 <- gmt[sites,]
        hap2 <- hapBlock2[1:nSNP,]
    } else {
        raw2 <- gmt[rownames(hapBlock2),]
        hap2 <- hapBlock2
    }
    
    raw <- rbind(raw1, raw2)
    
    hap11 <- rbind(hap1, hap2)
    
    hap2$hap <- flipFun(hap2$hap)
    hap22 <- rbind(hap1, hap2)
    
    sum11 <- sum(apply(raw, 2, function(v) cvCountFun(hap11$hap, v)))
    sum22 <- sum(apply(raw, 2, function(v) cvCountFun(hap22$hap, v)))
    
    if (sum11 > sum22) {
        hapBlock2$hap <- flipFun(hapBlock2$hap)
    }
    
    hapBlock <- rbind(hapBlock1, hapBlock2)
    
    message (paste('Number of crossovers given haplotype 1/2: ', 
        sum11, '/', sum22, sep=''))
    
    return (hapBlock)
}



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