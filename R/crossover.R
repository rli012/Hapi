##' @title Indentify crossovers in gamete cells
##' @description Indentify crossovers in gamete cells
##' @param hap a dataframe of the two haplotypes
##' @param gmt a dataframe of genotype data of gamete cells
##' @param hmm a list containing probabilities of a HMM. Default is \code{NULL}
##' @return a dataframe containing crossover information in each gamete cell
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- sample(c('A','T'),500, replace=TRUE)
##' alt <- sample(c('C','G'),500, replace=TRUE)
##' 
##' hap <- data.frame(hap1=ref, hap2=alt, stringsAsFactors = FALSE)
##' rownames(hap) <- seq_len(500)
##' 
##' gmt <- data.frame(gmt1=ref, gmt2=alt, gmt3=ref,
##'     gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##'     stringsAsFactors = FALSE)
##'     
##' cvOutput <- hapiIdentifyCV(hap=hap, gmt=gmt)
hapiIdentifyCV <- function(hap, gmt, hmm=NULL) {
    
    if (is.null(hmm)) {
        hmm = initHMM(States=c("F","M"), 
            Symbols=c("f","m"), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    if (is.vector(gmt)) {
        cvPos <- hmmCVFun(hap[,1], gmt, hmm)
        if (nrow(cvPos)==0) {
            return (NULL)
        } else {
            cvCoord <- lapply(cvPos, function(v) rownames(hap)[v])
            cvCoord <- do.call(cbind, cvCoord)
            
            cv <- data.frame(cvCoord, stringsAsFactors = FALSE)
            colnames(cv) <- c('start', 'end')
            
            cv$start <- as.numeric(cv$start)
            cv$end <- as.numeric(cv$end)
            
            cv$pos <- as.integer((cv$start+cv$end)/2)
            cv$res <- cv$end-cv$start-1
            return (cv)
        }
    } else {
        cv <- c()
        cvPos <- apply(gmt, 2, function(v) hmmCVFun(hap[,1], v, hmm))
        
        for (i in 1:length(cvPos)) {
            if (nrow(cvPos[[i]])==0) {
                next
            } else {
                cvCoord <- lapply(cvPos[[i]], function(v) rownames(hap)[v])
                cvCoord <- do.call(cbind, cvCoord)
                cvInfo <- c()
                for (j in 1:nrow(cvCoord)) {
                    cvInfo <- rbind(cvInfo, c(names(cvPos)[i], cvCoord[j,]))
                }
                cv <- rbind(cv, cvInfo)
            }
        }
        
        cv <- data.frame(cv, stringsAsFactors = FALSE)
        colnames(cv) <- c('gmt', 'start', 'end')
        
        cv$start <- as.numeric(cv$start)
        cv$end <- as.numeric(cv$end)
        
        cv$pos <- as.integer((cv$start+cv$end)/2)
        cv$res <- cv$end-cv$start-1
        
        return (cv)
    }
}



##' @title Histogram of crossover resolution
##' @description Histogram of crossover resolution
##' @param cv a dataframe of crossover information
##' @return a histogram
##' @export
##' @author Ruidong Li
##' @examples
##' data(crossover)
##' hapiCVResolution(cv=crossover)
hapiCVResolution <- function(cv) {
    ggplot(cv, aes(res/1000)) +
        geom_histogram(binwidth=10, colour="black", fill="white") +
        scale_y_continuous(sec.axis = sec_axis(~./40, 
            name = "Cumulative percentage (%)"), limits = c(0, 40)) +
        xlab('Resolution of crossovers (kb)') + 
        ylab('Number of crossover events') +
        #geom_step(aes(y=..y..*40),stat="ecdf")
        geom_line(aes(y=..y..*40),stat = "ecdf", color='red') +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
            axis.text = element_text(size=14),
            axis.title = element_text(size=16),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour='black'),
            panel.background = element_blank()) +
        theme(axis.title.y.right = element_text(angle = 90, hjust=0.5))
    
}




##' @title Histogram of crossover distance
##' @description Histogram of crossover distance
##' @param cv a dataframe of crossover information
##' @return a histogram
##' @export
##' @author Ruidong Li
##' @examples
##' data(crossover)
##' hapiCVDistance(cv=crossover)
hapiCVDistance <- function(cv) {
    
    cvDist <- c()
    for (gmt in unique(cv$gmt)) {
        gmtCV <- cv[cv$gmt==gmt,]
        for (chr in unique(gmtCV$chr)) {
            cvNew <- gmtCV[gmtCV$chr==chr,]
            #print (cvNew)
            if (nrow(cvNew)>=2) {
                distance <- diff(cvNew$pos)
                cvDist <- c(cvDist, distance)
            }
        }
    }
    cvDist <- data.frame(num=1:length(cvDist), cvDist)
    
    ### distribution
    ggplot(cvDist, aes(cvDist/1000000)) + 
        geom_histogram(aes(y=..count..), binwidth=10,
            colour="black", fill="white") + #ylim(0,15) +
        #geom_density(aes(y=10*..count..), color='red') +
        geom_line(aes(y=10*..count..), stat="density", color='red') +
        xlab('Distance between neighboring crossovers (Mb)') + 
        ylab('Event counts') +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
            axis.text = element_text(size=14),
            axis.title = element_text(size=16),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
}


##' @title Visualization of crossover map
##' @description Visualization of crossover map
##' @param cv a dataframe of crossover information
##' @param chr a dataframe of chromosome information, including length, 
##' and centrometric regions
##' @param step a numeric value of genomic interval in Mb. 
##' Default is \code{5} 
##' @param gap a dataframe of unassembled regions with the first column is
##' chromosme, the second column is start position, and third column is the 
##' end position of the gap. Default is gap for hg19. 
##' If no gap region is provided, use \code{gap=NULL}
##' @param x.limits a numeric value of limits on x axis
##' @param y.breaks a vector of positions to show labels on y axis. 
##' Default is \code{NULL}
##' @param y.labels a vector of labels on the y axis. Default is \code{NULL}
##' @return a plot of crossover map on all the chromosomes
##' @export
##' @author Ruidong Li
##' @examples
##' data(crossover)
##' \dontrun{hapiCVMap(cv=crossover)}
hapiCVMap <- function(cv, chr=hg19, step=5, gap=gap.hg19, 
    x.limits=6, y.breaks=NULL, y.labels=NULL) {
    
    nChr <- nrow(chr)
    chrSize <- chr$length
    
    cvPlot <- c()
    
    for (j in 1:nChr) {
        nCV <- c()
        x <- seq(2.5,chrSize[j]/1000000,step)
        
        for (i in x) {
            nCV <- c(nCV, sum(cv$pos[cv$chr==j]/1000000 > i-2.5 & 
                cv$pos[cv$chr==j]/1000000 <= i-2.5+step))
        }
        
        cvDa <- data.frame(x=c(0,x,chrSize[j]/1000000),y=c(0,nCV,0), 
            chr=rep(j,length(x)+2))
        cvPlot <- rbind(cvPlot, cvDa)
    }
    
    maxChr <- max(chrSize)
    
    chrLength <- chr$length
    cenStart <- chr$cenStart
    cenEnd <- chr$cenEnd
    
    xMin <- 1:nChr-0.25
    xMax <- 1:nChr+0.25
    yMin <- rep(0,nChr)
    yMax <- chrLength
    
    boxDa <- data.frame(xMin,xMax,yMin,yMax)
    
    #####
    cenXMin <- xMin
    cenXMax <- xMax
    
    cenYMin <- cenStart
    cenYMax <- cenEnd
    
    cenDa <- data.frame(cenXMin,cenXMax,cenYMin,cenYMax)
    
    cenX <- 1:nChr
    cenY <- (cenStart+cenEnd)/2
    
    cenDa <- data.frame(cenX, cenY)
    
    chrPlot <- data.frame(cenYMin = cenYMin/1000000,
        cenYMax = cenYMax/1000000,
        cenXMin=rep(0,nChr), cenXMax=rep(1.2,nChr),
        yMin = yMin/1000000,
        yMax = yMax/1000000,
        xMin=rep(0,nChr), xMax=rep(1.2,nChr),
        chr = 1:nChr)
    
    cvSitePlot <- data.frame(chr=cv$chr, pos=cv$pos/1000000, 
        y=rep(0.7,nrow(cv)))

    cvSitePlot$chr <- factor(cvSitePlot$chr, 
        labels = paste('chr',1:nChr,sep=''))
    cvPlot$chr <- factor(cvPlot$chr, 
        labels = paste('chr',1:nChr,sep=''))
    chrPlot$chr <- factor(chrPlot$chr, 
        labels = paste('chr',1:nChr,sep=''))
    
    if (is.null(y.breaks)) {
        y.breaks = c(0,50,100,150,200)
    }
    if (is.null(y.labels)) {
        y.labels = c(0,'50Mb','100Mb','150Mb','200Mb')
    }
    

    if (is.null(gap)) {
        ggplot()+
            geom_line(data = cvPlot,
                aes(x = x, y = y+1.25), colour = 'blue', 
                stat="identity", size=0.6) + 
            geom_rect(data=chrPlot,aes(xmin=yMin, xmax=yMax, 
                ymin=xMin, ymax=xMax), color='limegreen',
                fill='limegreen',size=0.2) + 
            geom_rect(data=chrPlot, aes(xmin=cenYMin, xmax=cenYMax, 
                ymin=cenXMin, ymax=cenXMax), color='black',
                fill='black',size=0.2) +
            geom_point(data=cvSitePlot, aes(x=pos,y=y), 
                color='red',fill='red',size=2.5, shape=4, stroke =0.8) +
            coord_flip() +  facet_grid(.~chr, scales = "free", space='free') +
            
            #geom_point(data=cenDa, aes(x=cenX[1], y=cenY[1]), 
            #size=3, shape=1, stroke=2) +
            scale_x_reverse(breaks = y.breaks, limits = c(maxChr/1000000,0),
                            labels = y.labels) +
            scale_y_continuous(limits = c(0,x.limits+1.25), 
                breaks=0:x.limits+1.25, labels=0:x.limits) +
            #scale_color_manual(values = c('green','hotpink'),
            #na.value='white', labels=c('H1','H2')) +
            theme_bw()+
            theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.title=element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size=12)) + 
            #legend.position = 'none',axis.title.y=element_blank(), 
            xlab('Genomic position')+ylab('')  #+ 
        
    } else {
        gapYMin <- gap[,2]
        gapYMax <- gap[,3]
        
        gapPlot <- data.frame(gapYMin = gapYMin/1000000,
            gapYMax = gapYMax/1000000,
            gapXMin=rep(0,nChr), gapXMax=rep(1.2,nChr),
            chr = gap[,1])
        
        gapPlot$chr <- factor(gapPlot$chr, 
            labels = paste('chr',gapPlot$chr,sep=''))
        
        ggplot()+
            geom_line(data = cvPlot,
                aes(x = x, y = y+1.25), colour = 'blue', 
                stat="identity", size=0.6) + 
            geom_rect(data=chrPlot,aes(xmin=yMin, xmax=yMax, 
                ymin=xMin, ymax=xMax), 
                color='limegreen',fill='limegreen',size=0.2) + 
            geom_rect(data=chrPlot, aes(xmin=cenYMin, xmax=cenYMax, 
                ymin=cenXMin, ymax=cenXMax), 
                color='black',fill='black',size=0.2) +
            geom_rect(data=gapPlot, aes(xmin=gapYMin, xmax=gapYMax, 
                ymin=gapXMin, ymax=gapXMax), 
                color='black',fill='white',size=0.2) +
            geom_point(data=cvSitePlot, aes(x=pos,y=y), 
                color='red',fill='red',size=2.5, shape=4, stroke =0.8) +
            coord_flip() +  facet_grid(.~chr, scales = "free", space='free') +
            
            #geom_point(data=cenDa, aes(x=cenX[1], y=cenY[1]), 
            #size=3, shape=1, stroke=2) +
            scale_x_reverse(breaks = y.breaks, limits = c(maxChr/1000000,0),
                labels = y.labels) +
            scale_y_continuous(limits = c(0,x.limits+1.25), 
                breaks=0:x.limits+1.25, labels=0:x.limits) +
            #scale_color_manual(values = c('green','hotpink'),
            #na.value='white', labels=c('H1','H2')) +
            theme_bw()+
            theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.title=element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size=12),
                strip.text = element_text(size=12)) + 
            #legend.position = 'none',axis.title.y=element_blank(), 
            xlab('Genomic position')+ylab('')  #+ 
    }
}
