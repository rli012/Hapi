##' @title Visualization of haplotypes in a single gamete cell
##' @description Visualization of haplotypes in a single gamete cell
##' @param hap a dataframe of all the phased hetSNPs in all chromosomes
##' @param chr a dataframe of chromosome information, including length, 
##' and centrometric regions
##' @param hap.color a vector of colors for the two haplotypes. 
##' Default is \code{c('deepskyblue2','darkorange2')}
##' @param centromere.fill a character of the color for the centromeres. 
##' Default is \code{'black'}
##' @param x.breaks a vector of positions to show labels on x axis. 
##' Default is \code{NULL}
##' @param x.labels a vector of labels on the x axis. 
##' Default is \code{NULL}
##' @param y.breaks a vector of positions to show labels on y axis. 
##' Default is \code{NULL}
##' @param y.labels a vector of labels on the y axis. 
##' Default is \code{NULL}
##' @import ggplot2
##' @return a plot of haplotypes in a single gamete cell
##' @export
##' @author Ruidong Li
##' @examples
##' data(gamete11)
##' \dontrun{hapiGameteView(hap=gamete11)}
hapiGameteView <- function(hap, chr=hg19, 
    hap.color=c('deepskyblue2','darkorange2'),
    centromere.fill='black', x.breaks=NULL,
    x.labels=NULL, y.breaks=NULL, y.labels=NULL) {
    
    nChr <- nrow(chr)
    nSNPs <- sapply(1:nChr, function(x) sum(hap$chr==x))

    ypos <- hap$pos
    
    xstart <- rep(1:nChr-0.25, nSNPs)
    xend <- rep(1:nChr+0.25, nSNPs)

    genoStatus <- hap$hap

    chrSize <- chr$length
    
    maxChr <- max(chrSize)
    
    htmpGeno <- data.frame(ypos, genoStatus, xstart, xend)
    
    chrLength <- chr$length
    cenStart <- chr$cenStart
    cenEnd <- chr$cenEnd
    
    xMin <- 1:nChr-0.25
    xMax <- 1:nChr+0.25
    yMin <- rep(0,nChr)
    yMax <- chrLength
    
    boxDa <- data.frame(xMin,xMax,yMin,yMax)
    
    
    cenXMin <- xMin
    cenXMax <- xMax
    
    cenYMin <- cenStart
    cenYMax <- cenEnd
    
    cenDa <- data.frame(cenXMin,cenXMax,cenYMin,cenYMax)
    
    ####
    cenX <- 1:nChr
    cenY <- (cenStart+cenEnd)/2
    
    cenDa <- data.frame(cenX, cenY)
    
    
    if (is.null(x.breaks)) {
        x.breaks = 1:nChr
    }
    if (is.null(x.labels)) {
        x.labels = 1:nChr
    }
    
    if (is.null(y.breaks)) {
        y.breaks = c(0,5,10,15,20)*10000000
    }
    if (is.null(y.labels)) {
        y.labels = c(0,'50Mb','100Mb','150Mb','200Mb')
    }

    ggplot() +
        geom_segment(data=htmpGeno, 
            aes(x=xstart, y=ypos, yend = htmpGeno$ypos, xend = htmpGeno$xend, 
            col=as.factor(genoStatus)), size=1) +
        geom_rect(data=cenDa,aes(xmin=cenXMin, xmax=cenXMax, 
            ymin=cenYMin, ymax=cenYMax), 
            color=NA,fill=centromere.fill,size=0.3) +
        geom_rect(data=boxDa,aes(xmin=xMin, xmax=xMax, ymin=yMin, ymax=yMax), 
            color='black',fill=NA,size=0.3) +
        
        #geom_point(data=cenDa, aes(x=cenX, y=cenY), size=3, 
        #shape=1, stroke=2) +
        scale_y_reverse(breaks = y.breaks, limits = c(maxChr,0),
            labels = y.labels) +
        #scale_color_manual(values = c('green','hotpink'),
        #na.value='white', labels=c('H1','H2')) +
        scale_color_manual(values = hap.color,na.value='white', 
            labels=c('H1','H2')) +
        theme_bw()+
        theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.title=element_blank(), 
            legend.position = 'none', 
            axis.title=element_text(size=16)) +
        #xlab('Chromosome')+ylab('Genomic position')  + 
        #coord_flip()+
        #scale_x_continuous(breaks=1:22,labels=paste('Chr',1:22, sep='')) +
        scale_x_continuous(breaks = x.breaks,labels = x.labels) +
        theme(axis.text = element_text(size=14)) +
        xlab('Chromosome') + ylab('Genomic Position')
}








