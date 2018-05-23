###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}


##############################################################
##############################################################
##  Final figures for paper







##############################################################
##############################################################
##  Exploratory figures

###########################################################################################
#' Deterministic Simulations: Equilibrium frequencies of Autosomal vs. X-linked inversions
#'  -- Equal selection through each sex (no sex-specific selection for Autosomal model)
#'  -- Additive fitness effects at SA locus
#' @title Eq. frequencies of Auto vs. X-linked inversions
#' @author Colin Olito
#' @export
autoVsXEqFreqs  <-  function() {
    
    # Import data for plotting
    auto  <-  read.csv("./output/data/simResults/determRecSim-Autosomal_r0.01.csv", header=TRUE)
    xchr  <-  read.csv("./output/data/simResults/determRecSim-X-linked_r0.01.csv", header=TRUE)

    # Calculate equilibrium frequencies for X-linked model 
    xEQ.freqs  <-  (2*xchr[,7:11] + xchr[,12:16])/3
    colnames (xEQ.freqs)  <-  c("f1","f2","f3","f4","f5")
    xchr       <-  cbind(xchr,xEQ.freqs)
    XAutodiff  <-  xchr$f5 - auto$x5

    # Factor levels for subsetting  
    hs  <-  unique(auto$h)
    ms  <-  unique(auto$m)

    # Color scheme
    COLS     <-  c(transparentColor('dodgerblue', opacity=1),
                   transparentColor('darkolivegreen', opacity=1),
                   transparentColor('tomato', opacity=1),
                   transparentColor('dodgerblue2', opacity=1))
    COLS.bg  <-  c(transparentColor('dodgerblue', opacity=0.5),
                   transparentColor('darkolivegreen', opacity=0.5),
                   transparentColor('tomato', opacity=0.5),
                   transparentColor('dodgerblue2', opacity=0.5))

    # Set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)


## Row 1: Equilibrium frequencies
    # Panel 1: m = 0.01 
    par(omi=c(0.5, 0.5, 0.5, 0.5), mar = c(4,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(auto$s)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        lines(auto$x5[auto$h==hs[1] & auto$m==ms[1]] ~ auto$s[auto$h==hs[1] & auto$m==ms[1]], col=COLS[1], lwd=3, lty=1)
        lines(auto$x5[auto$h==hs[2] & auto$m==ms[1]] ~ auto$s[auto$h==hs[2] & auto$m==ms[1]], col=COLS[2], lwd=3, lty=1)
        lines(auto$x5[auto$h==hs[3] & auto$m==ms[1]] ~ auto$s[auto$h==hs[3] & auto$m==ms[1]], col=COLS[3], lwd=3, lty=1)
        # X-linked frequencies
        lines(xchr$f5[xchr$h==hs[1] & xchr$mf==ms[1]] ~ xchr$sf[xchr$h==hs[1] & xchr$mf==ms[1]], col=COLS[1], lwd=3, lty=2)
        lines(xchr$f5[xchr$h==hs[2] & xchr$mf==ms[1]] ~ xchr$sf[xchr$h==hs[2] & xchr$mf==ms[1]], col=COLS[2], lwd=3, lty=2)
        lines(xchr$f5[xchr$h==hs[3] & xchr$mf==ms[1]] ~ xchr$sf[xchr$h==hs[3] & xchr$mf==ms[1]], col=COLS[3], lwd=3, lty=2)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.45,  1.1,   expression(paste(italic(m)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.655, 1.11,  substitute(m,list(m=ms[1])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Eq. inversion frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
              x       =  usr[2]*0.825,
              y       =  usr[4]*0.5,
              legend  =  c(
                          expression(paste(italic(h)~"="~0)),
                          expression(paste(italic(h)~"="~1/2)),
                          expression(paste(italic(h)~"="~1))),
#              pch     =  c(21,21,21),
#              pt.bg   =  c(COLS[1],COLS[2],COLS[3]),
              col     =  c(COLS[1],COLS[2],COLS[3]),
              lty     =  c(1,1,1),
              lwd     =  c(3,3,3),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.2,
              legend  =  c(
                          expression(paste("X-linked")),
                          expression(paste("Autosomal"))),
#              pch     =  c(21,21,21),
#              pt.bg   =  c(COLS[1],COLS[2],COLS[3]),
              col     =  1,
              lty     =  c(2,1),
              lwd     =  c(3,3),
              seg.len  =  3,
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )

    # Panel 2: m = 0.05     
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(auto$s)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        lines(auto$x5[auto$h==hs[1] & auto$m==ms[2]] ~ auto$s[auto$h==hs[1] & auto$m==ms[2]], col=COLS[1], lwd=3, lty=1)
        lines(auto$x5[auto$h==hs[2] & auto$m==ms[2]] ~ auto$s[auto$h==hs[2] & auto$m==ms[2]], col=COLS[2], lwd=3, lty=1)
        lines(auto$x5[auto$h==hs[3] & auto$m==ms[2]] ~ auto$s[auto$h==hs[3] & auto$m==ms[2]], col=COLS[3], lwd=3, lty=1)
        # X-linked frequencies
        lines(xchr$f5[xchr$h==hs[1] & xchr$mf==ms[2]] ~ xchr$sf[xchr$h==hs[1] & xchr$mf==ms[2]], col=COLS[1], lwd=3, lty=2)
        lines(xchr$f5[xchr$h==hs[2] & xchr$mf==ms[2]] ~ xchr$sf[xchr$h==hs[2] & xchr$mf==ms[2]], col=COLS[2], lwd=3, lty=2)
        lines(xchr$f5[xchr$h==hs[3] & xchr$mf==ms[2]] ~ xchr$sf[xchr$h==hs[3] & xchr$mf==ms[2]], col=COLS[3], lwd=3, lty=2)
        # Axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.45,  1.1,   expression(paste(italic(m)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.655, 1.11,  substitute(m,list(m=ms[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)


## Row 2: Difference in X vs. Autosome inversion frequencies
    # Panel 3: m = 0.01 
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(auto$s)), ylim = c(0,max(XAutodiff)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Difference in X vs. Autosomal inversion frequencies
        lines(XAutodiff[auto$h==hs[1] & auto$m==ms[1]] ~ auto$s[auto$h==hs[1] & auto$m==ms[1]], col=COLS[1], lwd=3, lty=1)
        lines(XAutodiff[auto$h==hs[2] & auto$m==ms[1]] ~ auto$s[auto$h==hs[2] & auto$m==ms[1]], col=COLS[2], lwd=3, lty=1)
        lines(XAutodiff[auto$h==hs[3] & auto$m==ms[1]] ~ auto$s[auto$h==hs[3] & auto$m==ms[1]], col=COLS[3], lwd=3, lty=1)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(Delta," Eq. inversion frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,   expression(paste(italic(s))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Panel 4: m = 0.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(auto$s)), ylim = c(0,max(XAutodiff)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Difference in X vs. Autosomal inversion frequencies
        lines(XAutodiff[auto$h==hs[1] & auto$m==ms[2]] ~ auto$s[auto$h==hs[1] & auto$m==ms[2]], col=COLS[1], lwd=3, lty=1)
        lines(XAutodiff[auto$h==hs[2] & auto$m==ms[2]] ~ auto$s[auto$h==hs[2] & auto$m==ms[2]], col=COLS[2], lwd=3, lty=1)
        lines(XAutodiff[auto$h==hs[3] & auto$m==ms[2]] ~ auto$s[auto$h==hs[3] & auto$m==ms[2]], col=COLS[3], lwd=3, lty=1)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,   expression(paste(italic(s))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
} 




###########################################################################################
#' Fig. 4 -- Enrichment of fixed vs. polymorphic inversions
#' @title Fig. 4
#' @author Colin Olito
#' @export
Fig4  <-  function() {

  # Import data from supplementary tables S1 & S2
  fixDat   <-  read.csv("./data/TableS1-plotting.csv", header=TRUE)
  polyDat  <-  read.csv("./data/TableS2-plotting.csv", header=TRUE)

  # Exclude species w/ fewer than 10 known inversions
  fixDat   <-  fixDat[(fixDat$nAinv+fixDat$nXinv) >= 10,]
  polyDat  <-  polyDat[(polyDat$nAinv+polyDat$nXinv) >= 10,]

  # Make plot
  par(omi=rep(0.5, 4), mar = c(3.5,3.5,2,2), bty='o', xaxt='s', yaxt='s')
  plot(NA, axes=FALSE, type='n', main='',xlim = c(0,3.5), ylim = c(0,6), ylab='', xlab='', cex.lab=1.2)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
  plotGrid(lineCol='grey80')
  box()
  # Equilibrium frequencies for different 
  hist(fixDat$enrichX, breaks=15, col=transparentColor('tomato', opacity=0.5), add=TRUE)
  hist(polyDat$enrichX, breaks=7, col=transparentColor('dodgerblue', opacity=0.5), add=TRUE)
  abline(v=1, lty=2, lwd=3)
  # axes
  axis(1, las=1)
  axis(2, las=1)
  proportionalLabel(0.5, -0.15, expression(paste("Enrichment on the X")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
  proportionalLabel(-0.15, 0.5, expression(paste(Count)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
#  proportionalLabel(0.75, 0.9, expression(paste(Inversion~type)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
  #legend
  legend(
        x       =  usr[2]*0.975,
        y       =  usr[4]*0.975,
        legend  =  c(
                    expression(paste("Fixed")),
                    expression(paste("Polymorphic"))),
        fill    =  c(transparentColor('tomato', opacity=0.5),
                     transparentColor('dodgerblue', opacity=0.5)),
        cex     =  1.2,
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA
        )
}



###########################################################################################
#' Deterministic Simulations: Equilibrium frequencies of Autosomal vs. X-linked inversions
#'  -- Equal selection through each sex (no sex-specific selection for Autosomal model)
#'  -- Additive fitness effects at SA locus
#' @title Eq. frequencies of Auto vs. X-linked inversions
#' @author Colin Olito
#' @export
propEstSuppFigs  <-  function() {

    # Import data for plotting
      aDataN30kH0    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/SexSpecFig_sexspM_sexspS_sexspR_N30000_h0_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      xDataN30kH0    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/XlinkedFig_sexspM_sexspecS_N30000_h0_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      aDataN30kH0.5  <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/SexSpecFig_sexspM_sexspS_sexspR_N30000_h0.5_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      xDataN30kH0.5  <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/XlinkedFig_sexspM_sexspecS_N30000_h0.5_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      aDataN30kH1    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/SexSpecFig_sexspM_sexspS_sexspR_N30000_h1_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      xDataN30kH1    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/XlinkedFig_sexspM_sexspecS_N30000_h1_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)

      aDataH0    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/SexSpecFig_sexspM_sexspS_sexspR_N5e+05_h0_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      xDataH0    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/XlinkedFig_sexspM_sexspecS_N5e+05_h0_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      aDataH0.5  <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/SexSpecFig_sexspM_sexspS_sexspR_N5e+05_h0.5_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
#      xDataH0.5  <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/XlinkedFig_sexspM_sexspecS_N5e+05_h0.5_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      aDataH1    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/SexSpecFig_sexspM_sexspS_sexspR_N5e+05_h1_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)
      xDataH1    <-  read.csv(file = "/Volumes/VERBATIM\ HD/XvAutosomeSuppData/XlinkedFig_sexspM_sexspecS_N5e+05_h1_n100_u1e-05_sDel_none_nReps1e+06.csv", header=TRUE)

    # indexes for plotting
    mfs   <-  unique(aDataH0$mf)
    mms   <-  unique(aDataH0$mm)
    sfs   <-  unique(aDataH0$sf)
    sms   <-  unique(aDataH0$sm)
    rfs   <-  unique(aDataH0$rf)
    rms   <-  unique(aDataH0$rm)
    mfsX  <-  unique(xDataH0$mf)
    mmsX  <-  unique(xDataH0$mm)
    sfsX  <-  unique(xDataH0$sf)
    smsX  <-  unique(xDataH0$sm)
    rsX   <-  unique(xDataH0$rf)

    # Color scheme
    COLS  <-  c('grey70',2)

    # Set plot layout
    layout.mat <- matrix(
                        c(1,1,2,2,3,3,10,11,11,12,12,13,13,
                          1,1,2,2,3,3,10,11,11,12,12,13,13,
                          4,4,5,5,6,6,10,14,14,15,15,16,16,
                          4,4,5,5,6,6,10,14,14,15,15,16,16,
                          7,7,8,8,9,9,10,17,17,18,18,19,19,
                          7,7,8,8,9,9,10,17,17,18,18,19,19
                          ), 
                        nrow=6, ncol=13, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)
#    layout.show(layout)

## Row 1: Equal and Sex-Specific Migration Rates
    # Panel 1: h_j = 0
    aDat   <-  aDataN30kH0[aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$sf == aDataN30kH0$sm & aDataN30kH0$rf == aDataN30kH0$rm,][-c(12,13),]
    xDat   <-  xDataN30kH0[xDataN30kH0$mf == xDataN30kH0$mm & xDataN30kH0$sf == xDataN30kH0$sm,]
    aDat2  <-  aDataN30kH0[aDataN30kH0$mf >  aDataN30kH0$mm & aDataN30kH0$sf == aDataN30kH0$sm & aDataN30kH0$rf == aDataN30kH0$rm,][-c(12,13),]
    xDat2  <-  xDataN30kH0[xDataN30kH0$mf >  xDataN30kH0$mm & xDataN30kH0$sf == xDataN30kH0$sm,]
    aDat3  <-  aDataN30kH0[aDataN30kH1$mf <  aDataN30kH0$mm & aDataN30kH0$sf == aDataN30kH0$sm & aDataN30kH0$rf == aDataN30kH0$rm,][-c(12,13),]
    xDat3  <-  xDataN30kH0[xDataN30kH0$mf <  xDataN30kH0$mm & xDataN30kH0$sf == xDataN30kH0$sm,]
#mar = c(4,4,0.5,0.5)
    par(omi=c(0.5, 1, 0.5, 0.5), mar = c(4,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH0$PropEst, xDataN30kH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(s[f])," = ",italic(s[m]),", ",italic(r[f])," = ",italic(r[m]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.05,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(h[j])," = ",0)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.655, 1.11,  substitute(m,list(m=ms[1])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(italic(Pi[j])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # Panel 2: h_j = 1/2
    aDat   <-  aDataN30kH0.5[aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$sf == aDataN30kH0.5$sm & aDataN30kH0.5$rf == aDataN30kH0.5$rm,][-c(12,13),]
    xDat   <-  xDataN30kH0.5[xDataN30kH0.5$mf == xDataN30kH0.5$mm & xDataN30kH0.5$sf == xDataN30kH0.5$sm,]
    aDat2  <-  aDataN30kH0.5[aDataN30kH0.5$mf >  aDataN30kH0.5$mm & aDataN30kH0.5$sf == aDataN30kH0.5$sm & aDataN30kH0.5$rf == aDataN30kH0.5$rm,][-c(12,13),]
    xDat2  <-  xDataN30kH0.5[xDataN30kH0.5$mf >  xDataN30kH0.5$mm & xDataN30kH0.5$sf == xDataN30kH0.5$sm,]
    aDat3  <-  aDataN30kH0.5[aDataN30kH0.5$mf <  aDataN30kH0.5$mm & aDataN30kH0.5$sf == aDataN30kH0.5$sm & aDataN30kH0.5$rf == aDataN30kH0.5$rm,][-c(12,13),]
    xDat3  <-  xDataN30kH0.5[xDataN30kH0.5$mf <  xDataN30kH0.5$mm & xDataN30kH0.5$sf == xDataN30kH0.5$sm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH0.5$PropEst, xDataN30kH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.4,   expression(paste(italic(N), " = ", 3%*%10^4)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.05,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(h[j])," = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Panel 3: h_j = 1
    aDat   <-  aDataN30kH1[aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$sf == aDataN30kH1$sm & aDataN30kH1$rf == aDataN30kH1$rm,][-c(12,13),]
    xDat   <-  xDataN30kH1[xDataN30kH1$mf == xDataN30kH1$mm & xDataN30kH1$sf == xDataN30kH1$sm,]
    aDat2  <-  aDataN30kH1[aDataN30kH1$mf >  aDataN30kH1$mm & aDataN30kH1$sf == aDataN30kH1$sm & aDataN30kH1$rf == aDataN30kH1$rm,][-c(12,13),]
    xDat2  <-  xDataN30kH1[xDataN30kH1$mf >  xDataN30kH1$mm & xDataN30kH1$sf == xDataN30kH1$sm,]
    aDat3  <-  aDataN30kH1[aDataN30kH1$mf <  aDataN30kH1$mm & aDataN30kH1$sf == aDataN30kH1$sm & aDataN30kH1$rf == aDataN30kH1$rm,][-c(12,13),]
    xDat3  <-  xDataN30kH1[xDataN30kH1$mf <  xDataN30kH1$mm & xDataN30kH1$sf == xDataN30kH1$sm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH1$PropEst, xDataN30kH1$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(h[j])," = ",1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.05,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*0.42,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              pch     =  c(21,21),
              pt.bg   =  c('grey70',2),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.24,
              legend  =  c(
                          expression(paste(italic(m[f])," = ",italic(m[m]))),
                          expression(paste(italic(m[f])>italic(m[m]))),
                          expression(paste(italic(m[f])<  italic(m[m])))),
              pch     =  c(21,22,23),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )


## Row 2: Equal and Sex-Specific Selection
    # Panel 4: h_j = 0
    aDat   <-  aDataN30kH0[aDataN30kH0$sf == aDataN30kH0$sm & aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$rf == aDataN30kH0$rm,][-c(12,13),]
    xDat   <-  xDataN30kH0[xDataN30kH0$sf == xDataN30kH0$sm & xDataN30kH0$mf == xDataN30kH0$mm,]
    aDat2  <-  aDataN30kH0[aDataN30kH0$sf >  aDataN30kH0$sm & aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$rf == aDataN30kH0$rm,][-c(12,13),]
    xDat2  <-  xDataN30kH0[xDataN30kH0$sf >  xDataN30kH0$sm & xDataN30kH0$mf == xDataN30kH0$mm,]
    aDat3  <-  aDataN30kH0[aDataN30kH0$sf <  aDataN30kH0$sm & aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$rf == aDataN30kH0$rm,][-c(12,13),]
    xDat3  <-  xDataN30kH0[xDataN30kH0$sf <  xDataN30kH0$sm & xDataN30kH0$mf == xDataN30kH0$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH0$PropEst, xDataN30kH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(m[f])," = ",italic(m[m]),", ",italic(r[f])," = ",italic(r[m]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.05,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(italic(Pi[j])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        


    # Panel 5: h_j = 1/2
    aDat   <-  aDataN30kH0.5[aDataN30kH0.5$sf == aDataN30kH0.5$sm & aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$rf == aDataN30kH0.5$rm,][-c(12,13),]
    xDat   <-  xDataN30kH0.5[xDataN30kH0.5$sf == xDataN30kH0.5$sm & xDataN30kH0.5$mf == xDataN30kH0.5$mm,]
    aDat2  <-  aDataN30kH0.5[aDataN30kH0.5$sf >  aDataN30kH0.5$sm & aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$rf == aDataN30kH0.5$rm,][-c(12,13),]
    xDat2  <-  xDataN30kH0.5[xDataN30kH0.5$sf >  xDataN30kH0.5$sm & xDataN30kH0.5$mf == xDataN30kH0.5$mm,]
    aDat3  <-  aDataN30kH0.5[aDataN30kH0.5$sf <  aDataN30kH0.5$sm & aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$rf == aDataN30kH0.5$rm,][-c(12,13),]
    xDat3  <-  xDataN30kH0.5[xDataN30kH0.5$sf <  xDataN30kH0.5$sm & xDataN30kH0.5$mf == xDataN30kH0.5$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH0.5$PropEst, xDataN30kH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Panel 6: h_j = 1
    aDat   <-  aDataN30kH1[aDataN30kH1$sf == aDataN30kH1$sm & aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$rf == aDataN30kH1$rm,][-c(12,13),]
    xDat   <-  xDataN30kH1[xDataN30kH1$sf == xDataN30kH1$sm & xDataN30kH1$mf == xDataN30kH1$mm,]
    aDat2  <-  aDataN30kH1[aDataN30kH1$sf >  aDataN30kH1$sm & aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$rf == aDataN30kH1$rm,][-c(12,13),]
    xDat2  <-  xDataN30kH1[xDataN30kH1$sf >  xDataN30kH1$sm & xDataN30kH1$mf == xDataN30kH1$mm,]
    aDat3  <-  aDataN30kH1[aDataN30kH1$sf <  aDataN30kH1$sm & aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$rf == aDataN30kH1$rm,][-c(12,13),]
    xDat3  <-  xDataN30kH1[xDataN30kH1$sf <  xDataN30kH1$sm & xDataN30kH1$mf == xDataN30kH1$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH1$PropEst, xDataN30kH1$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*0.42,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              pch     =  c(21,21),
              pt.bg   =  c('grey70',2),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.24,
              legend  =  c(
                          expression(paste(italic(s[f])," = ",italic(s[m]))),
                          expression(paste(italic(s[f])>italic(s[m]))),
                          expression(paste(italic(s[f])<italic(s[m])))),
              pch     =  c(21,22,23),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )


## Row 3: Equal and Sex-Specific Recombination
    # Panel 7: h_j = 0
    aDat   <-  aDataN30kH0[aDataN30kH0$rf == aDataN30kH0$rm & aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$sf == aDataN30kH0$sm,][-c(12,13),]
    aDat2  <-  aDataN30kH0[aDataN30kH0$rf >  aDataN30kH0$rm & aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$sf == aDataN30kH0$sm,][-c(12,13),]
    aDat3  <-  aDataN30kH0[aDataN30kH0$rf <  aDataN30kH0$rm & aDataN30kH0$mf == aDataN30kH0$mm & aDataN30kH0$sf == aDataN30kH0$sm,][-c(12,13),]
    xDat   <-  xDataN30kH0[xDataN30kH0$sf == xDataN30kH0$sm & xDataN30kH0$mf == xDataN30kH0$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH0$PropEst, xDataN30kH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rm, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(m[f])," = ",italic(m[m]),", ",italic(s[f])," = ",italic(s[m]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.05,  1.075, expression(paste(bold(G))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(italic(Pi[j])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(0.5,  -0.3,   expression(paste("Recombination rate ",italic((r[j])))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


    # Panel 8: h_j = 1/2
    aDat   <-  aDataN30kH0.5[aDataN30kH0.5$rf == aDataN30kH0.5$rm & aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$sf == aDataN30kH0.5$sm,][-c(12,13),]
    aDat2  <-  aDataN30kH0.5[aDataN30kH0.5$rf >  aDataN30kH0.5$rm & aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$sf == aDataN30kH0.5$sm,][-c(12,13),]
    aDat3  <-  aDataN30kH0.5[aDataN30kH0.5$rf <  aDataN30kH0.5$rm & aDataN30kH0.5$mf == aDataN30kH0.5$mm & aDataN30kH0.5$sf == aDataN30kH0.5$sm,][-c(12,13),]
    xDat   <-  xDataN30kH0.5[xDataN30kH0.5$sf == xDataN30kH0.5$sm & xDataN30kH0.5$mf == xDataN30kH0.5$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH0.5$PropEst, xDataN30kH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rm, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(H))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste("Recombination rate ",italic((r[j])))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

    # Panel 9: h_j = 1
    aDat   <-  aDataN30kH1[aDataN30kH1$rf == aDataN30kH1$rm & aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$sf == aDataN30kH1$sm,][-c(12,13),]
    aDat2  <-  aDataN30kH1[aDataN30kH1$rf >  aDataN30kH1$rm & aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$sf == aDataN30kH1$sm,][-c(12,13),]
    aDat3  <-  aDataN30kH1[aDataN30kH1$rf <  aDataN30kH1$rm & aDataN30kH1$mf == aDataN30kH1$mm & aDataN30kH1$sf == aDataN30kH1$sm,][-c(12,13),]
    xDat   <-  xDataN30kH1[xDataN30kH1$sf == xDataN30kH1$sm & xDataN30kH1$mf == xDataN30kH1$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH1$PropEst, xDataN30kH1$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rm, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(I))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste("Recombination rate ",italic((r[j])))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*0.42,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              pch     =  c(21,21),
              pt.bg   =  c('grey70',2),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.24,
              legend  =  c(
                          expression(paste(italic(r[f])," = ",italic(r[m]))),
                          expression(paste(italic(r[f])," > ",italic(r[m]))),
                          expression(paste(italic(r[f])," < ",italic(r[m])))),
              pch     =  c(21,22,23),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )

#########################################
plot.new()
#     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataN30kH1$PropEst, xDataN30kH1$PropEst)), ylab='', xlab='', cex.lab=1.2, xpd=NA)
#########################################

## Row 1: Equal and Sex-Specific Migration Rates
    # Panel 1: h_j = 0
    aDat   <-  aDataH0[aDataH0$mf == aDataH0$mm & aDataH0$sf == aDataH0$sm & aDataH0$rf == aDataH0$rm,][-c(12,13),]
    xDat   <-  xDataH0[xDataH0$mf == xDataH0$mm & xDataH0$sf == xDataH0$sm,]
    aDat2  <-  aDataH0[aDataH0$mf >  aDataH0$mm & aDataH0$sf == aDataH0$sm & aDataH0$rf == aDataH0$rm,][-c(12,13),]
    xDat2  <-  xDataH0[xDataH0$mf >  xDataH0$mm & xDataH0$sf == xDataH0$sm,]
    aDat3  <-  aDataH0[aDataH1$mf <  aDataH0$mm & aDataH0$sf == aDataH0$sm & aDataH0$rf == aDataH0$rm,][-c(12,13),]
    xDat3  <-  xDataH0[xDataH0$mf <  xDataH0$mm & xDataH0$sf == xDataH0$sm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH0$PropEst, xDataH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(s[f])," = ",italic(s[m]),", ",italic(r[f])," = ",italic(r[m]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.05,  1.075, expression(paste(bold(J))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(h[j])," = ",0)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.655, 1.11,  substitute(m,list(m=ms[1])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(italic(Pi[j])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # Panel 2: h_j = 1/2
    aDat   <-  aDataH0.5[aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$sf == aDataH0.5$sm & aDataH0.5$rf == aDataH0.5$rm,][-c(12,13),]
    xDat   <-  xDataH0[xDataH0$mf == xDataH0$mm & xDataH0$sf == xDataH0$sm,]
    aDat2  <-  aDataH0.5[aDataH0.5$mf >  aDataH0.5$mm & aDataH0.5$sf == aDataH0.5$sm & aDataH0.5$rf == aDataH0.5$rm,][-c(12,13),]
    xDat2  <-  xDataH0[xDataH0$mf >  xDataH0$mm & xDataH0$sf == xDataH0$sm,]
    aDat3  <-  aDataH0.5[aDataH0.5$mf <  aDataH0.5$mm & aDataH0.5$sf == aDataH0.5$sm & aDataH0.5$rf == aDataH0.5$rm,][-c(12,13),]
    xDat3  <-  xDataH0[xDataH0$mf <  xDataH0$mm & xDataH0$sf == xDataH0$sm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH0.5$PropEst, xDataH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
#        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
#        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
#        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.4,   expression(paste(italic(N), " = ", 5%*%10^5)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.05,  1.075, expression(paste(bold(K))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(h[j])," = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Panel 3: h_j = 1
    aDat   <-  aDataH1[aDataH1$mf == aDataH1$mm & aDataH1$sf == aDataH1$sm & aDataH1$rf == aDataH1$rm,][-c(12,13),]
    xDat   <-  xDataH1[xDataH1$mf == xDataH1$mm & xDataH1$sf == xDataH1$sm,]
    aDat2  <-  aDataH1[aDataH1$mf >  aDataH1$mm & aDataH1$sf == aDataH1$sm & aDataH1$rf == aDataH1$rm,][-c(12,13),]
    xDat2  <-  xDataH1[xDataH1$mf >  xDataH1$mm & xDataH1$sf == xDataH1$sm,]
    aDat3  <-  aDataH1[aDataH1$mf <  aDataH1$mm & aDataH1$sf == aDataH1$sm & aDataH1$rf == aDataH1$rm,][-c(12,13),]
    xDat3  <-  xDataH1[xDataH1$mf <  xDataH1$mm & xDataH1$sf == xDataH1$sm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH1$PropEst, xDataH1$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(h[j])," = ",1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.05,  1.075, expression(paste(bold(L))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*0.42,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              pch     =  c(21,21),
              pt.bg   =  c('grey70',2),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.24,
              legend  =  c(
                          expression(paste(italic(m[f])," = ",italic(m[m]))),
                          expression(paste(italic(m[f])>italic(m[m]))),
                          expression(paste(italic(m[f])<  italic(m[m])))),
              pch     =  c(21,22,23),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )


## Row 2: Equal and Sex-Specific Selection
    # Panel 4: h_j = 0
    aDat   <-  aDataH0[aDataH0$sf == aDataH0$sm & aDataH0$mf == aDataH0$mm & aDataH0$rf == aDataH0$rm,][-c(12,13),]
    xDat   <-  xDataH0[xDataH0$sf == xDataH0$sm & xDataH0$mf == xDataH0$mm,]
    aDat2  <-  aDataH0[aDataH0$sf >  aDataH0$sm & aDataH0$mf == aDataH0$mm & aDataH0$rf == aDataH0$rm,][-c(12,13),]
    xDat2  <-  xDataH0[xDataH0$sf >  xDataH0$sm & xDataH0$mf == xDataH0$mm,]
    aDat3  <-  aDataH0[aDataH0$sf <  aDataH0$sm & aDataH0$mf == aDataH0$mm & aDataH0$rf == aDataH0$rm,][-c(12,13),]
    xDat3  <-  xDataH0[xDataH0$sf <  xDataH0$sm & xDataH0$mf == xDataH0$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH0$PropEst, xDataH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(m[f])," = ",italic(m[m]),", ",italic(r[f])," = ",italic(r[m]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.05,  1.075, expression(paste(bold(M))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(italic(Pi[j])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # Panel 5: h_j = 1/2
    aDat   <-  aDataH0.5[aDataH0.5$sf == aDataH0.5$sm & aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$rf == aDataH0.5$rm,][-c(12,13),]
    xDat   <-  xDataH0[xDataH0$sf == xDataH0$sm & xDataH0$mf == xDataH0$mm,]
    aDat2  <-  aDataH0.5[aDataH0.5$sf >  aDataH0.5$sm & aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$rf == aDataH0.5$rm,][-c(12,13),]
    xDat2  <-  xDataH0[xDataH0$sf >  xDataH0$sm & xDataH0$mf == xDataH0$mm,]
    aDat3  <-  aDataH0.5[aDataH0.5$sf <  aDataH0.5$sm & aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$rf == aDataH0.5$rm,][-c(12,13),]
    xDat3  <-  xDataH0[xDataH0$sf <  xDataH0$sm & xDataH0$mf == xDataH0$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH0.5$PropEst, xDataH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
#        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
#        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
#        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(N))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Panel 6: h_j = 1
    aDat   <-  aDataH1[aDataH1$sf == aDataH1$sm & aDataH1$mf == aDataH1$mm & aDataH1$rf == aDataH1$rm,][-c(12,13),]
    xDat   <-  xDataH1[xDataH1$sf == xDataH1$sm & xDataH1$mf == xDataH1$mm,]
    aDat2  <-  aDataH1[aDataH1$sf >  aDataH1$sm & aDataH1$mf == aDataH1$mm & aDataH1$rf == aDataH1$rm,][-c(12,13),]
    xDat2  <-  xDataH1[xDataH1$sf >  xDataH1$sm & xDataH1$mf == xDataH1$mm,]
    aDat3  <-  aDataH1[aDataH1$sf <  aDataH1$sm & aDataH1$mf == aDataH1$mm & aDataH1$rf == aDataH1$rm,][-c(12,13),]
    xDat3  <-  xDataH1[xDataH1$sf <  xDataH1$sm & xDataH1$mf == xDataH1$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH1$PropEst, xDataH1$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rf, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(O))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*0.42,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              pch     =  c(21,21),
              pt.bg   =  c('grey70',2),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.24,
              legend  =  c(
                          expression(paste(italic(s[f])," = ",italic(s[m]))),
                          expression(paste(italic(s[f])>italic(s[m]))),
                          expression(paste(italic(s[f])<italic(s[m])))),
              pch     =  c(21,22,23),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )



## Row 3: Equal and Sex-Specific Recombination
    # Panel 4: h_j = 0
    aDat   <-  aDataH0[aDataH0$rf == aDataH0$rm & aDataH0$mf == aDataH0$mm & aDataH0$sf == aDataH0$sm,][-c(12,13),]
    aDat2  <-  aDataH0[aDataH0$rf >  aDataH0$rm & aDataH0$mf == aDataH0$mm & aDataH0$sf == aDataH0$sm,][-c(12,13),]
    aDat3  <-  aDataH0[aDataH0$rf <  aDataH0$rm & aDataH0$mf == aDataH0$mm & aDataH0$sf == aDataH0$sm,][-c(12,13),]
    xDat   <-  xDataH0[xDataH0$sf == xDataH0$sm & xDataH0$mf == xDataH0$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH0$PropEst, xDataH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rm, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(m[f])," = ",italic(m[m]),", ",italic(s[f])," = ",italic(s[m]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.05,  1.075, expression(paste(bold(P))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(italic(Pi[j])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(0.5,  -0.3,   expression(paste("Recombination rate ",italic((r[j])))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


    # Panel 5: h_j = 1/2
    aDat   <-  aDataH0.5[aDataH0.5$rf == aDataH0.5$rm & aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$sf == aDataH0.5$sm,][-c(12,13),]
    aDat2  <-  aDataH0.5[aDataH0.5$rf >  aDataH0.5$rm & aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$sf == aDataH0.5$sm,][-c(12,13),]
    aDat3  <-  aDataH0.5[aDataH0.5$rf <  aDataH0.5$rm & aDataH0.5$mf == aDataH0.5$mm & aDataH0.5$sf == aDataH0.5$sm,][-c(12,13),]
    xDat   <-  xDataH0[xDataH0$sf == xDataH0$sm & xDataH0$mf == xDataH0$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH0.5$PropEst, xDataH0$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
#        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
#        points(PropEst ~ rf, pch = 22, col=1, bg = '2', data=xDat2)
        points(PropEst ~ rm, pch = 23, col=1, bg = 'grey70', data=aDat3)
#        points(PropEst ~ rf, pch = 23, col=1, bg = '2', data=xDat3)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(Q))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste("Recombination rate ",italic((r[j])))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

    # Panel 6: h_j = 1
    aDat   <-  aDataH1[aDataH1$rf == aDataH1$rm & aDataH1$mf == aDataH1$mm & aDataH1$sf == aDataH1$sm,][-c(12,13),]
    aDat2  <-  aDataH1[aDataH1$rf >  aDataH1$rm & aDataH1$mf == aDataH1$mm & aDataH1$sf == aDataH1$sm,][-c(12,13),]
    aDat3  <-  aDataH1[aDataH1$rf <  aDataH1$rm & aDataH1$mf == aDataH1$mm & aDataH1$sf == aDataH1$sm,][-c(12,13),]
    xDat   <-  xDataH1[xDataH1$sf == xDataH1$sm & xDataH1$mf == xDataH1$mm,]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,max(aDataH1$PropEst, xDataH1$PropEst)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Autosomal frequencies
        abline(h=(2*aDat$mf[1]), lwd=3, lty=2)
        points(PropEst ~ rf, pch = 21, col=1, bg = 'grey70', data=aDat)
        points(PropEst ~ rf, pch = 22, col=1, bg = 'grey70', data=aDat2)
        points(PropEst ~ rm, pch = 23, col=1, bg = 'grey70', data=aDat3)
        points(PropEst ~ rf, pch = 21, col=1, bg = '2', data=xDat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(R))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste("Recombination rate ",italic((r[j])))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
              x       =  usr[2],
              y       =  usr[4]*0.42,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              pch     =  c(21,21),
              pt.bg   =  c('grey70',2),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.24,
              legend  =  c(
                          expression(paste(italic(r[f])," = ",italic(r[m]))),
                          expression(paste(italic(r[f])," > ",italic(r[m]))),
                          expression(paste(italic(r[f])," < ",italic(r[m])))),
              pch     =  c(21,22,23),
              col     =  c(1),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )

} 
