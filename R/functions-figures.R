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