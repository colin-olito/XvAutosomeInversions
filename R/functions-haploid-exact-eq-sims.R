## Dependencies ##
#source('./R/functions-plots.R')


#' Simulations to validate weak selection/migration approximations for haploid two-locus model
#' of migration-seleciton balance
#' 
#' @title Simulations to validate weak selection/migration approximations for haploid two-locus model
#' @param m Migration rate
#' @param m Recombination rate
#' @param generations Timescale for the simulations
#' @param initFreqs Initial genotypic frequencies - must be a numeric vector of length 4 that sums to 1
#' @return Returns a data frame with columns: "m", "r", "s", "approx.a", "freq.a", "freq.b", "LD"
#' @seealso 
#' @export
#' @author Colin Olito & Tim Connallon
#'
hapApproxCompare  <-  function(m = 0.01, r = 0.5, generations = 10000, initFreqs = rep(0.25,4)) {

  # Parameters
  sels  <-  seq(from = 0.01, to = 0.2, by = 0.005)
  
  # Data frame for results
  res  <-  data.frame("m" = rep(m,times=length(sels)), 
                      "r" = rep(r,times=length(sels)), 
                      "s" = sels, 
                      "approx.a" = rep(0, times=length(sels)), 
                      "freq.a" = rep(0, times=length(sels)), 
                      "freq.b" = rep(0, times=length(sels)), 
                      "LD" = rep(0, times=length(sels)))
  
  # Simulation Loop
  for(j in 1:length(sels)) {
  
    #parameters
    s = sels[j]
  
    #initial haplotype frequencies
    x1 = 0.25
    x2 = 0.25
    x3 = 0.25
    x4 = 0.25
  
    for(i in 1:generations){
      # Tim's recursions
      W = x1 + x2*(1 + s) + x3*(1 + s) + x4*(1 + s)^2
      D.adults = ((1 - m)*(1 + s)^2)*((x1*x4 - x2*x3)*(1 - m)/W + m*x4)/W
      #frequency next generation
      y1 = x1*(1 - m)/W + m - r*D.adults
      y2 = x2*(1 + s)*(1 - m)/W + r*D.adults
      y3 = x3*(1 + s)*(1 - m)/W + r*D.adults
      y4 = x4*((1 - m)*(1 + s)^2)/W - r*D.adults
      x1 = y1
      x2 = y2
      x3 = y3
      x4 = y4
    }
  
  #approx.a = m/s
  res$approx.a[j] =  ((s + m*(1 + s)) - sqrt((s - m*(1 + s))^2))/(2*s)
  res$freq.a[j]   =  x1 + x3
  res$freq.b[j]   =  x1 + x2
  res$LD[j]       =  (x1*x4 - x2*x3)
  }

  # Return results data frame
  res
}




#' Plotting functions for summary figures for simulations to validate weak selection/migration approximations for haploid two-locus model
#' of migration-seleciton balance
#' 
#' @title Plottinf functions for simulations to validate weak selection/migration approximations for haploid two-locus model
#' @param m Migration rate
#' @param generations Timescale for the simulations
#' @param initFreqs Initial genotypic frequencies - must be a numeric vector of length 4 that sums to 1
#' @param makePlots Switch - should simulation automatically make summary plots?
#' @return Returns a data frame with columns sapprox.a
#' @seealso 
#' @export
#' @author Colin Olito & Tim Connallon
#'
hapApproxComparePlots  <-  function(data=res) {

    # Recover parameters
    m     <-  res$m[1]
    r     <-  res$r[1]
    sels  <-  res$s

    # Plot Layout
    par(mfrow=c(1,2), omi=rep(0.5,4), mai = rep(0.75,4))

    # Plot allele frequencies 
      # Plot of allele frequencies
      plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(data$s)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
      usr  <-  par('usr')
      rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
      plotGrid(lineCol='grey80')
      box()
      # Plot points/lines
      points(freq.a  ~ sels, pch=21, col='dodgerblue4', bg=transparentColor('dodgerblue', opacity=0.7), data=res)
      lines(approx.a ~ sels, lwd=2, data=res)
      # axes
      axis(1, las=1)
      axis(2, las=1)
        proportionalLabel(0.24, 1.1, expression(paste(italic(m)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.425, 1.11, substitute(m,list(m=rounded(m,precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.67, 1.1, expression(paste(italic(r)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.825, 1.11, substitute(r,list(r=rounded(r,precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
      proportionalLabel(0.5, -0.25, expression(paste(italic(s))), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
      proportionalLabel(-0.25, 0.5, expression(paste("Allele Frequency")), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=90)
        legend(
            x       =  usr[2]*0.99,
            y       =  usr[4],
            legend  =  c(
                        expression(paste(hat(italic(p))," approx.")),
                        expression(paste("Exact"))),
            lty     =  c(1, NA),
            lwd     =  c(2,1),
            pch     =  c(NA, 21),
            pt.bg   =  c(NA, 'dodgerblue'),
            col     =  c('black','dodgerblue4'),
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)


      # Plot of LD
      plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(sels)), ylim = c(0,max(res$LD)*1.05), ylab='', xlab='', cex.lab=1.2)
      usr  <-  par('usr')
      rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
      plotGrid(lineCol='grey80')
      box()
      # Plot points/lines
      lines(LD ~ sels, lwd=3, col='dodgerblue2', data=res)
      # axes
      axis(1, las=1)
      axis(2, las=1)
        proportionalLabel(0.24, 1.1, expression(paste(italic(m)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.425, 1.11, substitute(m,list(m=rounded(m,precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.67, 1.1, expression(paste(italic(r)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.825, 1.11, substitute(r,list(r=rounded(r,precision=2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
      proportionalLabel(0.5, -0.25, expression(paste(italic(s))), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
      proportionalLabel(-0.3, 0.5, expression(paste(italic(LD))), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=90)
}
