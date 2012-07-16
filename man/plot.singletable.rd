\name{plot.singletable}
\alias{plot.singletable}
\title{Plot Method for \code{singletable} objects}
\description{
   Produces various plots for single table analysis.
}
\usage{
    \method{plot}{singletable}(x,type=type,file=NULL,select=c(1,2),xlab=NULL,ylab=NULL,addline=NULL,xlim=NULL,ylim=NULL,...)
}

\arguments{  
   \item{x}{an object inheriting from class \code{singletable}.}
   \item{type}{a chracter string specifying the type of plots to
     produce. Options are \code{sidebyside} and \code{overlap}. See details}
   \item{file}{a character string specifying the filename as which the plots
     are saved. Default is NULL to view on screen. See details.}
   \item{select}{a numeric value or vector specifying which distribution
     should be plotted. \code{select=1} is posterior distribution. \code{select=2} is
     prior distribution. \code{select=c(1,2)} is both. Default is
     c(1,2). This argument is only used when \code{type="sidebyside"}.}
   \item{xlab}{a character string specifying the x-axis label in the
     plot. Default is the name of the measure of association}
   \item{ylab}{a character string specifying the x-axis label in the plot. Default is "Density"}
   \item{addline}{a numeric value specifying the x-value for a vertical
     reference line at \code{x=addline}. Default is NULL}
   \item{xlim, ylim}{a numeric vectors of length 2 specifying the lower
     and upper limits of the axes}
   \item{...}{Other arguments can be passed to plot function} 
 }
 
\details{   
     If \code{type="sidebyside"}, the posterior distribution of measure
     and the prior distribution are drawn side by side in two plots. If
     \code{type="overlap"}, the posterior distribution of measure and
     the prior distribution are overlaid in one plot. 
     
     If \code{file=NULL}, the plots will be displayed on screen. Or
     else, the plots will be saved as "./mmeta/code{file}.eps", where
     "./" denotes current working directory.
}


\references{
Chen, Y., Luo, S., (2011a). A few remarks on 'Statistical distribution of
the difference of two proportions'. Statistics in Medicine 30, 1913-1915.

Chen, Y., Chu, H., Luo, S., Nie, L., and Chen, S. (2011b). Bayesian
analysis on meta-analysis of case-control studies accounting for
within-study correlation. Statistical Methods in Medical Research,
Published online on Dec 4, 2011, PMID: 22143403. doi: 10.1177/0962280211430889

Chen, Y., Luo, S., Chu, H., Su, X., and Nie, L. (2012a). An empirical
Bayes method for multivariate meta-analysis with an application in
clinical trials. in press at Communication in Statistics:
Theory and Methods.
}



\seealso{\code{\link{singletable}}}
\examples{
# Inference under Jeffreys prior distribution
#single.OR.Jeffreys <- singletable(a1=0.5, b1=0.5, a2=0.5,
#                                  b2=0.5, y1=40, n1=96, y2=49, n2=109,
#                                  model="Independent",
#                                  measure="OR", method="exact")
#summary(single.OR.Jeffreys)

# Inference under Laplace prior distribution
#single.OR.Laplace <- singletable(a1=1, b1=1, a2=1, b2=1,
#                                 y1=40, n1=96, y2=49, n2=109,
#                                 model="Independent", measure="OR",
#                                 method="exact")
# Inference under Sarmanov prior distribution with positive correlation
#single.OR.Sar1 <- singletable(a1=0.5, b1=0.5, a2=0.5, b2=0.5,
#                              rho=0.5, y1=40, n1=96, y2=49, n2=109,
#                              model="Sarmanov",
#                              measure="OR", method="exact")
# Inference under Sarmanov prior distribution with negative correlation
#single.OR.Sar2 <- singletable(a1=0.5, b1=0.5, a2=0.5, b2=0.5,
#                              rho=-0.5, y1=40, n1=96, y2=49, n2=109,
#                              model="Sarmanov",
#                              measure="OR", method="exact")
# generate a 2X2 panel plot
#par(mfrow=c(2,2))
#plot(single.OR.Jeffreys, type="overlap", xlim=c(0.5, 2),
#    main="Jefferys Prior")
#plot(single.OR.Laplace, type="overlap", xlim=c(0.5, 2),
 #    main="Laplace Prior")
#plot(single.OR.Sar1, type="overlap", xlim=c(0.5, 2),
#     main=expression(paste("Sarmanov Prior ",rho," = 0.5")))
#plot(single.OR.Sar2, type="overlap", xlim=c(0.5, 2),
#     main=expression(paste("Sarmanov Prior ",rho," = -0.5")))
}


\keyword{singletable}