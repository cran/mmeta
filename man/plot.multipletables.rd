\name{plot.multipletables}
\alias{plot.multipletables}
\title{Plot Method for \code{multipletables} objects}
\description{
   Produces a variety of plots for multiple tables analysis
}
\usage{
    \method{plot}{multipletables}(x,type=NULL,select=NULL,file=NULL, xlim=NULL,ylim=NULL,
                                xlabel=NULL,mar=c(5, 10, 4, 9),xlog=TRUE,
                                addline=NULL,xlab=NULL,ylab=NULL,ciShow=TRUE,...)
}

\arguments{
   \item{x}{an object inheriting from class \code{multipletables}.}
   \item{type}{a chracter string specifying the type of plots to
     produce. Options are \code{sidebyside}, \code{overlap}, and \code{forest}. See details}
   \item{select}{a numeric value or vector specifying which studies to
     be plotted. By default (when \code{NULL}), all of the studies will be plotted.}  
   \item{xlab}{a character string specifying the x-axis label in the
     plot. Default is the name of the measure of association}
   \item{ylab}{a character string specifying the x-axis label in the plot. Default is "Density"}
   \item{file}{a character string specifying the filename as which the plots
     are saved. By default (when \code{NULL}), the plots are displayed on screen. See details.}
   \item{xlim, ylim}{a numeric vectors of length 2 specifying the lower
     and upper limits of the axes. By default (when \code{NULL}), \code{xlim} and
     \code{ylim} are computed. For forest plots, if the lower bound of
     any measure is smaller than \code{xlim[1]} or the upper bound of
     any measure is larger than \code{xlim[2]}, arrows will be used at
     the limits to denote the bounds exceed the specified ranges.}
   \item{xlabel}{a numeric vector specifying at which tick-marks are to
     be drawn. By default (when \code{NULL}), tickmark locations are
     computed.}
   \item{addline}{a numeric value specifying the x-value for a vertical
     reference line at \code{x=addline}. Default is \code{NULL}}
   \item{xlog}{a logical value indicating whether a logarithmic scale
     should be used for x-axis. Default is \code{TRUE} for measures
     \code{OR} and \code{RR} and \code{FALSE} for measure \code{RD}.} 
   \item{mar}{A numerical vector of 4 values which control the space (in the number of lines)
     between the axes and the border of the graph of the form
     \code{c(bottom, left, top, right)} the default values are \code{c(5.1, 4.1, 4.1, 2.1)}.}
    \item{ciShow}{a logical value; if \code{TRUE} (default), the true
      credible intervals numbers will display at the right hand side of
      the forest plot.}
    \item{...}{Other arguments can be passed to plot function}   
}

\details{   
     If \code{type="sidebyside"}, the posterior distributions of all
     study-specific measure are displayed side by side in 4-panel plots
     with study names.
     
     If \code{type="overlap"}, the posterior distributions of all
     study-specific measure are displayed in one graph. To clarity, it
     is advisable to specify a few studies by \code{select} argument. 

     If \code{type="forest")}, a forest plot of all study-specific and
       overall measure with 95\% credible/confidence intervals are
       plotted.
       
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

\author{Xiao Su <Xiao.Su@uth.tmc.edu>}

\seealso{\code{\link{multipletables}}
         \code{\link{summary.multipletables}}
         \code{\link{xtable.multipletables}}
       }

\examples{
#library(mmeta)
# Analyze the dataset colorectal to conduct exact inference of the odds ratios
#data(colorectal)
#multiple.OR <- multipletables(data=colorectal, measure="OR", model="Sarmanov", method="exact")
# Generate the forest plot with 95\% CIs of study-specific odds ratios
#and 95\% CI of overall odds ratio
#plot(multiple.OR, type="forest", addline=1)
# Plot the posterior density functions of some target studies in an overlaying manner
#plot(multiple.OR, type="overlap", select=c(4,14,16,20))
# Plot the posterior density functions of some target studies in a
#side-by-side manner 
#plot(multiple.OR, type="sidebyside", select=c(4,14,16,20), ylim=c(0,2.7), xlim=c(0.5,1.5))


# Analyze the dataset withdrawal to conduct inference of the relative risks
#data(withdrawal)
#multiple.RR <- multipletables(data=withdrawal, measure="RR",model="Sarmanov")
#plot(multiple.RR, type="forest", addline=1)
#plot(multiple.RR, type="overlap", select=c(3,8,14,16))
#plot(multiple.RR, type="sidebyside", select=c(3,8,14,16), ylim=c(0,1.2),
#xlim=c(0.4,3))

# Analyze the dataset withdrawal to conduct inference of the risk differences
#data(withdrawal)
#multiple.RD <- multipletables(data=withdrawal, measure="RD", model="Sarmanov")
#summary(multiple.RD)
#plot(multiple.RD, type="forest", addline=0)
#plot(multiple.RD, type="overlap", select=c(3,8,14,16))
#plot(multiple.RD, type="sidebyside", select=c(3,8,14,16))
#plot(multiple.RD, type="sidebyside", select=c(3,8,14,16),
#     ylim=c(0,6), xlim=c(-0.2,0.4))
}


\keyword{methods}