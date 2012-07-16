\name{xtable.multipletables}
\alias{xtable.multipletables}
\title{Create export tables in Latex code of objects \code{multipletables}}
\description{
    Function converting an object inheriting from class
    \code{multipletables} to an \code{xtable} object, which can then be
    printed as LaTeX table.
  }
  
\usage{
    \method{xtable}{multipletables}(x,caption = NULL, label = NULL, align = NULL,
    digits = NULL, display = NULL,...)
}

\arguments{
  \item{x}{an object inheriting from class \code{multipletables}.}
  \item{caption}{the \code{caption} argument of \code{xtable()}.}
  \item{label}{the \code{label} argument of \code{xtable()}.}
  \item{align}{the \code{align} argument of \code{xtable()}.}
  \item{digits}{the \code{digits} argument of \code{xtable()}.}
  \item{display}{the \code{display} argument of \code{xtable()}.}
   \item{...}{additional arguments; currently none is used.}
 }


\value{
  An object of class \code{xtable} which inherits the \code{data.frame} class and contains 
  several additional attributes specifying the table formatting options.
}


\note{
  This latex table generated in \code{R 2.13.0} by \code{xtable 1.5-6} package.
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



\seealso{\code{\link{multipletables}}
         \code{\link{plot.multipletables}}
         \code{\link{summary.multipletables}}
	 \code{\link{xtable}}
       }

       
\examples{
#library(mmeta)

# Analyze the dataset colorectal to conduct exact inference of the odds ratios
#data(colorectal)
#multiple.OR <- multipletables(data=colorectal, measure="OR", model="Sarmanov", method="exact")
## Convert R object to an xtable object to be printed as a Latex or HTML table
#multiple.OR.table <- xtable(multiple.OR)
#print(multiple.OR.table)
#print(multiple.OR.table, type="html")

# Analyze the dataset withdrawal to conduct inference of the relative risks
#data(withdrawal)
#multiple.RR <- multipletables(data=withdrawal, measure="RR",
#                              model="Sarmanov")
#multiple.RR.table <- xtable(multiple.RR)
#print(multiple.RR.table)
#print(multiple.RR.table, type="html")

# Analyze the dataset withdrawal to conduct inference of the risk differences
#data(withdrawal)
#multiple.RD <- multipletables(data=withdrawal, measure="RD",
#                              model="Sarmanov")
#multiple.RD.table <- xtable(multiple.RD)
#print(multiple.RD.table)
#print(multiple.RD.table, type="html")
}


\keyword{ latex }