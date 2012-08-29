\name{bspline}
\alias{bspline}


\title{
  B-spline estimator wrapper function for 'mis_calc'
}

\description{
  This is a wrapper function for the B-spline estimator program
  "mis_calc" (Daub 2004) obtained from Carsten Daub (oncaphillis@snafu.de).
}

\usage{
  bspline(expmat)
}

\arguments{
  \item{expmat}{
     gene expression matrix, where rows define genes and columns samples
   }
   
}

\details{
Mutual information estimator
}

\value{
The function 'bspline' returns a mutual information matrix, where rows
and columns are defined by genes.

}

\references{
Daub CO, Steuer R, Selbig J, Kloska S., Estimating mutual information
using B-spline functions--an improved similarity measure for analysing
gene expression data. BMC Bioinformatics. 2004 Aug 31;5:118.

de Matos Simoes R, Emmert-Streib F. Bagging statistical network inference from
large-scale gene expression data. PLoS One. 2012;7(3):e33624. Epub 2012 Mar 30.

}
\author{
de Matos Simoes R, Emmert-Streib F.
}

\note{
 The function requires the use of a linux system and the installation of "mis_calc" obtained from Carsten Daub (oncaphillis@snafu.de).
}

%\seealso{

%}

\examples{
data(expmat)
mim=bspline(expmat)

}


%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }