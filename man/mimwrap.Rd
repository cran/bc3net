\name{mimwrap}
\alias{mimwrap}

\title{Wrapper function for mutual information matrix estimators}

\description{
  Mutual information matrix estimation wrapper function for various
  mutual information estimators (minet package, b-spline).
 }

\usage{
  mimwrap(dataset, estimator="pearson", disc="none")
}

\arguments{
   \item{dataset}{
      Data gene expression matrix where rows denote genes (features) and
      columns samples.
    }  

  \item{estimator}{
    minet package (continuous estimators)
            "pearson", "spearman", "kendall", "spearman"

    minet package (discrete estimators)
            "mi.empirical", "mi.mm","mi.sg","mi.shrink"

    c3net gaussian estimator (pearson)
            "gaussian"

  }  
  \item{disc}{
	    only required for discrete estimators (see minet package)
            "equalwidth" (default), "globalequalwidth" , "equalfreq"
	  }
	  
}

\value{
  \code{mimwrap} returns a symmetric mutual information matrix for
  various mutual information estimators.
}
 
\details{
      A mutual information matrix is estimated from a gene expression data set
    }

\references{
  % add minet/misc

  Patrick E Meyer, Frederic Lafitte and Gianluca Bontempi, minet: A
  R/Bioconductor Package for Inferring Large Transcriptional Networks
  Using Mutual Information, BMC Bioinformatics 2008, 9:461
  
  Carsten O. Daub, Ralf Steuer, Joachim Selbig, and Sebastian Kloska,
  Estimating mutual information using B-spline functions - an improved similarity measure for analysing gene expression data, BMC Bioinformatics. 2004; 5: 118

  de Matos Simoes R, Emmert-Streib F., Bagging statistical network
  inference from large-scale gene expression data., PLoS
  One. 2012;7(3):e33624. Epub 2012 Mar 30.
  
}

%\seealso{\code{\link{c3}}, \code{\link{copula}} }

\examples{
 
data(expmat)
mim <- mimwrap(expmat)

}

\keyword{misc}




