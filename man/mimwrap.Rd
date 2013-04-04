\name{mimwrap}
\alias{mimwrap}

\title{Wrapper function for mutual information matrix estimators}

\description{
  Mutual information matrix estimation wrapper function for various
  mutual information estimators. Depends on infotheo package for mutual information estimators on discrete variables.
 }

\usage{
  mimwrap(dataset, estimator="pearson", disc="equalwidth", bins = sqrt(ncol(dataset)))
}

\arguments{
   \item{dataset}{
      Data gene expression matrix where rows denote genes (features) and
      columns samples.
    }  

  \item{estimator}{
    estimators for continuous variables
            "pearson" (default), "spearman", "kendall", "spearman"

    estimators for discrete variables (infotheo package)
            "emp", "mm","sg","shrink"

  }  
  \item{disc}{
	    only required for discrete estimators (see infotheo package)
            "equalwidth" (default), "globalequalwidth" , "equalfreq"
	  }
  \item{bins}{
            number of bins for the descretize function (infotheo), 
            default sqrt(ncol(dataset))
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




