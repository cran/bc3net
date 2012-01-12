\name{makenull}
\alias{makenull}

\title{
Generate a mutual information null distribution from a gene expression matrix

}
\description{

In order to determine the statistical significance of the mutual information values between genes we test for each pair of genes the following null hypothesis.

H_0^I: The mutual information between gene i and j is zero.

Because we are using a nonparametric test we need to obtain the
 corresponding null distribution for H_0^I from a randomization of the
 data. 

'makenull' performs a randomization with the formulated null hypothesis
by permuting the sample and gene labels for all genes of the entire
expression matrix at once. The vector of the mutual information null
distribution is obtained from repeated randomizations for a given
number of iterations. 

}
\usage{
makenull(dataset, nullit=NA, estimator="pearson", disc="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
gene expression data matrix   
}
  \item{nullit}{

  nullit defines the size of the generated null distribution vector
  used for hypothesis testing of significant edges inferred by
  c3net. The null distribution of mutual information is generated from
  sample and gene label randomization.  

  default number of iterations:
  nullit=ceiling(10^5/(((genes*genes)/2)-genes))
  genes: number of genes

}
\item{estimator}{
  
    minet package (continuous estimators)
            "pearson", "spearman", "kendall", "spearman"

    minet package (discrete estimators)
            "mi.empirical", "mi.mm","mi.sg","mi.shrink"

  }
  
  \item{disc}{
	    only required for discrete estimators (minet package)
            "equalfreq", "equalwidth"
	  }
	  

}
%\details{
%Mutual information null distribution
%}

\value{
  The function 'makenull' generates a numeric vector of the mutual information
  null distribution estimated from a given number of randomized
  datasets.
}
\references{
  de Matos Simoes R, Emmert-Streib F., Bagging statistical network
  inference from large-scale gene expression data., PLoS
  One. 2012;7(3):e33624. Epub 2012 Mar 30.

  
}
\author{
Ricardo de Matos Simoes  <r.dematossimoes@qub.ac.uk>
Frank Emmert-Streib <f.emmert-streib@qub.ac.uk>
}

%\note{

%}


\seealso{
c3mtc bc3net
}
\examples{
 data(expmat)
 null=makenull(expmat, nullit=5)
}

\keyword{internal}
