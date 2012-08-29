\name{bc3net-package}

\alias{bc3net-package}

\docType{package}

\title{
BC3NET Gene Regulatory network Inference
}

\description{
The basic idea of BC3NET is to generate from one dataset D_s,
consisting of s samples, an ensemble of B independent bootstrap datasets
D_k by sampling from D(s) with replacement by using a non-parametric
bootstrap (Efron 1993). Then, for each generated data set D_k in the
ensemble, a network G^b_k is inferred by using C3NET (Altay 2010a). From the ensemble of networks G^b_k we construct one weighted network 
G^b_w which is used to determine the statistical significance of the
connection between gene pairs. This results in the final binary,
undirected network G.

A base component of BC3NET is the inference method C3NET introduced in
Altay (2010a), which we present in the following in a modified form to obtain a more efficient implementation. Briefly, C3NET consists of three main steps. First, mutual information values among all gene pairs are estimated. Second, an extremal selection strategy is applied allowing each of the p genes in a given dataset to contribute at most one edge to the inferred network. That means we need to test only p different hypotheses and not p(p-1)/2. This potential edge corresponds to the hypothesis test that needs to be conducted for each of the p genes. Third, a multiple testing procedure is applied to control the type one error. In the above described context, this results in a network G^b_k.
  
}
\details{
\tabular{ll}{
Package: \tab bc3net\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2012-01-12\cr
License: \tab GPL (>=2)\cr
}

 
  bc3net.R c3mtc.R makenull.R mimwrap.R getpval.R mat2igraph.R

}
\author{
Ricardo de Matos Simoes  <r.dematossimoes@qub.ac.uk>
Frank Emmert-Streib <f.emmert-streib@qub.ac.uk>
Maintainer: Ricardo de Matos Simoes <r.dematossimoes@qub.ac.uk>
}

\references{
  de Matos Simoes R, Emmert-Streib F., Bagging statistical network
  inference from large-scale gene expression data., PLoS One. 2012;7(3):e33624. Epub 2012 Mar 30.
}


\keyword{ package }
\seealso{
C3NET, MINET, INFOTHEO
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

data(expmat)
bnet=bc3net(expmat)

data(expmat)
cnet=c3mtc(expmat)

}
