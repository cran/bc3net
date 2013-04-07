\name{bc3net}
\alias{bc3net}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bc3net gene regulatory network inference}

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
Altay (2010a), which we present in the following in a modified form to
obtain a more efficient implementation. Briefly, C3NET consists of three
main steps. First, mutual information values among all gene pairs are
estimated. Second, an extremal selection strategy is applied allowing
each of the p genes in a given dataset to contribute at most one edge to
the inferred network. That means we need to test only p different
hypotheses and not p(p-1)/2. This potential edge corresponds to the
hypothesis test that needs to be conducted for each of the p
genes. Third, a multiple testing procedure is applied to control the
type one error. In the above described context, this results in a
network G^b_k.

}

\usage{
bc3net(dataset, boot=100, estimator="pearson", disc="equalwidth", mtc1=TRUE,
alpha1=0.05, nullit=NA, null=c(), adj1="bonferroni", mtc2=TRUE,
alpha2=0.05, adj2="bonferroni",
weighted=TRUE, igraph=TRUE, verbose=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      gene expression dataset where rows define genes and columns samples
  }
  \item{boot}{
      default 100 bootstrap datasets are generated to infer an ensemble
      of c3net gene regulatory networks
  }
  \item{estimator}{
    estimators for continuous variables
            "pearson", "spearman", "kendall", "spearman"

    estimators for discrete variables
            "emp", "mm","sg","shrink"
  }
  \item{disc}{
	    required for discrete estimators, method for discretize function (see infotheo package)
            "equalwidth" (default),  "equalfreq", "globalequalwidth" 
  }
  \item{nullit}{
     nullit defines the size of the generated null distribution vector
     used for hypothesis testing of significant edges inferred by
     c3net. The null distribution of mutual information is generated from
     sample and gene label randomization.
     
     number of iterations, where the default is defined by

     nullit=ceiling(10^5/(((genes*genes)/2)-genes))

     genes: number of genes 
   }
   \item{null}{
     assign alternatively an external null distribution vector
   }
  \item{mtc1}{
    consider multiple hypothesis testing for edges inferred by c3net
  }
  \item{alpha1}{ significance level for mtc1 }
  \item{adj1}{ if mtc1==TRUE default multiple hypothesis
     testing procedure for c3net inferred edges using "bonferroni"
     (default)

     alternatively use "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr",
     "none" (see ?p.adjust())
  }

  \item{mtc2}{
    Consider multiple hypothesis testing for edges inferred by bc3net. A
    binomial test is performed for each gene pair with an ensemble
    consensus rate >0 consider multiple hypothesis testing for edges
    inferred by bc3net
  }
  \item{alpha2}{
    significance level for mtc2
  }
  \item{adj2}{
    Consider multiple hypothesis testing for edges inferred by bc3net.
    if mtc2==TRUE  "bonferroni" is used as multiple hypothesis testing procedure.
    alternatively use "holm", "hochberg", "hommel", "bonferroni", "BH",
    "BY","fdr", "none"}
  \item{weighted}{ A weighted network is returned, where the
    weights denote the ensemble consensus rate of bc3net.}
  \item{igraph}{ A bc3net igraph object is returned.}
  \item{verbose}{ Return processing information of running procedures.}
}
\details{
BC3NET Gene regulatory network inference
}
\value{
  'bc3net' returns a gene regulatory network formated as adjacency
  matrix, as weighted matrix where the edge weights are defined by the
  corresponding mutual information values or as undirected weighted or
  unweighted igraph object.

}
\references{
  Altay G, Emmert-Streib F. Inferring the conservative causal core of gene
regulatory networks. BMC Syst Biol. 2010 Sep 28;4:132.

 de Matos Simoes R, Emmert-Streib F. Bagging statistical network inference from
large-scale gene expression data. PLoS One. 2012;7(3):e33624. Epub 2012 Mar 30.

 de Matos Simoes R, Emmert-Streib F. Influence of statistical estimators of
mutual information and data heterogeneity on the inference of gene regulatory
networks. PLoS One. 2011;6(12):e29279. Epub 2011 Dec 29.

}

\author{
de Matos Simoes R, Emmert-Streib F.
}

%\note{

%}

\seealso{
C3NET c3mtc 
}

\examples{
 data(expmat)
 bnet=bc3net(expmat)
}


\keyword{bc3net}
\keyword{c3mtc}
