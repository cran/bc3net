\name{c3mtc}
\alias{c3mtc}

\title{
'c3mtc' gene regulatory network inference using c3net with multiple testing correction procedure
}

\description{
We present in the following the inference method C3NET introduced in
Altay (2010a) in a modified form to obtain a more efficient
implementation. Briefly, C3NET consists of three main steps. First,
mutual information values among all gene pairs are estimated. Second, an
extremal selection strategy is applied allowing each of the p genes in a
given dataset to contribute at most one edge to the inferred
network. That means we need to test only p different hypotheses and not
p(p-1)/2. This potential edge corresponds to the hypothesis test that
needs to be conducted for each of the p genes. Third, a multiple testing
procedure is applied to control the type one error.

In order to determine the statistical significance of the mutual information values between genes we test for each pair of genes the following null hypothesis.

H_0^I: The mutual information between gene i and j is zero.

Because we are using a nonparametric test we need to obtain the
 corresponding null distribution for H_0^I from a randomization of the
 data. 

The formulated null hypothesis is performed by permuting the sample and
gene labels for all genes of the entire expression matrix at once. The
vector of the mutual information null distribution is obtained from
repeated randomizations for a given number of iterations. 

}

\usage{

c3mtc(dataset, null=NULL, mtc=TRUE, adj="bonferroni", alpha=0.05, nullit=NA,
estimator="pearson", disc="none", adjacency=FALSE, igraph=TRUE)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{dataset}{
      gene expression dataset where rows define genes and columns samples
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

      c3net gaussian estimator (pearson)
            "gaussian"

      bspline requires installation of "mis_calc"
            "bspline"}
  
   \item{disc}{
	    only required for discrete estimators (minet package)
            "equalfreq", "equalwidth"
   }
   \item{mtc}{
     consider multiple hypothesis testing for edges inferred by c3net
   }
   \item{adj}{ if mtc==TRUE default multiple hypothesis
     testing procedure for c3net inferred edges using "bonferroni"
     (default)

     alternatively use "holm", "hochberg", "hommel", "bonferroni", "BH",
   "BY","fdr", "none" (see ?p.adjust())
   }
   \item{alpha}{
     significance level for mtc after multiple hypothesis
     testing correction
   }
   \item{adjacency}{
     return an adjacency matrix
   }
   \item{igraph}{
     return igraph object
   }
   \item{null}{
     If NULL a null distribution vector is generated from a
     sample label and gene label permutation of the gene expression
     matrix. For the ensemble inference of one dataset an external null
     distribution vector is suggested for decreasing running time.
   }
}

%\details{
  
%}

\value{
  'c3mtc' returns a gene regulatory network formated as adjacency
  matrix, as weighted matrix where the edge weights are defined by the
  corresponding mutual information values or as undirected weighted or
  unweighted igraph object.
}

\references{
  Altay G, Emmert-Streib F. Inferring the conservative causal core of gene
regulatory networks. BMC Syst Biol. 2010 Sep 28;4:132. PubMed PMID: 20920161;
PubMed Central PMCID: PMC2955605.

 de Matos Simoes R, Emmert-Streib F. Bagging statistical network inference from
large-scale gene expression data. PLoS One. 2012;7(3):e33624. Epub 2012 Mar 30.
PubMed PMID: 22479422; PubMed Central PMCID: PMC3316596.

 de Matos Simoes R, Emmert-Streib F. Influence of statistical estimators of
mutual information and data heterogeneity on the inference of gene regulatory
networks. PLoS One. 2011;6(12):e29279. Epub 2011 Dec 29. PubMed PMID: 22242113;
PubMed Central PMCID: PMC3248437.
}

\author{
  de Matos Simoes R, Emmert-Streib F.
}

% \note{
% }

\seealso{
\code{\link{c3}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
data(expmat)
net=c3mtc(expmat)

}


\keyword{c3net}
\keyword{multiple}
\keyword{hypothesis}
\keyword{testing}
