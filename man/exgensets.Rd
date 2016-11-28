\name{exgensets}
\alias{exgensets}
\docType{data}

\title{Example gene sets defined from the CPDB database}

\description{
  This list data object contains 25 gene sets defined from the 
  CPDB database and is used for testing purposes.
}
\usage{
data(exgensets)
}

\format{
  A list object with 25 character vectors.
}

\details{
  An list object of 25 pathway gene sets from CPDB with gene symbols.
}

\source{
http://consensuspathdb.org/
}

\references{
A. Kamburov, U. Stelzl, H. Lehrach and R. Herwig, The ConsensusPathDB interaction database: 2013 update, Nucleic Acids Research, Volume 41, Issue D1, Pp. D793-D800
 }
 
\examples{
  data(exgensets)

  # The gene set collection from CPDB can be formatted by:
  # "CPDB_pathways_genes.tab" is available in the download section in http://consensuspathdb.org/
  # cpdb = readLines("CPDB_pathways_genes.tab", warn = FALSE)
  # cpdb = lapply(cpdb, function(x) strsplit(x, "\t")[[1]])
  # names(cpdb) = sapply(cpdb, function(x) paste(x[3], x[2], x[1], sep=":") )       
  # cpdb = lapply(cpdb, function(x) x[-c(1:3)])
  # cpdb = cpdb[-length(cpdb)]
  # cpdb = lapply(cpdb, function(x) strsplit(x, ",")[[1]])
}


\keyword{datasets}
