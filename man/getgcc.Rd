\name{getgcc}
\alias{getgcc}

\title{Extracting the giant connected component from an igraph object}

\description{
 The function uses igraph::clusters and igraph::induced.subgraph to extract
the giant component of an igraph object.
 }

\usage{
  getgcc(net)
}

\arguments{
   \item{net}{
      An igraph object.
    }	  
}

\value{
  \code{getgcc} returns an igraph object with the giant connected 
  component of an igraph object.
}
 
\details{
     In a connected component of an undirected graph all pairs of vertices u and v are reachable by a path. The giant connected component describes the largest connected component of a graph object.
 }

\references{
    G. Csardi and T. Nepusz, The igraph software package for complex network research,
    InterJournal, Complex Systems,1695, 2006, http://igraph.org  
}


\examples{
# The example network "exanet" contains 814 genes with 1103 edges
data(exanet)

# The giant connected component of the example 
# network contains 382 genes with 989 edges
getgcc(exanet)

}

\keyword{giant connected component}




