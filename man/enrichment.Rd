\name{enrichment} 

\alias{enrichment} 

%- Also NEED an '\alias' for EACH other topic documented here. 

\title{ 

Function that performs a functional enrichment analysis based on a one-sided Fisher's exact teset (hypergeometric test).

} 

\description{ 

For a given set of candidate genes, reference genes and a list object of gene sets (e.g. Gene Ontology terms or gene sets from pathways) a one-sided Fisher's exact test ("greater") is performed for each gene set in the collection.

} 

\usage{  
   enrichment(genes, reference, genesets, adj = "fdr", verbose = FALSE)  
} 

\arguments{ 

  \item{genes}{ 
     A character vector of gene identifiers used as a candidate gene list that is assessed by the functional enrichment analysis. The candidate gene list is a subset of the reference gene list.
} 

  \item{reference}{ 

  A character vector of gene identifiers used as the reference gene list. Note all candidate genes are included in the reference gene list.

} 

  \item{genesets}{ 

     A named list object of a collection of gene sets. The identifiers used for the candidate and reference genes need to match the identifier types used for the gene sets. An example list of gene sets is given in data(exgensets) showing an example list object of gene sets from pathways with gene symbols.

$`Reactome:REACT_115566:Cell Cycle`
  [1] "APITD1"           "TAOK1"            "CDC23"           
  [4] "RB1"              "PRKCA"            "HIST1H4J"        
  [7] "MCM10"            "PPP1CC"           "NUP153"          
...

$`Reactome:REACT_152:Cell Cycle, Mitotic`
  [1] "APITD1"           "TAOK1"            "CDKN2C"          
  [4] "RB1"              "PRKCA"            "MCM10"           
  [7] "HIST1H2BH"        "NUP153"           "TUBGCP3"         
 [10] "APEX1"            "RPA2"             "PRKACA"          
...

} 


  \item{adj}{
    The default value is <fdr> (False discovery rate using the Benjamini-Hochberg approach). 
    Multiple testing correction based on the function stats::p.adjust() with available options for 
    "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" and "none"
  } 
 

  \item{verbose}{ 
    The default value is <FALSE>. If this option is set <TRUE> the number and name of the 
    gene sets during their processing is reported.
  } 
 

} 

\details{ 
  The enrichment analysis is based on a one-sided Fisher's exact test.
} 


\value{ 
The function returns a data.frame object with the columns 

"TermID" the given name of a gene set i from the named gene set collection list object S.
"genes"  the number of candidate genes present in the gene set i
"all"    the number of all genes present in the gene set i
"pval"   the nominal p-value from a one-sided fisher's exact test
"padj"   the adjusted p-value to consider for multiple testing

} 

\references{ 

de Matos Simoes R, Tripathi S, Emmert-Streib F. Organizational structure and the periphery of the gene regulatory network in B-cell lymphoma. BMC Syst Biol. 2012; 6:38.

Inference and Analysis of Gene Regulatory Networks in R: Applications in Biology, Medicine, and Chemistry, DOI: 10.1002/9783527694365.ch10
In book: Computational Network Analysis with R, 2016, pp.289-306
} 

\author{ 

Ricardo de Matos Simoes

} 

\examples{

# In the following candidate genes are defined from a
# giant connected component of an example igraph network 
# where we use the remaining genes of a given network 
# as a reference list.

data(exanet)
data(exgensets)

candidate=V(getgcc(exanet))$name
reference=V(exanet)$name

# hypergeometric test is performed for 
# each defined set of genes in the 
# list object exgensets
tab.hypg=enrichment(candidate, reference, exgensets, verbose=TRUE) 

}


