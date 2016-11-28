\name{gpea} 

\alias{gpea} 


\title{ 

Gene pair enrichment analysis (GPEA)

} 

\description{ 

When a network G contains n interactions, of which k interactions are among genes from the given gene set S, then a p-value for the enrichment of gene pairs of this gene set S can be calculated based on a e.g., one-sided Fisher's exact test. For p genes there is a total of N=p(p-1)/2 different gene pairs (clique graph) with the assumption that all genes within a gene set are associated to each other. If there are pS genes for a particular gene set (S) then the total number of gene pairs for this gene set is mS=pS(pS-1)/2. 
} 

\usage{ 


gpea(gnet, genesets, verbose = TRUE, cmax = 1000, cmin = 3,
                 adj = "bonferroni")
} 


\arguments{ 

  \item{gnet}{ 
 igraph object (e.g., inferred from bc3net) of a given network where the gene identifiers [V(net)$names] correspond to the provided gene identifiers in the reference gene sets.

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

  \item{verbose}{ 
    The default value is <FALSE>. If this option is set <TRUE> the number and name of the 
    gene sets during their processing is reported.
  } 


  \item{cmax}{ 
   All provided genesets with more than cmax genes will be excluded from the analysis (default cmax=1000).
  } 

  \item{cmin}{ 
   All provided genesets with less than cmin genes will be excluded from the analysis
   (default cmin>=3).
} 
  \item{adj}{ 
    The default value is <fdr> (False discovery rate using the Benjamini-Hochberg approach). 
    Multiple testing correction based on the function stats::p.adjust() with available options for 
    "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" and "none"
 } 
} 

\details{ 
  The enrichment analysis is based on a one-sided Fisher's exact test.
} 

\value{ 
The function returns a data.frame object with the columns 

"TermID" the given name of a gene set i from the named gene set collection list object S.
"edges"  the number of connected gene pairs present a given geneset
"genes"  the number of candidate genes present in the gene set i
"all"    the number of all genes present in the gene set i
"pval"   the nominal p-value from a one-sided fisher's exact test
"padj"   the adjusted p-value to consider for multiple testing

} 

\references{ 

 Inference and Analysis of Gene Regulatory Networks in R: Applications in Biology, Medicine, and Chemistry, DOI: 10.1002/9783527694365.ch10
In book: Computational Network Analysis with R, 2016, pp.289-306

Urothelial cancer gene regulatory networks inferred from large-scale RNAseq, Bead and Oligo gene expression data, BMC Syst Biol. 2015; 9: 21. 

} 

\author{ 

Ricardo de Matos Simoes

} 

\seealso{ 

See Also as \code{\link{enrichment}}

} 

\examples{

data(exanet)
data(exgensets) ## example gene sets from the CPDB database (http://www.consensuspathdb.org)

res = gpea(exanet, exgensets, cmax=1000, cmin=2)


}



