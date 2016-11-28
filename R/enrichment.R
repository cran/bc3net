enrichment=function (genes, reference, genesets, adj = "fdr", verbose=FALSE) 
{
    tab = lapply(1:length(genesets), function(i) {
        if(verbose==TRUE){
           cat("processing term", i, names(genesets)[i], "\n")
        }
        
        reference = reference[!reference %in% genes]
        RinSet = sum(reference %in% genesets[[i]])
        RninSet = length(reference) - RinSet
        GinSet = sum(genes %in% genesets[[i]])
        GninSet = length(genes) - GinSet
        fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2, 
            ncol = 2, byrow = F)
        colnames(fmat) = c("inSet", "ninSet")
        rownames(fmat) = c("genes", "reference")
        fish = fisher.test(fmat, alternative = "greater")
        pval = fish$p.value
        inSet = RinSet + GinSet
        res = c(GinSet, inSet, pval)
        res
    })
    rtab = do.call("rbind", tab)
    rtab = data.frame(as.vector(names(genesets)), rtab)
    rtab = rtab[order(rtab[, 4]), ]
    colnames(rtab) = c("TermID", "genes", "all", "pval")    
    padj = p.adjust(rtab[, 4], method = adj)
    tab.out = data.frame(rtab, padj)
   
    return(tab.out)
}
