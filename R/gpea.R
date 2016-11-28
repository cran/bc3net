gpea <- function(gnet, genesets, verbose=TRUE, cmax=1000, cmin=3, adj="bonferroni"){

    
  genes = V(gnet)$name  
  genesets = lapply(genesets, function(x) x[x%in%genes])
  genesets = genesets[sapply(genesets, function(x) length(x)>=cmin & length(x)<cmax)]


  # considering only genes defined in the genesets
  allgenes = unique(unlist(genesets))

  ref.net = induced.subgraph(gnet, allgenes)
  v.ref = vcount(ref.net)
  e.ref = ecount(ref.net)

  ref.clique = (vcount(ref.net)*(vcount(ref.net)-1))/2

  res=matrix(0, ncol=4, nrow=length(genesets))
  rownames(res)=names(genesets)
  colnames(res)=c("edges", "genes", "pval", "padj")

 
   for(i in 1:nrow(res)){

      
       if(verbose==TRUE){
           cat(i," of ",nrow(res),"\n")
       }
       
       term = genesets[[i]]
       sub.net = induced.subgraph(gnet, term)

       v.sub=vcount(sub.net)
       e.sub=ecount(sub.net)

       e.clique=(v.sub*(v.sub-1))/2

       
     # sub of genesets[[i]] in GRN | sub not in GRN (clique)  
     # not in genesets[[i]] in GRN | not in genesets[[i]] (clique)

     # Example: #######################################
     #     170 | 2,680
     #  59,933 | 54,561,831
     ############################## p=8.15e-227
     
     mat = matrix(c(e.sub, e.ref-e.sub,
         e.clique-e.sub, ref.clique-e.clique-e.ref+e.sub),
         ncol=2, nrow=2)
       
     pval = fisher.test(mat, alternative="greater")$p.value

     res[i,1:3]=c(e.sub, v.sub, pval)
   }

   # res=res[res[,2]<=cmax & res[,2]>cmin,]
   res=res[order(res[,3]),]
   res[,4]=p.adjust(res[,3], method=adj)  

   
   tab=data.frame(TermID=rownames(res),res)

   return(tab)
}
