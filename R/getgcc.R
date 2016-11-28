getgcc <- function(net){
  
  netclust = clusters(net)
  gcc = V(net)$name[netclust$membership==which.max(netclust$csize)]
  gccnet = induced.subgraph(net, gcc)

  return(gccnet)
}


