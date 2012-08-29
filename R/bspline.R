bspline <-
function(expmat)
{

  require(igraph)  
 
  expmat=data.frame(rownames(expmat),expmat)

  pid=Sys.getpid()
  file=paste("/tmp/",pid,sep="")
  
  # tmp files for mis_calc
  out_tmp = paste(file,".tmp",sep="")
  out_bspline = paste(file,".bspline",sep="")
  write.table(expmat,file=out_tmp,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  
  # run mis_calc
  system(paste("mis_calc --t=4 10 3 --no-qt  <",out_tmp," >",out_bspline," 2>/dev/null",sep=""))

  mim=(file.info(out_bspline)[1]!=0)*1

  if(mim!=0){
   mim=read.csv(out_bspline,sep="\t",header=TRUE)
   net=graph.data.frame(mim,directed=FALSE)
   E(net)$weight=mim[,3]
   mim=get.adjacency(net,attr="weight")
  }

  # removing tmp files

  system(paste("rm",out_tmp))
  system(paste("rm",out_bspline))

return(mim)
}
