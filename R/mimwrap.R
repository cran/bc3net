
mimwrap <-
function (dataset,estimator="pearson",disc="none") {

   minetest=c("pearson", "spearman","kendall","gauss","mi.empirical","mi.mm",
   "mi.shrink","mi.sg")

   if(estimator %in% minetest){
      mim=build.mim(t(dataset), estimator = estimator, disc = disc)
      mim[is.na(mim)] = 0
   } else if(estimator=="bspline"){
      mim=bspline(dataset) 
   } else if(estimator=="cor"){
      mim=cor(t(dataset))
   }
   
return(mim)
}

