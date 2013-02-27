
mimwrap <-
function (dataset,estimator="pearson",disc="none") {

   mi1=c("pearson", "spearman","kendall","gauss")
   mi2=c("mi.empirical","mi.mm","mi.shrink","mi.sg")
   
   if(estimator %in% mi1){
      mim=build.mim(t(dataset), estimator = estimator, disc = disc)
      mim[is.na(mim)] = 0 
   } else if(estimator %in% mi2){
      if(disc=="none") {
        disc="equalwidth"
      }
      mim=build.mim(t(dataset), estimator = estimator, disc = disc)
      mim[is.na(mim)] = 0     
   } else if(estimator=="cor"){
      mim=cor(t(dataset))
   }
   
return(mim)
}

