alpha <- function(col, alpha){
  rgb(t(col2rgb(col)), alpha=alpha, maxColorValue = 255)
}

tercileBrewerColorRamp <- function(ncolors){
  ncolors <- ncolors-2
  tcols <- tercileColor()
  tmp <- brewer.pal(ncolors, tcols[1])
  rval <- data.frame(low=c(tmp[1],tmp[1],tmp), stringsAsFactors=F)
  tmp <- brewer.pal(ncolors, tcols[2])
  rval$middle <- c(tmp[1],tmp[1],tmp)
  tmp <- brewer.pal(ncolors, tcols[3])
  rval$high <- c(tmp[1],tmp[1],tmp)
  return(rval)
}

tercileColor <- function(){
  #c("red", "yellow", "blue")
  c("Blues", "Greys", "Reds")
}

tercileColorRamp <- function(ncolors){
  alphas <- seq(0,255,len=ncolors)
  tcols <- tercileColor()
  t.low <- alpha(tcols[1],alphas)
  t.mid <- alpha(tcols[2],alphas)
  t.hi <- alpha(tcols[3],alphas)
  return(data.frame(t.low,t.mid,t.hi))

}
