alpha <- function(col, alpha){
  alpha <- ifelse(is.na(alpha),0,alpha) # set NA to fully transparent
  alpha <- ifelse(alpha<0,0,alpha) # set negative to fully transparent
  rgb(t(col2rgb(col)), alpha=alpha, maxColorValue = 255)
}

deunshape <- function(x){
  dim(x) <- attr(x,"shape")
  attr(x,"shape") <- NULL
  return(x)
}

fldsubsample <- function(dsRobj, one.out.of = 2){
  dsRobj.dimNames <- attr(dsRobj$Data, "dimensions")
  lon.dim <- grep("lon", dsRobj.dimNames)
  lat.dim <- grep("lat", dsRobj.dimNames)
  nx <- length(dsRobj$xyCoords$x)
  ny <- length(dsRobj$xyCoords$y)
  filter.x <- seq(1,nx,by=one.out.of)
  filter.y <- seq(1,ny,by=one.out.of)
  rval <- dsRobj
  rval$xyCoords$x <- dsRobj$xyCoords$x[filter.x]
  rval$xyCoords$y <- dsRobj$xyCoords$y[filter.y]
  attr(rval$xyCoords, "resX") <- attr(dsRobj$xyCoords, "resX")*one.out.of
  attr(rval$xyCoords, "resY") <- attr(dsRobj$xyCoords, "resY")*one.out.of
  rval$Data <- extract(dsRobj$Data,
                       indices=list(filter.x,filter.y),
                       dims=c(lon.dim, lat.dim)
  )
  attr(rval$Data, "dimensions") <- attr(dsRobj$Data, "dimensions")
  return(rval)
}

summary.dsRgrid <- function(x){
  cat(sprintf("Variable: %s\n", x$Variable$varName))
  cat("Shape: ", dim(x$Data), "(", attr(x$Data, "dimensions"), ")\n")
  nx <- length(x$xyCoords$x)
  ny <- length(x$xyCoords$y)
  nt <- length(x$Dates$start)
  cat(sprintf("X range: %6.2f to %6.2f (%d values)\n", min(x$xyCoords$x), max(x$xyCoords$x),nx))
  cat(sprintf("Y range: %6.2f to %6.2f (%d values)\n", min(x$xyCoords$y), max(x$xyCoords$y), ny))
  available.months <- sort(unique(as.integer(substr(x$Dates$start, 6,7))))
  if (identical(available.months, 1:12)) {
    available.months <- "All year"
  } else {
    available.months <- sprintf("Season: %s", paste(available.months, collapse=' '))
  }
  cat(sprintf("time range: %s to %s (%d values, %s)\n", x$Dates$start[1], x$Dates$end[nt], nt, available.months))
  if ("member" %in% attr(x$Data, "dimensions"))
    member.dim <- grep("member", attr(x$Data, "dimensions"))
  cat(sprintf("Members: %d\n", dim(x$Data)[member.dim]))
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

unshape <- function(x, MAR=c(1)){
  attr(x, "shape") <- dim(x)
  ndim <- length(dim(x))
  dim(x) <- c(dim(x)[MAR],prod(dim(x)[-MAR]))
  return(x)
}

yearmean <- function(dsRobj, yrOfSeasonEnd = TRUE){
  dsRobj.dimNames <- attr(dsRobj$Data, "dimensions")
  lon.dim <- grep("lon", dsRobj.dimNames)
  lat.dim <- grep("lat", dsRobj.dimNames)
  member.dim <- grep("member", dsRobj.dimNames)
  time.dim <- grep("time", dsRobj.dimNames)
  n.mem <- dim(dsRobj$Data)[member.dim]
  if (yrOfSeasonEnd) {
    yrs <- getYearsAsINDEX(dsRobj)
  } else {
    yrs <- as.integer(substr(dsRobj$Dates$start, 0,4))
  }
  rval <- dsRobj
  margin <- c(lat.dim, lon.dim, member.dim)
  rval$Data <- apply(dsRobj$Data, MARGIN = margin,
    FUN = function(x) {tapply(x, INDEX = yrs, FUN = mean, na.rm = TRUE)}
  )
  # Recover the original dimension order
  newdims <- c("time", dsRobj.dimNames[margin])
  rval$Data <- aperm(rval$Data, match(dsRobj.dimNames, newdims))
  attr(rval$Data, "dimensions") <- dsRobj.dimNames
  rval$Dates$start <- tapply(dsRobj$Dates$start, INDEX=yrs, FUN=min)
  rval$Dates$end <- tapply(dsRobj$Dates$end, INDEX=yrs, FUN=max)
  return(rval)
}