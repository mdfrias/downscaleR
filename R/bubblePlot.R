#' @title Bubble plot for visualization of the skill of an ensemble forecast prediction
#' 
#' @description Bubble plot for the visualization of the skill of an ensemble forecast prediction. It provides a
#'  spatially-explicit representation of the skill, resolution and reliability of a probabilistic predictive system in
#'  a single map.
#' 
#' @param mm.obj A multi-member object with predictions, either a field or a multi-member station object as a result of
#' downscaling of a forecast using station data. See details.
#' @param obs The benchmarking observations for forecast verification. 
#' @param select.year Year within the whole verification period to display the results for.
#' @param score Logical indicating if the relative operating characteristic skill score (ROCSS) should be included. See 
#'  details. Default is TRUE
#' @param size.as.probability Logical indicating if the tercile probabilities (magnitude proportional to bubble radius) 
#'  are drawn in the plot. See details. Default is TRUE.
#' @param pie. Logical flag indicating if pie charts should be plot. Default is FALSE.
#' 
#' @importFrom scales alpha
#' @importFrom verification roc.area
#' @importFrom mapplots draw.pie
#' @importFrom mapplots add.pie
#' 
#' @export
#' 
#' @details 
#' 
#' For each member, the daily predictions are averaged to obtain a single seasonal forecast. The corresponding terciles 
#' for each ensemble member are then computed for the analysis period. Thus, each particular grid point, member and season,
#' are categorized into three categories (above, between or below), according to their respective climatological 
#' terciles. Then, a probabilistic forecast is computed year by year by considering the number of members falling 
#' within each category. For instance, probabilities below 1/3 are very low, indicating that a minority of the members 
#' falls in the tercile. Conversely, probabilities above 2/3 indicate a high level of member agreement (more than 66\% of members
#' falling in the same tercile). Color represents the tercile with the highest probability for the selected year. The bubble size 
#' indicates the probability of that tercile. This option is not plotted if the size.as.probability argument is FALSE.
#' 
#' Finally, the ROC Skill Score (ROCSS) is computed. For each tercile, it provides a quantitative measure of the forecast skill,
#' and it is commonly used to evaluate the performance of probabilistic systems (Joliffe and Stephenson 2003). The value of 
#' this score ranges from 1 (perfect forecast system) to -1 (perfectly bad forecast system). A value zero indicates no skill 
#' compared with a random prediction. The transparency of the bubble is associated to the ROCSS (negative values are
#' plotted with x).  This option is not plotted if the score argument is FALSE.
#' 
#' 
#' @note The computation of climatological terciles requires a representative period to obtain meaningful results.
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}, M.D. Frias and J. Fernandez 
#' 
#' @family visualization
#' 
#' @references
#'  Jolliffe, I. T. and Stephenson, D. B. 2003. Forecast Verification: A Practitioner's Guide in 
#'  Atmospheric Science, Wiley, NY
#'  


bubblePlot <- function(mm.obj, obs, select.year, score = TRUE, size.as.probability = TRUE, pie = FALSE, only.at = NULL) {      
      mm.dimNames <- attr(mm.obj$Data, "dimensions")
      obs.dimNames <- attr(obs$Data, "dimensions")
      if (!("member" %in% mm.dimNames)) {
            stop("The input data for 'multimember' is not a multimember field")
      }
      if ("member" %in% obs.dimNames) {
            stop("The verifying observations can't be a multimember")
      }
      if ("var" %in% mm.dimNames | "var" %in% obs.dimNames) {
            stop("Multifields are not a valid input")
      }
#       if (!("lon" %in% obs.dimNames & "lat" %in% obs.dimNames)) {
#             stop("The observed reference must have 'lon' and 'lat' dimensions")
#       }
      if (!identical(c("time", "lat", "lon"), obs.dimNames)) {
            stop("The observed reference must be a 3D array of the form [time,lat,lon]")
      }
      x.mm <- mm.obj$xyCoords$x; nx <- length(x.mm)
      y.mm <- mm.obj$xyCoords$y; ny <- length(y.mm)
      yrs <- getYearsAsINDEX(mm.obj)
      yy <- unique(yrs)
      if (!select.year %in% yy) {
            stop("Target year outside temporal data range")
      }
      iyear <- which(yy[1]:yy[length(yy)] == select.year)
      # Why not... ?
      # iyear <- which(yy == select.year)
      lon.dim <- grep("lon", mm.dimNames)
      lat.dim <- grep("lat", mm.dimNames)
      member.dim <- grep("member", mm.dimNames)
      time.dim <- grep("time", mm.dimNames)
      n.mem <- dim(mm.obj$Data)[member.dim]

      # Computation of terciles and exceedance probabilities
      # yearmean changes the data dimension. time.dim is in the first dimension!!

      y.mean <- apply(mm.obj$Data, MARGIN = c(lat.dim, lon.dim, member.dim), FUN = function(x) {
            tapply(x, INDEX = yrs, FUN = mean, na.rm = TRUE)
      })
      terciles <- apply(y.mean, MARGIN=c(2, 3, 4), FUN=quantile, c(1/3, 2/3), na.rm = TRUE)
      # Compute the probability for each tercile
      prob <- array(dim = c(3,dim(y.mean)[1:3]))
      for (i in 1:(dim(y.mean)[1]) ){
            prob[3, i, , ] <- apply(y.mean[i, , , ] > terciles[2, , , ], MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) / n.mem
            prob[1, i, , ] <- apply(y.mean[i, , , ] < terciles[1, , , ], MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) / n.mem
            prob[2, i, , ] <- 1 - prob[3, i, , ] - prob[1, i, , ]
      }
      prob <- ifelse(prob<0,0,prob) # round-off errors lead to some tiny neg. values
      # Tercile for the maximum probability from the terciles
      t.max.prob <- apply(prob, MARGIN = c(2, 3, 4), FUN = which.max)
      # ... as index matrix for the target year
      idxmat.max.prob <- cbind(c(t.max.prob[iyear,,]), 1:prod(dim(t.max.prob[iyear,,])))
      # Probability of the most likely tercile
      max.prob <- apply(prob, MARGIN = c(2, 3, 4), FUN = max)
      # Terciles for the observations
      obs.y.mean <- apply(obs$Data, MARGIN = c(2, 3), FUN = function(x) {
            tapply(x, INDEX = yrs, FUN = mean, na.rm = TRUE)
      })
      obs.terciles <- apply(obs.y.mean, MARGIN = c(2, 3), FUN = quantile, c(1/3, 2/3), na.rm = TRUE)
      obs.t <- array(dim = dim(obs.y.mean))
      for (i in 1:(dim(obs.y.mean)[1]) ) {
            obs.t[i, , ] <- (obs.y.mean[i, , ] > obs.terciles[1, , ]) +
                            (obs.y.mean[i, , ] > obs.terciles[2, , ]) + 1
      }
      # Filter points with observations in model data 
      # Select a year and remove the NaN cases detected for the observations. 
      v.max.prob <- as.vector(max.prob[iyear, , ])
      v.t.max.prob <- as.vector(t.max.prob[iyear, , ])
      v.valid <- complete.cases(as.vector(obs.t[iyear, , ]))
      ve.max.prob <- v.max.prob[v.valid]
      v.prob <- array(dim = c(sum(v.valid),3))
      for (i in 1:3){
        v.prob[,i] <- as.vector(prob[i,iyear,,])[v.valid]
      }
      if (!size.as.probability){
            ve.max.prob <- rep(1, length(ve.max.prob))
      }
      df <- data.frame(max.prob = ve.max.prob, t.max.prob = v.t.max.prob[v.valid])
      df$color <- "black"
      t.colors <- c("blue", "darkgrey", "red")
      df$color[df$t.max.prob == 3] <- t.colors[3]
      df$color[df$t.max.prob == 2] <- t.colors[2]
      df$color[df$t.max.prob == 1] <- t.colors[1]
      yx <- as.matrix(expand.grid(y.mm, x.mm))
      nn.yx <- yx[v.valid, ]
      if (score) { # Compute ROCSS for all terciles
        rocss.fun <- function(o,p){
          if (length(unique(o))==2){
            suppressWarnings(roc.area(o,p)$A*2-1)
          } # ROCSS cannot be computed if all obs are equal.
          else NA
        }
        rocss <- array(dim=dim(prob)[-2]) # no year dim
            for (i.tercile in 1:3){
                  rocss[i.tercile, , ] <- apply(
                       array(c(obs.t==i.tercile, prob[i.tercile, , ,]), dim=c(dim(obs.t),2)),
                       MAR=c(2,3),
                       FUN=function(x){rocss.fun(x[,1],x[,2])})
            }
            if (!pie) { # Select the rocss for the tercile with max prob.
              rocss <- unshape(rocss)
              max.rocss <- rocss[idxmat.max.prob]
              rocss <- deunshape(rocss)
              dim(max.rocss) <- dim(rocss)[-1]
              v.score <- c(max.rocss)[v.valid]
              pos.val <- v.score >= 0
              neg.val <- v.score < 0
            }
      }
      # Bubble plot
      mons <- unique(months(as.POSIXlt(obs$Dates$start), abbreviate = T))
      title <- sprintf("%s, %s to %s, %d", obs$Variable$varName, mons[1],last(mons), year.target)
      par(bg = "white", mar = c(3, 3, 2, 1))
      plot(0, xlim=range(x.mm), ylim=range(y.mm), type="n")
      mtext(title, side=3, line=0.5, at=min(x.mm), adj=0, cex=1.2)
      symb.size <- (df$max.prob-0.33) * 4
      if (pie){
          dx <- diff(x.mm[1:2])
          dy <- diff(y.mm[1:2])
          radius <- min(dx,dy)/2*0.8
          if (score){
            rocss <- unshape(rocss)
            col <- array(dim=dim(rocss))
            t.colors[2] <- gray(0.2)
            for (i.tercile in 1:3) {col[i.tercile,] <- alpha(t.colors[i.tercile],255*rocss[i.tercile,])}
            score.valid <- ! is.na(rocss[1,v.valid])
          } else {
            col <- matrix(rep(t.colors, dim(nn.yx)[1]), nrow = 3)
            score.valid <- rep(TRUE, sum(v.valid))
          }
          if (!is.null(only.at)) {
            ns <- nrow(only.at$Stations$LonLatCoords)
            score.valid <- rep(FALSE, length(score.valid))
            for (i.station in 1:ns){
              score.valid[which.min((only.at$Stations$LonLatCoords[i.station,1] - nn.yx[,2])^2 +
                                   (only.at$Stations$LonLatCoords[i.station,2] - nn.yx[,1])^2)] <- TRUE
            }
          }
          for (i.loc in which(score.valid)){
              add.pie(v.prob[i.loc,], nn.yx[i.loc, 2], nn.yx[i.loc, 1], col=col[,i.loc],
                      radius=radius, init.angle=90, clockwise = F, border="lightgray", labels=NA
              )  
          }
          # Highlight those whose ROCSS cannot be computed due to constant obs conditions (e.g. always dry)
          if (!is.null(only.at)) {score.valid <- rep(TRUE, sum(v.valid))}
          for (i.loc in which(!score.valid)){
            add.pie(v.prob[i.loc,], nn.yx[i.loc, 2], nn.yx[i.loc, 1], col=NA,
                    radius=radius, init.angle=90, clockwise = F, border="green", labels=NA
            )  
          }
      } else {
        if (score) {
            points(nn.yx[pos.val, 2], nn.yx[pos.val, 1], cex = symb.size[pos.val], col = alpha(df$color[pos.val], 255*v.score[pos.val]), pch = 16, xlab = "", ylab = "")
        } else {
            points(nn.yx[ , 2], nn.yx[ , 1], cex=symb.size, col = df$color, pch = 16, xlab = "", ylab = "")
        }
      }
      #points(nn.yx[neg.val, 2], nn.yx[neg.val, 1], pch=4, cex=0.5) # To add negative values
      #flat.val <- is.na(v.score)
      # 
      #points(nn.yx[flat.val, 2], nn.yx[flat.val, 1], pch=5, cex=0.5) # To add obs constant values
      world(add = TRUE, interior = T)      
      world(add = TRUE, interior = F, lwd=3)    
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend('topright', c("Below", "Normal", "Above"), pch=c(19, 19, 19), col = t.colors, horiz = T, inset = c(0, 0), xpd = TRUE, bty = "n")      
}
# End
