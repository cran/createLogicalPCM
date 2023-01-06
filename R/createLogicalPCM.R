#' @title Simulated Logical Pairwise Comparison Matrix for Analytic Hierarchy Process
#'
#' @description Creates a logical pairwise comparison matrix for the Analytic Hierarchy Process such as would be created by a rational decision maker based on a relative vector of preferences for the alternatives involved. Choices of the pairwise comparison ratios are from the Fundamental Scale and simulate a reasonable degree of error. The algorithm is modified from a paper by Bose, A [2022], \doi{https://doi.org/10.1002/mcda.1784}
#' @param ord The desired order of the Pairwise Comparison Matrix
#' @param prefVec The preference vector of length as the order of the input matrix
#' ' @return A Logical Pairwise Comparison Matrix
#' @importFrom stats runif
#' @examples
#' lPCM <- createLogicalPCM(3,c(1,2,3));
#' lPCM <- createLogicalPCM(5,c(0.25,0.4,0.1,0.05,0.2));
#' @export
createLogicalPCM <- function(ord, prefVec=rep(NA,ord)) {
  if (is.na(ord)) stop("The first parameter is mandatory")
  if (!is.numeric(ord) || ord %% 1 != 0) stop("The first parameter has to be an integer")
  if (!all(is.na(prefVec)) && !is.numeric(prefVec)) stop("The second parameter has to be a numeric vector")
  if (!all(is.na(prefVec)) && length(prefVec)!=ord) stop("The length of the second parameter has to be the same as the first parameter")
  opt <- c(1/(9:2),1:9)

  if (is.na(prefVec[1]))
    prefVec <- runif(ord)

  mperfect <- outer(prefVec, prefVec, "/")

  bestMat <- bestM(mperfect)
  # bm <- abs(Re(eigen(bestMat)$vectors[,1]))/sum(abs(Re(eigen(bestMat)$vectors[,1])))
  # bevl <- abs(Re(eigen(bestMat)$values[1])) -  nrow(bestMat)

  PCM <- list(ord=ord, orgVec=prefVec, ppcm=bestMat)

  # creating a logical PCM from a consistent PCM (bestM)
  m <-  PCM$ppcm
  ord <- PCM$ord
  # now creating a logical PCM
  for (r in 1:(ord-1)) {
    for (c in (r+1):ord) {
      m1 <- which.min(abs(opt-m[r,c]))
      m2 <- which.min(abs(opt[-m1]-m[r,c]))
      m3 <- which.min(abs(opt[-c(m1,m2)]-m[r,c]))
      # random choice from the nearest 3
      allChoices <- choices <- c(m1, m2, m3)
      if (m[r,c] >= 1) {
        choices <- allChoices[opt[allChoices] >= 1]
      } else if (m[r,c] < 1) {
        choices <- allChoices[opt[allChoices] <= 1]
      }
      m[r,c] <- sample(opt[choices],1)
      m[c,r] <- 1/m[r,c]
    }
  }
  e <- Re(eigen(m)$vectors[,1])/sum(Re(eigen(m)$vectors[,1]))

  mdev <- matrix(0,ord, ord)
  for (i in 1:ord)   mdev[,i] <- abs(rank(-m[,i])-rank(-e))

  PCM$pcm <- m
  # PCM$evc <- e
  # ev <- Re(eigen(m)$values[1]) - ord
  # ev <- ifelse(abs(ev)<eps,0,ev)
  # PCM$evl <- ev
  #return(list(logicalPCM=PCM$pcm, bestPCM=PCM$ppcm))
  return(logicalPCM=PCM$pcm)
}


bestM <- function(pcm) {
  tSc <- c(1/(9:1),2:9)
  p <- pcm
  o <- nrow(pcm)
  bestMatrix <- diag(o)
  ep <- abs(Re(eigen(p)$vectors[,1]))
  for (r in 1:(o-1))
    for (c in (r+1):o) {
      #b <- closest(ep[r]/ep[c], 9)    # fails for sl9T[[119]]$pcm[1,3]
      b <- tSc[which.min(abs(ep[r]/ep[c]-tSc))[1]]
      bestMatrix[r, c] <- b
      bestMatrix[c, r] <- 1/b
    }
  return(bestMatrix)
}

