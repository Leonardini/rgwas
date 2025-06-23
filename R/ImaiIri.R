#' Turn an input 2-dimensional vector into a unit vector in the same direction
#' @noRd
normalizeVector = function(vector) {
  len <- pracma::hypot(vector[1], vector[2])
  if (dplyr::near(len, 0)) {
    return(vector)
  }
  normVector <- vector/sqrt(len)
  normVector
}

#' Determine the (sign of the) angle by which we need to turn to get from point 1 to point 3 via point 2 (ignore magnitudes!)
#' The return value is positive if the turn is a right turn, negative if the turn is a left turn, 0 if the points are collinear
#' @noRd
getAngle = function(threePointSet) {
  stopifnot(all(dim(threePointSet) == c(3, 2)))
  threePointSet <- threePointSet %>%
    as.matrix
  point1 <- threePointSet[1, ]
  point2 <- threePointSet[2, ]
  point3 <- threePointSet[3, ]
  if (comparePoints(point1, point2) || comparePoints(point2, point3)) { # degenerate situation!
    return(0)
  }
  vector1 <- normalizeVector(point2 - point1)
  vector2 <- normalizeVector(point2 - point3)
  angle <- vector1[1] * vector2[2] - vector2[1] * vector1[2]
  angle
}

#' Remove all but the first and the last point of any group of three or more consecutive input points that are collinear.
#' Note that this function is idempotent (at least in exact arithmetic), i.e. running it twice is the same as running it once.
#' @noRd
removeCollinear = function(pointSet, tol = TOL) {
  N <- nrow(pointSet)
  if (N <= 2) {
    return(pointSet)
  }
  badInds <- rep(FALSE, N)
  for (ind in 2:(N - 1)) {
    if (abs(getAngle(pointSet[(ind - 1):(ind + 1), ])) < tol) {
      badInds[ind] <- TRUE
    }
  }
  keepPoints <- which(!badInds)
  pointSet <- pointSet %>%
    dplyr::slice(keepPoints)
  pointSet
}

#' Find the slope and y-intercept of a line given by a 2 x 2 matrix, whose rows are points and columns are x and y coordinates
#' Note: special cases are defined in the dependent functions for horizontal lines (slope = 0) and vertical lines (slope = Inf)
#' @noRd
findSlopeAndIntercept = function(endPoints) {
  stopifnot(nrow(endPoints) == 2)
  endPoints <- endPoints %>%
    as.matrix
  delta <- endPoints[2, ] - endPoints[1, ]
  slope <- delta[2]/delta[1]
  if (is.finite(slope)) {
    intercept <- endPoints[1, 2] - slope * endPoints[1, 1]
  } else if (is.na(slope)) {
    stop("The two endpoints are the same!")
  } else {
    intercept <- endPoints[1, 1]
  }
  output <- c(slope, intercept)
  output
}

#' Find the intersection between two lines, each specified by a pair of points (2 x 2 matrices; row = point, columns = coords)
#' Note: does not check that the intersection point lies between the starting and ending points of each segment!
#' @noRd
findIntersection = function(endPoints1, endPoints2) {
  line1 <- findSlopeAndIntercept(endPoints1)
  line2 <- findSlopeAndIntercept(endPoints2)
  if (is.finite(line1[1]) && is.finite(line2[1])) {
    delta <- line2 - line1
    if (dplyr::near(delta[1], 0)) {
      stop("The lines are parallel or near parallel!")
    }
    x0 <- -delta[2]/delta[1]
    y0 <- line1[1] * x0 + line1[2]
  } else if (is.finite(line1[1])) {
    x0 <- line2[2]
    y0 <- line1[1] * x0 + line1[2]
  } else if (is.finite(line2[1])) {
    x0 <- line1[2]
    y0 <- line2[1] * x0 + line2[2]
  } else {
    stop("The lines are both vertical!")
  }
  output <- c(x0, y0)
  output
}

#' Compare two points numerically, up to a specified number of dimensions Dim (which can be smaller than the full dimension)
#' Return TRUE if they are numerically identical, FALSE otherwise. Do not check that the points have the same dimensions!
#' @noRd
comparePoints = function(point1, point2, Dim = 2) {
  comp <- TRUE
  for (ind in 1:Dim) {
    if (!dplyr::near(point1[ind], point2[ind], tol = TOL)) {
      comp <- FALSE
    }
  }
  comp
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
computeIndex <- function(index, sign, N) {
  allSigns <- c("positive", "negative", "computed", "other")
  return(index + (which(allSigns == sign) - 1) * N)
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
getIndex <- function(point) {
  return(point %>% dplyr::pull(rowid))
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
getPoint <- function(pointSet, index) {
  return(pointSet %>% dplyr::slice(index))
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
getPrevPoint <- function(pointSet, prevMap, point) {
  return(getPoint(pointSet, prevMap[getIndex(point)]))
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
getNextPoint <- function(pointSet, nextMap, point) {
  return(getPoint(pointSet, nextMap[getIndex(point)]))
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
setPoint <- function(pointSet, index, point) {
  pointSet[index, ] <- point
  pointSet
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
setPrevPoint <- function(prevMap, curPoint, prevPoint) {
  prevMap[getIndex(curPoint)] <- getIndex(prevPoint)
  prevMap
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
setNextPoint <- function(nextMap, curPoint, nextPoint) {
  nextMap[getIndex(curPoint)] <- getIndex(nextPoint)
  nextMap
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
makePoint <- function(point, index) {
  output <- tibble::tibble(Total = point[1], Sum = point[2], rowid = index)
  output
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
getSignedAngle <- function(point1, point2, point3, sign) {
  myPoints <- dplyr::bind_rows(point1, point2, point3) %>%
    dplyr::select(-rowid) %>%
    as.matrix
  myAngle <- getAngle(myPoints) * SIGNS[sign]
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
checkAcute <- function(point1, point2, point3, sign, tol = TOL) {
  myAngle <- getSignedAngle(point1, point2, point3, sign)
  return(myAngle > tol)
}

#' Auxiliary function for the Imai-Iri algorithm
#' @noRd
checkObtuse <- function(point1, point2, point3, sign, tol = TOL) {
  myAngle <- getSignedAngle(point1, point2, point3, sign)
  return(myAngle < -tol)
}

#' The Imai-Iri algorithm for finding the shortest piecewise linear path lying between pointsP and pointsM.
#' Assume without checking that the x-coordinates of pointsP and pointsM are identical, and ordered from smallest to largest.
#' Also assume without checking that the y-coordinates of pointsP are pointwise larger than those of pointsM.
#' This implementation is a corrected version of the pseudocode in Sabine Neubauer's student thesis (University of Karlsruhe).
#' @noRd
ImaiIriAlgorithm = function(pointsP, pointsM) {
  n <- nrow(pointsP)
  stopifnot(nrow(pointsM) == n)
  stopifnot(ncol(pointsP) == 2)
  stopifnot(ncol(pointsM) == 2)
  if (n <= 2) {
    return(pointsM)
  }
  # Global numbering scheme: the positive points are 1:n; the negative points are (n+1):2n; the computed Q points are (2n+1):4n
  allQ <- matrix(NA, n, 2, dimnames = list(c(), colnames(pointsP))) %>%
    tibble::as_tibble()
  allPoints <- dplyr::bind_rows(pointsP, pointsM, allQ, allQ) %>%
    tibble::rowid_to_column() %>%
    dplyr::select(-rowid, rowid) %>%
    dplyr::mutate_all(as.double)
  nextMap <- rep(NA, nrow(allPoints))
  names(nextMap) <- 1:nrow(allPoints)
  prevMap <- rep(NA, nrow(allPoints))
  names(prevMap) <- 1:nrow(allPoints)
  allP <- vector("list")
  allL <- vector("list")
  allR <- vector("list")
  for (sign in c("positive", "negative")) {
    nextMap <- nextMap %>%
      setNextPoint(getPoint(allPoints, computeIndex(1, sign, n)), getPoint(allPoints, computeIndex(2, sign, n)))
    prevMap <- prevMap %>%
      setPrevPoint(getPoint(allPoints, computeIndex(2, sign, n)), getPoint(allPoints, computeIndex(1, sign, n)))
    allP[[sign]] <- getPoint(allPoints, computeIndex(1, sign, n))
    allL[[sign]] <- getPoint(allPoints, computeIndex(1, sign, n))
    allR[[sign]] <- getPoint(allPoints, computeIndex(2, sign, n))
  }
  j <- 1
  print(paste("There are", n, "boundary points to process"))
  for (i in 3:n) {
    if (i %% 500 == 0) {
      print(paste("Processed", i, "boundary points so far"))
    }
    nextWindow <- FALSE
    for (sign in c("positive", "negative")) {
      curPrev  <- getPoint(allPoints, computeIndex(i - 1, sign, n))
      curPoint <- getPoint(allPoints, computeIndex(i, sign, n))
      curP     <- allP[[sign]]
      while (!comparePoints(curPrev, curP) && checkObtuse(curPoint, curPrev, getPrevPoint(allPoints, prevMap, curPrev), sign)) {
        curPrev <- getPrevPoint(allPoints, prevMap, curPrev)
      }
      nextMap <- nextMap %>%
        setNextPoint(curPrev, curPoint)
      prevMap <- prevMap %>%
        setPrevPoint(curPoint, curPrev)
    }
    for (index in 1:2) {
      signStar <- ifelse(index == 1, "positive", "negative")
      signDiam <- ifelse(index == 1, "negative", "positive")
      curPoint <- getPoint(allPoints, computeIndex(i, signStar, n))
      curL <- allL[[signStar]]
      curR <- allR[[signDiam]]
      if (!nextWindow && checkAcute(curPoint, curL, curR, signStar)) {
        indQ <- computeIndex(j, "computed", n)
        myIntersect <- makePoint(findIntersection(dplyr::bind_rows(curL, curR), dplyr::bind_rows(allP[[signStar]], allP[[signDiam]])), indQ)
        allPoints <- allPoints %>%
          setPoint(indQ, myIntersect)
        allP[[signDiam]][1,] <- allR[[signDiam]][1,]
        prevP <- getPoint(allPoints, computeIndex(i - 1, signStar, n))
        curP  <- getPoint(allPoints, computeIndex(i, signStar, n))
        indR <- computeIndex(j, "other", n)
        otherIntersect <- makePoint(findIntersection(dplyr::bind_rows(curL, curR), dplyr::bind_rows(prevP, curP)), indR)
        allPoints <- allPoints %>%
          setPoint(indR, otherIntersect)
        allP[[signStar]][1,] <- otherIntersect
        j <- j %>%
          magrittr::add(1)
        nextMap <- nextMap %>%
          setNextPoint(otherIntersect, curP)
        prevMap <- prevMap %>%
          setPrevPoint(curP, otherIntersect)
        allR[[signStar]][1,] <- getPoint(allPoints, computeIndex(i, signStar, n))
        allR[[signDiam]][1,] <- getPoint(allPoints, computeIndex(i, signDiam, n))
        allL <- allP
        while (checkAcute(allL[[signDiam]], allR[[signStar]], getNextPoint(allPoints, nextMap, allL[[signDiam]]), signDiam)) {
          allL[[signDiam]] <- getNextPoint(allPoints, nextMap, allL[[signDiam]])
        }
        nextWindow <- TRUE
      }
    }
    if (!nextWindow) {
      for (index in 1:2) {
        signStar <- ifelse(index == 1, "positive", "negative")
        signDiam <- ifelse(index == 1, "negative", "positive")
        curP <- getPoint(allPoints, computeIndex(i, signStar, n))
        if (checkAcute(curP, allL[[signDiam]], allR[[signStar]], signStar)) {
          allR[[signStar]][1,] <- getPoint(allPoints, computeIndex(i, signStar, n))
          while (checkAcute(curP, allL[[signDiam]], getNextPoint(allPoints, nextMap, allL[[signDiam]]), signStar)) {
            allL[[signDiam]][1,] <- getNextPoint(allPoints, nextMap, allL[[signDiam]])
          }
        }
      }
    }
  }
  m <- j + 1
  if (comparePoints(allL[["positive"]], allP[["positive"]])) {
    supportLine <- rbind(allL[["positive"]], allR[["negative"]])
  } else {
    supportLine <- rbind(allL[["negative"]], allR[["positive"]])
  }
  indQ <- computeIndex(m - 1, "computed", n)
  allPoints <- allPoints %>%
    setPoint(indQ, makePoint(findIntersection(supportLine, rbind(allP[[1]], allP[[2]])), indQ))
  indP <- computeIndex(n, "positive", n)
  indM <- computeIndex(n, "negative", n)
  indR <- computeIndex(m, "computed", n)
  allPoints <- allPoints %>%
    setPoint(indR, makePoint(findIntersection(supportLine, rbind(getPoint(allPoints, indP), getPoint(allPoints, indM))), indR))
  outputInds <- computeIndex(1, "computed", n):indR
  output <- allPoints %>%
    dplyr::slice(outputInds) %>%
    dplyr::select(-rowid)
  output
}
