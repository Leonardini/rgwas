#' Find a minimal piecewise linear function that ensures that all integer points below it lie on or below the input-defined one.
#' Return the breakpoints of this function, which are guaranteed to be as few as possible in number (see the original paper).
#' The width specifies the width of the "tunnel" used; for best results, it should be close to 1, but strictly smaller than it.
#' Use the standard Imai-Iri algorithm of given width.
callImaiIri = function(extremePoints, width = WIDTH, CPP = TRUE, fname = "Test.csv") {
  stopifnot(width <= 1)
  if (nrow(extremePoints) == 0) {
    return(list(breakpoints = NA, inputSize = 0, input = extremePoints))
  }
  firstPoint <- extremePoints %>%
    dplyr::slice(1) %>%
    dplyr::select(Total, Sum)
  extremePoints %<>%
    dplyr::group_by(Total) %>%
    dplyr::summarize(Sum = max(Sum), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(Total, Sum)
  if (CPP) {
    points <- extremePoints %>%
      dplyr::mutate_at("Sum", ~{magrittr::add(., 1/2)})
    if (nrow(points) <= 2) {
      breakpoints <- points
    } else {
      execPath <- system.file("bin", "plfoptq", package = "rgwas")
      if (file.exists(execPath)) {
        readr::write_csv(points, fname)
        breakpoints <- system2(execPath, args = c(fname, width/2), stdout = TRUE) %>%
          stringr::str_split_fixed(",", 2) %>%
          magrittr::set_colnames(c("Total", "Sum")) %>%
          tibble::as_tibble() %>%
          dplyr::mutate_all(as.numeric)
      } else {
        print("The compiled version of the C++ code is not available; falling back on the pure R option.")
        CPP = FALSE
      }
    }
  }
  if (!CPP) {
    extremePoints %<>%
      removeCollinear %>%
      dplyr::mutate_at("Sum", ~{magrittr::add(., 1/2)})
  }
  pointsP <- extremePoints
  pointsM <- extremePoints
  pointsP %<>%
    dplyr::mutate(Sum = Sum + width/2)
  pointsM %<>%
    dplyr::mutate(Sum = Sum - width/2)
  if (!CPP) {
    breakpoints <- ImaiIriAlgorithm(pointsP, pointsM)
  }
  # if (nrow(breakpoints) > 0) { # check that the required bounds are satisfied by the output
  #   print("Verification of correctness")
  #   pointsL <- pointsM %>%
  #     dplyr::mutate(Sum = round(Sum))
  #   pointsU <- pointsP %>%
  #     dplyr::mutate(Sum = round(Sum))
  #   stopifnot(checkPosition(breakpoints, pointsL, below = TRUE))
  #   stopifnot(checkPosition(breakpoints, pointsU, below = FALSE))
  # }
  if (nrow(breakpoints) > 0) {
    if (min(breakpoints$Sum) != firstPoint$Sum) {
      breakpoints %<>%
        dplyr::bind_rows(firstPoint) %>%
        dplyr::arrange(Total, Sum)
    }
  }
  output <- list(breakpoints = breakpoints, inputSize = nrow(pointsP), input = pointsP)
  output
}

#' Check to ensure that the input endpoints define a boundary lying entirely above the PWLF defined by the checkpoints
#' If below = TRUE, check that the boundary lies below the PWLF instead (abbreviation used: PWLF = piecewise linear function).
#' If reportCoincident = FALSE, also return a list of all checkpoints that coincide exactly with one of the endpoints.
#' Return TRUE if the condition is fulfilled; otherwise print out the details of the earliest found violation, return FALSE.
checkPosition = function(endpoints, checkpoints, below = TRUE, reportCoincident = FALSE) {
  checkpoints %<>%
    as.matrix
  endpoints %<>%
    as.matrix
  stopifnot(all(dplyr::near(range(endpoints[,1]), range(checkpoints[,1]), tol = TOL)))
  correct <- TRUE
  if (nrow(endpoints) == 1) {
    return((below && endpoints[,2] >= checkpoints[,2]) || (!below && endpoints[,2] <= checkpoints[,2]))
  }
  pos <- 1
  curLeft  <- endpoints[pos,]
  curRight <- endpoints[pos + 1,]
  for (ind in seq_len(nrow(checkpoints))) {
    curPoint <- checkpoints[ind,]
    while (!(curPoint[1] - curLeft[1] >= -TOL && curPoint[1] - curRight[1] <= -TOL) && pos < nrow(endpoints) - 1) {
      pos %<>%
        magrittr::add(1)
      curLeft  <- endpoints[pos,]
      curRight <- endpoints[pos + 1,]
    }
    curAngle <- getAngle(rbind(curLeft, curPoint, curRight))
    if ((below && (curAngle > TOL)) || ((!below) && (curAngle < -TOL))) {
      correct <- FALSE
      print(paste0("Verification failed at point number ", ind, " with coordinates (", curPoint[1], ", ", curPoint[2], ")"))
      print(paste0("It should lie ", ifelse(below, "above", "below"), " the segment between points ", pos, " and ", pos + 1))
      print(paste0("Their coordinates are: (", curLeft[1], ", ", curLeft[2], ") and (", curRight[1], ", ", curRight[2], ")"))
      break
    }
  }
  correct
}
