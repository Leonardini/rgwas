### Finds a short piecewise linear function that ensures that all integer points below it lie on or below the input-defined one.
### Returns the breakpoints of this function, which are guaranteed to be as few as possible in number (see the original paper).
### if dryRun = TRUE, only returns the extreme points to be used, pruned for collinearity; otherwise computes the best boundary.
### If verify = TRUE, also check that the input points lie below (above) the output points but no more than 1 below (above) them.
### The width specifies the width of the "tunnel" used; for best results, it should be close to 1, but strictly smaller than it.
### Uses the standard Imai-Iri algorithm of given width (the restricted = TRUE option, now obsolete, is found in LegacyCode.R).
callImaiIri = function(extremePoints, verify = FALSE, width = WIDTH, dryRun = FALSE, CPP = TRUE, fname = "Test.csv") {
  stopifnot(width <= 1)
  if (nrow(extremePoints) == 0) {
    return(list(breakpoints = NA, inputSize = 0, input = extremePoints))
  }
  firstPoint <- extremePoints %>%
    slice(1) %>%
    select(Total, Sum)
  extremePoints %<>%
    group_by(Total) %>%
    summarize(Sum = max(Sum), .groups = "keep") %>%
    ungroup %>%
    select(Total, Sum)
  # ADDED ON JULY 23, 2020; REMOVED ON DECEMBER 30, 2020
  # lastPoint <- extremePoints %>%
  #   slice(nrow(extremePoints))
  # if (lastPoint$Total != lastPoint$Sum) {
  #   lastPoint %<>%
  #     mutate(Total = Sum)
  #   extremePoints %<>%
  #     bind_rows(lastPoint)
  # }
  if (CPP) {
    points <- extremePoints %>%
      mutate_at("Sum", ~{add(., 1/2)})
    if (nrow(points) <= 2) {
      breakpoints <- points
    } else {
      write_csv(points, fname)
      cmdName <- paste0(ifelse(HPC, "/hpc/grid/wip_cmg_systems-immunology/chindl/", "/Users/lchindelevitch/"), 
                        paste0(ifelse(HPC, "", "Downloads/"), "ReverseGWAS/"), "rgwas/plfoptq")
      breakpoints <- system2(cmdName, args = c(fname, width/2), stdout = TRUE) %>%
        str_split_fixed(",", 2) %>%
        set_colnames(c("Total", "Sum")) %>%
        as_tibble() %>%
        mutate_all(as.numeric)
    }
  } else {
    extremePoints %<>%
      removeCollinear %>%
      mutate_at("Sum", ~{add(., 1/2)})
  }
  pointsP <- extremePoints
  if (!dryRun) {
    pointsM <- extremePoints
    pointsP %<>%
      mutate(Sum = Sum + width/2)
    pointsM %<>%
      mutate(Sum = Sum - width/2)
    if (!CPP) {
      breakpoints <- ImaiIriAlgorithm(pointsP, pointsM)
    }
  } else {
    breakpoints <- extremePoints
  }
  if (verify && !dryRun && nrow(breakpoints) > 0) { ### check that the required bounds are satisfied by the output
    print("Verification of correctness")
    pointsL <- pointsM %>%
      mutate(Sum = round(Sum))
    pointsU <- pointsP %>%
      mutate(Sum = round(Sum))
    stopifnot(checkPosition(breakpoints, pointsL, below = TRUE))
    stopifnot(checkPosition(breakpoints, pointsU, below = FALSE))
  }
  if (nrow(breakpoints) > 0) {
    if (min(breakpoints$Sum) != firstPoint$Sum) {
      breakpoints %<>%
        bind_rows(firstPoint) %>%
        arrange(Total, Sum)
    }
  }
  output <- list(breakpoints = breakpoints, inputSize = nrow(pointsP), input = pointsP)
  output
}

### Checking function to ensure that the input endpoints define a boundary lying entirely above the PWLF defined by the checkpoints
### If below = TRUE, check that the boundary lies below the PWLF instead (abbreviation used: PWLF = piecewise linear function).
### If reportCoincident = FALSE, also returns a list of all checkpoints that coincide exactly with one of the endpoints.
### Returns TRUE if the condition is fulfilled, otherwise prints out the details of the earliest found violation, returns FALSE.
checkPosition = function(endpoints, checkpoints, below = TRUE, reportCoincident = FALSE) {
  checkpoints %<>%
    as.matrix
  endpoints %<>%
    as.matrix
  stopifnot(all(near(range(endpoints[,1]), range(checkpoints[,1]), tol = TOL)))
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
        add(1)
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
