
exactTestDoubleTail <- function (y1, y2, dispersion = 0, big.count = 900) 
{
  ntags <- NROW(y1)
  n1 <- NCOL(y1)
  n2 <- NCOL(y2)
  if (n1 > 1) 
    s1 <- round(rowSums(y1))
  else s1 <- round(y1)
  if (n2 > 1) 
    s2 <- round(rowSums(y2))
  else s2 <- round(y2)
  if (length(dispersion) == 1) 
    dispersion <- rep(dispersion, ntags)
  s <- s1 + s2  # sum of rowSums
  mu <- s/(n1 + n2)  # normalize by total # of columns
  mu1 <- n1 * mu  # Same as n1/(n1+n2) * (sum of rowSums)
  mu2 <- n2 * mu
  pvals <- rep(1, ntags)
  names(pvals) <- names(y1)
  pois <- dispersion <= 0
  if (any(pois)) 
    pvals[pois] <- binomTest(s1[pois], s2[pois], p = n1/(n1 + 
                                                           n2))
  big <- s1 > big.count & s2 > big.count
  if (any(big)) {
    y1 <- as.matrix(y1)
    y2 <- as.matrix(y2)
    pvals[big] <- exactTestBetaApprox(y1[big, , drop = FALSE], 
                                      y2[big, , drop = FALSE], dispersion[big])
  }
  p.bot <- size1 <- size2 <- rep(0, ntags)
  left <- s1 < mu1 & !pois & !big
  if (any(left)) {
    p.bot[left] <- dnbinom(s[left], size = (n1 + n2)/dispersion[left], 
                           mu = s[left])
    size1[left] <- n1/dispersion[left]
    size2[left] <- n2/dispersion[left]
    for (g in which(left)) {
      x <- 0:s1[g]
      p.top <- dnbinom(x, size = size1[g], mu = mu1[g]) * 
        dnbinom(s[g] - x, size = size2[g], mu = mu2[g])
      pvals[g] <- 2 * sum(p.top)
    }
    pvals[left] <- pvals[left]/p.bot[left]
  }
  right <- s1 > mu1 & !pois & !big
  if (any(right)) {
    p.bot[right] <- dnbinom(s[right], size = (n1 + n2)/dispersion[right], 
                            mu = s[right])
    size1[right] <- n1/dispersion[right]
    size2[right] <- n2/dispersion[right]
    for (g in which(right)) {
      x <- s1[g]:s[g]
      p.top <- dnbinom(x, size = size1[g], mu = mu1[g]) * 
        dnbinom(s[g] - x, size = size2[g], mu = mu2[g])
      pvals[g] <- 2 * sum(p.top)
    }
    pvals[right] <- pvals[right]/p.bot[right]
  }
  pmin(pvals, 1)
}



exactTest <- function (object, pair = 1:2, dispersion = "auto", rejection.region = "doubletail", 
                       big.count = 900, prior.count = 0.125) 
{
  if (!is(object, "DGEList")) 
    stop("Currently only supports DGEList objects as the object argument.")
  if (length(pair) != 2) 
    stop("Pair must be of length 2.")
  rejection.region <- match.arg(rejection.region, c("doubletail", 
                                                    "deviance", "smallp"))
  group <- as.factor(object$samples$group)
  levs.group <- levels(group)
  if (is.numeric(pair)) 
    pair <- levs.group[pair]
  else pair <- as.character(pair)
  if (!all(pair %in% levs.group)) 
    stop("At least one element of given pair is not a group.\n Groups are: ", 
         paste(levs.group, collapse = " "))
  if (is.null(dispersion)) 
    dispersion <- "auto"
  if (is.character(dispersion)) {
    dispersion <- match.arg(dispersion, c("auto", "common", 
                                          "trended", "tagwise"))
    dispersion <- switch(dispersion, common = object$common.dispersion, 
                         trended = object$trended.dispersion, tagwise = object$tagwise.dispersion, 
                         auto = getDispersion(object))
    if (is.null(dispersion)) 
      stop("specified dispersion not found in object")
    if (is.na(dispersion[1])) 
      stop("dispersion is NA")
  }
  ldisp <- length(dispersion)
  ntags <- nrow(object$counts)
  if (ldisp != 1 && ldisp != ntags) 
    stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
  if (ldisp == 1) 
    dispersion <- rep(dispersion, ntags)
  group <- as.character(group)
  j <- group %in% pair
  y <- object$counts[, j, drop = FALSE]
  lib.size <- object$samples$lib.size[j]
  norm.factors <- object$samples$norm.factors[j]
  group <- group[j]
  if (is.null(rownames(y))) 
    rownames(y) <- paste("tag", 1:ntags, sep = ".")
  lib.size <- lib.size * norm.factors
  offset <- log(lib.size)
  lib.size.average <- exp(mean(offset))
  prior.count <- prior.count * lib.size/mean(lib.size)
  offset.aug <- log(lib.size + 2 * prior.count)
  j1 <- group == pair[1]
  n1 <- sum(j1)
  if (n1 == 0) 
    stop("No libraries for", pair[1])
  y1 <- y[, j1, drop = FALSE]
  abundance1 <- mglmOneGroup(y1 + matrix(prior.count[j1], 
                                         ntags, n1, byrow = TRUE), offset = offset.aug[j1], dispersion = dispersion)
  j2 <- group == pair[2]
  n2 <- sum(j2)
  if (n1 == 0) 
    stop("No libraries for", pair[2])
  y2 <- y[, j2, drop = FALSE]
  abundance2 <- mglmOneGroup(y2 + matrix(prior.count[j2], 
                                         ntags, n2, byrow = TRUE), offset = offset.aug[j2], dispersion = dispersion)
  logFC <- (abundance2 - abundance1)/log(2)
  abundance <- mglmOneGroup(y, dispersion = dispersion, offset = offset)
  e <- exp(abundance)
  input.mean <- matrix(e, ntags, n1)
  output.mean <- input.mean * lib.size.average
  input.mean <- t(t(input.mean) * lib.size[j1])
  y1 <- q2qnbinom(y1, input.mean = input.mean, output.mean = output.mean, 
                  dispersion = dispersion)
  input.mean <- matrix(e, ntags, n2)
  output.mean <- input.mean * lib.size.average
  input.mean <- t(t(input.mean) * lib.size[j2])
  y2 <- q2qnbinom(y2, input.mean = input.mean, output.mean = output.mean, 
                  dispersion = dispersion)
  exact.pvals <- switch(rejection.region, doubletail = exactTestDoubleTail(y1, 
                                                                           y2, dispersion = dispersion, big.count = big.count), 
                        deviance = exactTestByDeviance(y1, y2, dispersion = dispersion, 
                                                       big.count = big.count), smallp = exactTestBySmallP(y1, 
                                                                                                          y2, dispersion = dispersion, big.count = big.count))
  AveLogCPM <- object$AveLogCPM
  if (is.null(AveLogCPM)) 
    AveLogCPM <- aveLogCPM(object)
  de.out <- data.frame(logFC = logFC, logCPM = AveLogCPM, 
                       PValue = exact.pvals)
  rn <- rownames(object$counts)
  if (!is.null(rn)) 
    rownames(de.out) <- make.unique(rn)
  new("DGEExact", list(table = de.out, comparison = pair, 
                       genes = object$genes))
}



et <- exactTest(y)


