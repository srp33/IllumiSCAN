scanNorm <- function(signalExprData, signalProbeSequences, signalPValueData=NULL, controlExprData=NULL, convThreshold=0.5, intervalN=10000, binsize=500, nbins=25, maxIt=100, asUPC=FALSE, numCores=1, verbose=FALSE)
{
  #############################################
  # Check parameters
  #############################################

  if (!is.matrix(signalExprData))
    stop("signalExprData must be a matrix.")
  if (!is.vector(signalProbeSequences))
    stop("signalProbeSequences must be a vector")
  if (nrow(signalExprData) != length(signalProbeSequences))
    stop("The dimensions of signalExprData and signalProbeSequences are incompatible.")
  if (!is.null(signalPValueData))
  {
    if (!is.matrix(signalPValueData))
      stop("signalPValueData must be a matrix.")

    if (all(dim(signalExprData) != dim(signalPValueData)))
      stop("The dimensions of signalExprData and signalPValueData are incompatible.")
  }

  if (!is.null(controlExprData))
  {
    if (!is.matrix(controlExprData))
      stop("controlExprData must be a matrix.")

    if (ncol(signalExprData) != ncol(controlExprData))
      stop("The dimensions of signalExprData and controlExprData are incompatible.")
  }
  
  #############################################
  # Perform background correction, if necessary
  #############################################

  exprData <- doLog2(signalExprData)
  
  if (all(signalExprData > 0)) { # No background subtraction has been performed
    if (is.null(controlExprData)) {
      if (!is.null(signalPValueData))
        exprData <- nec(x = doLog2(signalExprData), detection.p = signalPValueData)
    } else {
        status <- c(rep("regular", nrow(signalExprData)), rep("negative", nrow(controlExprData)))

        combinedExprData <- rbind(signalExprData, controlExprData)

        exprData <- nec(x = doLog2(combinedExprData), status = status)
        exprData <- exprData[1:nrow(signalExprData),,drop=FALSE]

##      exprData <- doLog2(rbind(signalExprData, controlExprData))
    }
  }

  #############################################
  # Normalize
  #############################################

  mx = buildDesignMatrix(signalProbeSequences)

  numSamples <- ncol(exprData)
  
  if (numSamples > 1 & numCores > 1)
  {
    cl <- makeCluster(numCores, outfile="")
    registerDoParallel(cl)
  }
  
  if (numSamples == 1)
  {
    normData <- as.matrix(scanNormVector(exprData[,1], mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose))
  }
  else
  {
    normData <- foreach(i = 1:ncol(exprData), .combine = cbind, .export=c("scanNormVector", "sampleProbeIndices", "EM_vMix", "mybeta", "assign_bin", "vsig", "vresp", "dn", "vbeta", "sig")) %dopar%
    {
      scanNormVector(exprData[,i], mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose)
    }
  }

  if (numSamples > 1 & numCores > 1)
    stopCluster(cl)
  
  rownames(normData) <- rownames(exprData)
  colnames(normData) <- colnames(exprData)
  
  ##########################
  normData <- normData[rownames(signalExprData),,drop=FALSE]

  return(normData)
}

scanNormVector <- function(my, mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose)
{
  # Add a tiny amount of random noise
  set.seed(0)
  noise = rnorm(length(my)) / 10000000
  my = my + noise

  nGroups = floor(length(my) / binsize)
  samplingProbeIndices = sampleProbeIndices(total=length(my), intervalN=intervalN, verbose=verbose)

  mixResult = EM_vMix(y=my[samplingProbeIndices], X=mx[samplingProbeIndices,], nbins=nbins, convThreshold=convThreshold, maxIt=maxIt, verbose=verbose)

  m1 = mx %*% mixResult$b1
  m2 = mx %*% mixResult$b2

  index = order(m1)
  y_norm = rep(0, length(my))
  for (i in 1:nGroups)
  {
    tmp = index[(binsize * i):min(binsize * i + binsize, length(my))]
    tmpSd = as.vector(sig(y=my[tmp], m=m1[tmp], verbose=verbose))
    y_norm[tmp] = ((my[tmp] - m1[tmp]) / tmpSd)
  }

  bin = assign_bin(y=m1, nbins=nbins, verbose=verbose)
  gam = vresp(y=my, X=mx, bin=bin, p=mixResult$p, b1=mixResult$b1, s1=mixResult$s1, b2=mixResult$b2, s2=mixResult$s2, verbose=verbose)[,2]

  y_norm = round(y_norm, 8)
  gam = round(gam, 8)

  if (asUPC)
  {
    return(gam)
  } else {
    return(y_norm)
  }
}

doLog2 <- function(x)
{
  # This is a semi-crude way of checking whether the values were not previously log-transformed
  if (max(x) > 100)
    x <- log2(x)
  
  return(x)
}

buildDesignMatrix = function(seqs, verbose=TRUE)
{
  mx = sequenceDesignMatrix(seqs)

  numA = apply(mx[,which(grepl("^A_", colnames(mx)))], 1, sum)
  numC = apply(mx[,which(grepl("^C_", colnames(mx)))], 1, sum)
  numG = apply(mx[,which(grepl("^G_", colnames(mx)))], 1, sum)
  numT = 60 - (numA + numC + numG)

  mx = cbind(numT, mx, numA^2, numC^2, numG^2, numT^2)
  #mx = cbind(numA, numC, numG, numA^2, numC^2, numG^2, numT^2)
  #mx = cbind(numA, numC, numG)
  mx = apply(mx, 2, as.integer)  

  return(mx)
}

sampleProbeIndices = function(total, intervalN, verbose=TRUE)
{
  interval = floor(total / intervalN)
  if (interval <= 1)
    interval = 1

  seq(1, total, interval)
}

EM_vMix = function(y, X, nbins, convThreshold=.01, maxIt=100, verbose=TRUE)
{
  if (verbose)
    message("Starting EM")

  quan = sort(y)[floor(0.5 * length(y)) - 1]
  gam = cbind(as.integer(y <= quan), as.integer(y > quan))

  p = apply(gam, 2, mean)

  b1 = mybeta(y=y, X=X, gam=gam[,1], verbose=verbose)
  b2 = mybeta(y=y, X=X, gam=gam[,2], verbose=verbose)
  bin = assign_bin(y=y, nbins=nbins, verbose=verbose)
  s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
  s2 = vsig(y=y, X=X, b=b2, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)

  theta_old=c(p, b1, s1, b2, s2)

  it = 0
  conv = 1000000

  while (conv > convThreshold & it < maxIt)
  {
    # Expectation Step:
    gam = vresp(y=y, X=X, bin=bin, p=p, b1=b1, s1=s1, b2=b2, s2=s2, verbose=verbose)

    #M-Step
    p = apply(gam, 2, mean)
    b1 = vbeta(y=y, X=X, bin=bin, gam=gam[,1], s2=s1, prof=TRUE, verbose=verbose)
    bin = assign_bin(y=(X %*% b1), nbins=nbins, verbose=verbose)
    b2 = vbeta(y=y, X=X, bin=bin, gam=gam[,2], s2=s2, prof=FALSE, verbose=verbose)
    s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
    s2 = vsig(y=y, X=X, b=b2, gam=gam[,2], bin=bin, nbins=nbins, verbose=verbose)

    theta = c(p, b1, s1, b2, s2)
    conv = max(abs(theta - theta_old) / theta_old)
    theta_old = theta
    it = it + 1

    if (verbose)
      message("Attempting to converge...iteration ", it, ", c = ", round(conv, 6))
  }

  if (verbose)
  {
    if (it == maxIt)
    {
      message("Reached convergence limit...", it, " iterations. Proportion of background probes: ", round(p[1], 6))
    } else {
      message("Converged in ", it, " iterations. Proportion of background probes: ", round(p[1], 6))
    }
  }

  list(p=p, b1=b1, b2=b2, s1=s1, s2=s2, bin=bin)
}

mybeta = function(y, X, gam, verbose=TRUE)
{
  sqgam = sqrt(gam)
  Xw = sqgam * X
  yw = sqgam * y

  z = t(Xw) %*% Xw
  a = solve(z)

  b = a %*% t(Xw)
  as.numeric(b %*% yw)
}

assign_bin = function(y, nbins, verbose=TRUE)
{
  quans = sort(y)[floor(length(y) * 1:nbins / nbins)]
  bins = sapply(y, function(x) { sum(x>quans) }) + 1

  if (length(table(bins)) != nbins)
  {
    if (verbose)
      message("The values were not separated into enough bins, so a tiny amount of noise will be added to make this possible.")

    set.seed(1)
    noise = rnorm(length(y)) / 10000000
    bins = assign_bin(y + noise, nbins, verbose)
  }

  bins
}

vsig = function(y, X, b, gam, bin, nbins, verbose=TRUE)
{
  s2 = NULL

  for (i in 1:nbins)
  {
    ystar = y[bin==i]
    Xstar = X[bin==i,]
    gamstar = gam[bin==i] + .01
    resid = as.numeric(ystar - Xstar %*% b)

    s2 = c(s2, ((resid * gamstar) %*% resid) / sum(gamstar))
  }

  s2
}

vresp = function(y, X, bin, p, b1, s1, b2, s2, verbose=TRUE)
{
  vars0 = s1[bin]
  L0 = dn(y=y, m=(X %*% b1), s2=vars0, verbose=verbose)
  vars1 = s2[bin]
  L1 = dn(y=y, m=(X %*% b2), s2=vars1, verbose=verbose)

  gam1 = p[1] * L0 / (p[1] * L0 + p[2] * L1)
  gam2 = 1 - gam1
  cbind(gam1, gam2)
}

dn = function(y, m, s2, verbose=TRUE)
{
  1 / (sqrt(2 * pi * s2)) * exp(-1 / (2 * s2) * (y - m)^2)
}

vbeta = function(y, X, bin, gam, s2, prof, verbose=TRUE)
{
  vars = sqrt(s2[bin])
  sqgam = sqrt(gam)
  vars_sqgam = vars * sqgam

  Xw = 1 / vars * sqgam * X
  yw = 1 / vars * sqgam * y

  tXw = t(Xw)
  tXwXw = tXw %*% Xw
  stXwXw = solve(tXwXw)
  stXwXwtXw = stXwXw %*% tXw
  result = stXwXwtXw %*% yw

  result
}

sig = function(y, m, verbose=TRUE)
{
  resid = y - m
  sqrt((resid %*% resid) / length(y))
}