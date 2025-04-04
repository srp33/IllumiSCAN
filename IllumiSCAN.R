# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# for (packageName in c("httr", "limma", "oligo", "doParallel" "illuminaHumanv4.db")) {
#   BiocManager::install(packageName)
# }

library(doParallel)
library(httr)
library(limma)
library(oligo)

normalizeBeadChipData <- function(filePath, species="Human", platformVersion="4", probeIDColumn="ID_REF", exprColumnPattern="GSM", detectionPValueColumnPattern="Detection Pval", numCores=1, verbose=FALSE) {
  # TODO: Make sure species and platform are valid.

  platformPackagePrefix <- paste0("illumina", species, "v", platformVersion)
  platformPackageName <- paste0(platformPackagePrefix, ".db")
  message(paste0("Loading package ", platformPackageName))
  
  library(platformPackageName, character.only=TRUE)
  
  getRefInfo <- function(platformPackagePrefix, suffix) {
    info <- get(paste0(platformPackagePrefix, suffix))
    info <- info[mappedkeys(info)] # Returns the subset of mapped keys.
    return(as.data.frame(info))
  }

  # https://bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
  probeSequenceRef <- getRefInfo(platformPackagePrefix, "PROBESEQUENCE")
  # probeReporterRef <- getRefInfo(platformPackagePrefix, "REPORTERGROUPID") # Control, housekeeping, permuted probes
  probeQualityRef <- getRefInfo(platformPackagePrefix, "PROBEQUALITY")
  probeGeneRef <- getRefInfo(platformPackagePrefix, "ENSEMBLREANNOTATED")
  
  # Remove low-quality probes
  perfectProbes = which(grepl("Perfect", probeQualityRef$ProbeQuality))
  goodProbes = which(grepl("Good", probeQualityRef$ProbeQuality))
  probesToKeep = probeQualityRef$IlluminaID[c(perfectProbes, goodProbes)]

  # Keep probes mapped to genes
  probesToKeep = intersect(probesToKeep, probeGeneRef$IlluminaID)

  nonNormData <- read.ilmn(filePath, probeid = probeIDColumn, expr = exprColumnPattern, other.columns = detectionPValueColumnPattern)
  
  exprData <- nonNormData$E
  detectionPValues <- nonNormData$other$`Detection Pval`
  
  colnames(exprData) <- paste0(exprColumnPattern, colnames(exprData))
  colnames(detectionPValues) <- paste0(exprColumnPattern, colnames(detectionPValues))

  probesToKeep <- intersect(probesToKeep, rownames(exprData))
  exprData <- exprData[probesToKeep,]
  detectionPValues <- detectionPValues[probesToKeep,]
  
  rownames(probeSequenceRef) = probeSequenceRef$IlluminaID
  probeSequences = probeSequenceRef[probesToKeep, 2]
  
  normData <- backgroundCorrect(exprData, signalPValueData=detectionPValues)
  
  normData <- scanNorm(exprData, probeSequences, numCores=numCores, verbose=verbose)

  rownames(probeGeneRef) = probeGeneRef$IlluminaID
  probeGenes = probeGeneRef[probesToKeep, 2]
  # TODO: Add an argument to summarize at the gene level and then implement the logic for this.

  return(normData)
}

getNonNormalizedDataFromGEO <- function(gseID) {
  nonNormalizedURL = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gseID, "&format=file&file=", gseID, "%5Fnon%5Fnormalized%2Etxt%2Egz")
  
  downloadDirPath = paste0(tempdir(), "/", gseID)
  
  if (!dir.exists(downloadDirPath)) {
    dir.create(downloadDirPath)
  }
  
  downloadFilePath = paste0(downloadDirPath, "/data.txt.gz")
  
  if (file.exists(downloadFilePath)) {
    return(downloadFilePath)
  }
  
  response <- HEAD(nonNormalizedURL)
  
  if (status_code(response) == 200) {
    download.file(nonNormalizedURL, destfile = downloadFilePath, mode = "wb")
    
    return(downloadFilePath)
  } else {
    stop(paste0("A non-normalized GEO file could not be found using the standard URL for ", gseID, ". You will need to provide the URL."))
  }
}

backgroundCorrect <- function(signalExprData, controlExprData=NULL, signalPValueData=NULL)
{
  #############################################
  # Check parameters
  #############################################
  
  if (!is.matrix(signalExprData))
    stop("signalExprData must be a matrix.")
  
  if (is.null(signalPValueData)) {
    if (is.null(controlExprData)) {
      stop("If signalPValueData is not provided, controlExprData must be provided.")
    } else {
      if (!is.matrix(controlExprData))
        stop("controlExprData must be a matrix.")
      
      if (ncol(signalExprData) != ncol(controlExprData))
        stop("The dimensions of signalExprData and controlExprData are incompatible.")
    }
  } else {
    if (!is.matrix(signalPValueData))
      stop("signalPValueData must be a matrix.")
    
    if (all(dim(signalExprData) != dim(signalPValueData)))
      stop("The dimensions of signalExprData and signalPValueData are incompatible.")
  }
  
  if (any(signalExprData < 0)) {
    stop("It appears that the values in signalExprData have already been background corrected. Please use raw, probe-level expression intensities.")
  }
  
  #############################################
  # Perform background correction (limma).
  #   It models the observed signal as a combination of: Observed_Signal = True_Signal + Background_Noise
  #   It removes non-specific fluorescence, helping improve the accuracy of low-intensity probe measurements.
  #   https://academic.oup.com/nar/article/38/22/e204/1049223
  #############################################

  status <- rep("regular", nrow(signalExprData))
  
  if (is.null(controlExprData)) { # We have detection p-values only.
    exprData <- nec(x = signalExprData, status = status, detection.p = signalPValueData)
  } else { # We have control values.
    status <- c(status, rep("negative", nrow(controlExprData)))

    exprData <- nec(x = rbind(signalExprData, controlExprData), status = status)
    # signalCorrectedData <- exprData[1:nrow(signalExprData),]
    # controlCorrectedData <- exprData[(nrow(signalExprData) + 1):(nrow(signalExprData) + nrow(controlExprData)),]
  }
  
  return(exprData)
}

# intervalN: Interval for probe sampling, Controls how many probes are selected for estimating model parameters.
#            Instead of using every single probe (which is computationally expensive), the function samples probes at intervals.
#            Smaller values → More probes are used, leading to better accuracy but slower computation.
#            Larger values → Fewer probes are used, making computation faster but potentially less precise.
# binsize:   Size of bins for intensity normalization
#            In SCAN normalization, probes are grouped into bins based on their intensity.
#            binsize controls how many probes go into each bin.
#            Smaller bins = More precise normalization but slower computation.
#            Larger bins = Smoother adjustments but may miss finer details.
# nbins:     Number of bins
#            Determines how many bins are used for normalizing intensity distributions.
#            Instead of normalizing each probe individually, the algorithm divides probes into nbins groups based on their signal intensity.
#            More bins = Finer-grained correction but slower computation.
#            Fewer bins = More general correction, which may not capture subtle effects.
# maxIt:     Maximum number of iterations
#            Controls how many times the EM algorithm is allowed to run before stopping.
#            If the algorithm has not converged after maxIt iterations, it forces an early stop.
#            Higher values allow EM to run longer, improving accuracy.
#            Lower values can speed up processing but may result in incomplete convergence.
# asUPC:     Whether to return values as Universal exPression Codes (probabilistic indicators of expression).
# numCores:  The number of CPU cores to use when processing the data.
# verbose:   Whether to display verbose output when processing the data.

scanNorm <- function(exprData, probeSequences, convThreshold=0.5, intervalN=10000, binsize=500, nbins=25, maxIt=100, asUPC=FALSE, numCores=1, verbose=FALSE)
{
  #############################################
  # Check parameters
  #############################################

  if (!is.matrix(exprData))
    stop("exprData must be a matrix.")
  if (!is.vector(probeSequences))
    stop("probeSequences must be a vector")
  if (nrow(exprData) != length(probeSequences))
    stop("The dimensions of exprData and probeSequences must be identical.")
  
  #TODO: Make sure other parameters are in valid range.
  
  #############################################
  # SCAN normalize
  #############################################
  
  mx = buildDesignMatrix(probeSequences)
  
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
    normData <- foreach(i = 1:ncol(exprData), .combine = cbind, .export=c("scanNormVector", "doLog2", "sampleProbeIndices", "EM_vMix", "mybeta", "assign_bin", "vsig", "vresp", "dn", "vbeta", "sig")) %dopar%
    {
      scanNormVector(exprData[,i], mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose)
    }
  }
  
  if (numSamples > 1 & numCores > 1)
    stopCluster(cl)
  
  rownames(normData) <- rownames(exprData)
  colnames(normData) <- colnames(exprData)
  
  ##########################
  normData <- normData[rownames(exprData),,drop=FALSE]
  
  return(normData)
}

scanNormVector <- function(my, mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose)
{
  # Log transform the data, if needed.
  my = doLog2(my)

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
  y_norm = 2^y_norm # Reverse log2 transformation.
  
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
  # This is a semi-crude way of checking whether the values were not previously log-transformed.
  if (max(x) > 30)
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
  #mx = cbind(numT, mx, as.integer(numA^2), as.integer(numC^2), as.integer(numG^2), as.integer(numT^2))
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