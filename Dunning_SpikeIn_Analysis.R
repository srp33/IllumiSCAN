# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("broom", "doParallel", "limma", "oligo", "tidyverse"))

library(broom)
library(doParallel)
library(limma)
library(oligo)
library(tidyverse)

readProbeData <- function(filePath, signalColumnSuffix, pValueColumnSuffix)
{
  # Read probe data
  probeData <- suppressMessages(read_tsv(filePath))

  # Keep columns of interest
  exprData <- probeData[,c(2, which(grepl(paste(signalColumnSuffix, "$", sep=""), colnames(probeData)))), drop=FALSE]
  colnames(exprData) <- sub(paste(".", signalColumnSuffix, sep=""), "", colnames(exprData))
  
  pValueData <- probeData[,c(2, which(grepl(paste(pValueColumnSuffix, "$", sep=""), colnames(probeData)))), drop=FALSE]
  colnames(pValueData) <- sub(paste(".", pValueColumnSuffix, sep=""), "", colnames(pValueData))

  return(list(exprData=exprData, pValueData=pValueData))
}

# From Ringo/Bioconductor package
compute.gc <- function(probe.sequences, digits=2)
{
  stopifnot(is.character(probe.sequences))

  splitted.seqs <- strsplit(toupper(probe.sequences),split="")
  round(sapply(splitted.seqs, function(x) length(grep("[GC]",x)))/
  listLen(splitted.seqs), digits=digits)
}

tmpDir <- tempdir()
zipUrl <- "https://github.com/markdunning/statistical-issues-illumina-microarray/raw/master/spike_beadstudio_output.zip"
zipFilePath <- paste(tmpDir, "/spike_beadstudio_output.zip", sep="")
if (!file.exists(zipFilePath)) {
  download.file(zipUrl, zipFilePath)
}

sampleDataFilePath <- unz(zipFilePath, "SampleProbeProfile.txt")
controlDataFilePath <- unz(zipFilePath, "ControlProbeProfile.txt")

sampleAnnotationFilePath <- "https://raw.githubusercontent.com/markdunning/statistical-issues-illumina-microarray/master/NewAnnotationM1.txt"
controlAnnotationFilePath <- "https://raw.githubusercontent.com/markdunning/statistical-issues-illumina-microarray/master/ControlSequences.txt"
spikeInProfileFilePath <- "https://raw.githubusercontent.com/markdunning/statistical-issues-illumina-microarray/master/spikeins_profile.csv"
targetsFilePath <- "https://github.com/markdunning/statistical-issues-illumina-microarray/raw/master/spike_targets.txt"

sampleProbeData <- readProbeData(sampleDataFilePath, "AVG_Signal", "Detection Pval")
controlProbeData <- readProbeData(controlDataFilePath, "AVG_Signal", "Detection Pval") # This also contains the spike-in data

sampleAnnotationData <- suppressMessages(read_tsv(sampleAnnotationFilePath)) %>%
  dplyr::select(ProbeId0, Sequence) %>%
  dplyr::rename(ProbeID=ProbeId0)

controlAnnotationData <- suppressMessages(read_tsv(controlAnnotationFilePath)) %>%
  dplyr::rename(Sequence=`Seq(5'to3')`, ProbeID=Common_Olicode_Name, ProbeType=Reporter_Group_Name) %>%
  dplyr::select(Sequence, ProbeID, ProbeType) %>%
  dplyr::filter(ProbeType=="negative") %>%
  dplyr::select(-ProbeType) %>%
  distinct()

# There are a few control probes (biotin, housekeeping genes, etc.) that are present
# in the control annotations and sample annotations. But expression values for these
# are only present in sampleProbeData, so we need to account for this. We will use
# the annotations in controlAnnotationData for these probes because this file
# indicates what type of control probes they are.
overlappingProbes <- intersect(sampleAnnotationData$ProbeID, controlAnnotationData$ProbeID)
sampleAnnotationData <- dplyr::filter(sampleAnnotationData, !ProbeID %in% overlappingProbes)

spikeInAnnotationData <- suppressMessages(read_csv(spikeInProfileFilePath)) %>%
  dplyr::select(ProbeID, Sequence, TargetID)
spikeInAnnotationData$TargetID <- sapply(spikeInAnnotationData$TargetID, function(x) {strsplit(x, "_")[[1]][1]})
spikeInAnnotationData$TargetID <- factor(spikeInAnnotationData$TargetID)

spikeInProbeIDs <- as.character(pull(spikeInAnnotationData, ProbeID))

allAnnotationData <- rbind(sampleAnnotationData, dplyr::select(controlAnnotationData, ProbeID, Sequence), dplyr::select(spikeInAnnotationData, ProbeID, Sequence)) %>%
  distinct()

signalExprData <- dplyr::filter(sampleProbeData$exprData, ProbeID %in% allAnnotationData$ProbeID)
controlExprData <- dplyr::filter(controlProbeData$exprData, ProbeID %in% allAnnotationData$ProbeID)
signalPValueData <- dplyr::filter(sampleProbeData$pValueData, ProbeID %in% allAnnotationData$ProbeID)
controlPValueData <- dplyr::filter(controlProbeData$pValueData, ProbeID %in% allAnnotationData$ProbeID)

# Move the spike-in data from the control data to the signal data
signalExprData <- rbind(signalExprData, dplyr::filter(controlExprData, ProbeID %in% spikeInProbeIDs))
controlExprData <- dplyr::filter(controlExprData, !(ProbeID %in% spikeInProbeIDs))
signalPValueData <- rbind(signalPValueData, dplyr::filter(controlPValueData, ProbeID %in% spikeInProbeIDs))

signalExprData <- inner_join(allAnnotationData, signalExprData, by="ProbeID")
controlExprData <- inner_join(allAnnotationData, controlExprData, by="ProbeID")
signalPValueData <- inner_join(allAnnotationData, signalPValueData, by="ProbeID")

signalProbeIDs <- pull(signalExprData, ProbeID)
signalProbeSequences <- pull(signalExprData, Sequence)
signalExprData <- as.matrix(dplyr::select(signalExprData, -ProbeID, -Sequence))
signalPValueData <- as.matrix(dplyr::select(signalPValueData, -ProbeID, -Sequence))
rownames(signalExprData) <- signalProbeIDs
rownames(signalPValueData) <- signalProbeIDs

controlProbeIDs <- pull(controlExprData, ProbeID)
controlProbeSequences <- pull(controlExprData, Sequence)
controlExprData <- as.matrix(dplyr::select(controlExprData, -ProbeID, -Sequence))
rownames(controlExprData) <- controlProbeIDs

signalProbeSequenceGC <- compute.gc(signalProbeSequences, digits=2)
signalProbeSequenceGC <- round(signalProbeSequenceGC * 2, 1) / 2

targetData <- suppressMessages(read_delim(targetsFilePath, delim=" "))
targetData <- targetData %>% dplyr::select(ArrayNo, SpikeConc) %>% dplyr::rename(SampleID=ArrayNo)
targetData$SampleID <- sapply(targetData$SampleID, function(x) {paste(strsplit(x, "_")[[1]][1:2], collapse="_")})
spikeLevels <- c(0e+00, 1e-02, 3e-02, 1e-01, 3e-01, 1e+00, 3e+00, 1e+01, 3e+01, 1e+02, 3e+02, 1e+03)
#targetData$SpikeConc <- factor(targetData$SpikeConc, levels=spikeLevels)
targetData <- distinct(targetData)

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

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
  # Add some tiny random noise
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

#####################################################################

restructureConsistencyData <- function(exprData, outFilePath) {
  normSpikeInData <- exprData[spikeInProbeIDs,]
  normSpikeInData <- bind_cols(tibble(ProbeID=spikeInProbeIDs), normSpikeInData)
  normSpikeInData <- pivot_longer(normSpikeInData, -ProbeID, names_to = "SampleID", values_to = "Value")
  normSpikeInData <- inner_join(normSpikeInData, targetData, by="SampleID")
  
  evalData <- spikeInAnnotationData %>%
    dplyr::select(ProbeID, TargetID)
  evalData$ProbeID <- factor(as.character(evalData$ProbeID))
  evalData <- inner_join(evalData, normSpikeInData, by="ProbeID")
  
  write_tsv(evalData, outFilePath)
  return(evalData)
}

plotConsistency <- function(evalData, metric, description, outFilePath) {
  metric <- bquote(R^2 == .(round(metric, 2)))
  
  p <- mutate(evalData, SpikeConc = factor(SpikeConc, levels = spikeLevels)) %>%
    ggplot(aes(x = SpikeConc, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.08) +
    annotate("text", x = 10.5, y = 0.5, label = metric, size = 3, color = "red") +
    ggtitle(description) +
    theme_bw()
  
  print(p)
  
  ggsave(outFilePath, width = 8, height = 6)
}

calcConsistencyMetric <- function(evalData)
{
  evalFit <- lm(Value~SpikeConc, data=evalData)

  return(summary(evalFit)$r.squared)
}

# convThresholdOptions <- c(10)
# intervalNOptions <- c(1000)
# binsizeOptions <- c(50)

convThresholdOptions <- c(0.01, 0.1, 1, 10)
intervalNOptions <- c(1000, 10000, 20000, 50000)
binsizeOptions <- c(50, 500, 5000)

#TODO: Plot the summarized results? Or just examine it as a table?

numCores = 4
verbose = FALSE

paramCombos <- expand.grid(convThresholdOptions, intervalNOptions, binsizeOptions)
colnames(paramCombos) <- c("convThreshold", "intervalN", "binsize")

paramTuningOutFilePath <- "Param_Tuning_Results.tsv"

if (!file.exists(paramTuningOutFilePath))
{
  comparisonResults <- NULL

  for (i in 1:nrow(paramCombos))
  {
    print(paste("Executing parameter combination when using control data ", i, "...", sep=""))
    convThreshold <- paramCombos[i,1]
    intervalN <- paramCombos[i,2]
    binsize <- paramCombos[i,3]
  
    normData <- scanNorm(signalExprData, signalProbeSequences, controlExprData = controlExprData, convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, numCores=numCores, verbose=verbose)

    outFilePrefix <- paste0("Figures/", convThreshold, "_", intervalN, "_", binsize, "_NoControls")

    evalData <- restructureConsistencyData(normData, paste0(outFilePrefix, ".tsv"))
    metric <- calcConsistencyMetric(evalData)

    description = paste0("convThreshold = ", convThreshold, "; intervalN = ", intervalN, "; binsize = ", binsize, "; no controls")
    plotConsistency(evalData, metric, description, paste0(outFilePrefix, ".pdf"))

    comparisonResults <- rbind(comparisonResults, c("No controls", convThreshold, intervalN, binsize, metric))

    print(paste("Executing parameter combination when using detection p-values", i, "...", sep=""))
  
    normData <- scanNorm(signalExprData, signalProbeSequences, signalPValueData = signalPValueData, convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, numCores=numCores, verbose=verbose)

    outFilePrefix <- paste0("Figures/", convThreshold, "_", intervalN, "_", binsize, "_DetectionPValues")

    evalData <- restructureConsistencyData(normData, paste0(outFilePrefix, ".tsv"))
    metric <- calcConsistencyMetric(evalData)
    
    description = paste0("convThreshold = ", convThreshold, "; intervalN = ", intervalN, "; binsize = ", binsize, "; detection p-values")
    plotConsistency(evalData, metric, description, paste0(outFilePrefix, ".pdf"))
    
    comparisonResults <- rbind(comparisonResults, c("Detection p-values", convThreshold, intervalN, binsize, metric))

    print(paste("Executing parameter combination with no controls or detection p-values", i, "...", sep=""))
  
    normData <- scanNorm(signalExprData, signalProbeSequences, convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, numCores=numCores, verbose=verbose)

    outFilePrefix <- paste0("Figures/", convThreshold, "_", intervalN, "_", binsize, "_ExpressionOnly")
    
    evalData <- restructureConsistencyData(normData, paste0(outFilePrefix, ".tsv"))
    metric <- calcConsistencyMetric(evalData)
    
    description = paste0("convThreshold = ", convThreshold, "; intervalN = ", intervalN, "; binsize = ", binsize, "; expression only")
    plotConsistency(evalData, metric, description, paste0(outFilePrefix, ".pdf"))
    
    comparisonResults <- rbind(comparisonResults, c("Expression only", convThreshold, intervalN, binsize, metric))
  }

  colnames(comparisonResults) <- c("Input data", "convThreshold", "intervalN", "binsize", "Metric")
  comparisonResults <- as_tibble(comparisonResults)
  write_tsv(arrange(comparisonResults, desc(Metric)), paramTuningOutFilePath)
}