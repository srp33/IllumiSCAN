# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("broom", "doParallel", "limma", "oligo", "tidyverse"))

library(broom)
library(doParallel)
library(limma)
library(oligo)
library(tidyverse)

#####################################################################
# Defining functions
#####################################################################

source("IllumiSCAN.R")

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

plotGC <- function(exprData, outFilePath) {
  meanEx <- apply(exprData, 1, mean)
  r = cor(meanEx, signalProbeSequenceGC, method = "spearman")
  r = paste0("r = ", round(r, 3))
  
  data <- tibble(GC = factor(signalProbeSequenceGC), Mean_Expression = meanEx)

  p <- ggplot(data, aes(x = GC, y = Mean_Expression)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color = "blue", alpha = 0.01, size = 0.5) +
    annotate("text", x = Inf, y = Inf, label = r, hjust = 1.5, vjust = 3, color = "red") +
    xlab("Proportion G/C nucleotides") +
    ylab("Expression level") +
    theme_bw(base_size = 14)
  
  print(p)
  
  ggsave(outFilePath, width = 8, height = 6)
}

parseConsistencyData <- function(exprData) {
  normSpikeInData <- exprData[spikeInProbeIDs,]
  normSpikeInData <- bind_cols(tibble(ProbeID=spikeInProbeIDs), normSpikeInData)
  normSpikeInData <- pivot_longer(normSpikeInData, -ProbeID, names_to = "SampleID", values_to = "Value")
  normSpikeInData <- inner_join(normSpikeInData, targetData, by="SampleID")
  
  evalData <- spikeInAnnotationData %>%
    dplyr::select(ProbeID, TargetID)
  evalData$ProbeID <- factor(as.character(evalData$ProbeID))
  evalData <- inner_join(evalData, normSpikeInData, by="ProbeID")

  return(evalData)
}

plotConsistency <- function(evalData, metric, description, outFilePath) {
  metric <- bquote(R^2 == .(round(metric, 3)))
  
  p <- mutate(evalData, SpikeConc = factor(SpikeConc, levels = spikeLevels)) %>%
    ggplot(aes(x = SpikeConc, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.08) +
    annotate("text", x = 10.5, y = 1, label = metric, size = 5, color = "red") +
    ggtitle(description) +
    theme_bw(base_size = 14)
  
  print(p)
  
  ggsave(outFilePath, width = 8, height = 6)
}

calcConsistencyMetric <- function(evalData)
{
  evalFit <- lm(Value~SpikeConc, data=evalData)
  
  return(summary(evalFit)$r.squared)
}

#####################################################################
# Retrieving spike-in data and configuring the experiment
#####################################################################

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

comparisonResults <- NULL

#####################################################################
# Save and plot the non-normalized data and detection p-values
#####################################################################

outFilePrefix <- "Figures/NonNormalizedExpression"
outFilePath <- paste0(outFilePrefix, ".tsv.gz")
if (!file.exists(outFilePath))
  write_tsv(signalExprData %>% as.data.frame() %>% rownames_to_column("ProbeID"), outFilePath)

if (!file.exists(paramTuningOutFilePath)) {
  plotGC(log2(signalExprData), paste0(outFilePrefix, "_GC.pdf")) # There is a large dynamic range, so log2-transform the values.

  evalData <- parseConsistencyData(signalExprData) %>%
    mutate(Value = log2(Value)) # There is a large dynamic range, so log2-transform the values.
  metric <- calcConsistencyMetric(evalData)

  plotConsistency(evalData, metric, "Non-normalized expression values", paste0(outFilePrefix, "_Eval.pdf"))
  comparisonResults <- rbind(comparisonResults, c("Non-normalized expression values", NA, NA, NA, metric))
}

outFilePrefix <- "Figures/DetectionPValues"
outFilePath <- paste0(outFilePrefix, ".tsv.gz")
if (!file.exists(outFilePath))
  write_tsv(signalPValueData %>% as.data.frame() %>% rownames_to_column("ProbeID"), outFilePath)

if (!file.exists(paramTuningOutFilePath)) {
  plotGC(log2(signalPValueData), paste0(outFilePrefix, "_GC.pdf"))
  
  evalData <- parseConsistencyData(signalPValueData)
  metric <- calcConsistencyMetric(evalData)
  
  plotConsistency(evalData, metric, "Detection p-values", paste0(outFilePrefix, "_Eval.pdf"))
  comparisonResults <- rbind(comparisonResults, c("Non-normalized expression values", NA, NA, NA, metric))
}

#####################################################################
# Parameter combination evaluation
#####################################################################

convThresholdOptions <- c(0.01, 0.1, 1, 10)
intervalNOptions <- c(1000, 10000, 20000, 50000)
binsizeOptions <- c(50, 500, 5000)

numCores = 4
verbose = FALSE

paramCombos <- expand.grid(convThresholdOptions, intervalNOptions, binsizeOptions)
colnames(paramCombos) <- c("convThreshold", "intervalN", "binsize")

paramTuningOutFilePath <- "Param_Tuning_Results.tsv"

if (!file.exists(paramTuningOutFilePath))
{
  for (i in 1:nrow(paramCombos))
  {
    print(paste("Executing parameter combination when using control data ", i, "...", sep=""))
    convThreshold <- paramCombos[i,1]
    intervalN <- paramCombos[i,2]
    binsize <- paramCombos[i,3]

    outFilePrefix <- paste0("Figures/", convThreshold, "_", intervalN, "_", binsize, "_WithControls")
    normalizedFilePath <- paste0(outFilePrefix, "_NormData.tsv.gz")

    if (file.exists(normalizedFilePath)) {
      normData <- read_tsv(normalizedFilePath)
      probeIDs <- pull(normData, ProbeID)
      normData <- select(normData, -ProbeID)
      normData <- as.matrix(normData)
      rownames(normData) <- probeIDs
    } else {
      normData <- scanNorm(signalExprData, signalProbeSequences, controlExprData = controlExprData, convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, numCores=numCores, verbose=verbose)
      write_tsv(normData %>% as.data.frame() %>% rownames_to_column("ProbeID"), normalizedFilePath)
    }
    
    plotGC(normData, paste0(outFilePrefix, "_GC.pdf"))

    evalData <- parseConsistencyData(normData)
    metric <- calcConsistencyMetric(evalData)

    description = paste0("convThreshold = ", convThreshold, "; intervalN = ", intervalN, "; binsize = ", binsize, "; with controls")
    plotConsistency(evalData, metric, description, paste0(outFilePrefix, "_Eval.pdf"))

    comparisonResults <- rbind(comparisonResults, c("With controls", convThreshold, intervalN, binsize, metric))

    print(paste("Executing parameter combination when using detection p-values", i, "...", sep=""))
    
    outFilePrefix <- paste0("Figures/", convThreshold, "_", intervalN, "_", binsize, "_DetectionPValues")
    normalizedFilePath <- paste0(outFilePrefix, "_NormData.tsv.gz")
    
    if (file.exists(normalizedFilePath)) {
      normData <- read_tsv(normalizedFilePath)
      probeIDs <- pull(normData, ProbeID)
      normData <- select(normData, -ProbeID)
      normData <- as.matrix(normData)
      rownames(normData) <- probeIDs
    } else {
      normData <- scanNorm(signalExprData, signalProbeSequences, signalPValueData = signalPValueData, convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, numCores=numCores, verbose=verbose)
      write_tsv(normData %>% as.data.frame() %>% rownames_to_column("ProbeID"), normalizedFilePath)
    }
    
    plotGC(normData, paste0(outFilePrefix, "_GC.pdf"))

    evalData <- parseConsistencyData(normData)
    metric <- calcConsistencyMetric(evalData)
    
    description = paste0("convThreshold = ", convThreshold, "; intervalN = ", intervalN, "; binsize = ", binsize, "; detection p-values")
    plotConsistency(evalData, metric, description, paste0(outFilePrefix, "_Eval.pdf"))
    
    comparisonResults <- rbind(comparisonResults, c("Detection p-values", convThreshold, intervalN, binsize, metric))

    print(paste("Executing parameter combination with no controls or detection p-values", i, "...", sep=""))
    
    outFilePrefix <- paste0("Figures/", convThreshold, "_", intervalN, "_", binsize, "_ExpressionOnly")
    normalizedFilePath <- paste0(outFilePrefix, "_NormData.tsv.gz")
    
    if (file.exists(normalizedFilePath)) {
      normData <- read_tsv(normalizedFilePath)
      probeIDs <- pull(normData, ProbeID)
      normData <- select(normData, -ProbeID)
      normData <- as.matrix(normData)
      rownames(normData) <- probeIDs
    } else {
      normData <- scanNorm(signalExprData, signalProbeSequences, convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, numCores=numCores, verbose=verbose)
      write_tsv(normData %>% as.data.frame() %>% rownames_to_column("ProbeID"), normalizedFilePath)
    }
    
    plotGC(normData, paste0(outFilePrefix, "_GC.pdf"))
    
    evalData <- parseConsistencyData(normData)
    metric <- calcConsistencyMetric(evalData)
    
    description = paste0("convThreshold = ", convThreshold, "; intervalN = ", intervalN, "; binsize = ", binsize, "; expression only")
    plotConsistency(evalData, metric, description, paste0(outFilePrefix, "_Eval.pdf"))
    
    comparisonResults <- rbind(comparisonResults, c("Expression only", convThreshold, intervalN, binsize, metric))
    #stop("got here")
  }

  colnames(comparisonResults) <- c("Input data", "convThreshold", "intervalN", "binsize", "Metric")
  comparisonResults <- as_tibble(comparisonResults)
  write_tsv(arrange(comparisonResults, desc(Metric)), paramTuningOutFilePath)
}

#TODO:
#  PlotGC bias before and after.
#  Plot the summarized results? Or just examine it as a table?
#  Create a scatterplot of the normalized data for the top-performing parameters vs. the bottom-performing params.