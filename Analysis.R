# install.packages("ppcor")
# install.packages("tidyverse")

library(ppcor)
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
  
  splitted.seqs <- strsplit(toupper(probe.sequences), split="")
  round(sapply(splitted.seqs, function(x) length(grep("[GC]", x))) / listLen(splitted.seqs), digits=digits)
}

averageAcrossReplicates <- function(data, log2Transform = FALSE) {
  data %>%
    pivot_longer(-ProbeID, names_to = "SampleID", values_to = "Value") %>%
    mutate(Value = if (log2Transform) log2(Value) else Value) %>%
    inner_join(targetData, by = "SampleID") %>%
    group_by(ProbeID, SpikeConc) %>%
    dplyr::summarize(Value = mean(Value), .groups = "drop") %>%
    return()
}

plotGC <- function(exprData, signalProbeMeta, outFilePath, ylabSuffix="") {
  mean_rho <- inner_join(exprData, signalProbeMeta, by = "ProbeID") %>%
    group_by(SpikeConc) %>%
    dplyr::summarize(rho = cor(GC_Proportion, Value, method = "spearman")) %>%
    pull(rho) %>%
    mean()
  
  set.seed(0)
  
  p <- inner_join(exprData, signalProbeMeta, by = "ProbeID") %>%
    ggplot(aes(x = GC_Bin, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = SpikeConc), alpha = 0.01, size = 0.5) +
    annotate("text", x = Inf, y = Inf, label = paste0("Mean rho = ", round(mean_rho, 3)), hjust = 1.5, vjust = 3, color = "red", size = 4) +
    xlab("Proportion G/C nucleotides") +
    ylab(paste0("Expression value", ylabSuffix)) +
    theme_bw(base_size = 14)
  
  print(p) %>%
    suppressWarnings()
  
  ggsave(outFilePath, width = 8, height = 6)
  
  return(mean_rho)
}

plotConcentrations <- function(normData, outFilePath, ylabSuffix="") {
  evalData <- inner_join(normData, spikeInAnnotationData)
  
  mean_rho <- dplyr::group_by(evalData, ProbeID) %>%
    dplyr::summarize(rho = cor(SpikeConc, Value, method = "spearman")) %>%
    pull(rho) %>%
    mean()
  
  set.seed(0)
  
  p <- mutate(evalData, SpikeConc = factor(SpikeConc, levels = spikeLevels)) %>%
    ggplot(aes(x = SpikeConc, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.08) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.35, vjust = 2, label = paste0("Mean rho = ", round(mean_rho, 3)), color = "red", size = 4) +
    xlab("Spike-in concentration") +
    ylab(paste0("Expression value", ylabSuffix)) +
    theme_bw(base_size = 14)
  
  print(p) %>%
    suppressWarnings()
  
  ggsave(outFilePath, width = 8, height = 6)
  
  return(mean_rho)
}

getStats <- function(normData) {
  x <- pull(normData, Value)
  return(c(min(x), mean(x), median(x), max(x), sd(x)))
}

#####################################################################
# Retrieve spike-in data and configure the experiment
#####################################################################

tmpDir <- "tmp"
zipUrl <- "https://github.com/markdunning/statistical-issues-illumina-microarray/raw/master/spike_beadstudio_output.zip"
zipFilePath <- paste(tmpDir, "/spike_beadstudio_output.zip", sep="")
if (!file.exists(zipFilePath)) {
  download.file(zipUrl, zipFilePath)
}

sampleDataFilePath <- unz(zipFilePath, "SampleProbeProfile.txt")
controlDataFilePath <- unz(zipFilePath, "ControlProbeProfile.txt")

#TODO: Change this to use read.ilmn?
sampleProbeData <- readProbeData(sampleDataFilePath, "AVG_Signal", "Detection Pval")
controlProbeData <- readProbeData(controlDataFilePath, "AVG_Signal", "Detection Pval") # This also contains the spike-in data

sampleAnnotationURL <- "https://raw.githubusercontent.com/markdunning/statistical-issues-illumina-microarray/master/NewAnnotationM1.txt"
controlAnnotationURL <- "https://raw.githubusercontent.com/markdunning/statistical-issues-illumina-microarray/master/ControlSequences.txt"
spikeInProfileURL <- "https://raw.githubusercontent.com/markdunning/statistical-issues-illumina-microarray/master/spikeins_profile.csv"
targetsURL <- "https://github.com/markdunning/statistical-issues-illumina-microarray/raw/master/spike_targets.txt"

sampleAnnotationFilePath <- "tmp/NewAnnotationM1.txt"
controlAnnotationFilePath <- "tmp/ControlSequences.txt"
spikeInProfileFilePath <- "tmp/spikeins_profile.csv"
targetsFilePath <- "tmp/spike_targets.txt"

if (!file.exists(targetsFilePath)) {
  download.file(sampleAnnotationURL, sampleAnnotationFilePath)
  download.file(controlAnnotationURL, controlAnnotationFilePath)
  download.file(spikeInProfileURL, spikeInProfileFilePath)
  download.file(targetsURL, targetsFilePath)
}

sampleAnnotationData <- suppressMessages(read_tsv(sampleAnnotationFilePath)) %>%
  dplyr::select(ProbeId0, Sequence) %>%
  dplyr::rename(ProbeID=ProbeId0)

controlAnnotationData <- suppressWarnings(suppressMessages(read_tsv(controlAnnotationFilePath))) %>%
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

signalExprMeta <- dplyr::select(signalExprData, ProbeID, Sequence) %>%
  # mutate(ProbeID = as.character(ProbeID)) %>%
  mutate(GC_Proportion = compute.gc(Sequence, digits=2)) %>%
  mutate(GC_Bin = round(GC_Proportion * 2, 1) / 2) %>%
  mutate(GC_Bin = factor(GC_Bin))

signalExprMatrix <- as.matrix(dplyr::select(signalExprData, -ProbeID, -Sequence))
signalPValueMatrix <- as.matrix(dplyr::select(signalPValueData, -ProbeID, -Sequence))
rownames(signalExprMatrix) <- pull(signalExprData, ProbeID)
rownames(signalPValueMatrix) <- pull(signalExprData, ProbeID)

controlProbeIDs <- pull(controlExprData, ProbeID)
controlProbeSequences <- pull(controlExprData, Sequence)
controlExprMatrix <- as.matrix(dplyr::select(controlExprData, -ProbeID, -Sequence))
rownames(controlExprMatrix) <- controlProbeIDs

targetData <- suppressMessages(read_delim(targetsFilePath, delim=" "))
targetData <- targetData %>% dplyr::select(ArrayNo, SpikeConc) %>% dplyr::rename(SampleID=ArrayNo)
targetData$SampleID <- sapply(targetData$SampleID, function(x) {paste(strsplit(x, "_")[[1]][1:2], collapse="_")})
spikeLevels <- c(0e+00, 1e-02, 3e-02, 1e-01, 3e-01, 1e+00, 3e+00, 1e+01, 3e+01, 1e+02, 3e+02, 1e+03)
#targetData$SpikeConc <- factor(targetData$SpikeConc, levels=spikeLevels)
targetData <- distinct(targetData)

#####################################################################
# Parameter combination evaluation
#####################################################################

dir.create("Outputs", showWarnings = FALSE, recursive = TRUE)

comparisonResults <- NULL
paramTuningOutFilePath <- "Param_Tuning_Results.tsv"

# backgroundCorrectOptions <- c("none")
# normalizationOptions <- c("none")
# scanOptions <- c("none")

backgroundCorrectOptions <- c("none", "controls", "detectionP")
normalizationOptions <- c("none", "quantile", "vsn")
scanOptions <- c("none", "SCAN", "UPC", "SCANnocontrols")

# convThresholdOptions <- c(0.01, 0.5, 1, 5)
# intervalNOptions <- c(1000, 10000, 50000)
# binsizeOptions <- c(50, 500, 5000)

numCores = 4
verbose = FALSE

paramCombos <- expand_grid(backgroundCorrectOptions, normalizationOptions, scanOptions)
colnames(paramCombos) <- c("backgroundCorrectOption", "normalizationOption", "scanOption")

if (!file.exists(paramTuningOutFilePath))
{
  for (i in 1:nrow(paramCombos))
  {
    print(paste("Executing parameter combination when using control data ", i, "...", sep=""))
    paramCombo = as.vector(unlist(paramCombos[i,]))
    backgroundCorrectOption <- paramCombo[1]
    normalizationOption <- paramCombo[2]
    scanOption <- paramCombo[3]
    
    outFilePrefix <- paste0("Outputs/", paste(paramCombo, collapse="_"))
    normalizedFilePath <- paste0(outFilePrefix, "_Data.tsv.gz")
    
    probeSequences <- signalExprData$Sequence
    
    if (!file.exists(normalizedFilePath)) {
      if (backgroundCorrectOption == "none") {
        normData <- signalExprMatrix
      } else {
        if (backgroundCorrectOption == "controls") {
          normData <- backgroundCorrect(signalExprMatrix, controlExprData=controlExprMatrix)
          probeSequences <- c(probeSequences, controlProbeSequences)
        } else {
          normData <- backgroundCorrect(signalExprMatrix, signalPValueData=signalPValueMatrix)
        }
      }
      
      # This is supposed to be done with raw intensities (I assume background correction is okay).    
      if (normalizationOption == "vsn") {
        normData <- normalizeVSN(normData)
        normData <- 2^normData # Reverse the log2 transformation
      }
      
      if (scanOption %in% c("SCAN", "UPC", "SCANnocontrols")) {
        # It is possible backgroundCorrection uses controls but SCAN does not.
        # Also, when backgroundCorrection = "none", we should get approximately
        #   the same result whether scanOption is "SCAN" or "SCANnocontrols".
        if (scanOption == "SCANnocontrols") {
          normData2 <- normData[1:nrow(signalExprMatrix),]
          probeSequences2 <- signalExprData$Sequence
        } else {
          normData2 <- normData
          probeSequences <- probeSequences
        }
        
        normData <- scanNorm(normData2, probeSequences2, asUPC = scanOption == "UPC", verbose = TRUE)
        #convThreshold = convThreshold, intervalN = intervalN, binsize = binsize, normalizationType=normalizationType, quantileNormalize=quantileNormalize, numCores=numCores, verbose=verbose
      }
      
      if (normalizationOption == "quantile") {
        normData <- normalizeBetweenArrays(normData, method = "quantile")
      }
      
      # Remove control probes if they are there.
      normData <- normData[1:nrow(signalExprMatrix),]
      
      if (scanOption != "UPC") {
        normData <- log2(normData)
      }
      
      normData <- as.data.frame(normData) %>%
        rownames_to_column("ProbeID") %>%
        mutate(ProbeID = as.integer(str_replace(ProbeID, "^X", ""))) %>%
        averageAcrossReplicates()
      
      write_tsv(normData, normalizedFilePath)
    }
    
    normData <- read_tsv(normalizedFilePath)
    
    gc_rho <- plotGC(normData, signalExprMeta, paste0(outFilePrefix, "_GC.pdf"))
    conc_rho <- plotConcentrations(normData, paste0(outFilePrefix, "_Spike.pdf"))
    stats <- getStats(normData)
    
    comparisonResults <- rbind(comparisonResults, c(backgroundCorrectOption, normalizationOption, scanOption, stats, gc_rho, conc_rho))
  }
  
  colnames(comparisonResults) <- c("backgroundCorrection", "normalization", "scan", "min", "mean", "median", "max", "sd", "GC rho", "Spike-in concentration rho")
  comparisonResults <- as_tibble(comparisonResults) %>%
    mutate(min = as.numeric(min),
           max = as.numeric(max),
           sd = as.numeric(sd),
           `GC rho` = as.numeric(`GC rho`),
           `Spike-in concentration rho` = as.numeric(`Spike-in concentration rho`))
  
  write_tsv(comparisonResults, paramTuningOutFilePath)
}

#####################################################################
# Evaluate "real-world" datasets
#####################################################################

calculateReplicateScore <- function(eSet, groupList) {
  probeGCProportion <- compute.gc(fData(eSet)$Probe_Sequence)
  
  # Compute mean correlation between two sample sets
  getMeanCorrelation <- function(samples1, samples2) {
    corVals <- c()
    for (s1 in samples1) {
      for (s2 in samples2) {
        if (s1 != s2) {
          corMethod = "spearman"
          
          print(s1)
          corValue <- pcor.test(exprs(eSet)[,s1], exprs(eSet)[,s2], probeGCProportion, method=corMethod)$estimate
          corVals <- c(corVals, corValue)
        }
      }
    }
    
    return(mean(corVals))
  }
  
  result <- NULL
  groupNames <- names(groupList)
  
  for (i in seq_along(groupNames)) {
    for (j in seq_along(groupNames)) {
      groupsSorted <- sort(c(groupNames[i], groupNames[j]))
      
      group1 <- groupsSorted[1]
      group2 <- groupsSorted[2]
      
      if (!is.null(result)) {
        alreadyDone <- dplyr::filter(result, Group1 == group1 & Group2 == group2) %>%
          nrow() %>%
          `>`(0)
        
        if (alreadyDone) {
          next
        }
      }
      
      samples1 <- groupList[[group1]]
      samples2 <- groupList[[group2]]
      
      meanCorr <- getMeanCorrelation(samples1, samples2)
      comparisonType <- if (group1 == group2) "within-group" else "between-group"
      
      row <- tibble(Group1 = group1, Group2 = group2, Mean_Correlation = meanCorr, Comparison_Type = comparisonType)
      if (is.null(result)) {
        result <- row
      } else {
        result <- bind_rows(result, row)
      }
    }
  }
  
  return(result)
}

############################################################
gseID <- "GSE31909"
############################################################

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=FALSE,
                                           useSCAN=FALSE,
                                           quantileNormalize=FALSE,
                                           log2Transform=FALSE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=FALSE,
                                           useSCAN=FALSE,
                                           quantileNormalize=FALSE,
                                           log2Transform=TRUE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           quantileNormalize=FALSE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           quantileNormalize=TRUE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=FALSE,
                                           useSCAN=TRUE, scanConvThreshold=0.5,
                                           quantileNormalize=FALSE,
                                           numCores=4,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=FALSE,
                                           useSCAN=TRUE, scanConvThreshold=0.5,
                                           quantileNormalize=TRUE,
                                           numCores=4,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=TRUE,
                                           useSCAN=TRUE, scanConvThreshold=0.5,
                                           quantileNormalize=TRUE,
                                           numCores=4,
                                           verbose=TRUE)

groups <- list(
  HEMn = c("GSM790975", "GSM790976", "GSM790977"),
  HEMa = c("GSM790978", "GSM790979", "GSM790980"),
  SKMEL28 = c("GSM790981", "GSM790982", "GSM790983"),
  LOXIMVI = c("GSM790984", "GSM790985", "GSM790986")
)

# To assess within-group consistency while accounting for overall gene expression variability, we calculated the Within-to-Total Variability Ratio (WTVR) for each gene. For each normalization method, we first computed the standard deviation of expression values for each gene within each biological group (i.e., set of replicates). We then averaged the within-group standard deviations across the four groups to obtain a single within-group variability estimate per gene. To contextualize this variability, we divided the averaged within-group standard deviation by the standard deviation of the same gene across all samples. This ratio reflects the proportion of a gene’s total variability that is attributable to variation within biological replicates. Lower WTVR values indicate greater consistency among replicates relative to the gene’s overall variation. We used the distribution of WTVR values across all genes to compare the effectiveness of each normalization method in reducing technical noise while preserving biological structure.
calculateRelativeVariability <- function(exprMatrix, groups) {
  allSdRatios <- c()

  for (group in names(groups)) {
    groupSamples <- groups[[group]]

    overallSDs <- apply(exprMatrix, 1, sd)
    groupSDs <- apply(exprMatrix[,groupSamples], 1, sd)
    
    sdRatios <- groupSDs / overallSDs
    sdRatios <- sdRatios[!is.nan(sdRatios)]
    
    allSdRatios <- c(allSdRatios, mean(sdRatios))
  }

  tibble(Group = names(groups), WTVR = allSdRatios) %>%
    return()
}

# Calculate within-group consistency
metrics <- calculateRelativeVariability(exprs(normalized), groups)

# Calculate within- and between group consistency
#metrics <- calculateReplicateScore(normalized, groups)

#TODO: Calculate GC correlation across all genes

############################################################
gseID <- "GSE43692"
############################################################

normalized1 <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           verbose=TRUE)

normalized2 <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=FALSE,
                                           useSCAN=TRUE, scanConvThreshold=0.5,
                                           numCores=4,
                                           verbose=TRUE)

# This is a good one to test having multiple non-normalized files with less-conventional names.
#   However, it probably doesn't make sense for a benchmark comparison because the values are not fully unnormalized.

############################################################
gseID <- "GSE169568"
############################################################

# This is a large, biomarker dataset (n = 204).

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=FALSE,
                                           useSCAN=TRUE, scanConvThreshold = 0.5,
                                           numCores=4,
                                           verbose=TRUE)

############################################################
gseID <- "GSE201405"
############################################################

# This one has IDAT files (n = 36) and no non-normalized files.
normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           verbose=TRUE)

############################################################
gseID <- "GSE224309"
############################################################

# This one has non-normalized data in Excel (only) (n = 16).

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           nonNormalizedFilePattern="GSE224309_Raw_Data.xlsx",
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           verbose=TRUE)

#####################################################################
# Plot the summarized results of parameter evaluation
#####################################################################

comparisonResults <- read_tsv(paramTuningOutFilePath) %>%
  mutate(GC_Rank = rank(abs(`GC rho`))) %>%
  mutate(Spike_Rank = rank(-`Spike-in concentration rho`)) %>%
  mutate(Combined_Rank = GC_Rank + Spike_Rank) %>%
  arrange(Combined_Rank)

stop("got to here...................")

# filter(comparisonResults, !is.na(convThreshold)) %>%
#   ggplot(aes(x = as.factor(convThreshold), y = `GC rho`)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#     xlab("SCAN convergence threshold parameter") +
#     ylab("Mean rho statistic") +
#     ggtitle("Correlation between GC content and expression values per spike-in concentration") +
#     theme_bw(base_size = 14)
# 
# filter(comparisonResults, !is.na(intervalN)) %>%
#   ggplot(aes(x = as.factor(intervalN), y = `GC rho`)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#   xlab("SCAN intervalN parameter") +
#   ylab("Mean rho statistic") +
#   ggtitle("Correlation between GC content and expression values per spike-in concentration") +
#   theme_bw(base_size = 14)
# 
# filter(comparisonResults, !is.na(binsize)) %>%
#   ggplot(aes(x = as.factor(binsize), y = `GC rho`)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#   xlab("SCAN binsize parameter") +
#   ylab("Mean rho statistic") +
#   ggtitle("Correlation between GC content and expression values per spike-in concentration") +
#   theme_bw(base_size = 14)
# 
# filter(comparisonResults, !is.na(convThreshold)) %>%
#   ggplot(aes(x = as.factor(convThreshold), y = `Spike-in concentration rho`)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#     xlab("SCAN convergence threshold parameter") +
#     ylab("Mean rho statistic") +
#     ggtitle("Correlation between spike-in concentration and expression values per spike-in probe") +
#     theme_bw(base_size = 14)
# 
# filter(comparisonResults, !is.na(intervalN)) %>%
#   ggplot(aes(x = as.factor(intervalN), y = `Spike-in concentration rho`)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#   xlab("SCAN intervalN parameter") +
#   ylab("Mean rho statistic") +
#   ggtitle("Correlation between spike-in concentration and expression values per spike-in probe") +
#   theme_bw(base_size = 14)
# 
# filter(comparisonResults, !is.na(binsize)) %>%
#   ggplot(aes(x = as.factor(binsize), y = `Spike-in concentration rho`)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#   xlab("SCAN binsize parameter") +
#   ylab("Mean rho statistic") +
#   ggtitle("Correlation between spike-in concentration and expression values per spike-in probe") +
#   theme_bw(base_size = 14)
# 
# baselineResults = filter(comparisonResults, `Input data` == "Non-normalized expression values")
# benchmarkResults = filter(comparisonResults, `Input data` != "Non-normalized expression values") %>%
#   filter(`Input data` != "Detection p-values")
# 
# ggplot(benchmarkResults, aes(x = `Input data`, y = Combined_Rank)) +
#   geom_boxplot() +
#   geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
#   geom_hline(yintercept = pull(baselineResults, Combined_Rank), color = "blue", linetype = "dashed") +
#   xlab("Input data type") +
#   ylab("Combined rank (lower is better)") +
#   theme_bw()

# TODO: Assess consistency of replicates: GSE31909, GSE43692, GSE169568
#         https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=500&platform=10558&zsort=date
#       Add quantile normalization?
# TODO: Support IDAT files in GEO.
# TODO: Use arrayQualityMetrics to assess quality and add quality findings to the phenoData.
# TODO: For now, we just picked one of the parameter combinations below. Need to compare them?
# TODO: Rework the above graphs.
# TODO: Try VSN with "Normalization against an existing reference dataset"? (see docs for vsn package)
#         Or at least make a note of it.
# TODO: Work with Down Syndrome data from Process_Illumina_Microarray.R.

#####################################################################
# Plot different versions of the data against each other.
#####################################################################

# nonNormalized = read_tsv("Outputs/NonNormalizedExpression.tsv.gz") %>%
#   averageAcrossReplicates(log2Transform = TRUE) %>%
#   dplyr::rename(`Not normalized` = Value)
# 
# normWithControls = read_tsv("Outputs/0.01_1000_50_WithControls_NormData.tsv.gz") %>%
#   averageAcrossReplicates() %>%
#   dplyr::rename(`Normalized with controls` = Value)
# 
# normDetectionP = read_tsv("Outputs/0.01_1000_50_DetectionPValues_NormData.tsv.gz") %>%
#   averageAcrossReplicates() %>%
#   dplyr::rename(`Normalized with detection p-values` = Value)
# 
# normExpressionOnly = read_tsv("Outputs/0.01_1000_50_ExpressionOnly_NormData.tsv.gz") %>%
#   averageAcrossReplicates() %>%
#   dplyr::rename(`Normalized with expression data only` = Value)
# 
# plot_data = inner_join(nonNormalized, normWithControls) %>%
#   inner_join(normDetectionP) %>%
#   inner_join(normExpressionOnly)
# 
# cor(pull(plot_data, `Not normalized`), pull(plot_data, `Normalized with controls`), method="spearman")
# cor(pull(plot_data, `Normalized with controls`), pull(plot_data, `Normalized with detection p-values`), method="spearman")
# cor(pull(plot_data, `Normalized with controls`), pull(plot_data, `Normalized with expression data only`), method="spearman")
# cor(pull(plot_data, `Normalized with detection p-values`), pull(plot_data, `Normalized with expression data only`), method="spearman")
# 
# plot_data %>%
#   # ggplot(aes(x = `Not normalized`, y = `Normalized with controls`)) +
#   # ggplot(aes(x = `Normalized with controls`, y = `Normalized with detection p-values`)) +
#   ggplot(aes(x = `Normalized with controls`, y = `Normalized with expression data only`)) +
#     # geom_point(alpha = 0.05, size = 0.5) +
#     geom_hex(bins = 100) +
#     scale_fill_viridis_c() +
#    # geom_smooth(method  = "lm", color = "red", size = 0.5) +
#     facet_wrap(vars(SpikeConc)) +
#     theme_bw()