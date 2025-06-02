# install.packages("ppcor")
# install.packages("tidyverse")

library(ppcor)
library(tidyverse)

#####################################################################
# Defining functions
#####################################################################

source("IllumiSCAN.R")

# From Ringo/Bioconductor package
compute.gc <- function(probe.sequences, digits=2)
{
  stopifnot(is.character(probe.sequences))
  
  splitted.seqs <- strsplit(toupper(probe.sequences), split="")
  round(sapply(splitted.seqs, function(x) length(grep("[GC]", x))) / listLen(splitted.seqs), digits=digits)
}

averageAcrossReplicates <- function(data) {
  data %>%
    pivot_longer(-ProbeID, names_to = "SampleID", values_to = "Value") %>%
    inner_join(targetData, by = "SampleID") %>%
    group_by(ProbeID, SpikeConc) %>%
    dplyr::summarize(Value = mean(Value), .groups = "drop") %>%
    return()
}

plotGC <- function(exprData, signalProbeMeta, outFilePath, ylabSuffix="") {
  # We remove the spike-in probes so we can focus on GC bias in non-signal probes.
  data <- anti_join(exprData, spikeInAnnotationData, by = "ProbeID")

  mean_rho <- inner_join(data, signalProbeMeta, by = "ProbeID") %>%
    group_by(SpikeConc) %>%
    dplyr::summarize(rho = cor(GC_Proportion, Value, method = "spearman")) %>%
    pull(rho) %>%
    mean()
  
  set.seed(0)
  
  p <- inner_join(data, signalProbeMeta, by = "ProbeID") %>%
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

getDescriptiveStats <- function(normData) {
  x <- pull(normData, Value)
  return(c(min(x), mean(x), median(x), max(x), sd(x)))
}

# We calculate a consistency metric for each probe and then average across the probes.
# We only the highest spike-in concentration to evaluate normalization consistency,
# since these genes have the strongest signal and are less affected by detection noise.
calcConsistencyAcrossReplicates <- function(data) {
  spikeInData <- pivot_longer(data, -ProbeID, names_to = "SampleID", values_to = "Value") %>%
    inner_join(targetData, by = "SampleID") %>%
    inner_join(spikeInAnnotationData)

  # For each target and concentration, look across all replicates for both probes.
  # Calculate a dispersion metric and average across all targets and concentrations.
  # We expect that SCAN may perform better here because we are correcting for sequence composition.
  crossTargetConsistency <- group_by(spikeInData, TargetID, SpikeConc) %>%
    dplyr::summarize(SD = calcNormalizedStandardDeviation(Value)) %>%
    pull(SD) %>%
    mean()
  
  # Calculate a dispersion metric for cross-sample consistency for probes that have the same sequence.
  # We expect that SCAN may perform worse here because correcting for sequence composition is not relevant.
  crossSampleConsistency <- group_by(spikeInData, ProbeID, SpikeConc) %>%
    dplyr::summarize(SD = calcNormalizedStandardDeviation(Value)) %>%
    pull(SD) %>%
    mean()

  # Calculate a dispersion metric for within-sample consistency for probes with the same target.
  # We expect that SCAN may perform better here because we are correcting for sequence composition.
  withinSampleConsistency <- group_by(spikeInData, SampleID, TargetID, SpikeConc) %>%
    dplyr::summarize(SD = calcNormalizedStandardDeviation(Value)) %>%
    pull(SD) %>%
    mean()

  return(c(crossTargetConsistency, crossSampleConsistency, withinSampleConsistency))
}

calcNormalizedStandardDeviation <- function(x, numTrim = 0) {
  rmse <- sqrt(mean((x - mean(x))^2))

  y <- sort(x)
  range_x <- y[length(y) - numTrim] - y[1 + numTrim] # Find the range but ignore one on each extreme.

  if (range_x == 0) {
    return(0)
  }
  
  nSD <- rmse / range_x
  return(nSD)
}

calcDiffExpressionMetrics <- function(data, useSpikeIns, lowConc, highConc) {
  data <- pivot_longer(data, -ProbeID, names_to = "SampleID", values_to = "Value") %>%
    inner_join(targetData, by = "SampleID")
  
  if (useSpikeIns) {
    data <- inner_join(data, spikeInAnnotationData, by = "ProbeID")
  } else {
    data <- anti_join(data, spikeInAnnotationData, by = "ProbeID")
  }

  data <- dplyr::filter(data, SpikeConc %in% c(lowConc, highConc)) %>%
    pivot_wider(id_cols = ProbeID, names_from = SampleID, values_from = Value) %>%
    as.data.frame()

  rownames(data) <- data$ProbeID
  data$ProbeID <- NULL

  highSamples <- filter(targetData, SpikeConc == highConc) %>%
    pull(SampleID)
  
  cls <- as.integer(colnames(data) %in% highSamples)
  design <- model.matrix(~ cls)
  fit <- lmFit(data, design)
  fit2 <- eBayes(fit)
  results <- topTable(fit2, coef = 2, number = Inf)

  numSignificant <- dplyr::filter(results, adj.P.Val < 0.05) %>%
    nrow()
  propSignificant <- numSignificant / nrow(results)

  return(propSignificant)
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

sampleDataFilePath <- unzip(zipFilePath, files = "SampleProbeProfile.txt", exdir = tempdir())
controlDataFilePath <- unzip(zipFilePath, files = "ControlProbeProfile.txt", exdir = tempdir())

sampleProbeData <- read.ilmn(sampleDataFilePath)
controlProbeData <- read.ilmn(controlDataFilePath) # This also contains the spike-in data

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
  dplyr::filter(is.na(OtherHits)) %>%
  dplyr::filter(`SpliceJunction?` == "No") %>%
  dplyr::filter(`PerfectMatch?` == "Yes") %>%
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
  dplyr::select(ProbeID, Sequence, TargetID) %>%
  filter(!str_starts(TargetID, "neo")) # There is only one probe sequence for this target.
spikeInAnnotationData$TargetID <- sapply(spikeInAnnotationData$TargetID, function(x) {strsplit(x, "_")[[1]][1]})
spikeInAnnotationData$TargetID <- factor(spikeInAnnotationData$TargetID)
spikeInAnnotationData <- mutate(spikeInAnnotationData, GC_Proportion = compute.gc(Sequence, digits=2)) %>%
  mutate(GC_Bin = round(GC_Proportion * 2, 1) / 2) %>%
  mutate(GC_Bin = factor(GC_Bin))

allAnnotationData <- bind_rows(sampleAnnotationData, dplyr::select(controlAnnotationData, ProbeID, Sequence)) %>%
  distinct() %>%
  mutate(GC_Proportion = compute.gc(Sequence, digits=2)) %>%
  mutate(GC_Bin = round(GC_Proportion * 2, 1) / 2) %>%
  mutate(GC_Bin = factor(GC_Bin)) %>%
  bind_rows(dplyr::select(spikeInAnnotationData, ProbeID, Sequence, GC_Proportion, GC_Bin))

# Find which probes are in common between the annotations and data.
signalProbeIDs <- sort(intersect(allAnnotationData$ProbeID, rownames(sampleProbeData$E)))
controlProbeIDs <- sort(intersect(allAnnotationData$ProbeID, rownames(controlProbeData$E)))
spikeInProbeIDs <- sort(as.character(pull(spikeInAnnotationData, ProbeID)))

# Move the spike-in data from the control data to the signal data.
signalProbeIDs <- c(signalProbeIDs, spikeInProbeIDs)
controlProbeIDs <- setdiff(controlProbeIDs, spikeInProbeIDs) %>%
  setdiff(signalProbeIDs)

# Combine all the data.
allExprData <- rbind(sampleProbeData$E, controlProbeData$E)
allPValueData <- rbind(sampleProbeData$other$Detection, controlProbeData$other$Detection)

# Keep data only for the relevant probes.
signalExprData <- allExprData[signalProbeIDs,]
controlExprData <- allExprData[controlProbeIDs,]
signalPValueData <- allPValueData[signalProbeIDs,]
# controlPValueData <- allPValueData[controlProbeIDs,]

# Get the metadata together
signalExprMeta <- dplyr::filter(allAnnotationData, ProbeID %in% rownames(signalExprData))
signalExprMeta <- signalExprMeta[match(rownames(signalExprData), signalExprMeta$ProbeID),]
controlExprMeta <- dplyr::filter(allAnnotationData, ProbeID %in% rownames(controlExprData))
controlExprMeta <- controlExprMeta[match(rownames(controlExprData), controlExprMeta$ProbeID),]

targetData <- suppressMessages(read_delim(targetsFilePath, delim=" "))
targetData <- targetData %>%
  dplyr::select(ArrayNo, SpikeConc) %>%
  dplyr::rename(SampleID=ArrayNo)
targetData$SampleID <- sapply(targetData$SampleID, function(x) {paste(strsplit(x, "_")[[1]][1:2], collapse="_")})
spikeLevels <- c(0e+00, 1e-02, 3e-02, 1e-01, 3e-01, 1e+00, 3e+00, 1e+01, 3e+01, 1e+02, 3e+02, 1e+03)
#targetData$SpikeConc <- factor(targetData$SpikeConc, levels=spikeLevels)
targetData <- distinct(targetData)

#####################################################################
# Parameter combination evaluation
#####################################################################

testParamCombos1 <- function(paramCombos, paramTuningOutFilePath, numCores = 4) {
  outDirPath <- dirname(paramTuningOutFilePath)
  
  if (!file.exists(paramTuningOutFilePath))
  {
    dir.create(outDirPath, showWarnings = FALSE, recursive = TRUE)
    colnames(paramCombos) <- c("backgroundCorrectOption", "normalizationOption", "scanOption")
    comparisonResults <- NULL
    
    for (i in 1:nrow(paramCombos))
    {
      print(paste("Executing parameter combination ", i, "...", sep=""))
      paramCombo = as.vector(unlist(paramCombos[i,]))

      backgroundCorrectOption <- paramCombo[1]
      normalizationOption <- paramCombo[2]
      scanOption <- paramCombo[3]
      
      outFilePrefix <- paste0(outDirPath, "/", paste(paramCombo, collapse="____"))
      normalizedFilePath <- paste0(outFilePrefix, "_Data.tsv.gz")
      
      if (!file.exists(normalizedFilePath)) {
        normData <- normalizeBeadChipData(signalExprData,
                                          signalExprMeta$Sequence,
                                          controlExprData=controlExprData,
                                          controlProbeSequences=controlExprMeta$Sequence,
                                          detectionPValues=signalPValueData,
                                          correctBackgroundType=backgroundCorrectOption,
                                          vsnNormalize=(normalizationOption=="vsn"),
                                          quantileNormalize=(normalizationOption=="quantile"),
                                          scanNormalize=(scanOption=="SCAN"),
                                          numCores = numCores,
                                          verbose = TRUE)

        normData <- as.data.frame(normData) %>%
          rownames_to_column("ProbeID") %>%
          mutate(ProbeID = as.integer(str_replace(ProbeID, "^X", "")))

        write_tsv(normData, normalizedFilePath)
      }
      
      normData <- read_tsv(normalizedFilePath)
      
      normDataAveraged <- averageAcrossReplicates(normData)

      stats <- getDescriptiveStats(normDataAveraged)
      consistencyMetrics <- calcConsistencyAcrossReplicates(normData)
      notSpikedDiffExprMetric <- calcDiffExpressionMetrics(normData, useSpikeIns = FALSE, lowConc = 0, highConc = 1000)
      lowSpikedDiffExprMetric <- calcDiffExpressionMetrics(normData, useSpikeIns = TRUE, lowConc = 0.03, highConc = 0.1)
      gc_rho <- plotGC(normDataAveraged, signalExprMeta, paste0(outFilePrefix, "_GC.pdf"))
      conc_rho <- plotConcentrations(normDataAveraged, paste0(outFilePrefix, "_Spike.pdf"))

      comparisonResults <- rbind(comparisonResults, c(backgroundCorrectOption, normalizationOption, scanOption, stats, consistencyMetrics, notSpikedDiffExprMetric, lowSpikedDiffExprMetric, gc_rho, conc_rho))
    }
    

    colnames(comparisonResults) <- c("backgroundCorrection", "normalization", "scan", "min", "mean", "median", "max", "sd", "CrossTargetConsistency", "CrossSampleConsistency", "WithinSampleConsistency", "NoSpikeInDiffExprProportion", "LowSpikeInDiffExprProportion", "GC rho", "Spike-in concentration rho")

    comparisonResults <- as_tibble(comparisonResults) %>%
      mutate(min = as.numeric(min),
             max = as.numeric(max),
             sd = as.numeric(sd),
             CrossTargetConsistency = as.numeric(CrossTargetConsistency),
             CrossSampleConsistency = as.numeric(CrossSampleConsistency),
             WithinSampleConsistency = as.numeric(WithinSampleConsistency),
             NoSpikeInDiffExprProportion = as.numeric(NoSpikeInDiffExprProportion),
             LowSpikeInDiffExprProportion = as.numeric(LowSpikeInDiffExprProportion),
             `GC rho` = as.numeric(`GC rho`),
             `Spike-in concentration rho` = as.numeric(`Spike-in concentration rho`))
    
    write_tsv(comparisonResults, paramTuningOutFilePath)
  }
  
  return(read_tsv(paramTuningOutFilePath))
}

# This round of tuning looks for patterns across different types of input data and ways of normalizing.

backgroundCorrectOptions <- c("none", "controls", "detectionP")
normalizationOptions <- c("none", "quantile", "vsn")
scanOptions <- c("none", "SCAN")

paramCombos <- expand_grid(backgroundCorrectOptions, normalizationOptions, scanOptions)

paramTuning1FilePath <- "Param_Tuning_1/Param_Tuning_Summary.tsv"

comparisonResults <- testParamCombos1(paramCombos, paramTuning1FilePath)

comparisonResults <- read_tsv(paramTuning1FilePath) %>%
  mutate(GC_Rank = rank(abs(`GC rho`))) %>%
  mutate(Spike_Rank = rank(-`Spike-in concentration rho`)) %>%
  mutate(CrossTargetConsistency_Rank = rank(CrossTargetConsistency)) %>%
  mutate(CrossSampleConsistency_Rank = rank(CrossSampleConsistency)) %>%
  mutate(WithinSampleConsistency_Rank = rank(WithinSampleConsistency)) %>%
  mutate(NoSpikeInDiffExprProportion_Rank = rank(NoSpikeInDiffExprProportion)) %>%
  mutate(LowSpikeInDiffExprProportion_Rank = rank(-LowSpikeInDiffExprProportion)) %>%
  mutate(Combined_Rank = GC_Rank + Spike_Rank + CrossTargetConsistency_Rank + CrossSampleConsistency_Rank + WithinSampleConsistency_Rank + NoSpikeInDiffExprProportion_Rank + LowSpikeInDiffExprProportion_Rank) %>%
  arrange(Combined_Rank)

# This round of tuning looks for the best hyperparameter combination for SCAN.

testParamCombos2 <- function(paramCombos, paramTuningOutFilePath, numCores = 4) {
  outDirPath <- dirname(paramTuningOutFilePath)
  
  if (!file.exists(paramTuningOutFilePath))
  {
    dir.create(outDirPath, showWarnings = FALSE, recursive = TRUE)
    colnames(scanParamCombos) <- c("convThreshold", "intervalN", "binsize", "nbins", "controls")
    comparisonResults <- NULL
    
    for (i in 1:nrow(paramCombos))
    {
      print(paste("Executing parameter combination ", i, "...", sep=""))

      paramCombo = as.vector(unlist(paramCombos[i,]))
      convThreshold <- as.numeric(paramCombo[1])
      intervalN <- as.numeric(paramCombo[2])
      binsize <- as.numeric(paramCombo[3])
      nbins <- as.numeric(paramCombo[4])
      controls <- paramCombo[5]
      
      outFilePrefix <- paste0(outDirPath, "/", paste(paramCombo, collapse="____"))
      normalizedFilePath <- paste0(outFilePrefix, "_Data.tsv.gz")
      
      if (!file.exists(normalizedFilePath)) {
        if (controls == "controls") {
          controlExprDataTest <- controlExprData
          controlProbeSequencesTest <- controlExprMeta$Sequence
        } else {
          controlExprDataTest <- NULL
          controlProbeSequencesTest <- NULL
        }

        normData <- normalizeBeadChipData(signalExprData,
                                          signalExprMeta$Sequence,
                                          controlExprData = controlExprDataTest,
                                          controlProbeSequences = controlProbeSequencesTest,
                                          detectionPValues = signalPValueData,
                                          vsnNormalize = TRUE,
                                          scanNormalize = TRUE,
                                          scanConvThreshold = convThreshold,
                                          scanIntervalN = intervalN,
                                          scanBinsize = binsize,
                                          scanNbins = nbins,
                                          numCores = numCores,
                                          verbose = TRUE)

        normData <- as.data.frame(normData) %>%
          rownames_to_column("ProbeID") %>%
          mutate(ProbeID = as.integer(str_replace(ProbeID, "^X", ""))) %>%
          write_tsv(normalizedFilePath)
      }

      normData <- read_tsv(normalizedFilePath)

      normDataAveraged <- averageAcrossReplicates(normData)

      stats <- getDescriptiveStats(normDataAveraged)
      consistencyMetrics <- calcConsistencyAcrossReplicates(normData)
      notSpikedDiffExprMetric <- calcDiffExpressionMetrics(normData, useSpikeIns = FALSE, lowConc = 0, highConc = 1000)
      lowSpikedDiffExprMetric <- calcDiffExpressionMetrics(normData, useSpikeIns = TRUE, lowConc = 0.03, highConc = 0.1)
      gc_rho <- plotGC(normDataAveraged, signalExprMeta, paste0(outFilePrefix, "_GC.pdf"))
      conc_rho <- plotConcentrations(normDataAveraged, paste0(outFilePrefix, "_Spike.pdf"))

      comparisonResults <- rbind(comparisonResults, c(convThreshold, intervalN, binsize, nbins, controls, stats, consistencyMetrics, notSpikedDiffExprMetric, lowSpikedDiffExprMetric, gc_rho, conc_rho))
    }

    colnames(comparisonResults) <- c("convThreshold", "intervalN", "binsize", "nbins", "controls", "min", "mean", "median", "max", "sd", "CrossTargetConsistency", "CrossSampleConsistency", "WithinSampleConsistency", "NoSpikeInDiffExprProportion", "LowSpikeInDiffExprProportion", "GC rho", "Spike-in concentration rho")

    comparisonResults <- as_tibble(comparisonResults) %>%
      mutate(min = as.numeric(min),
             max = as.numeric(max),
             sd = as.numeric(sd),
             CrossTargetConsistency = as.numeric(CrossTargetConsistency),
             CrossSampleConsistency = as.numeric(CrossSampleConsistency),
             WithinSampleConsistency = as.numeric(WithinSampleConsistency),
             NoSpikeInDiffExprProportion = as.numeric(NoSpikeInDiffExprProportion),
             LowSpikeInDiffExprProportion = as.numeric(LowSpikeInDiffExprProportion),
             `GC rho` = as.numeric(`GC rho`),
             `Spike-in concentration rho` = as.numeric(`Spike-in concentration rho`))

    write_tsv(comparisonResults, paramTuningOutFilePath)
  }

  return(read_tsv(paramTuningOutFilePath))
}

convThresholdOptions <- c(0.01, 0.1, 0.5, 1)
intervalNOptions <- c(1000, 5000, 10000, 50000)
binsizeOptions <- c(50, 500, 5000)
nbinsOptions <- c(10, 25, 50)
controlsOptions <- c("controls", "nocontrols")

paramTuning2FilePath <- "Param_Tuning_2/Param_Tuning_Summary.tsv"
scanParamCombos <- expand_grid(convThresholdOptions, intervalNOptions, binsizeOptions, nbinsOptions, controlsOptions)
comparisonResults <- testParamCombos2(scanParamCombos, paramTuning2FilePath)

comparisonResults <- read_tsv(paramTuning2FilePath) %>%
  mutate(GC_Rank = rank(abs(`GC rho`))) %>%
  mutate(Spike_Rank = rank(-`Spike-in concentration rho`)) %>%
  mutate(CrossTargetConsistency_Rank = rank(CrossTargetConsistency)) %>%
  mutate(CrossSampleConsistency_Rank = rank(CrossSampleConsistency)) %>%
  mutate(WithinSampleConsistency_Rank = rank(WithinSampleConsistency)) %>%
  mutate(NoSpikeInDiffExprProportion_Rank = rank(-NoSpikeInDiffExprProportion)) %>%
  mutate(LowSpikeInDiffExprProportion_Rank = rank(-LowSpikeInDiffExprProportion)) %>%
  mutate(Combined_Rank = GC_Rank + Spike_Rank + CrossTargetConsistency_Rank + CrossSampleConsistency_Rank + WithinSampleConsistency_Rank + NoSpikeInDiffExprProportion_Rank + LowSpikeInDiffExprProportion_Rank) %>%
  arrange(Combined_Rank)

#####################################################################
# Plot the summarized results of parameter evaluation
#####################################################################

ggplot(comparisonResults, aes(x = as.factor(convThreshold), y = Combined_Rank)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
  xlab("SCAN convergence threshold parameter") +
  ylab("Combined rank (lower is better)") +
  theme_bw(base_size = 14)

ggplot(comparisonResults, aes(x = as.factor(intervalN), y = Combined_Rank)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
  xlab("SCAN intervalN parameter") +
  ylab("Combined rank (lower is better)") +
  theme_bw(base_size = 14)

ggplot(comparisonResults, aes(x = as.factor(binsize), y = Combined_Rank)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
  xlab("SCAN binsize parameter") +
  ylab("Combined rank (lower is better)") +
  theme_bw(base_size = 14)

ggplot(comparisonResults, aes(x = as.factor(nbins), y = Combined_Rank)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
  xlab("SCAN nbins parameter") +
  ylab("Combined rank (lower is better)") +
  theme_bw(base_size = 14)

ggplot(comparisonResults, aes(x = as.factor(controls), y = Combined_Rank)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "red", alpha = 0.5, size = 1.5) +
  xlab("SCAN controls parameter") +
  ylab("Combined rank (lower is better)") +
  theme_bw(base_size = 14)