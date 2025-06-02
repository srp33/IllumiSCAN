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
  # We remove the spike-in probes so we can focus just on GC bias in non-signal probes.
  data <- anti_join(exprData, spikeInAnnotationData, by = "ProbeID")

  mean_rho <-
    inner_join(data, signalProbeMeta, by = "ProbeID") %>%
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

getStats <- function(normData) {
  x <- pull(normData, Value)
  return(c(min(x), mean(x), median(x), max(x), sd(x)))
}

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
                                           exprColumnPattern="GSM",
                                           adjustBackground=FALSE,
                                           useSCAN=FALSE,
                                           quantileNormalize=FALSE,
                                           log2Transform=FALSE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="GSM",
                                           adjustBackground=FALSE,
                                           useSCAN=FALSE,
                                           quantileNormalize=FALSE,
                                           log2Transform=TRUE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="GSM",
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           quantileNormalize=FALSE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="GSM",
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           quantileNormalize=TRUE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="GSM",
                                           adjustBackground=FALSE,
                                           useSCAN=TRUE, scanConvThreshold=0.5,
                                           quantileNormalize=FALSE,
                                           numCores=4,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="GSM",
                                           adjustBackground=FALSE,
                                           useSCAN=TRUE, scanConvThreshold=0.5,
                                           quantileNormalize=TRUE,
                                           numCores=4,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="GSM",
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

############################################################
gseID <- "GSE43692"
############################################################

normalized1 <- normalizeBeadChipDataFromGEO(gseID,
                                            exprColumnPattern="Value ",
                                            adjustBackground=TRUE,
                                            useSCAN=FALSE,
                                            verbose=TRUE)

normalized2 <- normalizeBeadChipDataFromGEO(gseID,
                                            exprColumnPattern="Value ",
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
                                           exprColumnPattern="SAMPLE ",
                                           adjustBackground=TRUE,
                                           useSCAN=FALSE,
                                           verbose=TRUE)

normalized <- normalizeBeadChipDataFromGEO(gseID,
                                           exprColumnPattern="SAMPLE ",
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

# TODO: Use the annotations to look for control probes when we do SCAN normalization. Is this status in the metadata?
#      Then discard these probes after normalization?
# TODO: Use arrayQualityMetrics to assess quality and add quality findings to the phenoData.
#       https://github.com/ebi-gene-expression-group/microarray-import/blob/develop/bin/arrayQC.R
# TODO: Work with Down Syndrome data from Process_Illumina_Microarray.R.
# TODO: Try VSN with "Normalization against an existing reference dataset"? (see docs for vsn package)
#         Or at least make a note of it.
# TODO: For the paper: "Given these observations and previous work for Affymetrix arrays, it would seem that more sophisticated methods than background normalisation are needed to account for sequence-specific hybridisation effects." (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-85)