# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# for (packageName in c("doParallel", "GEOquery", "httr", "limma", "oligo", "illuminaHumanv4.db", "readxl")) {
#   BiocManager::install(packageName)
# }

library(doParallel)
library(GEOquery)
library(httr)
library(limma)
library(oligo)
library(readxl)

# normalizeBeadChipDataFromGEO <- function(gseID, nonNormalizedDataFilePaths=NULL, nonNormalizedFileSuffix="non_normalized", nonNormalizedFileFieldDelimiter="\t", probeIDColumn=NULL, exprColumnPattern=NULL, detectionPValueColumnPattern=NULL, adjustBackground=TRUE, useSCAN=TRUE, scanConvThreshold=0.5, numCores=1, verbose=FALSE) {
normalizeBeadChipDataFromGEO <- function(gseID, nonNormalizedFilePattern="non.*normalized.*\\.txt\\.gz$", adjustBackground=TRUE, useSCAN=TRUE, scanConvThreshold=0.5, quantileNormalize=TRUE, log2Transform=TRUE, numCores=1, verbose=FALSE) {
  supplementaryFilePaths <- getSupplementaryFilesFromGEO(gseID)
  nonNormalizedFilePaths <- supplementaryFilePaths[grepl(nonNormalizedFilePattern, supplementaryFilePaths, ignore.case = TRUE)]

  if (length(nonNormalizedFilePaths) > 0) {
    # Convert Excel files to TSV, when necessary.
    for (i in 1:length(nonNormalizedFilePaths)) {
      f <- nonNormalizedFilePaths[i]

      if (grepl("\\.xlsx?$", f)) {
        stop("got here")
        excelData <- read_excel(f)
        write.table(excelData, paste0(f, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, na = "\"\"")
        nonNormalizedFilePaths[i] <- paste0(f, ".tsv")
        unlink(f)
      }
    }

    nonNormList <- readNonNormalizedData(nonNormalizedFilePaths)
  } else {
    nonNormList <- retrieveFromIDAT(gseID, supplementaryFilePaths)
  }

  annotationPackagePrefix <- getAnnotationPackagePrefixFromGEO(gseID)

  exprMatrix <- normalizeBeadChipData(nonNormList,
                                      annotationPackagePrefix,
                                      adjustBackground=adjustBackground,
                                      useSCAN=useSCAN, scanConvThreshold=scanConvThreshold,
                                      numCores=numCores,
                                      verbose=verbose)

  phenoData <- getPhenoDataFromGEO(gseID)

  gse <- getGEO(gseID, GSEMatrix=TRUE)
  experimentData <- experimentData(gse[[1]])
  
  featureData <- getFeatureData(exprMatrix, annotationPackagePrefix)
  featureData <- new("AnnotatedDataFrame", data = featureData)

  # Column names must match between sampleData and exprMatrix.
  colnames(exprMatrix) <- rownames(phenoData)
    
  # Row names must match between featureData and exprMatrix.
  exprMatrix <- exprMatrix[rownames(featureData), ]
  
  if (quantileNormalize) {
    exprMatrix <- normalizeBetweenArrays(exprMatrix, method = "quantile")
  }
  
  if (log2Transform) {
    exprMatrix <- doLog2(exprMatrix)
  }

  # https://www.bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf
  return(ExpressionSet(assayData = exprMatrix,
                       phenoData = phenoData, 
                       experimentData = experimentData,
                       featureData = featureData,
                       annotation = annotationPackagePrefix))
}

getSupplementaryFilesFromGEO <- function(gseID) {
  downloadDirPath <- paste0(tempdir(), "/", gseID)
  
  if (dir.exists(downloadDirPath)) {
    filePaths <- list.files(path = downloadDirPath, full.names = TRUE)
    
    if (length(filePaths) > 0) {
      return(filePaths)
    } else {
      unlink(downloadDirPath, recursive = TRUE, force = TRUE)
    }
  }
  
  getGEOSuppFiles(gseID, baseDir = tempdir())
  
  return(list.files(path = downloadDirPath, full.names = TRUE))
}

readNonNormalizedData <- function(filePaths, probeIDColumn=NULL, detectionPValueColumnPattern=NULL, exprColumnPattern=NULL, verbose=FALSE) {
  # Infer the file field delimiter
  delimiter <- inferFileDelimiter(filePaths[1])
  
  # Read the first line to find the column names
  originalColnames <- getColumnNames(filePaths, delimiter)

  # Auto-detect the probeid column.
  if (is.null(probeIDColumn)) {
    if (verbose) {
      message(paste0("Auto-detecting probeIDColumn."))
    }
    
    probeIDColumn <- originalColnames[1]
  }

  # Auto-detect detectionPValueColumnPattern.
  if (is.null(detectionPValueColumnPattern)) {
    if (verbose) {
      message(paste0("Auto-detecting detectionPValueColumnPattern."))
    }
    
    candidatePValueColnames <- originalColnames[grep("detection", originalColnames, ignore.case = TRUE)]

    detectionPValueColumnPattern <- findLongestCommonPrefix(candidatePValueColnames)

    if (detectionPValueColumnPattern == "") {
      detectionPValueColumnPattern <- findLongestCommonSuffix(candidatePValueColnames)
    }

    if (detectionPValueColumnPattern == "") {
      message(paste0("No value was specified for the detectionPValueColumnPattern parameter, and a pattern could not be auto-detected."))
      stop()
    }
  }
  
  # Auto-detect exprColumnPattern.
  if (is.null(exprColumnPattern)) {
    if (verbose) {
      message(paste0("Auto-detecting exprColumnPattern."))
    }

    candidateExprColnames <- originalColnames[grep("signal", originalColnames, ignore.case = TRUE)]
    
    if (length(candidateExprColnames) == 0) {
      detectionPIndices <- grep(detectionPValueColumnPattern, originalColnames, ignore.case = TRUE)
      candidateExprColnames <- detectionPIndices - 1
    }
    
    exprColumnPattern <- findLongestCommonPrefix(candidateExprColnames)
    
    if (exprColumnPattern == "") {
      exprColumnPattern <- findLongestCommonSuffix(candidateExprColnames)
    }
    
    if (exprColumnPattern == "") {
      message(paste0("No value was specified for the exprColumnPattern parameter, and a pattern could not be auto-detected."))
      stop()
    }
  }

print(probeIDColumn)
print(candidatePValueColnames)
print(candidateExprColnames)
print(detectionPValueColumnPattern)
print(exprColumnPattern)
print(filePaths)
print(delimiter)

filePaths = "/var/folders/5y/bv281l6d0m90njln4ywknhjm0000gq/T//Rtmp5yHeDY/GSE224309/GSE224309_Raw_Data.xlsx.tsv"
probeIDColumn = "PROBE_ID"
exprColumnPattern = "AVG_Signal"
detectionPValueColumnPattern = "Detection Pval"
delimiter = "\t"

  tmpFile <- read.table(filePaths[1], header = TRUE, sep = delimiter, fill = TRUE, check.names = FALSE)
  tmpFile <- read.columns(filePaths[1], required.col=NULL, text.to.search="", sep="\t", quote="\"", skip=0, fill=TRUE, blank.lines.skip=TRUE, comment.char="", allowEscapes=FALSE)
  #tmpFile <- read.table(filePaths[1], header = TRUE, sep = delimiter, check.names = FALSE)
  nonNormList <- read.ilmn(filePaths, probeid = probeIDColumn, expr = exprColumnPattern, other.columns = detectionPValueColumnPattern, sep=delimiter, text.to.search="", quote="\"")
print("got here")
stop()

  # When there are unusual delimiters, the read.ilmn function can't handle it well, and the number of columns is zero.
  # In this scenario, we read the file and then save it with tabs as delimiters. We also have to fix the column
  # names when R adds suffixes to them.
  if (ncol(nonNormList$E) == 0) {
    filePaths2 <- c()

    for (filePath in filePaths) {
      tmpFilePath <- paste0(tempfile(), ".tsv")

      raw <- read.table(filePath, header = TRUE, sep = "", quote = "\"", stringsAsFactors = FALSE, check.names = FALSE)
      colnames(raw) <- originalColnames
      stopifnot(length(originalColnames) == ncol(raw))

      # Save the data to temp file.
      write.table(raw, file = tmpFilePath, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = TRUE)

      filePaths2 <- c(filePaths2, tmpFilePath)
    }

    nonNormList <- read.ilmn(filePaths2, probeid = probeIDColumn, expr = exprColumnPattern, other.columns = detectionPValueColumnPattern, sep="\t")
  }

  return(nonNormList)
}

retrieveFromIDAT <- function(gseID, supplementaryFilePaths) {
  tarFilePath <- supplementaryFilePaths[grepl("_RAW\\.tar$", supplementaryFilePaths)]
  
  if (length(tarFilePath) == 0) {
    return(c())
  }

  contentsDirPath <- paste0(tempdir(), "/", gseID, "__files")
  untar(tarFilePath, exdir = contentsDirPath)

  for (f in list.files(path = contentsDirPath, pattern = "\\.idat\\.gz$", full.names = TRUE)) {
    gunzip(f, overwrite = TRUE)
  }

  idatFilePaths <- list.files(path = contentsDirPath, pattern = "\\.idat$", full.names = TRUE)

  if (length(idatFilePaths) == 0) {
    message(paste0("No IDAT files were found in the supplementary archive for ", gseID, "."))
    return(c())
  }

  bgxFilePath <- list.files(path = contentsDirPath, pattern = "_B\\.txt\\.gz", full.names = TRUE)
  bgxFilePath <- bgxFilePath[1] # There is often more than one, but we can just take the first one.

  rgSet <- read.idat(idatFilePaths, bgxfile = bgxFilePath)

  rownames(rgSet$E) <- rgSet$genes$Probe_Id

  rgSet$other$`Detection Pval` <- detectionPValues(rgSet)
  rownames(rgSet$other$`Detection Pval`) <- rgSet$genes$Probe_Id
  
  return(rgSet)
}

normalizeBeadChipData <- function(nonNormList, annotationPackagePrefix, adjustBackground=TRUE, useSCAN=TRUE, scanConvThreshold=0.5, numCores=1, verbose=FALSE) {
  platformPackageName <- paste0(annotationPackagePrefix, ".db")
  message(paste0("Loading package ", platformPackageName))
  library(platformPackageName, character.only=TRUE)
  
  probeSequenceRef <- getRefInfo(annotationPackagePrefix, "PROBESEQUENCE")
  probeQualityRef <- getRefInfo(annotationPackagePrefix, "PROBEQUALITY")
  
  # Remove low-quality probes
  perfectProbes = which(grepl("Perfect", probeQualityRef$ProbeQuality))
  goodProbes = which(grepl("Good", probeQualityRef$ProbeQuality))
  probesToKeep = probeQualityRef$IlluminaID[c(perfectProbes, goodProbes)]
  
  exprData <- nonNormList$E
  detectionPValues <- nonNormList$other$`Detection Pval`
  
  probesToKeep <- intersect(probesToKeep, rownames(exprData))
  exprData <- exprData[probesToKeep,]
  detectionPValues <- detectionPValues[probesToKeep,]
  
  if (adjustBackground) {
    if (min(exprData) < 0) {
      if (max(exprData) > 30) {
        message(paste0("There are negative values in the data, and there are relatively large positive values. This suggests that a background adjustment has already been performed and a log2 transformation has not. To avoid adding noise, we will not perform a background adjustment."))
      } else {
        message(paste0("There are negative values in the data, and there are relatively small positive values. This suggests that a background adjustment and/or log2 transformation have been performed. To avoid adding noise, we will not perform a background adjustment."))
      }
    } else {
      exprData <- backgroundCorrect(exprData, signalPValueData=detectionPValues)
    }
  }
  
  if (useSCAN) {
    if (any(exprData < 0, na.rm = TRUE)) {
      message("Negative values detected. Before SCAN can be applied, they must be shifted.")
      
      shiftAmount <- abs(min(exprData, na.rm = TRUE)) + 1
      exprData <- exprData + shiftAmount
    }
    
    if (max(exprData, na.rm = TRUE) < 30) {
      message("It appears the data have been log-transformed. Before SCAN can be applied, Negative values detected: before SCAN can be applied, we must reverse the log transformation.")
      exprData <- 2^exprData
    }
    
    rownames(probeSequenceRef) <- probeSequenceRef$IlluminaID
    probeSequences = probeSequenceRef[probesToKeep, 2]
    
    exprData <- scanNorm(exprData, probeSequences, convThreshold=scanConvThreshold, numCores=numCores, verbose=verbose)
  }
    
  return(exprData)
}

getFeatureData <- function(exprMatrix, annotationPackagePrefix) {
  platformPackageName <- paste0(annotationPackagePrefix, ".db")
  message(paste0("Loading package ", platformPackageName))
  library(platformPackageName, character.only=TRUE)

  # https://bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
  probeSequenceRef <- getRefInfo(annotationPackagePrefix, "PROBESEQUENCE", "Probe_Sequence")

  probeReporterRef <- getRefInfo(annotationPackagePrefix, "REPORTERGROUPNAME", "Probe_Reporter_Type")
  # Add rows for "regular" probes.
  regularProbes <- setdiff(rownames(exprMatrix), probeReporterRef[,1])
  probeReporterRef <- rbind(probeReporterRef, data.frame(Illumina_ID = regularProbes, Probe_Reporter_Type = "regular"))

  probeQualityRef <- getRefInfo(annotationPackagePrefix, "PROBEQUALITY", "Probe_Quality")
  probeChrRef <- getRefInfo(annotationPackagePrefix, "CHR", "Chromosome")
  probeEnsemblRef <- getRefInfo(annotationPackagePrefix, "ENSEMBLREANNOTATED", "Ensembl_Gene_ID")
  probeEntrezRef <- getRefInfo(annotationPackagePrefix, "ENTREZREANNOTATED", "Entrez_Gene_ID")
  probeSymbolRef <- getRefInfo(annotationPackagePrefix, "SYMBOLREANNOTATED", "Gene_Symbol")
  probeNuRef <- getRefInfo(annotationPackagePrefix, "NUID", "nuID")

  featureData <- merge(probeSequenceRef, probeReporterRef, by = "Illumina_ID", sort = FALSE)
  featureData <- merge(featureData, probeQualityRef, by = "Illumina_ID", sort = FALSE)
  featureData <- merge(featureData, probeChrRef, by = "Illumina_ID", sort = FALSE)
  featureData <- merge(featureData, probeEnsemblRef, by = "Illumina_ID", sort = FALSE)
  featureData <- merge(featureData, probeEntrezRef, by = "Illumina_ID", sort = FALSE)
  featureData <- merge(featureData, probeSymbolRef, by = "Illumina_ID", sort = FALSE)
  featureData <- merge(featureData, probeNuRef, by = "Illumina_ID", sort = FALSE)
  
  # In some cases, there are two chromosome values for the same probe. This
  # code collapses those and separates them with commas.
  featureData <- aggregate(. ~ Illumina_ID, data = featureData, FUN = function(x) {
    if (is.character(x)) {
      paste(sort(unique(x)), collapse = ",")
    } else {
      mean(x[1])
    }
  })

  featureData <- featureData[featureData$Illumina_ID %in% rownames(exprMatrix), ]
  rownames(featureData) <- featureData$Illumina_ID
  featureData$Illumina_ID <- NULL

  return(featureData)
}

inferFileDelimiter <- function(filePath, nLines = 5) {
  # Candidate delimiters
  delimiters <- c("," = ",", "\t" = "\t", ";" = ";", "|" = "|", ":" = ":", " " = " ")
  
  # Read first few lines
  lines <- readLines(filePath, n = nLines)
  
  # Clean up quotes (in case values are quoted)
  lines_cleaned <- gsub('(^"|"$)', "", lines)
  lines_cleaned <- gsub('"[[:space:]]+"', " ", lines_cleaned)
  
  # Count delimiter occurrences
  counts <- sapply(delimiters, function(delim) {
    mean(sapply(lines_cleaned, function(line) {
      length(strsplit(line, delim, fixed = TRUE)[[1]]) - 1
    }))
  })
  
  total_counts <- sum(counts)
  
  # Check if tab or comma make up at least 25% of all delimiters
  tab_fraction <- counts["\t"] / total_counts
  comma_fraction <- counts[","] / total_counts

  if (!is.na(tab_fraction) && tab_fraction >= 0.25) {
    return("\t")
  } else if (!is.na(comma_fraction) && comma_fraction >= 0.25) {
    return(",")
  } else {
    # Otherwise, pick the delimiter with the highest average count
    return(names(which.max(counts)))
  }
}

getColumnNames <- function(filePaths, fileFieldDelimiter) {
  allColNames <- character(0)  # Start with an empty vector
  
  for (i in seq_along(filePaths)) {
    # Read the first line
    headerItems <- read.table(filePaths[i], header = FALSE, sep = fileFieldDelimiter, 
                              nrows = 1, stringsAsFactors = FALSE, check.names = FALSE)
    colNames <- unlist(headerItems[1, ], use.names = FALSE)
    
    # Remove empty strings and NAs and trim whitespace
    colNames <- trimws(colNames[!is.na(colNames) & colNames != ""])
    
    if (i == 1) {
      allColNames <- colNames
    } else {
      # Drop the first column from all subsequent files
      allColNames <- c(allColNames, colNames[-1])
    }
  }
  
  return(allColNames)
}

getRefInfo <- function(annotationPackagePrefix, suffix, secondColumnName=NULL) {
  info <- get(paste0(annotationPackagePrefix, suffix))
  info <- info[mappedkeys(info)] # Returns the subset of mapped keys.

  df <- as.data.frame(info)
  
  if (!is.null(secondColumnName)) {
    colnames(df) <- c("Illumina_ID", secondColumnName)
  }

  return(df)
}

findLongestCommonPrefix <- function(x) {
  if (length(x) == 0)
    return("")

  # Split all strings into characters
  chars <- strsplit(x, "")
  
  # Get the minimum length across strings
  max_prefix_len <- min(sapply(chars, length))
  
  prefix <- character()
  
  for (i in seq_len(max_prefix_len)) {
    ith_chars <- sapply(chars, `[[`, i)
    if (length(unique(ith_chars)) == 1) {
      prefix <- c(prefix, ith_chars[1])
    } else {
      break
    }
  }
  
  paste0(prefix, collapse = "")
}

findLongestCommonSuffix <- function(x) {
  if (length(x) == 0)
    return("")
  
  # Split all strings into characters and reverse them
  chars <- lapply(strsplit(x, ""), rev)
  
  # Get the minimum length across strings
  max_suffix_len <- min(sapply(chars, length))
  
  suffix <- character()
  
  for (i in seq_len(max_suffix_len)) {
    ith_chars <- sapply(chars, `[[`, i)
    if (length(unique(ith_chars)) == 1) {
      suffix <- c(suffix, ith_chars[1])
    } else {
      break
    }
  }
  
  # Re-reverse the suffix to return it in the correct order
  paste0(rev(suffix), collapse = "")
}

# getNonNormalizedDataFromGEO <- function(gseID, fileSuffix="non_normalized") {
#   downloadDirPath = paste0(tempdir(), "/", gseID)
# 
#   if (!dir.exists(downloadDirPath)) {
#     dir.create(downloadDirPath)
#   }
#   
#   nonNormalizedURL = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gseID, "&format=file&file=", gseID, "_", fileSuffix, ".txt.gz")
# 
#   downloadFilePath = paste0(downloadDirPath, "/", fileSuffix, ".txt.gz")
#   
#   if (file.exists(downloadFilePath)) {
#     return(downloadFilePath)
#   }
#   
#   response <- HEAD(nonNormalizedURL)
#   
#   if (status_code(response) == 200) {
#     download.file(nonNormalizedURL, destfile = downloadFilePath, mode = "wb")
#     
#     return(downloadFilePath)
#   } else {
#     stop(paste0("A non-normalized GEO file could not be found using the standard URL (", nonNormalizedURL,") for ", gseID, ". You will need to provide a custom suffix."))
#   }
# }

# getIDATDataFromGEO <- function(gseID) {
#   
#   
#   tarFilePath <- list.files(path = paste0(tempdir(), "/", gseID), pattern = "\\.tar$|\\.tar\\.gz$", full.names = TRUE)
# 
#   # if (!file.exists(tarFilePath)) {
#   #   getGEOSuppFiles(gseID, baseDir = tempdir())
#   # }
# 
#   contentsDirPath <- paste0(tempdir(), "/", gseID, "_files")
#   untar(tarFilePath, exdir = contentsDirPath)
# 
#   idatFilePaths <- list.files(path = contentsDirPath, pattern = "\\.idat\\.gz$", full.names = TRUE)
# 
#   if (length(idatFilePaths) == 0) {
#     message(paste0("No gzipped IDAT files were found in the supplementary archive for ", gseID, "."))
#     stop()
#   }
#   
#   for (f in idatFilePaths) {
#     gunzip(f, overwrite = TRUE)
#   }
# 
#   idatFilePaths <- list.files(path = contentsDirPath, pattern = "\\.idat$", full.names = TRUE)
#   
#   if (length(idatFilePaths) == 0) {
#     message(paste0("No IDAT files were found in the supplementary archive after decompression for ", gseID, "."))
#     stop()
#   }
#   
#   bgxFilePath <- list.files(path = contentsDirPath, pattern = "_B\\.txt\\.gz", full.names = TRUE)
#   bgxFilePath <- bgxFilePath[1] # There is often more than one, but we can just take the first one.
# 
#   rgSet <- read.idat(idatFilePaths, bgxfile = bgxFilePath)
#   rawExprMatrix <- rgSet$E
#   rownames(rawExprMatrix) <- rgSet$genes$Probe_Id
#   
#   return(rawExprMatrix)
# }

getPhenoDataFromGEO <- function(gseID) {
  metadataESet <- getGEO(gseID, GSEMatrix = TRUE)

  phenotypeData <- pData(metadataESet[[1]])

  colsToExclude <- c()
  for (i in 1:ncol(phenotypeData)) {
    colValues <- phenotypeData[,i]
    colValues <- colValues[!is.na(colValues)]

    if (length(unique(colValues)) == 1) {
      colsToExclude <- c(colsToExclude, i)
    }
  }

  phenotypeData <- phenotypeData[, -colsToExclude, drop = FALSE]

  return(new("AnnotatedDataFrame", data = phenotypeData))
}

getAnnotationPackagePrefixFromGEO <- function(gseID, gplID=NULL) {
  # Get the GEO series
  gse <- getGEO(gseID, GSEMatrix = FALSE)

  if (is.null(gplID)) {
    # Extract platform(s)
    gplIDs <- Meta(gse)$platform_id
    
    if (length(gplIDs) > 1) {
      stop(paste0("Two platforms were found for ", gseID, ", so you must specify a value for gplID."))
    }
    
    gplID <- gplIDs[1]
  }

  # https://bioconductor.org/packages/release/data/annotation/
  # Platform / annotation package mapping
  platformDF <- read.table("BeadChip_Platforms.tsv", header = TRUE, sep = "\t", quote = "\"", 
                           stringsAsFactors = FALSE, check.names = FALSE)

  return(platformDF[which(platformDF$Accession == gplID),]$AnnotationPackagePrefix)
}

backgroundCorrect <- function(signalExprData, controlExprData=NULL, signalPValueData=NULL) {
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

scanNorm <- function(exprData, probeSequences, convThreshold=0.5, intervalN=10000, binsize=500, nbins=25, maxIt=100, asUPC=FALSE, numCores=1, verbose=FALSE) {
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
    normData <- as.matrix(scanNormVector(1, exprData[,1], mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose))
  }
  else
  {
    if (numCores > 1) {
      normData <- foreach(i = 1:ncol(exprData), .combine = cbind, .export=c("scanNormVector", "doLog2", "sampleProbeIndices", "EM_vMix", "mybeta", "assign_bin", "vsig", "vresp", "dn", "vbeta", "sig")) %dopar%
      {
        scanNormVector(i, exprData[,i], mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose)
      }
    } else {
      normData <- NULL

      for (i in 1:ncol(exprData)) {
        sampleNormData <- scanNormVector(i, exprData[,i], mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose)
        normData <- cbind(normData, sampleNormData)
      }
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

scanNormVector <- function(sampleIndex, my, mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose) {
  # Make an adjustment to accommodate any zero values (should happen rarely).
  my <- my + 1

  # Log transform the data, if needed.
  my = doLog2(my)

  # Add a tiny amount of random noise
  set.seed(0)
  noise = rnorm(length(my)) / 10000000
  my = my + noise
  
  nGroups = floor(length(my) / binsize)
  samplingProbeIndices = sampleProbeIndices(total=length(my), intervalN=intervalN, verbose=verbose)
  
  mixResult = EM_vMix(sampleIndex, y=my[samplingProbeIndices], X=mx[samplingProbeIndices,], nbins=nbins, convThreshold=convThreshold, maxIt=maxIt, verbose=verbose)
  
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

doLog2 <- function(x) {
  # This is a semi-crude way of checking whether the values were not previously log-transformed.
  if (max(x) > 30) {
    if (any(x <= 0)) {
      min_val <- min(x, na.rm = TRUE)
      shift <- abs(min_val) + 1  # small positive shift
      x <- x + shift
    }
    
    x <- log2(x)
  }

  return(x)
}

# reverseLog2 <- function(x) {
#   
# }

buildDesignMatrix = function(seqs, verbose=FALSE) {
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

sampleProbeIndices = function(total, intervalN, verbose=FALSE) {
  interval = floor(total / intervalN)
  if (interval <= 1)
    interval = 1
  
  seq(1, total, interval)
}

EM_vMix = function(sampleIndex, y, X, nbins, convThreshold=.01, maxIt=100, verbose=FALSE) {
  if (verbose)
    message(paste0("Starting EM for sample ", sampleIndex))

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
      message(paste0("Attempting to converge for sample ", sampleIndex, "...iteration ", it, ", c = ", round(conv, 6)))
  }
  
  if (verbose)
  {
    if (it == maxIt)
    {
      message(paste0("Reached convergence limit for sample ", sampleIndex, "...", it, " iterations. Proportion of background probes: ", round(p[1], 6)))
    } else {
      message(paste0("Converged for sample ", sampleIndex, " in ", it, " iterations. Proportion of background probes: ", round(p[1], 6)))
    }
  }
  
  list(p=p, b1=b1, b2=b2, s1=s1, s2=s2, bin=bin)
}

mybeta = function(y, X, gam, verbose=FALSE) {
  sqgam = sqrt(gam)
  Xw = sqgam * X
  yw = sqgam * y
  
  z = t(Xw) %*% Xw
  a = solve(z)
  
  b = a %*% t(Xw)
  as.numeric(b %*% yw)
}

assign_bin = function(y, nbins, verbose=FALSE) {
  quans = sort(y)[floor(length(y) * 1:nbins / nbins)]
  bins = sapply(y, function(x) { sum(x>quans) }) + 1
  
  # if (length(table(bins)) != nbins)
  # {
  #   if (verbose)
  #     message("The values were not separated into enough bins, so a tiny amount of noise will be added to make this possible.")
  #   
  #   set.seed(1)
  #   noise = rnorm(length(y)) / 10000000
  #   bins = assign_bin(y + noise, nbins, verbose)
  # }
  
  bins
}

vsig = function(y, X, b, gam, bin, nbins, verbose=FALSE) {
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

vresp = function(y, X, bin, p, b1, s1, b2, s2, verbose=FALSE) {
  vars0 = s1[bin]
  L0 = dn(y=y, m=(X %*% b1), s2=vars0, verbose=verbose)
  vars1 = s2[bin]
  L1 = dn(y=y, m=(X %*% b2), s2=vars1, verbose=verbose)
  
  gam1 = p[1] * L0 / (p[1] * L0 + p[2] * L1)
  gam2 = 1 - gam1
  cbind(gam1, gam2)
}

dn = function(y, m, s2, verbose=FALSE) {
  1 / (sqrt(2 * pi * s2)) * exp(-1 / (2 * s2) * (y - m)^2)
}

vbeta = function(y, X, bin, gam, s2, prof, verbose=FALSE) {
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

sig = function(y, m, verbose=FALSE) {
  resid = y - m
  sqrt((resid %*% resid) / length(y))
}