# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# for (packageName in c("doParallel", "GEOquery", "httr", "limma", "oligo", "illuminaHumanv4.db", "readxl", "vsn")) {
#   BiocManager::install(packageName)
# }

library(doParallel)
library(GEOquery)
library(httr)
library(limma)
library(oligo)
library(readxl)
library(vsn)

normalizeBeadChipDataFromGEO <- function(gseID,
                                         nonNormalizedFilePattern="non.*normalized.*\\.txt\\.gz$",
                                         exprColumnPattern=NULL,
                                         correctBackgroundType=NULL,
                                         scanNormalize=TRUE, scanConvThreshold=0.5, scanIntervalN=50000, scanBinsize=500, scanNbins=25, scanMaxIt=100, scanAsUPC=FALSE,
                                         vsnNormalize=TRUE,
                                         quantileNormalize=FALSE,
                                         log2Transform=TRUE,
                                         numCores=1,
                                         verbose=FALSE) {
  #TODO: Make sure parameters are in valid range.
  supplementaryFilePaths <- getSupplementaryFilesFromGEO(gseID)
  nonNormalizedFilePaths <- supplementaryFilePaths[grepl(nonNormalizedFilePattern, supplementaryFilePaths, ignore.case = TRUE)]

  if (length(nonNormalizedFilePaths) > 0) {
    # Convert Excel files to TSV, when necessary.
    for (i in 1:length(nonNormalizedFilePaths)) {
      f <- nonNormalizedFilePaths[i]

      if (grepl("\\.xlsx?$", f)) {
        excelData <- read_excel(f)
        write.table(excelData, paste0(f, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, na = "\"\"")
        nonNormalizedFilePaths[i] <- paste0(f, ".tsv")
        unlink(f)
      }
    }

    nonNormList <- readNonNormalizedData(nonNormalizedFilePaths, exprColumnPattern = exprColumnPattern)
  } else {
    nonNormList <- retrieveFromIDAT(gseID, supplementaryFilePaths)
  }

  annotationPackagePrefix <- getAnnotationPackagePrefixFromGEO(gseID)

  #TODO: Invoke extractDataUsingAnnotations before invoking this function.
  exprMatrix <- normalizeBeadChipData(nonNormList,
                                      annotationPackagePrefix,
                                      correctBackgroundType=correctBackgroundType,
                                      scanNormalize=scanNormalize, scanConvThreshold=scanConvThreshold, scanIntervalN=scanIntervalN, scanBinsize=scanBinsize, scanNbins=scanNbins, scanMaxIt=scanMaxIt, scanAsUPC=scanAsUPC,
                                      vsnNormalize=vsnNormalize,
                                      quantileNormalize=quantileNormalize,
                                      log2Transform=log2Transform,
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

    candidateExprColnames <- originalColnames[grep("avg_signal", originalColnames, ignore.case = TRUE)]

    if (length(candidateExprColnames) == 0) {
      message(paste0("No value was specified for the exprColumnPattern parameter, and a pattern could not be auto-detected. You can view the downloaded file at ", filePaths[1], "."))
      stop()
      # detectionPIndices <- grep(detectionPValueColumnPattern, originalColnames, ignore.case = TRUE)
      # candidateExprColnames <- detectionPIndices - 1
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

  nonNormList <- read.ilmn(filePaths, probeid = paste0("^ *", probeIDColumn, " *$"), expr = exprColumnPattern, other.columns = detectionPValueColumnPattern, sep=delimiter)
    
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

  dp <- detectionPValues(rgSet)
  rgSet$other <- NULL
  rgSet$other$`Detection Pval` <- dp
  rownames(rgSet$other$`Detection Pval`) <- rgSet$genes$Probe_Id

  return(rgSet)
}

extractDataUsingAnnotations <- function(nonNormList,
                                        annotationPackagePrefix,
                                        verbose=FALSE) {
  platformPackageName <- paste0(annotationPackagePrefix, ".db")
  message(paste0("Loading package ", platformPackageName))
  library(platformPackageName, character.only=TRUE)
  
  probeSequenceRef <- getRefInfo(annotationPackagePrefix, "PROBESEQUENCE")
  probeQualityRef <- getRefInfo(annotationPackagePrefix, "PROBEQUALITY")

  # TODO: Use this to find which are control probes and separate based on that.
  # probeReporterRef <- getRefInfo(annotationPackagePrefix, "REPORTERGROUPNAME", "Probe_Reporter_Type")
  # # Add rows for "regular" probes.
  # regularProbes <- setdiff(rownames(exprMatrix), probeReporterRef[,1])
  # probeReporterRef <- rbind(probeReporterRef, data.frame(Illumina_ID = regularProbes, Probe_Reporter_Type = "regular"))
  
  # Remove low-quality probes
  perfectProbes = which(grepl("Perfect", probeQualityRef$ProbeQuality))
  goodProbes = which(grepl("Good", probeQualityRef$ProbeQuality))
  probesToKeep = probeQualityRef$IlluminaID[c(perfectProbes, goodProbes)]
  
  exprData <- nonNormList$E
  #detectionPValues <- nonNormList$other$`Detection Pval`
  detectionPValues <- nonNormList$other[[1]] # Sometimes the names are not consistent, but detection p-values should be the only thing in "other".
  
  probesToKeep <- intersect(probesToKeep, rownames(exprData))
  exprData <- exprData[probesToKeep,]
  detectionPValues <- detectionPValues[probesToKeep,]

  # rownames(probeSequenceRef) <- probeSequenceRef$IlluminaID
  # probeSequences = probeSequenceRef[probesToKeep, 2]
  
  return(exprData)
}

normalizeBeadChipData <- function(signalExprData,
                                  signalProbeSequences,
                                  controlExprData=NULL,
                                  controlProbeSequences=NULL,
                                  detectionPValues=NULL,
                                  correctBackgroundType=NULL,
                                  vsnNormalize=FALSE,
                                  quantileNormalize=FALSE,
                                  scanNormalize=TRUE, scanConvThreshold=0.5, scanIntervalN=50000, scanBinsize=500, scanNbins=25, scanMaxIt=100, scanAsUPC=FALSE,
                                  log2Transform=TRUE,
                                  numCores=1,
                                  verbose=FALSE) {
  #############################################
  # Check parameters
  #############################################
  
  if (!is.matrix(signalExprData))
    stop("signalExprData must be a matrix.")
  
  if (is.null(detectionPValues)) {
    if (is.null(controlExprData)) {
      stop("If detectionPValues is not provided, controlExprData must be provided.")
    } else {
      if (!is.matrix(controlExprData))
        stop("controlExprData must be a matrix.")
      
      if (ncol(signalExprData) != ncol(controlExprData))
        stop("The dimensions of signalExprData and controlExprData are incompatible.")
    }
  } else {
    if (!is.matrix(detectionPValues))
      stop("detectionPValues must be a matrix.")
    
    if (all(dim(signalExprData) != dim(detectionPValues)))
      stop("The dimensions of signalExprData and detectionPValues are incompatible.")
  }

  if (is.null(signalProbeSequences)) {
    stop("A value was not provided for signalProbeSequences.")
  } else {
    if (!is.vector(signalProbeSequences))
      stop("signalProbeSequences must be a vector.")
    if (nrow(signalExprData) != length(signalProbeSequences))
      stop("The dimensions of signalExprData and signalProbeSequences must be identical.")
  }
  
  if (!is.null(controlExprData)) {
    if (is.null(controlProbeSequences))
      stop("A value was provided for controlExprData but not for controlProbeSequences.")
    if (!is.vector(controlProbeSequences))
      stop("controlProbeSequences must be a vector.")
    if (nrow(controlExprData) != length(controlProbeSequences))
      stop("The dimensions of controlExprData and controlProbeSequences must be identical.")
  }

  # Make sure all parameters are in valid range.
  if (!(correctBackgroundType %in% c(NULL, "controls", "detectionP")))
    stop("Invalid value for correctBackgroundType.")
  if (!is.logical(vsnNormalize))
    stop("vsnNormalize must be a logical value.")
  if (!is.logical(quantileNormalize))
    stop("quantileNormalize must be a logical value.")
  if (!is.logical(scanNormalize))
    stop("scanNormalize must be a logical value.")
  if (!is.logical(log2Transform))
    stop("log2Transform must be a logical value.")
  if (!is.numeric(scanConvThreshold) | scanConvThreshold < 0.01 | scanConvThreshold > 100)
    stop("Invalid value for scanConvThreshold.")
  if (!is.numeric(scanIntervalN) | !(as.integer(scanIntervalN) = scanIntervalN) | scanIntervalN < 1000 | scanIntervalN > 50000)
    stop("Invalid value for scanIntervalN.")
  if (!is.numeric(scanBinsize) | !(as.integer(scanBinsize) = scanBinsize) | scanBinsize < 50 | scanBinsize > 5000)
    stop("Invalid value for scanBinsize.")
  if (!is.numeric(scanNbins) | !(as.integer(scanNbins) = scanNbins) | scanNbins < 10 | scanNbins > 50)
    stop("Invalid value for scanNbins.")
  if (!is.numeric(scanMaxIt) | !(as.integer(scanMaxIt) = scanMaxIt) | scanMaxIt < 10 | scanMaxIt > 10000)
    stop("Invalid value for scanMaxIt.")
  if (!is.logical(scanAsUPC))
    stop("scanAsUPC must be a logical value.")
  if (!is.numeric(numCores) | !(as.integer(numCores) = numCores) | numCores < 1 | numCores > 10000)
    stop("Invalid value for numCores")
  if (!is.logical(verbose))
    stop("verbose must be a logical value.")
  
  if (!is.null(correctBackgroundType)) {
    if (min(signalExprData) < 0) {
      if (max(signalExprData) > 30) {
        message(paste0("There are negative values in the data, and there are relatively large positive values. This suggests that a background adjustment has already been performed and a log2 transformation has not. To avoid adding noise, we will not perform a background correction."))
      } else {
        message(paste0("There are negative values in the data, and there are relatively small positive values. This suggests that a background adjustment and/or log2 transformation have been performed. To avoid adding noise, we will not perform a background correction."))
      }
    } else {
      if (correctBackgroundType == "controls") {
        if (is.null(controlExprData)) {
          stop("A NULL value was specified for controlExprData, but 'controls' was specified for correctBackgroundType.")
        }

        correctedData <- performBackgroundCorrection(signalExprData, controlExprData=controlExprData)
        signalExprData <- correctedData[1:nrow(signalExprData),]
        controlExprData <- correctedData[(nrow(signalExprData) + 1):nrow(correctedData),]
      } else {
        if (correctBackgroundType == "detectionP") {
          if (is.null(detectionPValues)) {
            stop("A NULL value was specified for detectionPValues, but 'detectionP' was specified for correctBackgroundType.")
          }

          signalExprData <- performBackgroundCorrection(signalExprData, detectionPValues=detectionPValues)
        }
      }
    }
  }

  if (vsnNormalize) {
    # We are going directly to the vsn2 function so that background correction is not also applied.
    # We are not including the control probes in the VSN normalization to avoid distorting the signal.
    fit = vsn2(signalExprData)
    signalExprData = predict(fit, newdata=signalExprData)
    signalExprData <- 2^signalExprData # Reverse the log2 transformation.
  }

  if (quantileNormalize) {
    exprData <- signalExprData
    
    if (!is.null(controlExprData)) {
      exprData <- rbind(exprData, controlExprData)
    }

    exprData <- normalizeBetweenArrays(exprData, method = "quantile")
    signalExprData <- exprData[1:nrow(signalExprData),]
    
    if (!is.null(controlExprData)) {
      controlExprData <- exprData[(nrow(signalExprData) + 1):nrow(exprData),]
    }
  }

  if (scanNormalize) {
    exprData <- signalExprData
    probeSequences <- signalProbeSequences
    
    if (!is.null(controlExprData)) {
      exprData <- rbind(exprData, controlExprData)
      probeSequences <- c(probeSequences, controlProbeSequences)
    }
    
    if (any(exprData < 0, na.rm = TRUE)) {
      message("Negative values detected. Before SCAN can be applied, they must be shifted.")

      shiftAmount <- abs(min(exprData, na.rm = TRUE)) + 1
      exprData <- exprData + shiftAmount
    }

    if (max(exprData, na.rm = TRUE) < 30) {
      message("It appears the data have been log-transformed. Before SCAN can be applied, Negative values detected: before SCAN can be applied, we must reverse the log transformation.")
      exprData <- 2^exprData
    }

    exprData <- scanNorm(exprData, probeSequences, scanConvThreshold, scanIntervalN, scanBinsize, scanNbins, scanMaxIt, scanAsUPC, numCores, verbose)

    signalExprData <- exprData[1:nrow(signalExprData),]
  }

  if (log2Transform) {
    signalExprData <- doLog2(signalExprData)
  }
    
  return(signalExprData)
}

getFeatureData <- function(exprMatrix, annotationPackagePrefix) {
  platformPackageName <- paste0(annotationPackagePrefix, ".db")
  message(paste0("Loading package ", platformPackageName))
  library(platformPackageName, character.only=TRUE)

  # https://bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
  probeSequenceRef <- getRefInfo(annotationPackagePrefix, "PROBESEQUENCE", "Probe_Sequence")

  # TODO: Remove this info from here because it is used elsewhere for filtering?
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

performBackgroundCorrection <- function(signalExprData, controlExprData=NULL, detectionPValues=NULL) {
  #############################################
  # Perform background correction (limma).
  #   It models the observed signal as a combination of: Observed_Signal = True_Signal + Background_Noise
  #   It removes non-specific fluorescence, helping improve the accuracy of low-intensity probe measurements.
  #   https://academic.oup.com/nar/article/38/22/e204/1049223
  #############################################

  status <- rep("regular", nrow(signalExprData))

  if (is.null(controlExprData)) { # We have detection p-values only.
    exprData <- nec(x = signalExprData, status = status, detection.p = detectionPValues)
  } else { # We have control values.
    status <- c(status, rep("negative", nrow(controlExprData)))
    exprData <- nec(x = rbind(signalExprData, controlExprData), status = status)
  }
  
  return(exprData)
}

#FYI: Don't call this function directly.
scanNorm <- function(exprData, probeSequences, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, numCores, verbose) {
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

  return(normData)
}

# This function performs per-sample normalization of gene expression microarray data
# using a two-component mixture model (via EM), adjusting for GC bias or similar probe-specific effects.
# It returns either:
#   - A normalized expression vector (after bin-wise standardization and back-transform), or
#   - A UPC-like value (probability of being "active") if `asUPC = TRUE`.
#
# Inputs:
# - sampleIndex: index of the sample (for logging/debugging purposes).
# - my: raw expression vector for the current sample (1 value per probe).
# - mx: design matrix containing features per probe (e.g., GC content).
# - convThreshold: convergence threshold for the EM algorithm.
# - intervalN: number of evenly spaced probes to sample for EM fitting.
# - binsize: number of probes per bin for normalization.
# - nbins: number of bins to use for estimating variance profiles.
# - maxIt: maximum iterations for the EM algorithm.
# - asUPC: if TRUE, return posterior probability of being "active" (from EM); if FALSE, return normalized values.
# - verbose: if TRUE, print diagnostic messages.
scanNormVector <- function(sampleIndex, my, mx, convThreshold, intervalN, binsize, nbins, maxIt, asUPC, verbose) {
  # Add +1 to all expression values to avoid log(0); should be rare.
  my <- my + 1
  
  # Apply log2 transformation to the expression data.
  my <- doLog2(my)
  
  # Add small noise to break ties and avoid instability in downstream steps.
  set.seed(0)
  noise <- rnorm(length(my)) / 1e7
  my <- my + noise
  
  # Decide how many bins will be used for normalization.
  nGroups <- floor(length(my) / binsize)
  
  # Select a subset of probes evenly across the array for model fitting.
  samplingProbeIndices <- sampleProbeIndices(total = length(my), intervalN = intervalN, verbose = verbose)

  # Fit a 2-component mixture model using the sampled probes.
  mixResult <- EM_vMix(
    sampleIndex,
    y = my[samplingProbeIndices],
    X = mx[samplingProbeIndices, ],
    nbins,
    convThreshold,
    maxIt,
    verbose
  )
  
  # Compute the fitted values from both components for all probes.
  m1 <- mx %*% mixResult$b1  # Background component
  m2 <- mx %*% mixResult$b2  # Signal component (unused here, but computed)
  
  # Initialize normalized output vector.
  index <- order(m1)
  y_norm <- rep(0, length(my))
  
  # Standardize expression values within bins based on fitted background model.
  for (i in 1:nGroups) {
    start_idx <- (binsize * (i - 1)) + 1
    end_idx <- min(binsize * i, length(my))
    tmp <- index[start_idx:end_idx]
    
    # Estimate standard deviation in the bin based on residuals from m1.
    tmpSd <- as.vector(sig(y = my[tmp], m = m1[tmp], verbose = verbose))
    
    # Standardize expression relative to background model.
    y_norm[tmp] <- (my[tmp] - m1[tmp]) / tmpSd
  }
  
  # Recompute bin assignments based on m1.
  bin <- assign_bin(y = m1, nbins = nbins, verbose = verbose)

  # Compute posterior probabilities (gam) from EM model for all probes.
  # Only keep the second column, which represents "signal" probability.
  gam <- vresp(
    y = my,
    X = mx,
    bin = bin,
    p = mixResult$p,
    b1 = mixResult$b1,
    s1 = mixResult$s1,
    b2 = mixResult$b2,
    s2 = mixResult$s2,
    verbose = verbose
  )[, 2]

  y_norm <- round(y_norm, 8)
  y_norm <- 2^y_norm  # Reverse the log2 transformation.
  
  gam <- round(gam, 8)
  
  # Return either UPC-like probabilities or normalized values.
  if (asUPC) {
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

buildDesignMatrix <- function(seqs, verbose = FALSE) {
  # Generate the initial one-hot encoded design matrix from sequences
  # This matrix encodes the presence of A, C, G, and T at each position
  mx <- sequenceDesignMatrix(seqs)
  
  # Estimate the sequence length based on the number of one-hot encoded columns
  # Assumes four columns per base position (A_, C_, G_, T_)
  seq_length <- ncol(mx) / 4
  
  # Compute the total counts of each nucleotide across each sequence
  numA <- rowSums(mx[, grepl("^A_", colnames(mx))])
  numC <- rowSums(mx[, grepl("^C_", colnames(mx))])
  numG <- rowSums(mx[, grepl("^G_", colnames(mx))])
  
  # Derive the count of T by subtracting A + C + G from total length
  numT <- seq_length - (numA + numC + numG)
  
  # Combine features into a final matrix:
  # - numT as a separate feature
  # - the full one-hot encoded matrix (mx)
  # - squared counts of A, C, G, and T to capture potential nonlinear effects
  features <- cbind(
    numT,
    mx,
    numA^2,
    numC^2,
    numG^2,
    numT^2
  )
  
  # Convert all matrix entries to integer type for efficiency and compatibility
  storage.mode(features) <- "integer"
  
  # If requested, print summary information about the resulting matrix
  if (verbose) {
    message("Built design matrix with ", nrow(features), " sequences and ", ncol(features), " features.")
  }
  
  # Return the complete feature matrix for modeling
  return(features)
}

# This function returns a set of evenly spaced probe indices for sampling.
# It is useful in scenarios (e.g., mixture model initialization)
# where you want to work with a subset of probes that are representative of the
# full dataset, rather than using all probes, for efficiency.
sampleProbeIndices <- function(total, intervalN, verbose = FALSE) {
  # Compute the sampling interval: how many steps to skip between indices
  interval <- floor(total / intervalN)
  
  # Ensure the interval is at least 1 to avoid infinite loops or zero-length sequences
  if (interval <= 1)
    interval <- 1
  
  # Generate a sequence of indices from 1 to total using the computed interval
  return(seq(1, total, interval))
}

# EM_vMix:
# This function performs Expectation-Maximization (EM) for a two-component mixture model,
# where each component models gene expression as a linear regression with bin-specific variance.
# The function estimates mixture proportions (p), regression coefficients (b1, b2),
# and bin-specific variances (s1, s2) for each component.
EM_vMix <- function(sampleIndex, y, X, nbins, convThreshold, maxIt, verbose) {
  if (verbose)
    message(paste0("Starting EM for sample ", sampleIndex))

  # Initialize responsibilities using a median-based split.
  quan <- sort(y)[floor(0.5 * length(y)) - 1]
  gam <- cbind(as.integer(y <= quan), as.integer(y > quan))

  # Initialize mixture proportions.
  p <- apply(gam, 2, mean)
  
  # Initialize regression parameters using weighted least squares.
  b1 <- mybeta(y = y, X = X, gam = gam[,1], verbose = verbose)
  b2 <- mybeta(y = y, X = X, gam = gam[,2], verbose = verbose)

  # Assign bins to observations for variance modeling.
  bin <- assign_bin(y = y, nbins = nbins, verbose = verbose)

  # Initialize bin-specific variances for both components.
  s1 <- vsig(y = y, X = X, b = b1, gam = gam[,1], bin = bin, nbins = nbins, verbose = verbose)
  s2 <- vsig(y = y, X = X, b = b2, gam = gam[,2], bin = bin, nbins = nbins, verbose = verbose)
  
  # Store initial parameter vector for convergence checking.
  theta_old <- c(p, b1, s1, b2, s2)

  it <- 0
  conv <- 1e6  # Arbitrarily large initial value to start the loop

  # EM loop: iterate until convergence or maximum iterations reached.
  while (conv > convThreshold & it < maxIt) {
    # E-step: compute responsibilities (posterior probabilities).
    gam <- vresp(y = y, X = X, bin = bin, p = p, b1 = b1, s1 = s1, b2 = b2, s2 = s2, verbose = verbose)

    # M-step: update parameters using responsibilities.
    p <- apply(gam, 2, mean)  # Update mixture proportions.
    
    # Update regression coefficients.
    b1 <- vbeta(y = y, X = X, bin = bin, gam = gam[,1], s2 = s1, prof = TRUE, verbose = verbose)
    bin <- assign_bin(y = (X %*% b1), nbins = nbins, verbose = verbose)  # Re-bin based on updated prediction
    b2 <- vbeta(y = y, X = X, bin = bin, gam = gam[,2], s2 = s2, prof = FALSE, verbose = verbose)
    
    # Update bin-specific variances.
    s1 <- vsig(y = y, X = X, b = b1, gam = gam[,1], bin = bin, nbins = nbins, verbose = verbose)
    s2 <- vsig(y = y, X = X, b = b2, gam = gam[,2], bin = bin, nbins = nbins, verbose = verbose)

    # Convergence check.
    theta <- c(p, b1, s1, b2, s2)
    conv <- max(abs(theta - theta_old) / (abs(theta_old) + 1e-8))  # Add small constant to avoid divide-by-zero
    theta_old <- theta
    it <- it + 1
    
    # Optional progress message.
    if (verbose)
      message(paste0("Attempting to converge for sample ", sampleIndex, "...iteration ", it, ", c = ", round(conv, 6)))
  }
  
  # Final status message.
  if (verbose) {
    if (it == maxIt) {
      message(paste0("Reached convergence limit for sample ", sampleIndex, "...", it, " iterations. Proportion of background probes: ", round(p[1], 6)))
    } else {
      message(paste0("Converged for sample ", sampleIndex, " in ", it, " iterations. Proportion of background probes: ", round(p[1], 6)))
    }
  }
  
  # Return all final parameter estimates.
  return(list(p = p, b1 = b1, b2 = b2, s1 = s1, s2 = s2, bin = bin))
}

mybeta = function(y, X, gam, verbose = FALSE) {
  # Take the square root of the weights (gam) for weighted least squares.
  # This allows expressing the weighted regression as a transformed ordinary least squares.
  sqgam = sqrt(gam)
  
  # Apply weights to each row of X and y.
  # Equivalent to pre-multiplying both sides of the regression by the square root of the weight matrix.
  Xw = sqgam * X  # Element-wise multiplication; rows of X are weighted.
  yw = sqgam * y  # Response vector also weighted.
  
  # Compute the weighted normal equations: (XᵗWX)⁻¹.
  z = t(Xw) %*% Xw
  
  # Solve for the inverse (XᵗWX)⁻¹.
  a = solve(z)
  
  # Compute the weighted least squares estimator:
  # β̂ = (XᵗWX)⁻¹ XᵗWy
  b = a %*% t(Xw)
  
  # Multiply with weighted response to get final coefficient estimates.
  return(as.numeric(b %*% yw))  # Returns β̂ as a numeric vector
}

assign_bin <- function(y, nbins, verbose = FALSE) {
  # Compute quantile thresholds to divide the values into nbins equal-sized bins.
  # Use type = 1 for consistency with floor-based quantile slicing.
  quans <- quantile(y, probs = seq(0, 1, length.out = nbins + 1), type = 1)
  
  # Assign each value in y to a bin using vectorized method.
  # This returns integers from 1 to nbins.
  bins <- findInterval(y, quans, rightmost.closed = TRUE, all.inside = TRUE)
  
  # If bin assignment resulted in fewer than nbins unique bins, add small noise and retry.
  if (length(unique(bins)) < nbins) {
    if (verbose) {
      message("Fewer than ", nbins, " bins assigned. Adding small noise to break ties.")
    }
    set.seed(1)
    noise <- rnorm(length(y), mean = 0, sd = 1e-7)
    y_noisy <- y + noise
    
    # Recompute quantiles and bins using noisy values.
    quans <- quantile(y_noisy, probs = seq(0, 1, length.out = nbins + 1), type = 1)
    bins <- findInterval(y_noisy, quans, rightmost.closed = TRUE, all.inside = TRUE)
  }
  
  # Return bin assignments.
  return(bins)
}

vsig <- function(y, X, b, gam, bin, nbins, verbose = FALSE) {
  # Initialize empty vector to store variances for each bin.
  s2 <- NULL
  
  # Loop over each bin.
  for (i in 1:nbins) {
    # Extract subset of y, X, and gam corresponding to bin i.
    ystar <- y[bin == i]              # Response values in bin i.
    Xstar <- X[bin == i, ]            # Design matrix rows in bin i.
    gamstar <- gam[bin == i] + 0.01   # Add small constant to gam to avoid zero weight.
    
    # Compute residuals: difference between observed and predicted values.
    resid <- as.numeric(ystar - Xstar %*% b)
    
    # Compute weighted residual variance for bin i.
    # Formula: (∑ wᵢ * rᵢ²) / ∑ wᵢ
    var_i <- ((resid * gamstar) %*% resid) / sum(gamstar)
    
    # Append result to output vector.
    s2 <- c(s2, var_i)
  }
  
  # Return vector of bin-specific variances.
  return(s2)
}

vresp <- function(y, X, bin, p, b1, s1, b2, s2, verbose = FALSE) {
  # Get bin-specific variances for each observation under component 1.
  # s1 is a vector of length nbins, so we map bin indices to variance values.
  vars0 <- s1[bin]
  
  # Compute likelihood under component 1 (e.g., background).
  # dn() computes the probability density assuming Gaussian noise.
  L0 <- dn(y = y, m = (X %*% b1), s2 = vars0, verbose = verbose)
  
  # Repeat for component 2 (e.g., signal).
  vars1 <- s2[bin]
  L1 <- dn(y = y, m = (X %*% b2), s2 = vars1, verbose = verbose)
  
  # Compute posterior probabilities using Bayes’ Rule.
  # gam1 = P(component 1 | y)
  gam1 <- p[1] * L0 / (p[1] * L0 + p[2] * L1)
  
  # gam2 = P(component 2 | y).
  gam2 <- 1 - gam1
  
  # Return a matrix of responsibilities (rows = observations, cols = components).
  return(cbind(gam1, gam2))
}

dn <- function(y, m, s2, verbose = FALSE) {
  # Compute the normalizing constant of the Gaussian PDF.
  # This is 1 / sqrt(2πσ²) for each observation.
  norm_const <- 1 / sqrt(2 * pi * s2)
  
  # Compute the exponent part: (y - m)^2 / (2σ²).
  # This represents the squared deviation scaled by the variance.
  exponent <- -1 / (2 * s2) * (y - m)^2
  
  # Combine parts to compute the PDF for each y[i].
  # This gives the likelihood of observing y[i] given mean m[i] and variance s2[i].
  return(norm_const * exp(exponent))
}

vbeta <- function(y, X, bin, gam, s2, prof, verbose = FALSE) {
  # Get standard deviation per observation based on bin assignments.
  vars <- sqrt(s2[bin])  # Standard deviation vector for each sample.
  
  # Compute the square root of responsibilities (for weighting).
  sqgam <- sqrt(gam)
  
  # Combine weights from variance and responsibility.
  # This is the weight applied to each observation.
  vars_sqgam <- vars * sqgam
  
  # Apply weights to design matrix and response vector.
  # These represent a transformation for weighted least squares.
  Xw <- (sqgam / vars) * X  # Weighted design matrix
  yw <- (sqgam / vars) * y  # Weighted response vector
  
  # Compute the normal equations for weighted least squares.
  tXw <- t(Xw)  # Transpose
  tXwXw <- tXw %*% Xw  # XᵗWX
  stXwXw <- solve(tXwXw)  # Inverse of XᵗWX
  stXwXwtXw <- stXwXw %*% tXw  # (XᵗWX)⁻¹ Xᵗ
  
  # Estimate regression coefficients β.
  result <- stXwXwtXw %*% yw
  
  # Return β as a vector.
  return(result)
}

sig <- function(y, m, verbose = FALSE) {
  # Compute residuals between observed and predicted values
  resid <- y - m
  
  # Compute standard deviation of residuals.
  # This is √(∑(resid²) / n), the root mean square error.
  return(sqrt((resid %*% resid) / length(y)))
}