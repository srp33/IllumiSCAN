#!/bin/bash

mkdir -p tmp Output

imageName="srp33/illumiscan"
numCores=3

docker build --rm -t $imageName -f Dockerfile .

function normalize {
  gseID="$1"
  probeIDColumn="$2"
  exprColumnPattern="$3"
  detectionPValueColumnPattern="$4"
  platform="$5"
  numCores="$6"
  url="$7"

  docker run --rm \
      -v $(pwd)/tmp:/tmp \
      -v $(pwd)/Output:/Output \
      $imageName \
      Rscript --vanilla normalizeBeadChip.R "$gseID" "$probeIDColumn" "$exprColumnPattern" "$detectionPValueColumnPattern" "$platform" "$numCores" "$url"
}

normalize GSE22427 ID_REF LV- "" illuminaHumanv4 $numCores "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE22427&format=file&file=GSE22427%5Fnon%2Dnormalized%2Etxt%2Egz"
#normalize GSE56457 "" SEQC_ILM_ "" illuminaHumanv4 $numCores "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE56457&format=file&file=GSE56457%5FSEQC%5FILM%5FGEOSub%2Etxt%2Egz"
#normalize GSE85452 "" "" "Detection Pval" illuminaHumanv4 $numCores "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE85452&format=file&file=GSE85452%5FNon%2Dnormalized%5Fdata%2Etxt%2Egz"
exit

normalize GSE67057 ID_REF "_[A-Z]$" "" illuminaHumanv2 $numCores "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67057&format=file&file=GSE67057%5Fnon%2Dnormalized%5Fdata%2Etxt%2Egz"
normalize GSE60102 ID_REF "^17\d\d$" "" illuminaHumanv2 $numCores "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60102&format=file&file=GSE60102%5Fnon%2Dnormalized%2Etxt%2Egz"
normalize GSE14847 "" "" "" illuminaHumanv2 $numCores "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14847&format=file&file=GSE14847%2Etxt%2Egz"
