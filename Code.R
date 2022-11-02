## Microarray Analysis ##

## Loading Library ##
library(GEOquery)
library(Biobase)
library(dplyr)
library(pheatmap)


## Get GEO Soft file ##
my_id <- "GSE84437"
gse <- getGEO(my_id)

## if it Failed ##
# Download it manually and assign it using "download.file" ##
gse <- getGEO(filename = "GSE84437_series_matrix.txt.gz", GSEMatrix=TRUE)
show(gse)
length(gse)

gse <- gse[[1]] ## to analyze one dataset
## if two dataset we can use gse <- gse[[2]] ##

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

## Check the normalisation and scales used ##
summary(exprs(gse)) ## retrieve expression values
exprs(gse) <- log2(exprs(gse))

## Clinical data ##
sampleInfo <- pData(gse)
sampleInfo

## Pick out factors ##
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1)
sampleInfo <- rename(sampleInfo, group = source_name_ch1, patient=characteristics_ch1)

## Correlation Heatmap ##
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)





