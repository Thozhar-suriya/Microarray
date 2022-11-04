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

## Differential Gene Expression ##
library(limma)
design <- model.matrix(~0+sampleInfo$group)
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Normal","Tumour")

summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

contrasts <- makeContrasts(Tumour - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable(fit2, coef=1)




