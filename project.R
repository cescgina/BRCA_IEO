library(SummarizedExperiment)
library(geneplotter)
library(edgeR)
se <- readRDS(file = "seBRCA.rds")
## Visualizing data
mcols(se) # Visualize genes
metadata(se)$objectCreationDate
mcols(colData(se), use.names = TRUE) #Phenotypic information
colData(se)[1:5, 1:5]
rowRanges(se) #Feature data
names(se) <- rowRanges(se)$symbol
##
dge <- DGEList(counts = assays(se)$counts, genes = mcols(se))
names(dge) #Fields in DGEList
assays(se)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)# Calculate CPM and store it within se
assays(se)$logCPM[1:5, 1:5]
logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)

## Library sizes
sampledepth <- round(dge$sample$lib.size / 1e6, digits = 1)
names(sampledepth) <- substr(colnames(se), 6, 12)
types <- se$type
dge$samples$type <- types
ordering <- order(sampledepth)
sorteddepth <- sampledepth[ordering]
ordtype <- types[ordering]
# barplot(sorteddepth, col = c("red","blue")[as.integer(ordtype)], space = 0.1)
# plot(sorteddepth,ordtype, col=c("red","blue")[as.integer(ordtype)])
# legend("bottomleft", c("Normal", "Tumor"), fill = c("red", "blue"))
boxplot(sorteddepth ~ ordtype, col = c("red", "blue"), ylab = "Millions of reads")

## Distribution of expression levels
normal <- rownames(dge$samples[dge$samples$type == "normal",])
tumor <- rownames(dge$samples[dge$samples$type == "tumor",])
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(logCPM[,tumor])), xlab = "log2 CPM", 
             legend = NULL, main = "Tumor samples", cex.axis = 1.2, 
             cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(logCPM[,normal])), xlab = "log2 CPM", 
             legend = NULL, main = "Normal samples", cex.axis = 1.2, 
             cex.lab = 1.5, las = 1)

# ##Multidimensional plot
# Comented because it is very slow
# Also, given the huge amount of information is not very informative, maybe it
# is useful when working with a subset of the data
# plotMDS(dge, col = c("red", "blue")[as.integer(dge$samples$group)], cex = 0.7)
# legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)

##Filtering
# Filtering will most likely be necessary before proceeding with the analysis
cut <- 1
par(mfrow = c(1,1))
hist(rowMeans(logCPM), main = "")
abline(v = cut, col = "red", lwd = 2)
# Filter by logCPM > than 3, this is done solely with a testing purpose, cut-off is 
# arbritary and probably makes no sense
dgefilter <- dge[rowMeans(logCPM) > cut, ]

##TMM normalization
dge <- calcNormFactors(dge)
head(dge$samples$norm.factors)
assays(se)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)

##Quantile normalization
#TODO

## MA plot (ommited for now)
#par(mfrow = c(9, 3), mar = c(4, 5, 3, 1))
setmp <- se[, se$type == "normal"]
dgetmp <- dge[, se$type == "normal"]
for (i in 1:ncol(setmp)) {
  # For each sample, draw a scatter plot comparing the average logCPM values for 
  # each gene versus the difference of the logCPM values for each genes for that
  # sample and the average values
  A <- rowMeans(assays(setmp)$logCPM)
  M <- assays(setmp)$logCPM[, i] - A
  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
  smoothScatter(A, M, main = samplename, las = 1)
  abline(h = 0, col = "blue", lwd = 2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col = "red", lwd = 2)
}

#par(mfrow=c(22, 3), mar=c(4, 5, 3, 1))
setmp <- se[, se$type == "tumor"]
dgetmp <- dge[, se$type == "tumor"]
for (i in 1:ncol(setmp)) {
  A <- rowMeans(assays(setmp)$logCPM)
  M <- assays(setmp)$logCPM[, i] - A
  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
  smoothScatter(A, M, main = samplename, las = 1)
  abline(h = 0, col = "blue", lwd = 2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col = "red", lwd = 2)
}

## Batch identification
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
# tissue source sites (centers, for example 3C -> Columbia University)
tss <- substr(colnames(se), 6, 7)
table(tss)
# Sequencing center(07 -> University of North Carolina)
center <- substr(colnames(se), 27, 28)
table(center)
# plate: Order of plate in a sequence of 96-well plates 
plate <- substr(colnames(se), 22, 25)
table(plate)
# analyte: molecular specimen extracted from a portion. Usually of a particular molecular type; for example, total RNA, or whole genome amplified DNA. R -> RNA
portionanalyte <- substr(colnames(se), 18, 20)
table(portionanalyte)
# vial: portion of a sample
samplevial <- substr(colnames(se), 14, 16)
table(samplevial)

table(data.frame(TYPE = se$type, TSS = tss))
