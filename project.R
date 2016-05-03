library(SummarizedExperiment)
library(geneplotter)
library(edgeR)

se <- readRDS(file = "seBRCA.rds")

## Visualizing data
mcols(se) # Visualize genes
metadata(se)$objectCreationDate
mcols(colData(se), use.names = TRUE) #Phenotypic information
rowRanges(se) #Feature data
names(se) <- rowRanges(se)$symbol

## Subsetting, only samples that we have data for normal and tumor for same patient
codesnorm <- colData(se)[colData(se)$type == "normal",]$bcr_patient_barcode
codestum <- colData(se)[colData(se)$type == "tumor",]$bcr_patient_barcode
commoncodes <- codesnorm[codesnorm %in% codestum]
se <- se[,colData(se)$bcr_patient_barcode %in% commoncodes]

##
dge <- DGEList(counts = assays(se)$counts, genes = mcols(se))
names(dge) #Fields in DGEList
assays(se)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)# Calculate CPM and store it within se
assays(se)$logCPM[1:5, 1:5]
logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)

## Library sizes
# Library size is the sum of all counts for a sample
sampledepth <- round(dge$sample$lib.size / 1e6, digits = 1)
names(sampledepth) <- substr(colnames(se), 6, 12)
types <- se$type
dge$samples$type <- types
ordering <- order(sampledepth)
sorteddepth <- sampledepth[ordering]
ordtype <- types[ordering]
 barplot(sorteddepth, col = c("red","blue")[as.integer(ordtype)], space = 0.1, ylab = "Millions of reads", xlab = "Sample")
# plot(sorteddepth,ordtype, col=c("red","blue")[as.integer(ordtype)])
 legend("topleft", c("Normal", "Tumor"), fill = c("red", "blue"))
# boxplot(sorteddepth ~ ordtype, col = c("red", "blue"), ylab = "Millions of reads")

## Distribution of expression levels
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(logCPM[,se$type == "tumor"])), xlab = "log2 CPM", 
             legend = NULL, main = "Tumor samples", cex.axis = 1.2, 
             cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(logCPM[,se$type == "normal"])), xlab = "log2 CPM", 
             legend = NULL, main = "Normal samples", cex.axis = 1.2, 
             cex.lab = 1.5, las = 1)

##Filtering
# Filtering will most likely be necessary before proceeding with the analysis
avgexp <- rowMeans(logCPM)
cut <- 1
par(mfrow = c(1,1))
hist(avgexp, main = "")
abline(v = cut, col = "red", lwd = 2)
# Filter by logCPM > than 3, this is done solely with a testing purpose, cut-off is 
# arbritary and probably makes no sense
mask <- avgexp > 1
se <- se[mask, ]
dge <- dge[mask, ]

##TMM normalization
dge <- calcNormFactors(dge)
head(dge$samples$norm.factors)
assays(se)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)

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
  a <- rowmeans(assays(setmp)$logcpm)
  m <- assays(setmp)$logcpm[, i] - a
  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
  smoothscatter(a, m, main = samplename, las = 1)
  abline(h = 0, col = "blue", lwd = 2)
  lo <- lowess(m ~ a)
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

## Hirearchal clustering
logCPM <- cpm(dge, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se)
outcome <- paste(substr(colnames(se), 9, 12), as.character(se$type), sep="-")
names(outcome) <- colnames(se)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))


plotMDS(dge, labels=outcome, col=batch)
legend("bottomleft", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)), inset=0.05)

## Differential expression
library(sva)
mod <- model.matrix(~ se$type, colData(se))
mod0 <- model.matrix(~ 1, colData(se))
pv <- f.pvalue(assays(se)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.01)
hist(pv,main="", las=1)
# surrogate variables
sv <- sva(assays(se)$logCPM, mod, mod0)

modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(se)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="fdr") < 0.01)
hist(pvsv,main="", las=1)
