library(SummarizedExperiment)
library(geneplotter)
library(edgeR)

se <- readRDS(file = "data/seBRCA.rds")

## Visualizing data
#mcols(se) # Visualize genes
#metadata(se)$objectCreationDate
#mcols(colData(se), use.names = TRUE) #Phenotypic information
#rowRanges(se) #Feature data
names(se) <- rowRanges(se)$symbol

## Subsetting, only samples that we have data for normal and tumor for same patient
codesnorm <- colData(se)[colData(se)$type == "normal",]$bcr_patient_barcode
codestum <- colData(se)[colData(se)$type == "tumor",]$bcr_patient_barcode
commoncodes <- codesnorm[codesnorm %in% codestum]
codefilter <- colData(se)$bcr_patient_barcode %in% commoncodes
se <- se[,!codefilter]
#filtse <- se[,codefilter]

##
dge <- DGEList(counts = assays(se)$counts, genes = mcols(se))
names(dge) #Fields in DGEList
logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)
assays(se)$logCPM <- logCPM# Calculate CPM and store it within se
# assays(se)$logCPM[1:5, 1:5]

## Library sizes
# Library size is the sum of all counts for a sample
sampledepth <- round(dge$sample$lib.size / 1e6, digits = 1)
names(sampledepth) <- substr(colnames(se), 6, 12)
types <- se$type
dge$samples$type <- types
ordering <- order(sampledepth)
sorteddepth <- sampledepth[ordering]
ordtype <- types[ordering]
filtdepth<-sorteddepth[sorteddepth>40]
filttype<-ordtype[sorteddepth>40]
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
# head(dge$samples$norm.factors)
assays(se)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)

### MA plot (ommited for now)
##par(mfrow = c(9, 3), mar = c(4, 5, 3, 1))
#setmp <- se[, se$type == "normal"]
#dgetmp <- dge[, se$type == "normal"]
#for (i in 1:ncol(setmp)) {
#  # For each sample, draw a scatter plot comparing the average logCPM values for 
#  # each gene versus the difference of the logCPM values for each genes for that
#  # sample and the average values
#  A <- rowMeans(assays(setmp)$logCPM)
#  M <- assays(setmp)$logCPM[, i] - A
#  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
#  smoothScatter(A, M, main = samplename, las = 1)
#  abline(h = 0, col = "blue", lwd = 2)
#  lo <- lowess(M ~ A)
#  lines(lo$x, lo$y, col = "red", lwd = 2)
#}
#
##par(mfrow=c(22, 3), mar=c(4, 5, 3, 1))
#
#setmp <- se[, se$type == "tumor"]
#dgetmp <- dge[, se$type == "tumor"]
#for (i in 1:ncol(setmp)) {
#  a <- rowmeans(assays(setmp)$logcpm)
#  m <- assays(setmp)$logcpm[, i] - a
#  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
#  smoothscatter(a, m, main = samplename, las = 1)
#  abline(h = 0, col = "blue", lwd = 2)
#  lo <- lowess(m ~ a)
#  lines(lo$x, lo$y, col = "red", lwd = 2)
#}
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

batches = c("tss", "plate", "portionanalyte", "samplevial")
## Hirearchal clustering

batch_analysis <- function(batch_e){
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  d <- as.dist(1-cor(logCPM, method="spearman"))
  sampleClustering <- hclust(d)
  batch <- as.integer(factor(get(batch_e)))
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
  png(paste(batch_e,"_dend.png",sep=""),height=7, width=14, units="in", res=100)
  par(cex = 0.6)
  plot(sampleDendrogram, main = "Hierarchical clustering of samples")
  legend("topright", paste(levels(factor(get(batch_e)))), fill=sort(unique(batch)))
  dev.off() 
  png(paste(batch_e,"_mds.png",sep=""),height=7, width=14, units="in", res=100)
  par(cex = 0.8)
  plotMDS(dge, labels=outcome, col=batch)
  legend("bottomleft", paste(levels(factor(get(batch_e)))),
         fill=sort(unique(batch)), inset=0.05, horiz=TRUE, cex=0.8)
  dev.off()
}

logCPM <- cpm(dge, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se)
outcome <- substr(se$type,1,2)
names(outcome) <- colnames(se)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
dev.off()
png("dendogram.png", height=10, width=10, units="in", res=100)
par(cex = 0.6)
plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste(levels(factor(tss))), fill=sort(unique(batch)))
dev.off()
png("mds.png", height=10, width=10, units="in", res=100)
par(cex = 0.8)
plotMDS(dge, labels=outcome, col=batch)
legend("bottomleft", paste(levels(factor(tss))),
       fill=sort(unique(batch)), inset=0.05, horiz=TRUE, cex=0.8)
dev.off()
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

### Linear regression DE
adjust_model <- function(model, genesmd, se){
  fit <- lmFit(assays(se)$logCPM, design)
  fit <- eBayes(fit)
  FDRcutoff <- 0.1
  res <- decideTests(fit, p.value = FDRcutoff)
  summary(res)
  # Adding choromosome information
  genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se))), symbol = names(rowRanges(se)), stringsAsFactors = FALSE)
  fit$genes <- genesmd
  tt <- topTable(fit, coef = 2, n = Inf)
  return(c(fit, tt))
}
design <- model.matrix(~ type, data = colData(se))
# First fit
fit <- lmFit(assays(se)$logCPM, design)
fit <- eBayes(fit)
FDRcutoff <- 0.1
res <- decideTests(fit, p.value = FDRcutoff)
summary(res)
# Adding choromosome information
genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se))), symbol = names(rowRanges(se)), stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt)
# Show chromosome distribution of DE genes
sort(table(tt$chr[tt$adj.P.Val < FDRcutoff]), decreasing = TRUE)
# 
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(tt$P.Value, xlab = "Raw P-values", main = "", las = 1)
qqt(fit$t[, 2], df = fit$df.prior + fit$df.residual, main = "", pch = ".", cex = 3)
abline(0, 1, lwd = 2)
# Second fit
v <- voom(dge, design, plot=TRUE)
mod0 <- model.matrix(~concentration, colData(se))
sv <- sva(v$E, mod = design, mod0 = mod0)
design <- cbind(design, sv$sv)
colnames(design) <- c(colnames(design)[1:2], paste0("SV", 1:sv$n))
fit2 <- lmFit(v, design)
fit2 <- eBayes(fit2)
res2 <- decideTests(fit2, p.value = FDRcutoff)
summary(res2)
fit2$genes <- genesmd
tt2 <- topTable(fit2, coef = 2, n = Inf)
head(tt2)
DEgenes2 <- rownames(tt2)[tt2$adj.P.Val < FDRcutoff]
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(tt2$P.Value, xlab = "Raw P-values", main = "", las = 1)
qqt(fit2$t[, 2], df = fit2$df.prior + fit2$df.residual, main = "", pch = ".", cex = 3)
abline(0, 1, lwd = 2)
dev.off()
volcanoplot(fit2, coef = 2, highlight = 7, fit$genes$symbol, main = "Model 1", las = 1)
top7 <- order(fit2$lods[, 2], decreasing = TRUE)[1:7]
limma::plotMA(fit2, coef = 2, status = rownames(fit2$lods) %in% DEgenes2, legend = FALSE,
       main = "Model 1", hl.pch = 46, hl.cex = 4, bg.pch = 46, bg.cex = 3, las = 1)
text(fit2$Amean[top7], fit2$coef[top7, 2], fit2$genes$symbol[top7], cex = 0.5, pos = 4)
#Identify potential source of bias
n <- 100
phenot <- as.data.frame(colData(se)[1, 1:n])
rownames(phenot) <- NULL
phenot
# Ages
ages <- colData(se)[,c("age_at_diagnosis")]
ages <- as.numeric(levels(ages))[ages]
hist(ages)
# Tumor type
types <- colData(se)[,"histological_type"]
types <- na.omit(types)
types <- as.character(types)
table(types)
mat <- as.data.frame(table(types))
rownames(mat)<- c("IDC", "ILC", "MdC", "MH", "MuC","Other")
barplot(mat$Freq, names.arg = rownames(mat))

## Model with known covariates
model_design <- model.matrix( ~ type + age_at_diagnosis, data = colData(se))
fit4 <- lmFit(assays(se)$logCPM, model_design)
fit4 <- eBayes(fit)
FDRcutoff <- 0.1
res4 <- decideTests(fit, p.value = FDRcutoff)
summary(res4)
# Adding choromosome information
fit4$genes <- genesmd
tt4 <- topTable(fit4, coef = 2, n = Inf)

### Functional annotation
library(org.Hs.eg.db)
