# data exploration and analysis of GSE15471 (matching pairs)
# author: Isabel Wassink
# date: March 2023

# load libraries
library(GEOquery)
library(limma)
library(umap)
library(ggfortify)
library(dplyr)
library(EnhancedVolcano)
library(rWikiPathways)
library(enrichplot)
library(RColorBrewer)
library(RCy3)

# load series and platform data from GEO
gset <- getGEO("GSE15471", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("X1X1X1111111111111111111111111111111111X0X0X000000",
               "0000000000000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("tumor","normal"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

### EXPERIMENTAL DESCRIPTION ###
# read experimental annotations (separate file in repository)
desc <- read.delim("Description_file_GSE15471.txt",as.is=TRUE)
desc$FactorValue <- as.factor(desc$FactorValue)
desc$Patient <- factor(desc$Patient, levels = c("30162", "40728", "41027", "30057", "30068", 
                                                "30277", "30308", "30364", "30582", "30617",
                                                "40645", "40656", "40726", "40730", "40741", 
                                                "40836", "40843", "40875", "40892", "40899",
                                                "40975", "40977", "51084", "51091", "51176", 
                                                "51292", "51294", "51308", "51315", "51572", 
                                                "51628", "51677", "51681", "51721", "51722", 
                                                "51783"))

# replace FactorValue by Group as column names
colnames(desc)[2] <- "Group"

# create short sample names as row names
rownames(desc) <- desc$Name

# replace colnames in gset by SampleNames from desc
# check whether they are in the same order before doing so
sum(gsub("X","",colnames(gset))!=desc$SourceName)==0 #FALSE
# replace long column names in gset by short ones
colnames(gset) <- rownames(desc)
colnames(exprs(gset))

### QC PLOTS ###
# PCA plot 
dat <- as.matrix(gset)
dat.turned <- t(dat)

dat.pca <- prcomp(dat.turned, center = TRUE, scale. = TRUE)
summary(dat.pca)

dat.pca.plot <- autoplot(dat.pca)
dat.pca.plot

# dendrogram
dd <- dist(scale(dat.turned), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

plot(hc, hang = -1, cex = 0.6)

### STATISTICAL MODELING
# set up design matrix
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# create a factor variable for the replicate_numbers
replicate_numbers2 <- as.factor(desc$Patient)

# estimate the extra correlation between paired measurements
corfit <- duplicateCorrelation(gset,design,block=replicate_numbers2)
corfit$consensus

# fit linear model
fit <- lmFit(gset,design,block=replicate_numbers2,correlation=corfit$consensus)
fit <- eBayes(fit)

# filter probe duplicates
o <- order(fit$Amean, decreasing = TRUE)
dup <- duplicated(fit$genes$Gene.ID [o])
fit.unique <- fit[o,][!dup,]

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit.unique, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.ID"))

# filter out genes with multiple gene IDs
tT <- tT %>% filter(!grepl('///', Gene.ID))

# filter out genes without gene ID
tT <- tT[!tT$Gene.ID == "",]

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

EnhancedVolcano(tT, title="Tumor-Normal", lab = rownames(tT), labSize = 3, 
                x = "logFC", y = "adj.P.Val", pCutoff = 0.05, FCcutoff = 1)

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE15471", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE15471", "/", annotation(gset), "value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE15471")

### ANALYSIS ### 
deg <- unique(tT[tT$adj.P.Val < 0.05 & abs(tT$logFC) > 1, c(1:8)])
down <- tT[tT$logFC < -1 & tT$adj.P.Val < 0.05, c(1:8)]
up <- tT[tT$logFC > 1 & tT$adj.P.Val < 0.05, c(1:8)]

## enrichment analysis
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens",format = "gmt",date = "20230210")

wp2gene <- readPathwayGMT(gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

bkgd.genes <- unique(tT[,c(1-9)])

# WP - DEG
ewp.DEG <- clusterProfiler::enricher(
  deg$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.DEG.res <- as.data.frame(ewp.DEG) 

# GO - DEG
ego.DEG <- clusterProfiler::enrichGO(
  deg$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.DEG.res <- as.data.frame(ego.DEG)

# WP - down
ewp.down <- clusterProfiler::enricher(
  down$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.down.res <- as.data.frame(ewp.down) 

# GO - down
ego.down <- clusterProfiler::enrichGO(
  down$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.down.res <- as.data.frame(ego.down)

# WP - up
ewp.up <- clusterProfiler::enricher(
  up$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.up.res <- as.data.frame(ewp.up) 

# GO - up
ego.up <- clusterProfiler::enrichGO(
  up$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.up.res <- as.data.frame(ego.up)

### visualization of GO enrichment 
max.overlaps=Inf

edo <- pairwise_termsim(ego.DEG)
p <- emapplot(edo)
p

### visualization of WP enrichment
library(multienrichjam)

# filter out irrelevant results
ewp.DEG.filt <- ewp.DEG.res[!ewp.DEG.res$ID %in% c("WP5115","WP4816", "WP4217", "WP5218", "WP4786", 
                                                   "WP3668","WP2431", "WP2272", "WP2118", "WP4808", 
                                                   "WP2572", "WP2865", "WP4197"), ]
ewp.DEG.filt <- enrichDF2enrichResult(ewp.DEG.filt)
ewp.DEG.filt <- pairwise_termsim(ewp.DEG.filt)

# tree plot
label_func <- function(names) {
  labels <- c("Immune-related pathways", "Complement system pathways", "Cell signaling pathways",
              "Cell signaling/focal adhesion pathways", "Nutrient metabolism pathways", 
              "Energy metabolism pathways")
  return(labels[1:length(names)])
}

tree <- treeplot(ewp.DEG.filt, color = "p.adjust", nCluster = 6, 
                 showCategory <- nrow(ewp.DEG.filt), cex_category = 0.6, 
                 label_params = list(cex = 0.005), label_format = label_func)
print(tree)

png("treeplot.png", width=5100, height=4500, res=300)
print(tree)
dev.off()

### PATHWAY VISUALIZATION ### 
# data visualization on pathways
data.values <- c(-1,0,1)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))

# WP5285 - Immune infiltration in pancreatic cancer
RCy3::commandsRun('wikipathways import-as-pathway id=WP5285') 
toggleGraphicsDetails()
RCy3::mapTableColumn("Ensembl", "Homo sapiens", "Ensembl", "Entrez Gene")
loadTableData(tT, data.key.column = "Gene.symbol", table.key.column = "name")

RCy3::setNodeLabelMapping("name")

RCy3::copyVisualStyle("WikiPathways", "pathway_viz")

RCy3::setNodeCustomHeatMapChart("logFC", slot = 2,
                                style.name = "pathway_viz", colors = node.colors)

setNodeLabelColorDefault('#000000', "pathway_viz")
setNodeBorderColorDefault('#000000', "pathway_viz")
setNodeColorDefault("#FFFFFF", "pathway_viz")

RCy3::setVisualStyle("pathway_viz")

# WP5078 - T cell modulation in pancreatic cancer
RCy3::commandsRun('wikipathways import-as-pathway id=WP5078') 
toggleGraphicsDetails()
RCy3::mapTableColumn("Ensembl", "Homo sapiens", "Ensembl", "Entrez Gene")
loadTableData(tT, data.key.column = "Gene.symbol", table.key.column = "name")

RCy3::setNodeLabelMapping("name")

RCy3::copyVisualStyle("WikiPathways", "pathway_viz")

RCy3::setNodeCustomHeatMapChart("logFC", slot = 2,
                                style.name = "pathway_viz", colors = node.colors)

setNodeLabelColorDefault('#000000', "pathway_viz")
setNodeBorderColorDefault('#000000', "pathway_viz")
setNodeColorDefault("#FFFFFF", "pathway_viz")

RCy3::setVisualStyle("pathway_viz")