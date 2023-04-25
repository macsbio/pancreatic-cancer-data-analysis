# data exploration and analysis of GSE19650
# author: Isabel Wassink
# date: Nov 2022

# load libraries
library(GEOquery)
library(limma)
library(umap)
library(data.table)
library(magrittr)
library(dplyr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(rWikiPathways)
library(clusterProfiler)
library(ggvenn)
library(gplots)
library(RColorBrewer)
library(WGCNA)
library(RCy3)
library(rstudioapi)
library(readr)
library(enrichplot)
library(ggnewscale)
library(plotly)

# load series and platform data from GEO
gset <- getGEO("GSE19650", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "3333333000000111111222"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("IPMA","IPMC","Cancer","Normal"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# filter probe duplicates
o <- order(fit$Amean, decreasing = TRUE)
dup <- duplicated(fit$genes$Gene.ID [o])
fit.unique <- fit[o,][!dup,]

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[4],sep=""), paste(groups[2],"-",groups[4],sep=""), paste(groups[3],"-",groups[4],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit.unique, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))

################################################################
# visualize and quality control of the data

# build histogram of P-values for all genes. Normal test
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

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
ct <- 1        # choose contrast of interest
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

##############
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE19650", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE19650", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=9", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE19650")

################################################################
## IPMA vs Normal

# set up contrasts of interest and recalculate model coefficients
cts.IPMA <- c(paste(groups[1],"-",groups[4],sep=""))
cont.matrix.IPMA <- makeContrasts(contrasts=cts.IPMA, levels=design)
fit2.IPMA <- contrasts.fit(fit.unique, cont.matrix.IPMA)

# compute statistics
fit2.IPMA <- eBayes(fit2.IPMA, 0.01)
tT.IPMA <- topTable(fit2.IPMA, adjust="fdr", sort.by="B", number=Inf)

tT.IPMA <- subset(tT.IPMA, select=c("ID","Gene.title", "Gene.symbol", "Gene.ID",
                                    "GO.Function.ID", "GO.Process.ID", "GO.Component.ID",
                                    "logFC", "P.Value", "adj.P.Val"))

dT.IPMA <- decideTests(fit2.IPMA, adjust.method="fdr", p.value=0.05)
summary(dT.IPMA)

# volcano plot 
colnames(tT.IPMA) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2.IPMA, coef=ct, main=colnames(fit2.IPMA)[ct], pch=20,
            highlight=length(which(dT.IPMA[,ct]!=0)), names=rep('+', nrow(fit2.IPMA)))

EnhancedVolcano(tT.IPMA, title="IPMA-Normal", lab = tT.IPMA$Gene.symbol, labSize = 3, x = "logFC", y = "adj.P.Val", pCutoff = 0.05, FCcutoff = 1)

################################################################
## IPMC vs Normal
cts.IPMC <- c(paste(groups[2],"-",groups[4],sep=""))
cont.matrix.IPMC <- makeContrasts(contrasts=cts.IPMC, levels=design)
fit2.IPMC <- contrasts.fit(fit.unique, cont.matrix.IPMC)

fit2.IPMC <- eBayes(fit2.IPMC, 0.01)
tT.IPMC <- topTable(fit2.IPMC, adjust="fdr", sort.by="B", number=Inf)

tT.IPMC <- subset(tT.IPMC, select=c("ID","Gene.title", "Gene.symbol", "Gene.ID",
                                    "GO.Function.ID", "GO.Process.ID", "GO.Component.ID",
                                    "logFC", "P.Value", "adj.P.Val"))

dT.IPMC <- decideTests(fit2.IPMC, adjust.method="fdr", p.value=0.05)
summary(dT.IPMC)

# volcano plot 
colnames(tT.IPMC) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2.IPMC, coef=ct, main=colnames(fit2.IPMC)[ct], pch=20,
            highlight=length(which(dT.IPMC[,ct]!=0)), names=rep('+', nrow(fit2.IPMC)))

EnhancedVolcano(tT.IPMC, title="IPMC-Normal", lab = tT.IPMA$Gene.symbol, labSize = 3, x = "logFC", y = "adj.P.Val", pCutoff = 0.05, FCcutoff = 1)

################################################################
## Cancer vs Normal
cts.Can <- c(paste(groups[3],"-",groups[4],sep=""))
cont.matrix.Can <- makeContrasts(contrasts=cts.Can, levels=design)
fit2.Can <- contrasts.fit(fit.unique, cont.matrix.Can)

fit2.Can <- eBayes(fit2.Can, 0.01)
tT.Can <- topTable(fit2.Can, adjust="fdr", sort.by="B", number=Inf)

tT.Can <- subset(tT.Can, select=c("ID","Gene.title", "Gene.symbol", "Gene.ID",
                                  "GO.Function.ID", "GO.Process.ID", "GO.Component.ID",
                                  "logFC", "P.Value", "adj.P.Val"))

dT.Can <- decideTests(fit2.Can, adjust.method="fdr", p.value=0.05)
summary(dT.Can)

# volcano plot 
colnames(tT.Can) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2.Can, coef=ct, main=colnames(fit2.Can)[ct], pch=20,
            highlight=length(which(dT.Can[,ct]!=0)), names=rep('+', nrow(fit2.Can)))

EnhancedVolcano(tT.Can, title="Can-Normal", lab = tT.IPMA$Gene.symbol, labSize = 3, x = "logFC", y = "adj.P.Val", pCutoff = 0.05, FCcutoff = 1)

################################################################
# change column names to allow merging 
colnames(tT.IPMA)[8] <- "logFC.IPMA"
colnames(tT.IPMA)[9] <- "P.Value.IPMA"
colnames(tT.IPMA)[10] <- "adj.P.val.IPMA"

colnames(tT.IPMC)[8] <- "logFC.IPMC"
colnames(tT.IPMC)[9] <- "P.Value.IPMC"
colnames(tT.IPMC)[10] <- "adj.P.val.IPMC"

colnames(tT.Can)[8] <- "logFC.Can"
colnames(tT.Can)[9] <- "P.Value.Can"
colnames(tT.Can)[10] <- "adj.P.val.Can"

# merge data frames
rownames(tT.IPMA) <- NULL
rownames(tT.IPMC) <- NULL
rownames(tT.Can) <- NULL

data.all <- merge(tT.IPMA, tT.IPMC[,c(4,8,9,10)], by = 'Gene.ID', all = TRUE)
data.all <- merge(data.all, tT.Can[,c(4,8,9,10)],  by = 'Gene.ID', all = TRUE)

# filter out genes with multiple gene IDs
data.all <- data.all %>% filter(!grepl('///', Gene.ID))

# filter out genes without gene ID
data.all <- data.all[!data.all$Gene.ID == "",]

# remove irrelevant columns
data.all = subset(data.all, select = c(-ID, -Gene.title, -GO.Function.ID, -GO.Component.ID))

################################################################
# venn diagrams
# all differentially expressed genes
IPMA.deg <- unique(data.all[data.all$adj.P.val.IPMA < 0.05 & abs(data.all$logFC.IPMA) > 1, c(1,2)])
IPMC.deg <- unique(data.all[data.all$adj.P.val.IPMC < 0.05 & abs(data.all$logFC.IPMC) > 1, c(1,2)])
Can.deg <- unique(data.all[data.all$adj.P.val.Can < 0.05 & abs(data.all$logFC.Can) > 1, c(1,2)])

deg <- list(IPMA = IPMA.deg$Gene.ID,
            IPMC = IPMC.deg$Gene.ID,
            Can = Can.deg$Gene.ID)

ggvenn(deg, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), stroke_size = 0.5,
       set_name_size = 4, show_percentage = FALSE)

# up-regulated genes
IPMA.up <- unique(data.all[data.all$logFC.IPMA > 1 & data.all$adj.P.val.IPMA < 0.05, c(1,2)])
IPMC.up <- unique(data.all[data.all$logFC.IPMC > 1 & data.all$adj.P.val.IPMC < 0.05, c(1,2)])
Can.up <- unique(data.all[data.all$logFC.Can > 1 & data.all$adj.P.val.Can < 0.05, c(1,2)])

up <- list(IPMA = IPMA.up$Gene.ID,
           IPMC = IPMC.up$Gene.ID,
           Can = Can.up$Gene.ID)

ggvenn(up, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), stroke_size = 0.5,
       set_name_size = 4, show_percentage = FALSE)

# down-regulated genes
IPMA.down <- unique(data.all[data.all$logFC.IPMA < 1 & data.all$adj.P.val.IPMA < 0.05, c(1,2)])
IPMC.down <- unique(data.all[data.all$logFC.IPMC < 1 & data.all$adj.P.val.IPMC < 0.05, c(1,2)])
Can.down <- unique(data.all[data.all$logFC.Can < 1 & data.all$adj.P.val.Can < 0.05, c(1,2)])

down <- list(IPMA = IPMA.down$Gene.ID,
             IPMC = IPMC.down$Gene.ID,
             Can = Can.down$Gene.ID)

ggvenn(down, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), stroke_size = 0.5,
       set_name_size = 4, show_percentage = FALSE)

## Sankey plot
# calculate values for Sankey input 
IPMA.notDE <- as.data.frame(outersect(data.all[,2], IPMA.deg[,2]))
IPMC.notDE <- as.data.frame(outersect(data.all[,2], IPMC.deg[,2]))
Can.notDE <- as.data.frame(outersect(data.all[,2], Can.deg[,2]))

ipma_up_ipmc_up <- as.data.frame(intersect(IPMA.up[,2], IPMC.up[,2]))
ipma_up_ipmc_down <- as.data.frame(intersect(IPMA.up[,2], IPMC.down[,2]))
ipma_up_ipmc_not <- as.data.frame(intersect(IPMA.up[,2], IPMC.notDE[,1]))

ipma_down_ipmc_up <- as.data.frame(intersect(IPMA.down[,2], IPMC.up[,2]))
ipma_down_ipmc_down <- as.data.frame(intersect(IPMA.down[,2], IPMC.down[,2]))
ipma_down_ipmc_not <- as.data.frame(intersect(IPMA.down[,2], IPMC.notDE[,1]))

ipma_not_ipmc_up <- as.data.frame(intersect(IPMA.notDE[,1], IPMC.up[,2]))
ipma_not_ipmc_down <- as.data.frame(intersect(IPMA.notDE[,1], IPMC.down[,2]))
ipma_not_ipmc_not <- as.data.frame(intersect(IPMA.notDE[,1], IPMC.notDE[,1]))

ipmc_up_can_up <- as.data.frame(intersect(IPMC.up[,2], Can.up[,2]))
ipmc_up_can_down <- as.data.frame(intersect(IPMC.up[,2], Can.down[,2]))
ipmc_up_can_not <- as.data.frame(intersect(IPMC.up[,2], Can.notDE[,1]))

ipmc_down_can_up <- as.data.frame(intersect(IPMC.down[,2], Can.up[,2]))
ipmc_down_can_down <- as.data.frame(intersect(IPMC.down[,2], Can.down[,2]))
ipmc_down_can_not <- as.data.frame(intersect(IPMC.down[,2], Can.notDE[,1]))

ipmc_not_can_up <- as.data.frame(intersect(IPMC.notDE[,1], Can.up[,2]))
ipmc_not_can_down <- as.data.frame(intersect(IPMC.notDE[,1], Can.down[,2]))
ipmc_not_can_not <- as.data.frame(intersect(IPMC.notDE[,1], Can.notDE[,1]))

# plot 
fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  node = list(
    label = c("IPMA_up", "IPMA_down","IPMA_not", "IPMC_up", "IPMC_down", "IPMC_not", "Cancer_up", "Cancer_down","Cancer_not"),
    color = c("red","blue", "white","red","blue", "white","red","blue", "white"),
    pad = 15,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),
    target = c(3,4,5,3,4,5,3,4,5,6,7,8,6,7,8,6,7,8),
    value = c(744,82,394,73,1955,1323,538,1547,15822,789,23,493,124,1734,1533,549,1340,15373)
  )
)

fig <- fig %>% layout(
  title = "Gene expression - disease progression",
  font = list(
    size = 10
  )
)

fig

################################################################
# filtering of immune genes
immunegenes <- data.all[data.all$GO.Process.ID %like% "GO:0002376" | 
                          data.all$GO.Process.ID %like% "GO:0006955" |
                          data.all$GO.Process.ID %like% "GO:0006954",]

# heatmap of immune genes 
rnames <- immunegenes[,1]                                   # assign labels in column 1 to "rnames"
mat_data <- data.matrix(immunegenes[, c(4, 7, 10)])         # transform logFC columns into a matrix
rownames(immunegenes) <- rnames  

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

colors <- seq(from=-3, to=3, length.out=300)
map <- heatmap.2(mat_data,
                 col = my_palette,
                 breaks = colors,
                 trace = "none", 
                 density.info = "none",
                 Colv = FALSE,
                 cexRow = 0.4,
                 cexCol = 0.8,
                 keysize = 0.5,
                 key.par = list(cex=0.4)
)

# use if new plot doesn't work after error:
# dev.off()

################################################################
# data visualization on pathways
data.values <- c(-1,0,1)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))

# WP5285 - Immune infiltration in pancreatic cancer
RCy3::commandsRun('wikipathways import-as-pathway id=WP5285') 
toggleGraphicsDetails()
RCy3::mapTableColumn("Ensembl", "Homo sapiens", "Ensembl", "Entrez Gene")
loadTableData(data.all, data.key.column = "Gene.ID", table.key.column = "Entrez Gene")

RCy3::setNodeLabelMapping("name", "pathway_viz")

RCy3::copyVisualStyle("WikiPathways", "pathway_viz")

RCy3::setNodeCustomHeatMapChart(c("logFC.IPMA", "logFC.IPMC", "logFC.Can"), slot = 2,
                                style.name = "pathway_viz", colors = node.colors)

setNodeLabelColorDefault('#000000', "pathway_viz")
setNodeBorderColorDefault('#000000', "pathway_viz")
setNodeColorDefault("#FFFFFF", "pathway_viz")

RCy3::setVisualStyle("pathway_viz")

# WP5078 - T cell modulation in pancreatic cancer
RCy3::commandsRun('wikipathways import-as-pathway id=WP5078') 
toggleGraphicsDetails()
RCy3::mapTableColumn("Ensembl", "Homo sapiens", "Ensembl", "Entrez Gene")
loadTableData(data.all, data.key.column = "Gene.ID", table.key.column = "Entrez Gene")

RCy3::setNodeLabelMapping("name", "pathway_viz")

RCy3::copyVisualStyle("WikiPathways", "pathway_viz")

RCy3::setNodeCustomHeatMapChart(c("logFC.IPMA", "logFC.IPMC", "logFC.Can"), slot = 2,
                                style.name = "pathway_viz", colors = node.colors)

setNodeLabelColorDefault('#000000', "pathway_viz")
setNodeBorderColorDefault('#000000', "pathway_viz")
setNodeColorDefault("#FFFFFF", "pathway_viz")

RCy3::setVisualStyle("pathway_viz")

################################################################
# overrepresentation analysis
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens",format = "gmt",date = "20221210")

wp2gene <- readPathwayGMT(gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

bkgd.genes <- unique(data.all[,c(1,2)])

## IPMA DEG
# WP
ewp.IPMA <- clusterProfiler::enricher(
  IPMA.deg$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.IPMA.res <- as.data.frame(ewp.IPMA) 

length(ewp.IPMA@universe)         # number of genes measured in pathways
length(IPMA.deg$Gene.ID[IPMA.deg$Gene.ID %in% unique(wp2gene$gene)]) # number of DEG in pathways

num.pathways.IPMA <- dim(ewp.IPMA.res)[1]

# GO
ego.IPMA <- clusterProfiler::enrichGO(
  IPMA.deg$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.IPMA.res <- as.data.frame(ego.IPMA) 

## IPMC DEG
# WP
ewp.IPMC <- clusterProfiler::enricher(
  IPMC.deg$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.IPMC.res <- as.data.frame(ewp.IPMC) 

length(ewp.IPMC@universe)        
length(IPMC.deg$Gene.ID[IPMC.deg$Gene.ID %in% unique(wp2gene$gene)]) 

num.pathways.IPMC <- dim(ewp.IPMC.res)[1]

# GO
ego.IPMC <- clusterProfiler::enrichGO(
  IPMC.deg$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.IPMC.res <- as.data.frame(ego.IPMC) 

## Can DEG
# WP
ewp.Can <- clusterProfiler::enricher(
  Can.deg$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.Can.res <- as.data.frame(ewp.Can) 

length(ewp.Can@universe)        
length(Can.deg$Gene.ID[Can.deg$Gene.ID %in% unique(wp2gene$gene)]) 

num.pathways.Can <- dim(ewp.Can.res)[1]

# GO
ego.Can <- clusterProfiler::enrichGO(
  Can.deg$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.Can.res <- as.data.frame(ego.Can) 

## venn diagram of DEG pathways
pathways.deg <- list(IPMA = ewp.IPMA.res$ID, 
                     IPMC = ewp.IPMC.res$ID,
                     Can = ewp.Can.res$ID)

ggvenn(pathways.deg, show_percentage = FALSE)

## venn diagram of DEG GO
GO.deg <- list(IPMA = ego.IPMA.res$ID, 
               IPMC = ego.IPMC.res$ID,
               Can = ego.Can.res$ID)

ggvenn(GO.deg, show_percentage = FALSE)

## visualization of DEG GO enrichment 
max.overlaps=Inf

edo.IPMA <- pairwise_termsim(ego.IPMA)
p.IPMA <- emapplot(edo.IPMA)
p.IPMA

edo.IPMC <- pairwise_termsim(ego.IPMC)
p.IPMC <- emapplot(edo.IPMC)
p.IPMC

edo.Can <- pairwise_termsim(ego.Can)
p.Can <- emapplot(edo.Can)
p.Can

cowplot::plot_grid(p.IPMA, p.IPMC, p.Can, ncol=2, labels=LETTERS[1:3])

## table of DEG for online GO analysis
write.table(IPMA.deg[,2], file="deg.ipma.tsv", quote=F, sep="\t", row.names = F, col.names = F)

write.table(IPMC.deg[,2], file="deg.ipmc.tsv", quote=F, sep="\t", row.names = F, col.names = F)

write.table(Can.deg[,2], file="deg.can.tsv", quote=F, sep="\t", row.names = F, col.names = F)

#################################
## IPMA up
# WP
ewp.IPMA.up <- clusterProfiler::enricher(
  IPMA.up$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.IPMA.up.res <- as.data.frame(ewp.IPMA.up) 

length(ewp.IPMA.up@universe)         
length(IPMA.up$Gene.ID[IPMA.up$Gene.ID %in% unique(wp2gene$gene)])

num.pathways.IPMA.up <- dim(ewp.IPMA.up.res)[1]

# GO
ego.IPMA.up <- clusterProfiler::enrichGO(
  IPMA.up$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.IPMA.up.res <- as.data.frame(ego.IPMA.up) 

## IPMC up
# WP
ewp.IPMC.up <- clusterProfiler::enricher(
  IPMC.up$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.IPMC.up.res <- as.data.frame(ewp.IPMC.up) 

length(ewp.IPMC.up@universe)         
length(IPMC.up$Gene.ID[IPMC.up$Gene.ID %in% unique(wp2gene$gene)])

num.pathways.IPMC.up <- dim(ewp.IPMC.up.res)[1]

# GO
ego.IPMC.up <- clusterProfiler::enrichGO(
  IPMC.up$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.IPMC.up.res <- as.data.frame(ego.IPMC.up) 

## Can up
# WP
ewp.Can.up <- clusterProfiler::enricher(
  Can.up$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.Can.up.res <- as.data.frame(ewp.Can.up) 

length(ewp.Can.up@universe)         
length(Can.up$Gene.ID[Can.up$Gene.ID %in% unique(wp2gene$gene)])

num.pathways.Can.up <- dim(ewp.Can.up.res)[1]

# GO
ego.Can.up <- clusterProfiler::enrichGO(
  Can.up$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.Can.up.res <- as.data.frame(ego.Can.up) 

## venn diagram of up pathways
pathways.up <- list(IPMA = ewp.IPMA.up.res$ID, 
                    IPMC = ewp.IPMC.up.res$ID,
                    Can = ewp.Can.up.res$ID)

ggvenn(pathways.up, show_percentage = FALSE)

# venn diagram of up GO
GO.up <- list(IPMA = ego.IPMA.up.res$ID, 
              IPMC = ego.IPMC.up.res$ID,
              Can = ego.Can.up.res$ID)

ggvenn(GO.up, show_percentage = FALSE)

#################################
## IPMA down
# WP
ewp.IPMA.down <- clusterProfiler::enricher(
  IPMA.down$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.IPMA.down.res <- as.data.frame(ewp.IPMA.down) 

length(ewp.IPMA.down@universe)         
length(IPMA.down$Gene.ID[IPMA.down$Gene.ID %in% unique(wp2gene$gene)])

num.pathways.IPMA.down <- dim(ewp.IPMA.down.res)[1]

# GO
ego.IPMA.down <- clusterProfiler::enrichGO(
  IPMA.down$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.IPMA.down.res <- as.data.frame(ego.IPMA.down) 

## IPMC down
# WP
ewp.IPMC.down <- clusterProfiler::enricher(
  IPMC.down$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.IPMC.down.res <- as.data.frame(ewp.IPMC.down) 

length(ewp.IPMC.down@universe)         
length(IPMC.down$Gene.ID[IPMC.down$Gene.ID %in% unique(wp2gene$gene)])

num.pathways.IPMC.down <- dim(ewp.IPMC.down.res)[1]

# GO
ego.IPMC.down <- clusterProfiler::enrichGO(
  IPMC.down$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.IPMC.down.res <- as.data.frame(ego.IPMC.down) 

## Can down
# WP
ewp.Can.down <- clusterProfiler::enricher(
  Can.down$Gene.ID,
  universe = bkgd.genes$Gene.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

ewp.Can.down.res <- as.data.frame(ewp.Can.down) 

length(ewp.Can.down@universe)         
length(Can.down$Gene.ID[Can.down$Gene.ID %in% unique(wp2gene$gene)])

num.pathways.Can.down <- dim(ewp.Can.down.res)[1]

# GO
ego.Can.down <- clusterProfiler::enrichGO(
  Can.down$Gene.ID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

ego.Can.down.res <- as.data.frame(ego.Can.down) 

## venn diagram of down pathways
pathways.down <- list(IPMA = ewp.IPMA.down.res$ID, 
                      IPMC = ewp.IPMC.down.res$ID,
                      Can = ewp.Can.down.res$ID)

ggvenn(pathways.down, show_percentage = FALSE)

## venn diagram of down GO
GO.down <- list(IPMA = ego.IPMA.down.res$ID, 
                IPMC = ego.IPMC.down.res$ID,
                Can = ego.Can.down.res$ID)

ggvenn(GO.down, show_percentage = FALSE)

#################################
# WP pathway overlap visualization
data.values.wp <-c(-1,0,1) 
node.colors.wp <- c(rev(brewer.pal(length(data.values), "RdBu")))

# choose genes for analysis (IPMA/IPMC/CAN, deg/up/down)
x <- ewp.Can.res
y <- Can.deg
title <- "Cancer-DEG"

# visualization 
pwy <- unique(x[,c(1,2)])
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'
edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")
genes <- unique(y)
colnames(genes) <- c("id","label")
genes <- transform(genes, id = as.character(id))
genes$type <- 'gene'
edges <- unique(edges[edges$target %in% genes$id,])
genes <- genes[genes$id %in% edges$target,]
nodes <- dplyr::bind_rows(genes, pwy)
rownames(nodes) <- NULL
createNetworkFromDataFrames(nodes=nodes,edges=edges,title= title, collection="PathwayGeneCrosstalk")
loadTableData(data.all, data.key.column = "Gene.ID", table.key.column = "id")

# set the correct visual style (IPMA/IPMC/Can)

# set visual style - IPMA
RCy3::copyVisualStyle("default","IPMA.vis")
RCy3::setNodeLabelMapping("label", style.name="IPMA.vis")
RCy3::lockNodeDimensions(TRUE, style.name="IPMA.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="IPMA.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "IPMA.vis")
setNodeColorMapping("logFC.IPMA", data.values.wp, node.colors.wp, default.color = "#99FF99", style.name = "IPMA.vis")
RCy3::setVisualStyle("IPMA.vis")
RCy3::toggleGraphicsDetails()

# set visual style - IPMC
RCy3::copyVisualStyle("default","IPMC.vis")
RCy3::setNodeLabelMapping("label", style.name="IPMC.vis")
RCy3::lockNodeDimensions(TRUE, style.name="IPMC.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="IPMC.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "IPMC.vis")
setNodeColorMapping("logFC.IPMC", data.values.wp, node.colors.wp, default.color = "#99FF99", style.name = "IPMC.vis")
RCy3::setVisualStyle("IPMC.vis")
RCy3::toggleGraphicsDetails()

# set visual style - Can
RCy3::copyVisualStyle("default","Can.vis")
RCy3::setNodeLabelMapping("label", style.name="Can.vis")
RCy3::lockNodeDimensions(TRUE, style.name="Can.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="Can.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "Can.vis")
setNodeColorMapping("logFC.Can", data.values.wp, node.colors.wp, default.color = "#99FF99", style.name = "Can.vis")
RCy3::setVisualStyle("Can.vis")
RCy3::toggleGraphicsDetails()