#Data input, cleaning and pre-processing

##Loading Data
getwd();
workingDir="~/Desktop/ppanini_data/gene_network/subnetworks";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors=FALSE);

filename="stool/bifid/norm_stool_bifid_table.csv"
inputData = read.csv(filename)
dim(inputData);
names(inputData);
datExpr0 = as.data.frame(t(inputData));
names(datExpr0) = inputData[,1];
datExpr0 = datExpr0[2:39,] ##custom

gsg = goodSamplesGenes(datExpr0, verbose=3);
gsg$allOK

if (!gsg$allOK)
{
	if (sum(!gsg$goodGenes)>0)
		printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse=", ")));
	if (sum(!gsg$goodSamples)>0)
		printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")));
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#Removing samples: SRS023850
sampleTree = hclust(dist(datExpr0), method="average")

pdf(paste(filename, "outliers.pdf", sep=""), width=12, height=9)
par(cex=0.6);
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(h=5e5, col="red")
dev.off()
#Removing samples: SRS017497
clust = cutreeStatic(sampleTree, cutHeight=5e5, minSize=10) #custom
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

gsg = goodSamplesGenes(datExpr, verbose=3);
gsg$allOK
if (!gsg$allOK)
{
	if (sum(!gsg$goodGenes)>0)
		printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse=", ")));
	if (sum(!gsg$goodSamples)>0)
		printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")));
	datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
powers = c(c(1:10), seq(from=12, to=50, by=2))
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose=5)
pdf(paste(filename, "softhresholdfit.pdf", sep=""), width=9, height=5)
par(mfrow=c(1,2));
cex1=0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold(power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main=paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h ylab="Scale Free Topology"labels=powers,cex=cex1,	col="red");
abline(h=0.9, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",main=paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

net = blockwiseModules(data.matrix(datExpr), power=6, TOMType="unsigned", minModuleSize=30, reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE, pamRespectsDendro=FALSE, saveTOMs=TRUE, saveTOMFileBase= paste(filename, "PFTOM", sep=""), verbose=3)
pdf(paste(filename, "_dendroAndcolorsModulecolors.pdf", sep=""), width=12, height=9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
					"Module colors",
					dendroLabels=FALSE,
					hang=0.03,
					addGuide=TRUE,
					guideHang=0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file=paste(filename, "-autonetworkconstruction.RData", sep=""))

TOM = TOMsimilarityFromExpr(data.matrix(datExpr), power =6)

unique(moduleColors)
modules = c("turquoise")
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM, edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse='-'),".txt", sep=""), nodeFile=paste("CytoscapeInput-nodes", paste(modules, collapse="-"),".txt", sep=""), weighted=TRUE, threshold=0.02, nodeNames = modProbes, nodeAttr=moduleColors[inModule]);


nSelect = 159
set.seed(10)
dissTOM = 1-TOM
select = sample(nGenes, size=nSelect);
selectTOM = dissTOM[select, select];
selectTree =hclust(as.dist(selectTOM), method="average")
selectColors = moduleColors[select];
pdf(paste(filename, "_NetworkheatmapPLOTselect.pdf", sep=""), width=9, height=9);
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main="Network heatmap plot, selected genes")
dev.off()

# ##Loading clinical trait data_files
# clinical_filename="data_files/FemaleLiver-Data/ClinicalTraits.csv" #custom
# traitData = read.csv(clinical_filename)
# dim(traitData)
# names(traitData)

# allTraits = traitData[, -c(31, 16)]; #custom
# allTraits = allTraits[, c(2, 11:36)]; #custom
# dim(allTraits)
# names(allTraits)

# femaleSamples = rownames(datExpr);
# traitRows = match(femaleSamples, allTraits$Mice) #custom
# datTraits = allTraits[traitRows, -1]
# rownames(datTraits) = allTraits[traitRows, -1];

# collectGarbage();
# # Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method="average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits, signed=FALSE)
# # Plot the sample dendrogram and the colors underneath.
# pdf(paste(filename, "_dendroAndColorsforTraitsHeatmap.pdf", sep=""), width=12, height=9)
# plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits), main="Sample dendrogram and trait heatmap")
# dev.off()
# save(datExpr, datTraits, file=paste(filename, "-dataInput.RData", sep=""))