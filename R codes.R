#Install & Load the WGCNA package********************************************************
install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("EnrichmentBrowser")
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("bestNormalize")
install.packages("xlsx")
BiocManager::install("GEOquery")
BiocManager::install("genefilter")
BiocManager::install("Organism.dplyr")
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library("Organism.dplyr")
library("GO.db")
library('BiocManager')
library("Biobase")
library("EnrichmentBrowser")
library("matrixStats")
library("Hmisc") 
library("lattice")
library("survival")
library("Formula")
library("ggplot2")
library("splines") 
library("foreach")
library("doParallel") 
library("fastcluster")
library("dynamicTreeCut")
library("survival")
library("xlsx")
library("Biobase")
library("EnrichmentBrowser")
library(readxl)
library(WGCNA)
library(bestNormalize)
library("GEOquery")
library(foreign)
library(genefilter)
library(utils)
options(stringsAsFactors = FALSE)
# Current working directory
getwd();
workingDir = "D:\\NewWGCNA";
setwd(workingDir); 

#Reading Data*******************************************************************
D1=data.frame(read_xlsx("D:\\NewWGCNA\\Data_all.xlsx"))
id=colnames(D1)[-1]; group=c(rep(0,30),rep(1,44))
gnames=D1[,1]
D2=t(D1[,-1])
colnames(D2)=gnames
rownames(D2)=id
x=D2


####X=t(read.table("D\\Data.txt", sep='\t', header=T,fill=T))
#Filtering Data*****************************************************************
FX1<-varFilter(as.matrix(t(x)), var.func=IQR, var.cutoff=0.95, filterByQuantile=TRUE);dim(FX1)
FX=data.frame(t(FX1));dim(FX)
write.table(FX,"FX.txt")

FX=read.table("D:\\NewWGCNA\\FX.txt")
dim(FX)
#=====================================================================================
#  Is Data Good Enough?
#=====================================================================================
gsg = goodSamplesGenes(FX, verbose = 3);
gsg$allOK


#=====================================================================================
#  Hierarchical Clustering
#=====================================================================================
sampleTree = hclust(dist(FX), method = "average");
sizeGrWindow(12,9)# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot the sample dendrogram and the colors underneath.
MyTraitData=group
# Re-cluster samples
sampleTree2 = hclust(dist(FX), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(MyTraitData, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(MyTraitData), 
                    main = "Sample dendrogram and trait heatmap")
#***************************************************************************
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()

#=====================================================================================
#  Choose a set of soft-thresholding powers
#=====================================================================================
powers = c(c(1:20))
# Call the network topology analysis function
sft = pickSoftThreshold(FX, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net= blockwiseModules(FX, power = 4,
                      TOMType = "unsigned", minModuleSize = 30,
                      reassignThreshold = 0, mergeCutHeight = 0.25,
                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                      verbose = 0)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
softPower = 4;
adjacency = adjacency(FX, power = softPower);
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



# Calculate eigengenes
MEList = moduleEigengenes(FX, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(FX, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")


write.xlsx(table(moduleColors),"D:\\NewWGCNA\\NetColours.xlsx" )
ModuleCL=cbind(moduleLabels, moduleColors)
write.xlsx(as.matrix(MEs),"D:\\NewWGCNA\\NetME.xlsx")


# Calculate eigengenes
MEList = moduleEigengenes(FX, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")



#*********************************************************************************
#Relating modules to external clinical traits
#*********************************************************************
nGenes = ncol(FX);
nSamples = nrow(FX)
MEs0 = moduleEigengenes(FX, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, MyTraitData, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
CorPval=cbind(moduleTraitCor,moduleTraitPvalue);colnames(CorPval)=c("Correlation","P-Value")
write.xlsx(CorPval,"moduleTraitCorrPvalue.xlsx")

sizeGrWindow(50,20)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Group",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


geneModuleMembership=as.data.frame(cor(FX, MEs, use = "p"))
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
write.xlsx(data.frame(moduleColors,moduleLabels,geneModuleMembership),"GM_cor.xlsx")
write.xlsx(data.frame(moduleColors,moduleLabels,MMPvalue),"GM_p.xlsx")
cors=c(0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95)
modss=names(geneModuleMembership)[-14]
corout=matrix(0,length(cors),13)
for( i in 1:8){
        for( j in 1:8){
        corout[i,j]=length(geneModuleMembership[,j][geneModuleMembership[,j]>cors[i]])
}}
dimnames(corout)=list(cors,modss)

#Hub Genes
HG_unsigned=chooseTopHubInEachModule(FX,moduleColors,omitColors = "grey",power = 4,type="unsigned") 
HG_signed=chooseTopHubInEachModule(FX,moduleColors,omitColors = "grey",power = 4,type="signed") 
HG=data.frame(HG_unsigned,HG_signed)
write.xlsx(HG,"HG.xlsx")


intModules = c("brown", "red", "blue",
               "turquoise","yellow","green",
               "black","pink")
for (module in modss)
{
        modGenes = (moduleColors==module)
        modLLIDs = probes[modGenes];
        fileName = paste("Probes-", module, ".xlsx");
        write.xlsx(as.data.frame(modLLIDs),file = fileName)
}
#=====================================================================================
#
#  Code chunk 3 PPPPPPProblem ???!!!!
#
#=====================================================================================


GOenr = GOenrichmentAnalysis(moduleColors, probes, organism = "human", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab









#*************************
#Intra-modular analysis: identifying genes with high GS and MM
#************************


#Do it for turquoise, blue, brown, yellow, green, red



modNames = substring(names(MEs), 3)
for(i in 1:13){
        mypath <- file.path("D:",paste("myplot_", modNames[i], ".jpg", sep = ""))
        module = modNames[i]; 
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        sizeGrWindow(10, 10)
        geneTraitSignificance = as.data.frame(cor(FX, group, use = "p"));
        jpeg(file=mypath)
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Group",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        dev.off()
}


#*****************************
#Visualizing the gene network
#*****************************
plotTOM = dissTOM^4;
diag(plotTOM) = NA;
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, main = "Network heatmap plot, all genes")


#VisAnt
modNames = substring(names(MEs), 3)
#"magenta"     "greenyellow" "salmon"      "green"       "red"         "tan"         "blue"        "brown"       "turquoise"  
#"yellow"      "pink"        "black"       "purple"      "grey" 
kIM = intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE) 
write.xlsx(kIM,"KIM.xlsx")
KIM1=read.xlsx("D:\\NewWGCNA\\KIM.xlsx",sheetIndex = 1)
dim(KIM1)


module ="turquoise" ;probes = names(FX);inModule = (moduleColors==module);modProbes = probes[inModule];
dissmodTOM = dissTOM[inModule, inModule];modTOM =1- dissmodTOM
dimnames(modTOM) = list(modProbes, modProbes);IMConn = softConnectivity(FX[, modProbes]);
vis = exportNetworkToVisANT(modTOM,file = paste("VisANTInput-", "turquoise", ".txt", sep=""),weighted = TRUE,threshold = 0 )
gre1=data.frame(cbind(vis$from,vis$to,vis$weight))
colnames(gre1)=c("Gene 1","Gene 2","Connection weight")
gre2=data.frame(cbind(substring(gre1$`Gene 1`, 5),substring(gre1$`Gene 2`, 5),gre1$`Connection weight`))
write.xlsx(gre2,"D:\\NewWGCNA\\Module_Weights_turquoise.xlsx") 

kIM = intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE) 





# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


