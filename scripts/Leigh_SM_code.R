
# Seurat guided analysis of the axolotl limb
#This markdown takes cell by gene matrices from each individual time point sampled, performs quality control, filters out low quality cells, #and then determines clusters and cluster-enriched markers. All mar

# Homeostatic limb
#This section analyses the homeostatic limb, which includes two samples. 

#First, load the required packages. 
library(Seurat)
library(dplyr)


#Now, load the cell by gene matrix that includes the two homeostatic limb samples. This cell by gene matrix has three other samples, that we #will not use right now. So we need to select out only samples 1 and 2. We will use intact to as shorthand for homeostatic limbs. 

#load in data
inDrops3.data = read.table('intact_and_contralateral.repGene', header = T, row.names = 1, sep = '\t')
#pull out samples 1 and 2, which are the intact limb samples
inDrops3.intact = inDrops3.data[,grep('^S[12]_', colnames(inDrops3.data))]

#Remove data matrix with extra samples
rm(inDrops3.data) 
#Create Seurat object and make sparse
seurat_inDrops3.intact = CreateSeuratObject(inDrops3.intact, project = 'inDrops3.intact', min.cells = 8, min.genes = 200)
seurat_inDrops3.intact = MakeSparse(seurat_inDrops3.intact)

#While we did some filtering above, we need to perform further quality control to ensure that the cells we are working with aren't apoptotic #or have a dearth of genes. First, we need to identify the mitochondrial genes present in this matrix. The axolotl mitochondrial genome can #be found here: https://www.ncbi.nlm.nih.gov/nuccore/AJ584639. Remember that the genes are written as protein names when greping for #mitochondrial genes. 

#find mitochonrial genes in matrix. The protein name should be used and changed for each gene within the mitochondrial genome.
grep(pattern = "*CYB_*", x = rownames(x = seurat_inDrops3.intact@data), value = TRUE)
#list of all mitochondrial genes in this intact matrix
mito.genes.intact <- c("c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1025234_g1_i1^sp|O63796|NU2M_ANACA", "c1068681_g4_i1^sp|Q9ZZM3|NU5M_SALSA^sp|P82013|VDAC2_MELGA^Porin_3", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.intact <- Matrix::colSums(seurat_inDrops3.intact@raw.data[mito.genes.intact, ])/Matrix::colSums(seurat_inDrops3.intact@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_inDrops3.intact <- AddMetaData(object = seurat_inDrops3.intact, metadata = percent.mito.intact, col.name = "percent.mito")

#Now perform quality control on matrix by filtering out cells with high percent mitochondrial RNA and low and high number of genes. We filter #out cells that by visualize inspection appear to have relatively high mitochondrial RNA content or high or low number of genes. These #numbers can be modified to be more or less inclusive. 

#visualize number of genes, unique molecular identifiers (UMI), and percent mitochondrial RNA
VlnPlot(object = seurat_inDrops3.intact, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#filter out cells
seurat_inDrops3.intact <- FilterCells(object = seurat_inDrops3.intact, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(4000, 0.125))

#normalize data
seurat_inDrops3.intact <- NormalizeData(seurat_inDrops3.intact, normalization.method= "LogNormalize", scale.factor= 10000)

#find variable genes
seurat_inDrops3.intact <- FindVariableGenes(object = seurat_inDrops3.intact, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)

#scale data and regress out nUMI and percent.mito
seurat_inDrops3.intact <- ScaleData(seurat_inDrops3.intact, vars.to.regress = c('nUMI', 'percent.mito'))


#Next, we perform linear dimensional reduction and visualize the results in a few different ways. 

seurat_inDrops3.intact <- RunPCA(object = seurat_inDrops3.intact, pc.genes = seurat_inDrops3.intact@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#visualize results
PrintPCA(object = seurat_inDrops3.intact, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = seurat_inDrops3.intact, pcs.use = 1:2)
PCAPlot(object = seurat_inDrops3.intact, dim.1 = 1, dim.2 = 2)
seurat_inDrops3.intact <- ProjectPCA(object = seurat_inDrops3.intact, do.print = FALSE)
PCHeatmap(object = seurat_inDrops3.intact, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seurat_inDrops3.intact, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = seurat_inDrops3.intact, pc.use = 13:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

#plot standard deviations to chose PCs to use in downstream analysis, here we chose 18
PCElbowPlot(object = seurat_inDrops3.intact)


#Now we can identify cell populations within the homeostatic limb, vizualize the resulting populations using tSNE, and subsequently find markers that define these different populations 

#find clusters using first 18 PCs
seurat_inDrops3.intact <- FindClusters(object = seurat_inDrops3.intact, reduction.type = "pca", dims.use = 1:18, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#run non-linear dimensional reduction
seurat_inDrops3.intact <- RunTSNE(object = seurat_inDrops3.intact, dims.use = 1:18, do.fast = TRUE)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
seurat_inDrops3.intact <- BuildClusterTree(seurat_inDrops3.intact, do.reorder=TRUE, reorder.numeric=TRUE)

#visualize tSNE 
set.seed(5)
TSNEPlot(object = seurat_inDrops3.intact, do.label = T)

#visulize tSNE based on sample to determine how similar the two samples are to one another
TSNEPlot(object = seurat_inDrops3.intact, group.by = 'orig.ident')

#assess nodes
node.scores <- AssessNodes(seurat_inDrops3.intact)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores


#merge first 7 nodes
#select nodes to merge
nodes.merge <- node.scores[1:7, ]
nodes.to.merge <- sort(x = nodes.merge$node)

#create a new Seurat object in which we will merge our selected nodes
merged <- seurat_inDrops3.intact
#merge nodes
for (n in nodes.to.merge) {merged <- MergeNode(object = merged, node.use = n)}


#re-visualize the tSNE after we have merged the non-distinct nodes
set.seed(5)
TSNEPlot(merged, do.label = TRUE)


#determine differentially expressed genes for each population

#find markers for each population
all.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#write DE results to table for inspection
write.table(all.markers, 'intact.only.markers.txt', sep = '\t')

#save Rdata
save.image('intact.Rdata')

#end session
q()

#Wound healing samples
#This section analyses three samples collected during wound healing (3 days post amputation) 

#First, load the required packages. 
library(Seurat)
library(dplyr)

#load in data
data = read.table('wound_healing_and_medium_bud.repGene', header=T, row.names=1, sep='\t')
#pull out samples N4, N5, and N6 which are the wound healing stage limb samples
N4_N5_N6 = data[,grep('^N[456]', colnames(data))]

#Remove data matrix with extra samples
rm(data)
#Create Seurat object and make sparse
seurat_3dpa = CreateSeuratObject(N4_N5_N6, project = '3dpa', min.cells = 8, min.genes = 200)
seurat_3dpa = MakeSparse(seurat_3dpa)

#While we did some filtering above, we need to perform further quality control to ensure that the cells we are working with aren't apoptotic #or have a dearth of genes. First, we need to identify the mitochondrial genes present in this matrix. The axolotl mitochondrial genome can #be found here: https://www.ncbi.nlm.nih.gov/nuccore/AJ584639. Remember that the genes are written as protein names when greping for #mitochondrial genes. 

#find mitochonrial genes in matrix. The protein name should be used and changed for each gene within the mitochondrial genome.
grep(pattern = "*CYB_*", x = rownames(x = seurat_3dpa.intact@data), value = TRUE)
#list of all mitochondrial genes in this wound healing matrix
mito.genes.3dpa <- c("c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c1057599_g1_i2^sp|P00018|CYC_DRONO^Cytochrome_CBB3", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c1053715_g3_i1^sp|O03539|COX1_NOTPE", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c289614_g1_i1^sp|P05503|COX1_RAT", "c959712_g1_i1^sp|P00416|COX3_MOUSE", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1025234_g1_i1^sp|O63796|NU2M_ANACA", "c1068681_g4_i4^sp|Q9ZZM3|NU5M_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.3dpa <- Matrix::colSums(seurat_3dpa@raw.data[mito.genes.3dpa, ])/Matrix::colSums(seurat_3dpa@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_3dpa <- AddMetaData(object = seurat_3dpa, metadata = percent.mito.3dpa, col.name = "percent.mito")

#Now perform quality control on matrix by filtering out cells with high percent mitochondrial RNA and low and high number of genes. We filter #out cells that by visualize inspection appear to have relatively high mitochondrial RNA content or high or low number of genes. These #numbers can be modified to be more or less inclusive. 

#visualize number of genes, unique molecular identifiers (UMI), and percent mitochondrial RNA
VlnPlot(object = seurat_3dpa, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#filter out cells
seurat_3dpa <- FilterCells(object = seurat_3dpa, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(5000, 0.15))

#normalize data
seurat_3dpa <- NormalizeData(object = seurat_3dpa, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes
seurat_3dpa <- FindVariableGenes(object = seurat_3dpa, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#scale data and regress out nUMI and percent.mito
seurat_3dpa <- ScaleData(object = seurat_3dpa, vars.to.regress = c("nUMI", "percent.mito"))


#Next, we perform linear dimensional reduction and visualize the results in a few different ways. 
seurat_3dpa <- RunPCA(object = seurat_3dpa, pc.genes = seurat_3dpa@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#visualize results
VizPCA(object = seurat_3dpa, pcs.use = 1:2)
PCAPlot(object = seurat_3dpa, dim.1 = 1, dim.2 = 2)
seurat_3dpa <- ProjectPCA(object = seurat_3dpa, do.print = FALSE)
PCHeatmap(object = seurat_3dpa, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seurat_3dpa, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = seurat_3dpa)

#plot standard deviations to chose PCs to use in downstream analysis, here we chose 19
seurat_3dpa <- FindClusters(object = seurat_3dpa, reduction.type = "pca", dims.use = 1:19, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#Now we can identify cell populations within the homeostatic limb, vizualize the resulting populations using tSNE, and subsequently find markers that define these different populations 

#find clusters using first 18 PCs
seurat_inDrops3.intact <- FindClusters(object = seurat_inDrops3.intact, reduction.type = "pca", dims.use = 1:18, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#run non-linear dimensional reduction
seurat_3dpa <- RunTSNE(object = seurat_3dpa, dims.use = 1:19, do.fast = TRUE)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
seurat_3dpa <- BuildClusterTree(seurat_3dpa, do.reorder=TRUE, reorder.numeric=TRUE)

#visualize tSNE 
set.seed(5)
TSNEPlot(object = seurat_3dpa)
#visulize tSNE based on sample to determine how similar the two samples are to one another
TSNEPlot(object = seurat_3dpa, group.by = "orig.ident")

#assess nodes
node.scores <- AssessNodes(seurat_3dpa)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores


#merge first 5 nodes
#select nodes to merge
nodes.merge=node.scores[1:5,]
nodes.to.merge <- sort(x = nodes.merge$node)

#create a new Seurat object in which we will merge our selected nodes
merged <- seurat_3dpa
#merge nodes
for (n in nodes.to.merge) {merged <- MergeNode(object = merged, node.use = n)}


#re-visualize the tSNE after we have merged the non-distinct nodes
set.seed(5)
TSNEPlot(merged, do.label = TRUE)


#determine differentially expressed genes for each population

#find markers for each population
all.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#write DE results to table for inspection
write.table(all.markers, 'wound.healing.markers.txt', sep = '\t')

#save Rdata
save.image('wound_healing.Rdata')

#end session
q()

#Early-bud stage blastema samples
#This section analyses two samples collected during early bud stage (14 days post amputation) 

#First, load the required packages. 
library(Seurat)
library(dplyr)

#load in data
inDrops4.d1 = read.table('early_and_medium_bud.repGene', header = T, row.names = 1, sep = '\t')
#pull out samples S3 and S5 which are the early-bud blastema samples
S3_S5 = inDrops4.d1[,grep('^S[35]_', colnames(inDrops4.d1))]

#Remove data matrix with extra samples
rm(inDrops4.d1)
#Create Seurat object and make sparse
seurat_14dpa = CreateSeuratObject(S3_S5, project = '14dpa', min.cells = 8, min.genes = 200)
seurat_14dpa = MakeSparse(seurat_14dpa)

#While we did some filtering above, we need to perform further quality control to ensure that the cells we are working with aren't apoptotic #or have a dearth of genes. First, we need to identify the mitochondrial genes present in this matrix. The axolotl mitochondrial genome can #be found here: https://www.ncbi.nlm.nih.gov/nuccore/AJ584639. Remember that the genes are written as protein names when greping for #mitochondrial genes. 

#list of all mitochondrial genes in this early-bud blastema matrix
mito.genes.14dpa <- c("c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1068681_g4_i1^sp|Q9ZZM3|NU5M_SALSA^sp|P82013|VDAC2_MELGA^Porin_3", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.14dpa <- Matrix::colSums(seurat_14dpa@raw.data[mito.genes.14dpa, ])/Matrix::colSums(seurat_14dpa@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_14dpa <- AddMetaData(object = seurat_14dpa, metadata = percent.mito.14dpa, col.name = "percent.mito")

#Now perform quality control on matrix by filtering out cells with high percent mitochondrial RNA and low and high number of genes. We filter #out cells that by visualize inspection appear to have relatively high mitochondrial RNA content or high or low number of genes. These #numbers can be modified to be more or less inclusive. 

#visualize number of genes, unique molecular identifiers (UMI), and percent mitochondrial RNA
VlnPlot(object = seurat_14dpa, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#filter out cells
seurat_14dpa <- FilterCells(object = seurat_14dpa, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(5000, 0.10))

#normalize data
seurat_14dpa <- NormalizeData(object = seurat_14dpa, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes
seurat_14dpa <- FindVariableGenes(object = seurat_14dpa, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#scale data and regress out nUMI and percent.mito
seurat_14dpa <- ScaleData(object = seurat_14dpa, vars.to.regress = c("nUMI", "percent.mito"))


#Next, we perform linear dimensional reduction and visualize the results in a few different ways. 
seurat_14dpa <- RunPCA(object = seurat_14dpa, pc.genes = seurat_14dpa@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#visualize results
VizPCA(object = seurat_14dpa, pcs.use = 1:2)
PCAPlot(object = seurat_14dpa, dim.1 = 1, dim.2 = 2)
seurat_14dpa <- ProjectPCA(object = seurat_14dpa, do.print = FALSE)
PCHeatmap(object = seurat_14dpa, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seurat_14dpa, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = seurat_14dpa)

#plot standard deviations to chose PCs to use in downstream analysis, here we chose 19
seurat_14dpa <- FindClusters(object = seurat_14dpa, reduction.type = "pca", dims.use = 1:19, resolution = 1, print.output = 0, save.SNN = TRUE)


#Now we can identify cell populations within the homeostatic limb, vizualize the resulting populations using tSNE, and subsequently find markers that define these different populations 

#find clusters using first 18 PCs
seurat_inDrops3.intact <- FindClusters(object = seurat_inDrops3.intact, reduction.type = "pca", dims.use = 1:18, resolution = 1.5, print.output = 0, save.SNN = TRUE)

#run non-linear dimensional reduction
seurat_14dpa <- RunTSNE(object = seurat_14dpa, dims.use = 1:19, do.fast = TRUE)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
seurat_14dpa <- BuildClusterTree(seurat_14dpa, do.reorder=TRUE, reorder.numeric=TRUE)

#visualize tSNE 
set.seed(5)
TSNEPlot(object = seurat_14dpa)
#visulize tSNE based on sample to determine how similar the two samples are to one another
TSNEPlot(object = seurat_14dpa, group.by = "orig.ident")

#assess nodes
node.scores <- AssessNodes(seurat_14dpa)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores


#merge first nodes
#select nodes to merge
nodes.merge=node.scores[1:1,]
nodes.to.merge <- sort(x = nodes.merge$node)

#create a new Seurat object in which we will merge our selected nodes
merged <- seurat_14dpa
#merge nodes
for (n in nodes.to.merge) {merged <- MergeNode(object = merged, node.use = n)}


#re-visualize the tSNE after we have merged the non-distinct nodes
set.seed(5)
TSNEPlot(merged, do.label = TRUE)


#determine differentially expressed genes for each population

#find markers for each population
all.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#write DE results to table for inspection
write.table(all.markers, 'early_bud_blastema.markers.txt', sep = '\t')

#save Rdata
save.image('early_bud_blastema.Rdata')

#end session
q()

#Medium-bud stage blastema samples
#This section analyses six samples collected during medium bud stage (23 days post amputation) 

#First, load the required packages. 
library(Seurat)
library(dplyr)

#load in first three medium-bud samples
data.d1 = read.table("early_and_medium_bud.repGene", header=T, row.names=1, sep='\t')
#sample S1, S2, and S4 are medium-bud samples collected on day 1
S1_S2_S4 = data.d1[,grep("^S[124]_",colnames(data.d1))]

#create Seurat object and make sparse
seurat_S1_S2_S4 = CreateSeuratObject(raw.data = S1_S2_S4, project = "23dpa_d1", min.cells = 8, min.genes = 200)
seurat_S1_S2_S4 <- MakeSparse(seurat_S1_S2_S4) 

#list mito genes medium-bud day 1 matrix
mito.genes.i4.d1 <- c("c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1057599_g1_i1^sp|P00018|CYC_DRONO^Cytochrome_CBB3", "c1127119_g1_i1^sp|P81280|CYC_ALLMI", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c959712_g1_i1^sp|P00416|COX3_MOUSE", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1060846_g2_i1^sp|P03921|NU5M_MOUSE", "c1068681_g4_i1^sp|Q9ZZM3|NU5M_SALSA^sp|P82013|VDAC2_MELGA^Porin_3")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.i4.d1 <- Matrix::colSums(seurat_S1_S2_S4@raw.data[mito.genes.i4.d1, ])/Matrix::colSums(seurat_S1_S2_S4@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_S1_S2_S4 <- AddMetaData(object = seurat_S1_S2_S4, metadata = percent.mito.i4.d1, col.name = "percent.mito") 

#remove raw data
rm(data.d1, S1_S2_S4)

#load in second set of three medium-bud samples collected on day 2
data.d2 = read.table('wound_healing_and_medium_bud.repGene', header=T, row.names=1, sep='\t')
#samples N1, N2, and N3 are medium-bud stage blastema samples collected on day 2
N1_N2_N3 = data.d2[,grep("^N[123]",colnames(data.d2))]

#create seurat object and make sparse
seurat_N1_N2_N3 = CreateSeuratObject(raw.data = N1_N2_N3, project = "23dpa_d2", min.cells = 8, min.genes = 200)
seurat_N1_N2_N3 <- MakeSparse(seurat_N1_N2_N3)

#list mito genes from day 2 matrix
mito.genes.i4.d2 <- c("c786641_g1_i1^sp|Q9B205|CYB_CAICR", "c1084180_g1_i1^sp|Q8LWP6|CYB_RANSI", "c1043008_g2_i1^sp|Q9B205|CYB_CAICR", "c1060846_g1_i1^sp|Q8WA47|CYB_MUSMA", "c1084180_g3_i1^sp|Q8LWP6|CYB_RANSI", "c1027109_g1_i1^sp|Q35920|ATP6_SALSA", "c1088733_g1_i1^sp|Q9ZZM6|COX1_SALSA", "c220469_g1_i1^sp|P00397|COX1_MOUSE", "c1451851_g1_i1^sp|Q9ZXY2|COX1_PAPHA", "c289614_g1_i1^sp|P05503|COX1_RAT", "c959712_g1_i1^sp|P00416|COX3_MOUSE", "c1049442_g1_i1^sp|Q96133|COX3_CARAU", "c1083417_g1_i2^sp|P00419|COX3_XENLA", "c934922_g1_i1^sp|Q9ZXX8|COX3_PAPHA", "c1083535_g6_i1^sp|Q4JQI7|NU1M_TETNG", "c1060846_g2_i1^sp|P03921|NU5M_MOUSE", "c1068681_g4_i4^sp|Q9ZZM3|NU5M_SALSA")

#calculate the percentage mitochondrial RNA for each cell
percent.mito.i4.d2 <- Matrix::colSums(seurat_N1_N2_N3@raw.data[mito.genes.i4.d2, ])/Matrix::colSums(seurat_N1_N2_N3@raw.data)
#add the percent mitochondrial content of each cell to the Seurat object
seurat_N1_N2_N3 <- AddMetaData(object = seurat_N1_N2_N3, metadata = percent.mito.i4.d2, col.name = "percent.mito")

#save each sample info in a metadata sample slot
seurat_N1_N2_N3 <- StashIdent(seurat_N1_N2_N3, save.name = 'sample')
seurat_S1_S2_S4 <- StashIdent(seurat_S1_S2_S4, save.name = 'sample')

#custom filters
seurat_S1_S2_S4 <- FilterCells(object = seurat_S1_S2_S4, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(6000, 0.1))
seurat_N1_N2_N3 <- FilterCells(object = seurat_N1_N2_N3, subset.names = c("nGene", "percent.mito"), low.thresholds = c(850, -Inf), high.thresholds = c(6000, 0.1))

#normalize data
seurat_N1_N2_N3 <- NormalizeData(object = seurat_N1_N2_N3, normalization.method = "LogNormalize", scale.factor= 10000)
seurat_S1_S2_S4 <- NormalizeData(object = seurat_S1_S2_S4, normalization.method = "LogNormalize", scale.factor= 10000)

#find variable genes
seurat_N1_N2_N3 <- FindVariableGenes(object = seurat_N1_N2_N3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
seurat_S1_S2_S4 <- FindVariableGenes(object = seurat_S1_S2_S4, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#scale data and regress out nUMI and percent.mito
seurat_N1_N2_N3 <- ScaleData(object = seurat_N1_N2_N3, vars.to.regress = c('nUMI', 'percent.mito'))
seurat_S1_S2_S4 <- ScaleData(object = seurat_S1_S2_S4, vars.to.regress = c('nUMI', 'percent.mito'))

#add day collected information into a metadata column
seurat_S1_S2_S4@meta.data$day <- "d1"
seurat_N1_N2_N3@meta.data$day <- "d2"

#find highly variable genes
hvg.d1 <- rownames(x = head(x = seurat_S1_S2_S4@hvg.info, n = 2000))
hvg.d2 <- rownames(x = head(x = seurat_N1_N2_N3@hvg.info, n = 2000))
hvg.union <- union(x = hvg.d1, y = hvg.d2)

#run canonical correlation analysis and save resulting matrix as a Seurat object called combined
combined <- RunCCA(object = seurat_S1_S2_S4, object2= seurat_N1_N2_N3, genes.use = hvg.union, num.cc = 30)

#determine cca dimensions to use in downstream analysis
DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 10:19, do.balanced = TRUE)
DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, dim.use = 20:29, do.balanced = TRUE)

#we chose first 25 dimensions
combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "day", dims.align = 1:25)

#perform integrated analysis on all cells
combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:25, do.fast = T)
combined <- FindClusters(combined, reduction.type = "cca.aligned", resolution = 1, dims.use = 1:25)

# Build a phylogenetic tree to see how cells are related while simultaneously renaming and reordering cluster names according to their #position on the tree. This will be important to determine when deciding whether similar populations should be merged. 
combined <- BuildClusterTree(combined, do.reorder=TRUE, reorder.numeric=TRUE)

#assess nodes
node.scores <- AssessNodes(combined)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
node.scores

#merge first 2 nodes
nodes.merge <- node.scores[1:2, ] 
nodes.to.merge <- sort(x = nodes.merge$node)
combined.merged <- combined 
for (n in nodes.to.merge) {combined.merged <- MergeNode(object = combined.merged, node.use = n) }

#visualize tSNE by sample, day collected, and assigned cluster 
TSNEPlot(combined.merged, do.return = T, pt.size = 0.5, group.by = "sample")
TSNEPlot(combined.merged, do.return = T, pt.size = 0.5, group.by = "day")
TSNEPlot(combined.merged, do.label = T, do.return = T, pt.size = 0.5)

#find DE genes for each population
all.markers.merged <- FindAllMarkers(object = combined.merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, logfc.threshold = 0.35, max.cells.per.ident = 2000)

#write DE results to table for inspection
write.table(all.markers.merged, 'all.markers.medium.bud.txt', sep = '\t') 

save.image('medium_bud_blastema.RData')

q()

#Pseudotime analysis of blastema-resident cells

#load required packages
library(URD)
library(Seurat)
library(dplyr)
#load medium-bud data
load('medium_bud_blastema.RData')

#pull out all non-immune blastema cells
dpa23_FDB <- SubsetData(combined.merged, ident.use=c('5','6','10', '13', '14','15', '16','17'), subset.name = 'day', accept.value = 'd1')
#remove medium-bud Seurat object
rm(combined.merged)

#add time point metadata to Seurat object
dpa23_FDB@meta.data$time <- "medium-bud"

#load in early-bud data
load('early_bud_blastema.Rdata')

#pull out all non-immune blastema cells
dpa14_blastema <- SubsetData(merged, ident.use= c('3', '4', '5','6','11'))
#add time point metadata to Seurat object
dpa14_blastema@meta.data$time <- "early-bud"

#remove excess data
rm(seurat_14dpa, percent.mito.14dpa, mito.genes.14dpa, merged)

#load in wound healing data
load('wound_healing.Rdata')
#pull out all non-immune blastema cells
dpa3_blastema <- SubsetData(merged, ident.use = c('2','21','24','25'))

#add time point metadata to Seurat object
dpa3_blastema@meta.data$time <- 'wound_healing'

#remove excess data
rm(seurat_3dpa, mito.genes.3dpa, percent.mito.3dpa, node.scores, merged, all.markers, nodes.merge, nodes.to.merge, n)

#combine above Seurat object so we can make a raw data file and metadata file with only the cells we want to put into URD
dpa14.23 <- MergeSeurat(dpa14_blastema, dpa23_FDB)
combined <- MergeSeurat(dpa14.23,dpa3_blastema)

#pull out raw data
raw.data <- as.matrix(combined@raw.data)

#pull out metadata
meta.data <- as.matrix(combined@meta.data)

#remove excess Seurat objects
rm(dpa14.23, dpa3_blastema, dpa14_blastema, dpa23_FDB)

#create URD object with all of the cells we selected above
blastema <- createURD(count.data = raw.data, meta = meta.data, min.cells=3, min.counts=3)

#now remove raw data and Seurat objects
rm(raw.data, meta.data, combined)

#add time point info to URD stage slot
blastema@group.ids$stage <- as.character(blastema@meta[rownames(blastema@group.ids),"time"])

#find variable genes for each time point
stages <- sort(unique(blastema@group.ids$stage))
var.3dpa <- findVariableGenes(blastema, cells.fit=cellsInCluster(blastema, "stage", 'wound_healing'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.14dpa <- findVariableGenes(blastema, cells.fit=cellsInCluster(blastema, "stage", 'early-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.23dpa <- findVariableGenes(blastema, cells.fit=cellsInCluster(blastema, "stage", 'medium-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)

#combine into one 
var.genes <- sort(unique(unlist(var.3dpa, var.14dpa, var.23dpa)))

#add to URD object
blastema@var.genes <- var.genes

#calculate PCA
blastema <- calcPCA(blastema, mp.factor = 2)

pcSDPlot(blastema)

# Calculate tSNE
set.seed(19)
blastema <- calcTsne(object = blastema)

#visualize tSNE by time point
plotDim(blastema, "time", plot.title = "tSNE: DPA")

#create an URD object with only medium-bud cells
blastema.23dpa <- urdSubset(blastema, cells.keep=cellsInCluster(blastema, "stage", "medium-bud"))

#use variable genes that were calculated only for medium-bud
blastema.23dpa@var.genes <- var.23dpa

#calculate PCA and tSNE
blastema.23dpa <- calcPCA(blastema.23dpa, mp.factor = 1.5)
pcSDPlot(blastema.23dpa)
set.seed(20)
blastema.23dpa <- calcTsne(blastema.23dpa)

#perform graph based clustering
blastema.23dpa <- graphClustering(blastema.23dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
blastema.23dpa <- graphClustering(blastema.23dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(blastema.23dpa, c, legend=F))}

#we chose to go forward with Infomap-100
clusters <- unique(blastema.23dpa@group.ids$'Infomap-100')

#find markers for these populations to get an idea of what they are
pr.markers <- lapply(clusters, function(c) markersAUCPR(blastema.23dpa, clust.1 = c, clustering = 'Infomap-100', genes.use= blastema.23dpa@var.genes))

#you can export each table like below, or if R studio look at top markers in viewer
#note that cluster number and number in pr.markers[[]] are not necessarily equal!
#cluster number can be found in headings of the marker lists
write.table(pr.markers[[1]], 'clus1.markers.txt', sep = '\t')

#looking at cluster genes to see if any cells we don't want to include may have made it through prior clustering
head(pr.markers[[3]], 20)

#we will not use clusters 7 (Myeloid cells), 10 (Wound epidermis), 13 (Erythrocytes) so these aren't included in trajectory analysis
#get cell names for the cells we will use in URD
dpa23_good_cells <- cellsInCluster(blastema.23dpa, clustering = 'Infomap-100', cluster = c('1','2','3','4','5','6','8','9','11','12'))

#let's clean up the other two time points since it's clear that WE, erythrocytes, etc. may be contaminating the clusters
#create URD object with just early-bud sample
blastema.14dpa <- urdSubset(blastema, cells.keep=cellsInCluster(blastema, "stage", "early-bud"))

#use variable genes found for early-bud
blastema.14dpa@var.genes <- var.14dpa

#calculate PCA and tSNE
blastema.14dpa <- calcPCA(blastema.14dpa, mp.factor = 1.5)
pcSDPlot(blastema.14dpa)
set.seed(20)
blastema.14dpa <- calcTsne(blastema.14dpa)

#perform graph-based clustering
blastema.14dpa <- graphClustering(blastema.14dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
blastema.14dpa <- graphClustering(blastema.14dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(blastema.14dpa, c, legend=F))}

#iwe chose to go forward with Infomap-50 based clustering
clusters <- unique(blastema.14dpa@group.ids$'Infomap-50')

#calculate marker genes
pr.markers_14dpa <- lapply(clusters, function(c) markersAUCPR(blastema.14dpa, clust.1 = c, clustering = 'Infomap-50', genes.use= blastema.14dpa@var.genes))

#inspect genes that define subsets, remembering that cluster number is found in column names
#either or export to table or inspect in R/Rstudio
write.table(pr.markers_14dpa[[1]], '081218.clus1.txt', sep = '\t')
#looking at cluster genes
head(pr.markers_14dpa[[10]])

#need to remove clusters 1 (WE), 7 (WE), 11(Erythrocytes), 13(Immune) which would be 217 cells that likely wont participate in this
#get names of cells we will use in URD
dpa14_good_cells <- cellsInCluster(blastema.14dpa, clustering = 'Infomap-50', cluster = c('2','3','4','5','6','8','9','10','12','14','15','16')) 

#finally, let's clean the wound healing cells
#create a subsetted object of cells from wound healing
blastema.3dpa <- urdSubset(blastema, cells.keep=cellsInCluster(blastema, "stage", "wound_healing"))

#use wound healing variable genes
blastema.3dpa@var.genes <- var.14dpa

#calculate PCA and tSNE
blastema.3dpa <- calcPCA(blastema.3dpa, mp.factor = 1.5)

pcSDPlot(blastema.3dpa)
set.seed(20)
blastema.3dpa <- calcTsne(blastema.3dpa)

#graph-based clustering
blastema.3dpa <- graphClustering(blastema.3dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
blastema.3dpa <- graphClustering(blastema.3dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(blastema.3dpa, c, legend=F))}


#we chose to go forward with Infomap-60
clusters <- unique(blastema.3dpa@group.ids$'Infomap-60')

#find marker genes
pr.markers_3dpa <- lapply(clusters, function(c) markersAUCPR(blastema.3dpa, clust.1 = c, clustering = 'Infomap-60', genes.use= blastema.3dpa@var.genes))

#inspect marker genes to determine what to remove
head(pr.markers_3dpa[[8]])

#we will remove 3 (WE) and 12 (erythrocytes)
dpa3_good_cells <- cellsInCluster(blastema.3dpa, clustering = 'Infomap-60', cluster = c('1','2','4','5','6','7','8','9','10','11','13'))

#collect all cells to use in URD
paste(dpa23_good_cells, dpa14_good_cells) -> cells_to_use
paste(cells_to_use, dpa3_good_cells) -> cells_for_urd


#make an URD object with the cells we want to use
cleaned <- urdSubset(blastema, cells.keep= c(dpa3_good_cells, dpa23_good_cells, dpa14_good_cells))

#re-calculate variable genes for these
stages <- sort(unique(cleaned@group.ids$stage))
var.3dpa <- findVariableGenes(cleaned, cells.fit=cellsInCluster(cleaned, "stage", 'wound_healing'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.14dpa <- findVariableGenes(cleaned, cells.fit=cellsInCluster(cleaned, "stage", 'early-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.23dpa <- findVariableGenes(cleaned, cells.fit=cellsInCluster(cleaned, "stage", 'medium-bud'), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, do.plot = T)
var.genes <- sort(unique(unlist(var.3dpa, var.14dpa, var.23dpa)))
cleaned@var.genes <- var.genes

#calculate PCA
cleaned <- calcPCA(cleaned, mp.factor = 2)
pcSDPlot(cleaned)

#calculate tSNE, see Figure 7a
set.seed(19)
cleaned <- calcTsne(object = cleaned)
plotDim(cleaned, "time", plot.title = "tSNE by time point", legend = F, point.size = 2)

#calculate diffision map, allowing destinty to determine sigma (this value was determined to be 19.538, which we used)
cleaned <- calcDM(cleaned, knn = 54, sigma = 19.538)

#visualize dim arrays
plotDimArray(cleaned, reduction.use = "dm", dims.to.plot = 1:16, outer.title = "Diffusion Map (Sigma 19.538, 54 NNs): DPA", label="stage", plot.title="", legend=F)

#tsne with transitions
plotDim(cleaned, "time", transitions.plot = 10000, plot.title="DPA (with transitions)")

#use cells from wound healing as root
root.cells <- cellsInCluster(cleaned, "stage", "wound_healing")

#run 'flood' simulations
cleaned.floods <- floodPseudotime(cleaned, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

#process simulations into a pseudotime
cleaned <- floodPseudotimeProcess(cleaned, cleaned.floods, floods.name="pseudotime")

#check for adequate number of simulations (should reach asymptote)
pseudotimePlotStabilityOverall(cleaned)

#visualize tSNE with pseudotime overlaid
plotDim(cleaned, "pseudotime")

#plot pseudtime at each time point, Figure 7b
plotDists(cleaned, "pseudotime", "time", plot.title="Pseudotime by time point", legend = F)

#create URD object with just medium-bud cells
cleaned.23dpa <- urdSubset(cleaned, cells.keep=cellsInCluster(cleaned, "stage", "medium-bud"))

#use medium-bud variable genes
cleaned.23dpa@var.genes <- var.23dpa

# Calculate PCA and tSNE
cleaned.23dpa <- calcPCA(cleaned.23dpa, mp.factor = 1.5)
pcSDPlot(cleaned.23dpa)
set.seed(20)
cleaned.23dpa <- calcTsne(cleaned.23dpa)

#perform graph-based clustering
cleaned.23dpa <- graphClustering(cleaned.23dpa, num.nn = c(30, 40, 50, 70), method = 'Louvain', do.jaccard=T)
cleaned.23dpa <- graphClustering(cleaned.23dpa, num.nn = c(50, 60, 70, 80, 100), method = 'Infomap', do.jaccard=T)
clusterings <- c(paste0('Louvain-',c(30, 40, 50, 70)), paste0('Infomap-', c(50, 60, 70, 80, 100)))
for (c in clusterings) {plot(plotDim(cleaned.23dpa, c, legend=T))}

#tSNE plot Figure 7c
plotDim(cleaned.23dpa, 'Infomap-100', plot.title = "Medium-bud blastema populations", legend = F, point.size = 2)

#we chose to go forward with Infomap-100
clusters <- unique(cleaned.23dpa@group.ids$'Infomap-100')

#find marker genes, see Supplementary Data files for marker lists
pr.markers_23dpa_cleaned <- lapply(clusters, function(c) markersAUCPR(cleaned.23dpa, clust.1 = c, clustering = 'Infomap-100', genes.use= cleaned.23dpa@var.genes))
#inspect markers in R
head(pr.markers_23dpa_cleaned[[3]], 20)

#distal blastema population has some WE markers, let's remove WE cells
plotDot(cleaned.23dpa, genes = c('c1069858_g1_i1^sp|Q90X25|HXA13_CHICK^HoxA13_N','c1070920_g1_i4^sp|Q2VL56|PAX9_SAGOE^PAX^Tm2','c1020768_g1_i2^sp|Q9H2S6|TNMD_HUMAN^BRICHOS^Tm1','c1083312_g1_i2^sp|P70390|SHOX2_MOUSE','c1081900_g1_i4^sp|P25815|S100P_HUMAN','c1091168_g1_i2^sp|Q66S13|NATT4_THANI^DUF946^sigP'), clustering = 'Infomap-100')
distal.score <- apply(cleaned.23dpa@logupx.data[c('c1069858_g1_i1^sp|Q90X25|HXA13_CHICK^HoxA13_N','c1070920_g1_i4^sp|Q2VL56|PAX9_SAGOE^PAX^Tm2','c1020768_g1_i2^sp|Q9H2S6|TNMD_HUMAN^BRICHOS^Tm1','c1083312_g1_i2^sp|P70390|SHOX2_MOUSE'), cellsInCluster(cleaned.23dpa, 'Infomap-100', '1')], 2, sum.of.logs)
new.distal <- names(which(distal.score > 0))
remove.distal <- names(which(distal.score <= 0))

#add names for all pops
i100.n <- length(unique(cleaned.23dpa@group.ids$'Infomap-100'))
i100.cluster.assignments <- data.frame(clusters = 1:i100.n, name = rep(NA, i100.n), tip = rep(NA, i100.n), row.names = 1:i100.n)
i100.cluster.assignments['1', 'name'] <- "Distal Blastema"
i100.cluster.assignments['2', 'name'] <- "FAPs"
i100.cluster.assignments['3', 'name'] <- "Synovial Fibroblasts"
i100.cluster.assignments['4', 'name'] <- "Cartilage"
i100.cluster.assignments['5', 'name'] <- "Osteoblast-like"
i100.cluster.assignments['6', 'name'] <- "Joint"
i100.cluster.assignments['7', 'name'] <- "Schwann" 
i100.cluster.assignments['8', 'name'] <- "Endothelial"
i100.cluster.assignments['9', 'name'] <- "Myogenic"
i100.cluster.assignments['10', 'name'] <- "Pericytes" 


cluster.assignments <- i100.cluster.assignments

cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name), ]
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)


cleaned@group.ids$'23dpa-Infomap-100' <- NA
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Infomap-100'] <- cleaned.23dpa@group.ids$'Infomap-100'
cleaned@group.ids$'23dpa-Cluster' <- NA
cleaned.23dpa@group.ids$clusters.23dpa.name <- NA
cleaned.23dpa@group.ids$clusters.23dpa.num <- NA

#so now for cluster 1 I'll just give it the cleaned new.distal 
cleaned.23dpa@group.ids[new.distal, "clusters.23dpa.name"] <- cluster.assignments[1, "name"] 
cleaned.23dpa@group.ids[new.distal, "clusters.23dpa.num"] <- as.character(cluster.assignments[1, "cluster.new"])

cells_2 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '2') 
cleaned.23dpa@group.ids[cells_2, "clusters.23dpa.name"] <- cluster.assignments[2, "name"] 
cleaned.23dpa@group.ids[cells_2, "clusters.23dpa.num"] <- as.character(cluster.assignments[2, "cluster.new"])

cells_3 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '3') 
cleaned.23dpa@group.ids[cells_3, "clusters.23dpa.name"] <- cluster.assignments[3, "name"] 
cleaned.23dpa@group.ids[cells_3, "clusters.23dpa.num"] <- as.character(cluster.assignments[3, "cluster.new"])

cells_4 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '4') 
cleaned.23dpa@group.ids[cells_4, "clusters.23dpa.name"] <- cluster.assignments[4, "name"] 
cleaned.23dpa@group.ids[cells_4, "clusters.23dpa.num"] <- as.character(cluster.assignments[4, "cluster.new"])

cells_5 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '5') 
cleaned.23dpa@group.ids[cells_5, "clusters.23dpa.name"] <- cluster.assignments[5, "name"] 
cleaned.23dpa@group.ids[cells_5, "clusters.23dpa.num"] <- as.character(cluster.assignments[5, "cluster.new"])

cells_6 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '6') 
cleaned.23dpa@group.ids[cells_6, "clusters.23dpa.name"] <- cluster.assignments[6, "name"] 
cleaned.23dpa@group.ids[cells_6, "clusters.23dpa.num"] <- as.character(cluster.assignments[6, "cluster.new"])

cells_7 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '7') 
cleaned.23dpa@group.ids[cells_7, "clusters.23dpa.name"] <- cluster.assignments[7, "name"] 
cleaned.23dpa@group.ids[cells_7, "clusters.23dpa.num"] <- as.character(cluster.assignments[7, "cluster.new"])

cells_8 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '8') 
cleaned.23dpa@group.ids[cells_8, "clusters.23dpa.name"] <- cluster.assignments[8, "name"] 
cleaned.23dpa@group.ids[cells_8, "clusters.23dpa.num"] <- as.character(cluster.assignments[8, "cluster.new"])

cells_9 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '9') 
cleaned.23dpa@group.ids[cells_9, "clusters.23dpa.name"] <- cluster.assignments[9, "name"] 
cleaned.23dpa@group.ids[cells_9, "clusters.23dpa.num"] <- as.character(cluster.assignments[9, "cluster.new"])

cells_10 <- cellsInCluster(cleaned.23dpa, clustering = 'Infomap-100', cluster = '10') 
cleaned.23dpa@group.ids[cells_10, "clusters.23dpa.name"] <- cluster.assignments[10, "name"] 
cleaned.23dpa@group.ids[cells_10, "clusters.23dpa.num"] <- as.character(cluster.assignments[10, "cluster.new"])


cleaned@group.ids$'23dpa-Infomap-100' <- NA
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Infomap-100'] <- cleaned.23dpa@group.ids$'Infomap-100'
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Cluster'] <- cleaned.23dpa@group.ids$clusters.23dpa.name
cleaned@group.ids$'23dpa-Cluster-Num'<- NA
cleaned@group.ids[rownames(cleaned.23dpa@group.ids), '23dpa-Cluster-Num'] <- cleaned.23dpa@group.ids$clusters.23dpa.num


cleaned@group.ids[rownames(cleaned.23dpa@group.ids), "tip.clusters"] <- cleaned.23dpa@group.ids$clusters.23dpa.num

#determine potential for terminal populations
potential <- clusterTipPotential(cleaned, 'pseudotime', 'tip.clusters', name.store = 'tip.potential')
potential 

#tips = 3, 4, 5, 6, 7, 8, 9, 10, so everything but "distal blastema" and "FAPs"

only.tips <- cellsInCluster(cleaned, clustering= '23dpa-Cluster-Num', cluster = c('3','4','5','6','7','8','9','10'))
tips <- urdSubset(cleaned, cells.keep= only.tips)
cleaned@group.ids[rownames(tips@group.ids), "real.tip.clusters"] <- tips@group.ids$'23dpa-Cluster-Num'

#determine parameters of the logistic used to bias the transition probabilities
cleaned.ptlogistic <- pseudotimeDetermineLogistic(cleaned, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)

#bias the transition matrix acording to pseudotime
cleaned.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(cleaned, "pseudotime", logistic.params=cleaned.ptlogistic))

#simulate the biased random walks from each tip
cleaned.walks <- simulateRandomWalksFromTips(cleaned, tip.group.id="real.tip.clusters", root.cells=root.cells, transition.matrix = cleaned.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)

#process the biased random walks into visitation frequencies
cleaned <- processRandomWalksFromTips(cleaned, cleaned.walks, verbose = F)

#color only tip clusters on tSNE
plotDim(cleaned, "real.tip.clusters", plot.title="Cells in each tip")

#load tip clusters into tree
cleaned.tree <- loadTipCells(cleaned, "real.tip.clusters")

#build tree
cleaned.tree <- buildTree(cleaned.tree, pseudotime = "pseudotime", tips.use=c('3','4','5','6','7','8','9','10'), divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

#rename clusters
cleaned.tree <- nameSegments(cleaned.tree, segments=c('3','4','5','6','7','8','9','10'), segment.names = c("Synovial Fibros", "Cartilage", "Osteoblast-like", "Joint","Schwann", "Endothelial", "Myogenic", "Pericyte"), short.names = c("Synovial Fibros", "Cartilage",  "Osteoblast-like", "Joint","Schwann", "Endothelial", "Myogenic", "Pericyte"))

#plot tree with time point info overlaid
plotTree(cleaned.tree, "stage", title="DPA")

#plot tree with medium-bud blastema colors overlaid
plotTree(cleaned.tree, 'tip.clusters', title = '23dpa_clusters', cell.alpha = 0.5, cell.size = 2, tree.alpha = 0.5, tree.size = .25)

joint.markers <- aucprTestAlongTree(cleaned.tree, pseudotime="pseudotime", tips='Joint', log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, only.return.global=F, must.beat.sibs=0.6, report.debug=T)

synovial.markers <- aucprTestAlongTree(cleaned.tree, pseudotime="pseudotime", tips='Synovial Fibros', log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, only.return.global=F, must.beat.sibs=0.6, report.debug=T)

#visuzalize population specific markers
#osteoblast-like

p1<- plotTree(cleaned.tree, 'c862122_g2_i1^sp|Q6DJ00|OSTCN_XENTR^sp|P40147|OSTCN_XENLA^Gla^sigP', title = 'OSTCN')
p2 <- plotTree(cleaned.tree, 'c1034953_g1_i2^sp|A1YQ92|ODAM_MACMU', title = 'ODAM')

#example of graph in Supplementary Fig7
plot_grid(p1,p2)

#can continue to overlay markers like shown above

#save
save.image('blastema.URD.RData')

#quit
q()


##Pseudotime of epidermis with Monocle

#homeostatic epidermis
library(monocle)
library(Seurat)

load("intact.Rdata")

#add IDs to metadata in preparation of importCDS
merged@meta.data$stashed.id <- merged@ident

pdf('intact_timecourse_tsne.pdf')
tsne.plot(merged, do.label = T)
dev.off()

WE_filtered <- SubsetData(merged, ident.remove = c(1,2,4,5,6,8,11,12))

pdf('intact_WE_timecourse_tsne.pdf')
tsne.plot(WE_filtered, do.label = T)
dev.off()

WE_filtered <- importCDS(WE_filtered)

WE_filtered <- detectGenes(WE_filtered, min_expr = .1)
expressed_genes <- row.names(subset(fData(WE_filtered),
                                    num_cells_expressed >= 3))

WE_filtered <- WE_filtered[expressed_genes,]

head(pData(WE_filtered))

WE_filtered <- estimateSizeFactors(WE_filtered)
WE_filtered <- detectGenes(WE_filtered, min_expr = .1)


head(fData(WE_filtered))
head(pData(WE_filtered))

pdf('intact_WE_variance.pdf')
plot_pc_variance_explained(WE_filtered, return_all = FALSE)
dev.off()

WE_filtered <- reduceDimension(WE_filtered,
                               max_components = 2,
                               norm_method = 'log',
                               num_dim = 8,
                               cores = 4,
                               residualModelFormulaStr = "~num_genes_expressed",
                               reduction_method = 'tSNE',
                               verbose = TRUE)

WE_filtered <- clusterCells(WE_filtered, verbose = TRUE)

WE_filtered$tree.ident <- as.character(WE_filtered$tree.ident)

pdf('intact_WE_clusters.pdf')
plot_cell_clusters(WE_filtered)
dev.off()

pdf('intact_WE_numGenes.pdf')
plot_cell_clusters(WE_filtered, color_by = "num_genes_expressed")
dev.off()

pdf('intact_WE_byseurat.pdf')
plot_cell_clusters(WE_filtered, color_by = "stashed.id")
dev.off()


pdf('intact_WE_rho_delta.pdf')
plot_rho_delta(WE_filtered, rho_threshold = 2, delta_threshold = 10)
dev.off()

clustering_DEG_genes <- differentialGeneTest(WE_filtered,
                                             fullModelFormulaStr = '~stashed.id',
                                             cores = 4)

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]
WE_filtered <- setOrderingFilter(WE_filtered, ordering_genes = my_ordering_genes)
WE_filtered <- reduceDimension(WE_filtered,
                               method = 'DDRTree',
                               fullModelFormulaStr = '~stashed.id', verbose = F, scaling = T, maxIter = 100, norm_method = 'log', max_components = 6, param.gamma = 20)

WE_filtered <- orderCells(WE_filtered)

pdf('intact_WE_pseudotime.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_byCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Cluster", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster_facet.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id") +
  facet_wrap(~stashed.id, nrow = 1)
dev.off()

save.image("082018_intact_WE_postDEG.Rdata")

###

load("082018_intact_WE_postDEG.Rdata")

# create vector of no's
my_vector <- rep('no', nrow(pData(WE_filtered)))

# change status to yes if the cell was in cluster 10
my_vector[pData(WE_filtered)$stashed.id == 10] <- rep('yes', sum(pData(WE_filtered)$stashed.id == 10))

# add vector to phenoData
pData(WE_filtered)$test <- my_vector
head(pData(WE_filtered))

clustering_DEG_genes <- differentialGeneTest(WE_filtered,
                                             fullModelFormulaStr = '~test',
                                             cores = 4)

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:150]
WE_filtered <- setOrderingFilter(WE_filtered, ordering_genes = my_ordering_genes)
WE_filtered <- reduceDimension(WE_filtered,
                               method = 'DDRTree',
                               fullModelFormulaStr = '~stashed.id', verbose = F, scaling = T, maxIter = 100, norm_method = 'log', max_components = 15, param.gamma = 20)

WE_filtered <- orderCells(WE_filtered, reverse = TRUE)

pdf('intact_WE_pseudotime.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_byCluster.pdf')
plot_cell_trajectory(WE_filtered, color_by = "Cluster", show_branch_points = FALSE)
dev.off()

pdf('intact_WE_bySeuratCluster_facet.pdf')
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE) +
  facet_wrap(~stashed.id, nrow = 1)
dev.off()

###

save.image("intact_WE_postDEG.Rdata")

rm(list=setdiff(ls(), "WE_filtered")) 

save.image("intact_WE_slim.Rdata")

png('intact_WE_facet.png', width = (7/4)*3, height = 2, units = 'in', res = 300)
plot_cell_trajectory(WE_filtered, color_by = "stashed.id", show_branch_points = FALSE) +
  facet_wrap(~stashed.id, nrow = 1)+
  theme(legend.position="none")
dev.off()

png('intact_WE_pseudotime.png', width = 2.1, height = 1.77, units = 'in', res = 300)
plot_cell_trajectory(WE_filtered, color_by = "Pseudotime", show_branch_points = FALSE)+
  theme(legend.position="none")
dev.off()

###Regenerating epidermis


library(monocle)
library(Seurat)

#load medium-bud blastema data (medium-bud and 23dpa are equivalent)
load("medium_bud_blastema.RData")

pdf('timecourse_tsne.pdf')
tsne.plot(combined, do.label = T)
dev.off()

combined <- importCDS(combined)

pData(combined)$"day1" <- pData(combined)$protocol == "d1" 

combined <- detectGenes(combined, min_expr = .1)
expressed_genes <- row.names(subset(fData(combined),
                                    num_cells_expressed >= 3))

#only want expressed genes and the day 1 subset of cells
combined_day1 <- combined[expressed_genes, pData(combined)$"day1"]

#only want WE populations
pData(combined_day1)$"WE" <- pData(combined_day1)$tree.ident == "1" |
  pData(combined_day1)$tree.ident == "2" |
  pData(combined_day1)$tree.ident == "3" |
  pData(combined_day1)$tree.ident == "4"

combined_day1_WE <- combined_day1[expressed_genes, pData(combined_day1)$"WE"]

head(pData(combined_day1_WE))

combined_day1_WE <- estimateSizeFactors(combined_day1_WE)
combined_day1_WE <- detectGenes(combined_day1_WE, min_expr = .1)


head(fData(combined_day1_WE))
head(pData(combined_day1_WE))

pdf('23dpa_1_WE_variance.pdf')
plot_pc_variance_explained(combined_day1_WE, return_all = FALSE)
dev.off()

combined_day1_WE <- reduceDimension(combined_day1_WE,
                                    max_components = 2,
                                    norm_method = 'log',
                                    num_dim = 8,
                                    cores = 4,
                                    residualModelFormulaStr = "~num_genes_expressed",
                                    reduction_method = 'tSNE',
                                    verbose = TRUE)

combined_day1_WE <- clusterCells(combined_day1_WE, verbose = TRUE)

#plot clusters
pdf('23dpa_day1_WE_clusters.pdf')
plot_cell_clusters(combined_day1_WE)
dev.off()

pdf('23dpa_day1_WE_numGenes.pdf')
plot_cell_clusters(combined_day1_WE, color_by = "num_genes_expressed")
dev.off()

pdf('23dpa_day1_WE_byseurat.pdf')
plot_cell_clusters(combined_day1_WE, color_by = "tree.ident")
dev.off()


pdf('WE_rho_delta.pdf')
plot_rho_delta(combined_day1_WE, rho_threshold = 2, delta_threshold = 10)
dev.off()

clustering_DEG_genes <- differentialGeneTest(combined_day1_WE,
                                             fullModelFormulaStr = '~tree.ident',
                                             cores = 4)

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:150]
combined_day1_WE <- setOrderingFilter(combined_day1_WE, ordering_genes = my_ordering_genes)
combined_day1_WE <- reduceDimension(combined_day1_WE,
                                    method = 'DDRTree',
                                    fullModelFormulaStr = '~tree.ident', verbose = F, scaling = T, norm_method = 'log', max_components =8, param.gamma = 20)

combined_day1_WE$tree.ident <- as.character(combined_day1_WE$tree.ident)

combined_day1_WE <- orderCells(combined_day1_WE)

#plot pseudotime graphs
pdf('23dpa_day1_WE_pseudotime.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

pdf('23dpa_day1_WE_bySeuratCluster.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "tree.ident", show_branch_points = FALSE)
dev.off()

pdf('23dpa_day1_WE_byState.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "Cluster", show_branch_points = FALSE)
dev.off()

pdf('23dpa_day1_WE_bySeuratCluster_facet.pdf')
plot_cell_trajectory(combined_day1_WE, color_by = "tree.ident", show_branch_points = FALSE) +
  facet_wrap(~tree.ident, nrow = 1)
dev.off()

save.image("23dpa_day1_WE_postDEG.Rdata")

rm(list=setdiff(ls(), "combined_day1_WE")) 

save.image("23dpa_day1_WE_slim.Rdata")

png('23dpa_day1_WE_facet.png', width = (7/5)*3, height = 2, units = 'in', res = 300)
plot_cell_trajectory(combined_day1_WE, color_by = "tree.ident", show_branch_points = FALSE) +
  facet_wrap(~tree.ident, nrow = 1)+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(-6,-3,0,2))
dev.off()

png('23dpa_day1_WE_pseudotime.png', width = (7/5)*1.15, height = 1.70, units = 'in', res = 300)
plot_cell_trajectory(combined_day1_WE, color_by = "Pseudotime", show_branch_points = FALSE) +
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(-6,-3,0,2))
dev.off()

png('23dpa_day1_WE_pseudotime_K17.png', width = 7, height = 2, units = 'in', res = 300)
genes_to_plot <- c('c1083200_g2_i3^sp|A1L595|K1C17_BOVIN^Filament')
plot_genes_in_pseudotime(combined_day1_WE[genes_to_plot,], color_by = "tree.ident")
dev.off()



