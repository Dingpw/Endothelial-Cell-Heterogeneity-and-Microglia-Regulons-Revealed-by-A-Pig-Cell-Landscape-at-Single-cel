################################
########-snRNA pipeline-########
################################
library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(getopt)
command=matrix(c("tissuelib","t",1,"character",
                 "tissue","T",1,"character",
                 "tissue_abb","a",1,"character",
                 "input_dir1","i",1,"character" ),byrow=T,ncol=4)
args=getopt(command)
tissuelib=args$tissuelib
tissue=args$tissue
tissue_abb=args$tissue_abb
input_dir1=args$input_dir1
tissuelib<-unlist(strsplit(tissuelib,"; "))

ribosome_genes<-read.table("/ref/pig_ribosome_gene.txt",header = F,sep = "\t",stringsAsFactors = F)
exc_genes<-read.table("/ref/pig_extracellular.txt",header = F,sep = "\t",stringsAsFactors = F)
cyto<-read.table("/ref/pig_cytoplasm.txt",header = F,sep = "\t",stringsAsFactors = F)
memb<-read.table("/ref/pig_membrane.txt",header = F,sep = "\t",stringsAsFactors = F)
lnc_gene<-read.table("/ref/Pig_lncGene.txt",header = T,sep = "\t",stringsAsFactors = F)
protein_coding_gene<-read.table("/ref/pig_prcGene.txt",header = T,sep = "\t",stringsAsFactors = F)
hk_gene<-c("B2M","HMBS","HPRT1")
tf_gene<-read.table("/ref/Sus_scrofa_TF.txt",header = T,sep = "\t",stringsAsFactors = F)

dir.create(tissue)
setwd(tissue)
sample.list<-read.table(input_dir1,header = F,stringsAsFactors = F)
Samples<-as.character(sample.list[,1])
#Samples <- c("PG-3_54","PG-4_62","PG-5_69","PG-LR-1")
i=0
idrop_LR<-list()
for (s in Samples) {
  i=i+1
  input_dir<-paste(s, "_outs/count_mtx.tsv",sep = "")
  seurat.data <- read.table(input_dir, header = T, sep="\t")
  rownames(seurat.data) <- seurat.data[,1]
  seurat.data <- seurat.data[,2:dim(seurat.data)[2]]
  
  seurat_data <- CreateSeuratObject(counts = seurat.data, 
                                    project = paste(tissue_abb,i,sep = ""), 
                                    min.cells = 3, 
                                    min.features = 0)
  idrop_LR[[i]]<-seurat_data
}
i<-length(idrop_LR)
#m_i=1
idrop_LR_object<-merge(idrop_LR[[1]],y=idrop_LR[[2]],add.cell.ids = paste(c("idrop_1", "idrop_2"),tissue_abb,sep = "_"), project = paste("idrop",tissue_abb,sep = ""))
m_i<-2
while(m_i<i){
  m_i<-m_i+1
  idrop_LR_object<-merge(idrop_LR_object,y=idrop_LR[[m_i]],add.cell.ids = c("",paste("idrop",m_i,tissue_abb,sep = "_")), project = paste("idrop",tissue_abb,sep = ""))
}
#features = c("ND2","COX1","COX2","ATP8", "ATP6","COX3","ND3","ND4L", "ND4","ND5", "ND6","CYTB")
idrop_LR_object[["percent.mt"]] <- PercentageFeatureSet(idrop_LR_object, pattern = "^MT-")

idrop_LR_object[["percent.rib"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),ribosome_genes$V2))
idrop_LR_object[["percent.exc"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),exc_genes$V2))
idrop_LR_object[["percent.cyt"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),cyto$V2))

idrop_LR_object[["percent.mem"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),memb$V2))
idrop_LR_object[["percent.pcg"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),protein_coding_gene$Gene.name))

idrop_LR_object[["percent.tf"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                        features=intersect(rownames(idrop_LR_object),c(tf_gene$Symbol,tf_gene$Ensembl)))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

idrop_LR_object <- CellCycleScoring(idrop_LR_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#filtering
idrop_LR_object<-subset(idrop_LR_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
#clustering
idrop_LR_object <- NormalizeData(idrop_LR_object, normalization.method = "LogNormalize", scale.factor = 10000)
idrop_LR_object <- FindVariableFeatures(idrop_LR_object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(idrop_LR_object)
idrop_LR_object <- ScaleData(idrop_LR_object, features = all.genes)
idrop_LR_object <- RunPCA(idrop_LR_object, features = VariableFeatures(object = idrop_LR_object))
pdf("PCADimPlot.pdf")
DimPlot(object = idrop_LR_object, group.by="orig.ident",
        reduction = "pca")
dev.off()
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
idrop_LR_object <- JackStraw(idrop_LR_object, num.replicate = 100)
idrop_LR_object <- ScoreJackStraw(idrop_LR_object, dims = 1:20)
n_PC<-max(which(idrop_LR_object@reductions[["pca"]]@jackstraw@overall.p.values[,2]<0.01))
n_PC

idrop_LR_object <- RunUMAP(idrop_LR_object, reduction = "pca", dims = 1:n_PC)
idrop_LR_object<-RunTSNE(idrop_LR_object, reduction = "pca", dims = 1:n_PC)

p1 <- DimPlot(idrop_LR_object, reduction = "umap", group.by = "orig.ident")
pdf("idrop_umap.pdf",height = 8,width = 10)
print(p1)
dev.off()

p1 <- DimPlot(idrop_LR_object, reduction = "tsne", group.by = "orig.ident")
pdf("idrop_tsne.pdf",height = 8,width = 10)
print(p1)
dev.off()

idrop_LR_object <- FindNeighbors(idrop_LR_object, dims = 1:n_PC)
idrop_LR_object <- FindClusters(idrop_LR_object, resolution = 1.0)

p1 <- DimPlot(idrop_LR_object, reduction = "umap", group.by = "orig.ident",label = TRUE)
p2 <- DimPlot(idrop_LR_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, 
              repel = TRUE)
p<-p1 + p2
pdf("clustering_umap.pdf",height = 8,width = 20)
print(p)
dev.off()
saveRDS(idrop_LR_object,paste(tissue_abb,"idrop_object.rds",sep = ""))

################################
########-scRNA pipeline-########
################################
library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(getopt)
command=matrix(c("tissuelib","t",1,"character",
                 "tissue","T",1,"character",
                 "tissue_abb","a",1,"character",
                 "input_dir1","i",1,"character" ),byrow=T,ncol=4)
args=getopt(command)
tissue=args$tissue
tissue_abb=args$tissue_abb
input_dir1=args$input_dir1
ribosome_genes<-read.table("/ref/pig_ribosome_gene.txt",header = F,sep = "\t",stringsAsFactors = F)
exc_genes<-read.table("/ref/pig_extracellular.txt",header = F,sep = "\t",stringsAsFactors = F)
cyto<-read.table("/ref/pig_cytoplasm.txt",header = F,sep = "\t",stringsAsFactors = F)
memb<-read.table("/ref/pig_membrane.txt",header = F,sep = "\t",stringsAsFactors = F)
lnc_gene<-read.table("/ref/Pig_lncGene.txt",header = T,sep = "\t",stringsAsFactors = F)
protein_coding_gene<-read.table("/ref/pig_prcGene.txt",header = T,sep = "\t",stringsAsFactors = F)
hk_gene<-c("B2M","HMBS","HPRT1")
tf_gene<-read.table("/ref/Sus_scrofa_TF.txt",header = T,sep = "\t",stringsAsFactors = F)

dir.create(tissue)
setwd(tissue)
seurat.data<-Read10X(data.dir = input_dir1)  
seurat_data <- CreateSeuratObject(counts = seurat.data, 
                                  project = tissue_abb, 
                                  min.cells = 3, 
                                  min.features = 0)
idrop_LR_object<-seurat_data
#features = c("ND2","COX1","COX2","ATP8", "ATP6","COX3","ND3","ND4L", "ND4","ND5", "ND6","CYTB")
idrop_LR_object[["percent.mt"]] <- PercentageFeatureSet(idrop_LR_object, pattern = "^MT-")

idrop_LR_object[["percent.rib"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),ribosome_genes$V2))
idrop_LR_object[["percent.exc"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),exc_genes$V2))
idrop_LR_object[["percent.cyt"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),cyto$V2))

idrop_LR_object[["percent.mem"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),memb$V2))
idrop_LR_object[["percent.pcg"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                         features=intersect(rownames(idrop_LR_object),protein_coding_gene$Gene.name))

idrop_LR_object[["percent.tf"]] <- PercentageFeatureSet(idrop_LR_object, 
                                                        features=intersect(rownames(idrop_LR_object),c(tf_gene$Symbol,tf_gene$Ensembl)))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

idrop_LR_object <- CellCycleScoring(idrop_LR_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#filtering
idrop_LR_object<-subset(idrop_LR_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
#clustering
idrop_LR_object <- NormalizeData(idrop_LR_object, normalization.method = "LogNormalize", scale.factor = 10000)
idrop_LR_object <- FindVariableFeatures(idrop_LR_object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(idrop_LR_object)
idrop_LR_object <- ScaleData(idrop_LR_object, features = all.genes)
idrop_LR_object <- RunPCA(idrop_LR_object, features = VariableFeatures(object = idrop_LR_object))
pdf("PCADimPlot.pdf")
DimPlot(object = idrop_LR_object, group.by="orig.ident",
        reduction = "pca")
dev.off()
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
idrop_LR_object <- JackStraw(idrop_LR_object, num.replicate = 100)
idrop_LR_object <- ScoreJackStraw(idrop_LR_object, dims = 1:20)
n_PC<-max(which(idrop_LR_object@reductions[["pca"]]@jackstraw@overall.p.values[,2]<0.01))
n_PC

idrop_LR_object <- RunUMAP(idrop_LR_object, reduction = "pca", dims = 1:n_PC)
idrop_LR_object<-RunTSNE(idrop_LR_object, reduction = "pca", dims = 1:n_PC)

p1 <- DimPlot(idrop_LR_object, reduction = "umap", group.by = "orig.ident")
pdf("10x_umap.pdf",height = 8,width = 10)
print(p1)
dev.off()

p1 <- DimPlot(idrop_LR_object, reduction = "tsne", group.by = "orig.ident")
pdf("10x_tsne.pdf",height = 8,width = 10)
print(p1)
dev.off()

idrop_LR_object <- FindNeighbors(idrop_LR_object, dims = 1:n_PC)
idrop_LR_object <- FindClusters(idrop_LR_object, resolution = 1.0)

p1 <- DimPlot(idrop_LR_object, reduction = "umap", group.by = "orig.ident",label = TRUE)
p2 <- DimPlot(idrop_LR_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, 
              repel = TRUE)
p<-p1 + p2
pdf("clustering_umap.pdf",height = 8,width = 20)
print(p)
dev.off()
saveRDS(idrop_LR_object,paste(tissue_abb,"10x_object.rds",sep = ""))

#========hclust===================
sample.cluster<-AverageExpression(idrop_LR_object)$RNA
hc = hclust(dist(t(sample.cluster)))
hcd = as.dendrogram(hc)
pdf("hclust_integrated.pdf")
plot(hcd)
dev.off()
#-----------Identify conserved cell type markers
DefaultAssay(object = idrop_LR_object) <- "RNA"

markers <- FindAllMarkers(object = idrop_LR_object, 
                          only.pos = T, 
                          min.pct = 0.5, 
                          logfc.threshold = 0.5)
write.table(markers,file="markers_0.5.txt",sep="\t",quote=F)

################################
####-Fig.1: merge all object-###
################################
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
files<-read.table("all_rds_20211024.txt",header=F,sep="\t")
files<-files[,c(2,3)]
object_all<-readRDS(as.character(files[1,2]))
object_all$Tissue<-as.character(files[1,1])
for(i in 2:nrow(files)){
  object<-readRDS(as.character(files[i,2]))
  object$Tissue<-as.character(files[i,1])
  object_all<-merge(x=object_all,y=object,project = "Pig")
}
object_all <- NormalizeData(object = object_all, normalization.method = "LogNormalize")
object_all<- FindVariableFeatures(object = object_all, mean.function = ExpMean, dispersion.function = LogVMR)
object_all <- ScaleData(object_all, verbose = FALSE)
object_all <- RunPCA(object_all, npcs=20, verbose=FALSE)
object_all <- RunUMAP(object_all, reduction = "pca", dims = 1:20)
object_all <- RunTSNE(object_all, reduction = "pca", dims = 1:20)
saveRDS(object_all,"merge_all.rds")

################################
######-Fig.S1: SCTransform-#####
################################
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(pheatmap)
ten<-readRDS("10x_spleen.rds")
idrop<-readRDS("idrop_spleen.rds")

list<-list(ten,idrop)
for (i in 1:length(list)) {
  list[[i]] <- SCTransform(list[[i]], verbose = FALSE)
}
features<-SelectIntegrationFeatures(object.list = list,nfeatures = 3000)
list<-PrepSCTIntegration(object.list=list,anchor.features=features)
anchors<-FindIntegrationAnchors(object.list = list,normalization.method = "SCT",anchor.features = features)
spleen<-IntegrateData(anchorset = anchors,normalization.method = "SCT")
spleen<-RunPCA(spleen)
spleen <- JackStraw(spleen, num.replicate = 100)
spleen <- ScoreJackStraw(spleen, dims = 1:20)
n_PC<-max(which(spleen@reductions[["pca"]]@jackstraw@overall.p.values[,2]<0.01))
n_PC
spleen<-RunUMAP(spleen,dims = 1:n_PC)
spleen<-RunTSNE(spleen,dims = 1:n_PC)

spleen <- FindNeighbors(spleen, dims = 1:n_PC)
spleen <- FindClusters(spleen, resolution = 1.0)
saveRDS(spleen,"integrated_spleen.sct.rds")

dir.create("plots_spleen")
setwd("plots_spleen/")

spleen$Platform<-as.character(spleen$orig.ident)
spleen$Platform[which(spleen$orig.ident %in% "spleen")]<-"10X"
spleen$Platform[which(spleen$orig.ident %in% c("SN1","SN2","SN3","SN4"))]<-"iDrop"
### plots
stat<-as.data.frame(table(spleen$Celltype,spleen$Platform))
stat$num <- ifelse(stat$Var2 == "10X", 
                   stat$Freq *1,
                   stat$Freq * -1)
pdf("stat_number.pdf")
ggplot(data = stat) + 
  geom_col(aes(x = factor(Var1), y = num,fill = Var2)) + 
  coord_flip() + theme_bw()
dev.off()

for (i in unique(stat$Var2)) {
  stat1<-stat[which(stat$Var2 == i),]
  label_value <- paste('(', round(stat1$Freq/sum(stat1$Freq) * 100, 1), '%)', sep = '')
  label <- paste(stat1$Var1, label_value, sep = '')
  p<-ggplot(data = stat1, mapping = aes(x = "Number", y = Freq, fill = Var1)) + 
    geom_bar(stat = 'identity', position = 'stack', width = 0.4)+
    coord_polar(theta = "y")+
    labs(x = '', y = '', title = '')+
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    scale_fill_discrete(labels = label)
  
  pdf(paste0("stat_pie_",i,".pdf"),width = 8,height = 6)
  print(p)
  dev.off()
}

write.table(stat,"stat.txt",quote = F,sep = "\t",row.names = F)

########Variable gene
variable_genes<-spleen@assays$integrated@var.features[1:3000]
Idents(spleen)<-"Celltype"
cluster.plat<-AverageExpression(spleen,add.ident = "Platform",features = variable_genes,assays = "integrated",slot = "scale.data")
ave<-cluster.plat$integrated
ten.ave<-ave[,grep("10X",colnames(ave))]
idrop.ave<-ave[,grep("iDrop",colnames(ave))]
genes<-variable_genes
ten.ave_1 <- ten.ave[genes,]
idrop.ave_1 <- idrop.ave[genes,]

R1 <- data.frame(matrix(NA, nrow = dim(ten.ave_1)[2], ncol = dim(idrop.ave_1)[2]))
P1 <- data.frame(matrix(NA, nrow = dim(ten.ave_1)[2], ncol = dim(idrop.ave_1)[2]))

colnames(R1) <- paste0("idrop_",colnames(idrop.ave_1))
rownames(R1) <- paste0("tenx_",colnames(ten.ave_1))
colnames(P1) <- paste0("idrop_",colnames(idrop.ave_1))
rownames(P1) <- paste0("tenx_",colnames(ten.ave_1))

for (i in colnames(ten.ave_1)) {
  tmp1 <- ten.ave_1[,i]
  for (j in colnames(idrop.ave_1)) {
    tmp2 <- idrop.ave_1[,j]
    MyR <- cor.test(tmp1,tmp2,method = "pearson")
    R1[paste0("tenx_",i),paste0("idrop_",j)] <- as.numeric(MyR$estimate)
    P1[paste0("tenx_",i),paste0("idrop_",j)] <- as.numeric(MyR$p.value)
  }
}

write.csv(R1, "tenx_vs_idrop_Pearson_Cor_integrated_variable_scaled.csv", quote = F, row.names = T)
write.csv(P1, "tenx_vs_idrop_Pearson_Pva_integrated_variable_scaled.csv", quote = F, row.names = T)

R1<-R1[,c(6,4,3,5,2,1)]
R1<-R1[c(3,9,2,5,4,8,7,6,1),]
pdf("Corr_VariableGene_all_integrated.pdf")
pheatmap(R1,display_numbers = T,cluster_rows = F,cluster_cols = F,drop_levels = F)
dev.off()

################################
#####-Fig.2: GO enrichment-#####
################################
retina_GO<-read.table("Retina_Human_GO.txt",header = T,sep = "\t",quote = "")
#### Retina
selected_go<-c(## Bip
  "axonogenesis","axonal transport","axo-dendritic transport","anterograde axonal transport",#"regulation of neuron projection development"
  ## Cone
  "photoreceptor cell cilium","9+0 non-motile cilium","photoreceptor outer segment",
  ## End
  "endothelium development","endothelial cell differentiation","endothelial cell development","vasculogenesis","endothelial cell migration","blood vessel endothelial cell migration",
  ## Microglia
  "macrophage differentiation","macrophage activation","microglial cell activation","glial cell activation",
  #"regulation of neuron death",
  ## Muller
  "negative regulation of neurogenesis","negative regulation of nervous system development","negative regulation of cell development","gliogenesis","neural retina development",#"negative regulation of neuron differentiation"
  ## RGC
  "regulation of neuron projection development","negative regulation of neuron differentiation","negative regulation of neuron projection development","regulation of axonogenesis","neuron recognition",
  ## Rod
  "detection of light stimulus","phototransduction","rhodopsin mediated signaling pathway","phototransduction, visible light","regulation of rhodopsin mediated signaling pathway","photoreceptor cell differentiation",
  ## T
  "T cell activation","T cell differentiation","lymphocyte differentiation","T cell proliferation"
)

selected<-retina_GO[which(retina_GO$Description %in% selected_go),]
selected$Description<-factor(selected$Description,levels = selected_go)
selected$Cluster<-factor(selected$Cluster,levels = rev(levels(selected$Cluster)))
p <- ggplot(data=selected,mapping=aes_string(x="Cluster",y="Description"))+
  geom_point(mapping=aes_string(size="Count",color="p.adjust"),show.legend = TRUE)+
  scale_color_gradient(low="red",high = "blue")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "Count")) +
  labs(x = "", y = "",color="p.adjust",size="Count")+
  theme_bw()+
  coord_flip()+
  theme(strip.text.x = element_text(angle=0,face="bold",size = 16),strip.text.y=element_text(angle=90,size = 16))+
  theme(axis.text.x = element_text(face = "bold",colour = "black", angle = 90, hjust = 1),axis.text.y = element_text(colour = "black"))
print(p)


kidney_GO<-read.table("Kidney_Human_GO.txt",header = T,sep = "\t",quote = "")
#### kidney
selected_go<-c(## Collec
  "regulation of peptidase activity","response to metal ion","sodium ion homeostasis","cellular response to drug","sodium ion transmembrane transport",
  ## DCT
  "positive regulation of sodium ion","positive regulation of sodium ion transmembrane transport","positive regulation of sodium ion transport","regulation of membrane potential","regulation of sodium ion transmembrane transporter activity",
  ##End
  "endothelial cell proliferation","regulation of endothelial cell proliferation","cell-matrix adhesion","endothelium development","regulation of angiogenesis",
  ##Fib
  "ficolin-1-rich granule lumen","ficolin-1-rich granule","myofibril","contractile fiber","myofilament",
  ##Loop
  "potassium ion homeostasis","potassium ion import across plasma membrane","potassium ion import","inorganic cation import across plasma membrane","cellular potassium ion homeostasis","chloride ion homeostasis","metanephric nephron tubule development","monovalent inorganic anion homeostasis",
  "anion:sodium symporter activity","cation:chloride symporter activity","anion:cation symporter activity",
  ## Podo
  "glomerulus development","glomerular basement membrane development","nephron development","renal system development","urogenital system development",
  "renal filtration cell differentiation","glomerular epithelial cell differentiation",
  ##PTcell
  "anion transmembrane transport","endocytic vesicle","solute:sodium symporter activity","solute:cation symporter activity","anion transmembrane transporter activity",
  ##T
  ## Ure
  "cellular monovalent inorganic cation homeostasis","monovalent inorganic cation homeostasis","excretion","negative regulation of proteolysis","proton transmembrane transport"
)

selected<-kidney_GO[which(kidney_GO$Description %in% selected_go),]
selected$Description<-factor(selected$Description,levels = selected_go)
selected$Cluster<-factor(selected$Cluster,levels = rev(levels(selected$Cluster)))
p <- ggplot(data=selected,mapping=aes_string(x="Cluster",y="Description"))+
  geom_point(mapping=aes_string(size="Count",color="p.adjust"),show.legend = TRUE)+
  scale_color_gradient(low="red",high = "blue")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "Count")) +
  labs(x = "", y = "",color="p.adjust",size="Count")+
  theme_bw()+
  coord_flip()+
  theme(strip.text.x = element_text(angle=0,face="bold",size = 16),strip.text.y=element_text(angle=90,size = 16))+
  theme(axis.text.x = element_text(face = "bold",colour = "black", angle = 90),axis.text.y = element_text(colour = "black"))
print(p)

################################
####-Fig.3&S3: ECs analysis-####
################################
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
### step1 filter
library(Seurat)
library(ggplot2)
library(cowplot)
all<-readRDS("merge_all.rds")
all$CD31<-all@assays$RNA@counts["PECAM1",]
all$CD45<-all@assays$RNA@counts["PTPRC",]
all$epi<-all@assays$RNA@counts["PECAM",]
end<-subset(all,subset = CD31 > 0 & CD45 == 0 &  epi == 0)

### step2 recluster
end_19cts<-subset(end,subset = Tissue %in% c("Liver","Heart","Kidney","Spleen","Adipose-S","Adipose-V",
                                             "Brain", "Cerebellum","AreaPostrema","OVoLT","Retina","Intestine","SubfomicalOrgan",
                                             "Lung","FrontalLobe","Hypothalamus","OccipitalLobe","ParietalLobe","TemporalLobe"))
end_tissues<-SplitObject(end_19cts,split.by = "Tissue")
Liver<-end_tissues$Liver
Heart<-end_tissues$Heart
Kidney<-end_tissues$Kidney
Spleen<-end_tissues$Spleen
Adipose_S<-end_tissues$`Adipose-S`
Adipose_V<-end_tissues$`Adipose-V`
Retina<-end_tissues$Retina
Intestine<-end_tissues$Intestine
Lung<-end_tissues$Lung
Brain<-merge(x = end_tissues$Brain, y = c(end_tissues$Cerebellum, end_tissues$AreaPostrema, end_tissues$OVoLT,
                                          end_tissues$SubfomicalOrgan, end_tissues$FrontalLobe, end_tissues$Hypothalamus,
                                          end_tissues$OccipitalLobe, end_tissues$ParietalLobe, end_tissues$TemporalLobe),
             add.cell.ids = c("Brain", "Cerebellum", "AreaPostrema",
                              "OVoLT","SubfornicalOrgan","FrontalLobe","Hypothalamus",
                              "OccipitalLobe","ParietalLobe","TemporalLobe"))
list<-list()
list[[1]]<-Liver
list[[2]]<-Heart
list[[3]]<-Kidney
list[[4]]<-Spleen
list[[5]]<-Adipose_S
list[[6]]<-Adipose_V
list[[7]]<-Retina
list[[8]]<-Intestine
list[[9]]<-Lung
list[[10]]<-Brain

rm(Liver,Heart,Kidney,Spleen,Adipose_S,Adipose_V,Retina,Intestine,Lung,Brain)
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = list)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features,k.filter = 100)
saveRDS(anchors,"anchors.rds")
all <- IntegrateData(anchorset = anchors)

all <- ScaleData(all, verbose = FALSE)
all<-RunPCA(all)

n_PC<-20
all<-RunUMAP(all,dims = 1:n_PC)
all<-RunTSNE(all,dims = 1:n_PC)

all <- FindNeighbors(all, dims = 1:n_PC)
all <- FindClusters(all, resolution = 1.0)

all$NewTissue<-as.character(all$Tissue)
all$NewTissue[which(all$Tissue == "Liver")] <- "Liver"
all$NewTissue[which(all$Tissue == "Heart")] <- "Heart"
all$NewTissue[which(all$Tissue == "Kidney")] <- "Kidney"
all$NewTissue[which(all$Tissue == "Spleen")] <- "Spleen"
all$NewTissue[which(all$Tissue == "Adipose_S")] <- "Adipose_S"
all$NewTissue[which(all$Tissue == "Adipose_V")] <- "Adipose_V"
all$NewTissue[which(all$Tissue == "Retina")] <- "Retina"
all$NewTissue[which(all$Tissue == "Intestine")] <- "Intestine"
all$NewTissue[which(all$Tissue == "Lung")] <- "Lung"
all$NewTissue[which(all$Tissue %in% c("Brain", "Cerebellum", "AreaPostrema",
                                      "OVoLT","SubfomicalOrgan","FrontalLobe","Hypothalamus",
                                      "OccipitalLobe","ParietalLobe","TemporalLobe"))] <- "Brain"
saveRDS(all,"19cts.Inte.recluster.umaptsne.rds")

### step3 plot
Endo<-readRDS("19cts.Inte.recluster.umaptsne.rds")
VlnPlot(Endo,features = c("PECAM1","PROX1","ICAM2","CDH5","VWF","TIE1","TEK","ACTA2","TAGLN","S100A4","CDH2","SERPINE1","FN1"),
        group.by = "NewTissue",pt.size = 0,slot = "data")

DotPlot(Endo,features = c("PECAM1","PROX1","ICAM2","CDH5","VWF","TIE1","TEK","ACTA2","TAGLN","S100A4","CDH2","SERPINE1","FN1"),
        group.by = "NewTissue",assay = "RNA",scale.by = "size",cols = c("lightgrey","red"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust = 0.5))
data<-p$data

ggplot(data=data,mapping=aes_string(x='id',y='features.plot'))+
  geom_point(mapping=aes_string(size='pct.exp',color='avg.exp.scaled'),show.legend = TRUE)+
  scale_color_gradient(low="blue",high = "red",breaks=NULL)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme_cowplot()+
  coord_flip()+theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust = 0.5))

### step4 filter
Endo$Fib<-Endo@assays$RNA@counts["COL1A1",]
Endo$Peri<-Endo@assays$RNA@counts["PDGFRB",]
Endo$Blood<-Endo@assays$RNA@counts["HBB",]
Endo_sub<-subset(Endo,subset = Fib %in% c(0) & Blood %in% c(0,1,2) & Peri %in% c(0))
Endo_sub.list<-SplitObject(Endo_sub,split.by = "NewTissue")
features <- SelectIntegrationFeatures(object.list = Endo_sub.list)
anchors <- FindIntegrationAnchors(object.list = Endo_sub.list, anchor.features = features,k.filter = 100)
all <- IntegrateData(anchorset = anchors)

all <- ScaleData(all, verbose = FALSE)
all<-RunPCA(all)
n_PC<-20
all<-RunUMAP(all,dims = 1:n_PC)
all<-RunTSNE(all,dims = 1:n_PC)

all <- FindNeighbors(all, dims = 1:n_PC)
all <- FindClusters(all, resolution = 1.0)
saveRDS(all,"Filtered.rds")

col_flg <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(all$seurat_clusters))))
col_flg1 <- colorRampPalette(brewer.pal(12,"Paired"))((length(levels(as.factor(all$NewTissue)))))

p.tsne1 <- DimPlot(object = all, reduction= "tsne", cols = col_flg, label = T, label.size = 5, group.by = "seurat_clusters") 
p.tsne2 <- DimPlot(object = all, reduction= "tsne", cols = col_flg, label = F, group.by = "seurat_clusters") 
pdf(paste0("01_all","_tsne_seurat_clusters.pdf"), w=16, h=6)
p.tsne1 + p.tsne2
dev.off()
p.tsne1 <- DimPlot(object = all, reduction= "tsne", cols = col_flg1, label = T, label.size = 5, group.by = "NewTissue") 
p.tsne2 <- DimPlot(object = all, reduction= "tsne", cols = col_flg1, label = F, group.by = "NewTissue") 
pdf(paste0("01_all","_tsne_NewTissue.pdf"), w=16, h=6)
p.tsne1 + p.tsne2
dev.off()

pdf("02_end_19cts_Heatmap.pdf",width = 10,height = 5)
DoHeatmap(object = all,features = c("PECAM1","EPCAM","PTPRC","COL1A1","PDGFRB"),group.by = "NewTissue",assay = "RNA",slot = "data",lines.width = 1,group.colors = col_flg3,draw.lines = F)+
  scale_fill_gradientn(colors = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

pdf("03_EndMT_markers_NewTissue.pdf",width = 10,height = 30)
DotPlot(all, features = c("PECAM1","CDH5","VWF","TIE1","TEK",
                          "ACTA2","TAGLN","SNAI1","SNAI2","TWIST1","S100A4"),
        group.by = "NewTissue",assay = "RNA",
        cols = c("lightgrey","red")) + coord_flip()
dev.off()
DefaultAssay(all)<-"RNA"
Idents(all)<-"seurat_clusters"
pdf("04_macrophage_Featureplots.pdf",width = 26,height = 18)
FeaturePlot(all, features = c("C1QA","C1QB","C1QC","CD68","MRC1","MARCO","HLA-DRA","AIF1","CD74","CTSS","CSF1R",
                              "MRC1","MSR1","CD36","CD68","CTSS","CTSD","CTSB","CTSZ","CST3","CSTB"),
            slot = "data",cols =c("lightgrey","red"),label = T,reduction = "tsne",min.cutoff = 2) 
dev.off()

GeneBoxplot<-function(Object = Object, Object_Mat = Object_Mat, Feature = Feature){
  require(ggplot2)
  TargetGeneMat <- Object_Mat[Feature,]
  TmpMat <- data.frame(Cluster = factor(x = as.character(Idents(Object)), levels = levels(Idents(Object))), Expression = Object_Mat[Feature,])
  col_flg <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(all$seurat_clusters))))
  p <- ggplot(TmpMat, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_boxplot(show.legend = F, outlier.size = 0, outlier.fill = NULL, outlier.shape = NA) +
    #geom_boxplot(show.legend = F, outlier.size = 1) +
    labs(y = Feature) +
    scale_fill_manual(values = col_flg)+
    theme_bw() +
    theme(axis.text.y = element_text(colour = "black"),
          axis.text.x = element_text(angle = 90, colour = "black"), 
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(size=1,colour="black")) 
  return(p)
}
for (i in c("C1QA","C1QB","C1QC","CD68","MRC1","MARCO","HLA-DRA","AIF1","CD74","CTSS","CSF1R",
            "MRC1","MSR1","CD36","CD68","CTSS","CTSD","CTSB","CTSZ","CST3","CSTB")) {
  pdf(paste0("04_macrophage_Boxplots_",i,".pdf"),width = 8,height = 3)
  print(GeneBoxplot(Object = all,Object_Mat = as.matrix(all@assays$RNA@data),Feature = i))
  dev.off()
}


all.markers<-FindAllMarkers(all,min.pct = 0.5,logfc.threshold = 0.5)
write.table(all.markers,file = "all.markers_0.5.txt",sep = "\t",quote = F)

######VocanoPlots
marker<-read.table("all.markers_0.5.txt",header = T,sep = "\t",quote = "")
dir.create("VocanoPlots")
setwd("VocanoPlots/")
library(ggpubr)
library(ggrepel)

for (i in as.character(unique(marker$cluster))) {
  print(i)
  m14<-marker[which(marker$cluster == i),]
  data<-m14
  data$change = ifelse(data$p_val_adj < 0.001 & abs(data$avg_logFC) >= 1, 
                       ifelse(data$avg_logFC> 1 ,'Up','Down'),
                       'Stable')
  data$gene<-as.character(data$gene)
  data$label=data$gene
  data$label[which(!(data$gene %in% c("C1QA","C1QB","C1QC","CD68","MRC1","MARCO","HLA-DRA","AIF1","CD74","CTSS","CSF1R")))]<-""
  pdf(paste0("Vo_C ",i,".pdf"),width = 10,height = 5)
  p<-ggplot(data = data, 
            aes(x = avg_logFC, 
                y = -log10(data$p_val_adj), 
                colour=change,
                label = data$gene)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("blue","grey","red"))+
    xlim(c(-4.5, 4.5)) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.001),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)",
         title=paste0("C",i," DEGs"))  +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())+
    geom_text_repel(data = data, aes(x = avg_logFC, 
                                     y = -log10(data$p_val_adj), 
                                     label = label),
                    size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE)
  print(p)
  dev.off()
}

##### step4 annotation

##### step5 adipose
all<-readRDS("Filtered.rds")
adipose_sub<-subset(all,subset = NewTissue %in% c("Adipose-S","Adipose-V"))
list<-SplitObject(adipose_sub,split.by = "Tissue")
features <- SelectIntegrationFeatures(object.list = list)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
adipose_sub <- IntegrateData(anchorset = anchors)

adipose_sub <- ScaleData(adipose_sub, verbose = FALSE)
adipose_sub<-RunPCA(adipose_sub)
n_PC<-20
adipose_sub<-RunUMAP(adipose_sub,dims = 1:n_PC)
adipose_sub<-RunTSNE(adipose_sub,dims = 1:n_PC)
adipose_sub <- FindNeighbors(adipose_sub, dims = 1:n_PC)
adipose_sub <- FindClusters(adipose_sub, resolution = 1.0)

setwd("Adipose")
saveRDS(adipose_sub,"Adipose.recluster.umaptsne.rds")

col_flg <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(adipose_sub$seurat_clusters))))

p.tsne1 <- DimPlot(object = adipose_sub, reduction= "tsne", cols = col_flg, label = T, label.size = 5, group.by = "seurat_clusters") 
p.tsne2 <- DimPlot(object = adipose_sub, reduction= "tsne", cols = col_flg, label = F, group.by = "seurat_clusters") 
pdf(paste0("01_adipose_sub","_tsne_seurat_clusters.pdf"), w=16, h=6)
p.tsne1 + p.tsne2
dev.off()

adipose_sub$celltype1<-as.character(adipose_sub$seurat_clusters)
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(6))]<-"Large artery ECs"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(1,9))]<-"Artery ECs"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(16))]<-"Capillary ECs(1)"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(18))]<-"Capillary ECs(2)"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(7,8))]<-"Capillary-venous ECs"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(2,3,4,5,14,15,17))]<-"Vein ECs"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(0,12))]<-"Large vein ECs"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(10))]<-"Lymphatic ECs"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(11))]<-"Scavenging ECs"

adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(13))]<-"EndMT cells"
adipose_sub$celltype1[which(adipose_sub$seurat_clusters %in% c(19))]<-"Proliferating ECs"

library(clusterProfiler)
library(org.Ss.eg.db)
library(org.Hs.eg.db)
genes<-as.character(adipose.celltype.markers$gene[which(adipose.celltype.markers$cluster == "EndMT cells" & adipose.celltype.markers$avg_logFC > 0)])
mygene<-select(org.Hs.eg.db,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL",keys=as.character(genes))
mygene<-mygene$ENTREZID
mygene<-na.omit(mygene)
ego<-enrichGO(OrgDb = "org.Hs.eg.db",gene=mygene,readable=TRUE,pvalueCutoff = 1)
ego<-data.frame(ego,stringsAsFactors = F)
write.table(ego,paste("Mesenchymal_cells","_GO_pig.txt",sep = ""),sep="\t",row.names = F,quote = F)

ekegg<-enrichKEGG(gene=mygene,pvalueCutoff=0.05,organism = "hsa")
ekegg<-data.frame(ego,stringsAsFactors = F)
write.table(ekegg,paste("Mesenchymal_cells","_KEGG.txt",sep = ""),sep="\t",row.names = F,quote = F)
saveRDS(adipose_sub,"Adipose.celltype1.rds")

marker<-FindAllMarkers(adipose_sub,group.by = "celltype1",assay = "RNA",slot = "data",features = features)

##### step 6 monocle2
library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)
adipose_sub<-readRDS("Adipose.celltype1.rds")
setwd("monocle2/EndMT/")
levels(adipose_sub)
endmt<-subset(adipose_sub,subset = celltype1 %in% "EndMT cells")

DefaultAssay(endmt)<-"RNA"
expr_matrix <- as(as.matrix(endmt@assays$RNA@counts), 'sparseMatrix')
p_data <- endmt@meta.data 
p_data$celltype <- endmt@active.ident  
f_data <- data.frame(gene_short_name = row.names(endmt),row.names = row.names(endmt))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds_endmt <- newCellDataSet(expr_matrix,
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())

cds_endmt <- estimateSizeFactors(cds_endmt)
cds_endmt <- estimateDispersions(cds_endmt)
cds_endmt <- detectGenes(cds_endmt, min_expr = 0.1) 
print(head(fData(cds_endmt)))#
expressed_genes <- row.names(subset(fData(cds_endmt),
                                    num_cells_expressed >= 10)) 


diff <- differentialGeneTest(cds_endmt[expressed_genes,],fullModelFormulaStr="~celltype",cores=1) 
head(diff)

deg <- subset(diff, qval < 0.01) #
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

ordergene <- rownames(deg) 
cds_endmt <- setOrderingFilter(cds_endmt, ordergene)  

pdf("train.ordergenes.pdf")
plot_ordering_genes(cds_endmt)
dev.off()

cds_endmt <- reduceDimension(cds_endmt, max_components = 2,
                             method = 'DDRTree')
cds_endmt <- orderCells(cds_endmt)

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(cds_endmt,color_by="Pseudotime", size=1,show_backbone=TRUE) 
dev.off()

pdf("train.monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds_endmt,color_by="celltype", size=1,show_backbone=TRUE)
dev.off()


pdf("train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds_endmt, color_by = "State",size=1,show_backbone=TRUE)
dev.off()

plot_cell_trajectory(cds_endmt, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_cell_trajectory(cds_endmt, color_by = "celltype") + facet_wrap("~State", nrow = 1)

p1 <- plot_cell_trajectory(cds_endmt, x = 1, y = 2, color_by = "celltype") + 
  theme(legend.position='none',panel.border = element_blank())  
p2 <- plot_complex_cell_trajectory(cds_endmt, x = 1, y = 2,
                                   color_by = "celltype")+
  #scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 
p1|p2

saveRDS(cds_endmt,"monocle2_endmtcells.rds")

keygenes <- c("PECAM1","VWF","ICAM1","CDH5","ICAM2",
              "ACTA2","TAGLN","VIM","DCN","CDH2",
              "TGFB2","TGFB3","TGFBR2","TGFBR3","NOTCH3",
              "S100A4","CD44","SNAI1","SNAI2","TWIST1","TWIST2","CNN1",
              "ZEB2","COL1A2","COL3A1",
              "ADIRF")
cds_subset <- cds_endmt[keygenes,]

p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 30)

p1 <- plot_genes_jitter(cds_endmt[keygenes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds_endmt[keygenes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(cds_endmt[keygenes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 20)

colnames(pData(cds_endmt))
pData(cds_endmt)$ACTA2 = log2( exprs(cds_endmt)['ACTA2',]+1)
p1=plot_cell_trajectory(cds_endmt, color_by = "ACTA2")  + scale_color_gsea()
pData(cds_endmt)$PECAM1 = log2(exprs(cds_endmt)['PECAM1',]+1)
p2=plot_cell_trajectory(cds_endmt, color_by = "PECAM1")    + scale_color_gsea()
library(patchwork)
p1+p2

#ordergene <- rownames(deg) 
Time_diff <- differentialGeneTest(cds_endmt, cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds_endmt[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

###### Figures
marker_genes1 <- row.names(subset(fData(cds_endmt),
                                  gene_short_name %in% c("PECAM1","VWF","ICAM1","CDH5","ICAM2",
                                                         "ACTA2","TAGLN","VIM","CDH2",
                                                         "TGFB2","TGFB3","TGFBR2","NOTCH3",
                                                         "S100A4","CD44","SNAI1","SNAI2","TWIST1","TWIST2","CNN1",
                                                         "ZEB2","FABP4")))
my_palette <- colorRampPalette(c("#2c7fb8", "#ffffd9", "#e34a33"))(n = 10000)
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 10000)
my_palette <- colorRampPalette(c("#2c7fb8", "white", "#e34a33"))(n = 10000)
pdf("01_markers_heatmap.pdf")
plot_pseudotime_heatmap(cds_endmt[marker_genes1,],
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = T,hmcols = my_palette,return_heatmap = T)
dev.off()

selected_genes <- c("PECAM1","VWF","ICAM1",
                    "ACTA2","TAGLN","VIM","CDH2","CNN1",
                    "TGFB2","TGFB3","NOTCH3",
                    "S100A4","SNAI1",
                    "ZEB2","FABP4")
pdf("02_markers_pseudotime.pdf")
plot_genes_in_pseudotime(cds_endmt[selected_genes,], color_by = "State",ncol = 3)
dev.off()

pdf("SNAI1.pdf",width = 3,height = 1.6)
plot_genes_in_pseudotime(cds_endmt["SNAI1",], color_by = "State",min_expr = 0)
dev.off()

################################
####-Fig.4&S4: ECs analysis-####
################################
### Aorta
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
raw.data<-Read10X(data.dir = "./raw data/EC-Aorta/")
Sample<-"EC_Aorta"
object <- CreateSeuratObject(counts = raw.data, min.cells = 3, min.features = 200, project = "EC_Lung")
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
Vlnplot <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dir.create("test")
pdf(paste0("01_Vlnplot_","EC_aorta",".pdf"), w=12, h=9)
Vlnplot
dev.off()

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("02_FeatureScatter_","EC_Aorta",".pdf"), w=10, h=5)
plot1 + plot2
dev.off()

############## Control data
data <- subset(x = object, subset = nFeature_RNA > 200 & percent.mt < 25 )
dim(data)
############## Double Finder
data <- NormalizeData(object = data, normalization.method = "LogNormalize")
data <- FindVariableFeatures(object = data, mean.function = ExpMean, dispersion.function = LogVMR)
data <- ScaleData(object = data, vars.to.regress = c("nCount_RNA"))
data <- RunPCA(data, verbose = FALSE)
pdf('03_ElbowPlot.pdf')
ElbowPlot(data, ndims=50)  ##change dims
dev.off()
data.dims = 20  ##change dims

data <- FindNeighbors(data, dims = 1:data.dims)   #dim
data <- FindClusters(data, resolution = 1)
data <- RunTSNE(data, reduction.use =  "pca", dims.use = 1:data.dims, do.fast = T)   #dim

sweep.res.list_SM <- paramSweep_v3(data, PCs = 1:data.dims)  #dim
sweep.staDMH_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.staDMH_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- data@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.075*length(data@meta.data$orig.ident))  #cell_load
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
data <- doubletFinder_v3(data, PCs = 1:data.dims, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)  #dim
data <- doubletFinder_v3(data, PCs = 1:data.dims, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  #dim
df.classifications <- paste("DF.classifications", pN_value, pK_value, nExp_poi.adj, sep="_")
data@meta.data$Doublet <- data@meta.data[ ,df.classifications]

#### plot
p.doublet <- DimPlot(data,reduction="tsne", group.by = "Doublet", label = TRUE) + ggtitle(Sample) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0("04_Doublet_","EC_Aorta",".pdf"), w=5, h=5)
p.doublet
dev.off()

data1 <- subset(data, subset = Doublet == 'Singlet')
p.doubletfeature <- VlnPlot(data1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pdf(paste0("05_Doublet_Feature_","EC_Aorta",".pdf"), w=10, h=5)
p.doubletfeature
dev.off()

EC <- data1
dim(EC)

##### normalization
EC <- NormalizeData(object = EC, normalization.method = "LogNormalize", scale.factor = 10000)
EC <- FindVariableFeatures(object = EC, mean.function = ExpMean, dispersion.function = LogVMR)

hv.genes <- head(VariableFeatures(EC), 2000)  ##???
length(VariableFeatures(EC))

##### Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EC), 10)
plot1 <- VariableFeaturePlot(EC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("06_VariableFeature.pdf")
plot2
dev.off()

##### scale data
#EC <- ScaleData(object = EC, genes.use = hv.genes, display.progress = FALSE, vars.to.regress = c("percent.mt","nCount_RNA"), do.par = TRUE, num.cores = 4)
EC <- ScaleData(object = EC)
EC <- RunPCA(object = EC, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
ElbowPlot(EC, ndims = 50)

##### Run TSNE and Find Cluster
EC <- FindNeighbors(EC, dims=1:data.dims)
EC <- FindClusters(EC, resolution = 1)
EC <- RunTSNE(object = EC, reduction = "pca", dims= 1:data.dims, reduction.name = "tSNE")  #dim
EC <- RunUMAP(object = EC, reduction = "pca", dims= 1:data.dims, reduction.name = "umap")  #dim

aorta<-EC
aorta$celltype1<-as.character(aorta$seurat_clusters)
aorta$celltype1[which(aorta$seurat_clusters %in% c(0,6))]<-"ECs(G1)"
aorta$celltype1[which(aorta$seurat_clusters %in% c(1,3,7))]<-"Proliferating ECs(S/G2M)"
aorta$celltype1[which(aorta$seurat_clusters %in% c(4))]<-"Intermediate ECs"
aorta$celltype1[which(aorta$seurat_clusters %in% c(2))]<-"Fibroblasts"
aorta$celltype1[which(aorta$seurat_clusters %in% c(5))]<-"Mesenchymal cells"
aorta$celltype1<-factor(aorta$celltype1,levels = c("ECs(G1)","Intermediate ECs","Proliferating ECs(S/G2M)","Mesenchymal cells","Fibroblasts"))

p.tsne1 <- DimPlot(object = aorta, reduction= "tsne", cols = colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(aorta$celltype1)))) , label = T, label.size = 5, group.by = "celltype1") 
p.tsne2 <- DimPlot(object = aorta, reduction= "tsne", cols = colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(aorta$celltype1)))) , label = F, group.by = "celltype1") 
pdf(paste0("05_aorta","_tsne_subtype.pdf"), w=16, h=6)
p.tsne1 + p.tsne2
dev.off()
saveRDS(aorta,file = "Aorta_celltype1.rds")


aorta<-readRDS("Aorta_celltype1.rds")
pdf("aorta_dotplot.pdf",width = 10)
DotPlot(aorta,features = rev(c("PECAM1","VCAM1","ICAM2",
                               "RPL7A","RPL3",
                               "TOP2A","PCNA","MCM6",
                               "ACTA2","TAGLN","S100A4",
                               "COL1A1","COL1A2")), 
        cols = c("lightgrey","red"),
        assay = "RNA",group.by = "celltype1")+coord_flip()
dev.off()


Idents(aorta)<-"celltype1"
sample.cluster<-AverageExpression(aorta)$RNA
hc = hclust(dist(t(sample.cluster)))
hcd = as.dendrogram(hc)
pdf("hclust_celltype1.pdf")
plot(hcd)
dev.off()

aorta<-readRDS("Aorta_celltype1.rds")
data<-data.frame(Phase = aorta$Phase,celltype1 = aorta$celltype1)
data<-as.data.frame(table(aorta$Phase,aorta$celltype1))
colnames(data)<-c("Stage","celltype","value")
data1<-c()
for (i in unique(data$celltype)) {
  tmp<-data[which(data$celltype %in% i),]
  sum<-sum(tmp$value)
  tmp$percentage<-tmp$value/sum
  data1<-rbind(data1,tmp)
}

data1$Stage<-factor(data1$Stage,levels = c("G1","S","G2M"))
pdf("aorta_cellcycle_order.pdf")
ggplot(data1, aes(x=celltype, y=percentage, fill=Stage))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = "Accent")
dev.off()

##### Lung
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
raw.data<-Read10X(data.dir = "./raw data/EC-Lung/")
Sample<-"EC_Lung"
object <- CreateSeuratObject(counts = raw.data, min.cells = 3, min.features = 200, project = "EC_Lung")
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
Vlnplot <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dir.create("test")
pdf(paste0("01_Vlnplot_","EC_lung",".pdf"), w=12, h=9)
Vlnplot
dev.off()

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("02_FeatureScatter_","EC_lung",".pdf"), w=10, h=5)
plot1 + plot2
dev.off()

############## Control data
data <- subset(x = object, subset = nFeature_RNA > 200 & percent.mt < 25 )
dim(data)
############## 
data <- NormalizeData(object = data, normalization.method = "LogNormalize")
data <- FindVariableFeatures(object = data, mean.function = ExpMean, dispersion.function = LogVMR)
data <- ScaleData(object = data, vars.to.regress = c("nCount_RNA"))
data <- RunPCA(data, verbose = FALSE)
pdf('03_ElbowPlot.pdf')
ElbowPlot(data, ndims=50)  ##change dims
dev.off()
data.dims = 20  ##change dims

data <- FindNeighbors(data, dims = 1:data.dims)   #dim
data <- FindClusters(data, resolution = 0.5)
data <- RunTSNE(data, reduction.use =  "pca", dims.use = 1:data.dims, do.fast = T)   #dim

lung<-data
lung$celltype1<-as.character(lung$seurat_clusters)
lung$celltype1[which(lung$seurat_clusters %in% c(0,3,5))]<-"ECs(G1)"
lung$celltype1[which(lung$seurat_clusters %in% c(9))]<-"Intermidiate ECs"
lung$celltype1[which(lung$seurat_clusters %in% c(1,8))]<-"Proliferating ECs(S)"
lung$celltype1[which(lung$seurat_clusters %in% c(2,4))]<-"Proliferating ECs(G2M)"
lung$celltype1[which(lung$seurat_clusters %in% c(6,7))]<-"Fibroblasts"
lung$celltype1<-factor(lung$celltype1,levels = c("ECs(G1)","Intermidiate ECs","Proliferating ECs(S)","Proliferating ECs(G2M)","Fibroblasts"))
saveRDS(lung,"lung_celltype1.rds")
p.tsne1 <- DimPlot(object = lung, reduction= "tsne", cols = colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(lung$celltype1)))) , label = T, label.size = 5, group.by = "celltype1") 
p.tsne2 <- DimPlot(object = lung, reduction= "tsne", cols = colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(lung$celltype1)))) , label = F, group.by = "celltype1") 
pdf(paste0("05_lung","_tsne_subtype.pdf"), w=16, h=6)
p.tsne1 + p.tsne2
dev.off()
pdf("lung_dotplot_new.pdf",width = 10)
DotPlot(lung,features = rev(c("PECAM1","VCAM1","ICAM2",
                              "RPL3","RPL7A",
                              "TOP2A","PCNA","MCM6",
                              "ACTA2","TAGLN","S100A4",
                              "COL1A1","COL1A2")), 
        cols = c("lightgrey","red"),
        assay = "RNA",group.by = "celltype1")+coord_flip()
dev.off()

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
lung <- CellCycleScoring(lung, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

lung<-readRDS("lung_celltype1.rds")
data<-data.frame(Phase = lung$Phase,celltype1 = lung$celltype1)
data<-as.data.frame(table(lung$Phase,lung$celltype1))
colnames(data)<-c("Stage","celltype","value")
data1<-c()
for (i in unique(data$celltype)) {
  tmp<-data[which(data$celltype %in% i),]
  sum<-sum(tmp$value)
  tmp$percentage<-tmp$value/sum
  data1<-rbind(data1,tmp)
}

data1$Stage<-factor(data1$Stage,levels = c("G1","S","G2M"))
pdf("lung_cellcycle_order.pdf")
ggplot(data1, aes(x=celltype, y=percentage, fill=Stage))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = "Accent")
dev.off()

################################
######-Fig.5&S5: CellChat-######
################################
# Cell communications analysis by CellChat (version 0.0.1)

# ================================ CellChat Start ================================

# Load the required libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#----------------- Part I: Data input & processing and initialization of CellChat object -----------------

# ----------- Pig Heart ----------- 
## Load data
Pig.Heart <- readRDS("Pig.Heart.celltype.rds")

## Create a CellChat object
Pig_cellchat <- createCellChat(data = Pig.Heart@assays$RNA@data)
identity <- data.frame(group = Pig$celltype, row.names = names(Pig$celltype))
identity$group <- factor(identity$group, levels = c("Cardiomyocytes","Endothelial cells","Fibroblasts","Myeloid cells","Lymphoid cells"))
unique(identity$group)

# ----------- Human Heart ----------- 
## Load data
Seurat_object <- readRDS("Human.Heart.celltype.rds")

## Create a CellChat object
Cellchat_object <- createCellChat(data = Seurat_object@assays$RNA@data)
identity <- data.frame(group = Seurat_object$celltype, row.names = names(Seurat_object$celltype))
identity$group <- factor(identity$group, levels = c("Cardiomyocytes","Endothelial cells","Fibroblasts","Myeloid cells","Lymphoid cells"))
unique(identity$group)

# ----------- Pig Liver ----------- 
## Load data
Seurat_object <- readRDS("Pig.Liver.celltype.rds")

## Create a CellChat object
Cellchat_object <- createCellChat(data = Seurat_object@assays$RNA@data)
identity <- data.frame(group = Seurat_object$celltype, row.names = names(Seurat_object$celltype))
identity$group <- factor(identity$group, levels = c("Hepatocytes","Endothelial cells","Kupffer cells","B cells","T/NK cells","Erythroid cells"))
unique(identity$group)

# ----------- Human Liver ----------- 
## Load data
Seurat_object <- readRDS("Human.Liver.celltype.rds")

## Create a CellChat object
Cellchat_object <- createCellChat(data = Seurat_object@assays$RNA@data)
identity <- data.frame(group = Seurat_object$celltype, row.names = names(Seurat_object$celltype))
identity$group <- factor(identity$group, levels = c("Hepatocytes","Endothelial cells","Kupffer cells","B cells","T/NK cells","Erythroid cells"))
unique(identity$group)

# ----------- Pig Kidney ----------- 
## Load data
Seurat_object <- readRDS("Pig.Kidney.celltype.rds")

## Create a CellChat object
Cellchat_object <- createCellChat(data = Seurat_object@assays$RNA@data)
identity <- data.frame(group = Seurat_object$celltype, row.names = names(Seurat_object$celltype))
identity$group <- factor(identity$group, levels = c("Epithelial cells","Podocytes","Proximal tubule cells","Collecting duct cells","Endothelial cells","Distal convoluted tubule cells"))
unique(identity$group)

# ----------- Human Kidney ----------- 
## Load data
Seurat_object <- readRDS("Human.Kidney.celltype.rds")

## Create a CellChat object
Cellchat_object <- createCellChat(data = Seurat_object@assays$RNA@data)
identity <- data.frame(group = Seurat_object$celltype, row.names = names(Seurat_object$celltype))
identity$group <- factor(identity$group, levels = c("Epithelial cells","Podocytes","Proximal tubule cells","Collecting duct cells","Endothelial cells","Distal convoluted tubule cells"))
unique(identity$group)


#----------------- Part II: Perform standered workflow of CellChat based on specific Cellchat object-----------------

## Add cell information into meta slot of the object
Cellchat_object <- addMeta(Cellchat_object, meta = identity, meta.name = "labels")
Cellchat_object <- setIdent(Cellchat_object, ident.use = "labels") # set "labels" as default cell identity
levels(Cellchat_object@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(Cellchat_object@idents)) # number of cells in each cell group

## Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human
showDatabaseCategory(CellChatDB)

## Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

## Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis

## Set the used database in the object
Cellchat_object@DB <- CellChatDB.use

## Preprocessing the expression data for cell-cell communication analysis
### subset the expression data of signaling genes for saving computation cost
Cellchat_object <- subsetData(Cellchat_object) 
future::plan("multiprocess", workers = 4) # do parallel
Cellchat_object <- identifyOverExpressedGenes(Cellchat_object)
Cellchat_object <- identifyOverExpressedInteractions(Cellchat_object)
### Project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
Cellchat_object <- projectData(Cellchat_object, PPI.human)


#----------------- Part III: Inference of cell-cell communication network ----------------- 

## Compute the communication probability and infer cellular communication network
Cellchat_object <- computeCommunProb(Cellchat_object) 

## Infer the cell-cell communication at a signaling pathway level
Cellchat_object <- computeCommunProbPathway(Cellchat_object)

## Calculate the aggregated cell-cell communication network
Cellchat_object <- aggregateNet(Cellchat_object)


#----------------- Part VI: Visualization and systems analysis of cell-cell communication network ----------------- 

pathways.show <- c("VEGF")  # Fig5A
pathways.show <- c("PDGF")  # Fig5B
pathways.show <- c("TGFb")  # Fig5(S5)B
pathways.show <- c("BMP")   # Fig5(S5)B

groupSize <- as.numeric(table(Cellchat_object@idents))

## Visualize each signaling pathway using circle plot (Fig5A, Fig5B, FigS5B)
netVisual_aggregate(Cellchat_object, signaling = pathways.show, layout = "circle", vertex.weight = groupSize)

# Compute the network centrality scores
Cellchat_object <- netAnalysis_computeCentrality(Cellchat_object, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

## Visualize the computed centrality scores using heatmap  (Fig5A, Fig5B, FigS5B)
netAnalysis_signalingRole_network(Cellchat_object, signaling = pathways.show, slot.name = "netP")

## Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway (FigS5A)
netAnalysis_contribution(Cellchat_object, signaling = pathways.show)

# ================================ CellChat End ================================

################################
#######-Fig.6&S6: GENIE3-#######
################################
### merge pig all brain
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)

files<-read.table("all_brain_rds.txt",header=F,sep="\t")
brain.list<-list()
for (i in 1:nrow(files)){
 object<-readRDS(as.character(files[i,3]))
 object$Tissue<-as.character(files[i,2])
 brain.list[[i]]<-object
}
features <- SelectIntegrationFeatures(object.list = brain.list)
brain.anchors <- FindIntegrationAnchors(object.list = brain.list, anchor.features = features)
brain.anchors<-readRDS("brain.anchors.rds")
brain.combined <- IntegrateData(anchorset = brain.anchors)
DefaultAssay(brain.combined) <- "integrated"
brain.combined <- ScaleData(brain.combined, verbose = FALSE)
brain.combined <- RunPCA(brain.combined, npcs=20, verbose=FALSE)
brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:20)
brain.combined <- RunTSNE(brain.combined, reduction = "pca", dims = 1:20)
brain.combined <- FindNeighbors(brain.combined, reduction = "pca", dims = 1:20)
brain.combined <- FindClusters(brain.combined, resolution = 1.0)
saveRDS(brain.combined,"integrated_all_brain.rds")
object<-brain.combined
col_flg <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(object$Celltype))))
col_flg2 <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(as.factor(object$Tissue))))

p.tsne1 <- DimPlot(object = object, reduction= "tsne", cols = col_flg, label = T, label.size = 5, group.by = "Celltype") 
p.tsne2 <- DimPlot(object = object, reduction= "tsne", cols = col_flg, label = F, group.by = "Celltype") 
pdf(paste0("01_all","_tsne_Celltype.pdf"), w=32, h=6)
p.tsne1 + p.tsne2
dev.off()
p.tsne1 <- DimPlot(object = object, reduction= "tsne", cols = col_flg2, label = T, label.size = 5, group.by = "Tissue") 
p.tsne2 <- DimPlot(object = object, reduction= "tsne", cols = col_flg2, label = F, group.by = "Tissue") 
pdf(paste0("01_all","_tsne_Tissue.pdf"), w=16, h=6)
p.tsne1 + p.tsne2
dev.off()

### run GENIE3
library(Seurat)
library(GENIE3)
TF<-read.table("/ref/Sus_scrofa_TF.txt",header = T,sep = "\t",stringsAsFactors = F)
PCGs<-read.table("/ref/pig_prcGene.txt",header = T,sep = "\t",stringsAsFactors = F)

### pig 10x brain
object<-readRDS("MG_Homo_10x_Brain.rds")
object
sub<-object[,sample(colnames(object),size=1000,replace=F)]
saveRDS(sub,"MG_Homo_10x_Brain_1000.rds")

object<-sub
print(Sys.time())
set.seed(123) # For reproducibility of results
exprMatr<-as.matrix(object@assays$RNA@data)
rownames(exprMatr)<-rownames(object@assays$RNA@data)
colnames(exprMatr)<-colnames(object@assays$RNA@data)
myfun<-function(x){length(which(x!=0))>=dim(exprMatr)[2]*0.05}
test<-apply(exprMatr, 1, myfun)
test<-unlist(test)
weightMat <- GENIE3(exprMatr, regulators=intersect(names(which(test==TRUE)),intersect(TF$Symbol,rownames(object))),targets=intersect(names(which(test==TRUE)),intersect(PCGs$Gene.name,rownames(object))),nCores=16, verbose=TRUE)
saveRDS(weightMat,paste0("MG","_weightMat.rds"))
linkList <- getLinkList(weightMat)
linkList <- getLinkList(weightMat, threshold=0.01)
write.table(linkList,paste0("MG","_linkList_0.01.txt"),row.names = F,sep = "\t",quote = F)
print(Sys.time())
### public data
obj_list<-readRDS("MG_pub_list.rds")
myfun<-function(x){length(which(x!=0))>=dim(exprMatr)[2]*0.05}
pub_mg<-read.table("./pub_stat_0708.txt",row.names = 1,header = T,sep = "\t",stringsAsFactors = F)
for (i in 1:length(obj_list)) {
  set.seed(123) 
  object<-obj_list[[i]]
  exprMatr<-as.matrix(object@assays$RNA@data)
  rownames(exprMatr)<-rownames(object@assays$RNA@data)
  colnames(exprMatr)<-colnames(object@assays$RNA@data)
  test<-apply(exprMatr, 1, myfun)
  test<-unlist(test)
  weightMat <- GENIE3(exprMatr, regulators=intersect(names(which(test==TRUE)),intersect(TF$Symbol,rownames(object))),targets=intersect(names(which(test==TRUE)),intersect(PCGs$Gene.name,rownames(object))),nCores=4, verbose=TRUE)
  saveRDS(weightMat,paste("pub_",pub_mg[i,3],"_weightMat.rds",sep = ""))
}

### get weight list
list<-read.table("list.txt",header=F,sep="\t")
for(i in 1:nrow(list)){
  weightMat<-readRDS(as.character(list[i,2]))
  linkList <- getLinkList(weightMat, threshold=0.01)
  write.table(linkList,paste0(as.character(list[i,1]),"_linkList_0.01.txt"),row.names = F,sep = "\t",quote = F)
}
### write table
list0.01<-list.files(path = "./03_GRNs",pattern = "0.01",full.names = T)
names<-gsub("_linkList_0.01.txt","",list0.01)
names<-gsub("./03_GRNs/","",names)
all<-c()
for (i in 1:length(list0.01)) {
  object<-read.table(list0.01[i],header = T,sep = "\t")
  object$Species<-names[i]
  all<-rbind(all,object)
}
write.table(all,"all_0.01.txt",row.names = F,sep = "\t",quote = F)
### plot
library(reshape2)
stat<-dcast(all,regulatoryGene~Species)
rownames(stat)<-stat$regulatoryGene
stat<-stat[,-1]
stat[stat > 0]<-1
stat<-stat[,-grep("07_Pig_1000cells",colnames(stat))]
stat$sum<-rowSums(stat)
write.table(stat,"stat_TF.txt",row.names=T,sep = "\t",quote = F)

stat1<-dcast(all,regulatoryGene+targetGene~Species)
rownames(stat1)<-paste(stat1$regulatoryGene,stat1$targetGene,sep = "_")
stat1<-stat1[,c(-1,-2)]
stat1[!is.na(stat1)]<-1
stat1[is.na(stat1)]<-0
write.table(stat1,"stat_TF_target.txt",row.names=T,sep = "\t",quote = F)
stat1<-read.table("stat_TF_target.txt",header = T,sep = "\t",row.names = 1,check.names = F)
stat1<-stat1[,-grep("07_Pig_1000cells",colnames(stat1))]
stat1$sum<-rowSums(stat1)
write.table(stat1,"stat_TF_target.txt",row.names=T,sep = "\t",quote = F)


conserved_TFs<-rownames(stat[which(stat$sum > 8),])
for (i in 1:length(list0.01)) {
  object<-read.table(list0.01[i],header = T,sep = "\t")
  object$Species<-names[i]
  object<-object[which(object$regulatoryGene == conserved_TFs),]
  relations <- data.frame(from=object$regulatoryGene,
                          to=object$targetGene,
                          #paris=f$pair,
                          weight=object$weight
                          #mode=f$mode
  )
  g <- graph_from_data_frame(relations, directed=TRUE,vertices=unique(union(object$targetGene,object$regulatoryGene)))
  #x<-length(unique(V(g)))
  #print(x)
  V(g)[as.character(unique(object$targetGene))]$color<-"#6CB92D"
  V(g)[as.character(unique(object$regulatoryGene))]$color<-"#3DBBC3"
  #V(g)$count<-count_source[V(g),2] 
  V(g)[as.character(unique(object$targetGene))]$size<-3
  V(g)[as.character(unique(object$regulatoryGene))]$size <- 7
  V(g)[as.character(unique(object$targetGene))]$label.cex<-0.001
  V(g)[as.character(unique(object$targetGene))]$label.color<-"#6CB92D"
  V(g)[as.character(unique(object$regulatoryGene))]$label.cex <- 0.4
  V(g)[as.character(unique(object$regulatoryGene))]$label.color <- "red"
  E(g)$arrow.size <- 0
  #E(g)$edge.color <- "gray80"
  E(g)$width <- E(g)$weight*0.01
  
  #l <- layout_as_star(g)
  #l <- layout_with_fr(g)
  l<-layout_nicely(g)
  pdf(paste0("species",names[i],"_MG_GRNs.pdf"),width = 6,height = 6)
  plot(g,layout=l,
       #vertex.label.color="red",
       #vertex.label.cex=0.6,
       vertex.frame.color="white",
       vertex.color=V(g)$color,
       vertex.shape="circle",
       vertex.size=V(g)$size,
       edge.width=E(g)$width,
       edge.arrow.width=0,
       edge.arrow.size=0,
       arrow.mode="0",
       edge.color="gray50",
       #mark.groups=unique(V(g)),
       #mark.border=0
  )
  #plot(p)
  dev.off()
}


conserved_TF_targets<-rownames(stat1[which(stat1$sum > 4),])
conserved_TF_targets_TF<-gsub("_.*","",conserved_TF_targets)
conserved_TF_targets_target<-gsub(".*_","",conserved_TF_targets)
for (i in 1:length(list0.01)) {
  object<-read.table(list0.01[i],header = T,sep = "\t")
  object$Species<-names[i]
  object$pairs<-paste(object$regulatoryGene,object$targetGene,sep = "_")
  #object<-object[which(object$regulatoryGene %in% conserved_TF_targets_TF & object$targetGene %in% conserved_TF_targets_target),]
  object<-object[which(object$pairs %in% conserved_TF_targets),]
  relations <- data.frame(from=object$regulatoryGene,
                          to=object$targetGene,
                          #paris=f$pair,
                          weight=object$weight
                          #mode=f$mode
  )
  g <- graph_from_data_frame(relations, directed=TRUE,vertices=unique(union(object$targetGene,object$regulatoryGene)))
  #x<-length(unique(V(g)))
  #print(x)
  V(g)[as.character(unique(object$targetGene))]$color<-"#6CB92D"
  V(g)[as.character(unique(object$regulatoryGene))]$color<-"#3DBBC3"
  #V(g)$count<-count_source[V(g),2] 
  V(g)[as.character(unique(object$targetGene))]$size<-3
  V(g)[as.character(unique(object$regulatoryGene))]$size <- 7
  V(g)[as.character(unique(object$targetGene))]$label.cex<-0.001
  V(g)[as.character(unique(object$targetGene))]$label.color<-"#6CB92D"
  V(g)[as.character(unique(object$regulatoryGene))]$label.cex <- 0.4
  V(g)[as.character(unique(object$regulatoryGene))]$label.color <- "red"
  E(g)$arrow.size <- 0
  #E(g)$edge.color <- "gray80"
  E(g)$width <- E(g)$weight*0.01
  
  #l <- layout_as_star(g)
  #l <- layout_with_fr(g)
  l<-layout_nicely(g)
  pdf(paste0("pairs_",names[i],"_MG_GRNs.pdf"),width = 6,height = 6)
  plot(g,layout=l,
       #vertex.label.color="red",
       #vertex.label.cex=0.6,
       vertex.frame.color=NA,
       vertex.color=V(g)$color,
       vertex.shape="circle",
       vertex.size=V(g)$size,
       edge.width=E(g)$width,
       edge.arrow.width=0,
       edge.arrow.size=0,
       arrow.mode="0",
       edge.color="gray50",
       #mark.groups=unique(V(g)),
       #mark.border=0
  )
  #plot(p)
  dev.off()
}

