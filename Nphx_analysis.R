# Nasopharynx data (with mm10 reference)
# 2023-02-20

wd = "~/Projects/Nasopharynx/Nphx"
setwd(wd)
source("~/Projects/Organotypic_LEC/pre_processing_SS3.R")
library(dplyr)
adult_matrix = make_count_matrix_from_fc_out(
  fc_dir = "./Adult/FC_OUT/combined.counts",
  gtf_dir = "~/Reference/ss3/gencode.vM16.primary_assembly.annotation.gtf"
)

aged_matrix = make_count_matrix_from_fc_out(
  fc_dir = "./Aged/FC_OUT/combined.counts",
  gtf_dir = "~/Reference/ss3/gencode.vM16.primary_assembly.annotation.gtf"
)

adult_plate_info = make_plate_info(num_plate = 16, sample_prefix = "Adult",count_matrix = adult_matrix)
aged_plate_info  = make_plate_info(num_plate = 12, sample_prefix = "Aged",count_matrix = aged_matrix)
# 
source("~/Projects/Seurat_mjyang_function.R")
adult = CreateSeuratObject(counts = adult_matrix, project = "Adult_nphx", min.cells = 3)
aged = CreateSeuratObject(counts = aged_matrix, project = "Aged_nphx", min.cells = 3)

adult$plate = adult_plate_info$plate_index_vec
aged$plate = aged_plate_info$plate_index_vec

adult$percent.mito = PercentageFeatureSet(adult, pattern = "^mt-")
aged$percent.mito = PercentageFeatureSet(aged, pattern = "^mt-")

adult$sample = "adult"
aged$sample = "aged"
np = merge(adult, aged)
np = subset(np, percent.mito < 10 & nFeature_RNA > 2000)

ribo_genes = grep(rownames(np) , pattern = "^Rpl|^Rps", value = T)
np = AddModuleScore(np, features = list(genes.dissoc), name = "dissoc")
np = AddModuleScore(np, features = list(ribo_genes), name = "ribo")

set.seed(7)
np = np %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst") %>% ScaleData(vars.to.regress = c("nCount_RNA", "percent.mito", "dissoc1","ribo1")) %>% RunPCA()
np = np %>% RunUMAP(reduction = "pca", dims = 1:17) %>% FindNeighbors(reduction ="pca", dims = 1:17) %>% FindClusters(resolution = 0.3)
DimPlot(np)

# Annotate
np$CellType = ""
np$CellType[WhichCells(np, idents = 0)] = "Col"
np$CellType[WhichCells(np, idents = 1)] = "Ptx3_Cap"
np$CellType[WhichCells(np, idents = 2)] = "Cap"
np$CellType[WhichCells(np, idents = 3)] = "Up_Valve"
np$CellType[WhichCells(np, idents = 4)] = "Down_Valve"

# Adult
adult = subset(np, sample %in% c("adult"))
adult = adult %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst") %>% ScaleData(vars.to.regress = c("nCount_RNA", "percent.mito", "dissoc1","ribo1")) %>% RunPCA()
adult = adult %>% RunUMAP(reduction = "pca", dims = 1:15, seed.use = 555) %>% FindNeighbors(reduction ="pca", dims = 1:15) %>% FindClusters(resolution = 0.3)
adult = RunUMAP(adult, reduction = "pca", dims = 1:15, seed.use = 77, min.dist = 0.4)
adult = FindClusters(adult , resolution = 0.4)
plot_grid(
  DimPlot(adult),
  DimPlot(adult, group.by = "CellType")
)

adult$CellType[WhichCells(adult, idents = 0)] = "Col"
adult$CellType[WhichCells(adult, idents = 1)] = "Ptx3_Cap"
adult$CellType[WhichCells(adult, idents = 2)] = "Up_Valve"
adult$CellType[WhichCells(adult, idents = 3)] = "Cap"
adult$CellType[adult_down_valve_cells] = "Down_Valve"
np$CellType[colnames(adult)] = adult$CellType

adult$CellType = factor(adult$CellType , levels = c("Ptx3_Cap", "Cap", "Col","Up_Valve", "Down_Valve"))
np$CellType = factor(np$CellType , levels = c("Ptx3_Cap", "Cap", "Col","Up_Valve", "Down_Valve"))

dplot_cols = c("#F75F55", "#ED8141", "#5bb300","#C77CFF", "#7997FF")


tiff("adult_nplp_dplot.tiff", width = 7, height = 7, units = "in", res = 300)
dimplot_fig(adult, pt.size = 5, cols = dplot_cols)
dev.off()

DimPlot(adult, group.by = "CellType", cols = dplot_cols)

vplot_genes = c("Cdh5", "Pecam1", "Prox1", "Lyve1", "Foxp2", "Efnb2", "Rgs3", "Plat", "Tgfa", "Mrc1", "Itih5", "Nrp2", "Ccl21a", "Reln", "Angpt2",
                "Aqp1", "Stab2", "Ptx3", "Fndc1", "Piezo2", "Adgrg3", "Mmrn1", "Fgl2", "Ptn", "Slco2b1", "Dlg1",
                "Cldn11", "Lama5", "Cd24a", "Esam", "F11r", "Itga9", "Gata2", "Gja4")
setwd("./Figures/")
dir.create("Vplots")
setwd("./Vplots")

Idents(adult) = "CellType"
for ( i in 1:length(vplot_genes)){
  Vplot_func(prefix = "Adult", object = adult, features = vplot_genes[i], assay = "RNA", cols = dplot_cols)
}

Idents(np) = "CellType"
tiff("Integrated_nplp_dplot.tiff", width = 7, height = 7, units = "in", res = 300)
dimplot_fig(np, pt.size = 5, cols = dplot_cols)
dev.off()

set.seed(1)
x = dimplot_fig(np, pt.size = 5, group.by = "sample", cols = c("#E2EBB5","#2380BA"))
x$data = x$data[sample(1:nrow(x$data)),]
tiff("Integrated_nplp_dplot_sample.tiff", width = 7, height = 7, units = "in", res = 300)
print(x)
dev.off()

make_barplot_by_samples(np, split = "sample", meta_to_use = "CellType", cols = dplot_cols) 

## Compare adult / aged by celltypes
mito_genes = grep(rownames(np), pattern = "^mt", value = T)

celltypes = names(table(np$CellType))
age_deg_list = list()

for (i in 1:length(celltypes)){
  celltype = celltypes[i]
  sub = subset(np, CellType %in% celltype)
  Idents(sub) = "sample"
  deg = FindMarkers(sub, ident.1 = "aged", logfc.threshold = 0.5, min.pct = 0.3, test.use = "MAST",
                    features = setdiff(rownames(np), c(genes.dissoc, ribo_genes, mito_genes))) # No only.pos to see number of overall degs
  # deg = deg %>% sig_markers(cutoff = )
  print(paste0(celltype, " age_degs_number:" , nrow(deg)))
  age_deg_list[[i]] = deg
  names(age_deg_list)[i] = celltype
}

for (i in 1:length(celltypes)){
  deg = age_deg_list[[i]]
  print(nrow(deg))
}

for (i in 1:length(celltypes)){
  deg = age_deg_list[[i]]
  deg_label = paste0(celltypes[i], "_age_degs.csv")
  deg %>% write.csv(deg_label)
}


## Adult heatmap
adult_markers = FindAllMarkers(adult, only.pos = T, test.use = "MAST", logfc.threshold = 0.3, min.pct = 0.3, features = setdiff(rownames(np), c(genes.dissoc, ribo_genes, mito_genes)))
adult_markers = adult_markers %>% sig_markers(cutoff = 0.05)
adult_markers_10 = adult_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
adult = ScaleData(adult, features = c(VariableFeatures(adult), adult_markers_10$gene),vars.to.regress = c("nCount_RNA", "percent.mito", "dissoc1","ribo1"))
DoHeatmap(adult, features = adult_markers_10$gene, disp.min = -1.5, disp.max = 1.5) + scale_fill_gradientn(colours = c("blue", "white", "red")) + NoLegend()
DoHeatmap(adult, features = adult_markers_10$gene, disp.min = -1.5, disp.max = 1.5) + scale_fill_gradientn(colours = c("blue", "white", "red"))

## 
Idents(np) = "sample"
age_markers = FindAllMarkers(np, only.pos = T, test.use = "MAST", logfc.threshold = 0.3, min.pct = 0.3, features = setdiff(rownames(np), c(genes.dissoc, ribo_genes, mito_genes))) %>%
              sig_markers(cutoff = 0.05)

age_markers$cluster %>% table
age_markers_30 = age_markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
np = ScaleData(np, features = c(VariableFeatures(np), age_markers_50$gene),vars.to.regress = c("nCount_RNA", "percent.mito", "dissoc1","ribo1"))
DoHeatmap(np, features = age_markers_30$gene, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colours = c("#558DB3","#F2F0AE","#B42E38")) + NoLegend()

age_markers %>% write.csv("age_markers_pseudobulk.csv")

age_markers
# genes to exclude in heatmap 
exclude_genes = c("Ddx3y", "Eif2s3y", "Malat1", "CT010467.1", "Gm10736", "Gm10095",  "Xist", "Gm1821", "2010111I01Rik")
exclude_genes %in% rownames(age_markers)

age_markers_heatmap = age_markers %>% filter(gene %ni% exclude_genes)
age_markers_heatmap %>% dim
age_markers_heatmap = age_markers_heatmap %>% group_by(cluster) %>% top_n(50, avg_log2FC)

age_markers_heatmap %>% dim
np = ScaleData(np, features = c(VariableFeatures(np), age_markers_heatmap$gene),vars.to.regress = c("nCount_RNA", "percent.mito", "dissoc1","ribo1"))
DoHeatmap(np, features = age_markers_heatmap$gene, disp.min = -1.5, disp.max = 1.5) + scale_fill_gradientn(colours = c("#558DB3","#F2F0AE","#B42E38"))
age_markers_heatmap %>% write.csv("heatmap_genes.csv")
save.image("Nphx_analysis.Rdata")

# Make adult vs aged violin plots 
load('Nphx_analysis.Rdata')
aged_vplot_genes = c("Fn1", "Mcl1", "Il33" ,"H2-Q6", "Tppp3", "Tgfbi", "Irf7", "Ifitm2", "Ifitm3", "Zbp1", "mhc_genes1")
mhc_genes = grep(rownames(np), pattern = "^H2-", value = T) # -> make as module score
mhc_genes %>% length

np = AddModuleScore(np, features = list(mhc_genes), name= "mhc_genes")

setwd("Figures/")
dir.create("aged_vplot") ;setwd("aged_vplot")

for ( i in 1:length(aged_vplot_genes)){
  Vplot_func(prefix = "adult_vs_aged", object = np, features = aged_vplot_genes[i], assay = "RNA", cols = c("#E2EBB5","#2380BA"))
}

# Irf7, Ifitm2, Ifitm3 bottom rows in heatmap
age_heatmap_genes = age_markers_heatmap$gene
age_heatmap_genes2 = c(age_heatmap_genes[-which(age_heatmap_genes %in% c("Irf7", "Ifitm2", "Ifitm3"))], c("Irf7", "Ifitm2", "Ifitm3"))
DoHeatmap(np, features = age_heatmap_genes2, disp.min = -1.5, disp.max = 1.5) + scale_fill_gradientn(colours = c("#558DB3","#F2F0AE","#B42E38")) + NoLegend()

# Preparation for GEO
library(Seurat)
np_list = SplitObject(np, split.by = "sample")
getdata(np_list[[1]]) %>% write.csv("Adult_nphx_normalized_expression_matrix.csv")
getdata(np_list[[2]]) %>% write.csv("Aged_nphx_normalized_expression_matrix.csv")



#
wilcox.test(x = np_list[[1]]$mhc_genes1, y = np_list[[2]]$mhc_genes1, paired = F)


