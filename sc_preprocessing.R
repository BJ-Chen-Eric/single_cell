seu_preprocess <- function(seurat, specy='human',pre_process=T, mt_p=10 , max=''
                           , decomposition=F, dim=1:20, res=0.05, nHVG=3000)  {
  t1 <- Sys.time()
  library(Seurat)
  library(dplyr)
  if(isTRUE(pre_process))  {
    mito <- '^MT-'
    if(specy=='mouse')  {mito <- '^mt-'}
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    # A <- as.data.frame(seurat@assays$RNA@counts@Dimnames[[1]])
    # grep(pattern = '^IGHA', x = seurat@assays$RNA@counts@Dimnames[[1]], value = T)
    print(VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    # plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
    # plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    # plot1 + plot2
    if(identical('', max))  {
      max <- min(boxplot.stats(seurat$nCount_RNA)$out)
    }
    seurat <- subset(seurat, subset = nFeature_RNA > 200 & nCount_RNA <max & percent.mt < mt_p)
    dim(seurat)
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    # seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    # top10 <- head(VariableFeatures(seurat), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(seurat)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # plot1 + plot2
    seurat <- ScaleData(seurat, features = rownames(seurat))
    t2 <- Sys.time()
    print(paste('preprocess time is ', t2-t1, sep = ''))
  }
  
  A <- as.data.frame(rownames(seurat))
  row <- ''
  if(identical(specy, 'human'))  {
    gene_gift <- c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
      ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
      , mito,'^Z[0-9]+')
  }else{
    gene_gift <- c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
                   ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
                   , '^mt-','^Z[0-9]+')
  }
  for(i in gene_gift)  { #, '^IG[HKL]'
    row <- c(row, grep(rownames(seurat), pattern = i))
  }
  row <- row[-c(1)]
  seurat <- seurat[-c(as.numeric(row)),]
  
  # VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
  # DimPlot(seurat, reduction = "pca")
  # 
  # seurat <- JackStraw(seurat, num.replicate = 100)
  # seurat <- ScoreJackStraw(seurat, dims = 1:20)
  # ElbowPlot(seurat)
  if(decomposition==T)  {
    t1 <- Sys.time()
    # seurat <- seurat[VariableFeatures(all), ]
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = nHVG)
    seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
    t2 <- Sys.time()
    print(paste('PCA time is ', t2-t1, sep = ''))
    
    
    t1 <- Sys.time()
    seurat <- FindNeighbors(seurat, dims = dim)
    seurat <- FindClusters(seurat, resolution = res)
    
    seurat <- RunTSNE(seurat, dims = dim, check_duplicates = FALSE)
    seurat <- RunUMAP(seurat, dims = dim)
    t2 <- Sys.time()
    print(paste('de-dimension time is ', t2-t1, sep = ''))
  }
  return(seurat)
}
