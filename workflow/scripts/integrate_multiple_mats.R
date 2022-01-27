suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(argparse))
suppressMessages(library(glmGamPoi))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(miQC))
suppressMessages(library(flexmix))
suppressMessages(library(SCNPrep))
set.seed(1)

parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--data',
                    nargs = "+",
                    help = 'Path to count matrixes')
parser$add_argument('--annot',
                    type = "character",
                    help = 'annotation path')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')
parser$add_argument('--sample_id',
                    type = "character",
                    help = 'Path to output directory')
## SET VARIABLES

args <- parser$parse_args()

## FUNCTIONS


add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mito_log10']] <- log10(data[['percent.mito']] + 1)
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']] + 1)
  data[['scaled_mito']] <- scale(percent.mito)
  data[['scaled_nCount_RNA']] <- scale(data[['nCount_RNA_log10']])
  attr(data$scaled_nCount_RNA, "scaled:center") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:scale") <- NULL
  attr(data$scaled_mito, "scaled:center") <- NULL
  attr(data$scaled_mito, "scaled:scale") <- NULL
  data
}

## GATHERING DATA TOGETHER


options(future.globals.maxSize = 20000 * 1024^2)


get_object <- function(counts, annot) {
  data <- fread(counts) %>%
    tibble::column_to_rownames('V1')
  obj <- CreateSeuratObject(data)
  df <- fread(annot) %>%
    filter(cellID %in% colnames(obj)) %>% 
    dplyr::select(-cellID) %>% 
    as.data.frame()
  obj@meta.data <- cbind(obj@meta.data, df)
  obj <- add_metadata(obj)
  Seurat::SplitObject(obj, split.by = 'patient')
}

print(args)
# 
# whole <- unlist(lapply(args$data, function(x) get_object(x, args$annot)))
# whole <- Reduce(function(x,y) merge(x,y,add.cell.ids = c(x@project.name,y@project.name)) , whole)
# whole <- SplitObject(whole, split.by = 'patient')


setwd(args$out_dir)

## Number of cells before
# 
# cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])
# whole <- whole[cells.before > 100]
# cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])


## NORMALIZATION

# whole <- sapply(whole, function(x) SCTransform(
#   x,
#   ncells=min(100000, ncol(x)),
#   vars.to.regress = c("percent.mito"),
#   method = "glmGamPoi",
#   verbose = T,
#   conserve.memory = T
# ))

## INTEGRATION
# 
# whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)
# 
# whole <- lapply(X = whole, FUN = function(x) {
#   x <- RunPCA(x, features = whole.features)
# })
# 
# whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features,
#                             verbose = FALSE)
# whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT",
#                                         anchor.features = whole.features, verbose = FALSE, reduction = 'rpca')
# whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
#                                   verbose = FALSE)
# 
# ## PCA
# gc()
# 
# whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)
# 
# ## UMAP
# 
# whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)
# 
# 
# ## CLUSTERING
# 
# whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
# whole.integrated <- FindClusters(object = whole.integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1))
# 
# ## SAVING: DATASET
# 
# save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = "object.RData")
# 
# ## FINDING ANS SAVING MARKERS
# 
# analyze_object <- function(object, ident) {
#   Idents(object) <- object[[ident]]
#   if (length(levels(object)) == 1) {
#     return(message(sprintf('%s: since only one cluster was identified, markers can not be found', ident)))
#   }
#   out_dir <- paste0('markers/', ident)
#   dir.create(out_dir, recursive = T)
#   whole.markers <- FindAllMarkers(object = object,
#                                   assay='SCT',
#                                   only.pos = TRUE,
#                                   min.pct = 0.10,
#                                   test.use = 'wilcox',
#                                   max.cells.per.ident = 3e3,
#                                   random.seed = 42)
#   write.table(whole.markers, paste(out_dir, "markers.tsv", sep = '/'), sep="\t", quote=F, row.names=F)
# }
# 
# sapply(c('meta.cluster', grep('snn_res', colnames(whole.integrated@meta.data), value = T)),
#        function(ident) analyze_object(object = whole.integrated, ident = ident))
# 
# 
# ## Number of cells after
# 
# cells.after <- sapply(whole, function(x) length(colnames(x = x)))
# cells.diff <- cells.before-cells.after
# rbind(cells.before, cells.after, cells.diff)

## Conversion

load('object.RData')

markers_files <- list.files(pattern = 'markers.tsv', recursive = T, full.names = T)
markers <- lapply(markers_files, function(x) fread(x) %>% dplyr::rename(avg_logFC = 'avg_log2FC'))
names(markers) <- gsub('/|markers', '', stringr::str_extract(markers_files, 'markers/*.*/'))

whole.integrated@meta.data <- whole.integrated@meta.data %>%
  mutate_if(is.character, as.factor) %>% 
  mutate(clusters_labels = meta.cluster) %>% 
  dplyr::select(-c('meta.cluster'))

migrateSeuratObject(whole.integrated,
                    assay="SCT",
                    species='hs',
                    outdir = args$sample_id,
                    public = T,
                    curated = T,
                    markers = markers,
                    generateMarkers = F,
                    generateGMTS = F,
                    name=args$sample_id,
                    token=args$sample_id,
                    description=sprintf('Pan-cancer single-cell landscape of tumor-infiltrating T cells [%s].', args$sample_id),
                    link='https://www.science.org/doi/10.1126/science.abe6474')