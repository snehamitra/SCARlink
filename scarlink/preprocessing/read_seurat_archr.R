## Generates  coassay_matrix.h5
## Update file paths of variables archr_out, seurat_out, and out_dir
## ArchR object should have TileMatrix (tile size: 500)
## The cell names in Seurat and ArchR objects need to be the same
## The script also saves the LSI matrix, KNN graph, and cell_info
 
library(ArchR)
library(Seurat)
library(oncomix)
library(rhdf5)
library(parallel)

### FUNCTIONS
get_gene_tile_matrix <- function(scatac.object, scrna.object){

    # tile size set to 500
    tile_size <- 500

    # get tile matrix for variable genes
    geneRegions <- getGenes(scatac.object)
    seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
    geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]
    geneUpstream = 250000 
    geneDownstream = 250000
    geneRegions <- trim(extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream))
    geneRegions <- geneRegions[geneRegions$symbol %in% scrna.object$RNA@var.features, ]


    tm <- getMatrixFromProject(
        ArchRProj = scatac.object,
	useMatrix = "TileMatrix",
	binarize = FALSE)
    # reorder cells in tile matrix to match cell order in Seurat object
    tileGR <- rowData(tm)
    tileGR$end <- tileGR$start + tile_size
    tileGR <- makeGRangesFromDataFrame(tileGR)
    tile.overlaps <- findOverlaps(tileGR, geneRegions)

    tm.filtered <- tm[queryHits(tile.overlaps), ]
    rowData(tm.filtered)$strand <- geneRegions[subjectHits(tile.overlaps)]$strand
    rowData(tm.filtered)$geneRange <- geneRegions[subjectHits(tile.overlaps)]$ranges
    rowData(tm.filtered)$gene_id <- geneRegions[subjectHits(tile.overlaps)]$gene_id
    rowData(tm.filtered)$symbol <- geneRegions[subjectHits(tile.overlaps)]$symbol

    tm.filtered <- tm.filtered[, match(colnames(scrna.object), colnames(tm))]
    rowData(tm.filtered)$end <- rowData(tm.filtered)$start + tile_size
    return(tm.filtered)
}

write_hdf5 <- function(path, tile.matrix, gene)
{
    print(gene)
    path <- path.expand(path) # protect against tilde's.
    # h5createFile(path)
    group <- gene
    h5createGroup(path, group)

    tile.matrix.gene <- tile.matrix[rowData(tile.matrix)$symbol == gene, ]
    h5write(as.data.frame(rowData(tile.matrix.gene)), file=path, name=paste0(group, "/tile_info"))

    # Saving matrix information.
    x <- tile.matrix.gene@assays@data$TileMatrix
    h5write(x@x, file=path, name=paste0(group, "/data"))
    h5write(dim(x), file=path, name=paste0(group, "/shape"))
    h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
    h5write(x@p, file=path, name=paste0(group, "/indptr"))

    return(NULL)
}

write_hdf5_rna <- function(path, scrna, genes)
{
    path <- path.expand(path) # protect against tilde's.
    # h5createFile(path)
    group <- "gene_expression"
    h5createGroup(path, group)

    df = as.data.frame(scrna$RNA@meta.features[genes, ])
    df['gene_name'] = rownames(df)
    h5write(df, file=path, name=paste0(group, "/gene_info"))

    # Saving matrix information.
    x <- scrna$RNA@data[genes, ]
    h5write(x@x, file=path, name=paste0(group, "/data"))
    h5write(dim(x), file=path, name=paste0(group, "/shape"))
    h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
    h5write(x@p, file=path, name=paste0(group, "/indptr"))

    return(NULL)
}

write_files <- function(archr_out, seurat_out, out_dir){
    ### Load Seurat and ArchR objects
    scatac.object <- loadArchRProject(archr_out)
    scrna.object <- readRDS(seurat_out)

    ### Create output directory if it doesn't exist
    if(!file.exists(out_dir)){
	dir.create(out_dir, showWarnings = TRUE, recursive=TRUE)
    }

    scatac.object <- scatac.object[match(colnames(scrna.object), scatac.object$cellNames), ]

    tm.filtered <- get_gene_tile_matrix(scatac.object, scrna.object)
    cell.info <- cbind(scrna.object@meta.data, as.data.frame(scatac.object@cellColData))
    cell.info <- cell.info[!duplicated(as.list(cell.info))]
    selected.genes <- unique(rowData(tm.filtered)$symbol)

    h5file <- paste(out_dir, 'coassay_matrix.h5', sep = '/')
    h5createFile(h5file)
    out <- lapply(selected.genes, function(x) write_hdf5(h5file, tm.filtered, x))
    write_hdf5_rna(h5file, scrna.object, selected.genes)

    cell.info['cell_name'] = unlist(as.vector(lapply(rownames(cell.info), function(x) strsplit(x, '')[[1]][2])))
    h5write(cell.info, file=h5file, name="cell_info")


    ## Save LSI, cell_info, and KNN as csv files
    lsi <- getReducedDims(ArchRProj = scatac.object, reducedDims = "IterativeLSI")
    write.table(lsi, file = paste(out_dir, 'scatac_LSI.csv', sep=""), row.names = FALSE, quote = FALSE, sep='\t')

    write.table(cell.info, file=paste(out_dir, 'cell_info.txt', sep=""), row.names=FALSE, quote=FALSE, sep='\t')

    ## Create KNN graph and save
    nbrs <- nabor::knn(lsi, lsi, k = 50)$nn.idx
    dists <- nabor::knn(lsi, lsi, k = 50)$nn.dists
    nbrs <- nbrs - 1
    write.table(nbrs, file=paste(out_dir, "knn-50-scatac.csv", sep=""), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
}