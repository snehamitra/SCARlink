## Generates  coassay_matrix.h5
## Update file paths of variables archr_out, seurat_out, and out_dir
## ArchR object should have non-binarized TileMatrix 
## The cell names in Seurat and ArchR objects need to be the same
## The script also saves the LSI matrix, KNN graph, and cell_info
 
library(ArchR)
library(Seurat)
library(rhdf5)
library(parallel)
library(BiocParallel)

### FUNCTIONS
get_gene_tile_matrix <- function(scatac.object, scrna.object, window_size=250000){

    # get tile matrix for variable genes
    geneRegions <- getGenes(scatac.object)
    seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
    geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]
    geneUpstream = window_size 
    geneDownstream = window_size
    geneRegions <- trim(extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream))
    geneRegions <- geneRegions[geneRegions$symbol %in% VariableFeatures(scrna.object), ]

    tm <- getMatrixFromProject(
        ArchRProj = scatac.object,
	useMatrix = "TileMatrix",
	binarize = FALSE)
    # reorder cells in tile matrix to match cell order in Seurat object
    tileGR <- rowData(tm)
    tile_size <- tileGR$start[2] - tileGR$start[1]
    print(paste("Tile size:", tile_size))
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

split_write_hdf5 <- function(dirname, selected.genes, tm.filtered)
{
    h5file <- tempfile(pattern="tmp_coassay_matrix", fileext=".h5", tmpdir=path.expand(dirname))
    h5createFile(h5file)    
    out <- lapply(selected.genes, function(x) write_hdf5(h5file, tm.filtered, x))
    return(h5file)
}

gather_tmp_h5 <- function(output_file, tmp_files) {

  # create the output file
  fid <- H5Fcreate(name = output_file)
  on.exit(H5Fclose(fid))
  
  ## iterate over the temp files and copy the named dataset into our new file
  for(i in tmp_files) {
    fid2 <- H5Fopen(i)
    genes <- h5ls(fid2, recursive=FALSE)$name
    for(gene in genes){
    	H5Ocopy(fid2, gene, h5loc_dest=fid, name_dest=gene)
    }
    H5Fclose(fid2)
  }  
}

write_hdf5 <- function(path, tile.matrix, gene)
{
    path <- path.expand(path) # protect against tilde's.

    group <- gene
    h5createGroup(path, group)

    cond <- rowData(tile.matrix)$symbol == gene

    h5write(as.data.frame(rowData(tile.matrix[cond, ])), file=path, name=paste0(group, "/tile_info"))

    # Saving matrix information.
    h5write(tile.matrix[cond, ]@assays@data$TileMatrix@x, file=path, name=paste0(group, "/data"))
    h5write(dim(tile.matrix[cond, ]@assays@data$TileMatrix), file=path, name=paste0(group, "/shape"))
    h5write(tile.matrix[cond, ]@assays@data$TileMatrix@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
    h5write(tile.matrix[cond, ]@assays@data$TileMatrix@p, file=path, name=paste0(group, "/indptr"))
    
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
    h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed
    h5write(x@p, file=path, name=paste0(group, "/indptr"))

    return(NULL)
}

write_files <- function(archr_out, seurat_out, out_dir, window_size, ncores, scale){
    ### Load Seurat and ArchR objects
    scatac.object <- loadArchRProject(archr_out)
    scrna.object <- readRDS(seurat_out)

    ### Normalize and scale values by  median or 10,000
    if(scale=='median'){
	scale.factor = median(colSums(scrna.object[['RNA']]@counts))
    }
    else{
	scale.factor = 10000
    }
    
    scrna.object <- NormalizeData(scrna.object, normalization.method = "RC", scale.factor=scale.factor)

    ### Create output directory if it doesn't exist
    if(!file.exists(out_dir)){
	dir.create(out_dir, showWarnings = TRUE, recursive=TRUE)
    }

    ### Make sure same cells are in both scrna.object and scatac.object
    scatac.object <- scatac.object[scatac.object$cellNames %in% colnames(scrna.object), ]
    scrna.object <- scrna.object[, colnames(scrna.object) %in% scatac.object$cellNames, ]

    n_cells <- dim(scrna.object)[2]

    scatac.object <- scatac.object[match(colnames(scrna.object), scatac.object$cellNames), ]
    if(n_cells==0){
	stop("No cells in common between scATAC and scRNA objects. Please check the cell names")
    }
    print(paste(n_cells, "cells in common between scATAC and scRNA objects"))

    tm.filtered <- get_gene_tile_matrix(scatac.object, scrna.object, window_size=window_size)
    cell.info <- cbind(scrna.object@meta.data, as.data.frame(scatac.object@cellColData))
    cell.info <- cell.info[!duplicated(as.list(cell.info))]
    selected.genes <- unique(rowData(tm.filtered)$symbol)

    if(ncores > 1){
    	bpparam <- MulticoreParam(workers=ncores) 
    }
    else{
	bpparam <- BiocParallel::SerialParam()
    }

    # create chunks
    max_ix <- min(length(selected.genes), max(trunc(length(selected.genes)/50), 50))
    
    tmp_files <- bplapply(1:max_ix, FUN = function(i) {
    	   selected.genes.subset <- selected.genes[seq(i, length(selected.genes), max_ix)]

    	   split_write_hdf5(out_dir, selected.genes.subset, tm.filtered[rowData(tm.filtered)$symbol %in% selected.genes.subset, ])
           }, 
           BPPARAM=bpparam) 
    
    o <- h5closeAll()	    

    tmp_files <- as.vector(unlist(tmp_files))

    h5file <- paste(out_dir, 'coassay_matrix.h5', sep = '/')
    
    # h5createFile(h5file)

    # gathering output in one file
    gather_tmp_h5(h5file, tmp_files)

    file.remove(tmp_files)

    write_hdf5_rna(h5file, scrna.object, selected.genes)

    # cell.info['cell_name'] = unlist(as.vector(lapply(rownames(cell.info), function(x) strsplit(x, '')[[1]][2])))
    cell.info['cell_name'] = rownames(cell.info)
    h5write(cell.info, file=h5file, name="cell_info")


    ## Save LSI, cell_info, and KNN as csv files
    lsi <- getReducedDims(ArchRProj = scatac.object, reducedDims = "IterativeLSI")
    lsi <- lsi[match(colnames(scrna.object), rownames(lsi)), ]

    write.table(lsi, file = paste(out_dir, 'scatac_LSI.csv', sep=""), row.names = FALSE, quote = FALSE, sep='\t')

    write.table(cell.info, file=paste(out_dir, 'cell_info.txt', sep=""), row.names=FALSE, quote=FALSE, sep='\t')

    ## Create KNN graph and save
    nbrs <- nabor::knn(lsi, lsi, k = 50)$nn.idx
    dists <- nabor::knn(lsi, lsi, k = 50)$nn.dists
    nbrs <- nbrs - 1
    write.table(nbrs, file=paste(out_dir, "knn-50-scatac.csv", sep=""), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

    ## Save variable genes list
    # scrna.object <- FindVariableFeatures(scrna.object, selection.method = "vst", nfeatures = 5000)	
    write.table(VariableFeatures(scrna.object), file=paste(out_dir, "hvg.txt", sep=""), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}