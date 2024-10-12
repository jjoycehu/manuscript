# Seurat Tutorial
setwd("/Users/joycehu/Dev/hcla")
rm(list = ls())

# libraries
library(Seurat)
library(Azimuth)
library(glue)
library(dplyr)
library(hdf5r) # if conversion is needed

# path to reference dataset; refer to these two links and download the files for "Human - Lung v2 (HLCA)"
# https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format#required-files-for-scrna-seq-queries
# https://azimuth.hubmapconsortium.org/

ref_path <- 'seurat/reference'


###########################################
##### if your query dataset is .rds: #####
###########################################

# read in query dataset
query_data <- readRDS('data/hlca_core.rds')
annotations <- RunAzimuth(query=query_data, reference=ref_path, annotation.levels='ann_finest_level')


###########################################
##### if your query dataset is .h5ad: #####
###########################################

# set path for both query dataset
# before conversion, remove all columns from .obs besides the original label column
query_path <- 'data/cellref/cellref_set1.h5ad' 

# you can try running Azimuth without any modifications, but if error arises, you have to manually load the query object
annotations <- RunAzimuth(query=query_path, reference=ref_path)

######### if you get an error: ##############

# manually convert .h5ad object to Seurat object 
# this function is from Seurat but some changes were made to accomodate the specific datasets used
LoadH5AD <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Loading H5AD files requires hdf5r", call. = FALSE)
  }
  adata <- hdf5r::H5File$new(filename = path, mode = 'r')
  on.exit(expr = adata$close_all())
  Exists <- function(name) {
    name <- unlist(x = strsplit(x = name[1], split = '/', fixed = TRUE))
    hpath <- character(length = 1L)
    exists <- TRUE
    for (i in seq_along(along.with = name)) {
      hpath <- paste(hpath, name[i], sep = '/')
      exists <- adata$exists(name = hpath)
      if (isFALSE(x = exists)) {
        break
      }
    }
    return(exists)
  }
  GetIndex <- function(md) {
    return(
      if (adata[[md]]$attr_exists(attr_name = '_index')) {
        h5attr(x = adata[[md]], which = '_index')
      } else if (adata[[md]]$exists(name = '_index')) {
        '_index'
      } else if (adata[[md]]$exists(name = 'index')) {
        'index'
      } else {
        stop("Cannot find the rownames for '", md, "'", call. = FALSE)
      }
    )
  }
  GetRowNames <- function(md) {
    return(adata[[md]][[GetIndex(md = md)]][])
  }
  LoadMetadata <- function(md) {
    factor.cols <- if (adata[[md]]$exists(name = '__categories')) {
      names(x = adata[[md]][['__categories']])
    } else {
      NULL
    }
    index <- GetIndex(md = md)
    col.names <- names(x = adata[[md]])
    if (adata[[md]]$attr_exists(attr_name = 'column-order')) {
      tryCatch(
        expr = {
          col.order <- hdf5r::h5attr(x = adata[[md]], which = 'column-order')
          col.names <- c(
            intersect(x = col.order, y = col.names),
            setdiff(x = col.names, y = col.order)
          )
        },
        error = function(...) {
          return(invisible(x = NULL))
        }
      )
    }
    col.names <- col.names[!col.names %in% c('__categories', index)]
    df <- sapply(
      X = col.names,
      FUN = function(i) {
        x <- adata[[md]][[i]][['codes']][]
        if (i %in% factor.cols) {
          x <- factor(x = x, levels = adata[[md]][['categories']][[i]][])
        }
        return(x)
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    return(as.data.frame(x = df, row.names = GetRowNames(md = md)))
  }
  if (Exists(name = 'raw/X')) {
    md <- 'raw/var'
    x <- adata[['raw/X']]
  } else if (Exists(name = 'X')) {
    md <- 'var'
    x <- adata[['X']]
  } else {
    stop("Cannot find counts matrix", call. = FALSE)
  }
  # check different possible attributes to try and get matrix shape
  if (isTRUE(x$attr_exists(attr_name = 'h5sparse_shape'))) {
    mtx.shape <- h5attr(x, 'h5sparse_shape')
  } else if (isTRUE(x$attr_exists(attr_name = 'shape'))) {
    mtx.shape <- h5attr(x, 'shape')
  } else {
    warning("Could not determine matrix shape")
  }
  # check different attributes to try and deduce matrix type
  if (isTRUE(x$attr_exists(attr_name = 'h5sparse_format'))) {
    mtx.type <- h5attr(x, 'h5sparse_format')
  } else if (isTRUE(x$attr_exists(attr_name = 'encoding-type'))) {
    mtx.type <- substr(h5attr(x, 'encoding-type'), 0, 3)
  } else {
    mtx.type <- 'csr' # assume matrix is csr
    warning("Could not determine matrix format")
  }
  if (mtx.type != 'csr') {
    p <- as.integer(x[['indptr']][])
    i <- as.integer(x[['indices']][])
    data <- as.double(x[['data']][])
    # csc -> csr
    converted.mtx <- csc_tocsr(
      n_row = as.integer(mtx.shape[1]),
      n_col = as.integer(mtx.shape[2]),
      Ap = p,
      Ai = i,
      Ax = data
    )
    # csr -> dgC
    counts <- new(
      Class = 'dgCMatrix',
      p = converted.mtx$p,
      i = converted.mtx$i,
      x = converted.mtx$x,
      Dim = c(mtx.shape[2], mtx.shape[1])
    )
  } else {
    # x must be a CSR matrix
    counts <- as.matrix(x = x)
  }
  metadata <- LoadMetadata(md = 'obs') # gather additional cell-level metadata
  rownames <- GetRowNames(md = md)
  colnames <- rownames(metadata)
  rownames(x = counts) <- rownames
  colnames(x = counts) <- colnames
  options(Seurat.object.assay.calcn = TRUE)
  object <- CreateSeuratObject(counts = counts)
  if (ncol(x = metadata)) {
    object <- AddMetaData(object = object, metadata = metadata)
  }
  object <- subset(object, subset = nCount_RNA > 0)
  return(object)
}

# load query object using the function
query_obj <- LoadH5AD(path=query_path)

# run azimuth
annotations <- RunAzimuth(query=query_obj, reference=ref_path)

###########################################
##########   save output   ################
###########################################

# save output
file_name <- glue("seurat/output/prediction.rds")
saveRDS(annotations, file = file_name)

# save data frame
df <- annotations@meta.data
df <- select(df, predicted.ann_finest_level, ann_finest_level, predicted.ann_finest_level.score)

filename = glue("seurat/output/prediction.csv")
write.csv(df, file = filename, row.names = TRUE)






