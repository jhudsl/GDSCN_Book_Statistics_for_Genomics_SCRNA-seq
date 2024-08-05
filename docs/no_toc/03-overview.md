# OVERVIEW

## Scope

There are many excellent programs available for scRNA-seq analysis. In this lesson, we will be using the `scran` and `scater` packages available from Bioconductor. We do not cover .fastq processing, nor do we explore the [Seurat](https://satijalab.org/seurat/) and [scanpy](https://scanpy.readthedocs.io/en/stable/index.html) pipelines.

## Bioconductor reference book

The lessons and activities are adapted from the [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/) book by Robert Amezquita, Aaron Lun, Stephanie Hicks, and Raphael Gottardo.

## SingleCellExperiment

Bioconductor stores data from single-cell experiments in a special class (the `SingleCellExperiment` class). The rows in this data class represent features of interest (such as genes, transcripts, or genomic regions), while the columns represent individual cells. This data class is the main data structure used in the `scran` and `scater` packages. In addition to storing count data, it can also store the results of dimensionality reduction methods as well as data for additional, alternative features (like spike-in transcripts).

`SingleCellExperiment` objects are created using a command of the same name.


```r
# AnVIL::install(c("SingleCellExperiment"))
library(SingleCellExperiment)

#create an empty data frame to use as an example 
counts <- matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10)
pretend.cell.labels <- sample(letters, ncol(counts), replace = TRUE)
counts <- matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10)
pretend.cell.labels <- sample(letters, ncol(counts), replace = TRUE)
pretend.gene.lengths <- sample(10000, nrow(counts))

#create SingleCellExperiment object
sce <- SingleCellExperiment(list(counts = counts),
    colData=DataFrame(label = pretend.cell.labels),
    rowData=DataFrame(length = pretend.gene.lengths),
    metadata=list(study = "GSE111111"))
```

## Obtaining Data

For these exercises, we are using a dataset of cells from a mouse brain from [Zeisel et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25700174/). This dataset contains about 3000 brain cells of a wide variety of types (such as neurons, glia, and oligodendrocytes). The library was prepared using a UMI-based (unique molecular identifiers) protocol, and the expression data you are working with is UMI counts that have been mapped to each gene. You can learn more about UMI protocols and why they are useful [here](https://dnatech.genomecenter.ucdavis.edu/faqs/what-are-umis-and-why-are-they-used-in-high-throughput-sequencing/).

This dataset has already had low-quality cells removed prior to publication, but we will still run QC steps to verify that the remaining cells are of good quality.

The dataset is stored as a `SingleCellExperiment` object in the `scRNAseq` package. Here we retrieve the dataset, as well as do a little bit of extra formatting to merge together redundant rows (which are the result of alternative genomic coordinates for the same gene). 


```r
# AnVIL::install(c("scRNAseq", "scater"))
library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

#library(scater)
#sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))
```
