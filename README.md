# Benchmarking single cell transcriptome matching methods for incremental growth of cell atlases

## Abstract
The advancement of single cell technologies has driven significant progress in constructing a multiscale, pan-organ Human Reference Atlas for healthy human cells. Many multi-faceted cell atlases for different organs, species, and diseases now exist, though challenges remain in harmonizing cell types and unifying nomenclature among respective cell atlases. Multiple machine learning and artificial intelligence methods, including models pre-trained on large-scale cell atlas datasets, are publicly available for single cell community users to computationally map their cell clusters to the cell atlases. In this study, seven computational tools for cell type matching and label transfer – Azimuth, CellTypist, CellHint, FR-Match, scArches, scPred, and singleR – were benchmarked in ten organ systems.  Using the healthy human lung as an exemplar organ, when matching two well-annotated atlases – the Human Lung Cell Atlas (HLCA) and the LungMAP Single-Cell Reference (CellRef), large variations in the matching accuracy were observed, especially in rare cell types, underlining the need for a consensus strategy using a selective set of computational methods.  In the subsequent meta-analysis, top benchmarked methods were used to incrementally integrate 61 cell types from HLCA and 48 from CellRef, resulting in a meta-atlas of 41 matched, 20 HLCA-, and 7 CellRef-specific cell types.  Similar approach revealed 25 matched cell types existed in two independent kidney atlases.  Generalizability of the benchmarking performances were further demonstrated in a variety of organ systems. In summary, this study reveals the complementing strengths of benchmarked methods and presents a framework for incremental growth of cell types in cell atlases.

<div align="center">
  <img width = "55%" src="https://github.com/jjoycehu/manuscript/blob/main/figures/Figure1.png">
</div>

**Figure 1: Study Overview**  
**(A)** Unsupervised and supervised cell type matching between query and reference datasets.   
**(B)** Cell-based and cluster-based cell type matching methods.   
**(C)** A sustainable workflow for continued integration of reference atlases requires a strategy that identifies new
cell types from new datasets to be integrated to the reference atlas while preserving the cell
type memberships of existing cells, to support the incremental growth of the knowledgebase.

## Directory
* `Supplementary_results/`: contains supplementary results for the study.
* `tutorials/`: contains code used in this study documented in tutorial notebooks.
* `figures/`: contains scripts for all figures in the main manuscript.

## Citation
* Preprint: [https://pmc.ncbi.nlm.nih.gov/articles/PMC12190912/](https://www.biorxiv.org/content/10.1101/2025.04.10.648034)
