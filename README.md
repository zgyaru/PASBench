# PASBench


## Introduction
Single-cell RNA sequencing (scRNA-seq) analysis enables researchers to uncover more refined and novel cell clusters, which have greatly advanced our understanding of cellular states. There were many state-of-art computational tools developed for clustering cells, identifying marker genes, and visualizing scRNA-seq data. However, biological interpretation of the clustering results remains a big challenge. Pathway activity scores (PASs) analysis has been applied to transform the gene-level data into explainable gene sets representing biological processes or pathways to uncover the potential mechanism of cell heterogeneity. To the best of our knowledge, there were no systematic benchmark studies to evaluate the performance of these unsupervised PAS transformation algorithms.

This reposity contain seven PAS tools and provide a shiny APP to explore results

### Installation
```
devtools::install_github("ZhangBuDiu/PASBench")
```

### Running with an example
```
counts = PASBench::load_counts()
pas_score = PASBench::calculate_PAS(counts, tool='pagoda2',species='human',pathway='kegg')
vis_oj = PASBench::prepare_vis()
PASBench::PAS_vis(vis_oj)
```

### Results
a user-interactive visulation page
