# PASBench


## Introduction
Single-cell RNA sequencing (scRNA-seq) analysis enables researchers to uncover more refined and novel cell clusters, which have greatly advanced our understanding of cellular states. There were many state-of-art computational tools developed for clustering cells, identifying marker genes, and visualizing scRNA-seq data. However, biological interpretation of the clustering results remains a big challenge. Pathway activity scores (PASs) analysis has been applied to transform the gene-level data into explainable gene sets representing biological processes or pathways to uncover the potential mechanism of cell heterogeneity. To the best of our knowledge, there were no systematic benchmark studies to evaluate the performance of these unsupervised PAS transformation algorithms.

This reposity contain seven PAS tools and evaluation metrics

### Evaluation scheme

<div align=center style="border:5px solid #000"><img  src="https://github.com/ZhangBuDiu/PASBench/blob/master/pic/workflow.jpg"/> </div>




## Detials
### pathway
#### `human`
| Name | Detials  | Number of gene sets |
| - | :-: | -: |
| hallmarker | Hallmark gene sets | 50 |
| CGP | genetic and chemical perturbations | 3297 |
|biocarta | BioCarta pathway database | 289 |
|kegg | KEGG pathway database | 186 |
|PID | PID pathway database | 196 |
|reactome | Reactome pathway database | 1532 |
|TFT | transcriptional factor targets | 1137 |
|CGN | cancer gene neighborhoods | 427 |
|CM | cancer models | 431|
|GO.bp | GO biological process | 7530 |
|GO.cc | Co cellular Component | 999 |
|GO.mf | GO molecular fucntion | 1663|
|OncoG | oncogenic signatures | 189 |
|Immu | immunologic signatures | 4872 |
|panther | protein annotation through evolutionary relationship | 94 |
|humancyc | human metabonomics | 127 |


#### `mouse`
|Name | Detials  | Number of gene sets|
|- | :-: | -: |
|kegg | KEGG pathway database | 259|
|PID | PID pathway database | 193
|panther | protein annotation through evolutionary relationship | 151|
|mousecyc | mouse metabonomics | 321|
|biocarta | BioCarta pathway database | 176|
|reactome | Reactome pathway database | 4342|
|TFT | transcriptional factor targets | 373|
|GO.bp | GO biological process | 8203|
|GO.cc | Co cellular Component | 1082|
|GO.mf | GO molecular fucntion | 3240|

