# PASBench


## Introduction
Single-cell RNA sequencing (scRNA-seq) analysis enables researchers to uncover more refined and novel cell clusters, which have greatly advanced our understanding of cellular states. There were many state-of-art computational tools developed for clustering cells, identifying marker genes, and visualizing scRNA-seq data. However, biological interpretation of the clustering results remains a big challenge. Pathway activity scores (PASs) analysis has been applied to transform the gene-level data into explainable gene sets representing biological processes or pathways to uncover the potential mechanism of cell heterogeneity. To the best of our knowledge, there were no systematic benchmark studies to evaluate the performance of these unsupervised PAS transformation algorithms.

This reposity contain seven PAS tools and provide a shiny APP to explore results

## Installation
```
devtools::install_github("ZhangBuDiu/PASBench")
```

## Running with an example

### Step1: 
calculating pathway activity score. 

User could select species and pathway database or assign customized pathways in GMT format
```
data("counts")
pas_score = calculate_PAS(counts, tool='pagoda2',species='human',pathway='kegg')
# pas_score = calculate_PAS(counts, gmt_file = 'path/to/file.gmt')
```
### Step2: 
preparing visualization object
```
# data("vis_oj.rda")
vis_oj = prepare_vis(pas_score)
```
### Setp3: 
visulizing through an interactive webpage. 

The input is a Seurat object which obtained from step2 or user customized.
```
PAS_vis(vis_oj)
```

## Results


PASBench provides a user-interactive and flexible webpage for visulization and exploration


#### Panel1: clustering analysis


<div align=center style="border:5px solid #000"><img  src="https://github.com/ZhangBuDiu/PASBench/blob/master/pic/clustering.png"/> </div>


#### Panel2: differential analysis (click GO!)


<div align=center style="border:5px solid #000"><img  src="https://github.com/ZhangBuDiu/PASBench/blob/master/pic/differential.png"/> </div>


#### Panel2: trajecoty analysis (click GO!)


<div align=center style="border:5px solid #000"><img  src="https://github.com/ZhangBuDiu/PASBench/blob/master/pic/trajectory.png"/> </div>

## Detials
### pathway
#### `human`
| Name | Detials  | Number of gene sets |
| - | :-: | -: |
| hallmarker | Hallmark gene sets | 50 |
| CGP | genetic and chemical perturbations | 3297 |
|biocarta | BioCarta pathway database | 289|
|kegg | KEGG pathway database | 186|
|PID | PID pathway database | 196|
|reactome | Reactome pathway database | 1532|
|TFT | transcriptional factor targets | 1137|
|CGN | cancer gene neighborhoods | 427|
|CM | cancer models | 431|
|GO.bp | GO biological process | 7530|
|GO.cc | Co cellular Component | 999|
|GO.mf | GO molecular fucntion | 1663|
|OncoG | oncogenic signatures | 189|
|Immu | immunologic signatures | 4872|
|panther | protein annotation through evolutionary relationship | 94|
|humancyc | human metabonomics | 127|


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



