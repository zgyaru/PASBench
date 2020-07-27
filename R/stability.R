library(reticulate)
use_python('/share/apps/anaconda3/bin/python3',required = T)

source('./tools.R')
source('./evaluation.R')
source('./utils.R')

suppressPackageStartupMessages(library(VISION))
suppressPackageStartupMessages(library(pagoda2))
suppressPackageStartupMessages(library(scde))
suppressPackageStartupMessages(library(AUCell))
suppressPackageStartupMessages(library(slalom))
suppressPackageStartupMessages(library(rRoma))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(SCnorm))
suppressPackageStartupMessages(library("optparse"))

#library(scDblFinder)
# clustering
suppressPackageStartupMessages(library(SC3))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(HDF5Array))


### best choice of preprocessing
tools = c("plage","pagoda2","GSVA",'AUCell','ssGSEA','zscore','Vision')
preprocess = vector(mode="list",length=length(tools))
names(preprocess) = tools
preprocess$plage = c('none','log')
preprocess$pagoda2 = c('none','none')
preprocess$AUCell = c('none','none')
preprocess$GSVA = c('freq','log')
preprocess$ssGSEA = c('freq','none')
preprocess$zscore = c('freq','log')
preprocess$Vision = c('freq','none')




#data_path = './dataSets/human/lung/li.rds'
#data_path = './dataSets/human/blood/zheng_5k.rds'

#out_folder = './results/stability/li/'
n_cores = 5
times = 5
ratio = 1:9
ratio = ratio/10

args = commandArgs(T)
dropType = args[1]
start = as.integer(args[2])
end = as.integer(args[3])
species = args[4]
organ = args[5]
data = args[6]
if(species == 'human'){
  gSets_path = './gmtFiles/human/c2.kegg.gmt'
}else{
  gSets_path = './gmtFiles/mouse/kegg.gmt'
}


data_path = file.path('./dataSets',species, organ,paste0(data,'.rds'))

start = as.integer(args[2])
end = as.integer(args[3])
out_folder = file.path('./results/stability/',data)
ratio = ratio[start:end]
print(out_folder)

dir.create(out_folder)

data = readRDS(data_path)
if( class(data) == 'matrix'){
  counts = data
}else{
  counts = getCounts(data)
}

print(dim(counts))
counts = counts[which(rownames(counts) != ''),]
counts = counts[unique(rownames(counts)),]
print(dim(counts))
#eval = cal_all_tools(counts, 
#                     gSets_path,
#                     tools = c('AUCell','pagoda2',
#                               'Vision',
#                               'GSVA','ssGSEA',
#                               'seurat','zscore','plage'),
#                     filter = 'none',
#                     normalize = 'none',
#                     impute = 'none',
#                     highVarible = 'none',
#                     n_cores = n_cores)
#saveRDS(eval,file.path(out_folder,'allGenes.rds'))
#print('all genes success')

for(r in ratio){
  res = vector(mode='list',length=times)
  print(r)
  for(i in 1:times){
    if(dropType == 'genes'){
      curr = dropoutGenes_v1(counts,r)
      curr_gSets_path = gSets_path
    }else if(dropType == 'gSets'){
      curr = counts
      curr_gSets_path = './gmtFiles/stability/tmp.gmt'
      shuffleSets(gSets_path,r,
                  rownames(counts),
                  curr_gSets_path)
    }
    
    #n_genes = colSums(curr>0)
    #print(mean(n_genes))
    eval = vector(mode='list',length = length(tools))
    names(eval) = tools
    for(tool in tools){
      print(tool)
      x = c('none','sctransform')
      eval[[tool]] = cal_all_tools(curr, 
                                   curr_gSets_path,
                                   tools = tool,
                                   filter = x[1],
                                   normalize = x[2],
                                   impute = 'none',
                                   highVarible = 'none',
                                   n_cores = n_cores,
                                   ifMem = F)
    }
    
    #eval$n_genes = n_genes
    res[[i]] = eval
  }
  saveRDS(res,file.path(out_folder,paste0('ratio_V1_',r,'.rds')))
}

