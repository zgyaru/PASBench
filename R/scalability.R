library(reticulate)
use_python('/share/apps/anaconda3/bin/python3',required = T)
invisible(ulimit::memory_limit(102400))
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
suppressPackageStartupMessages(library(progeny))
suppressPackageStartupMessages(library("optparse"))

#library(scDblFinder)
# clustering
suppressPackageStartupMessages(library(SC3))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(HDF5Array))



data_folder = './dataSets/scalability/smart-seq2'
data_folder = './dataSets/scalability/10x'
gSets_folder = './gmtFiles/scalability/'
n_cores = 5

#dats = c('20_1','20_5','10_5','10_10','10_20','20_10')
#gsts = c('gSets_50','gSets_200','gSets_1000')
gsts = c('gSets_1000')

#dats = c('10_5','10_10','20_5','20_10','20_20')
dats = c('20_5','20_10','20_20')
#gsts = c('gSets_50','gSets_200')

all_runs = vector(mode='list',length=length(dats)*length(gsts))
i=1
res_name = c()
for(d in dats){
  for(g in gsts){
    all_runs[[i]] = c(d,g)
    res_name = c(res_name, paste(d,gsub('gSets_','',g),sep='_'))
    i=i+1
  }
}
names(all_runs) = res_name

runs_name = function(x){
  paste(x[1],gsub('gSets_','',x[2]),sep='_')
}

parseRuns = function(x,
                     tools,
                     data_folder,
                     gSets_folder){
  source('./utils.R')
  source('./tools.R')
  counts = readRDS(file.path(data_folder,paste0(x[1],'.rds')))
  counts = counts[!duplicated(rownames(counts)),]
  counts = counts[unique(rownames(counts)),]
  gSets_path = file.path(gSets_folder,paste0(x[2],'.gmt'))
  eval = cal_all_tools(counts,
                       gSets_path,
                       tools = tools,
                       filter = 'none',
                       normalize = 'none',
                       impute = 'none',
                       highVarible = 'none',
                       n_cores = 5,
                       ifStoreScore = F)
  #eval$cell_type = cell_type
  saveRDS(eval, file.path('./results/scalability/10x/',paste0(x[1],x[2],'_tools7.rds')))
  #saveRDS(eval, file.path('./results/scalability/smart-seq2/',paste0(x[1],x[2],'_tools7.rds')))
  #eval
}





#gsts = c('gSets_1000')
#dats = c('20_20')
#gsts = c('gSets_50','gSets_200')

#dats = c('10_5','20_5')
#gsts = c('gSets_50')
#eval_names = c()
#res = vector(mode='list')
#i=1
tools = c('AUCell','pagoda2','Vision','GSVA','ssGSEA','zscore','plage')
#tools = c('AUCell','pagoda2')
#tools = c('fscLVM')

#cl = makeCluster(n_cores)
#res_name = unlist(clusterApply(cl,all_runs,runs_name))
#print(all_runs)
#print(res_name)
#res = clusterApply(cl, all_runs, parseRuns, 
#                   tools,data_folder,
#                   gSets_folder)
#res_name = unlist(apply(all_runs, runs_name))
print(res_name)
#i=1
#for(r in all_runs){
#  print(r)
#  res[[i]] = parseRuns(r,tools,
#                       data_folder,
#                       gSets_folder)
#}
res = lapply(all_runs, parseRuns, 
            tools=tools, 
            data_folder=data_folder, 
            gSets_folder=gSets_folder)
#print(names(res))
#names(res) = res_name
#saveRDS(res, './results/scalability/smart-seq2/20200627_tools7.rds')
#saveRDS(res, './results/scalability/10x/20200627_tools7.rds')




# for(data_name in dats){
#   counts = readRDS(file.path(data_folder,paste0(data_name,'.rds')))
#   counts = counts[!duplicated(rownames(counts)),]
#   print(dim(counts))
#   for(gSets_name in gsts){
#     gSets_path = file.path(gSets_folder,paste0(gSets_name,'.gmt'))
#     eval = cal_all_tools(counts, 
#                          gSets_path, 
#                          tools = tools,
#                          filter = 'none',
#                          normalize = 'none',
#                          impute = 'none',
#                          highVarible = 'none',
#                          n_cores = n_cores,
#                          ifStoreScore = F)
#     gg = gsub('gSets_', '', gSets_name)
#     eval_names = c(eval_names, paste(data_name, gg, sep='_'))
#     res[[i]] = eval
#     i=i+1
#   }
# }
# names(res) = eval_names
# #saveRDS(res,'/share/pub/zhangyr/single-RNA/R_proj/comparison/results/scalability/res_smartSeq_20200622_v1.rds')
# saveRDS(res,'/share/pub/zhangyr/single-RNA/R_proj/comparison/results/scalability/res_10x_20200622_v2.rds')
# 

