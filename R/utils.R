getVarib = function(sc){
  tryCatch({
    sc = Seurat::FindVariableFeatures(sc,selection.method = 'vst',verbose=F)
    return('vst')
  },error=function(e){
    tryCatch({
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'disp',verbose=F)
      return('disp')
    },error=function(e){
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'mvp',verbose=F)
      return('mvp')
    })
  })
}


del_geneSets_roma = function(gSets){
  res = vector(mode='list')
  for( i in 1:length(gSets) ){
    res[[i]] = list(Name=names(gSets[i]),
                    Genes=gSets[[i]])
  }
  res
}


del_geneSets_vision = function(gSets){
  env = new.env(parent=globalenv())
  invisible(lapply(1:length(gSets),function(i) {
    genes = gSets[[i]]
    name = names(gSets[i])
    assign(name, genes, envir = env)
  }))
  return(env)
}

getGMT = function(gmt_file_path){
  paths = readLines(gmt_file_path)
  gsets = list()
  for(i in 1:length(paths)){
    t = strsplit(paths[i],'\t')[[1]]
    genes = t[3:length(t)]
    genes = genes[which(genes != "")]
    gsets[[t[1]]] = genes
  }
  return (gsets)
}

toGMT = function(gSets,
                 gmt_file_path){
  if(file.exists(gmt_file_path)){
    file.remove(gmt_file_path)
  }
  aa = gsub('-|,|/| |\\(|\\)','_',names(gSets))
  aa = gsub('_+','_',aa)
  aa = gsub('_\\b','',aa)
  aa = gsub(',','',aa)
  names(gSets) = aa
  for(i in 1:length(gSets)){
    genes = paste(gSets[[i]],collapse ='\t')
    aa = paste(names(gSets)[i],names(gSets)[i],genes,sep='\t')
    write.table(aa, file=gmt_file_path, row.names=F,col.names = F,quote=F,append = T)
  }
}

normToOne = function(score){
  center = sweep(score, 2, apply(score, 2, min),'-')
  R = apply(score, 2, max) - apply(score, 2, min)
  sweep(center, 2, R, "/")
}


getAllGenes = function(gSets_path){
  genes = c()
  gSets = getGMT(gSets_path)
  for(ss in gSets){
    genes = c(genes,ss)
  }
  unique(genes)
}




dropoutGenes_v1 = function(counts,
                           ratio=0.1){
  ## random drop expressed genes
  res = apply(counts,2,function(x){
    index = which(x > 0)
    x[sample(index,round(length(index)*ratio))] = 0
    x
  })
  res
}

dropoutGenes_v2 = function(counts,
                           ratio=0.1){
  res = apply(counts,2,function(x){
    rr = rank(x, ties.method = "first")
    drop_len = round(sum(x>0)*ratio)
    c = sort(rr[x>0])[1:drop_len]
    x[which(rr %in% c)] = 0
    x
  })
  res
}

dropGsets = function(gSets,
                     ratio=0.1){
  res = sapply(gSets, function(x){
    sample(x,length(x)*(1-ratio))
  })
  res
}

len_gSets = function(gSets){
  res = sapply(gSets,length)
  unlist(res)
}


getCounts = function(data){
  tryCatch({
    return(data@assays$data[[1]])
  },error=function(e){
    return(data@assays@data[[1]])
  })
}

dropoutSets = function(gSets_path,
                       ratio = 0.1,
                       out_path){
  gSets = getGMT(gSets_path)
  gSets2 = dropGsets(gSets, ratio)
  toGMT(gSets2, out_path)
}

shufflGsets = function(gSets,
                       ratio,
                       all_genes){
  res = sapply(gSets, function(x){
    c(sample(x,length(x)*(1-ratio)),
      sample(all_genes, length(x)*ratio))
  })
  res
}

shuffleSets = function(gSets_path,
                       ratio,
                       all_genes,
                       out_path){
  gSets = getGMT(gSets_path)
  gSets2 = shufflGsets(gSets, ratio, all_genes)
  toGMT(gSets2, out_path)
}
