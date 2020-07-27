library(reticulate)
use_python('/share/apps/anaconda3/bin/python3',required = T)
library(igraph)
library(scran)
library(Seurat)

#library(ClusterR)
library(dplyr)
##### evaluation metrics

## from https://github.com/mlampros/ClusterR/blob/master/R/clustering_functions.R
cluster_metrics = function(true_labels, clusters){
  
  if (!is.vector(true_labels)) stop('true_labels should be a numeric vector')
  if (!is.vector(clusters)) stop('clusters should be a numeric vector')
  if (length(true_labels) != length(clusters)) stop('the length of the true_labels vector should equal the length of the clusters vector')
  
  true_labels = as.numeric(true_labels)
  clusters = as.numeric(clusters)
  
  tbl = table(clusters, true_labels)
  conv_df = as.data.frame.matrix(tbl)
  
  
  tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))
  tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))
  tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))
  fp = tp_plus_fp - tp
  fn = tp_plus_fn - tp
  tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn
  
  ## ARI
  prod_comb = (tp_plus_fp * tp_plus_fn) / gmp::asNumeric(gmp::chooseZ(length(true_labels), 2))
  mean_comb = (tp_plus_fp + tp_plus_fn) / 2.0
  ARI = round((tp - prod_comb) / (mean_comb - prod_comb), 4)
  
  ## purity
  tmp_pur = apply(conv_df, 1, max)
  res_purity = sum(tmp_pur)/length(true_labels)
  res_purity = round(res_purity, 4)
  
  ## F1
  prec = tp / (tp + fp)
  rec = tp / (tp + fn)
  F1 = round(2.0 * ((prec * rec) / (prec + rec)), 4)
  
  ## NMI
  NMI = 0.0
  unq_true = unique(true_labels)
  unq_clust = unique(clusters)
  if (length(unq_true) == 1 && length(unq_clust) == 1) {  
    NMI = 1.0
  }else{
    mutual_information = 0.0
    for (i in 1:nrow(conv_df)) {
      for (j in 1:ncol(conv_df)) {
        if (conv_df[i,j] > 0.0) {
          mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log2(as.numeric(gmp::as.bigz(as.numeric(sum(tbl)) * as.numeric(conv_df[i,j])) / gmp::as.bigz(as.numeric(sum(conv_df[i,])) * as.numeric(sum(conv_df[,j]))))))
        }
      }
    }
    entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))
    entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))
    NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))
  }
  
  res = c(ARI, NMI, F1, res_purity)
  names(res) = c("ARI","NMI", "F1","purity")
  res
}

cluster_GMM = function(data,
                       n_class){
  #sd_data = apply(data,1,sd)
  #data = data[which(sd_data>quantile(sd_data)[4]),]
  #print(dim(data))
  data = center_scale(t(data),mean_center=T,sd_scale=T)
  gmm = GMM(data,n_class,
            dist_mode = 'maha_dist',
            seed_mode = 'random_subset',
            km_iter = 10,
            em_iter = 10,
            verbose = F)
  pr = predict_GMM(data, gmm$centroids,gmm$covariance_matrices,gmm$weights)
  pr$cluster_labels
}

data_ARI = function(reduc_dat, 
                    true_labels){
  n_class = length(unique(true_labels))
  p = kmeans(dist(reduc_dat),n_class)
  cluster_m = cluster_metrics(true_labels,p$cluster)
  cluster_m['ARI']
}

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

clustering_louvain=function(data,
                            true_labels,
                            n_class){
  graph = scran::buildSNNGraph(data, transposed=T,k=5,d=NA)
  res = igraph::cluster_louvain(graph)$membership
  if(max(res) == n_class){
    clu = res
  }else{
    n_class = min(n_class, max(res))
    print(n_class)
    if (max(res) <= 3){
      hclu <- hclust(dist(data))
      clu <- cutree(hclu,n_class)
    } else {
      
      cc <- aggregate(data,list(res),mean)
      cc <- as.matrix(cc[,-1])
      hclu <- hclust(dist(cc))
      clu <- cutree(hclu,n_class)
      clu <- clu[res]
    }
  }
  mclust::adjustedRandIndex(clu, true_labels)
}

clustering_kmeans = function(data,
                             true_labels,
                             n_class){
  p = kmeans(dist(data),n_class)
  mclust::adjustedRandIndex(p$cluster, true_labels)
}


dimension_reduc_metric = function(data,
                                  true_labels,
                                  dimension_red = 'umap',
                                  dimension_type = 'seurat',
                                  selection.method= 'vst',
                                  normalization=F){
  
  true_labels = as.numeric(true_labels)
  perplexity = min(floor((ncol(data)-1)/3),20)
  n_pca = min(20,nrow(data))
  #dims = min(20,n_pca)
  #print(perplexity)
  if(dimension_type == 'seurat'){
    sc = CreateSeuratObject(data)
    #sc = Seurat::NormalizeData(sc,verbose=F)
    if(normalization){
      sc = Seurat::SCTransform(sc,verbose=FALSE)
    }else{
      sc = ScaleData(sc,verbose = F)
    }
    sc = FindVariableFeatures(sc,selection.method = getVarib(sc),verbose=F)
    #sc = Seurat::FindVariableFeatures(sc,selection.method=selection.method,verbose=F)
    sc = RunPCA(sc,verbose=F,npcs = n_pca)
    pca = as.data.frame.matrix(Embeddings(sc[['pca']]))
    if(dimension_red == 'tsne'){
      #max_ari = 0
      #max_reduc_dat = NA
      #for(n_com in c(2,3)){
      sc = RunTSNE(sc,perplexity = perplexity,dim.embed = 2,verbose=F,check_duplicates = FALSE)
      reduc_dat = as.data.frame.matrix(Embeddings(sc[['tsne']]))
      #curr_ari = data_ARI(reduc_dat, true_labels)
      #if(curr_ari > max_ari){
      #  max_ari = curr_ari
      #  max_reduc_dat = reduc_dat
      #  #print(c(dims,n_com))
      #}
      #}
      #reduc_dat = max_reduc_dat
      
    }else if(dimension_red == 'umap'){
      #max_ari = 0
      #max_reduc_dat = NA
      dims = 10
      n_com = 2
      #for(dims in c(5,10,15,20)){
      #  for(n_com in c(3,5,7,9)){
      #    tryCatch({
      #     if(dims>n_com){
      #print(c(dims,n_com))
      sc = RunUMAP(sc,dims = 1:dims, n.components = n_com,
                   verbose = F,umap.method='umap-learn')
      reduc_dat = as.data.frame.matrix(Embeddings(sc[['umap']]))
      #        curr_ari = data_ARI(reduc_dat, true_labels)
      #        if(curr_ari > max_ari){
      #          max_ari = curr_ari
      #          max_reduc_dat = reduc_dat
      #        }  
      #      }
      #    },error=function(e){
      #      
      #    })
      #  }
      #}
      #reduc_dat = max_reduc_dat
    }
  }else if(dimension_type == 'original'){
    data = center_scale(t(data),mean_center=T,sd_scale=T)
    pca_dat = stats::princomp(data)
    pca = pca_dat$score[,1:n_pca]
    
    if(nrow(data)>ncol(data)){
      #pca_dat = stats::princomp(data)
      if(dimension_red == 'tsne'){
        reduc_dat = Rtsne::Rtsne(pca,pca=F,dims = 2,perplexity = perplexity)$Y
      }else if(dimension_red == 'umap'){
        reduc_dat = umap::umap(pca_dat$scores[,1:n_pca])$data[,1:2]
      }
    }else{
      if(dimension_red == 'tsne'){
        perplexity = floor((nrow(data)-1)/3)
        reduc_dat = Rtsne::Rtsne(pca,pca=T,dims = 3,perplexity = perplexity)$Y
      }else if(dimension_red == 'umap'){
        reduc_dat = umap::umap(data)$data[,1:3]
      }
    }
  }
  
  
  
  #pca_dat = stats::princomp(data)
  #tsne_dat = Rtsne::Rtsne(data,pca=F,dims = 3,initial_dims=n_pca,
  #                        pca_center=T,pca_scale=T)$Y
  
  #umap_output = umap::umap(pca_dat$scores[,1:30], dims = 3)
  sil1 = cluster::silhouette(true_labels, dist(reduc_dat))
  #sil2 = cluster::silhouette(true_labels, dist(umap_output$layout))
  sil = group_by(as.data.frame.matrix(sil1),cluster)
  sil = as.data.frame(summarise(sil, mean_sil=mean(sil_width)))
  sil = colMeans(sil)[2]
  #sil=as.data.frame.matrix(sil1)
  #  sil=mean(sil$sil_width,na.rm=T)
  n_class = length(unique(true_labels))
  print(n_class)
  #p = kmeans(dist(reduc_dat),n_class)
  #cluster_m = cluster_metrics(true_labels,p$cluster)
  
  ari_lov = clustering_louvain(pca,
                               true_labels,
                               n_class)
  ari_km = clustering_kmeans(pca,
                             true_labels,
                             n_class)
  
  
  res = c(sil, ari_lov, ari_km)
  names(res) = c('mean_sil','ARI_LV','ARI_KM')
  res
}



