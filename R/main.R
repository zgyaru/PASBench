#' calculation pipeline
#'
#' parameters `counts`,  `gSets_path` and `gSets` are consistent across all tools
#'
#' other parameters see details of tools
#'
#'
#'
#' @param counts gene expresion matrix, rows are genes and cols are cells
#' @param gSets_path patwhays/gene sets in clasical GMT format
#' @param tool select PAS tools
#' @param filter whether filtering for genes expressed in less than 5% cells
#' @param normalize normalization method

#' @export
cal_all_tools = function(counts,
                         gSets_path,
                         #cells_label,
                         tools = c('AUCell','pagoda2',
                                   'fscLVM','Vision',
                                   'ROMA','GSVA','ssGSEA',
                                   'plage','zscore'),
                         filter = F,
                         normalize = c('log','CLR','RC','scran','sctransform','none'),
                         n_cores = 4,
                         rand_seed = 123){


  if(round(counts[which(counts>0)[1],1]) == counts[which(counts>0)[1],1]){
    mat_type ='counts'
  }else{
    mat_type = 'tpm'
  }

  ## filtering
  if(filter){
    num_cells = Matrix::rowSums(x = counts > 0)
    counts = counts[which(x = num_cells >= ncol(counts)*0.05),]
    cat(paste("After filtering, dimensions of count are:",dim(counts),collapse=' '))
    cat('\n')
  }

  ## normalization
  tryCatch({
    if(normalize == 'log'){
      counts = Seurat::NormalizeData(counts,normalization.method = 'LogNormalize',verbose=0)
      counts = Seurat::ScaleData(counts)
    }else if(normalize == 'CLR'){
      counts = Seurat::NormalizeData(counts,normalization.method = 'CLR',verbose=0)
      counts = Seurat::ScaleData(counts)
    }else if(normalize == 'RC'){
      counts = Seurat::NormalizeData(counts,normalization.method = 'RC',verbose=0)
      counts = Seurat::ScaleData(counts)
    }else if(normalize == 'scran'){
      genes = rownames(counts)
      cells = colnames(counts)
      if(mat_type == 'counts'){
        counts = apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
        sc = SingleCellExperiment::SingleCellExperiment(
          assays = list(counts = counts)
        )
        if(ncol(counts)>300){
          clusters = scran::quickCluster(sc, assay.type = "counts")
          sc = scran::computeSumFactors(sc, clusters=clusters, assay.type = "counts")
        }else{
          sc = scran::computeSumFactors(sc,assay.type = "counts")
        }

        sc = scater::normalize(sc,exprs_values  ='counts')
        counts = sc@assays$data$logcounts
        rownames(counts) = genes
        colnames(counts) = cells
        rm(genes,cells)
      }else{
        sc = SingleCellExperiment::SingleCellExperiment(
          assays = list(tpm = counts)
        )
        if(ncol(counts)>300){
          clusters = scran::quickCluster(sc, assay.type = "tpm")
          sc = scran::computeSumFactors(sc, clusters=clusters, assay.type = "tpm")
        }else{
          sc = scran::computeSumFactors(sc,assay.type = "tpm")
        }

        sc = scater::normalize(sc,exprs_values  ='tpm')
        counts = sc@assays$data$logcounts
        rownames(counts) = genes
        colnames(counts) = cells
        rm(genes,cells)
      }


    }else if(normalize == 'sctransform'){
      counts = counts[which(rownames(counts) != ''),]
      expr_ob = Seurat::CreateSeuratObject(counts=counts[,])
      expr_ob = Seurat::SCTransform(expr_ob,verbose=FALSE)
      counts = as.matrix(expr_ob@assays$SCT@data)
      rm(expr_ob)
    }else if(normalize == 'scnorm_9'){
      tryCatch({
        counts = na.omit(counts)
        DataNorm = SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K=9,NCores=n_cores)
        counts = DataNorm@assays$data$normcounts
        rm(DataNorm)
      },error=function(e){
        print("scnorm error")
        print(e)
        return("error")
      })
    }else if(normalize == 'scnorm_5'){
      tryCatch({
        counts = na.omit(counts)
        DataNorm = SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K=5,NCores=n_cores)
        counts = DataNorm@assays$data$normcounts
        rm(DataNorm)
      },error=function(e){
        print("scnorm error")
        print(e)
        return("error")
      })

    }},error=function(e){
      print("normalize error")
      return("error")
    })
  cat("normalize success\n")



  eval_tools = vector(mode="list")

  n_cores = n_cores
  gSets = getGMT(gSets_path)

  counts = as.matrix(counts)
  for(i in 1:length(tools)){
    tool = tools[i]
    t_start = Sys.time()
    gc()
    score = switch(tool,
                   AUCell = cal_AUCell(counts[,],
                                       gSets,
                                       n_cores),  #0,0.8
                   pagoda2 = cal_pagoda2(counts[,],
                                         gSets,
                                         n_cores=n_cores),  #-13,65
                   fscLVM = cal_fscLVM(counts[,],
                                       gSets_path,
                                       type = mat_type),
                   Vision = cal_vision(counts[,],
                                       gSets_path,
                                       n_cores),  #-0.6895973, 3.32
                   ROMA = cal_ROMA(counts[,],
                                   gSets,
                                   n_cores),   #-0.07, 0.4
                   GSVA = cal_gsva(counts[,],
                                   gSets,
                                   n_cores), #-0.98,0.94
                   ssGSEA = cal_ssgsea(counts[,],
                                       gSets,
                                       n_cores),#-0.5,0.5
                   plage = cal_plage(counts[,],
                                     gSets,
                                     n_cores),  #-1,1
                   zscore = cal_zscore(counts[,],
                                       gSets,
                                       n_cores),   #-4,29
                   seurat = cal_SeuratScore(counts[,],
                                            gSets),  #-1753,6958
                   scSigScore = cal_scSigScore(),
                   gene = counts)


    if(!ifStoreScore){
      if(class(score) == 'matrix'){
        score='NA'
      }
    }
    eval_tools[[i]] = score

  }
  names(eval_tools) = tools
  eval_tools
}
