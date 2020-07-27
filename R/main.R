
#' calculation Pathway Activity Score
#'
#' parameters `counts`,  `gSets_path` and `gSets` are consistent across all tools
#'
#' other parameters see details of tools
#'
#'
#'
#' @param counts gene expresion matrix, rows are genes and cols are cells
#' @param tool tool name
#' @param gmt_file pathways in GMT format
#' @param species species; human or mouse
#' @param pathway abbreviation for pathway database
#' @param filter whether filtering for genes expressed in less than 5 percent cells
#' @param normalize normalization method
#'
#'
#' @examples
#' counts = load_counts()
#' score = calculate_PAS(counts, 'AUCell', 'human','kegg')

#' @export
calculate_PAS = function(counts,
                         tool,
                         species = 'none',
                         pathway = 'none',
                         gmt_file = 'none',
                         filter = F,
                         normalize = 'sctransform',
                         n_cores = 3,
                         rand_seed = 123){

  if(gmt_file == 'none'){
    if(species == 'none' || pathway == 'none'){
      stop("please define species and pathway when do not assign gmt_file")
    }

    gmt_folder = system.file("gmtFiles", package = "PASBench")
    if(gmt_folder == ""){
      stop("could not find shiny directory, try-re-installing 'PASBench'.")
    }

    gSets_path = file.path(gmt_folder,species,paste0(pathway,'.gmt'))
  }else{
    if(species != 'none' || pathway != 'none'){
      warning(" 'gmt_file' is already present.
              Ignoring 'species' and 'pathway'.")
    }
    gSets_path = gmt_file
  }

  if(tool %in% c('pagoda2','vision')){
    warning(" 'log' transform will counter error for 'pagoda2' or 'vision'.
            force the 'normalize' to 'sctransform' or you could modify your code
            to call other normalization function.")
  }
  counts = na.omit(counts)
  counts = counts[which(rownames(counts) != ''),]


  score = cal_all_tools(counts,
                        gSets_path,
                        tools = tool,
                        filter = filter,
                        normalize = normalize,
                        n_cores = n_cores,
                        rand_seed = rand_seed
                        )

  score[[tool]]
}







#' calculation seven tools
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
#' @param filter whether filtering for genes expressed in less than 5 percent cells
#' @param normalize normalization method

cal_all_tools = function(counts,
                         gSets_path,
                         #cells_label,
                         tools = c('AUCell','pagoda2',
                                   'fscLVM','Vision',
                                   'ROMA','GSVA','ssGSEA',
                                   'plage','zscore'),
                         filter = F,
                         normalize = c('log','CLR','RC','scran','sctransform','none'),
                         n_cores = 3,
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
  }else{
    cat("Do not filter genes.\n")
  }

  ## normalization
  tryCatch({
    if(normalize == 'log'){
      counts = Seurat::NormalizeData(counts,normalization.method = 'LogNormalize',verbose=0)
      counts = Seurat::ScaleData(counts)
      cat("log nomalization success\n")
    }else if(normalize == 'CLR'){
      counts = Seurat::NormalizeData(counts,normalization.method = 'CLR',verbose=0)
      counts = Seurat::ScaleData(counts)
      cat("CLR nomalization success\n")
    }else if(normalize == 'RC'){
      counts = Seurat::NormalizeData(counts,normalization.method = 'RC',verbose=0)
      counts = Seurat::ScaleData(counts)
      cat("RC nomalization success\n")
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
      cat("scran nomalization success\n")
    }else if(normalize == 'sctransform'){

      expr_ob = Seurat::CreateSeuratObject(counts=counts[,])
      expr_ob = Seurat::SCTransform(expr_ob,verbose=FALSE)
      counts = as.matrix(expr_ob@assays$SCT@data)
      rm(expr_ob)
      cat("scTransform nomalization success\n")
    }else if(normalize == 'scnorm_9'){
      tryCatch({
        counts = na.omit(counts)
        DataNorm = SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K=9,NCores=n_cores)
        counts = DataNorm@assays$data$normcounts
        rm(DataNorm)
        cat("SCnorm nomalization success\n")
      },error=function(e){
        print("SCnorm error")
        print(e)
        return("error")
      })
    }else if(normalize == 'scnorm_5'){
      tryCatch({
        counts = na.omit(counts)
        DataNorm = SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K=5,NCores=n_cores)
        counts = DataNorm@assays$data$normcounts
        rm(DataNorm)
        cat("SCnorm nomalization success\n")
      },error=function(e){
        print("SCnorm error")
        print(e)
        return("error")
      })

    }else{
      cat("Do not normalize PAS.\n")
    }
    },error=function(e){
      print("normalize error")
      return("error")
    })




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


    eval_tools[[i]] = score

  }
  names(eval_tools) = tools
  eval_tools
}


#' preparation of visulization object
#'
#' @param pas_score a matrix of PAS
#'
#'
#' @export
prepare_vis = function(pas_score,
                       n_pcs = 10,
                       resolution = 0.5){
  if(is.character(pas_score) == 'character'){
    stop(" 'pas_score' is a character, cannot convert a character to seurat object ")
  }
  perplexity = min(floor((ncol(pas_score)-1)/3),30)
  n_pcs = min(n_pcs,nrow(pas_score)-1)

  sc = Seurat::CreateSeuratObject(pas_score)
  sc = Seurat::ScaleData(sc)
  sc = Seurat::FindVariableFeatures(sc,selection.method = getVarib(sc),verbose=F)
  sc = Seurat::RunPCA(sc,verbose=F,npcs = n_pcs)
  sc = Seurat::RunUMAP(sc,dims = 1:n_pcs, n.components = 2,
                       verbose = F)
  sc = Seurat::RunTSNE(sc,dims = 1:n_pcs, perplexity=perplexity, verbose = F)
  sc = Seurat::FindNeighbors(sc,dims = 1:n_pcs)
  sc = Seurat::FindClusters(sc,resolution=resolution)
  return(sc)
}




#' visualization of PAS
#'
#'
#'
#'
#' @param seurat_oj a seurat object
#'

#' @export
PAS_vis = function(seurat_oj){
  suppressPackageStartupMessages("")
  appDir = system.file("shiny", package = "PASBench")
  if(appDir == ""){
    stop("could not find shiny directory, try-re-installing 'PASBench'.")
  }

  .GlobalEnv$.sc_oj = seurat_oj
  .GlobalEnv$.pathways = rownames(Seurat::GetAssayData(seurat_oj))


  on.exit(rm(list = c(".sc_oj",".pathways")))

  shiny::runApp(appDir, launch.browser = T, display.mode = "normal")
}
