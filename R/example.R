
#' loading counts matrix

#' @export
load_counts = function(){

  folder = system.file("example", package = "PASBench")
  if(folder == ""){
    stop("could not find example directory, try-re-installing 'PASBench'.")
  }
  readRDS(file.path(folder, 'counts.rds'))
}

#' loading visualization object

#' @export
load_visOJ = function(){

  folder = system.file("example", package = "PASBench")
  if(folder == ""){
    stop("could not find example directory, try-re-installing 'PASBench'.")
  }
  readRDS(file.path(folder, 'vis_oj.rds'))
}
