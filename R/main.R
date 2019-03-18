
#' initDb - Initialization of the package so the GTEx networks can be used with
#' CoExpNets
#'
#' @param mandatory If this parameter is `TRUE` then the networks will be added no matter whether they were already there.
#'
#' @return No value
#' @export
#'
#' @examples
initDb = function(mandatory=F){
  the.dir = system.file("", "gtexv7", package = "CoExpGTExV7")
  nets = CoExpGTExV7::getGTExTissues()
  for(net in nets){
    CoExpNets::addNet(which.one="CoExpGTExV7",
           tissue=net,
           netfile=CoExpGTExV7::getGTExNet(net),
           ctfile=paste0(the.dir,"/",net,"_celltype.tsv"),
           gofile=paste0(the.dir,"/",net,"_gprof.tsv"),
           exprdatafile=paste0(the.dir,"/",net,".rds"),
           overwrite=mandatory)
  }
}

#' gteGTExNet - Accessing a network object directly
#'
#' @param tissue One of the labels that can be obtained by calling getGTExTissues() method
#' to refer to a specific network within the package
#'
#' @return The RDS full file name to the network
#' @export
#'
#' @examples
getGTExNet = function(tissue){

  the.dir = system.file("", "gtexv7", package = "CoExpGTExV7")
  files = list.files(the.dir)

  net.file = files[grep(paste0("net",tissue,"\\.\\d+\\.it\\.50\\.rds$"),files)]
  if(length(net.file) == 0)
    return(NULL)

  return(paste0(the.dir,"/",net.file))
}

#' Title Getting all GTEx available tissues
#'
#' @return
#' @export
#'
#' @examples
getGTExTissues = function(){

  the.dir = system.file("", "gtexv7", package = "CoExpGTExV7")
  files = list.files(the.dir)

  net.files = files[grep("net\\w+\\.\\d+\\.it\\.50\\.rds$",files)]
  net.files = gsub("net","",net.files)
  net.files = gsub(".\\d+\\.it\\.50\\.rds$","",net.files)

  return(net.files)
}

repareGeneNames = function(){
  tissues = CoExpGTExV7::getGTExTissues()
  for(tissue in tissues){
    file = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one="CoExpGTExV7",only.file = T)
    cat("Reading network for",tissue,"from", file,"\n")
    net = readRDS(file)
    oldnames = names(net$moduleColors)
    names(net$moduleColors) = unlist(lapply(stringr::str_split(names(net$moduleColors),"\\."),function(x){return(x[[1]])}))
    if(!identical(oldnames,names(net$moduleColors))){
      saveRDS(net,file)
      cat("Writting down network",tissue)
    }
      
  }
}


testNetworks = function(){
  tissues = CoExpGTExV7::getGTExTissues()
  for(tissue in tissues){
    file = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one="CoExpGTExV7",only.file = T)
    stopifnot(file.exists(file))
    cat("Network for",tissue,"has",length(readRDS(file)$moduleColors),"genes\n")
    file = CoExpNets::getExprDataFromTissue(tissue=tissue,which.one="CoExpGTExV7",only.file = T)
    stopifnot(file.exists(file))
    cat("Network for",tissue,"has",nrow(readRDS(file)),
        "samples\n")
    go =  CoExpNets::findGO(tissue=tissue,which.one="CoExpGTExV7")
    stopifnot(file.exists(go))
    cat("Network for",tissue,"has",nrow(read.delim(go)),
        "annotation signals\n")
    ct =  CoExpNets::findCT(tissue=tissue,which.one="CoExpGTExV7")
    stopifnot(file.exists(ct))
    cat("Network for",tissue,"has",nrow(read.delim(ct)),
        "cell type annotation signals\n")
  }
  cat("Test OK!!\n")
}

