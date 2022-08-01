#' DEG analaysis using Limma
#'
#' This function performs DEG analysis using
#'
#' @param motherMatrix row name should be probe or gene and col name should be sample names.
#' @param case pattern of case samples
#' @param control pattern of control samples
#' @return res DEG results.
#' @examples
#' mothermatrix colnames: treat_A, treat_B, treat_C, no_treat_A, no_treat_B, no_treat_C
#' LimmaDEG(mat, "treat", "no_treat")
#' @export
LimmaDEG <- function(motherMatrix  , case, control){
  x = motherMatrix
  y = case
  z = control
  HT <- x %>% as.data.frame() %>%  dplyr::select(contains(y)|contains(z))
  rownames(HT) = rownames(motherMatrix)

  cse <- c()
  for(k in 1:length(grep(y, colnames(x)))){
    cse[k] <- y
  }
  ctrl <- c()
  for(k in 1:length(grep(z, colnames(x)))){
    ctrl[k] <- z
  }

  HT_group <- c(paste("A",cse,sep = "_"), paste("B",ctrl,sep = "_")) %>% as.factor()
  design <- model.matrix(~-1+ HT_group)

  HT_group <- factor(x= c(cse, ctrl),levels = c(y,z))
  colnames(design) <- levels(HT_group)

  fit <- lmFit(HT, design)
  constrast.matrix<-makeContrasts(contrasts = sprintf("%s - %s", y, z), levels=HT_group)
  fit.cont<-contrasts.fit(fit,constrast.matrix)
  efit<- eBayes(fit.cont)
  result<-topTable(efit, n=nrow(efit))
  DEGs <- result
  DEGs[, 'id'] <- rownames(DEGs)
  return(DEGs)
}
