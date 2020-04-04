#' diann: Report processing and protein quantification for MS-based proteomics.
#' @description A set of functions for dealing with mass spectrometry-based proteomics analysis reports.
#' @section diann functions:
#' diann_load
#' diann_matrix
#' diann_maxlfq
#' diann_save
#'
#' @docType package
#' @name diann

library(data.table)

cast <- function(df, sample.header, id.header, quantity.header) {
  x <- melt.data.table(df, id.vars = c(sample.header, id.header), measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- dcast.data.table(x, as.formula(paste0(id.header,'~',sample.header)), value.var = "value") 
  piv[[1]] <- NULL
  as.matrix(piv)
}

cast_aggregate <- function(df, sample.header, id.header, quantity.header) {
  x <- melt.data.table(df, id.vars = c(sample.header, id.header), measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- dcast.data.table(x, as.formula(paste0(id.header,'~',sample.header)), value.var = "value", fun.aggregate = function(x) max(x, na.rm=TRUE)) 
  piv[[1]] <- NULL
  piv = as.matrix(piv)
  piv[is.infinite(piv)] <- NA
  piv
}

pivot <- function(df, sample.header, id.header, quantity.header) {
  x <- melt.data.table(df, id.vars = c(sample.header, id.header), measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- as.data.frame(dcast.data.table(x, as.formula(paste0(id.header,'~',sample.header)), value.var = "value")) 
  rownames(piv) <- piv[[1]]
  piv[[1]] <- NULL
  piv <- piv[order(rownames(piv)),]
  as.matrix(piv)
}

pivot_aggregate <- function(df, sample.header, id.header, quantity.header) {
  x <- melt.data.table(df, id.vars = c(sample.header, id.header), measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- as.data.frame(dcast.data.table(x, as.formula(paste0(id.header,'~',sample.header)), value.var = "value", fun.aggregate = function(x) max(x, na.rm=TRUE))) 
  rownames(piv) <- piv[[1]]
  piv[[1]] <- NULL
  piv <- piv[order(rownames(piv)),]
  piv = as.matrix(piv)
  piv[is.infinite(piv)] <- NA
  piv
}

#' Read a table as a data frame
#'
#' @param file file name
#' @export
diann_load <- function(file) {
  as.data.frame(fread(file, stringsAsFactors = FALSE))
}

#' Generate a matrix with quantities. Zero quantities will be replaced with NA values.
#'
#' @param x DIA-NN report
#' @param id.header protein/peptide/precursor Id column name
#' @param quantity.header quantity column name
#' @param proteotypic.only only proteotypic peptides and the respective proteins should be considered?
#' @param q precursor q-value threshold
#' @param protein.q uniquely identified protein q-value threshold
#' @param pg.q protein group q-value threshold
#' @param gg.q gene group q-value threshold
#' @export
diann_matrix <- function(x, id.header = "Precursor.Id", quantity.header = "Precursor.Normalised", proteotypic.only = F, q = 0.01, protein.q = 1.0, pg.q = 1.0, gg.q = 1.0) {
  
  df <- as.data.table(x)
  if (proteotypic.only) df <- df[which(df[['Proteotypic']] != 0),]
  df <- unique(df[which(df[[id.header]] != "" & df[[quantity.header]] > 0 & df[['Q.Value']] <= q & df[['Protein.Q.Value']] <= protein.q & df[['PG.Q.Value']] <= pg.q & df[['GG.Q.Value']] <= gg.q),c("File.Name", id.header, quantity.header),with=FALSE])
  
  is_duplicated = any(duplicated(paste0(df[["File.Name"]],":",df[[id.header]])))
  if (is_duplicated) {
    warning("Multiple quantities per id: the maximum of these will be calculated")
    pivot_aggregate(df,"File.Name",id.header,quantity.header)
  } else {
    pivot(df,"File.Name",id.header,quantity.header)
  }
}

#' Quantify using the MaxLFQ algorithm https://doi.org/10.1074/mcp.M113.031591.
#' This function can be used to calculate protein quantities from peptide/precursor/fragment quantities or, for example, to calculate peptide quantities from precursor/fragment quantities.
#'
#' @param x data frame or data table
#' @param sample.header sample Id column name
#' @param group.header column name corresponding to the group Id, e.g. protein Id
#' @param id.header precursor/peptide Id column name
#' @param quantity.header precursor/peptide quantity column name
#' @param margin quantities below exp(margin) might be treated as NA
#' @export
#' @examples
#' df <- data.frame(File.Name = c("A","A","A","B","B","B"),Protein.Names=rep("ALB",6),
#'       Precursor.Id=rep(c("PEPTIDE","EPTIDEP","PTIDEPE"),2),Precursor.Normalised=c(20,10,5,25,12,NA))
#' diann_maxlfq(df)
diann_maxlfq <- function(x, sample.header = "File.Name", group.header = "Protein.Names", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised", margin = -10.0) {
  
  df <- as.data.table(x)
  df <- unique(df[which(df[[group.header]] != ""),c(sample.header, group.header, id.header, quantity.header),with=FALSE])
  
  df[[sample.header]] <- as.character(df[[sample.header]])
  df[[group.header]] <- as.character(df[[group.header]])
  df[[id.header]] <- as.character(df[[id.header]])
  df[[quantity.header]] <- as.numeric(df[[quantity.header]])
  if (any(df[[quantity.header]] < 0, na.rm=T)) stop("Only non-negative quantities accepted")
  
  is_duplicated = any(duplicated(paste0(df[[sample.header]],":",df[[group.header]],":",df[[id.header]])))
  if (is_duplicated) warning("Multiple quantities per id: the maximum of these will be calculated")
  
  if (margin > -2.0) {
    margin <- -2.0
    warning("margin reset to -2.0")
  }
  if (margin < -20.0) {
    margin <- -20.0
    warning("margin reset to -20.0")
  }
  
  df[[quantity.header]][which(df[[quantity.header]] == 0)] <- NA 
  df[[quantity.header]] <- log(df[[quantity.header]])
  df[[quantity.header]][which(df[[quantity.header]] <= margin)] <- NA 
  df <- df[!is.na(df[[quantity.header]]),]
  
  proteins <- unique(df[[group.header]])
  m <- length(proteins)
  samples <- unique(df[[sample.header]])
  n <- length(samples)
  
  result <- matrix(NA, nrow = m, ncol = n)
  rownames(result) = proteins
  colnames(result) = samples
  
  for (i in 1:length(proteins)) {
    if (is_duplicated) {
      piv <- cast_aggregate(df[which(df[[group.header]] == proteins[i]),], sample.header, id.header, quantity.header) 
    } else {
      piv <- cast(df[which(df[[group.header]] == proteins[i]),], sample.header, id.header, quantity.header) 
    }
    if (nrow(piv)==1 | ncol(piv) == 1) {
      res <- piv[1,]
    } else {
      piv[is.na(piv)] <- -1000000.0
      ref = col_max(as.vector(piv), nrow(piv), ncol(piv))
      columns <- which(ref > margin)
      identified <- piv[,columns]
      if (ncol(identified) >= 2) {
        res <- maxlfq_solve(as.vector(identified), nrow(identified), ncol(identified), margin * 1.001)
        ref[columns] <- res
      } else res <- ref
      res[which(res <= margin)] <- NA
    }
    result[i,match(colnames(piv), samples)] <- res
  }
  exp(result)
}

#' Save a data frame to a tab-separated file
#'
#' @param x data frame or data table
#' @param file file name
#' @param row.names save row names?
#' @export
diann_save <- function(x, file, row.names=F) {
  fwrite(x, file, row.names=row.names, quote=F, sep='\t')
}

.onUnload <- function(libpath) {
  library.dynam.unload("diann", libpath)
}