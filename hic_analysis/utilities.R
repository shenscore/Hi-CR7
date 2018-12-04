library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
###
sourceCpp('/YOUR/PATH/TO/straw/R/straw-R.cpp')
###
library(Matrix)
library(wvtool)
library(magrittr)
library(methods) # avoid script error



refine_chrlen <- function(x){
  if(x > 1000000){
    x <- round(x/1000000,1)
    paste(x,'MB',sep = '')
  }else{
    x <- round(x/1000,0)
    paste(x,'KB',sep = '')
  }
}

refine_pair_chrlen <- function(a,b){
  x <- max(a,b)
  if(x > 1000000){
    a <- round(a/1000000,1)
    b <- round(b/1000000,1)
    a <- paste(a,'MB',sep = '')
    b <- paste(b,'MB',sep = '')
  }else{
    a <- round(a/1000,0)
    b <- round(b/1000,0)
    a <- paste(a,'KB',sep = '')
    b <- paste(b,'KB',sep = '')
  }
  return(c(a,b))
}

get_comp_hic_mat <- function(hic1,hic2,res,norm = 'KR',chr=NULL,region=NULL){
  #TODO scale normalization
  mat1  <- get_mat_from_hic(hic1,norm,chr,res)
  mat2  <- get_mat_from_hic(hic2,norm,chr,res)
  ratio <- sum(mat1)/sum(mat2)
  mat2  <- mat2 * ratio
  mat2  <- triu(mat2,k=1)
  mat1  <- tril(mat1,k=-1)

  sta   <- region[1]/res + 1
  end   <- region[2]/res + 1
  mat1  <- mat1[sta:end,sta:end]
  mat2  <- mat2[sta:end,sta:end]

  mat   <- mat1 + mat2
}

get_mat_from_hic <- function(hic_file,norm='KR',chr,res,
                             output = 'matrix'){
  # get string
  straw_string <- paste(norm, hic_file, chr, chr, 'BP', format(res, scientific = FALSE), sep=' ')
  # get matrix
  #res_int <- as.numeric(res)
  mat <- straw_R(straw_string)
  mat$counts[is.na(mat$counts)] <- 0
  # avoid 0 matrix
  if(max(mat$counts) == 0){
    straw_string <- paste('VC', hic_file, chr, chr, 'BP', format(res, scientific = FALSE), sep=' ')
    print("use VC normalization instead !")
    mat <- straw_R(straw_string)
    mat$counts[is.na(mat$counts)] <- 0
  }

  if(output == 'matrix'){
    mat <- with(mat,sparseMatrix(i=x/res + 1, j=y/res + 1, x=counts, symmetric = TRUE))
  }else if(output == 'dataframe'){
    # mat <- as.data.frame(cbind(mat$x/res + 1, mat$y/res + 1, mat$counts))
    mat <- with(mat,data.frame(x=x/res + 1, y=y/res + 1, counts=counts))
  }else {
    stop("can't recongnize argument: output")
  }
}


