#' @title Get permuted groupwise correlations and pairwise differential correlations. 
#' @description Takes input and methods and randomly permutes the data to do getCor as well as group-specific pairwiseDCor.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param inputMatB Optional, secondary input matrix that allows you to calculate correlation and differential correlation for the rows between inputMat and imputMatB.
#' @param compare Vector of two character strings, each corresponding to one name in the list of correlation matrices that should be compared.
#' @param impute A binary variable specifying whether values should be imputed if there are missing values. Note that the imputation is performed in the full input matrix (i.e., prior to subsetting) and uses k-nearest neighbors.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman".
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param nPerms Number of permutations to generate.
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes. Default = 0.99
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none"). This could be used if you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @return An array of permuted differences in z-scores calculated between conditions, with the third dimension corresponding to the number of permutations performed.
#' @export
getDCorPerm <- function(inputMat, design, compare, inputMatB = NULL, impute = FALSE,
	nPerms = 10, corrType = "pearson", corr_cutoff = 0.99, signType = "none"){

	secondMat = FALSE
	if(!is.null(inputMatB)){
		zPermMat = array(dim = c(nrow(inputMat), nrow(inputMatB), nPerms))
		secondMat = TRUE
	} else {
		zPermMat = array(dim = c(nrow(inputMat), nrow(inputMat), nPerms))
	}

	for(i in 1:nPerms){
		message("Calculating permutation number ", i, ".")
		inputMat_perm = inputMat[ , sample(ncol(inputMat)), drop = FALSE]
		if(secondMat){
			inputMatB_perm = inputMatB[ , sample(ncol(inputMatB)), drop = FALSE]
			corMats_res = getCors(inputMat_perm, design = design,
				inputMatB = inputMatB_perm, corrType = corrType, impute = impute)
		} else {
			corMats_res = getCors(inputMat_perm, design = design,
				corrType = corrType, impute = impute)
		}
		
		corMats_bootstrap1 <- list()
		corMats_bootstrap2 <- list()
		
		col_index1 <- which(colnames(as.matrix(design)) == compare[1])
		row_indices1 <- which(as.matrix(design)[, col_index1] == 1)
		
		col_index2 <- which(colnames(as.matrix(design)) == compare[2])
		row_indices2 <- which(as.matrix(design)[, col_index2] == 1)
		
		inputMat1 <- inputMat_perm %>%
		  `colnames<-`(NULL)
		inputMat2 <- inputMat_perm %>%
		  `colnames<-`(NULL)
		
		fisher_z <- function(r) {
		  r_capped <- pmin(pmax(r, -0.999), 0.999)
		  return(0.5 * log((1 + r_capped) / (1 - r_capped)))
		}
		
		for (j in 1:100) {
		  # Resample with replacement
		  
		  inputMat2[, row_indices1] <- inputMat1[,sample(row_indices1, size = length(row_indices1), replace = TRUE)]
		  inputMat2[, row_indices2] <- inputMat1[,sample(row_indices2, size = length(row_indices2), replace = TRUE)]
		  
		  # Store bootstrap z value
		  corMats_bootstrap <- getCors(inputMat2, design = design,
		                                  corrType = corrType, impute = impute)
		  transformed_1 <- matrix(sapply(corMats_bootstrap@corMatList[[1]]$corrs, fisher_z),
		                          nrow=nrow(corMats_bootstrap@corMatList[[1]]$corrs), 
		                          ncol=ncol(corMats_bootstrap@corMatList[[1]]$corrs)) %>%
		    `colnames<-`(rownames(corMats_bootstrap@corMatList[[1]]$corrs)) %>%
		    `rownames<-`(rownames(corMats_bootstrap@corMatList[[1]]$corrs))
		  
		  transformed_2 <- matrix(sapply(corMats_bootstrap@corMatList[[2]]$corrs, fisher_z),
		                          nrow=nrow(corMats_bootstrap@corMatList[[2]]$corrs), 
		                          ncol=ncol(corMats_bootstrap@corMatList[[2]]$corrs)) %>%
		    `colnames<-`(rownames(corMats_bootstrap@corMatList[[2]]$corrs)) %>%
		    `rownames<-`(rownames(corMats_bootstrap@corMatList[[2]]$corrs))
		  
		  corMats_bootstrap1[[j]] <- transformed_1
		  corMats_bootstrap2[[j]] <- transformed_2
		}
		
		Get_variance <- function(matrix_list,nrows){
		  # Calculate the variance for each (i, j) position
		  variance_matrix <- matrix(0,nrow = nrows, ncol = nrows)
		  for (m in 1:nrows) {
		    for (n in 1:nrows) {
		      # Extract the (i, j) elements from all matrices
		      elements <- sapply(matrix_list, function(x) x[m, n])
		      elements <- elements[is.finite(elements)]
		      
		      # Calculate the variance of the extracted elements
		      variance_matrix[m, n] <- var(elements, na.rm = T)
		    }
		  }
		  variance_matrix <- variance_matrix %>%
		    `colnames<-`(colnames(matrix_list[[1]])) %>%
		    `rownames<-`(rownames(matrix_list[[1]]))
		  return(variance_matrix)
		}
		
		Variance_1 <- Get_variance(corMats_bootstrap1,dim(corMats_bootstrap1[[1]])[1])
		Variance_2 <- Get_variance(corMats_bootstrap2,dim(corMats_bootstrap1[[1]])[1])
		
		Permut_Variance_lst <- list(
		  Variance_1 = Variance_1,
		  Variance_2 = Variance_2
		)

		dcPairs_res = pairwiseDCor(corMats_res, compare, corr_cutoff = corr_cutoff,
			secondMat = secondMat, signType = signType,Permut_Variance_lst = Permut_Variance_lst)
		zscores = slot(dcPairs_res, "ZDiff")
		zPermMat[ , , i] = zscores
	}

	return(zPermMat)

}
