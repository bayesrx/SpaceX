## N: Total number of locations
## G: Number of Genes
## Gene_expression_mat: Gene expression matrix (N X G) [Should be provided as dataframe]
## Spatial_locations: Spatial locations and cluster annotations. [Should be provided as dataframe]
## (x,y) coordinates on the first two columns and
## cluster annotations on the third column.
## Dependency: doParallel,PQLseq MSFA

#' SpaceX function
#'
#' DESCRIPTION
#'
#' @param Gene_expression_mat Gene expression matrix (N X G) [Should be provided as dataframe]
#' @param Spatial_locations Spatial locations and cluster annotations. [Should be provided as dataframe]
#' (x,y) coordinates on the first two columns and
#' cluster annotations on the third column.
#'
#' @return A LIST
#'

SpaceX <- function(Gene_expression_mat,Spatial_locations){

#### Global Parameters ######h
G <-dim(Gene_expression_mat)[2]
L <- length(unique(Spatial_locations[,3]))
Clusters <- unique(Spatial_locations[,3])
N_l <- numeric()
sigma1_sq_est <- matrix(0,G,L)
sigma2_sq_est <- matrix(0,G,L)
u <-list()
Z_est <- list() ##latent gene expression matrix

### Cluster sizes ####
for (l in 1:L) {
  pos <- which(Spatial_locations[,3] == Clusters[l])
  N_l[l] <- length(pos)
  Z_est[[l]] <- matrix(0,nrow = N_l[l],ncol = G)
}


### Poisson mixed model with PQLSEQ algorithm
for (l in 1:L) {

  pos <- which(Spatial_locations[,3] == Clusters[l])

  ### Rho estimation ###
  a <- dist(Spatial_locations[pos,-3])
  a_max <- log10(2*max(a))
  a_min <- log10(min(a)/2)
  a_seq <- seq(a_min, a_max, length.out = 10)
  rho_l <- 10^(a_seq[5])

  Y_mat <- as.matrix(Gene_expression_mat[pos,], rownames.force = F)
  colnames(Y_mat) <- NULL
  location <- Spatial_locations[pos,-3]

  cov_kernel_l <- matrix(0,N_l[l],N_l[l])
  for (i in 1:N_l[l]) {
    for (j in 1:i) {
      dist_loc <- (Spatial_locations[i,1] - Spatial_locations[j,1])^2 + (Spatial_locations[i,2] - Spatial_locations[j,2])^2
      cov_kernel_l[i,j] <- exp(-dist_loc/(2*rho_l^2))
      cov_kernel_l[j,i] <- cov_kernel_l[i,j]
    }}

  fit <- pqlseq_modified(RawCountDataSet=t(Y_mat),Phenotypes= rep(1,N_l[l]), RelatednessMatrix = cov_kernel_l,
                   fit.model="PMM", numCore = 1)

  j = 1 + (0:(G-1))*N_l[l]
  sigma1_sq_est[,l] <- fit$tau1[j]
  sigma2_sq_est[,l] <- fit$tau2[j]
  u[[l]] <- matrix(fit$residual, nrow = G, byrow = TRUE)

  ## Estimation of latent gene expression
  for (g in 1:G) {

    if(sigma1_sq_est[g,l]==0){
      V <- (sigma1_sq_est[g,l]+0.001)*cov_kernel_l + (sigma2_sq_est[g,l])*diag(N_l[l])
    }
    else{
      V <- (sigma1_sq_est[g,l])*cov_kernel_l + (sigma2_sq_est[g,l])*diag(N_l[l])
    }

    if(sigma2_sq_est[g,l]==0){
      Z_est[[l]][,g] <- ((sigma2_sq_est[g,l]+0.001)*solve(V))%*%u[[l]][g,]
    }
    else{
      Z_est[[l]][,g] <- ((sigma2_sq_est[g,l])*solve(V))%*%u[[l]][g,]
    }


  }
  print(l)
}

## Applying multi-study factor model on latent gene expression matrix
print("Multi-Study Factor Model")
fit_MSFA = sp_msfa(Z_est,  k = 10,  j_s = rep(10,L), trace = FALSE)

return(fit_MSFA)

## SigmaPhi :: Shared Covariance matrix
## SigmaLambda :: Cluster specific Covaraince matrices
}









