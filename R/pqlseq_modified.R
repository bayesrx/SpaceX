##################################################################################
##
## The pqlseq algorithm is developed in the following paper:
## Title  : Heritability Estimation and Differential Analysis with Generalized Linear Mixed Models in Large-Scale Genomic Sequencing Studies
## Authors: Shiquan Sun, Jiaqiang Zhu, and Xiang Zhou: Package: PQLseq.
## This is a modified version of the pqlseq algorithm particularly for the SpaceX package.
##
##################################################################################

#' Fit Generalized Linear Mixed Model with Known Kinship Matrices Through Penalized-quasi Likelihood
#'
#' Fit a generalized linear mixed model with a random intercept. The covariance matrix of the random intercept is proportional to a known kinship matrix. This is a modified version of the pqlseq algorithm particularly for the SpaceX package. For more details check the pqlseq function from PQLseq package.
#'
#' @param RawCountDataSet a data frame containing the read count.
#' @param Phenotypes a vector containing the predictor of interest.
#' @param Covariates a data frame containing the covariates subject to adjustment (Default = NULL).
#' @param RelatednessMatrix a known relationship matrix (e.g. kinship matrix in genetic studies). When supplied with a matrix, this matrix should be a positive semi-definite matrix with dimensions equal to the sample size in count data, and the order of subjects in this matrix should also match the order of subjects in count data. Currently there is no ID checking feature implemented, and it is the user's responsibility to match the orders.
#' @param LibSize a data frame containing the total read count. For possion mixed model, it will be calculated automatically if users do not provide. For binomial mixed model, it is required.
#' @param fit.model a description of the error distribution and link function to be used in the model. Either "PMM" for possion model, or "BMM" for binomial model (default = "PMM").
#' @param fit.method method of fitting the generalized linear mixed model, currently only "REML" version is available.
#' @param fit.maxiter a positive integer specifying the maximum number of iterations when fitting the generalized linear mixed model (default = 500).
#' @param fit.tol a positive number specifying tolerance, the difference threshold for parameter estimates below which iterations should be stopped (default = 1e-5).
#' @param numCore a positive integer specifying the number of cores for parallel computing (default = 1).
#' @param filtering a logical switch for RNAseq data. By default, for each gene, at least two individuals should have read counts greater than 5. Otherwise, the gene is filtered (default = TRUE).
#' @param verbose a logical switch for printing detailed information (parameter estimates in each iteration) for testing and debugging purpose (default = FALSE).
#' @param ... additional arguments that could be passed to glm.
#'
#' @return
#' \item{numIDV}{number of individuals with data being analyzed}
#' \item{beta}{the fixed effect parameter estimate for the predictor of interest.}
#' \item{se_beta}{the standard deviation of fixed effect.}
#' \item{pvalue}{P value for the fixed effect, based on the wald test.}
#' \item{h2}{heritability of the transformed rate.}
#' \item{sigma2}{total variance component.}
#' \item{overdisp}{dispersion parameter estimate.}
#' \item{converged}{a logical indicator for convergence.}
#'
#' @references Sun, S., Hood, M., Scott, L., Peng, Q., Mukherjee, S., Tung, J., and Zhou, X. (2017). Differential expression analysis for rnaseq using poisson mixed models. Nucleicacids research, 45(11), e106–e106.
#'

pqlseq_modified <- function(RawCountDataSet, Phenotypes, Covariates=NULL, RelatednessMatrix=NULL, LibSize=NULL,
                      fit.model="PMM", fit.method = "AI.REML", fit.maxiter=500, fit.tol=1e-5, numCore=1,
                      filtering=TRUE, verbose=FALSE, ...) {
  # specify the number of cores we want to use
  if(numCore > 1){
    if(numCore>detectCores()){warning("PQLseq:: the number of cores you're setting is larger than detected cores!");numCore = detectCores()-1}
  }

  registerDoParallel(numCore)

  # cl <- makeCluster(numCore)
  # registerDoParallel(cl,cores=numCore)
  # on.exit(stopCluster(cl))

  # filtering genes/sites
  #if (filtering & fit.model == "PMM"){
  #  unfilterIdx <- apply(RawCountDataSet, 1, function(x){length(x[x>5])>=2} )
  #  CountData   <- RawCountDataSet[unfilterIdx,]
  #}else{
  #  CountData   <- RawCountDataSet
  #}
  CountData   <- RawCountDataSet
  rm(RawCountDataSet)

  numVar <- dim(CountData)[1]
  numIDV <- dim(CountData)[2]

  # remove the intercept
  if(length(unique(Covariates[,1])) == 1){
    Covariates<- Covariates[,-1]
  }

  if(is.null(Covariates)){
    numCov <- 0
  }else{
    numCov     <- dim(Covariates)[2]
    Covariates <- as.matrix(Covariates)
  }

  cat(paste("## number of total individuals: ", numIDV,"\n"))
  cat(paste("## number of total genes/sites: ", numVar,"\n"))
  cat(paste("## number of adjusted covariates: ", numCov,"\n"))


  CountData  <- as.matrix(CountData)
  Phenotypes <- as.matrix(Phenotypes)


  # if(is.null(RelatednessMatrix)){
  #   stop("PQLseq::please input relatedness matrix!")
  # }else{
  #   RelatednessMatrix <- as.matrix(RelatednessMatrix)
  #   scalerM           <- diag(numIDV)-(rep(1,numIDV)%*%t(rep(1,numIDV)))/numIDV
  #   eig               <- eigen(RelatednessMatrix)
  #   eigval            <- eig$value
  #   eigvector         <- eig$vectors
  #   if(any(eigval<1e-10)){
  #     warning("PQLseq::the relatedness matrix is singular, it has been modified!")
  #     RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix,corr=T)$mat)
  #   }
  #   rm(scalerM)
  #   rm(eig)
  #   rm(eigval)
  #   rm(eigvector)
  # }

  RelatednessMatrix <- list(RelatednessMatrix, diag(numIDV))

  #***********************************#
  #       Poisson Mixed Model         #
  #***********************************#
  if(fit.model == "PMM"){
    cat("# fitting Poisson mixed model ... \n")
    if(is.null(LibSize)){
      LibSize <- apply(CountData, 2, sum)
      LibSize <- as.matrix(LibSize)
    }else{
      LibSize <- as.matrix(t(LibSize))
    }


    # do parallel using foreach function
    iVar   <- NULL
    resPMM <-foreach(iVar=1:numVar,.combine=rbind)%dopar%{
      numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
      if(numCov==0){
        model0 <- try(glm(formula = CountData[iVar,]~1 + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = CountData[iVar,]~1 + offset(log(LibSize)), na.action = na.omit)),
                       rownames(model.frame(formula = CountData[iVar,]~1 + offset(log(LibSize)), na.action = na.pass)))
      }else{
        model0 <- try(glm(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.omit)),
                       rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.pass)))
      }

      if(verbose) {cat(paste("NO. Gene = ",iVar,"\n"))}

      tmpRelatednessMatrix <- RelatednessMatrix
      if(class(tmpRelatednessMatrix) == "matrix") {
        tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
      }else {
        for(ik in seq_len(length(tmpRelatednessMatrix)) ) {tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]}
      }

      names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")

      if(class(model0)[1]!="try-error"){
        # t1 <- system.time(model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix)))
        model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix))
      }else{
        model1 <- NULL
      }

      if(!is.null(model1)&(class(model1)!="try-error")){
        if(verbose){cat(paste("PQLseq::PMM::tau = ", model1$theta,"\n"))}
        numAnalysis <- length(idx)
        beta        <- model1$coefficients[length(model1$coefficients)]
        alpha       <- model1$coefficients[1]
        se_beta     <- sqrt(diag(model1$cov)[length(model1$coefficients)] )
        pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
        sigma2      <- model1$theta[2]+model1$theta[3]
        h2          <- model1$theta[2]/(sigma2)
        tau1        <- model1$theta[2]
        tau2        <- model1$theta[3]
        residual    <- model1$residuals
        fitted_values <- model1$fitted.values
        converged   <- model1$converged
      }else{converged <- FALSE}

      res <- data.frame(numIDV = numAnalysis, beta = beta, alpha=alpha, se_beta = se_beta,
                        pvalue = pvalue, h2 = h2, sigma2 = sigma2,tau1=tau1,tau2=tau2,
                        fitted_values = fitted_values,
                        residual = residual, converged = converged)
    }# end for iVar, parallel
    rm(iVar)
    closeAllConnections()
    # if(nrow(showConnections())!=0){closeAllConnections()}

    rownames(resPMM) <- rownames(CountData)
    return(resPMM)
  }# end PMM
  #***********************************#
  #       Binomial Mixed Model        #
  #***********************************#
  if(fit.model == "BMM"){
    cat("# fitting binomial mixed model ... \n")
    if(is.null(LibSize)){
      stop("PQLseq::BMM::ERROR: please input the LibSize (total counts) file!!")
    }else{
      LibSize <- as.matrix(LibSize)
    }

    ratio               <- CountData/LibSize
    ratio[is.na(ratio)] <- 0
    flag                <- ratio>1.0
    sumflag             <- apply(flag,1, sum)
    idx                 <- which(sumflag>0)

    if (length(idx)>0){
      CountData <- CountData[-idx,]
      LibSize   <- LibSize[-idx,]
    }else{
      CountData <- CountData
      LibSize   <- LibSize
    }

    numVar <- dim(CountData)[1]
    numIDV <- dim(CountData)[2]
    iVar   <- NULL

    # do parallel
    resBMM <- foreach(iVar=1:numVar,.combine=rbind)%dopar%{
      numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
      if(verbose){cat(paste("NO. Gene/Site = ",iVar,"\n"))}
      if(sum(dim(LibSize)==dim(CountData)) != 2){
        stop("PQLseq::BMM::ERROR: the dimensions of read counts and total read counts do not match!")
      }

      LibSize <- as.matrix(LibSize)

      if(numCov == 0){
        model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,])
        idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.omit)),
                        rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.pass)))
      }else{
        model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,] )
        idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.omit)),
                        rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.pass)))
      }

      model0$numTotal <- LibSize[iVar,idx]
      model0$numSucc  <- CountData[iVar,idx]

      redflag <- FALSE
      for( ierr in c(2:dim(model.matrix(model0))[2])){
        if(length(unique(model.matrix(model0)[,ierr])) == 1){
          warning(paste("PQLseq::BMM::the ",ierr-1,"-th column of covariates are the same for gene/site ",rownames(CountData)[iVar],"!",sep = "") )
          redflag <- TRUE
        }
      }
      if(!redflag){

        tmpRelatednessMatrix <- RelatednessMatrix
        if(class(tmpRelatednessMatrix) == "matrix") {
          tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
        }else {
          for(ik in seq_len(length(tmpRelatednessMatrix)) ) {
            tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]
          }
        }
        names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")

        # t1 <- system.time(model1 <- try( PQLseq.fit(model0, tmpRelatednessMatrix) ))
        model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix))

        if(class(model1) != "try-error"&!is.null(model1)){
          if(verbose){cat(paste("PQLseq::BMM::tau = ", model1$theta,"\n"))}
          numAnalysis <- length(idx)
          beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
          se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )
          pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
          sigma2      <- model1$theta[2]+model1$theta[3]
          h2          <- model1$theta[2]/(sigma2)
          tau1        <- model1$theta[2]
          tau2        <- model1$theta[3]
          converged   <- model1$converged
        }else{converged <- FALSE}

        res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta,
                          pvalue = pvalue, h2 = h2, sigma2 = sigma2,
                          converged = converged)
      }# end for iVar, parallel

    }
    rm(iVar)

    # if(nrow(showConnections())!=0){closeAllConnections()}
    closeAllConnections()
    rownames(resBMM) <- rownames(CountData)
    return(resBMM)
  }# end BMM

}# end function PQLseq


##########################################################
#           	   PQLseq FIT FUNCTION					 #
##########################################################

PQLseq.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {

  names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
  # if((method.optim == "AI")&(!sum(model0$fitted.values<1e-5))) {
  if(method.optim == "AI") {
    fixtau.old 	<- rep(0, length(RelatednessMatrix)+1)
    # to use average information method to fit alternative model
    model1 		<- PQLseq.AI(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new 	<- 1*(model1$theta < 1.01 * tol)

    while(any(fixtau.new != fixtau.old)) {
      fixtau.old <- fixtau.new
      model1 	<- PQLseq.AI(model0, RelatednessMatrix, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(model1$theta < 1.01 * tol)
    }
  }else{
    model1 <- NULL
  }
  return(model1)
}

##########################################################
#       PQLseq FIT AVERAGE INFORMATION FUNCTION			 #
##########################################################

PQLseq.AI <- function(model0, RelatednessMatrix, tau = rep(0, length(RelatednessMatrix)+1), fixtau = rep(0, length(RelatednessMatrix)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {

  if(model0$family$family %in% c("binomial")){
    y <- model0$numSucc
  }else{
    y <- model0$y
  }
  numIDV <- length(y)
  offset <- model0$offset
  if(is.null(offset)) {offset <- rep(0, numIDV)}

  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  mu.eta <- family$mu.eta(eta)
  D <- mu.eta/sqrt(model0$family$variance(mu))

  if(family$family %in% c("binomial")){
    mu.eta <- model0$numTotal*mu.eta
    D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
    mu <- model0$numTotal*mu
  }

  Y <- eta - offset + (y - mu)/mu.eta
  X <- model.matrix(model0)
  alpha <- model0$coef

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  numK <- length(RelatednessMatrix)
  idxtau <- which(fixtau == 0)
  numK2 <- sum(fixtau == 0)

  ### this part needs to be changed for intercept only model same as spark (Satwik)
  if(numK2 > 0) {
    tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)

    H <- tau[1]*diag(1/D^2)
    for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

    Hinv 	<- chol2inv(chol(H))
    HinvX 	<- crossprod(Hinv, X)
    XHinvX 	<- crossprod(X, HinvX)

    P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))

    if(class(P) == "try-error"){
      stop("Error in P matrix calculation!")
    }

    PY <- crossprod(P, Y)
    tau0 <- tau
    for(ik in 1:numK2) {
      if(ik == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
      else {
        PAPY <- crossprod(P, crossprod(RelatednessMatrix[[idxtau[ik]-1]], PY))
        tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(P*RelatednessMatrix[[idxtau[ik]-1]]))/numIDV)
      }
    }
  }

  for (iter in seq_len(maxiter)) {
    alpha0 	<- alpha
    tau0 	<- tau
    model1 	<- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)

    tau <- as.numeric(model1$tau)
    cov <- as.matrix(model1$cov)
    alpha <- as.numeric(model1$alpha)
    eta <- as.numeric(model1$eta) + offset


    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(family$variance(mu))

    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
      mu <- model0$numTotal*mu
    }

    Y <- eta - offset + (y - mu)/mu.eta

    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
    if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {

      iter <- maxiter
      break
    }
  }

  converged <- ifelse(iter < maxiter, TRUE, FALSE)
  res <- y - mu
  P <- model1$P
  return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, P = P, residuals = res, cov = cov, converged = converged))
}# end function


# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

AI <- function(Yin, Xin, numKin, Phiin, Din, tauin, fixtauin, tolin) {
  .Call('_PQLseq_AI', PACKAGE = 'PQLseq', Yin, Xin, numKin, Phiin, Din, tauin, fixtauin, tolin)
}

rcpparma_hello_world <- function() {
  .Call('_PQLseq_rcpparma_hello_world', PACKAGE = 'PQLseq')
}

rcpparma_outerproduct <- function(x) {
  .Call('_PQLseq_rcpparma_outerproduct', PACKAGE = 'PQLseq', x)
}

rcpparma_innerproduct <- function(x) {
  .Call('_PQLseq_rcpparma_innerproduct', PACKAGE = 'PQLseq', x)
}

rcpparma_bothproducts <- function(x) {
  .Call('_PQLseq_rcpparma_bothproducts', PACKAGE = 'PQLseq', x)
}


#########################################
#             CODE END                  #
#########################################
