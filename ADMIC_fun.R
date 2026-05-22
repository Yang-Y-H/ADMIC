library(stats)
library(glasso)
# data is an abundance matrix where the rows correspond to samples and the columns correspond to taxa. 
# Y is a vector of trait
# C is a matrix of confounders. The default is NULL.
# c is the L1 regularization penalty parameter used to control the precision matrix, it determines the sparsity of the precision matrix. 
#    c=2: The network is extremely sparse, with no erroneous edges. c=1: Theoretical equilibrium point. c=0.5: Allows for some noise, resulting in a denser network. The default value is 1.
# pseudo.count is the non-zero pseudo-count used to replace 0; it is recommended to use values such as 0.5 or 1. The default is 1.
# r.test is the parameter range for the Box-Cox transformation (must be set between 0 and 1); the default range is 0.1 to 0.9, with a step size of 0.1.

ADMIC <- function(data, Y, C = NULL, 
                  pseudo.count = 1, c = 1,
                  r.test = seq(from = 0.1, to = 0.9, by = 0.1)){
  
  if (!is.numeric(data)) {
    stop("Data Matrix must be numeric.")
  }
  
  if (
    any(is.na(data)) ||
    any(!is.finite(data)) ||
    any(data < 0)
  ) {
    stop("There is an issue with the Data Matrix, please correct it.")
  }
  if (is.data.frame(data)) { data <- data.matrix(data) }
  if (is.data.frame(Y)) { Y <- data.matrix(Y) }
  if (is.data.frame(C)) { C <- data.matrix(C) }
  
  # filtering out low prevalence taxa
  taxaname <- colnames(data)
  prop.presence <- colMeans(data > 0)
  remove.taxa <- which(prop.presence < 0.05)
  if(length(remove.taxa) > 0) {
    cat("The prevalence of",length(remove.taxa)," taxa is below the threshold and will be filtered out.")
    data <- data[ ,-remove.taxa]
  }
  
  n.tax <- ncol(data)
  n.sam <- nrow(data)
  raw.data <- data
  data[which(data==0)] <- pseudo.count
  # step 1 --- select reference taxa
  ref <- reference.select(raw.data, Y, pseudo.count = pseudo.count)
  
  # step 2 --- precision matrix 
  geo.mean <- apply(data[,ref], 1, mean)
  trans.data.raw <- data/geo.mean
  rho <- c*sqrt(log(n.tax)/n.sam)
  Precision.mat <- bonobo.precision(t(data), rho = rho)
  
  # step 3 --- main test model
  n <- length(r.test)
  AIC0 <- numeric(n)
  AIC.eff <- numeric(n)
  P0 <- matrix(0, nrow = n, ncol = n.tax)
  P.eff <- matrix(0, nrow = n, ncol = n.tax)
  for(r in 1:n){
    trans.data <- ((trans.data.raw)^r.test[r]-1)/r.test[r]
    eff <- matrix(0, n.sam, n.tax)
    for(i in 1:n.tax){
      for(j in 1:n.sam){
        eff[j,i] <- sum(Precision.mat[[j]][-i,i]*trans.data[j,-i])/(-Precision.mat[[j]][i,i])
      }
    }
    aic0 <- numeric(n.tax)
    aic1 <- numeric(n.tax)
    for(i in 1:n.tax){
      X <- trans.data[,i]
      effectsize <- eff[,i]# interaction effect
      if(is.null(C)){
        df <- data.frame(X = X, Y = Y, effectsize = effectsize)
        formula_lm0 <- as.formula("X ~ Y")
        formula_lm1 <- as.formula("X ~ Y + effectsize")
      }else{
        C <- as.data.frame(C)
        df <- data.frame(X = X, Y = Y, effectsize = effectsize, C)
        formula_lm0 <- as.formula(
          paste("X ~ Y +", paste(colnames(C), collapse = " + "))
        )
        formula_lm1 <- as.formula(
          paste("X ~ Y + effectsize +", paste(colnames(C), collapse = " + "))
        )
      }
      
      #ADMIC-Basic
      model0 <- lm(formula_lm0, data = df)
      P0[r,i] <- summary(model0)[["coefficients"]][2,4]
      aic_raw0 <- AIC(model0)
      aic0[i] <- aic_raw0 - 2 * (r.test[r] - 1) * sum(log(trans.data.raw[,i]))
      
      #ADMIC
      model1 <- lm(formula_lm1, data = df)
      P.eff[r,i] <- summary(model1)[["coefficients"]][2,4]
      aic_raw1 <- AIC(model1)
      aic1[i] <- aic_raw1 - 2 * (r.test[r] - 1) * sum(log(trans.data.raw[,i]))
    }
    AIC0[r] <- mean(aic0)
    AIC.eff[r] <- mean(aic1)
  }
  maxrho0 <- which.min(AIC0)
  maxrho.eff <- which.min(AIC.eff)
  p0 <- P0[maxrho0,]
  p.eff <- P.eff[maxrho.eff,]
  q0 <- p.adjust(p0, method = "BH")
  q.eff <- p.adjust(p.eff, method = "BH")
  
  res <- list(feature = colnames(data), p.basic = p0, q.basic =q0, r.basic = r.test[maxrho0],
              p.admic = p.eff, q.admic = q.eff, r.admic = r.test[maxrho.eff])
  return(res)
}





