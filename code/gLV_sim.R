#---------load function---------
sumresult <- function(qvalue, causal.otus, noncausal.otus, fdr.target=0.01) {
  
  otu.detected = sort(which(qvalue < fdr.target))
  n.otu = length(otu.detected)
  
  if (n.otu > 0) {
    sen = sum(otu.detected %in% causal.otus)/length(causal.otus)
    fdr = n.otu - sum(otu.detected %in% causal.otus)
    fdr = fdr/n.otu
  } else {
    sen = 0
    fdr = 0
  }
  ord=order(qvalue)[1:20]
  out = list(otu.d=otu.detected,otu.top=ord, sen=sen,fdr=fdr)
  return(out)
}
controlRecall <- function(qvalue, causal.otus, noncausal.otus, target=1) {
  cut <- ceiling(target*length(causal.otus))
  fdr.target <- sort(qvalue[causal.otus])[cut]
  otu.detected <- sort(which(qvalue <= fdr.target))
  n.otu <- length(otu.detected)
  if (n.otu > 0) {
    fdr = n.otu - sum(otu.detected %in% causal.otus)
    fdr = fdr/n.otu
  } else {
    sen = 0
    fdr = 1
  }
  return(fdr)
}
add_poisson_noise <- function(x_glv, depth = 150) {
  # x_glv: matrix / data.frame / vector
  n <- nrow(x_glv)
  depth.sam <- rnorm(n, mean = depth, sd = depth/5)
  x <- as.matrix(x_glv)
  
  if (any(x < 0)) {
    stop("gLV abundance contains negative values")
  }
  
  lambda <- x * depth.sam
  x_poisson <- matrix(
    rpois(length(lambda), lambda = lambda),
    nrow = nrow(x),
    ncol = ncol(x)
  )
  
  colnames(x_poisson) <- colnames(x)
  rownames(x_poisson) <- rownames(x)
  
  return(x_poisson)
}
ANCOM_BC <- function(feature.table, grp.name, grp.ind, struc.zero, adj.method = "BH", 
                     tol.EM = 1e-5, max.iterNum = 100, perNum = 1000, alpha = 0.05){
  n.taxa.raw = nrow(feature.table)
  taxa.id.raw = rownames(feature.table)
  n.samp = ncol(feature.table)
  sample.id = colnames(feature.table)
  
  n.grp = length(grp.ind)
  n.samp.grp = sapply(grp.ind, length)
  
  ### 0. Discard taxa with structural zeros for the moment
  comp.taxa.pos = which(apply(struc.zero, 1, function(x) all(x == 0))) # position of complete taxa (no structural zeros)
  O = feature.table[comp.taxa.pos, ]
  n.taxa = nrow(O)
  taxa.id = rownames(O)
  n.samp = ncol(O)
  y = as.matrix(log(O+1))
  
  ### 1. Initial estimates of sampling fractions and mean absulute abundances
  mu = t(apply(y, 1, function(i) tapply(i, rep(1:n.grp, n.samp.grp), function(j)
    mean(j, na.rm = T))))
  d = colMeans(y - mu[, rep(1:n.grp, times = n.samp.grp)], na.rm = T)
  
  ## Iteration in case of missing values of y
  iterNum = 0
  epsilon = 100
  while (epsilon > tol.EM & iterNum < max.iterNum) {
    # Updating mu
    mu.new = t(apply(t(t(y) - d), 1, function(i) tapply(i, rep(1:n.grp, n.samp.grp), function(j)
      mean(j, na.rm = T))))
    
    # Updating d
    d.new = colMeans(y - mu.new[, rep(1:ncol(mu.new), times = n.samp.grp)], na.rm = T)
    
    # Iteration
    epsilon = sqrt(sum((mu.new - mu)^2) + sum((d.new - d)^2))
    iterNum = iterNum + 1
    
    mu = mu.new
    d = d.new
    if(is.nan(epsilon)){epsilon=0}
  }
  
  mu.var.each = (y-t(t(mu[, rep(1:ncol(mu), times = n.samp.grp)])+d))^2
  mu.var = t(apply(mu.var.each, 1, function(x) tapply(x, rep(1:n.grp, n.samp.grp), function(y)
    mean(y, na.rm = T))))
  sample.size = t(apply(y, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) length(y[!is.na(y)])))))
  mu.var = mu.var/sample.size
  
  ### 2. Estimate the bias (between-group difference of sampling fractions) by E-M algorithm
  bias.em.vec = rep(NA, n.grp - 1)
  bias.wls.vec = rep(NA, n.grp - 1)
  bias.var.vec = rep(NA, n.grp - 1)
  for (i in 1:(n.grp-1)) {
    Delta = mu[, 1] - mu[, 1+i]
    nu = rowSums(mu.var[, c(1, 1+i)])
    
    ## 2.1 Initials
    pi0_0 = 0.75
    pi1_0 = 0.125
    pi2_0 = 0.125
    delta_0 = mean(Delta[Delta >= quantile(Delta, 0.25, na.rm = T)&
                           Delta <= quantile(Delta, 0.75, na.rm = T)], na.rm = T)
    l1_0 = mean(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    l2_0 = mean(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    kappa1_0 = var(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
    kappa2_0 = var(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    if(is.na(kappa2_0)|kappa2_0 == 0) kappa2_0 = 1
    
    ## 2.2 Apply E-M algorithm
    # 2.21 Store all paras in vectors/matrices
    pi0.vec = c(pi0_0); pi1.vec = c(pi1_0); pi2.vec = c(pi2_0)
    delta.vec = c(delta_0); l1.vec = c(l1_0); l2.vec = c(l2_0)
    kappa1.vec = c(kappa1_0); kappa2.vec = c(kappa2_0)
    
    # 2.22 E-M iteration
    iterNum = 0
    epsilon = 100
    while (epsilon > tol.EM & iterNum < max.iterNum) {
      # print(iterNum)
      ## Current value of paras
      pi0 = pi0.vec[length(pi0.vec)]; pi1 = pi1.vec[length(pi1.vec)]; pi2 = pi2.vec[length(pi2.vec)]
      delta = delta.vec[length(delta.vec)]; 
      l1 = l1.vec[length(l1.vec)]; l2 = l2.vec[length(l2.vec)]
      kappa1 = kappa1.vec[length(kappa1.vec)]; kappa2 = kappa2.vec[length(kappa2.vec)]
      
      ## E-step
      pdf0 = sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta, sqrt(nu[i])))
      pdf1 = sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta + l1, sqrt(nu[i] + kappa1)))
      pdf2 = sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta + l2, sqrt(nu[i] + kappa2)))
      r0i = pi0*pdf0/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r0i[is.na(r0i)] = 0
      r1i = pi1*pdf1/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r1i[is.na(r1i)] = 0
      r2i = pi2*pdf2/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r2i[is.na(r2i)] = 0
      
      ## M-step
      pi0_new = mean(r0i, na.rm = T); pi1_new = mean(r1i, na.rm = T); pi2_new = mean(r2i, na.rm = T)
      delta_new = sum(r0i*Delta/nu + r1i*(Delta-l1)/(nu+kappa1) + r2i*(Delta-l2)/(nu+kappa2), na.rm = T)/
        sum(r0i/nu + r1i/(nu+kappa1) + r2i/(nu+kappa2), na.rm = T)
      l1_new = min(sum(r1i*(Delta-delta)/(nu+kappa1), na.rm = T)/sum(r1i/(nu+kappa1), na.rm = T), 0)
      l2_new = max(sum(r2i*(Delta-delta)/(nu+kappa2), na.rm = T)/sum(r2i/(nu+kappa2), na.rm = T), 0)
      
      # Nelder-Mead simplex algorithm for kappa1 and kappa2
      obj.kappa1 = function(x){
        log.pdf = log(sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+l1, sqrt(nu[i]+x))))
        log.pdf[is.infinite(log.pdf)] = 0
        -sum(r1i*log.pdf, na.rm = T)
      }
      kappa1_new = neldermead(x0 = kappa1, fn = obj.kappa1, lower = 0)$par
      
      obj.kappa2 = function(x){
        log.pdf = log(sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+l2, sqrt(nu[i]+x))))
        log.pdf[is.infinite(log.pdf)] = 0
        -sum(r2i*log.pdf, na.rm = T)
      }
      kappa2_new = neldermead(x0 = kappa2, fn = obj.kappa2, lower = 0)$par
      
      ## Merge to the paras vectors/matrices
      pi0.vec = c(pi0.vec, pi0_new); pi1.vec = c(pi1.vec, pi1_new); pi2.vec = c(pi2.vec, pi2_new)
      delta.vec = c(delta.vec, delta_new)
      l1.vec = c(l1.vec, l1_new); l2.vec = c(l2.vec, l2_new)
      kappa1.vec = c(kappa1.vec, kappa1_new); kappa2.vec = c(kappa2.vec, kappa2_new)
      
      ## Calculate the new epsilon
      epsilon = sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 + (delta_new-delta)^2+
                       (l1_new-l1)^2 + (l2_new-l2)^2 + (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
      iterNum = iterNum+1
      if(is.nan(epsilon)){epsilon=0}
    }
    # 2.23 Estimate the bias
    bias.em.vec[i] = delta.vec[length(delta.vec)]
    
    # 2.24 The WLS estimator of bias
    # Cluster 0
    C0 = which(Delta >= quantile(Delta, pi1_new, na.rm = T) & Delta < quantile(Delta, 1 - pi2_new, na.rm = T))
    # Cluster 1
    C1 = which(Delta < quantile(Delta, pi1_new, na.rm = T))
    # Cluster 2
    C2 = which(Delta >= quantile(Delta, 1 - pi2_new, na.rm = T))
    
    nu_temp = nu
    nu_temp[C1] = nu_temp[C1] + kappa1_new
    nu_temp[C2] = nu_temp[C2] + kappa2_new
    wls.deno = sum(1 / nu_temp)
    
    wls.nume = 1 / nu_temp
    wls.nume[C0] = (wls.nume * Delta)[C0]
    wls.nume[C1] = (wls.nume * (Delta - l1_new))[C1]
    wls.nume[C2] = (wls.nume * (Delta - l2_new))[C2];   
    wls.nume = sum(wls.nume)
    
    bias.wls.vec[i] = wls.nume / wls.deno
    
    # 2.25 Estimate the variance of bias  
    bias.var.vec[i] = 1 / wls.deno
    if (is.na(bias.var.vec[i])) bias.var.vec[i] = 0
  }
  bias.em.vec = c(0, bias.em.vec)
  bias.wls.vec = c(0, bias.wls.vec)
  
  ### 3. Final estimates of mean absolute abundane and sampling fractions
  mu.adj.comp = t(t(mu) + bias.em.vec)
  colnames(mu.adj.comp) = grp.name; rownames(mu.adj.comp) = taxa.id
  
  d.adj = d - rep(bias.em.vec, sapply(grp.ind, length))
  names(d.adj) = sample.id
  
  ### 4. Hypothesis testing
  W.numerator = matrix(apply(mu.adj.comp, 1, function(x) combn(x, 2, FUN = diff)), ncol = n.taxa)
  W.numerator = t(W.numerator)
  # Variance of estimated mean difference
  W.denominator1 = matrix(apply(mu.var, 1, function(x) combn(x, 2, FUN = sum)), ncol = n.taxa)
  # Variance of delta_hat
  if (length(bias.var.vec) < 2) {
    W.denominator2 = bias.var.vec
  }else {
    W.denominator2 = c(bias.var.vec, combn(bias.var.vec, 2, FUN = sum))
  }
  W.denominator = W.denominator1 + W.denominator2 + 2 * sqrt(W.denominator1 * W.denominator2)
  W.denominator = t(sqrt(W.denominator))
  grp.pair = combn(n.grp, 2)
  colnames(W.numerator) = sapply(1:ncol(grp.pair), function(x) 
    paste0("mean.difference (", grp.name[grp.pair[2, x]], " - ", grp.name[grp.pair[1, x]], ")"))
  colnames(W.denominator) = sapply(1:ncol(grp.pair), function(x) 
    paste0("se (", grp.name[grp.pair[2, x]], " - ", grp.name[grp.pair[1, x]], ")"))
  rownames(W.numerator) = taxa.id; rownames(W.denominator) = taxa.id
  
  if (length(grp.name) == 2) {
    ## Two-group comparison
    W = W.numerator/W.denominator
    p.val = sapply(W, function(x) 2*pnorm(abs(x), mean = 0, sd = 1, lower.tail = F))
    q.val = p.adjust(p.val, method = adj.method)
    q.val[is.na(q.val)] = 1
  } else {
    ## Multi-group comparison: Permutation test
    # Test statistics
    W.each = W.numerator/W.denominator
    W.each[is.na(W.each)] = 0 # Replace missing values with 0s
    W = apply(abs(W.each), 1, max)
    
    # Test statistics under null
    W.null.list = lapply(1:perNum, function(x) {
      set.seed(x)
      mu.adj.comp.null = matrix(rnorm(n.taxa * n.grp), nrow = n.taxa, ncol = n.grp) * sqrt(mu.var)
      W.numerator.null = matrix(apply(mu.adj.comp.null, 1, function(x) combn(x, 2, FUN = diff)), ncol = n.taxa)
      W.numerator.null = t(W.numerator.null)
      
      W.each.null = W.numerator.null/W.denominator
      W.each.null[is.na(W.each.null)] = 0
      W.null = apply(abs(W.each.null), 1, max)
      return(W.null)
    })
    W.null = Reduce('cbind', W.null.list)
    
    # Test results
    p.val = apply(W.null - W, 1, function(x) sum(x > 0)/perNum)
    q.val = p.adjust(p.val, method = adj.method)
    q.val[is.na(q.val)] = 1
  }
  
  W = matrix(W, ncol = 1)
  colnames(W) = "W"
  res.comp = data.frame(W.numerator, W.denominator, W = W, p.val, q.val, check.names = FALSE)
  
  ### 5. Combine results from structural zeros
  mu.adj = matrix(NA, nrow = n.taxa.raw, ncol = n.grp)
  colnames(mu.adj) = grp.name; rownames(mu.adj) = taxa.id.raw
  mu.adj[comp.taxa.pos, ] = mu.adj.comp
  
  if (length(comp.taxa.pos) < n.taxa.raw) {
    O.incomp = feature.table[-comp.taxa.pos, ]
    ind.incomp = struc.zero[-comp.taxa.pos, rep(1:n.grp, times = n.samp.grp)]
    y.incomp = log(O.incomp + 1)
    d.incomp = t(t(1 - ind.incomp) * d) # Sampling fractions for entries considered to be structural zeros are set to be 0s
    y.adj.incomp = y.incomp - d.incomp
    mu.incomp = t(apply(y.adj.incomp, 1, function(i) 
      tapply(i, rep(1:n.grp, n.samp.grp), function(j) mean(j, na.rm = T))))
    # In case of negative values for mean absolute abundances
    mu.adj.incomp = mu.incomp
    mu.adj.incomp[mu.adj.incomp == 0] = NA
    mu.adj.incomp = t(t(mu.adj.incomp) + abs(apply(mu.incomp, 2, min)))
    mu.adj.incomp[is.na(mu.adj.incomp)] = 0
  }else{
    mu.adj.incomp = NA
  }
  mu.adj[-comp.taxa.pos, ] = mu.adj.incomp
  colnames(mu.adj) = paste0("mean.absolute.abundance (", grp.name, ")")
  rownames(mu.adj) = taxa.id.raw
  
  ### 6. Outputs
  W.numerator = matrix(apply(mu.adj, 1, function(x) combn(x, 2, FUN = diff)), ncol = n.taxa.raw)
  W.numerator = t(W.numerator)
  W.denominator = matrix(0, ncol = ncol(W.numerator), nrow = nrow(W.numerator))
  
  res = data.frame(W.numerator, W.denominator, W = Inf, p.val = 0, q.val = 0, check.names = FALSE)
  res[comp.taxa.pos, ] = res.comp
  colnames(res) = colnames(res.comp); rownames(res) = taxa.id.raw
  res = res%>%mutate(diff.abn = ifelse(q.val < alpha, TRUE, FALSE))
  
  out = list(feature.table = feature.table, res = res, d = d.adj, mu = mu.adj, bias.em = bias.em.vec, bias.wls = bias.wls.vec)
  return(out)
}
reference.select <- function(data, Y, pseudo.count = 1){
  raw.data <- data
  data[which(data==0)] <- pseudo.count
  n.tax <- ncol(data)
  p.r <- matrix(0,n.tax,n.tax)
  q.r <- matrix(0,n.tax,n.tax)
  for(i in 1:n.tax){
    for(j in i:n.tax){
      dcd <- log(data[,i]/data[,j])
      lmr <- lm(dcd~Y)
      lmrs <- summary(lmr)
      pv1 <- lmrs$coefficients[2,4]
      if(!is.na(pv1)){
        p.r[i,j] <- pv1
        p.r[j,i] <- pv1
      }
    }
  }
  p.adj <- function(x){p.adjust(x, method = "BH")}
  q.r <- apply(p.r, 1, p.adj)
  q.r.med <- apply(q.r, 1, median)
  rankq <- rank(q.r.med)
  z.num <- rep(1,n.tax)
  for(j in 1:n.tax){z.num[j] <- sum(raw.data[,j]!=0)}
  rankp <- rank(z.num)
  rankj <- rankq+rankp
  target.cut <- c(0.2,0.1,0.05,0.01)
  ttt <- 1
  exclude <- which(q.r.med<=target.cut[1])
  while(length(exclude) == n.tax){
    ttt <- ttt+1
    exclude <- which(q.r.med<=target.cut[ttt])
    if(ttt==length(target.cut)){
      break
    }
  }
  if(length(exclude) >= round(n.tax*0.8)){
    ref <- order(rankj, decreasing = TRUE)[1:round(n.tax*0.2)]
    n.ref <- length(ref)
  }else{
    if(length(exclude)/n.tax <0.5){
      cut <- round(n.tax*0.5)
      ref.candi <- setdiff(order(rankj, decreasing = TRUE),exclude)
      ref <- ref.candi[1:cut]
    }else{
      ref <- setdiff(order(rankj, decreasing = TRUE),exclude)
    }
    n.ref <- length(ref)
  }
  return(ref)
}
bonobo.precision <- function(DATA,rho){
  n.sample <- ncol(DATA)
  g <- nrow(DATA)
  DATA[which(DATA == 0)] <- 1
  x <- log(DATA)
  x.mean <- rowMeans(x)
  X <- x-x.mean
  leaveoneCOV <- function(i){
    data <- t(X[,-i])
    return(cov(data))
  }
  S <- lapply(1:n.sample, leaveoneCOV)
  diagS <- matrix(0, nrow = g, ncol = n.sample)
  for(i in 1:n.sample){
    diagS[,i] <- diag(S[[i]])
  }
  eta <- rep(0, g)
  for(j in 1:g){
    eta[j] <- var(diagS[j,])*(n.sample-1)/n.sample
  }
  eta.sum <- sum(eta)
  sk <- rep(0, n.sample)
  for(i in 1:n.sample){
    sk[i] <- 2*sum(diagS[,i]^2)
  }
  vi_g_3 <- sk/eta.sum
  vi_g <- vi_g_3+3
  delta <- 1/vi_g
  
  sample.specific.cov <- function(i){
    Vi <- delta[i]*X[,i]%*%t(X[,i])+(1-delta[i])*S[[i]]
    diag(Vi)[which(diag(Vi)==0)] <- 1
    return(Vi)
  }
  V <- lapply(1:n.sample, sample.specific.cov)
  
  sample.specific.precision <- function(i){
    precision <- glasso(V[[i]], rho=rho)[["wi"]]
    return(precision)
  }
  Precision <- lapply(1:n.sample, sample.specific.precision)
  
  return(Precision)
}
ADMIC <- function(data, Y, C = NULL, 
                  pseudo.count = 1, c = 0.1,
                  r.test = seq(from = 0.25, to = 0.75, by = 0.25)){
  
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
  Precision.mat <- bonobo.precision(t(data),rho = rho)
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


gLV_sim <- function(data){
  library(dirmult)
  library(permute)
  library(lattice)
  library(vegan)
  library(Rcpp)
  library(RcppArmadillo)
  library(utils)
  library(metap)
  library(psych)
  library(parallel)
  library(LOCOM)
  library(ROCR)
  library(carData)
  library(car)
  library(plyr)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(MASS)
  library(survival)
  library(NADA)
  library(truncnorm)
  library(zCompositions)
  library(nloptr)
  library(nlme)
  library(ALDEx2)
  library(Maaslin2)
  library(GUniFrac)
  library(stats)
  library(statmod)
  library(ggrepel)
  library(foreach)
  library(modeest)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(Matrix)
  library(magrittr)
  library(glmnet)
  library(glasso)
  library(scalreg)
  library(lars)
  library(hdi)
  library(psych)
  library(WGCNA)
  library(deSolve)
  library(igraph)
  library(phyloseq)
  library(ADAPT)
  methodsnum <- 10
  SEN <- matrix(0,nrow=methodsnum,ncol=4)
  FDR <- matrix(0,nrow=methodsnum,ncol=4)
  AUPR <- matrix(0,nrow=methodsnum,ncol=4)
  controlFDR <- matrix(0,nrow=methodsnum*4,ncol=4)
  for(d in 1:4){
    # Setting network parameters
    n_nodes <- 10*d    # number of taxa
    prob <- 0.25       # connectivity parameter
    delta <- 0.2
    n0.sam <- 50
    # metadata
    Y <- c(rep(0,n0.sam),rep(1,n0.sam))
    control <- which(Y==0)
    case <- which(Y==1)
    n.sam <- length(Y)
    sam.ID <- paste0("Sample", 1:n.sam)
    meta.data <- data.frame(group = Y, row.names = sam.ID ,stringsAsFactors = FALSE)
    colnames(meta.data) <- "group"
    n.sam.control <- length(control)
    n.sam.case <- length(case)
    
    # ER random graph
    network <- sample_gnp(n = n_nodes, p = prob, directed = T)
    # Compute the adjacency matrix
    adj_matrix <- as_adjacency_matrix(network, sparse = F)
    nonzero.site <- which(adj_matrix == 1)
    N <- length(nonzero.site)
    a0 <- rnorm(N, mean = 0, sd = delta)
    A <- adj_matrix
    A[nonzero.site] <- a0
    diag(A) <- -2.5*delta
    
    # gLV Equation
    gLV_model <- function(t, state, parameters) {
      N <- state
      r <- parameters$r
      alpha <- parameters$alpha
      dNdt <- numeric(length(N))
      for (i in 1:length(N)) {
        dNdt[i] <- N[i] * (r[i] + sum(alpha[i, ] * N))
      }
      return(list(dNdt))
    }
    state <- runif(n_nodes)
    # Setting model parameters
    parameters <- list(
      r = runif(n_nodes, min = 0.8, max = 1.3),  # growth rate
      alpha = A  # interaction matrix
    )
    # time vector
    time <- seq(0, 50, by = 1)
    
    # Setting DA taxa
    r0 <- parameters[["r"]]
    r1 <- r0
    causal.taxa <- sample(c(1:n_nodes),round(n_nodes*0.2))
    dt0 <- rep(0, n_nodes)
    dt0[causal.taxa] <- 1
    noncausal.taxa <- setdiff(c(1:n_nodes),causal.taxa)
    beta <- runif(round(n_nodes*0.2), min = -0.3, max = 0.3)
    r1[causal.taxa] <- r1[causal.taxa]+beta
    
    # Generating data
    R0 <- matrix(0, nrow = n0.sam, ncol = n_nodes)
    R1 <- matrix(0, nrow = n0.sam, ncol = n_nodes)
    for(i in 1:n_nodes){
      R0[,i] <- rnorm(n0.sam, mean = r0[i], sd = 0.05)
      R1[,i] <- rnorm(n0.sam, mean = r1[i], sd = 0.05)
    }
    table <- matrix(0, nrow = 2*n0.sam, ncol = n_nodes)
    for(i in 1:n0.sam){
      state <- runif(n_nodes)
      simparameters <- list(
        r =  R0[i,],  
        alpha = A 
      )
      simout <- ode(y = state, times = time, func = gLV_model, parms = simparameters)
      table[i,] <- simout[51,-1]
    }
    for(i in 1:n0.sam){
      state <- runif(n_nodes)
      simparameters <- list(
        r =  R1[i,],  
        alpha = A  
      )
      simout <- ode(y = state, times = time, func = gLV_model, parms = simparameters)
      table[i+n0.sam,] <- simout[51,-1]
    }
    
    table[which(table < 1e-6)] <- 1e-6
    
    table.raw <- add_poisson_noise(table)
    
    # -----------------------------------
    #  Filtering otus
    # -----------------------------------
    prop.presence0 <- colMeans(table.raw[control,] > 0)
    prop.presence1 <- colMeans(table.raw[case,] > 0)
    otus.keep0 <- which(prop.presence0 >= 0.05)
    otus.keep1 <- which(prop.presence1 >= 0.05)
    otus.keep <- intersect(otus.keep0, otus.keep1)
    otu.table.filter <- table.raw[, otus.keep]
    causal.otus.filter <- which(otus.keep %in% causal.taxa)
    noncausal.otus.filter <- setdiff(1:length(otus.keep), causal.otus.filter)
    dt <- dt0[otus.keep]
    n.tax <- length(otus.keep)
    n.sam <- nrow(table.raw)
    
    ###########################
    # DAA methods analysis
    ###########################
    data <- otu.table.filter
    data[which(data==0)] <- 1
    taxaname <- paste0("Taxa", 1:n.tax)
    colnames(data) <- taxaname
    
    #######################
    # ADMIC & ADMIC-Basic
    ######################
    res.ADMIC <- ADMIC(otu.table.filter,Y)
    q0 <- res.ADMIC$q.basic
    fq0 <- 1-q0
    aupr.cor0 <- performance(prediction(fq0,dt),"aucpr")
    AUPR[9,d] <- unlist(slot(aupr.cor0,"y.values"))
    
    q.eff <- res.ADMIC$q.admic
    fq.eff <- 1-q.eff
    aupr.cor.eff <- performance(prediction(fq.eff,dt),"aucpr")
    AUPR[10,d] <- unlist(slot(aupr.cor.eff,"y.values"))
  
    
    ###############
    # Wilcoxn-clr
    ###############
    clr.data <- data
    geo.mean.clr <- apply(clr.data, 1, geometric.mean)
    for(j in 1:n.sam){
      clr.data[j,] <- clr.data[j,]/geo.mean.clr[j]
    }
    clr.data <- log(clr.data)
    
    p.clr <- rep(1,n.tax)
    for(j in 1:n.tax){
      p.clr[j] <- wilcox.test(clr.data[control,j], clr.data[case,j])$p.value
    }
    
    q.clr <- p.adjust(p.clr, method = "BH")
    no.clr <- which(q.clr>=0)
    fq.clr <- 1-q.clr[no.clr]
    dt.clr <- dt[no.clr]
    pred.clr <- performance(prediction(fq.clr,dt.clr),"aucpr")
    AUPR[1,d] <- unlist(slot(pred.clr,"y.values"))
    
    ###########
    # MaAsLin2
    ###########
    OTU <- data.frame(data)
    name <- rep("N",n.tax)
    name[causal.otus.filter] <- 'Y'
    colnames(OTU) <- name
    rownames(OTU) <- sam.ID
    
    res <- Maaslin2(
      input_data = OTU,
      input_metadata = meta.data,
      output = "demo_output",
      fixed_effects = 'group',
      plot_heatmap =FALSE,
      plot_scatter = FALSE)
    
    p.maaslin2 <- res$results[,c("feature","pval")]$pval
    q.maaslin2 <- res$results[,c("feature","qval")]$qval
    dt.maaslin2 <- res$results[,c("feature","qval")]$feature
    lnn <- length(dt.maaslin2)
    for(l in 1:lnn){
      if(grepl('Y',dt.maaslin2[l])==TRUE){
        dt.maaslin2[l] <- "Y"
      }else{dt.maaslin2[l] <- "N"}
      c.otu <- c(which(dt.maaslin2=="Y"))
      n.otu <- setdiff(1:lnn, c.otu)
    }
    dt.ma <- rep(0,length(q.maaslin2))
    dt.ma[c.otu] <- 1
    fq.maaslin2 <- 1-q.maaslin2
    aupr.maaslin2 <- performance(prediction(fq.maaslin2,dt.ma),"aucpr")
    AUPR[2,d] <- unlist(slot(aupr.maaslin2,"y.values"))
    
    ###########
    # ALDEx2
    ###########
    C <- NULL
    otu.table.aldex2 <- t(data)
    row.names(otu.table.aldex2) <- taxaname
    
    aldex.sample <- aldex.clr(otu.table.aldex2, as.character(Y), mc.samples=128, verbose=TRUE)
    res.aldex2 <- aldex.ttest(aldex.sample, paired.test=FALSE, hist.plot = FALSE)
    lll <- row.names(res.aldex2)
    q.aldex2 <- rep(1,n.tax)
    keeptaxa <- c(which(taxaname %in% lll))
    q.aldex2[keeptaxa] <- res.aldex2$wi.eBH
    fq.aldex2 <- 1-q.aldex2
    aupr.aldex2 <- performance(prediction(fq.aldex2,dt),"aucpr")
    AUPR[3,d] <- unlist(slot(aupr.aldex2,"y.values"))
    
    
    ###########
    # LinDA
    ###########
    data.linda<- t(data)
    
    res.linda <- linda(data.linda,meta.data,formula="~group",alpha=0.05)
    
    p.linda <- res.linda[["output"]][["group"]][["pval"]]
    q.linda <- res.linda[["output"]][["group"]][["padj"]]
    fq.linda <- 1-q.linda
    pred.linda <- performance(prediction(fq.linda,dt),"aucpr")
    AUPR[4,d] <- unlist(slot(pred.linda,"y.values"))
    
    ###########
    # ANCOMBC
    ###########
    otu.table.ancombc <- t(data)
    colnames(otu.table.ancombc) <- sam.ID
    meta.data <- data.frame(group = Y, row.names = sam.ID, stringsAsFactors = FALSE)
    formula <- "group"
    group <- factor(meta.data[, formula])
    group.name <- levels(group)
    grp.ind <- lapply(1:nlevels(group), function(i) which(group == group.name[i]))
    struc.zero <- matrix(0, nrow = n.tax, ncol = n.sam)
    res.ancom_bc <- ANCOM_BC(otu.table.ancombc,group.name ,grp.ind ,struc.zero)
    
    p.ancombc <- res.ancom_bc[["res"]][["p.val"]]
    q.ancombc <- res.ancom_bc[["res"]][["q.val"]]
    fq.ancombc <- 1-q.ancombc
    aupr.ancombc <- performance(prediction(fq.ancombc,dt),"aucpr")
    AUPR[5,d] <- unlist(slot(aupr.ancombc,"y.values"))
   
    ###########
    # LOCOM
    ###########
    res.locom <- locom(otu.table = data, Y = Y, C = NULL)
    p.locom <- res.locom$p.otu
    q.locom <- res.locom$q.otu
    fq.locom <- 1-q.locom
    aupr.locom <- performance(prediction(fq.locom[1,],dt),"aucpr")
    AUPR[6,d] <- unlist(slot(aupr.locom,"y.values"))
   
    ###########
    # ZicoSeq
    ###########
    taxaname <- colnames(data)
    comm <- t(otu.table.filter)
    row.names(comm) <- taxaname
    
    zico.obj <- ZicoSeq(meta.dat = meta.data, feature.dat = comm, 
                        grp.name = 'group', adj.name =NULL, feature.dat.type = "other",
                        # Filter to remove rare taxa
                        prev.filter = 0.1, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
                        # Winsorization to replace outliers
                        is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                        # Posterior sampling to impute zeros
                        is.post.sample = TRUE, post.sample.no = 25, 
                        # Multiple link functions to capture diverse taxon-covariate relation
                        link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75), 
                        stats.combine.func = max,
                        # Permutation-based multiple testing correction
                        perm.no = 99,  strata = NULL, 
                        # Reference-based multiple stage normalization
                        ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                        # Family-wise error rate control
                        is.fwer = FALSE,
                        verbose = TRUE, return.feature.dat = FALSE)
    p.zicoseq <- as.data.frame(zico.obj$p.raw)
    p.zi <- p.zicoseq$`zico.obj$p.raw`
    q.zicoseq <- as.data.frame(zico.obj$p.adj.fdr)
    q.zi <- q.zicoseq$`zico.obj$p.adj.fdr`
    lll <- row.names(q.zicoseq)
    q.zicoseq <- rep(1,n.tax)
    keeptaxa <- c(which(taxaname %in% lll))
    q.zicoseq[keeptaxa] <- q.zi
    fq.zicoseq <- 1-q.zicoseq
    pred.zicoseq <- performance(prediction(fq.zicoseq,dt),"aucpr")
    AUPR[7,d] <- unlist(slot(pred.zicoseq,"y.values"))
    
    ###########
    # ADAPT
    ##########
    data.adapt <- t(matrix(as.integer(t(data)),n.tax,n.sam))
    taxaname <-  paste0("taxon", 1:n.tax)
    colnames(data.adapt) <- taxaname
    Ya <- as.character(Y)
    Y.adapt <- as.data.frame(Ya)
    phy.data <- phyloseq(otu_table(data.adapt,taxa_are_rows = F),sample_data(Y.adapt))
    res.adapt <- adapt(phy.data, cond.var = "Ya")
    remaintaxa <- which(taxaname %in% res.adapt@details[["Taxa"]])
    q.adapt <- rep(1,n.tax)
    q.adapt[remaintaxa] <- res.adapt@details[["adjusted_pval"]]
    fq.adapt <- 1-q.adapt
    pred.adapt <- performance(prediction(fq.adapt,dt),"aucpr")
    AUPR[8,d] <- unlist(slot(pred.adapt,"y.values"))
    
    otu.clr <- sumresult(q.clr,causal.otus.filter,noncausal.otus.filter)
    SEN[1,d] <- otu.clr$sen
    FDR[1,d] <- otu.clr$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[1+kk,d] <- controlRecall(q.clr,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.maaslin2 <- sumresult(q.maaslin2,c.otu,n.otu)
    SEN[2,d] <- otu.maaslin2$sen
    FDR[2,d] <- otu.maaslin2$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[2+kk,d] <- controlRecall(q.maaslin2,c.otu,n.otu,target = (0.2*i))
    }
    
    otu.aldex2 <- sumresult(q.aldex2,causal.otus.filter,noncausal.otus.filter)
    SEN[3,d] <- otu.aldex2$sen
    FDR[3,d] <- otu.aldex2$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[3+kk,d] <- controlRecall(q.aldex2,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.linda=sumresult(q.linda,causal.otus.filter,noncausal.otus.filter)
    SEN[4,d] <- otu.linda$sen
    FDR[4,d] <- otu.linda$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[4+kk,d] <- controlRecall(q.linda,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.ancombc <- sumresult(q.ancombc,causal.otus.filter,noncausal.otus.filter)
    SEN[5,d] <- otu.ancombc$sen
    FDR[5,d] <- otu.ancombc$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[5+kk,d] <- controlRecall(q.ancombc,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.locom <- sumresult(q.locom,causal.otus.filter,noncausal.otus.filter)
    SEN[6,d] <- otu.locom$sen
    FDR[6,d] <- otu.locom$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[6+kk,d] <- controlRecall(q.locom,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.zicoseq <- sumresult(q.zicoseq,causal.otus.filter,noncausal.otus.filter)
    SEN[7,d] <- otu.zicoseq$sen
    FDR[7,d] <- otu.zicoseq$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[7+kk,d] <- controlRecall(q.zicoseq,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.adapt <- sumresult(q.adapt,causal.otus.filter,noncausal.otus.filter)
    SEN[8,d] <- otu.adapt$sen
    FDR[8,d] <- otu.adapt$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[8+kk,d] <- controlRecall(q.adapt,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.lm <- sumresult(q0,causal.otus.filter,noncausal.otus.filter)
    SEN[9,d] <- otu.lm$sen
    FDR[9,d] <- otu.lm$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[9+kk,d] <- controlRecall(q0,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
    otu.admic.eff <- sumresult(q.eff,causal.otus.filter,noncausal.otus.filter)
    SEN[10,d] <- otu.admic.eff$sen
    FDR[10,d] <- otu.admic.eff$fdr
    for(i in 1:4){
      kk <- (i-1)*methodsnum
      controlFDR[10+kk,d] <- controlRecall(q.eff,causal.otus.filter,noncausal.otus.filter,target = (0.2*i))
    }
    
   #-----------------------------------------------
    
  }
  
  myres <- list(SEN,FDR,AUPR,controlFDR)
  return(myres)
}


#----Multicore Operation-------
closeAllConnections()

cores <- 30
cl <- makeCluster(cores)
registerDoParallel(cl)

start <- Sys.time()
myresults <- foreach(i = 1:30, .errorhandling = "pass") %dopar% {
  gLV_sim(i)
}
end <- Sys.time()
time <- difftime(end,start,units = "min")
print(time)

stopCluster(cl)

sen <- 0
fdr <- 0
aupr <- 0
cfdr <- 0
n <- 0
for(aa in 1:30){
  if(class(myresults[[aa]][[1]])[1] == "matrix"){
    sen <- myresults[[aa]][[1]]+sen
    fdr <- myresults[[aa]][[2]]+fdr
    aupr <- myresults[[aa]][[3]]+aupr
    cfdr <- myresults[[aa]][[4]]+cfdr
    n <- n+1
  }
}
print(n)
SEN <- sen/n
FDR <- fdr/n
AUPR <- aupr/n
CFDR <- cfdr/n
controlFDR <- matrix(rowMeans(CFDR), nrow = 10, ncol = 4)
methodname <- c("wilcox-clr","MaAsLin2", "Aldex2","Linda","ANCOMBC","LOCOM","ZicoSeq","ADAPT","ADMIC-Basic","ADMIC")
rownames(AUPR) <- methodname
rownames(controlFDR) <- methodname
rownames(FDR) <- methodname
rownames(SEN) <- methodname

