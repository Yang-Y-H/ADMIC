library(deSolve)
library(ggplot2)
library(reshape2)
#------Figure 1B-----------------------
# Setting network parameters
n_nodes <- 4    
adj_matrix <- matrix(0,4,4)
adj_matrix[1,3] <- 1
adj_matrix[2,c(1,3)] <- 1
adj_matrix[3,c(1,4)] <- 1
A <- adj_matrix
A[which(A==1)] <- 0.2
diag(A) <- -0.4

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

# Setting model parameters
parameters <- list(
  r = rep(1, n_nodes),  # growth rate
  alpha = A  # interaction matrix
)
# Setting DA taxa
r0 <- parameters[["r"]]
r1 <- r0
r1[1] <- r1[1]-0.3
state <- rep(5, n_nodes)
time <- seq(0, 50, by = 1)
out <- ode(y = state, times = time, func = gLV_model, parms = parameters)
state1 <- out[51,-1]
parameters1 <- list(
  r =  r1,  # growth rate
  alpha = A  # interaction matrix
)
out1 <- ode(y = state1, times = time, func = gLV_model, parms = parameters1)

#----Absolute abundance-------
myout <- rbind(out[-51,],out1[-51,])[,-1]

data.single <- data.frame(
  x = rep(c(1:100), 4),
  y = as.vector(myout),
  group = rep(c("Taxon1","Taxon2","Taxon3","Taxon4"), each=100)
)
colors <-c("#1F78B4","#33A02C","#E31A1C" ,"#FF7F00")

ggplot(data.single, aes(x=x, y=y, color = factor(group), group=group)) +
  geom_line(linewidth=0.8) +
  ylim(2, 9) +
  labs(x="Time", y="Population") +
  geom_vline(xintercept = 51, color = "black", linetype = "dashed")+
  scale_color_manual(values=colors) +
  theme_bw() +
  theme(legend.title=element_blank())

#-----Relative abundance---------
ra.myout <- myout/rowSums(myout)

ra.data.single <- data.frame(
  x = rep(c(1:100), 4),
  y = as.vector(ra.myout),
  group = rep(c("Taxon1","Taxon2","Taxon3","Taxon4"), each=100)
)

ggplot(ra.data.single, aes(x=x, y=y, color = factor(group), group=group)) +
  geom_line(linewidth=0.8) +
  ylim(0.1, 0.4) +
  labs(x="Time", y="Relative Abundance") +
  geom_vline(xintercept = 51, color = "black", linetype = "dashed")+
  scale_color_manual(values=colors) +
  theme_bw() +
  theme(legend.title=element_blank())


#------Figure 1C-----------------------
#--load function---
bonobo.precision <- function(DATA,rho){
  ########################
  #DATA：abundance matrix,rows are taxa, columns are samples
  # Data Preprocessing
  n.sample <- ncol(DATA)
  g <- nrow(DATA)
  # Add pseudo count
  y <- DATA[1,]
  DATA <- DATA[-1,]
  DATA[which(DATA==0)] <- 1
  x.taxa <- DATA%>% log()
  x <- rbind(y, x.taxa)
  x.mean <- x %>% rowMeans()
  X <- x-x.mean
  # Compute the covariance matrix
  leaveoneCOV <- function(i){
    data <- X[,-i] %>% t()
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
    sk[i] <- 2* sum(diagS[,i]^2)
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
    precision <- glasso(V[[i]],rho=rho)[["wi"]]
    return(precision)
  }
  Precision <- lapply(1:n.sample, sample.specific.precision)
  
  return(Precision)
}
sumresult <- function(qvalue, causal.otus, noncausal.otus, fdr.target=0.05) {
  
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
  ord=order(qvalue)[1:2]
  out = c(sen=sen,fdr=fdr)
  return(out)
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
sensitivity.top <- function(top.order,target.order){
  n <- length(target.order)
  sen.top <- length(intersect(top.order[1:n], target.order))/n
  return(sen.top)
}

mysim <- function(data){
  methodsnum <- 8
  AUPR <- rep(0, methodsnum)
  Q.value <- matrix(0, nrow=4, ncol=methodsnum)
  sam.number <- 100
  n0.sam <- 50
  # Setting metadata
  Y <- c(rep(0,n0.sam),rep(1,n0.sam))
  control <- which(Y==0)
  case <- which(Y==1)
  n.sam <- length(Y)
  sam.ID <- paste0("Sample", 1:n.sam)
  meta.data <- data.frame(group = Y, row.names = sam.ID ,stringsAsFactors = FALSE)
  colnames(meta.data) <- "group"
  n.sam.control <- length(control)
  n.sam.case <- length(case)
  
  # Setting network parameters
  n_nodes <- 4  
  adj_matrix <- matrix(0,4,4)
  adj_matrix[1,3] <- 1
  adj_matrix[2,c(1,3)] <- 1
  adj_matrix[3,c(1,4)] <- 1
  A <- adj_matrix
  A[which(A==1)] <- 0.2
  diag(A) <- -0.4
  
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
  
  parameters <- list(
    r = runif(n_nodes, min = 1, max = 1),  # growth rate
    alpha = A  # interaction matrix
  )
  time <- seq(0, 50, by = 1)
  
  # Setting DA taxa
  r0 <- parameters[["r"]]
  r1 <- r0
  causal.taxa <- 1
  dt0 <- rep(0, n_nodes)
  dt0[causal.taxa] <- 1
  noncausal.taxa <- setdiff(c(1:n_nodes),causal.taxa)
  r1[causal.taxa] <- r1[causal.taxa]-0.3
  
  # Generating data
  R0 <- matrix(0, nrow = n0.sam, ncol = n_nodes)
  R1 <- matrix(0, nrow = n0.sam, ncol = n_nodes)
  for(i in 1:n_nodes){
    R0[,i] <- rnorm(n0.sam, mean = r0[i], sd = 0.1)
    R1[,i] <- rnorm(n0.sam, mean = r1[i], sd = 0.1)
  }
  table <- matrix(0, nrow = 2*n0.sam, ncol = n_nodes)
  for(i in 1:n0.sam){
    state <- rep(1, n_nodes)
    simparameters <- list(
      r =  R0[i,],  
      alpha = A  
    )
    simout <- ode(y = state, times = time, func = gLV_model, parms = simparameters)
    table[i,] <- simout[51,-1]
  }
  for(i in 1:n0.sam){
    state <- table[i,]
    simparameters <- list(
      r =  R1[i,],  
      alpha = A 
    )
    simout <- ode(y = state, times = time, func = gLV_model, parms = simparameters)
    table[i+n0.sam,] <- simout[51,-1]
  }
  
  table[which(table < 1e-6)] <- 1e-6
  
  table.raw <- add_poisson_noise(table)
  
  otus.keep <- c(1:n_nodes)
  otu.table.filter <- table.raw
  causal.otus.filter <- causal.taxa
  noncausal.otus.filter <- setdiff(1:length(otus.keep), causal.otus.filter)
  dt <- dt0[otus.keep]
  n.tax <- length(otus.keep)
  n.sam <- nrow(table.raw)
  
  data <- otu.table.filter
  data[which(data==0)] <- 1
  taxaname <- paste0("Taxa", 1:n.tax)
  colnames(data) <- taxaname
  
  ###########
  # Wilcoxn-clr
  ###########
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
  AUPR[1] <- unlist(slot(pred.clr,"y.values"))
  Q.value[,1] <- q.clr
  
  ###########
  # MaAsLin2
  ###########
  OTU <- data.frame(data)
  colnames(OTU) <- taxaname
  rownames(OTU) <- sam.ID
  
  res <- Maaslin2(
    input_data = OTU,
    input_metadata = meta.data,
    output = "demo_output",
    fixed_effects = 'group',
    plot_heatmap =FALSE,
    plot_scatter = FALSE)
  
  p.maas <- res$results[,c("feature","pval")]$pval
  q.maas <- res$results[,c("feature","qval")]$qval
  maas.order <- match(taxaname, res$results[,c("feature","qval")]$feature)
  p.maaslin2 <- p.maas[maas.order]
  q.maaslin2 <- q.maas[maas.order]
  fq.maaslin2 <- 1-q.maaslin2
  aupr.maaslin2 <- performance(prediction(fq.maaslin2,dt),"aucpr")
  AUPR[2] <- unlist(slot(aupr.maaslin2,"y.values"))
  Q.value[,2] <- q.maaslin2 
  
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
  AUPR[3] <- unlist(slot(aupr.aldex2,"y.values"))
  Q.value[,3] <- q.aldex2
  
  ###########
  # LinDA
  ###########
  data.linda<- t(data)
  
  res.linda <- linda(data.linda,meta.data,formula="~group",alpha=0.05)
  
  p.linda <- res.linda[["output"]][["group"]][["pval"]]
  q.linda <- res.linda[["output"]][["group"]][["padj"]]
  fq.linda <- 1-q.linda
  pred.linda <- performance(prediction(fq.linda,dt),"aucpr")
  AUPR[4] <- unlist(slot(pred.linda,"y.values"))
  Q.value[,4] <- q.linda
  
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
  AUPR[5] <- unlist(slot(aupr.ancombc,"y.values"))
  Q.value[,5] <- q.ancombc
  
  ###########
  # LOCOM
  ###########
  res.locom <- locom(otu.table = data, Y = Y, C = NULL)
  p.locom <- res.locom$p.otu
  q.locom <- res.locom$q.otu
  fq.locom <- 1-q.locom
  aupr.locom <- performance(prediction(fq.locom[1,],dt),"aucpr")
  AUPR[6] <- unlist(slot(aupr.locom,"y.values"))
  Q.value[,6] <- q.locom
  
  ###########
  # ZicoSeq
  ###########
  taxaname <- colnames(data)
  comm <- t(otu.table.filter)
  row.names(comm) <- taxaname
  meta.dat <- meta.data
  
  zico.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
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
  AUPR[7] <- unlist(slot(pred.zicoseq,"y.values"))
  Q.value[,7] <- q.zicoseq
  
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
  AUPR[8] <- unlist(slot(pred.adapt,"y.values"))
  Q.value[,8] <- q.adapt
  
  colnames(Q.value) <- c("Wilcoxon-CLR","MaAsLin2", "Aldex2","LinDA","ANCOMBC","LOCOM","ZicoSeq","ADAPT")
  FDR <- sapply(c(1:8), function(p) {
    sumresult(qvalue = Q.value[,p], causal.otus = causal.taxa, noncausal.otus = noncausal.taxa)
  })
  
  RANK <- apply(Q.value,2,rank)
  top.taxa <- apply(RANK,2, order)
  SEN <- apply(top.taxa, 2, sensitivity.top, target.order = causal.taxa)
  res <- list(FDR = FDR[2,], SEN)
  return(res)
}

#----Multicore Operation-------
library(doParallel)
cores <- 30
cl <- makeCluster(cores)
registerDoParallel(cl)

start=Sys.time()
AAmyresults <- foreach(i = 1:100, .errorhandling = "pass") %dopar% {
  mysim(i)
}
end=Sys.time()
time=difftime(end,start,units = "min")
print(time)

stopCluster(cl)


mysen=0
myfdr=0
n=0
for(aa in 1:100){
  if(class(AAmyresults[[aa]][[1]])[1]=="numeric"){
    mysen=AAmyresults[[aa]][[2]]+mysen
    myfdr=AAmyresults[[aa]][[1]]+myfdr
    n=n+1
  }
}
print(n)
SEN0=mysen/n
FDR0=myfdr/n


df1 <- data.frame(
  Method = factor(names(SEN0),
                  levels = names(SEN0)[order(FDR0-SEN0)]),
  FDR = FDR0,
  SEN = SEN0
)
ggplot(df1, aes(x = Method, group = 1)) +
  geom_line(aes(y = FDR, color = "FDR"), linewidth = 1) +
  geom_point(aes(y = FDR, color = "FDR", shape = "FDR"), size = 2) +
  geom_line(aes(y = SEN, color = "SEN"), linewidth = 1) +
  geom_point(aes(y = SEN, color = "SEN", shape = "SEN"), size = 2) +
  scale_y_continuous(
    name = "FDR",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    sec.axis = sec_axis(~ ., name = "Sensitivity of Top")
  ) +
  scale_color_manual(values = c("FDR" = "#D62728", "SEN" = "#1F78B4")) +
  scale_shape_manual(values = c("FDR" = 16, "SEN" = 17)) +
  labs(x = NULL, title = "gLV", color = NULL, shape = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    legend.position = "bottom"
  )





#------Figure 1D-----------------------
QVAL1 <- read.csv("data/T4+F1_q.csv",row.names = 1)
QVAL2 <- read.csv("data/B40-8+VD13_q.csv",row.names = 1)
sumresult <- function(qvalue, causal.otus, noncausal.otus, fdr.target=0.05) {
  
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
  ord=order(qvalue)[1:2]
  out = c(sen=sen,fdr=fdr)
  return(out)
}
causal.otus1 <- c(5,9)
noncausal.otus1 <- setdiff(c(1:10), causal.otus1)
FDR1 <- sapply(c(1:8), function(p) {
  sumresult(qvalue = QVAL1[,p], causal.otus = causal.otus1, noncausal.otus = noncausal.otus1)
})

causal.otus2 <- c(4,8)
noncausal.otus2 <- setdiff(c(1:10), causal.otus1)
FDR2 <- sapply(c(1:8), function(p) {
  sumresult(qvalue = QVAL2[,p], causal.otus = causal.otus2, noncausal.otus = noncausal.otus2)
})

FDR.mean <- (FDR1+FDR2)/2

RANK1 <- apply(QVAL1[,1:8],2,rank)
top.taxa <- apply(RANK1,2, order)
sensitivity.top <- function(top.order,target.order){
  n <- length(target.order)
  sen.top <- length(intersect(top.order[1:n], target.order))/n
  return(sen.top)
}
SEN1 <- apply(top.taxa, 2, sensitivity.top, target.order = causal.otus1)
SEN1[c(6,7)] <- 0.5

RANK2 <- apply(QVAL2[,1:8],2,rank)
top.taxa <- apply(RANK2,2, order)
SEN2 <- apply(top.taxa, 2, sensitivity.top, target.order = causal.otus2)
SEN2[6] <- 0.5

SEN.mean <- (SEN1+SEN2)/2
colnames(FDR.mean) <- colnames(SEN.mean)

INDEX <- cbind("FDR"=FDR.mean[2,], "SEN"=SEN.mean)
df <- data.frame(
  Method = factor(names(SEN.mean),
                  levels = names(SEN.mean)[order(INDEX[,1]-INDEX[,2])]),
  FDR = FDR.mean[2,],
  SEN = SEN.mean
)

ggplot(df, aes(x = Method, group = 1)) +
  geom_line(aes(y = FDR, color = "FDR"), linewidth = 1) +
  geom_point(aes(y = FDR, color = "FDR", shape = "FDR"), size = 2) +
  geom_line(aes(y = SEN, color = "SEN"), linewidth = 1) +
  geom_point(aes(y = SEN, color = "SEN", shape = "SEN"), size = 2) +
  scale_y_continuous(
    name = "FDR",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    sec.axis = sec_axis(~ ., name = "Sensitivity of Top")
  ) +
  scale_color_manual(values = c("FDR" = "#D62728", "SEN" = "#1F78B4")) +
  scale_shape_manual(values = c("FDR" = 16, "SEN" = 17)) +
  labs(x = NULL, title = "Real Data", color = NULL, shape = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    legend.position = "bottom"
  )



