detect.function <- function(data, Y){
  otu.table.filter <- data
  n.sam <- nrow(data)
  n.sam.control <- length(control.sample)
  n.sam.case <- length(case.sample)
  n.tax <- ncol(data)
  databind <- data.frame(Y = factor(Y), data)
  colnames(databind) <- c("group", taxaname)
  control <- which(Y==0)
  case <- which(Y==1)
  # Add Pseudo Count
  data[which(data==0)] <- 1
  
  
  pvalue <- matrix(0, n.tax, 10)
  row.names(pvalue) <- taxaname
  colnames(pvalue) <- c("Wilcoxon-CLR","MaAsLin2", "ALDEx2","LinDA",
                        "ANCOMBC","LOCOM","ZicoSeq","ADAPT","ADMIC-Basic","ADMIC")
  qvalue <- pvalue
  
  ###########
  # ADMIC
  ###########
  res.ADMIC <- ADMIC(otu.table.filter,Y)
  p0 <- res.ADMIC$p.basic
  p.eff <- res.ADMIC$p.admic
  q0 <- res.ADMIC$q.basic
  q.eff <- res.ADMIC$q.admic
  pvalue[,9] <- p0
  pvalue[,10] <- p.eff 
  qvalue[,9] <- q0
  qvalue[,10] <- q.eff 
  
  
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
  pvalue[,1] <- p.clr
  qvalue[,1] <- q.clr
  
  
  ###########
  # MaAsLin2
  ###########
  OTU=data.frame(otu.table.filter)
  rownames(OTU)=sam.ID
  
  res=Maaslin2(
    input_data = OTU,
    input_metadata = meta.data,
    output = "demo_output1",
    fixed_effects = 'group',
    plot_heatmap =FALSE,
    plot_scatter = FALSE)
  
  p.maas <- res$results[,c("feature","pval")]$pval
  q.maas <- res$results[,c("feature","qval")]$qval
  feature.maaslin2 <- res$results[,c("feature","qval")]$feature
  sx.maaslin2 <- match(taxaname, feature.maaslin2)
  p.maaslin2 <- p.maas[sx.maaslin2]
  q.maaslin2 <- q.maas[sx.maaslin2]
  pvalue[,2] <- p.maaslin2
  qvalue[,2] <- q.maaslin2
  
  ###########
  # ALDEx2
  ###########
  C <- NULL
  otu.table.aldex2 <- t(round(otu.table.filter))
  row.names(otu.table.aldex2) <- taxaname
  
  aldex.sample <- aldex.clr(otu.table.aldex2, as.character(Y), mc.samples=128, verbose=TRUE)
  res.aldex2 <- aldex.ttest(aldex.sample, paired.test=FALSE, hist.plot = FALSE)
  lll <- row.names(res.aldex2)
  sx.aldex2 <- match(taxaname, lll)
  p.aldex2 <- (res.aldex2$wi.ep)[sx.aldex2]
  q.aldex2 <- (res.aldex2$wi.eBH)[sx.aldex2]
  pvalue[,3] <- p.aldex2
  qvalue[,3] <- q.aldex2
  
  
  
  ###########
  # LinDA
  ###########
  data.linda<- t(otu.table.filter)
  res.linda <- linda(data.linda,meta.data,formula="~group",alpha=0.05)
  p.linda <- res.linda[["output"]][["group"]][["pvalue"]]
  q.linda <- res.linda[["output"]][["group"]][["padj"]]
  pvalue[,4] <- p.linda
  qvalue[,4] <- q.linda
  
  
  ###########
  # ANCOMBC
  ###########
  otu.table.ancombc <- t(otu.table.filter)
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
  pvalue[,5] <- p.ancombc
  qvalue[,5] <- q.ancombc
  
  ###########
  # LOCOM
  ###########
  res.locom <- locom(otu.table = otu.table.filter, Y = Y, C = NULL)
  locom.name <- colnames(res.locom$p.otu)
  sx.locom <- match(taxaname, locom.name)
  p.locom <- (res.locom$p.otu)[sx.locom]
  q.locom <- (res.locom$q.otu)[sx.locom]
  pvalue[,6] <- p.locom
  qvalue[,6] <- p.locom
  
  
  ###########
  # ZicoSeq
  ###########
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
  zicoseq.name <- row.names(q.zicoseq)
  sx.zicoseq <- match(taxaname, zicoseq.name)
  p.zicoseq <- p.zi[sx.zicoseq]
  q.zicoseq <- q.zi[sx.zicoseq]
  pvalue[,7] <- p.zicoseq
  qvalue[,7] <- q.zicoseq
  
  
  ###########
  # ADAPT
  ##########
  data.adapt <- t(matrix(as.integer(t(round(otu.table.filter))),n.tax,n.sam))
  colnames(data.adapt) <- taxaname
  Ya <- as.character(Y)
  Y.adapt <- as.data.frame(Ya)
  phy.data <- phyloseq(otu_table(data.adapt,taxa_are_rows = F),sample_data(Y.adapt))
  res.adapt <- adapt(phy.data, cond.var = "Ya")
  sx.adapt <- match(taxaname, res.adapt@details[["Taxa"]])
  p.adapt <- (res.adapt@details[["pval"]])[sx.adapt]
  q.adapt <- (res.adapt@details[["adjusted_pval"]])[sx.adapt]
  pvalue[,8] <- p.adapt
  qvalue[,8] <- q.adapt

  pvalue <- as.data.frame(pvalue)
  qvalue <- as.data.frame(qvalue)
  res <- list(pvalue, qvalue)
  
  return(res)
}