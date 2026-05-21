# data is an abundance matrix where the rows correspond to samples and the columns correspond to taxa. 
# Y is a vector of trait
# pseudo.count is the non-zero pseudo-count used to replace 0; it is recommended to use values such as 0.5 or 1. The default is 1.

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