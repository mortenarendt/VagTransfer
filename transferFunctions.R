gt <- function(x){
  tb <- table(x$countC>0,x$count01) 
  dtb <- tb %>% dim()
  tb <- data.frame(n1 = dtb[1],n2 = dtb[2],nZero = sum(tb==0))
  return(tb)
}

getFisher <- function(x,doglm = TRUE){
  x$or <- (x$n00*x$n11) / (x$n10*x$n01)
  side <- ifelse(x$or<1,'less','greater')
  nn <- x$n00 + x$n10 + x$n01 + x$n11
  tb <- x[,c('n00','n01','n10','n11')] %>% 
    unlist() %>%
    matrix(2) %>%
    fisher.test(alternative = side) %>%
    tidy
  colnames(tb) <- paste('Fisher',colnames(tb),sep = '_')
  
  # bias corrected
  or_biascorr <- (x$n00 + 0.5)*(x$n11+0.5) / (x$n01 + 0.5)*(x$n10+0.5)
  t_biascorr <- (nn*(abs(x$n00*x$n11 - x$n10*x$n01) - 0.5*nn)^2) / ((x$n00 + x$n01)*(x$n00 + x$n10)*(x$n10 + x$n11)*(x$n01 + x$n11))
  pv_biascorr <- 2*(1 - pnorm(t_biascorr))
  biascorr <- data.frame(or_biascorr,t_biascorr,pv_biascorr)

  tb <- cbind(tb,biascorr)
  
  # glm
  if (doglm){
    mom <- c(rep(0,x$n00), rep(0,x$n01), rep(1,x$n10),rep(1,x$n11))
    child <- c(rep(0,x$n00), rep(1,x$n01), rep(0,x$n10),rep(1,x$n11))
    mdl <- glm(child~mom,family = binomial('logit')) %>% 
      tidy %>%
      filter(term=='mom') %>%
      mutate(or = exp(estimate)) %>%
      select(-term)
    colnames(mdl) <- paste('Glm',colnames(mdl),sep = '_')
    tb <- cbind(tb,mdl)
  } 
  return(tb)
}

truncateZerosInf <- function(or,trc = 100){
  trcm <- 10^-log10(trc)
  ornew <- or
  ornew[(is.infinite(or) & or>1) | or>trc ] <- trc
  ornew[(or==0 & or<1) | or<trcm] <- trcm
  return(ornew)
}

getWeigtedRatio <- function(x,e = 0.001){
  or <- x$Fisher_estimatetr
  pv <- -log10(x$Fisher_p.value)
  np <- sum(or>=1)
  nn <- sum(or<1)

  # mass
  area <- log(or)*pv
  ratio <- (sum(area[area>=0]) + e)/(-sum(area[area<0]) + e)
  R <- data.frame(np,nn , ratio)
  return(R)
}

shuffle <- function(nestedfactor){
  id <- 1:length(nestedfactor)
  df <- data.frame(nestedfactor,id,idnew=id)
  unf <- unique(df$nestedfactor)
  unf <- unf[!is.na(unf)]
  
  for (i in unf){
    ic <- df$nestedfactor==i & !is.na(df$nestedfactor)
    ids <- df$id[ic]
    df$idnew[ic] <- sample(ids)
  }
  ic <- is.na(df$nestedfactor)
  ids <- df$id[ic]
  df$idnew[ic] <- sample(ids)
  idnew <- df$idnew
  return(idnew)
}

getWeigtedRatio2 <- function(x,e = 0.001){
  or <- x$or_biascorr
  pv <- -log10(x$pv_biascorr)
  np <- sum(or>=1)
  nn <- sum(or<1)
  
  # mass
  area <- log(or)*pv
  ratio <- (sum(area[area>0]) + e)/(-sum(area[area<=0]) + e)
  R <- data.frame(np,nn , ratio)
  return(R)
}


extractPV <- function(permSTAT,modelratio,trm=100){
  niter <- dim(permSTAT)[3]
  tb <- c()
  if (is.na(niter)){
        xx <- permSTAT %>%
          t() %>% 
          data.frame() %>%
          mutate(Fisher_estimatetr = truncateZerosInf(Fisher_estimate,trm)) %>%
          do(getWeigtedRatio(x = .))
        tb <- rbind(tb,xx)
        pv <- sum(tb$ratio>modelratio)
  } else {
    for (i in 1:niter){
      xx <- permSTAT[,,i] %>%
        data.frame() %>%
        mutate(Fisher_estimatetr = truncateZerosInf(Fisher_estimate,trm)) %>%
        do(getWeigtedRatio(x = .))
      tb <- rbind(tb,xx)
    }
    pv <- sum(tb$ratio>modelratio) / niter
  }
  
  
  # estimate null distribution for ratio
  lgratio <- log(tb$ratio)
  SElgratio <- sqrt(sum(lgratio^2)/length(lgratio))

  print(c(i,median(tb$ratio),modelratio))
  df <- data.frame(pv, SElgratio, permmedian = median(tb$ratio), modelratio)
  
  return(df)
}