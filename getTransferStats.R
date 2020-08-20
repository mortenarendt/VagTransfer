get2by2table <- function(o1,o2){
  n11 <- t((o1>0)+0) %*% (o2>0) %>% diag
  n00 <- t((o1==0)+0) %*% (o2==0) %>% diag
  n10 <- t((o1>0)+0) %*% (o2==0) %>% diag
  n01 <- t((o1==0)+0) %*% (o2>0) %>% diag
  STAT <- data.frame(n00,n10,n01,n11)
  return(STAT)
}

dyadmatch <- intersect(sample_data(phy1)$dyadnb,sample_data(phy2)$dyadnb)
otu1 <- subset_samples(phy1,sample_data(phy1)$dyadnb %in% dyadmatch)
otu2 <- subset_samples(phy2,sample_data(phy2)$dyadnb %in% dyadmatch) 
deptho1 <- sample_sums(otu1)
totAbuM <- sum(deptho1)
deptho2 <- sample_sums(otu2)
totAbuC <- sum(deptho2)
# remove non-testabel OTU's
o1 <- otu1 %>% otu_table()
o2 <- otu2 %>% otu_table()
s01 <- apply(o1==0,1,sum)
s02 <- apply(o2==0,1,sum)
nn <- dim(o2)[2]
icvar <- s01>0 & s02>0 & s01<nn & s02<nn
otu1 <- subset_taxa(otu1,icvar)
otu2 <- subset_taxa(otu2,icvar)
p <- sum(icvar)
print(c(nn,p))
######
o1 <- otu1 %>% otu_table() %>% t() 
o2 <- otu2 %>% otu_table() %>% t()
o1 <- o1[order(rownames(o1)),]
deptho1 <- deptho1[order(names(deptho1))] 
o2 <- o2[order(rownames(o2)),]
deptho2 <- deptho2[order(names(deptho2))] 



abuM <- apply(o1,2,sum)
abuC <- apply(o2,2,sum)
abuMrel <- o1 %>% sweep(1,deptho1,'/') %>% apply(2,mean)
abuCrel <- o2 %>% sweep(1,deptho2,'/') %>% apply(2,mean)

STAT <- get2by2table(o1,o2)
STAT <- data.frame(STAT,abuM,abuMrel,abuC,abuCrel)
STAT$otu <- rownames(STAT)

# make permutation testing
if (nperm>0){
  permSTAT <- array(dim = c(dim(STAT)[1],4,nperm))
  for (i in 1:nperm){
    stp <- get2by2table(o1[sample(dim(o1)[1]),],o2)
    permSTAT[,,i] <-  as.matrix(stp)
  }
  dimnames(permSTAT)[[1]] <- rownames(stp)
  dimnames(permSTAT)[[2]] <- colnames(stp)
  
  permSTATfisher <- array(dim = c(dim(STAT)[1],3,nperm))
  for (i in 1:nperm){
    print(i)
    aa <- permSTAT[,,i] %>% data.frame()
    for (j in 1:dim(permSTAT)[1]){
      stp2 <- getFisher(aa[j,],doglm = F)
      # permSTATfisher[j,,i] <- c(stp2$Fisher_estimate, stp2$Fisher_p.value)
      # permSTATfisher[j,,i] <- c(stp2$Fisher_estimate,stp2$or_biascorr, stp2$Gtest_p.value, stp2$Fisher_p.value)
      permSTATfisher[j,,i] <- c(stp2$Fisher_estimate,stp2$or_biascorr, stp2$Fisher_p.value)
    }
  }
  dimnames(permSTATfisher)[[1]] <- rownames(stp)
  # dimnames(permSTATfisher)[[2]] <- c('Fisher_estimate', 'Fisher_p.value')
  # dimnames(permSTATfisher)[[2]] <- c('Fisher_estimate','or_biascorr', 'Gtest_p.value', 'Fisher_p.value')
  dimnames(permSTATfisher)[[2]] <- c('Fisher_estimate','or_biascorr', 'Fisher_p.value')
}

STAT <- STAT %>%
  group_by(otu) %>%
  do(getFisher(x = .)) %>%
  ungroup %>%
  left_join(STAT,by = 'otu')

TAXtb <- data.frame(tax_table(phy1))
STAT <- merge(STAT,TAXtb,by.x = 'otu',by.y = 'row.names')
STAT$totAbuM <- totAbuM
STAT$totAbuC <- totAbuC  

STAT <- STAT %>%
  mutate(Fisher_estimatetr =  truncateZerosInf(Fisher_estimate),
         Glm_ortr =  truncateZerosInf(Glm_or))
