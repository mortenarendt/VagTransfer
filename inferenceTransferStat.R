
###### Some functions

get2by2table <- function(o1,o2){
  n11 <- t((o1>0)+0) %*% (o2>0) %>% diag
  n00 <- t((o1==0)+0) %*% (o2==0) %>% diag
  n10 <- t((o1>0)+0) %*% (o2==0) %>% diag
  n01 <- t((o1==0)+0) %*% (o2>0) %>% diag
  STAT <- data.frame(n00,n10,n01,n11)
  return(STAT)
}


getRatios <- function(o1,o2,IDfactor){
  unFac <- unique(IDfactor)
  wr <- c()
  
  for (i in unFac){
    # print(i)
    STAT <- get2by2table(o1[IDfactor==i,],o2[IDfactor==i,])
    STAT$otu <- rownames(STAT)
    ID <- 1:dim(STAT)[1]
    ID <- (STAT$n00 + STAT$n10) > 0 & 
      (STAT$n00 + STAT$n01) > 0 &
      (STAT$n11 + STAT$n01) > 0 &
      (STAT$n11 + STAT$n10) > 0 
    ID <- which(ID)
    
    stsdf <- c()
    for (j in ID){
      # print(j)
      df <- getFisher(STAT[j,],doglm = F)
      df$otu <- STAT$otu[j]
      stsdf <- rbind(stsdf,df)
    }
    
    STAT <- merge(STAT,stsdf, by = 'otu')
    STAT <- STAT %>%
      mutate(Fisher_estimatetr = truncateZerosInf(Fisher_estimate))
    wr <- rbind(wr,
                data.frame(delivery = i,getWeigtedRatio(STAT)))
  }
  return(wr)
}

### Match data
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

o1 <- otu1 %>% otu_table() %>% t() 
o2 <- otu2 %>% otu_table() %>% t()
metaData <- sample_data(otu1)
order_o1 <- order(rownames(o1))
o1 <- o1[order_o1,]
metaData <- metaData[order_o1,]
deptho1 <- deptho1[order(names(deptho1))] 
o2 <- o2[order(rownames(o2)),]
deptho2 <- deptho2[order(names(deptho2))] 

abuM <- apply(o1,2,sum)
abuC <- apply(o2,2,sum)
abuMrel <- o1 %>% sweep(1,deptho1,'/') %>% apply(2,mean)
abuCrel <- o2 %>% sweep(1,deptho2,'/') %>% apply(2,mean)

IDfactor <- metaData$DELIVERY
IDfactor[IDfactor!='Normal'] <- 'Sectio'

# model
wr <- getRatios(o1,o2,IDfactor)
wr <- wr %>% 
  mutate(ratioratio = ratio[delivery=='Normal'] / ratio[delivery=='Sectio']) 
colnames(wr)[2:5] <-paste('Model',colnames(wr)[2:5],sep = '_' )

# # permutation
# wrperm <- c()
# for (jj in 1:nperm){
#   wr2 <- getRatios(o1,o2,sample(IDfactor))
#   wrperm <- rbind(wrperm,data.frame(permutation = jj,wr2))
# }
# 

registerDoMC()
getDoParWorkers()

wrperm <-foreach (jj=1:nperm,.combine = rbind) %dopar% {
  wr2 <- getRatios(o1,o2,sample(IDfactor))
}
wrperm <- data.frame(permutation = sort(rep(1:nperm,2)),wrperm)

wrperm <- wrperm %>% 
  group_by(permutation) %>%
  mutate(Perm_ratioratio = ratio[delivery=='Normal'] / ratio[delivery=='Sectio']) %>%
  left_join(wr,by = 'delivery')

