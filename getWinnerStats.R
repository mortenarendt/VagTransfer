
getRankedPrc <- function(X1,X2){
  mX2 <- X2 %>% 
    gather(otu,Ccount,-c(dyadnb,Time,Type,DELIVERY,delivery,CST_w36)) %>%
    select(otu,Ccount,dyadnb)
  
  mX1 <- X1 %>% 
    gather(otu,Mcount,-c(dyadnb,Time,Type,DELIVERY,delivery,CST_w36)) %>%
    group_by(dyadnb) %>%
    arrange(desc(Mcount)) %>%
    mutate(rnk = 1:n()) %>%
    
    left_join(mX2,by = c('otu','dyadnb')) %>%
    filter(!is.na(Ccount)) %>%
    group_by(delivery,rnk) %>%
    summarise(n = n(), 
              nC = sum(Ccount>0),
              nC1 = sum(Ccount>1),
              nM = sum(Mcount>0)) %>%
    ungroup() %>%
    filter(nM==n) 
  return(mX1)
}

getPermStats <- function(x){
  # x <- mX1[1,]
  prcM <- x$prcModel
  prcPerm <- x[,8:dim(x)[2]]
  pv <- sum(prcPerm>prcM) / length(prcPerm)
  df <- data.frame(x[,1:7],niter = length(prcPerm), pv)
  return(df)
}



######## Set Data

sd2 <- data.frame(sample_data(phy2)) %>% mutate(delivery = DELIVERY, delivery = replace(delivery, delivery!='Normal','Sectio'))
otu2 <- t(otu_table(phy2))
otu2 <- otu2[,colnames(otu2) %in% colnames(X1)]
X2 <- data.frame(sd2,otu2)


# Model 
mX1 <- getRankedPrc(X1,X2) %>% mutate(prcModel = nC/n)

# registerDoMC()
# permutation
# for (i in 1:nperm){
RR <- foreach (i=1:nperm) %do% {
  # print(i)
  X2r <- data.frame(sd2[shuffle(sd2$delivery),],otu2) 
  mX1r <- getRankedPrc(X1,X2r) %>% 
    mutate(prcP = nC/n) %>%
    select(delivery,rnk,prcP) 
  # mX1 <- merge(mX1,mX1r,by = c('delivery','rnk'))
}
mX1r <- Reduce(function(x, y) {merge(x, y, by = c('delivery','rnk'))}, RR)
mX1 <- merge(mX1,mX1r,by = c('delivery','rnk'))

tb <- mX1 %>% 
  group_by(delivery,rnk) %>%
  do(getPermStats(x = .))


# ggplot(data =mX1, aes(rnk,prcModel,color = delivery, group = factor(delivery))) + 
#   geom_point() + 
#   geom_line() + 
#   geom_text(data = tb, aes(label = paste('p =' ,pv)), color = 'black')
# # scale_x_log10() + 
# # facet_wrap(~count_chld)  
# 
# 
# 
