# Kappa statistics for Fortschegger M et al comparing ddPCR and real time PCR for the detection of mutation in bone marrow. 

# Content based on reviewers' comments:
# Receiver Operating Characteristic (ROC) curve and the corresponding area under the curve (Az) comparing the two methods - STR-qPCR based and ddPCR based.  
# Three indexes of κ, positive agreement, and negative agreement.

##########################
# libraries and settings #
##########################

library("ggplot2")
library("reshape2")
library("cowplot")
library("pROC")
library("psych")

theme_set(theme_cowplot())

########
# data #
########

# these data are generated in the respective publication
dat <- read.table("data/STR_ddPCR.txt", header=TRUE)

############
# analysis #
############

# correlation plot
p1 <- ggplot(dat,aes(x=ddPCR, y=STR_PCR))+
  geom_point()+
  geom_smooth(method="lm")+
  NULL

cor(dat, method="spearman")
cor(dat, method="pearson") 

dffin <- NULL 
for (i in c(2,4,6,8)*10){
  cat <- (dat[,2] > i)*1
  roci <- roc(cat,dat[,1])
  auc <- roci$auc
  sens <- roci$sensitivities
  spec <- roci$specificities
  df <- data.frame(sensitivity=sens, specificity=spec,cutoff=rep(i,length(sens)))
  dffin <- rbind(dffin,df)
}

p2 <- ggplot(dffin,aes(x=1-specificity, y=sensitivity,colour=cutoff))+
  geom_line()+
  geom_abline(slope=1,colour="grey")+
  facet_wrap(~cutoff, nrow=2)

dffin <- NULL 
for (i in 1:99){
  cat <- (dat[,1] > i)*1
  roci <- roc(cat,dat[,2])
  auc <- as.numeric(roci$auc)
  vec <- c(cutoff=i,AUC=auc)
  dffin <- rbind(dffin, vec)
}
dffin <- data.frame(dffin)

p3 <- ggplot(dffin, aes(x=cutoff, y=AUC))+
  geom_line()+
  NULL


# Kappa statistics

cohen.kappa(dat)

# 
dffin <- NULL 
for (i in 1:99){
  cat <- (dat > i)*1
  kap <- cohen.kappa(cat)$kappa
  agree <- cohen.kappa(cat)$agree
  confid <- cohen.kappa(cat)$confid
  vec <- c(cutoff=i,Kappa=kap, 
           Agreement=(2*agree[1,1])/(2*agree[1,1]+agree[1,2]+agree[2,1]), 
           Disagreement=(2*agree[2,2])/(2*agree[2,2]+agree[1,2]+agree[2,1]), 
           Overall_agreement=(agree[1,1]+agree[2,2])/(agree[1,1]+agree[2,2]+agree[1,2]+agree[2,1]),
           lower=confid[1,1], upper=confid[1,3])
  dffin <- rbind(dffin, vec)
}
dffin <- data.frame(dffin)

# 95% confidence interval
p4 <- ggplot(dffin, aes(x=cutoff, y=Kappa))+
  geom_line()+
  geom_line(aes(x=cutoff, y=upper),colour="grey")+
  geom_line(aes(x=cutoff, y=lower),colour="grey")+
  NULL

p5 <- ggplot(dffin, aes(x=cutoff, y=Overall_agreement))+
  geom_line()+
  geom_line(aes(x=cutoff, y=Agreement),colour="blue")+
  geom_line(aes(x=cutoff, y=Disagreement),colour="grey")+
  NULL

pdf("plots/200124_kappa_stats.pdf", useDingbats = FALSE)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()


### calculation of positive and negative error
#             2a                    2d
# PA  =  ----------;    NA  =  ----------.    
#         2a + b + c            2d + b + c

#             a + d         a + d
# po  =  -------------  =  -----.    
#         a + b + c + d       N

# raters are defined as a two-by-two matrix of +/-


