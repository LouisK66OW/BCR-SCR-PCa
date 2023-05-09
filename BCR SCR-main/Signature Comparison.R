library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(survival)
library(compareC)

load('data/signature.input.RData')
model<- signature_all[!duplicated(signature_all$Model),]$Model
signature_name <-  signature_all[!duplicated(signature_all$Author),]$Author

#signature
load("data/PCamachine_signaturescore.Rdata")

#CancerMap
CancerMap_time <- CancerMap[,c(1:3)]
CancerMap_gene<- CancerMap[,-c(1:3)]
data_exp<- data.frame()
signature<- list()
CancerMap_RS<- list()
CancerMap_RS$RS <- CancerMapRS

for(i in 1:length(model)) {
  print(i)
  signature[[i]] <- signature_all  %>%
    filter( Model ==  model[i] ) %>%
    select(c(1,4,6)) %>%
    column_to_rownames("ENSEMBL")
  data_exp<- CancerMap_gene[ ,colnames(CancerMap_gene) %in% rownames(signature[[i]])]
  for (k in colnames(data_exp)){
    print(k)
    data_exp[,k]<- data_exp[,k]*signature[[i]][k,2]
    
  }
  CancerMap_time[,i+3] <- rowSums(data_exp)
  CancerMap_RS[[i+1]]<- CancerMap_time[,c(2,3,i+3)] 
  colnames(CancerMap_RS[[i+1]])[3] <- model[i]
}
colnames(CancerMap_time)[-c(1:3)] <- signature_name
names(CancerMap_RS)[-1] <- signature_name

#CIT
CIT_time <- CIT[,c(1:3)]
CIT_gene<- CIT[,-c(1:3)]
data_exp<- data.frame()
CIT_RS <- list()
signature<- list()
CIT_RS<- list()
CIT_RS$RS <- CITRS

for(i in 1:length(model)) {
  print(i)
  signature[[i]] <- signature_all  %>%
    filter( Model ==  model[i] ) %>%
    select(c(1,4,6)) %>%
    column_to_rownames("ENSEMBL")
  data_exp<- CIT_gene[ ,colnames(CIT_gene) %in% rownames(signature[[i]])]
  for (k in colnames(data_exp)){
    print(k)
    data_exp[,k]<- data_exp[,k]*signature[[i]][k,2]
    
  }
  CIT_time[,i+3] <- rowSums(data_exp)
  CIT_RS[[i+1]]<- CIT_time[,c(2,3,i+3)] 
  colnames(CIT_RS[[i+1]])[3] <- model[i]
}
colnames(CIT_time)[-c(1:3)] <- signature_name
names(CIT_RS)[-1] <- signature_name

#cohorts
##CancerMap
cc_CancerMap <- data.frame(Cindex = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
                           HR = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$coefficients[1,2])}),
                           P = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$coefficients[1,5])}),
                           se = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
  rownames_to_column('ID')

##CIT
cc_CIT <- data.frame(Cindex = sapply(CIT_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
                     HR = sapply(CIT_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$coefficients[1,2])}),
                     P = sapply(CIT_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$coefficients[1,5])}),
                     se = sapply(CIT_RS,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
  rownames_to_column('ID')


########################merge##############################
######cancerMap####
rs_cancerMap<- merge(CancerMapRS,CancerMap_time,by.x = 0,by.y = 1)
rs_cancerMap<- rs_cancerMap %>%
  select(-c(5,6))  %>%
  column_to_rownames("Row.names")
colnames(rs_cancerMap)[1:2] <- c("OS.time","OS")
rt<- rs_cancerMap
cancerMap_compareC_p <- data.frame(Var = colnames(rt[, 4:105]), pval = c(1:length(colnames(rt[, 4:105]))))
for (i in colnames(rt[, 4:105])) {
  p <- compareC(rt$OS.time, rt$OS, rt$RS, rt[,i])$pval
  cancerMap_compareC_p[which(cancerMap_compareC_p$Var == i), 2] <- p
}

######CIT####
rs_CIT<- merge(CITRS,CIT_time,by.x = 0,by.y = 1)
rs_CIT<- rs_CIT %>%
  select(-c(5,6))  %>%
  column_to_rownames("Row.names")
colnames(rs_CIT)[1:2] <- c("OS.time","OS")

rt<- rs_CIT
CIT_compareC_p <- data.frame(Var=colnames(rt[, 4:105]),pval=c(1:length(colnames(rt[, 4:105]))))
for (i in colnames(rt[, 4:105])) {
  p <- compareC(rt$OS.time, rt$OS, rt$RS, rt[, i])$pval
  CIT_compareC_p[which(CIT_compareC_p$Var == i), 2] <- p
}

