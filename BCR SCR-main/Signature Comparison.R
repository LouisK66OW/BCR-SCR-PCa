
rm(list = ls())

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)


library(dplyr)
library(tidyr)
library(tibble)
library(data.table)


setwd("")

path <- "filtersubset/"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x,row.names = 1,header = T)})  



for (i in 1:length(data)){
  data[[i]] <- as.data.frame(t(data[[i]]))
}



for (i in 1:length(data)){
  rownames(data[[i]]) <- data[[i]]$Row.names
  
  colnames(data[[i]])[1:3] <- c("ID","BCR.time","BCR")
  data[[i]]<- mutate( data[[i]],BCR.time = BCR.time/12)
  
}

###102 signature
signature_all <- read.csv("")


CancerMap<- data[[1]]
CIT<- data[[2]]
CPC<- data[[3]]
GSE54460<- data[[4]]
Stockholm<- data[[5]]
Taylor<- data[[6]]
TCGA<- data[[7]]


library(survival)


load("machine_score.Rdata")

rs<- rs3
names(rs)
CancerMapRS<- rs[["CancerMap"]]
CITRS<- rs[["CIT"]]
CPCRS<- rs[["CPC"]]
GSE54460RS<- rs[["GSE54460"]]
StockholmRS<- rs[["Stockholm"]]
TaylorRS<- rs[["Taylor"]]
TCGARS<- rs[["TCGA"]]


#CancerMap
CancerMap_time <- CancerMap[,c(1:3)]
CancerMap_gene<- CancerMap[,-c(1:3)]
data_exp<- data.frame()
signature<- list()
CancerMap_RS<- list()
CancerMap_RS$RS <- CancerMapRS

###3- Calculate the C-index
library(survival)

##CancerMap
cc_CancerMap <- data.frame(Cindex = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR)~., x))$concordance[1])}),
                           P = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR)~., x))$coefficients[1,5])}),
                           se = sapply(CancerMap_RS,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR)~.,x))$concordance[2])})) %>%
  rownames_to_column('ID')

################Average C-index of 7 datasets

cc_CancerMap_1<- cc_CancerMap 
colnames(cc_CancerMap_1)[2,5]<- c("CancerMap_cindex","CancerMap_se")

cc_CIT_1<- cc_CIT
colnames(cc_CIT_1)[2,5]<- c("CIT_cindex","CIT_se")

cc_CPC_1<- cc_CPC
colnames(cc_CPC_1)[2,5]<- c("CPC_cindex","CPC_se")

cc_GSE54460_1<- cc_GSE54460
colnames(cc_GSE54460_1)[2,5]<- c("GSE54460_cindex","GSE54460_se")

cc_Stockholm_1<- cc_Stockholm
colnames(cc_Stockholm_1)[2,5]<- c("Stockholm_cindex","Stockholm_se")

cc_Taylor_1<-cc_Taylor
colnames(cc_Taylor_1)[2,5]<- c("Taylor_cindex","Taylor_se")

cc_TCGA_1<-cc_TCGA
colnames(cc_TCGA_1)[2,5]<- c("TCGA_cindex","TCGA_se")


data1<- merge(cc_CancerMap_1,cc_CIT_1,by=1)
data2<- merge(data1,cc_CPC_1,by=1)
data3<- merge(data2,cc_GSE54460_1,by=1)
data4<- merge(data3,cc_Stockholm_1,by=1)
data5<- merge(data4,cc_Taylor_1,by=1)
data6<- merge(data5,cc_TCGA_1,by=1)
colnames(data6)
data6_cindex<- data6[,c(1,2,4,6,8,10,12,14)]
data6_se<- data6[,c(1,3,5,7,9,11,13)]
data6_cindex<- data6_cindex %>%
  mutate(meancindex = rowMeans(.[,-1]))

data6_se<- data6_se %>%
  mutate(meanse = rowMeans(.[,-1]))

data66<- merge(data6_cindex,data6_se,by = 1)

colnames(data66)
data66<- data66[,c(1,9,16)]
data66<- data66 %>%
  arrange(desc(meancindex))


###The second step is to calculate the significance of the difference in C-index and RS between the other models.
library(compareC)

########################Merging your own model risk score with external signatures.##############################
######TCGA####
rs_TCGA<- merge(TCGARS,TCGA_time,by.x = 0,by.y = 1)

rs_TCGA<- rs_TCGA %>%
  select(-c(5,6))  %>%
  column_to_rownames("Row.names")

colnames(rs_TCGA)[1:2] <- c("BCR.time","BCR")

rt<- rs_TCGA

TCGA_compareC_p <- data.frame(Var = colnames(rt[, 4:105]), pval = c(1:length(colnames(rt[, 4:105]))))
for (i in colnames(rt[, 4:105])) {
  p <- compareC(rt$BCR.time, rt$BCR, rt$RS, rt[,i])$pval
  TCGA_compareC_p[which(TCGA_compareC_p$Var == i), 2] <- p
}


cancerMap_compareC_p <- rbind(c("RS", 1), cancerMap_compareC_p)
CIT_compareC_p <- rbind(c("RS", 1), CIT_compareC_p)
CPC_compareC_p <- rbind(c("RS", 1), CPC_compareC_p)
GSE54460_compareC_p  <- rbind(c("RS", 1), GSE54460_compareC_p)
Stockholm_compareC_p <- rbind(c("RS", 1), Stockholm_compareC_p)
Taylor_compareC_p <- rbind(c("RS", 1), Taylor_compareC_p)
TCGA_compareC_p  <- rbind(c("RS", 1), TCGA_compareC_p)

##merge
all_cancerMap <- data.frame(cc_CancerMap , p = cancerMap_compareC_p[, 2])
all_CIT <- data.frame(cc_CIT, p = CIT_compareC_p[, 2])
all_CPC  <- data.frame(cc_CPC, p = CPC_compareC_p[, 2])
all_GSE54460 <- data.frame(cc_GSE54460, p = GSE54460_compareC_p[, 2])
all_Stockholm <- data.frame(cc_Stockholm, p = Stockholm_compareC_p[, 2])
all_Taylor <- data.frame(cc_Taylor, p = Taylor_compareC_p[, 2])
all_TCGA <- data.frame(cc_TCGA, p = TCGA_compareC_p[, 2])


