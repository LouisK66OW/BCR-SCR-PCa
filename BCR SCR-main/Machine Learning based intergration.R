rm(list = ls())
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(tidyverse)
library(CoxBoost)
library(survivalsvm)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(BART)

load('')

##rename
names(data) <- c("CancerMap","CIT","CPC","GSE54460","Stockholm","Taylor","TCGA")
data <- lapply(data,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})
data_table<- data.frame()
data_train <- data$TCGA
data_all <- data
data_name <- colnames(data_train)[-c(1:3)]
data_train_surv <- data_train[,c('BCR.time','BCR',data_name)]
data_all_surv <- lapply(data_all,function(x){x[,c('BCR.time','BCR',data_name)]})

nodesize_value <- 5
seed <- 1
################################################ 1-1.RSF ################################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~.,data = data_train_surv,
             ntree = 2000,nodesize = nodesize_value,    #此参数可调整
             mtry = 4,
             nsplit = 10,
             block.size = 10,
             splitrule = 'logrank',
             importance = TRUE,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- 'RSF'
data_table<- rbind(data_table,cc)

################################################ 2-1.Enet ################################################
x1 <- as.matrix(data_train_surv[,data_name])
x2 <- as.matrix(Surv(data_train_surv$BCR.time,data_train_surv$BCR))
for (alpha in seq(0.1,0.9,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$model<- paste0('Enet','[α=',alpha,']')
  data_table<- rbind(data_table,cc)
}

################################################ 3-1.lasso ################################################
x1 <- as.matrix(data_train_surv[,data_name])
x2 <- as.matrix(Surv(data_train_surv$BCR.time,data_train_surv$BCR))
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- "lasso"
data_table<- rbind(data_table,cc)

############################################### 3-2.lasso+RSF ################################################
x1 <- as.matrix(data_train_surv[,data_name])
x2 <- as.matrix(Surv(data_train_surv$BCR.time,data_train_surv$BCR))
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
coef<- coef(fit,s= fit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
rid=row.names(coef)[index]
data_train_surv2 <- data_train[,c('BCR.time','BCR',rid)]
data_all_surv2 <- lapply(data_all,function(x){x[,c('BCR.time','BCR',rid)]})
fit <- rfsrc(Surv(BCR.time,BCR)~.,data = data_train_surv2,
             ntree = 2000,nodesize = nodesize_value,        ##此处参数建议调整 
             mtry = 4,
             nsplit = 10,
             block.size = 10,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(data_all_surv2,function(x){cbind(x[,1:2],score=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<-  paste0('lasso',' + RSF')
data_table<- rbind(data_table,cc)

################################################ 4.Ridge ################################################
x1 <- as.matrix(data_train_surv[,data_name])
x2 <- as.matrix(Surv(data_train_surv$BCR.time,data_train_surv$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=0,nfolds = 10)
rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- "Ridge"
data_table<- rbind(data_table,cc)

################################################ 5.GBM ################################################
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~.,data = data_train_surv,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~.,data = data_train_surv,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- paste0('GBM')
data_table<- rbind(data_table,cc)

################################################ 6-1.StepCox ################################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(BCR.time,BCR)~.,data_train_surv),direction = direction)
  rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$model<- paste0('StepCox','[',direction,']')
  data_table<- rbind(data_table,cc)
}

############################################### 7-1.CoxBoost ################################################
set.seed(seed)
pen <- optimCoxBoostPenalty(data_train_surv[,'BCR.time'],data_train_surv[,'BCR'],as.matrix(data_train_surv[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(data_train_surv[,'BCR.time'],data_train_surv[,'BCR'],as.matrix(data_train_surv[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(data_train_surv[,'BCR.time'],data_train_surv[,'BCR'],as.matrix(data_train_surv[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- paste0('CoxBoost')
data_table<- rbind(data_table,cc)

################################################ 7-2.CoxBoost+RSF ################################################
set.seed(seed)
pen <- optimCoxBoostPenalty(data_train_surv[,'BCR.time'],data_train_surv[,'BCR'],as.matrix(data_train_surv[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(data_train_surv[,'BCR.time'],data_train_surv[,'BCR'],as.matrix(data_train_surv[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(data_train_surv[,'BCR.time'],data_train_surv[,'BCR'],as.matrix(data_train_surv[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)

summary(fit)
coef<- coef(fit)
index <- which(coef != 0)
actCoef <- coef[index]
rid=names(coef)[index]
data_train_surv2 <- data_train[,c('BCR.time','BCR',rid)]
data_all_surv2 <- lapply(data_all,function(x){x[,c('BCR.time','BCR',rid)]})
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~.,data = data_train_surv2,
             ntree = 2000,nodesize = nodesize_value,       ##此处参数建议调整
             mtry = 4,
             nsplit = 10,
             block.size = 10,
             splitrule = 'logrank',
             importance = TRUE,
             proximity = T,
             forest = T,
             seed = seed)

rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost','+ RSF')
data_table <- rbind(data_table,cc)

################################################ 8.plsRcox ################################################
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=data_train_surv[,data_name],time=data_train_surv$BCR.time,status=data_train_surv$BCR),nt=10,verbose = FALSE)
fit <- plsRcox(data_train_surv[,data_name],time=data_train_surv$BCR.time,event=data_train_surv$BCR,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- paste0('plsRcox')
data_table<- rbind(data_table,cc)

################################################ 9.superpc ################################################
data <- list(x=t(data_train_surv[,-c(1,2)]),y=data_train_surv$BCR.time,censoring.status=data_train_surv$BCR,featurenames=colnames(data_train_surv)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(data_all_surv,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$BCR.time,censoring.status=w$BCR,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],score=rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- paste0('SuperPC')
data_table<- rbind(data_table,cc)

################################################ 10.survivalsvm ################################################
fit = survivalsvm(Surv(BCR.time,BCR)~., data= data_train_surv, gamma.mu = 1)

rs <- lapply(data_all_surv,function(x){cbind(x[,1:2],score=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~score,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$model<- paste0('survival-SVM')
data_table<- rbind(data_table,cc)
data_table2 <- data_table
data_table2$model<- gsub('α','a',data_table2$Model)

dd2 <- pivot_wider(data_table2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)

dd2<- dd2 %>%
  column_to_rownames("Model")

########visualization########
library(ComplexHeatmap)
library(BART)
library(snowfall)
library(RColorBrewer)

Cindex_mat <- dd2
avg_Cindex <- apply(Cindex_mat, 1, mean)          
avg_Cindex <- sort(avg_Cindex, decreasing = T)    
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]     
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,stat = 'identity',
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 11, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          
                          show_annotation_name = F)

cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              col = c("#4195C1", "#FFFFFF", "#CB5746"), 
              rect_gp = gpar(col = "black", lwd = 1), 
              cluster_columns = FALSE, cluster_rows = FALSE, 
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 14),
              
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 11))
              }
)

pdf(file.path(fig.path, "machine_Cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 3.5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())
