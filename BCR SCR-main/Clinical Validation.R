options(stringsAsFactors = F)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(timeROC)

load('')

###KM
risk<- TCGA
rt <- risk %>% 
  filter(BCR.time >= 1) %>% 
  mutate(BCR.time = BCR.time/12)
colnames(risk)
res.cut <- surv_cutpoint(rt, 
                         time = "BCR.time",
                         event = "BCR", 
                         variables = "riskscore", 
                         minprop = 0.3
                         
)                      
summary(res.cut) 
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(BCR.time, BCR) ~riskscore, data = res.cat)
#plot
pp <- ggsurvplot(fit,
                 data = res.cat,
                 size = 1,
                 palette = c("#DC0000FF", "#00A087FF"),
                 conf.int = F,
                 xlab ="Time(years)",
                 ylab = "Biochemical Recurrence Free Survival",
                 risk.table = T,
                 tables.theme = theme_bw(),
                 pval.method=TRUE,
                 risk.table.col = "strata",
                 font.legend = 19,
                 font.title = 19,font.x = 19,font.y = 19,
                 pval = TRUE,legend.labs=c("high","low"),ncensor.plot=F,
                 font.tickslab = 15,
                 pval.size=6.5,
                 pval.method.size=6.5,
                 fontsize=4.5
                 
)

####ROC
timeroc<-timeROC(T=dat$BCR.time,                        
                 delta=dat$BCR, marker=dat[,4],
                 cause=1,weighting="marginal",
                 times=c(6,8,10),           
                 ROC = TRUE,
                 iid=TRUE)


pdf(file="TCGA_ROC.pdf", width=5, height=5)
plot(timeroc, 
     time=6, col="#DC0000FF", lwd=2, title = "")  
plot(timeroc,
     time=8, col="#00A087FF", add=TRUE, lwd=2)    
plot(timeroc,
     time=10, col="#F39B7FFF", add=TRUE, lwd=2)
#add legends
legend("bottomright",
       c(paste0("6-year = ",round(timeroc[["AUC"]][1],3)), 
         paste0("8-year = ",round(timeroc[["AUC"]][2],3)), 
         paste0("10-year = ",round(timeroc[["AUC"]][3],3))),
       col=c("#DC0000FF", "#00A087FF", "#F39B7FFF"),
       lty=1, lwd=2,bty = "n")   

dev.off()

####Cox analysis
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file=forestFile, width=8, height=3)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
  abline(v=1, col="black", lty=2, lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
  axis(1)
  dev.off()
}

indep=function(expFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)     
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     
  sameSample=intersect(row.names(cli),row.names(exp))
  exp=exp[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(exp, cli)
  
  #Univariate analysis
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(BCR.time, BCR) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="#00A087FF")
  
  #Multivariate analysis
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.2,]
  rt1=rt[,c("BCR.time", "BCR", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(BCR.time, BCR) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="#DC0000FF")
}

indep(expFile="TCGArisk.txt",
      cliFile="TCGAclinical.txt",
      uniOutFile="TCGAuniCox.txt",
      multiOutFile="TCGAmultiCox.txt",
      uniForest="TCGAuniForest.pdf",
      multiForest="TCGAmultiForest.pdf")

dev.off()





