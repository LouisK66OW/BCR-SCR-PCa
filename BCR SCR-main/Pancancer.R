Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

library(circlize)
library(ComplexHeatmap)

SampleInfo <- read.csv("")
anno<- openxlsx::read.xlsx("")
anno<- unique(anno$Cohort)

#color
Cohort.Col <- c("ACC" = "#F5CCDC", "PCPG" = "#FED9EA", "PAAD" = "#F6B5B6", "THCA" = "#EEDEED",
                "LGG" = "#957D84", "GBM" = "#D1A6A8", "UVM" = "#7B7C75", "MESO" = "#F7F7BC",
                "SKCM" = "#D0D36F", "SARC" = "#E5E6AB", "UCS" = "#88C4A8", "UCEC" = "#79D9E4",
                "OV" = "#D8F3FA", "CESC" = "#83C3E4", "BRCA" = "#91BBE5", "TGCT" = "#71A5D0",
                "BLCA" = "#A0C2DD", "DLBC" = "#ED5C64", 
                "THYM" = "#A55A5E", "KICH" = "#F7E897", "KIRP" = "#BDB279", "KIRC" = "#E3E0B1",
                "ESCA" = "#DECEE8", "STAD" = "#CAB4D8", "HNSC" = "#FBF7B8", "LIHC" = "#F0A363",
                "CHOL" = "#AA5F39", "LUAD" = "#979599", "LUSC" = "#CCD4D6", "READ" = "#98D38F",
                "COAD" = "#8BCF9D")

Class1.Col <- c("Secretory gland" = "#B1CE46",
                "Soft Tissue" = "#BB9727",
                "Female\nreproductive\norgan" = "#DF9D73",
                "Male genitourinary system" = "#38679A",
                "Hematopoietic\nand immune" = "#32B897",
                "Digestive\nsystem" = "#CCE6F1",
                "Hepatobiliary system" = "#D8383A")

Class2.Col <- c("Adrenal\ngland" = "#957D84", "Thyroid" = "#E6CEE0",
                "Brain" = "#866460", "Eye" = "#585D5B",
                "Skin" = "#BCBD30", "Uterus" = "#21C2CE",
                "Ovary" = "#B6EDF9", "Breast" = "#2878B5",
                "Kidney" = "#F3D266",
                "Stomach" = "#C6B0D6", "Head\nand neck" = "#BDB279",
                "Liver" = '#FCC88F', "Lung" = "#A9A8A8",
                "Colon" = "#6C8E59")

Col = c(Cohort.Col, Class1.Col, Class2.Col)
CS.Col = colorRamp2(breaks = c(-2,0,4), c("#F7FBFC","#4DBBD5B2", "#00A087B2"))
SectorColor <- function(sector, col, track.index, sampleInfo, currentCol, standardCol = "Cohort", ...){
  cohorts.to.plot = unique(SampleInfo[[standardCol]][SampleInfo[[currentCol]] %in% sector])
  st.degree = max(unlist(lapply(cohorts.to.plot, get.cell.meta.data, name = "cell.start.degree")))
  ed.degree = min(unlist(lapply(cohorts.to.plot, get.cell.meta.data, name = "cell.end.degree")))
  
  draw.sector(start.degree = st.degree, end.degree = ed.degree,
              rou1 = get.cell.meta.data("cell.top.radius", track.index = track.index),
              rou2 = get.cell.meta.data("cell.bottom.radius", track.index = track.index),
              col = col[sector], ...)
}

SectorLabel <- function(sector, col, track.index, sampleInfo, currentCol, standardCol = "Cohort", niceFacing = T, ...){

  plot.data = sampleInfo[sampleInfo[[currentCol]] %in% sector, ]
  plot.data = plot.data[round(0.5*nrow(plot.data)), ]
  circos.text(x = plot.data$ID, y = 0.5, labels = sector, sector.index = plot.data$Cohort, track.index = track.index, col = col, niceFacing = niceFacing, ...)
}

SampleInfo$Cohort <- factor(SampleInfo$Cohort, levels = anno) 
SampleInfo$ID <- 1:nrow(SampleInfo)

###
circos.par("cell.padding" = c(0.02, 0.00, 0.02, 0.00), "start.degree" = 90) 
circos.initialize(sectors = SampleInfo$Cohort, x = SampleInfo$ID) 

### 1
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, cell.padding = c(0, 1))
for (i in 1:nrow(SampleInfo)) circos.rect(xleft = SampleInfo$ID[i], xright = SampleInfo$ID[i]+1, ybottom = 0, ytop = 1, 
                                          sector.index = SampleInfo$Cohort[i], track.index = 1, border = NA, col = CS.Col(SampleInfo$BCR.score[i]))

### 2
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2)
for (i in unique(SampleInfo$Cohort)) SectorColor(sector = i, col = Col, track.index = 2, sampleInfo = SampleInfo, currentCol = "Cohort") 
for (i in unique(SampleInfo$CohortLabel)) SectorLabel(sector = i, col = "black", track.index = 2, sampleInfo = SampleInfo, currentCol = "CohortLabel") 

### 3
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2)
for (i in unique(SampleInfo$Class2)) SectorColor(sector = i, col = Col, track.index = 3, sampleInfo = SampleInfo, currentCol = "Class2") 
for (i in unique(SampleInfo$Class2Label)) SectorLabel(sector = i, col = "white", track.index = 3, sampleInfo = SampleInfo, currentCol = "Class2Label") 

### 4

circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.lty = "blank")
for (i in unique(SampleInfo$Class1)) SectorColor(sector = i, col = Col, track.index = 4, sampleInfo = SampleInfo, currentCol = "Class1", lty = "blank") 
for (i in unique(SampleInfo$Class1Label)) SectorLabel(sector = i, col = "black", track.index = 4, sampleInfo = SampleInfo, currentCol = "Class1Label", niceFacing = F, facing = "clockwise", pos = 2) 
circos.clear()


lgd = Legend(col_fun = CS.Col, title = "BCR.score")
draw(lgd, x = unit(0.87, "npc"), y = unit(0.14, "npc"))



