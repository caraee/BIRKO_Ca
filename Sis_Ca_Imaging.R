library("tidyverse")
library("here")
library("mclust")
library("extrafont")
library("RColorBrewer")
library("scales")
library("grid")
library("gtable")
library("data.table")
library("readxl")
library("janitor")
library("cramer")
library("pheatmap")

temp <- list.files(pattern="*.csv") # Find all files in your folder that are of type .txt

list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv, as.is=T, check.names=F), envir = .GlobalEnv)

lapply(names(temp),function(nm){
  names(get(nm))<-paste(nm,names(get(nm),sep="_"))
})

temp<-gsub("*.csv$", "", temp)

image_list<-lapply(temp,get)
names(image_list)<-temp
image_list<-lapply(image_list,function(df){
  n<-nrow(df)
  i<-459-n
  if(i>=0) {df[n+i,] <- NA}
  df
  })

cells<-bind_cols(image_list)
cells$time<-seq(0,458*10/60,10/60)
cells[cells==-1]<-NA

tests<-read.table("tests.txt",sep="\t",stringsAsFactors = F, 
                  header=T, check.names=F,as.is=T)
tests$Time<-seq(0,458*10/60,10/60)

animalTable<-read.table("animalTable.txt",sep="\t",stringsAsFactors = F, 
                  header=T, check.names=F,as.is=T)
animalTable$Image<-paste("Image",animalTable$Image,sep="")

cells_full<-cells
#tail(cells[,"Image36_Cell1"])
cells<-cells[1:429,]
#tests<-tests[1:429,]

#cells<-cells_full
#filter bad cells
#  which.min(cells_long$Ratio[cells_long$Image=="Image27"])
#  which.max(cells_long$Ratio[cells_long$Image=="Image27"])
# 
#  which.max(cells_long$Ratio[cells_long$Image=="Image29"&cells_long$time>70])
# 
# cells_long[cells_long$Image=="Image72"&
#               cells_long$time>70,][which.max(cells_long$Ratio[cells_long$Image=="Image72"&
#                                                                 cells_long$time>70]),]
# 
#  cells_long[cells_long$Image=="Image75"&
#               cells_long$Stimulus=="baseline",][which.max(cells_long$Ratio[cells_long$Image=="Image75"&
#                                                                         cells_long$Stimulus=="baseline"]),]
# 
# 
#  p<-ggplot(cells_long[cells_long$Image=="Image75",],aes(x=time,y=Ratio,group=Cell,colour=Stimulus))
#  p+geom_path(alpha=0.2)+
#    #facet_wrap("Image")+
#    labs(x="Time (minutes)",y=labelsY)+
#    scale_x_continuous(breaks = pretty_breaks(n=10))+
#    theme(axis.title = element_text(family = "Arial", color="black", size=18))+
#    theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
#    theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
#    theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
#    theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#          panel.background = element_blank(), axis.line = element_line(colour = "black"))

cells<-cells[,!grepl("Image32",names(cells))]
cells<-cells[,!grepl("Image33",names(cells))]
cells<-cells[,!grepl("Image34",names(cells))]
cells<-cells[,!grepl("Image35",names(cells))]
cells[,"Image75_Cell91"]<-NULL
cells[,"Image76_Cell23"]<-NULL
cells[,"Image76_Cell33"]<-NULL
cells[,"Image76_Cell95"]<-NULL
cells[,"Image27_Cell35"]<-NULL
cells[,"Image27_Cell11"]<-NULL
cells[,"Image29_Cell2"]<-NULL
cells[,"Image30_Cell36"]<-NULL
cells[,"Image56_Cell6"]<-NULL
cells[,"Image56_Cell84"]<-NULL
cells[,"Image71_Cell1"]<-NULL 
cells[,"Image72_Cell45"]<-NULL 
cells[,"Image72_Cell80"]<-NULL 
cells[,"Image75_Cell46"]<-NULL


cells_long<-gather(data=cells, key="Cell",value="Ratio",-time)
cells_long$Image<-str_split(cells_long$Cell,"_",simplify = T)[,1]

cells_long$Stimulus<-sapply(1:nrow(cells_long),function(i){
  image<-cells_long$Image[i]
  time<-cells_long$time[i]
  tests[tests$Time==time,image]
})

cells_long<-cells_long[!is.na(cells_long$Ratio),]
cells_long$Genotype<-sapply(cells_long$Image,function(i) unique(animalTable$Genotype[animalTable$Image==i]))
cells_long$AnimalID<-sapply(cells_long$Image,function(i) unique(animalTable$AnimalID[animalTable$Image==i]))
# 15 min pre perifusion before recording, 
# 10min - 3mM Glucose  60 reps
# 30min - 15mM Glucose 180 reps
# 15min - 3mM Glucose 90 reps
# 5min - 30mM KCL 30 reps
# 10min - 3mM Glucose 60 reps

labelsY=expression(paste("Cytosolic Ca"^"2+", " Fura2 (340/380)",sep=""))

p<-ggplot(cells_long,aes(x=time,y=Ratio,group=Cell,colour=Stimulus))
p+geom_path(alpha=0.2)+
  facet_wrap("Image")+
  labs(x="Time (minutes)",y=labelsY)+
  scale_x_continuous(breaks = pretty_breaks(n=10))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=12))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=12))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("traces_images.png",path="./Figures",width = 30, height = 20, units = "cm")
ggsave("traces_images.svg",path="./Figures",device="svg",
       width = 30*0.393701, height = 20*0.393701,units="in")

nonstandard_images<-paste("Image",c(29,30,31,56,62,63,64,65,66),sep="")

cells_map<-pivot_wider(cells_long[!cells_long$Image%in%nonstandard_images,],
                       id_cols=c(time,Stimulus),names_from = Cell,values_from = Ratio)
cells_map<-cells_map[order(cells_map$time),]
cells_mat<-t(as.matrix(cells_map[,3:ncol(cells_map)]))
colnames(cells_mat)<-paste("time",1:429,sep="")

cells_mat_norm <- as.matrix(scale(cells_mat,center=T,scale=T))

ann_colours=list(
  Stimulus=c("baseline"="#E41A1C","15G"="#377EB8",
             "3G"="#4DAF4A","30KCl"="#984EA3"),
  Genotype=c("No Cre"="#000000",
             "HET"="#95251f",
             "KO"="#f4b9c2")
)

annotation_col<-data.frame(Stimulus=cells_map$Stimulus)
rownames(annotation_col)<-paste("time",1:429,sep="")
annotation_row<-data.frame(sapply(rownames(cells_mat),function(c) unique(cells_long$Genotype[cells_long$Cell==c])))
names(annotation_row)<-"Genotype"
annotation_row$Cell<-rownames(annotation_row)
annotation_row<-annotation_row[order(annotation_row$Genotype),]
cells_mat<-cells_mat[annotation_row$Cell,]
annotation_row$Cell<-NULL

png(filename = "./Figures/traces_heatmap.png",width=20,height=20,units="cm",res=300)
pheatmap(cells_mat,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         na_col = "#ffffff",
         cluster_cols = F,
         cluster_rows = F,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_legend = T,
         annotation_colors = ann_colours,
         show_rownames = F, show_colnames = F)
dev.off()

di<-dist(t(cells[,-3524]))
mds<-as.data.frame(cmdscale(di))

mds$Genotype<-sapply(row.names(mds),function(nm){
  image<-str_split(nm,"_",simplify = T)[,1]
  animalTable$Genotype[animalTable$Image==image]})

mds$AnimalID<-sapply(row.names(mds),function(nm){
  image<-str_split(nm,"_",simplify = T)[,1]
  animalTable$AnimalID[animalTable$Image==image]})

mds$Genotype<-factor(mds$Genotype)

p<-ggplot(mds,aes(x=V1,y=V2, colour=Genotype, shape=AnimalID)) #label=row.names(mds),
p+geom_point(alpha=0.8)+
  labs(x="Dimension 1",y="Dimension 2")+
  #scale_colour_manual(name="Cell Type",values=c("#999999", "#00FF00"))+
  theme(plot.title = element_text(family = "Open Sans Light", color="black",  size=32, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("MDS_raw.png",path="./Figures",width = 18, height = 10, units = "cm")

find_peaks <- function (x, lo, m = 3){ #count peaks based on local change in sign
  shape <- diff(sign(diff(x)))
  pks <- sapply(which(shape < 0), function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks <- sapply(pks, function(i){
    if(abs(lo-x[i])/((lo+x[i])/2)>0.2) return(i) else return(numeric(0)) #check if peak is 20% above median level
  })
  pks<-unlist(pks)
  pks
}

cell_nms<-names(cells)[-3524]

getProfile<-function(nm){
  image<-str_split(nm,"_",simplify = T)[,1]
  test<-tests[,image]
  cell<-cells[,nm]
  min<-min(cell,na.rm=T)
  max<-max(cell,na.rm=T)
  lo <- median(cell[test=="baseline"],na.rm=T)
  md <- median(cell[test=="15G"],na.rm=T)
  hi <- median(cell[test=="30KCl"],na.rm=T) 
  glu <- (md-lo)/(hi-lo) # glucose response
  kcl <- max(cell[test=="30KCl"],na.rm=T)/(max-lo)
  osc_lo <- mad(cell[test=="baseline"],na.rm=T)/(max-lo)
  osc_md <- mad(cell[test=="15G"],na.rm=T)/(max-lo)
  osc_hi <- mad(cell[test=="30KCl"],na.rm=T)/(max-lo)
  pks_lo <- length(find_peaks(cell[test=="baseline"],min))
  pks_md <- length(find_peaks(cell[test=="15G"],min))
  pks_hi <- length(find_peaks(cell[test=="30KCl"],min))
  c("GlucResp"=glu,"GlucOsc"=osc_md,"KClOsc"=osc_hi,
    "BaseOsc"=osc_lo,"BasePeaks"=pks_lo,
    "GlucPeaks"=pks_md,"KClPeaks"=pks_hi,"KClResp"=kcl)
}

features<-t(sapply(names(cells)[-3524],getProfile))
BIC = mclustBIC(features, G=1:10)

png(filename = "./Figures/BIC.png",width=20,height=20,units="cm",res=300)
plot(BIC)
dev.off()

mod1 = Mclust(features, x=BIC)
summary(mod1,parameters=T)
mod1dr = MclustDR(mod1)
ann_colours<-c("#239100",
               "#002ab3",
               "#f9e400",
               "#490036",
               "#b7f6ff",
               "#bb5800",
               "#00a283",
               "#edb1ff",
               "#917d00",
               "#003034")

summary(mod1dr)

png(filename = "./Figures/pairs.png",width=20,height=20,units="cm",res=300)
plot(mod1dr, what = "pairs",colors=ann_colours)
dev.off()

plot(mod1dr,what="scatterplot")

png(filename = "./Figures/boundaries.png",width=20,height=20,units="cm",res=300)
plot(mod1dr, what = "boundaries",ngrid=200,colors=ann_colours)
dev.off()

plot(mod1dr, what="classification")

featuresNorm <- as.data.frame(scale(features,center=T,scale=T))

pcaData <- prcomp(featuresNorm,center = T,scale.=T)
pcaData_out <- as.data.frame(pcaData$x)
pcaData_out$Genotype<-sapply(row.names(featuresNorm),function(nm){
  image<-str_split(nm,"_",simplify = T)[,1]
  animalTable$Genotype[animalTable$Image==image]})
pcaData_out$AnimalID<-sapply(row.names(featuresNorm),function(nm){
  image<-str_split(nm,"_",simplify = T)[,1]
  animalTable$AnimalID[animalTable$Image==image]})

pcaData_out$Cluster<-mod1dr$classification


percentage <- round(pcaData$sdev / sum(pcaData$sdev) * 100, 2)
percentage <- paste( colnames(pcaData_out), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(pcaData_out, aes(PC1, PC2, shape=Genotype,colour=Cluster)) +
  geom_point(size=3,alpha = 0.6,fill=NA) +
  #geom_text(size=3)+
  xlab(percentage[1]) + ylab(percentage[2])+
  scale_color_manual(name="Cluster",values = ann_colours)+
  #scale_shape_manual(values=c(21,24), name="Cell Type")+
  #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: \n",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(plot.title = element_text(family = "Open Sans Light", color="black",  size=32, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA_clusters_genotype.png",path="./Figures",width = 8, height = 10, units = "cm")

ggplot(pcaData_out, aes(PC1, PC2, shape=AnimalID,colour=Cluster)) +
  geom_point(size=3,alpha = 0.6,fill=NA) +
  #geom_text(size=3)+
  xlab(percentage[1]) + ylab(percentage[2])+
  scale_color_manual(name="Cluster",values = ann_colours)+
  #scale_shape_manual(values=c(21,24), name="Cell Type")+
  #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: \n",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(plot.title = element_text(family = "Open Sans Light", color="black",  size=32, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA_clusters_animal.png",path="./Figures",width = 8, height = 10, units = "cm")

df_out_r <- as.data.frame(pcaData$rotation)
df_out_r$feature <- row.names(df_out_r)

df_out_r

labs<-c("KClOsc"="KCl\noscillation","KClPeaks"="KCl\npeaks",
        "KClResp"="Response to\nKCl","BaseOsc"="3 mM\nglucose\noscillation",
        "BasePeaks"="3 mM\nglucose peaks","GlucOsc"="15 mM glucose\noscillation",
        "GlucPeaks"="15 mM glucose\npeaks","GlucResp"="Response to\n15 mM glucose")

p<-ggplot(df_out_r,aes(x=PC1,y=PC2,label=labs))
p+#geom_point(show.legend = F) +
  geom_text(size=3,show.legend = F)+
  xlab(percentage[1]) + ylab(percentage[2])+
  #coord_fixed() +
  theme(plot.title = element_text(family = "Open Sans Light", color="black",  size=32, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("PCA_features.png",path="./Figures",width = 20, height = 10, units = "cm")


#cells_long<-gather(data=cells, key="Cell",value="Ratio",-time)
#cells_long$Image<-str_split(cells_long$Cell,"_",simplify = T)[,1]

# cells_long$Stimulus<-sapply(1:nrow(cells_long),function(i){
#   image<-cells_long$Image[i]
#   time<-cells_long$time[i]
#   tests[tests$Time==time,image]
# })
# 
# cells_long<-cells_long[!is.na(cells_long$Ratio),]
cells_long$Cluster<-sapply(cells_long$Cell,function(c) unique(pcaData_out$Cluster[row.names(pcaData_out)==c]))
cells_long$AnimalID<-sapply(cells_long$Cell,function(c) unique(pcaData_out$AnimalID[row.names(pcaData_out)==c]))
cells_long$AnimalID<-factor(cells_long$AnimalID,levels=c("M10B1","M10C2","M10B2","M10C3","M10C4","M13B1"))

set.seed(1)
images<-sapply(unique(cells_long$Image),function(i){
  l<-length(unique(cells_long$Cell[cells_long$Image==i]))
  sels<-sample(1:l,10,replace=F)
  unique(cells_long$Cell[cells_long$Image==i])[sels]
})
images<-gather(as.data.frame(images),value="Image")
images<-images$Image

labs2<-c("baseline"="3 mM glucose","glucose"="15 mM glucose","kcl"="30 mM KCl")
labelsY<-expression(paste("Cytosolic Ca"^"2+", " Fura2 (340/380)",sep=""))
cells_long$Genotype<-factor(cells_long$Genotype, levels=c("No Cre","HET","KO"))
labsGen<-c("HET"=expression(paste(italic("Insr"^"-/wt"),italic("Ins1"^"cre/wt"),"nTnG"^"+/wt",sep="")),
           "KO"=expression(paste(italic("Insr"^"-/-"),italic("Ins1"^"cre/wt"),"nTnG"^"+/wt",sep="")),
           "No Cre"=expression(paste(italic("Insr"^"-/wt"),italic("Ins1"^"wt/wt"),"nTnG"^"+/wt",sep="")))

p<-ggplot(cells_long[cells_long$Cell%in%images,],
          aes(x=time,y=Ratio,group=Cell,colour=Genotype))
p+geom_path(alpha=0.2, size=0.5)+
  facet_wrap("AnimalID")+
  labs(x="Time (minutes)",y=labelsY)+
  scale_x_continuous(breaks = pretty_breaks(n=10))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  scale_colour_manual(name="Genotype",values=c("#000000",
                                               "#95251f",
                                               "#f4b9c2"),
                      labels=labsGen)+
  scale_fill_manual(name="Genotype",values=c("#000000",
                                             "#95251f",
                                             "#f4b9c2"),
                    labels=labsGen)+
  scale_shape_manual(values=c(19,15,17),labels=labsGen)+
  guides(colour = guide_legend(override.aes = list(size = 1,alpha=1)))+
  theme(strip.text = element_text(family="Arial",color="black",size=14),
        strip.background = element_blank())+
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14,hjust = 0))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("traces_genotype_all.svg",path="./Figures",width = 25, height = 15, units = "cm",dpi = 450)

p<-ggplot(filter(cells_long,Cell%in%images,AnimalID%in%c("M10C2","M10C3","M10C4")),
          aes(x=time,y=Ratio,group=Cell,colour=Genotype))
p+geom_path(alpha=0.2, size=0.5)+
  facet_wrap("AnimalID")+
  labs(x="Time (minutes)",y=labelsY)+
  scale_x_continuous(breaks = pretty_breaks(n=10))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  scale_colour_manual(name="Genotype",values=c("#000000",
                                               "#95251f",
                                               "#f4b9c2"),
                      labels=labsGen)+
  scale_fill_manual(name="Genotype",values=c("#000000",
                                             "#95251f",
                                             "#f4b9c2"),
                    labels=labsGen)+
  scale_shape_manual(values=c(19,15,17),labels=labsGen)+
  guides(colour = guide_legend(override.aes = list(size = 1,alpha=1)))+
  theme(strip.text = element_blank(),
        strip.background = element_blank())+
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14,hjust = 0))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("traces_genotype_sub.svg",path="./Figures",width = 25, height = 7.5, units = "cm",dpi = 450)

p<-ggplot(cells_long[cells_long$Cell%in%images,],
          aes(x=time,y=Ratio,group=Cell))+
  geom_path(alpha=0.2, size=0.5)+
  facet_wrap("AnimalID")+
  labs(x="Time (minutes)",y=labelsY)+
  scale_x_continuous(breaks = pretty_breaks(n=10))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  # scale_colour_manual(name="Genotype",values=c("#000000",
  #                                              "#95251f",
  #                                              "#f4b9c2"),
  #                     labels=labsGen)+
  # scale_fill_manual(name="Genotype",values=c("#000000",
  #                                            "#95251f",
  #                                            "#f4b9c2"),
  #                   labels=labsGen)+
  scale_shape_manual(values=c(19,15,17),labels=labsGen)+
  #guides(colour = guide_legend(override.aes = list(size = 1,alpha=1)))+
  theme(strip.text = element_text(family="Arial",color=NA,size=14))+
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_blank())+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#95251f","#f4b9c2","#f4b9c2","#000000",
           "#000000",
           "#95251f")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.newpage()
grid.draw(g)

ggsave("traces_genotype_rect.png",path="./Figures",width = 25, height = 15, units = "cm",dpi = 450)


#labelsCells<-c(expression("EGFP"^"-"),expression("EGFP"^"+"))

p<-ggplot(cells_long,aes(x=time,y=Ratio,colour=Genotype,shape=Stimulus,group=Cell))+
  geom_line(alpha=0.6)+
  geom_point(aes(colour=Genotype),size=1,fill=NA,alpha=0.6)+
  facet_grid(Cluster~.)+
  labs(x="Time (minutes)",y=labelsY)+
  scale_colour_manual(name="Genotype",values=c("#000000",
                                               "#95251f",
                                               "#f4b9c2"),
                      labels=labsGen)+
  scale_fill_manual(name="Genotype",values=c("#000000",
                                             "#95251f",
                                             "#f4b9c2"),
                    labels=labsGen)+
  scale_shape_manual(values=c(19,15,17),labels=labsGen)+
  theme(strip.text = element_text(family="Arial",color="black",size=14))+
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14),
        legend.title = element_text(family = "Arial",colour="black",size=14))+
  #guides(shape = guide_legend(override.aes = list(size = 3)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

z <- ggplotGrob(p)
z <- gtable_add_padding(z, unit(0.5, "cm"))
z <- gtable_add_cols(z, unit(0.8, 'cm'), 8)
z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, fill = gray(0.85)),width=unit(0.9,"npc")),
                          textGrob("Cluster", rot=-90,gp = gpar(col = gray(0),
                                                                fontfamily="Arial",fontsize=14))),
                     6, 9, 26,9, name = paste(runif(2))) #4n+1
z <- gtable_add_cols(z, unit(1/8, "line"), 7)
grid.newpage()
grid.draw(z)

ggsave("traces_cluster.png",plot=grid.draw(z),
       path="./Figures",width = 20, height = 20, units = "cm")


# p<-ggplot(cells_long,aes(x=Time,y=Ratio,colour=Stimulus,shape=Cluster,group=Cell))+
#   geom_line(alpha=0.6)+
#   geom_point(aes(colour=Exposure),size=1.5,fill=NA,alpha=0.6)+
#   facet_wrap(~Name)+
#   labs(x="Time (minutes)",y=labelsY)+
#   scale_x_continuous(breaks = pretty_breaks(n=3))+
#   scale_colour_manual(values=brewer.pal(5,"Set2"),
#                       labels=c("0.5 mM glucose","11 mM glucose","Adrenaline"))+
#   #scale_shape_manual(values=c(16,3,17),labels=labs2)+
#   theme(strip.text = element_text(family="Arial",color="black",size=12))+
#   theme(axis.title = element_text(family = "Arial", color="black", size=12))+
#   theme(axis.text.x = element_text(family = "Arial",colour="black",size=12))+
#   theme(axis.text.y = element_text(family = "Arial",colour="black",size=12))+
#   theme(legend.text = element_text(family = "Arial",colour="black",size=12))+
#   theme(legend.title = element_text(family = "Arial",colour="black",size=12))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         panel.spacing = unit(0.25,"cm"))
# p
# 
# ggsave("traces_indiv.png",
#        path="./Figures",width = 25, height = 20, units = "cm")

featuresNorm2<-as.matrix(featuresNorm)
rc <- rainbow(nrow(featuresNorm2), start = 0, end = .3)
cc <- rainbow(ncol(featuresNorm2), start = 0, end = .3)
rowCols<-sapply(row.names(featuresNorm2),
                function(c) unique(cells_long$Genotype[cells_long$Cell==c]))
rowCols<-case_when(rowCols=="No Cre" ~ "#fb523a",
                   rowCols=="KO" ~ "#9e0093",
                   TRUE ~ "#f1a9c4")


png(filename = "./Figures/heatmap.png",width=20,height=20,units="cm",res=300)
heatmap(featuresNorm2,Rowv=NA,Colv=NA,col = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(20),
        RowSideColors = rowCols,
        margins=c(10,3),labCol=labs)
dev.off()

features_full<-as.data.frame(features)
features_full$Cluster<-pcaData_out$Cluster
features_full$Cell<-row.names(features_full)

write.table(features_full,
            file="features.txt",col.names=T,row.names=F,quote=F,sep="\t")

results<-data.frame(cell=row.names(pcaData_out),cluster=pcaData_out$Cluster)
results$Genotype<-sapply(results$cell,function(c) unique(pcaData_out$Genotype[row.names(pcaData_out)==c]))

cells_long$Factor<-interaction(cells_long$Stimulus,cells_long$Genotype,cells_long$Cluster)

resultsTab<-cells_long%>%group_by(Factor)%>%summarize(y0=min(Ratio,na.rm = T),y75=quantile(Ratio,0.75,na.rm=T), 
                                                      y25=quantile(Ratio,0.25,na.rm=T),
                                                      y50=median(Ratio,na.rm=T),y100=max(Ratio,na.rm = T))

resultsTab$Cluster<-sapply(resultsTab$Factor,function(f) unique(cells_long$Cluster[cells_long$Factor==f]))
resultsTab$Cluster<-factor(resultsTab$Cluster,levels=levels(cells_long$Cluster))
resultsTab$Genotype<-sapply(resultsTab$Factor,function(f) unique(cells_long$Genotype[cells_long$Factor==f]))
#resultsTab$cell_type<-factor(resultsTab$cell_type,levels=levels(cells_long$cell_type))
resultsTab$Stimulus<-sapply(resultsTab$Factor,function(f) unique(cells_long$Stimulus[cells_long$Factor==f]))
resultsTab$Exposure<-factor(resultsTab$Exposure,levels=levels(cells_long$Exposure))

labels<-c("baseline"="3 mM glucose","15G"="15 mM glucose","3G"="3 mM glucose",
          "30KCl"="30 mM KCl")
# 
# p<-ggplot(cells_long,aes(x=Genotype,colour=Genotype,fill=Genotype,group=Factor))
# p<-p+geom_boxplot(aes(x=Genotype,middle=y50,ymin = y0,lower=y25,upper=y75,ymax = y100, fill = Genotype),
#                   stat="identity",data=resultsTab,alpha=0.6)+
#   labs(x="",y=labelsY)+
#   geom_point(aes(y=Ratio),size=1,fill=NA,alpha=0.6,
#              position=position_dodge2(width=0.2,preserve="total"))+
#   facet_grid(Cluster~Stimulus,labeller=labeller(Stimulus=labels))+
#   scale_colour_manual(name="Genotype",values=c("#000000",
#                                                "#95251f",
#                                                "#f4b9c2"))+
#   scale_fill_manual(name="Genotype",values=c("#000000",
#                                              "#95251f",
#                                              "#f4b9c2"))+
#   scale_shape_manual(values=c(19,15,17),labels=labs2)+
#   theme(strip.text = element_text(family="Arial",color="black",size=14))+
#   theme(axis.title = element_text(family = "Arial", color="black", size=14))+
#   theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
#   theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
#   theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
#   theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# z <- ggplotGrob(p)
# z <- gtable_add_padding(z, unit(0.5, "cm"))
# z <- gtable_add_cols(z, unit(0.8, 'cm'), 12)
# z <- gtable_add_grob(z, 
#                      list(rectGrob(gp = gpar(col = NA, fill = gray(0.85)),width=unit(0.9,"npc")),
#                           textGrob("Cluster", rot=-90,gp = gpar(col = gray(0),
#                                                                 fontfamily="Arial",fontsize=14))),
#                      9, 13, 22,13, name = paste(runif(2))) #4n+1
# z <- gtable_add_rows(z, unit(0.8, 'cm'), 2)
# z <- gtable_add_grob(z, 
#                      list(rectGrob(gp = gpar(col = NA, fill = gray(0.85)),height=unit(2.5,"npc")),
#                           textGrob("Exposure", gp = gpar(col = gray(0),
#                                                          fontfamily="Arial",fontsize=14))),
#                      3, 6, 3,12, name = paste(runif(2)))
# 
# z <- gtable_add_cols(z, unit(2/8, "line"), 7)
# z <- gtable_add_rows(z, unit(2/8, "line"), 3)
# grid.newpage()
# grid.draw(z)
# 
# ggsave("summary_stats.png",plot=grid.draw(z),dpi=320,
#        path="./Figures",width = 20, height = 15, units = "cm")

features_df<-as.data.frame(features)
features_df$Cluster<-sapply(row.names(features_df),function(c) unique(cells_long$Cluster[cells_long$Cell==c]))
#features_df$Cluster<-factor(features_df$Cluster,levels=levels(cells_long$Cluster))
features_df$Genotype<-sapply(row.names(features_df),function(c) unique(cells_long$Genotype[cells_long$Cell==c]))
features_df$Cell<-row.names(features_df)

features_long<-gather(data=features_df, key="Feature",value="Value",-Genotype,-Cluster,-Cell)

Feature.labs <- c("Baseline Oscillation", "Number of Peaks\n(Baseline)",
               "High Glucose Oscillation", "Number of Peaks\n(High Glucose)",
               "Response to Glucose",
               "KCl Oscillation", "Number of Peaks\n(KCl)",
               "Response to KCl")
names(Feature.labs) <- c("BaseOsc","BasePeaks","GlucOsc","GlucPeaks","GlucResp",
                      "KClOsc","KClPeaks","KClResp")

p<-ggplot(features_long,aes(y=Value,x=Genotype,colour=Genotype,fill=Genotype))
p+geom_boxplot(alpha=0.6)+
  labs(x="",y="Value")+
  #geom_point(aes(y=Value),size=1,fill=NA,alpha=0.6,
  #           position=position_dodge2(width=0.2,preserve="total"))+
  facet_wrap(~Feature,scales = "free",  labeller = labeller(Feature = Feature.labs))+
  scale_colour_manual(name="Genotype",values=c("#000000",
                                               "#95251f",
                                               "#f4b9c2"),
                      labels=labsGen)+
  scale_fill_manual(name="Genotype",values=c("#000000",
                                             "#95251f",
                                             "#f4b9c2"),
                    labels=labsGen)+
  scale_shape_manual(values=c(19,15,17),labels=labsGen)+
  theme(strip.text = element_text(family="Arial",color="black",size=14),
        strip.background = element_blank())+
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14,hjust = 0),
        legend.title = element_text(family = "Arial",colour="black",size=14),
        legend.position = c(0.8,0.15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("summary_stats_all.svg",path="./Figures",width = 20, height = 20, units = "cm")

mod1<-lm(Value~Genotype:Feature,data=features_long)
an<-aov(mod1)
dat<-TukeyHSD(an)

my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)

write.csv(dat$`Genotype:Feature`,file = "./multcomps.csv")

my_comparisons <- list( c("No Cre", "HET"), c("HET", "KO"))


anno<-dat$`Genotype:Feature`[c("KO:GlucPeaks-HET:GlucPeaks",
                         "No Cre:GlucPeaks-HET:BasePeaks"),"p adj"]

# Make plot with custom x and y position of the bracket
ggplot(iris, aes(x = Species, y = Sepal.Width, fill = Petal.Width > 1)) +
  geom_boxplot(position = "dodge") +
  geom_signif(
    annotation = formatC(anno, digits = 1),
    y_position = 4.05, xmin = 2.2, xmax = 3,
    tip_length = c(0.2, 0.04)
  )

p<-ggplot(filter(features_long,Feature=="GlucPeaks"),
          aes(y=Value,x=Genotype,colour=Genotype,fill=Genotype))
p+geom_boxplot(alpha=0.6,show.legend = F)+
  labs(x="",y="Number of Peaks")+
  geom_signif(
    comparisons = my_comparisons,
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = 0.08,
    step_increase = 0.05,
    tip_length = 0.01
  )+
  #stat_compare_means(comparisons = my_comparisons)+
  #geom_point(aes(y=Value),size=1,fill=NA,alpha=0.6,
  #           position=position_dodge2(width=0.2,preserve="total"))+
  #facet_wrap(~Feature,scales = "free",  labeller = labeller(Feature = Feature.labs))+
  scale_colour_manual(name="Genotype",values=c("#000000",
                                               "#95251f",
                                               "#f4b9c2"),
                      labels=labsGen)+
  scale_fill_manual(name="Genotype",values=c("#000000",
                                             "#95251f",
                                             "#f4b9c2"),
                    labels=labsGen)+
  scale_shape_manual(values=c(19,15,17),labels=labsGen)+
  theme(strip.text = element_text(family="Arial",color="black",size=14),
        strip.background = element_blank())+
  theme(axis.title = element_text(family = "Arial", color="black", size=14))+
  theme(axis.text.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_blank(),
        legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("summary_stats_peaks.svg",path="./Figures",width = 7.5, height = 7.5, units = "cm")
