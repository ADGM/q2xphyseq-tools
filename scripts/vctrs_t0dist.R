library(usedist)
library(metagMisc)

vctrs_t0dist=function(vctrs) {

vctrs2=vctrs
vctrs2$px_tn=paste0(vctrs2$Patient_number,"_",vctrs2$Cycle_agg2)
vctrs2$px_cyc=paste0(vctrs2$Patient_number,"_",unfactor(vctrs2$Cycle_number))

vctrs2$Axis.1.wtd=vctrs2$Axis.1*pcoa.euc$values$Relative_eig[1]
vctrs2$Axis.2.wtd=vctrs2$Axis.2*pcoa.euc$values$Relative_eig[2]

vctrs2=vctrs2 %>% group_by(Patient_number)

coords.df=vctrs2 %>% ungroup() %>% select("Axis.1.wtd","Axis.2.wtd")
coords.df=as.data.frame(coords.df)
rownames(coords.df)=vctrs2$px_cyc

eucdist=dist(coords.df,"euclidean")

t0=which(vctrs2$Timepoint_bin=="Baseline")

list_t0=list()
list_tn=list()

for (i in seq_along(t0)) {
  
  px=vctrs2$Patient_number[t0[i]]
  
  list_t0[[i]]=which(vctrs2$Patient_number %in% px & vctrs2$Timepoint_bin=="Baseline")
  list_tn[[i]]=which(vctrs2$Patient_number %in% px & vctrs2$Timepoint_bin=="Timepoint")
  
}

dist.bypx=list()

for (i in seq_along(t0)) {
  
  dist.bypx[[i]]=dist_get(eucdist,vctrs2$px_cyc[list_t0[[i]]],vctrs2$px_cyc[list_tn[[i]]])
  
}

dist.df=data.frame(dist.t0tn=reshape2::melt(dist.bypx)$value,px_cyc=vctrs2$px_cyc[reshape2::melt(list_tn)$value],Patient_number=vctrs2$Patient_number[reshape2::melt(list_tn)$value],Cycle_number=vctrs2$Cycle_number[reshape2::melt(list_tn)$value],Cycle_agg2=vctrs2$Cycle_agg2[reshape2::melt(list_tn)$value])

dist.allt0=dist_subset(eucdist,vctrs2$px_cyc[reshape2::melt(list_t0)$value])
dist.allt0.df=dist2list(dist.allt0,tri=TRUE)

#dist.allt0.df=dist.allt0.df[duplicated(dist.allt0.df$value)==FALSE,]
dist.df=rbind(dist.df,data.frame(dist.t0tn=dist.allt0.df$value,px_cyc=dist.allt0.df$row,Patient_number=paste0(gsub("_.*","",dist.allt0.df$row),"_",gsub("_.*","",dist.allt0.df$col)),Cycle_number=rep("0",rep=nrow(dist.allt0.df)),Cycle_agg2=rep("0",nrow(dist.allt0.df))))


dist.df$px_tn=paste0(dist.df$Patient_number,"_",vctrs2$Cycle_number)

dist.df$Cycle_agg2=factor(dist.df$Cycle_agg2,levels=c("0","1-11","12-17","18-49"))
comp=list(c("0","1-11"),c("1-11","12-17"),c("12-17","18-49"))
comp2=list(c("0","1-11"),c("0","12-17"),c("0","18-49"))

dist.df$Cycle_number=as.numeric(unfactor(dist.df$Cycle_number))

dist.df
}

# ggplot(vctrs2[!is.na(vctrs2$dist.t0tn),],aes(Cycle_agg2,dist.t0tn,color=Cycle_agg2)) + geom_boxplot(outlier.alpha=0,width=0.2,alpha=0.60) + geom_jitter() + stat_compare_means(comparisons=comp) + scale_color_viridis_d(option="H",begin=0.1,end=0.9)


