library(usedist)
library(metagMisc)

vctrs_t0dist=function(vctrs,px.grp,tbin.grp,t0.var,tn.var,tn.grp,trt.grp) {

vctrs2=vctrs

vctrs2$Subject_tn=make.unique(paste0(vctrs2$Subject_number,"_",vctrs2$Timepoint_number))

vctrs2$Axis.1.wtd=vctrs2$Axis.1*pcoa$values$Relative_eig[1]
vctrs2$Axis.2.wtd=vctrs2$Axis.2*pcoa$values$Relative_eig[2]

vctrs2=vctrs2 %>% group_by(Subject_number)


coords.df=vctrs2 %>% ungroup() %>% select("Axis.1.wtd","Axis.2.wtd")
coords.df=as.data.frame(coords.df)
rownames(coords.df)=vctrs2$Subject_tn

eucdist=dist(coords.df,"euclidean")

t0=which(vctrs2[,tbin.grp]==t0.var)

list_t0=list()
list_tn=list()

for (i in seq_along(t0)) {
  
  px=vctrs2$Subject_number[t0[i]]
  
  list_t0[[i]]=which(vctrs2$Subject_number %in% px & vctrs2[,tbin.grp]=="Baseline")
  list_tn[[i]]=which(vctrs2$Subject_number %in% px & vctrs2[,tbin.grp]=="Timepoint")
  
}

dist.bypx=list()

for (i in seq_along(t0)) {
  
  dist.bypx[[i]]=dist_get(eucdist,vctrs2$Subject_tn[list_t0[[i]]],vctrs2$Subject_tn[list_tn[[i]]])
  
}

dist.df=data.frame(dist.t0tn=reshape2::melt(dist.bypx)$value,Subject_tn=vctrs2$Subject_tn[reshape2::melt(list_tn)$value],Subject_number=vctrs2$Subject_number[reshape2::melt(list_tn)$value],Timepoint_number=vctrs2[reshape2::melt(list_tn)$value,tn.grp],Group=vctrs2$Group[reshape2::melt(list_tn)$value])

dist.allt0=dist_subset(eucdist,vctrs2$Subject_tn[reshape2::melt(list_t0)$value])
dist.allt0.df=dist2list(dist.allt0,tri=TRUE)

#dist.allt0.df=dist.allt0.df[duplicated(dist.allt0.df$value)==FALSE,]
dist.df=rbind(dist.df,data.frame(dist.t0tn=dist.allt0.df$value,Subject_tn=dist.allt0.df$row,Subject_number=paste0(gsub("_.*","",dist.allt0.df$row),"_",gsub("_.*","",dist.allt0.df$col)),Timepoint_number=rep("0",rep=nrow(dist.allt0.df)),Group=rep("within-t0",rep=nrow(dist.allt0.df))))

dist.df
}

# ggplot(vctrs2[!is.na(vctrs2$dist.t0tn),],aes(Cycle_agg2,dist.t0tn,color=Cycle_agg2)) + geom_boxplot(outlier.alpha=0,width=0.2,alpha=0.60) + geom_jitter() + stat_compare_means(comparisons=comp) + scale_color_viridis_d(option="H",begin=0.1,end=0.9)


