physeq_aitchdiv=function(physeq,group) {
  
physeq.clrtransfo=logratio.transfo(as.matrix(t(otu_table(physeq))),logratio='CLR',offset=0.00001)
physeq.clrtransfoMAT=matrix(physeq.clrtransfo,nrow=nrow(physeq.clrtransfo),ncol=ncol(physeq.clrtransfo))
colnames(physeq.clrtransfoMAT)=colnames(physeq.clrtransfo)
rownames(physeq.clrtransfoMAT)=rownames(physeq.clrtransfo)

dist=vegdist(physeq.clrtransfoMAT,method="euclidean")
pcoa.euc=pcoa(dist)

physeq.clrtransfoMAT.mtd=as.data.frame(as(sample_data(physeq),"matrix"))

if (all.equal(rownames(pcoa.euc$vectors),rownames(physeq.clrtransfoMAT.mtd))) {
  mtd=physeq.clrtransfoMAT.mtd
  
  mtd$px_tn=paste0(mtd$RECIST_response_RRNR,"-",mtd$Patient_number,"_",mtd$Timepoint_bin)
  mtd$px_tn=gsub("baseline","t0",mtd$px_tn)
  mtd$px_tn=gsub("timepoint","tn",mtd$px_tn)

  #px_tn=as.data.frame(mtd$px_tn)
  #rownames(px_tn)=rownames(mtd)
  #colnames(px_tn)[1]="px_tn"

  vctrs=as.data.frame(pcoa.euc$vectors)
  
  vctrs=cbind(vctrs,mtd)
  vctrs$ps.rownames=rownames(mtd)
}
f=as.formula(paste0("cbind(Axis.1,Axis.2)~",group))
centroidsR=aggregate(f,vctrs,mean)
names(centroidsR)=c(group,"Axis.1.centroid","Axis.2.centroid")
vctrs=left_join(vctrs,centroidsR,group)

}