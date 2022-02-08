physeq_aitchord=function(physeq) {
  
physeq.clrtransfo=logratio.transfo(as.matrix(t(otu_table(physeq))),logratio='CLR',offset=0.00001)
physeq.clrtransfoMAT=matrix(physeq.clrtransfo,nrow=nrow(physeq.clrtransfo),ncol=ncol(physeq.clrtransfo))
colnames(physeq.clrtransfoMAT)=colnames(physeq.clrtransfo)
rownames(physeq.clrtransfoMAT)=rownames(physeq.clrtransfo)

dist=vegdist(physeq.clrtransfoMAT,method="euclidean")
pcoa.euc=pcoa(dist)

}