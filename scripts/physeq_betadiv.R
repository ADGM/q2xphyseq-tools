physeq_betadiv=function(physeq,method,group) {

ps.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=min(sample_sums(physeq)), replace=F)
#factor variables before declaring new variable for collating alphadiv vals

sd=min(sample_sums(physeq))
print_sd=paste0("Rarefied at ",sd," depth")

print(print_sd)

dist = phyloseq::distance(ps.rarefied, method=method)
ordination = ordinate(ps.rarefied, method="PCoA", distance=dist)

if (all.equal(rownames(ordination$vectors),rownames(sample_data(ps.rarefied)))) {
  vctrs=cbind(ordination$vectors,sample_data(ps.rarefied))
  vctrs$psrar.rownames=rownames(sample_data(ps.rarefied))
}

f=as.formula(paste0("cbind(Axis.1,Axis.2)~",group))
centroidsR=aggregate(f,vctrs,mean)
names(centroidsR)=c(group,"Axis.1.centroid","Axis.2.centroid")
vctrs=left_join(vctrs,centroidsR,group)

}