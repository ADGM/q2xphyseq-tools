deseq_annotRes=function(dds,var="Tissue_type",alpha=0.05) {
  
contrast=resultsNames(dds)[2]

g0=levels(colData(dds)[,var])[1]
g1=levels(colData(dds)[,var])[2]

res = results(dds, alpha=alpha,name=contrast)
#numerator #then denominator

#order hits by p-adj value
res = res[order(res$padj, na.last=NA), ]

#keep only significant hits in a df
sigtab = res[(res$padj < alpha), ]

if (nrow(sigtab)==0) {

  print("No significant hits. Return DESeq2 results table without annotation.")

sigtab=as.data.frame(sigtab)

print(datatable(sigtab, caption="DESeq2 results"))

}

else if (nrow(sigtab)>0) {

sigtab$value=rownames(sigtab)

#retrieve normalized counts fr deseq
normalized_counts=counts(dds, normalized=TRUE)

#subset normalized counts table to signifcant hits only
normalized_counts.toptax=normalized_counts[sigtab$value,]
meannormcount=rowMeans(normalized_counts.toptax)
sigtabxannot$mean.normcount=meannormcount

#get list of CRC samples by group from phyloseq metadata
mtd=sample_data(physeq)
g1spls=rownames(mtd[mtd[[var]]==g1,])
g0spls=rownames(mtd[mtd[[var]]==g0,])

#match CRC-samples-by-group with normalized counts table
normcts.g1=normalized_counts.toptax[,g1spls]
normcts.g1=as.data.frame(normcts.g1)
rownames(normcts.g1)=rownames(normalized_counts.toptax)

normcts.g0=normalized_counts.toptax[,g0spls]
normcts.g0=as.data.frame(normcts.g0)
rownames(normcts.g0)=rownames(normalized_counts.toptax)

#compute prevalence of each taxon by group (# of samples where taxon is >0 / total # of samples in the group)
g1=adply(normcts.g1, 1, function(x) {(nnzero(x)/ncol(normcts.g1))*100})
g1=g1$V1
g0=adply(normcts.g0, 1, function(x) {(nnzero(x)/ncol(normcts.g0))*100})
g0=g0$V1

#integrate prevalence info in grp
sigtabxannot$grp1=g1
sigtabxannot$grp0=g0

#reformat table to keep one column for prevalence group
sigtabxannot.melt=pivot_longer(sigtabxannot,starts_with("grp"),names_to="Prevalence_grp",values_to="Prevalence_percent")

#arrange df by prevalence group
sigtabxannot.melt=sigtabxannot.melt %>% arrange(Prevalence_grp) 

#compute mean norm counts of each taxon by group 
mean.grp1=as.data.frame(rowMeans(normcts.g1))
colnames(mean.grp1)="Mean_normcount"
mean.grp1$value=rownames(mean.grp1)
mean.grp1$Prevalence_grp="grp1"

mean.grp0=as.data.frame(rowMeans(normcts.g0))
colnames(mean.grp0)="Mean_normcount"
mean.grp0$value=rownames(mean.grp0)
mean.grp0$Prevalence_grp="grp0"

#make a column for mean normcount
Meanct=rbind(mean.grp0,mean.grp1)

if (all.equal(sigtabxannot.melt$value,Meanct$value)) {
  if (all.equal(sigtabxannot.melt$Prevalence_grp,Meanct$Prevalence_grp)) {
  sigtabxannot.melt=cbind(sigtabxannot.melt,Mean_normcount=Meanct$Mean_normcount)
  }
}  

#if there is need to scale norm counts
#sigtabxannot.melt$mean.normcount=sigtabxannot.melt$mean.normcount/2

#factor prevalence group to have grp1 first before grp0
sigtabxannot.melt$Prevalence_grp=factor(sigtabxannot.melt$Prevalence_grp,levels=c("grp1","grp0"))


sigtabxannot.melt.reord=sigtabxannot.melt[,c((7:ncol(sigtabxannot.melt)),(1:6))]

print(datatable(sigtabxannot.melt.reord,caption="DESeq results with per-group taxon stats"))
sigtabxannot.melt.reord

}

}
