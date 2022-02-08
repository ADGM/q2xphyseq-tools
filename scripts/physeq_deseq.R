physeq_deseq=function(physeq,formula,ref.var,ref.grp,comp.grp) {
    
#physeq=readRDS("physeq_838_mtd2109xx.rds")


#design.formula=("~ Timepoint_bin + Pathological_response_bin + Seq_batch")
#ref.var="Timepoint_bin"
#ref.grp="Baseline"
#comp.grp="Timepoint"

#spls are cols for counts table for deseq
cts=as(otu_table(physeq),"matrix")
mtd=as(sample_data(physeq),"matrix")


if (all.equal(colnames(cts),rownames(mtd))) {


dds = DESeqDataSetFromMatrix(countData = cts,
                              colData = mtd,
                              design = paste0(formula))

}

dds[[ref.var]]=relevel(dds[[ref.var]], ref = ref.grp)
#ref is denominator

dds = DESeq(dds, test="Wald", fitType="local", sfType="poscounts")


#saveRDS(dds,"16S-DESeq-BvTall_res.rds")


}

#physeq=readRDS("physeq_838_mtd2109xx.rds")

#subset samples
#physeq=subset_samples(physeq, Pathological_response_bin!="NA")

#design.formula=("~ Timepoint_bin + Pathological_response_bin + Seq_batch")
#ref.var="Timepoint_bin"
#ref.grp="Baseline"
#comp.grp="Timepoint"

#dds=readRDS("16S-DESeq-BvTall_res.rds")

#set alpha
#alpha = 0.05
#datatable(as.data.frame(alpha),caption="Cut-off level of significance (alpha)")

#set contrast
#contrast=resultsNames(dds)[3]
#datatable(as.data.frame(contrast),caption="Groups contrasted in DESeq results (grp1 vs grp0)")


#res = results(dds, alpha=alpha,name=contrast)
#numerator #then denominator

#order hits by p-adj value
#res = res[order(res$padj, na.last=NA), ]

#keep only significant hits in a df
#sigtab = res[(res$padj < alpha), ]
#sigtab$value=rownames(sigtab)

#get taxonomy table from physeq object 
#taxnames.physeq=as(tax_table(physeq),"matrix")
#taxnames.physeq=cbind(rownames(tax_table(physeq)),taxnames.physeq)
#taxnames.physeq=as_tibble(taxnames.physeq)
#names(taxnames.physeq)[1]="value"

#match ASVs in sig hits table with taxonomy table
#sigtabdf=as_tibble(sigtab)
#sigtabxannot=left_join(sigtabdf,taxnames.physeq)

#find lowest identifiable taxon for each hit
#sigtabxannot$Genussp=paste0(sigtabxannot$Genus," ",sigtabxannot$Species)

#collapse ambiguous taxa names at Genus-level to Family taxa
#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"NA NA"))]=sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"NA NA"))]

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^uncultured NA"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^uncultured NA"))]," ","uncultured")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^uncultured gut metagenome"))]=sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^uncultured gut metagenome"))]

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^uncultured bacterium"))]=sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^uncultured bacterium"))]

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^metagenome"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^metagenome"))]," metagenome")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^uncultured organism"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^uncultured organism"))]," uncultured")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^uncultured organism"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^uncultured organism"))]," uncultured")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^CAG-56"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^CAG-56"))]," CAG-56")
#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^UBA1819"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^UBA1819"))]," UBA1819")

#collapse NA taxa names at Genus-level to nearest taxonomic level without NA
#sigtabxannot$Genussp[is.na(sigtabxannot$Genussp)]=sigtabxannot$Order[is.na(sigtabxannot$Genussp)]
#sigtabxannot$Genussp[is.na(sigtabxannot$Genussp)]=sigtabxannot$Class[is.na(sigtabxannot$Genussp)]
#sigtabxannot$Genussp[is.na(sigtabxannot$Genussp)]=sigtabxannot$Phylum[is.na(sigtabxannot$Genussp)]
#sigtabxannot$Genussp[is.na(sigtabxannot$Genussp)]=sigtabxannot$Kingdom[is.na(sigtabxannot$Genussp)]

#remove unnecessary NAs and symbols
#sigtabxannot$Genussp=gsub(" NA$","",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("\\[","",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("\\]","",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("= Timone 84634 = DSM 17679 = JCM 13223","",sigtabxannot$Genussp)

#remove redundants
#sigtabxannot$Genussp=gsub(" uncultured bacterium"," uncultured",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Lachnoclostridium uncultured Firmicutes bacterium","Lachnoclostridium uncultured bacterium",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Bacillus Burkholderia cepacia","Burkholderia cepacia",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Ruminiclostridium 9 uncultured Flavonifractor sp.","Flavonifractor sp.",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Blautia Ruminococcus sp. Marseille-P328","Ruminococcus sp. Marseille-P328",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Alistipes uncultured Alistipes sp.","Alistipes uncultured",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Parabacteroides Parabacteroides ","Parabacteroides ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Bacteroides Bacteroides ","Bacteroides ",sigtabxannot$Genussp)
##sigtabxannot$Genussp=gsub("Enterobacter Enterobacter ","Enterobacter ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Alistipes Alistipes ","Alistipes ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Streptococcus Streptococcus ","Streptococcus ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub(" Selenomonas sp.","",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Sutturella  ","",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Campylobacter Campylobacter ","Campylobacter ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Lachnospiraceae FCS020 group metagenome","Lachnospiraceae FCS020 group",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Peptoniphilus Peptoniphilus sp.","Peptoniphilus sp.",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Parabacteroides uncultured Bacteroidaceae bacterium","Parabacteroides uncultured",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Bacteroides Bacteroidaceae bacterium DJF_B220","Bacteroides bacterium DJF_B220",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Lactobacillus Lactobacillus ","Lactobacillus ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Alistipes Faecalibacterium ","Faecalibacterium ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Bifidobacterium Bifidobacterium ","Bifidobacterium ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub("Megasphaera Megasphaera ","Megasphaera ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub( " Bacteroidetes bacterium","",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub( "Candidatus Nucleicultrix Candidatus Nucleicultrix ","Candidatus Nucleicultrix ",sigtabxannot$Genussp)
#sigtabxannot$Genussp=gsub( "Erysipelatoclostridium Massiliomicrobiota ","Massiliomicrobiota ",sigtabxannot$Genussp)

#remove typos
#sigtabxannot$Genussp=gsub( "Sutterella Sutturella ","Sutturella ",sigtabxannot$Genussp)

#collapse ambiguous Family-level taxa names to Order
#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^S5-A14a uncultured organism"))]=paste0(sigtabxannot$Order[which(str_detect(sigtabxannot$Genussp,"^S5-A14a uncultured organism"))]," S5-A14a uncultured")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^Family XIII AD3011 group"))]=paste0(sigtabxannot$Order[which(str_detect(sigtabxannot$Genussp,"^Family XIII AD3011 group"))]," Family XIII AD3011 group")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^Family XIII AD3011 group gut metagenome"))]=paste0(sigtabxannot$Order[which(str_detect(sigtabxannot$Genussp,"^Family XIII AD3011 group gut metagenome"))]," Family XIII AD3011 group")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^TM7 phylum sp. oral clone FR058"))]=paste0(sigtabxannot$Order[which(str_detect(sigtabxannot$Genussp,"^TM7 phylum sp. oral clone FR058"))]," TM7 phylum oral clone FR058")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^uncultured bacterium"))]=paste0(sigtabxannot$Order[which(str_detect(sigtabxannot$Genussp,"^uncultured bacterium"))]," uncultured bacterium")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^GCA-900066225"))]=paste0(sigtabxannot$Family[which(str_detect(sigtabxannot$Genussp,"^GCA-900066225"))]," GCA-900066225")

#sigtabxannot$Genussp[which(str_detect(sigtabxannot$Genussp,"^gut metagenome"))]=paste0(sigtabxannot$Order[which(str_detect(sigtabxannot$Genussp,"^gut metagenome"))]," gut metagenome")

#filter out irrelevant taxa
#sigtabxannot=sigtabxannot[sigtabxannot$Genussp!="Bacteria",]
#sigtabxannot=sigtabxannot[sigtabxannot$Kingdom!="Unassigned",]
#sigtabxannot=sigtabxannot[sigtabxannot$Kingdom!="Archaea",]
#sigtabxannot=sigtabxannot[sigtabxannot$Genussp!="Mitochondria",]
#sigtabxannot=sigtabxannot[sigtabxannot$Genussp!="Chloroplast",]

#make name unique if duplicated at lowest-identifiable-taxonomic level 
##(instead of aggregating at this level; ASVs *may* be genomically the same or different but aggregating genomically different ASVs will make counts wrong)

#sigtabxannot$Genussp=make.unique(sigtabxannot$Genussp)

#retrieve normalized counts fr deseq
#normalized_counts=counts(dds, normalized=TRUE)

#subset normalized counts table to signifcant hits only
#normalized_counts.toptax=normalized_counts[sigtabxannot$value,]
#meannormcount=rowMeans(normalized_counts.toptax)
#sigtabxannot$mean.normcount=meannormcount

#get list of CRC samples by group from phyloseq metadata
#g1spls=rownames(sample_data(physeq)[sample_data(physeq)[[ref.var]]==comp.grp,])
#g0spls=rownames(sample_data(physeq)[sample_data(physeq)[[ref.var]]==ref.grp,])

#match CRC-samples-by-group with normalized counts table
#normcts.g1=normalized_counts.toptax[,g1spls]
#normcts.g1=as.data.frame(normcts.g1)
#rownames(normcts.g1)=rownames(normalized_counts.toptax)

#normcts.g0=normalized_counts.toptax[,g0spls]
#normcts.g0=as.data.frame(normcts.g0)
#rownames(normcts.g0)=rownames(normalized_counts.toptax)

#compute prevalence of each taxon by group (# of samples where taxon is >0 / total # of samples in the group)
#g1=adply(normcts.g1, 1, function(x) {(nnzero(x)/ncol(normcts.g1))*100})
#g1=g1$V1
#g0=adply(normcts.g0, 1, function(x) {(nnzero(x)/ncol(normcts.g0))*100})
#g0=g0$V1

#integrate prevalence info in grp
#sigtabxannot$grp1=g1
#sigtabxannot$grp0=g0

#reformat table to keep one column for prevalence group
#sigtabxannot.melt=pivot_longer(sigtabxannot,starts_with("grp"),names_to="Prevalence_grp",values_to="Prevalence_percent")

#Arrange df by prevalence group
#sigtabxannot.melt=sigtabxannot.melt %>% arrange(Prevalence_grp) 

#compute mean norm counts of each taxon by group 
#mean.grp1=as.data.frame(rowMeans(normcts.g1))
#colnames(mean.grp1)="Mean_normcount"
#mean.grp1$value=rownames(mean.grp1)
#mean.grp1$Prevalence_grp="grp1"

#mean.grp0=as.data.frame(rowMeans(normcts.g0))
#colnames(mean.grp0)="Mean_normcount"
#mean.grp0$value=rownames(mean.grp0)
#mean.grp0$Prevalence_grp="grp0"

#Make a column for mean normcount
#Meanct=rbind(mean.grp0,mean.grp1)

#if (all.equal(sigtabxannot.melt$value,Meanct$value)) {
#  if (all.equal(sigtabxannot.melt$Prevalence_grp,Meanct$Prevalence_grp)) {
#  sigtabxannot.melt=cbind(sigtabxannot.melt,Mean_normcount=Meanct$Mean_normcount)
#  }
#}  

#factor prevalence group to have grp1 first before grp0
#sigtabxannot.melt$Prevalence_grp=factor(sigtabxannot.melt$Prevalence_grp,levels=c("grp1","grp0"))

#if there is need to scale norm counts
#sigtabxannot.melt$mean.normcount=sigtabxannot.melt$mean.normcount/2

#plot DESeq 
#p1=ggplot(sigtabxannot.melt, aes(x=reorder(Genussp,log2FoldChange),y=log2FoldChange, color=Phylum)) +
#  geom_point(size=2) + coord_flip() + 
#      theme(axis.text.y=element_text(size=6.5)) +
#  geom_segment(aes(x=Genussp, xend=Genussp, y=0, yend=log2FoldChange),color="gray",size=0.25) + xlab("Lowest identifiable taxon") + ylab("log2FoldChange timepoint/baseline)") + ggtitle(label=NULL,subtitle="DESeq2 results at p-adj<0.05")
#p1 + theme(axis.text.y=element_text(size=10))

#plot percent prevalence
#p2.fixed=ggplot(sigtabxannot.melt, aes(fill=Prevalence_grp, x=reorder(Genussp,log2FoldChange),y=Prevalence_percent)) + geom_bar(position="stack",stat="identity") + scale_y_continuous() + 
#  theme(axis.text.y=element_blank(),axis.title.y=element_blank()) +
#  coord_flip() + ylab("Percent prevalence across samples per group") + ggtitle(label=NULL,subtitle="Per-group taxon stats")
#p2.fixed + theme(axis.text.y=element_text(size=10))

#plot mean normalized count
#sigtabxannot.melt$Group=sigtabxannot.melt$Prevalence_grp

#p3=ggplot(sigtabxannot.melt, aes(fill=Group, x=reorder(Genussp,log2FoldChange),y=Mean_normcount)) + geom_point(aes(x=reorder(Genussp,log2FoldChange),y=mean.normcount),size=0.65,color="gray1") + 
#  geom_line(aes(x=reorder(Genussp,log2FoldChange),y=mean.normcount), group=1,color="gray") +
#  geom_point(aes(x=reorder(Genussp,log2FoldChange),y=Mean_normcount,color=Group),size=1) + 
#  theme(axis.text.y=element_text(size=6.5),axis.title.y=element_blank()) +
#  coord_flip() + ylab("mean DESeq-normalized count") + ggtitle(label=NULL,subtitle="Per-group taxon stats")
#p3

#arrange plots
#gg1=ggplot_gtable(ggplot_build(p1)) 
#gg2=ggplot_gtable(ggplot_build(p2.fixed)) 
#gg3=ggplot_gtable(ggplot_build(p3)) 

#grid.arrange(gg1,gg2,gg3,ncol=2,widths=c(10/14,4/14))
#p3
#grid.arrange(gg1,gg2,gg3,ncol=2,nrow=2,widths=c(10/14,4/14))

#print("filter for taxa with mean norm count above median")
#p1=ggplot(sigtabxannot.melt[sigtabxannot.melt$mean.normcount>median(sigtabxannot.melt$mean.normcount),], aes(x=reorder(Genussp,log2FoldChange),y=log2FoldChange, color=Phylum)) +
#  geom_point(size=2) + coord_flip() + 
#      theme(axis.text.y=element_text(size=8)) +
#  geom_segment(aes(x=Genussp, xend=Genussp, y=0, yend=log2FoldChange),color="gray",size=0.25) + xlab("Lowest identifiable taxon") + ylab("log2FoldChange timepoint/baseline)") + ggtitle(label=NULL,subtitle="DESeq2 results at p-adj<0.05")
#p1 + theme(axis.text.y=element_text(size=10))

#p2.fixed=ggplot(sigtabxannot.melt[sigtabxannot.melt$mean.normcount>median(sigtabxannot.melt$mean.normcount),], aes(fill=Prevalence_grp, x=reorder(Genussp,log2FoldChange),y=Prevalence_percent)) + geom_bar(position="stack",stat="identity") + scale_y_continuous() + 
#  theme(axis.text.y=element_blank(),axis.title.y=element_blank()) +
#  coord_flip() + ylab("Percent prevalence across samples per group") + ggtitle(label=NULL,subtitle="Per-group taxon stats")
#p2.fixed + theme(axis.text.y=element_text(size=10))

#p3=ggplot(sigtabxannot.melt[sigtabxannot.melt$mean.normcount>median(sigtabxannot.melt$mean.normcount),], aes(fill=Group, x=reorder(Genussp,log2FoldChange),y=Mean_normcount)) + geom_point(aes(x=reorder(Genussp,log2FoldChange),y=mean.normcount),size=0.65,color="gray1") + 
#  geom_line(aes(x=reorder(Genussp,log2FoldChange),y=mean.normcount), group=1,color="gray") +
#  geom_point(aes(x=reorder(Genussp,log2FoldChange),y=Mean_normcount,color=Group),size=1) + 
#  theme(axis.text.y=element_text(size=8),axis.title.y=element_blank()) +
#  coord_flip() + ylab("mean DESeq-normalized count") + ggtitle(label=NULL,subtitle="Per-group taxon stats")
#p3

#arrange plots
#gg1=ggplot_gtable(ggplot_build(p1)) 
#gg2=ggplot_gtable(ggplot_build(p2.fixed)) 
#gg3=ggplot_gtable(ggplot_build(p3)) 

#grid.arrange(gg1,gg2,gg3,ncol=2,widths=c(10/14,4/14))
#p3
#grid.arrange(gg1,gg2,gg3,ncol=2,nrow=2,widths=c(10/14,4/14))

#sigtabxannot.melt.reord=sigtabxannot.melt[,c((7:ncol(sigtabxannot.melt)),(1:6))]

#datatable(sigtabxannot.melt.reord,caption="DESeq results with per-group taxon stats")