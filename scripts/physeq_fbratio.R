physeq_fbratio=function(physeq,)
physeq=readRDS("physeq_838_mtd2109xx.rds")
#physeq=subset_samples(physeq,Study=="LN-838"|Study=="LN-MIO")
#physeq=subset_samples(physeq,!is.na(Pathological_response_bin))
#physeq=subset_samples(physeq, Therapy_cycle_n_timepoint_n=="0")

physeq=tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

physeq.firm=subset_taxa(physeq,Phylum=="Firmicutes")
physeq.bac=subset_taxa(physeq,Phylum=="Bacteroidetes")

fb.ratio=as(otu_table(physeq.firm),"matrix")/as(otu_table(physeq.bac),"matrix")
fb.ratio=as.data.frame(t(fb.ratio))
colnames(fb.ratio)="fb.ratio"

mtd=as(sample_data(physeq.bac),"matrix")
mtd=as.data.frame(mtd)

if (all.equal(rownames(mtd),rownames(fb.ratio))) {
  
  mtd=data.frame(mtd,fb.ratio=fb.ratio)
  
}

mtd$scaled.fb.ratio=scale(mtd$fb.ratio,center=TRUE,scale=TRUE)

#ggplot(mtd,aes(Pathological_response_bin,fb.ratio,color=Pathological_response_bin)) + geom_jitter() + geom_boxplot(alpha=0.5,width=0.15,outlier.alpha=0) + stat_compare_means() + ggtitle("F/B ratio at baseline, by pathological response") + ylab("Firmicutes/Bacteroidetes ratio")


#print("filter samples with FB ratio SD > 2")

#mtd.filt=mtd[mtd$scaled.fb.ratio<=2,]

#ggplot(mtd.filt,aes(Pathological_response_bin,fb.ratio,color=Pathological_response_bin)) + geom_jitter() + geom_boxplot(alpha=0.5,width=0.15,outlier.alpha=0) + stat_compare_means() + ggtitle("F/B ratio at baseline, by pathological response") + ylab("Firmicutes/Bacteroidetes ratio")

#datatable(mtd[,c("Subject_number","Therapy_cycle_n_timepoint_n","fb.ratio","scaled.fb.ratio","Pathological_response_bin")],caption="List of samples and corresponding FB ratio")