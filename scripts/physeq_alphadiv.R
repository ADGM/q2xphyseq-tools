physeq_alphadiv=function(physeq,group) {

ps.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=min(sample_sums(physeq)), replace=F)
#factor variables before declaring new variable for collating alphadiv vals

sd=min(sample_sums(physeq)
print_sd=paste0("Rarefied at "sd," depth"))

print(print_sd)

alpha.vals=estimate_richness(ps.rarefied, measures=c("Observed","Chao1","Shannon","Fisher","Simpson","InvSimpson"))
rownames(alpha.vals)=gsub("X","",rownames(alpha.vals))
alpha.vals=subset(alpha.vals, select = -c(se.chao1))


if (all.equal(rownames(alpha.vals),rownames(sample_data(ps.rarefied)))) {
  alpha.vals=cbind(sample_data(ps.rarefied),alpha.vals)
  alpha.vals=pivot_longer(alpha.vals,cols=((ncol(sample_data(ps.rarefied))+1):ncol(alpha.vals)),names_to="alphadiv_metric",values_to="value")
}

}

