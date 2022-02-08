physeq_dist=function(physeq,method) {

ps.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=min(sample_sums(physeq)), replace=F)
#factor variables before declaring new variable for collating alphadiv vals

sd=min(sample_sums(physeq))

dist = phyloseq::distance(ps.rarefied, method=method)

}