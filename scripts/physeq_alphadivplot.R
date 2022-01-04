physeq_alphadivplot=function(physeq,group) {

stats=alpha.vals %>%
  group_by(alphadiv_metric) %>%
  rstatix::wilcox_test(value ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

p=ggplot(subset(alpha.vals, !is.na(`group`)),aes(x=`group`,y=value,color=`group`)) + facet_wrap(. ~ alphadiv_metric, scales="free_y") + 
  #geom_violin(alpha=0.65,width=.85) + 
  geom_boxplot(width=0.15) + geom_jitter(size=0.8) + scale_x_discrete(guide = guide_axis(angle = 90)) + ggtitle("Alpha-diversity by `group`") + scale_color_viridis_d(direction = -1,option="C",begin=0.5,end=0.9)

print(datatable(stats,caption="Wilcoxon test by `group`"))

}
