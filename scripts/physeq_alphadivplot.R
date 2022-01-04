physeq_alphadivplot=function(physeq,group) {

group=sym(group)
value=sym("value")

p=ggplot(alpha.vals,aes(x=!!group,y=!!value,color=!!group)) + facet_wrap(. ~ alphadiv_metric, scales="free_y") + geom_boxplot(width=0.3,alpha=0.5,outlier.alpha=0) +
geom_jitter(size=1,alpha=.8) + scale_x_discrete(guide = guide_axis(angle = 90)) + ggtitle(paste0("Alpha-diversity, by ",group)) + scale_color_viridis_d(option="H",begin=0.1,end=0.9) 

}
