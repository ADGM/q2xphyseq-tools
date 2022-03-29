t0dist_plot=function(dist.df,group) {

group=sym(group)
value=sym("value")

p=ggplot(dist.df,aes(!!group,dist.t0tn,color=!!group)) + geom_boxplot(outlier.alpha=0,width=0.2,alpha=0.60) + geom_jitter() + scale_color_viridis_d(option="H",begin=0.1,end=0.9)


}
