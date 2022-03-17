t0dist_plot=function(dist.df,group) {

group=sym(group)
value=sym("value")

p=ggplot(dist.df,aes(Cycle_agg2,dist.t0tn,color=Cycle_agg2)) + geom_boxplot(outlier.alpha=0,width=0.2,alpha=0.60) + geom_jitter() + scale_color_viridis_d(option="H",begin=0.1,end=0.9)


}
