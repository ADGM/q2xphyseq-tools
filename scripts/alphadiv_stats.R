alphadiv_stats=function(df,plot,group) {

f=as.formula(paste0("value ~ ",group))

stats=df %>%
  group_by(alphadiv_metric) %>%
  rstatix::wilcox_test(f) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_xy_position(x=group)

datatable(stats,caption="Wilcoxon test by ",group)
}
