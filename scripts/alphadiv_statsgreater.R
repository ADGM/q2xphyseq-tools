alphadiv_stats.greater=function(df,group) {

f=as.formula(paste0("value ~ ",group))

stats=df %>%
  group_by(alphadiv_metric) %>%
  rstatix::wilcox_test(f,alternative="greater") %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_xy_position(x=group)

}
