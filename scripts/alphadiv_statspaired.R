alphadiv_statspaired=function(df,group) {

f=as.formula(paste0("value ~ ",group))

stats=df %>%
  group_by(alphadiv_metric) %>%
  rstatix::wilcox_test(f,paired=TRUE) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_xy_position(x=group)

}
