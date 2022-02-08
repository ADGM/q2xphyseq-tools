alphavals_stats=function(alpha.vals,group) {

f=as.formula(paste0("value ~ ",group))

stats=alpha.vals %>%
  group_by(alphadiv_metric) %>%
  rstatix::wilcox_test(f) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_xy_position(x=group)

  # Add 10% spaces between the p-value labels and the plot border
#p + stat_pvalue_manual(stats,size = 2.5,label = paste0("        p = ","{p.adj}"),remove.bracket = TRUE,y.position=0.5)

#datatable(stats,caption="Wilcoxon test by ",group)
}
