### Plotting functions 
library(viridis)
library(ggpubr)
library(rstatix)

# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
ROI_col = viridis(6)

# Violin plot functions
violin_plot_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, colour_palette) {
  # make violin plot
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none")
  
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p.adj")
  stat.test <- stat.test[stat.test$p.adj < 0.05, ]
  print(stat.test)
  
  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  
  #return object
  return(bxp)
}

# Boxplot functions
boxplot_2grp_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, grp_variable,colour_palette) {
  
  # format data
  data_table$g <- data_table[, grp_variable]
  data_table$g <- as.factor(data_table$g)
  data_table$group <- factor(rep(c("grp1", "grp2"), nrow(data_table)/2))
  
  # order factor levels
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  data_table$g <- data_table[, grp_variable]
  
  stat.test <- data_table %>%
    group_by(x) %>%
    t_test(y ~ g) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  stat.test <- stat.test[order(unlist(stat.test[,"statistic"])),]
  fact_lvls <- stat.test$x[order(unlist(stat.test[,"statistic"]))]
  print(fact_lvls)
  data_table[, x_variable] <- factor(data_table[, x_variable] , levels = fact_lvls)
  
  # plot
  bxp <- ggboxplot(
    data_table, x = x_variable, y = y_variable,
    color = grp_variable, palette = colour_palette
  ) + ylab(y_lab) + xlab(x_lab)
  
  # add stat tests
  bxp + geom_pwc(
    aes(group = g), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )
}
