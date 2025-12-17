# source("~/Desktop/P9/P9_TimeVaryingCointegration/Radiative_forcings.r")

###########
# tabel
ic_table <- expand.grid(
  m = 1:nrow(aic),
  r = 1:ncol(aic)
) %>%
  mutate(
    AIC  = as.vector(aic),
    BIC  = as.vector(bic),
    HQ   = as.vector(hann)
  )
ic_table

# print(
#   xtable(ic_table, digits = 4),
#   include.rownames = FALSE,
#   type = "latex"
# )

###########
# plot beta_t
plot_beta <- function(k = ncol(y), m = m, var_names = NULL) {
  # beta_block: vælg det relevante sæt af kolonner
  beta_block <- betat[, (k * (m - 1) + 1):(k * m)]
  
  # hvis var_names ikke er givet, lav default navne
  if (is.null(var_names)) {
    var_names <- paste0("Var", 1:k)
  }
  colnames(beta_block) <- var_names
  
  # years for beta
  ind <- seq(from = p + 2, by = 1, length.out = n - p - 1)
  years_beta <- merged_data$year[ind]
  
  # build tidy data.frame for ggplot
  df_beta <- data.frame(year = years_beta, beta_block)
  df_long <- tidyr::pivot_longer(
    df_beta,
    cols = all_of(var_names),
    names_to = "variable",
    values_to = "beta"
  )
  
  # enforce ordering in the facets
  df_long$variable <- factor(df_long$variable, levels = var_names)
  
  # plot
  p_beta <- ggplot(df_long, aes(x = year, y = beta)) +
    geom_line() +
    geom_point(size = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    labs(
      title = expression("Estimated time-varying cointegration coefficients " ~ beta[t]),
      x = "Year",
      y = expression(beta[t])
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(angle = 0, vjust = 0.5) 
    )
  
  return(p_beta)
}
m
r
plot_beta(k=k, m=m, var_names = c("GMST", "ERF"))
plot_beta(k=k, m=5, var_names = c("GMST", "ERF"))
