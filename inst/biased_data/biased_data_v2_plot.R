theme_set(theme_bw())
theme_update(axis.text.x = element_text(size = 20),
             axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 20, margin=margin(20,0,0,0)),
             axis.title.y = element_text(size = 20, angle = 90, margin = margin(0,20,0,0)),
             panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                             colour = "gray"),
             panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                             colour = "gray"),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15),
             title = element_text(size = 30),
             strip.text = element_text(size = 20),
             strip.background = element_rect(fill="white"),
             legend.position = "bottom")

my_palette <- c("black", rgb(1,0.7,0.7), rgb(0.8,0.5,0.5), rgb(0.7,0.7, 1), rgb(0.5, 0.5, 0.8), rgb(0, 0.3, 0.6))
attach(biased_dataset)
load(file = paste0("inst/biased_data/biased_data_v2_jointmodel_n1_", n1, ".RData"))
load(file = paste0("inst/biased_data/biased_data_v2_twostage_n1_", n1, ".RData"))
df_samples <- data.frame(theta1 = mcmc_samples[,1], theta2 = mcmc_samples[,2], dist = 'Standard')
df_samples <- rbind(df_samples,
                    data.frame(theta1 = cut_samples[,1], theta2 = cut_samples[,2], dist = 'Cut'))

samples_laplace <- fast_rmvnorm(1e5, jointmodel_map, cov_laplace)
df_samples <- rbind(df_samples,
                    data.frame(theta1 = samples_laplace[,1], theta2 = samples_laplace[,2], dist = 'Standard Laplace'))

cut_laplace_samples <- fast_rmvnorm(1e5, twostage_map, laplace_cut_var)
df_samples <- rbind(df_samples,
                    data.frame(theta1 = cut_laplace_samples[,1], theta2 = cut_laplace_samples[,2], dist = 'Cut-Laplace'))

df_samples <- rbind(df_samples,
                    data.frame(theta1 = as.numeric(pbmi_samples[,1]), theta2 = as.numeric(pbmi_samples[,2]), dist = 'PBMI'))

tail(df_samples)

df_map <- data.frame(theta1 = c(jointmodel_map[1], twostage_map[1]),
                     theta2 = c(jointmodel_map[2], twostage_map[2]),
                     dist = c("Standard", "Cut"))


unique(factor(df_samples$dist))
# ## plot 95% credible region associated with mcmc_samples using stat_ellipse
g_std_vs_cut <- ggplot(df_samples, aes(x = theta1, y = theta2, color = factor(dist), linetype = factor(dist))) +
  stat_ellipse(level = 0.95, type = 'norm', linewidth = 2) +
  scale_linetype_manual(name = "", values = c(1,3,2, 1,3)) +
  scale_color_manual(name = "", values = c(my_palette[2], my_palette[3], my_palette[6], my_palette[4], my_palette[5])) +
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  xlim(-0.4, 1.3) + ylim(-0.4, 1.3) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 1) +
  geom_point(data = df_map, aes(x = theta1, y = theta2, color = factor(dist)),
             inherit.aes = FALSE,
             size = 4, shape = 4, stroke = 2)
  # ggtitle(paste0("n_1 = ", n1)) + #, subtitle = 'Standard Posterior in blue vs Cut in red vs PBMI in dark blue') +
  # theme(legend.position = "bottom")
g_std_vs_cut



ggsave(paste0('inst/biased_data/biased_data_v2_std_vs_twostage_n1_', n1, '.pdf'), g_std_vs_cut, width = 8, height = 5)


## plot density for theta1
g_theta1 <- ggplot(df_samples, aes(x = theta1, color = factor(dist), linetype = factor(dist))) +
  geom_density(linewidth = 1) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[6], my_palette[4], my_palette[5])) +
  scale_linetype_manual(name = "Distribution", values = c(1,3,1,3,1)) +
  xlab(expression(theta[1])) + ylab("Density") +
  xlim(-0.4,1.3)
g_theta1

## plot density for theta2
g_theta2 <- ggplot(df_samples, aes(x = theta2, color = factor(dist), linetype = factor(dist))) +
  geom_density(linewidth = 1) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[6], my_palette[4], my_palette[5])) +
  scale_linetype_manual(name = "Distribution", values = c(1,3,1,3,1)) +
  xlab(expression(theta[2])) + ylab("Density") +
  xlim(-.4, 1.3)
g_theta2


detach(biased_dataset)

