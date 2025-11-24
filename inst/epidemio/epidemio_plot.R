
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


load("inst/epidemio/epidemio_joint_mcmc_samples.RData")
df_samples <- data.frame(theta2 = joint_mcmc_samples[,14:15], dist = "Standard")
load('inst/epidemio/epidemio_cut_mcmc_samples.RData')
df_samples <- rbind(df_samples,
              data.frame(theta2.1 = theta2s[,1], theta2.2 = theta2s[,2], dist = "Cut"))
tail(df_samples)
load('inst/epidemio/epidemio_cutlaplace.RData')
twostage_mean <- colMeans(theta2s)
cut_laplace_samples <- fast_rmvnorm(1e6, twostage_mean, asympvar[14:15,14:15])
df_samples <- rbind(df_samples,
                    data.frame(theta2.1 = cut_laplace_samples[,1], theta2.2 = cut_laplace_samples[,2], dist = 'Cut-Laplace'))

load('inst/epidemio/epidemio_pbmi.RData')
df_samples <- rbind(df_samples,
                    data.frame(theta2.1 = pbmi_theta2[,1], theta2.2 = pbmi_theta2[,2], dist = 'PBMI'))

unique(factor(df_samples$dist))
jointmodel_mean <- colMeans(joint_mcmc_samples[,14:15])
df_map <- data.frame(theta2.1 = c(jointmodel_mean[1], twostage_mean[1]),
                     theta2.2 = c(jointmodel_mean[2], twostage_mean[2]),
                     dist = c("Standard", "Cut"))


g_std_vs_cut <- ggplot(df_samples, aes(x = theta2.1, y = theta2.2, color = factor(dist), linetype = factor(dist))) +
  stat_ellipse(level = 0.95, type = 'norm', linewidth = 2)
g_std_vs_cut <- g_std_vs_cut + scale_linetype_manual(name = "Distribution", values = c(1,3,1,1)) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[6], my_palette[4])) +
  xlab(expression(theta[2.2])) + ylab(expression(theta[2.1])) +
  # xlim(-0.4, 1.3) + ylim(-0.4, 1.3) +
  # geom_vline(xintercept = 0) + geom_hline(yintercept = 1) +
  geom_point(data = df_map, aes(x = theta2.1, y = theta2.2, color = factor(dist)),
             inherit.aes = FALSE,
             size = 4, shape = 4, stroke = 2)
g_std_vs_cut

ggsave('inst/epidemio/epidemio_std_vs_twostage.pdf', g_std_vs_cut, width = 8, height = 5)

df_samples$dist %>% unique
## plot histogram of first marginal
g_hist_theta2_1 <- ggplot(df_samples, aes(x = theta2.1, color = factor(dist), linetype = factor(dist))) +
  geom_density(aes(y=after_stat(density)), position="identity") +
  scale_linetype_manual(name = "Distribution", values = c(1, 2, 1, 1)) +
  scale_color_manual(name = "Distribution", values = c("Cut" = my_palette[2], "Cut-Laplace" = my_palette[3], "Standard" = my_palette[4], "PBMI" = my_palette[6])) +
  xlab(expression(theta[2.1])) +
  ylab("Density")
g_hist_theta2_1
ggsave('inst/epidemio/epidemio_hist_theta2_1.pdf', g_hist_theta2_1, width = 8, height = 5)

g_hist_theta2_2 <- ggplot(df_samples, aes(x = theta2.2, color = factor(dist), linetype = factor(dist))) +
  geom_density(aes(y=after_stat(density)), position="identity") +
  scale_linetype_manual(name = "Distribution", values = c(1, 2, 1, 1)) +
  scale_color_manual(name = "Distribution", values = c("Cut" = my_palette[2], "Cut-Laplace" = my_palette[3], "Standard" = my_palette[4], "PBMI" = my_palette[6])) +
  xlab(expression(theta[2.2])) +
  ylab("Density")
g_hist_theta2_2
ggsave('inst/epidemio/epidemio_hist_theta2_2.pdf', g_hist_theta2_2, width = 8, height = 5)

## now plot histograms of theta2.1 and theta2.2 on the same plot using facet_wrap
df_samples_long <- df_samples %>%
  pivot_longer(cols = c(theta2.1, theta2.2),
               names_to = "parameter",
               values_to = "value")
library(latex2exp)
df_samples_long$parameter <- factor(df_samples_long$parameter,
                                    levels = c("theta2.1", "theta2.2"),
                                     labels = c(TeX("$\\theta_{2,1}$"), TeX("$\\theta_{2,2}$")))

g_hist_theta2_both <- ggplot(df_samples_long, aes(x = value, color = factor(dist), linetype = factor(dist))) +
  geom_density(aes(y=after_stat(density)), position="identity") +
  scale_linetype_manual(name = "Distribution", values = c(1, 2, 1, 1)) +
  scale_color_manual(name = "Distribution", values = c("Cut" = my_palette[2], "Cut-Laplace" = my_palette[3], "Standard" = my_palette[4], "PBMI" = my_palette[6])) +
  xlab("Parameter value") +
  ylab("Density") +
  facet_wrap(~parameter, scales = "free", ncol = 1,
             labeller = label_parsed)

g_hist_theta2_both
ggsave('inst/epidemio/epidemio_hist_theta2_both.pdf', g_hist_theta2_both, width = 8, height = 5)

theta1hat = thetahat[1:13]
theta2hat = thetahat[14:15]

df_plot_poissonreg <- data.frame(x = Npop_normalized + theta2hat[1] + theta1hat * theta2hat[2],
     y = log(ncases))
ggplot(df_plot_poissonreg, aes(x = x, y = y)) +
  geom_point() +
  xlab(expression(log(mu[j]))) +
  ylab(expression(log(ncases[j]))) +
  ggtitle('Poisson regression fit')

gdata <- ggplot(df_plot_poissonreg, aes(x = exp(x), y = exp(y))) +
  geom_point() +
  # xlab add mu hat
  xlab(expression(hat(mu))) +
  ylab("number of cancer cases")
gdata <- gdata +
  theme(axis.text.x = element_text(size = 20, hjust = .75))
gdata
ggsave('inst/epidemio/epidemio_data.pdf', gdata, width = 8, height = 5)
