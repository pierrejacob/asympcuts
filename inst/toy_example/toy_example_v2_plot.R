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

## Load Joint Model results
load(paste0("inst/toy_example/toy_example_v2_jointmodel_scenario", SCENARIO, ".RData"))
df_samples <- data.frame(theta1 = mcmc_samples[,1], theta2 = mcmc_samples[,2], dist = "Standard", laplace = FALSE)
samples_laplace <- mvtnorm::rmvnorm(1e6, jointmodel_map, cov_laplace)
df_samples <- rbind(df_samples,
                    data.frame(theta1 = samples_laplace[,1], theta2 = samples_laplace[,2], dist = "Standard Laplace", laplace = TRUE))
## Load Two Stage results
load(paste0("inst/toy_example/toy_example_v2_twostage_scenario", SCENARIO, ".RData"))
df_samples <- rbind(df_samples,
                    data.frame(theta1 = cut_samples[,1], theta2 = cut_samples[,2], dist = "Cut", laplace = FALSE))
## Comparison Cut versus Cut-Laplace versus PBMI
cut_laplace_samples <- mvtnorm::rmvnorm(1e6, twostage_map, laplace_cut_var)
df_samples <- rbind(df_samples, data.frame(theta1 = cut_laplace_samples[,1], theta2 = cut_laplace_samples[,2],
                                           dist = "Cut-Laplace", laplace = TRUE))

## create data frame with jointmodel_map and twostage_map
df_map <- data.frame(theta1 = c(jointmodel_map[1], twostage_map[1]),
                     theta2 = c(jointmodel_map[2], twostage_map[2]),
                     dist = c("Standard", "Cut"))

#### Plot Standard versus Modular Bayesian approaches, and their Laplace approximation ####

## plot 95% credible region associated with mcmc_samples using stat_ellipse
g_std_vs_cut <- ggplot(df_samples, aes(x = theta1, y = theta2, color = factor(dist), linetype = factor(dist))) +
  stat_ellipse(level = 0.95, type = 'norm', linewidth = 2) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[4], my_palette[5])) +
  scale_linetype_manual(name = "Distribution", values = c(1,3,1,3)) +
  geom_point(data = df_map, aes(x = theta1, y = theta2, color = factor(dist)),
             inherit.aes = FALSE,
             size = 4, shape = 4, stroke = 2) +
  xlim(-0.16,0.41) + ylim(1.27, 1.61) +
  xlab(expression(theta[1])) + ylab(expression(theta[2]))   # ggtitle(paste0("Scenario ", SCENARIO), subtitle = 'Standard Posterior (and Laplace approx.) in blue vs Cut (and Laplace approx.) in red')
g_std_vs_cut
ggsave(paste0('inst/toy_example/toy_example_v2_std_vs_cut_scenario', SCENARIO, '.pdf'), g_std_vs_cut, width = 8, height = 5)

## plot density for theta1
g_theta1 <- ggplot(df_samples, aes(x = theta1, color = factor(dist), linetype = factor(dist))) +
  geom_density(linewidth = 1) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[4], my_palette[5])) +
  scale_linetype_manual(name = "Distribution", values = c(1,3,1,3)) +
  xlab(expression(theta[1])) + ylab("Density") +
  xlim(-0.16,0.41)
g_theta1

## plot density for theta2
g_theta2 <- ggplot(df_samples, aes(x = theta2, color = factor(dist), linetype = factor(dist))) +
  geom_density(linewidth = 1) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[4], my_palette[5])) +
  scale_linetype_manual(name = "Distribution", values = c(1,3,1,3)) +
  xlab(expression(theta[2])) + ylab("Density") +
  xlim(1.27, 1.61)
g_theta2


#### Plot Cut Posterior versus PBMI approaches ####

df_samples <- data.frame(theta1 = cut_samples[,1], theta2 = cut_samples[,2], dist = "Cut")
df_samples <- rbind(df_samples, data.frame(theta1 = cut_laplace_samples[,1], theta2 = cut_laplace_samples[,2],
                                           dist = "Cut-Laplace"))
df_samples <- rbind(df_samples, data.frame(theta1 = pbmi_samples[,1], theta2 = pbmi_samples[,2],
                                           dist = "PBMI"))


g_cut_vs_pbmi <- ggplot(df_samples, aes(x = theta1, y = theta2, color = factor(dist), linetype = factor(dist))) +
  # geom_point(alpha = 0.3) +
  stat_ellipse(level = 0.95, type = 'norm', linewidth = 2) +
  scale_color_manual(name = "Distribution", values = c(my_palette[2], my_palette[3], my_palette[6])) +
  scale_linetype_manual(name = "Distribution", values = c(1,3,1)) +
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  xlim(-0.2,0.14) + ylim(1.4, 1.7)
  # ggtitle(paste0("Scenario ", SCENARIO), subtitle = 'Cut (and Laplace approx.) in red, PBMI in dark blue')

g_cut_vs_pbmi
ggsave(paste0('inst/toy_example/toy_example_v2_cut_vs_pbmi_scenario', SCENARIO, '.pdf'), g_cut_vs_pbmi, width = 8, height = 5)

### all in one plot
df_samples <- data.frame(theta1 = mcmc_samples[,1], theta2 = mcmc_samples[,2], dist = "Standard", laplace = FALSE)
df_samples <- rbind(df_samples,
                    data.frame(theta1 = samples_laplace[,1], theta2 = samples_laplace[,2], dist = "Standard Laplace", laplace = TRUE))
df_samples <- rbind(df_samples,
                    data.frame(theta1 = cut_samples[,1], theta2 = cut_samples[,2], dist = "Cut", laplace = FALSE))
df_samples <- rbind(df_samples, data.frame(theta1 = cut_laplace_samples[,1], theta2 = cut_laplace_samples[,2],
                                           dist = "Cut-Laplace", laplace = TRUE))

df_samples <- rbind(df_samples, data.frame(theta1 = pbmi_samples[,1], theta2 = pbmi_samples[,2],
                                           dist = "PBMI", laplace = FALSE))
df_samples$dist <- factor(df_samples$dist,
                                  levels = c("Cut", "Cut-Laplace", "Standard", "Standard Laplace", "PBMI"))
g_std_vs_cut <- ggplot(df_samples, aes(x = theta1, y = theta2, color = factor(dist), linetype = factor(dist))) +
  stat_ellipse(level = 0.95, type = 'norm', linewidth = 2) +
  scale_color_manual(name = "", values = c("Cut"=my_palette[2], "Cut-Laplace"=my_palette[3], "Standard"=my_palette[4], "Standard Laplace"=my_palette[5], "PBMI" = my_palette[6])) +
  scale_linetype_manual(name = "", values = c(1,3,1,3, 1)) +
  geom_point(data = df_map, aes(x = theta1, y = theta2, color = factor(dist)),
             inherit.aes = FALSE,
             size = 4, shape = 4, stroke = 2) +
  xlim(-0.16,0.41) + ylim(1.27, 1.61) +
  xlab(expression(theta[1])) + ylab(expression(theta[2]))   # ggtitle(paste0("Scenario ", SCENARIO), subtitle = 'Standard Posterior (and Laplace approx.) in blue vs Cut (and Laplace approx.) in red')
g_std_vs_cut
ggsave(paste0('inst/toy_example/toy_example_v2_std_vs_cut_scenario', SCENARIO, '.pdf'), g_std_vs_cut, width = 8, height = 5)


load(file = paste0("inst/toy_example/toy_example_v2_twostage_coverage_scenario", SCENARIO, ".RData"))

## Compute coverage

## For theta1:
cutlaplace_results <- c(credible_regions_laplace %>% filter(parameter == "theta1") %>%
  mutate(coverage = ifelse((lower <= theta1star) & (upper >= theta1star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage)) %>% pull())
## For theta2:
cutlaplace_results <- c(cutlaplace_results, credible_regions_laplace %>% filter(parameter == "theta2") %>%
  mutate(coverage = ifelse((lower <= theta2star) & (upper >= theta2star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage)) %>% pull())
## For theta1 + theta2:
cutlaplace_results <- c(cutlaplace_results, credible_regions_laplace %>% filter(parameter == "theta1p2") %>%
  mutate(coverage = ifelse((lower <= (theta1star + theta2star)) & (upper >= (theta1star + theta2star)), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage)) %>% pull())

cutlaplace_results

## For theta1:
pbmi_results <- c(credible_regions_pbmi %>% filter(parameter == "theta1") %>%
  mutate(coverage = ifelse((lower <= theta1star) & (upper >= theta1star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage)) %>% pull())
## For theta2:
pbmi_results <- c(pbmi_results, credible_regions_pbmi %>% filter(parameter == "theta2") %>%
  mutate(coverage = ifelse((lower <= theta2star) & (upper >= theta2star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage)) %>% pull())
## For theta1 + theta2:
pbmi_results <- c(pbmi_results, credible_regions_pbmi %>% filter(parameter == "theta1p2") %>%
  mutate(coverage = ifelse((lower <= (theta1star + theta2star)) & (upper >= (theta1star + theta2star)), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage)) %>% pull())
pbmi_results

## Create data frame to export as a table using knitr::kable
df_coverage <- data.frame(Parameter = c("$\\theta_1$", "$\\theta_2$", "$\\theta_1 + \\theta_2$"),
                          Cut_Laplace = cutlaplace_results,
                          PBMI = pbmi_results)
df_coverage <- df_coverage %>% setNames(c("Parameter","Cut-Laplace", "PBMI"))
coverage_table <- knitr::kable(df_coverage, digits = 2, row.names = NA, format = 'latex', escape = FALSE, booktabs = TRUE)
coverage_table

coverage_table %>%  cat(., file = paste0("inst/toy_example/toy_example_v2_coverage_scenario", SCENARIO, ".tex"))


