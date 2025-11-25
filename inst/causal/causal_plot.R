# Causal inference histogram
library(ggplot2)
library(dplyr)
library(MatchIt)
data("lalonde")

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


lalonde %>% summary()
lalonde %>% group_by(treat) %>% summarize(median_re75 = median(re75),
                                     median_re78 = median(re78),
                                     mean_re75 = mean(re75),
                                     mean_re78 = mean(re78),
                                     n075 = sum(re75 == 0),
                                     n = n())

df_histogram <- rbind(cbind(lalonde$re75, lalonde$treat,
                            rep(c('pre-intervention (1975)'), 614)),
                      cbind(lalonde$re78, lalonde$treat,
                            rep(c('post-intervention (1978)'), 614))) %>% as.data.frame()
colnames(df_histogram) <- c('income', 'treatment', 'year')
df_histogram$year <- factor(df_histogram$year,
                            levels = c('pre-intervention (1975)',
                                       'post-intervention (1978)'))
df_histogram$income <- df_histogram$income %>% as.character()  %>% as.numeric()

df_histogram %>% group_by(treatment) %>% summarize(n = n())


## Histogram of income pre-intervention, untreated on top and treated underneath
## + vertical segments for median income
ghist_upsidedown_pre <- ggplot() + geom_histogram(data=df_histogram %>% filter(year == "pre-intervention (1975)", treatment == 0), aes(x = income, y = ..density.., fill = treatment), alpha=0.4, bins = 50) +
  geom_histogram(data=df_histogram %>% filter(year == "pre-intervention (1975)", treatment == 1), aes(x = income, y = -..density.., fill = treatment), alpha=0.4, bins = 50) +
  ylim(-1.3e-3, 1.3e-3) +
  scale_fill_manual(name = "", values = c( '0' = '#56B4E9',
                                '1' = '#CC79A7'), labels = c('Control', 'Treated'))

m0 <- df_histogram %>% filter(year == "pre-intervention (1975)", treatment == 0) %>% summarize(m = median(income)) %>% pull(m)
m1 <- df_histogram %>% filter(year == "pre-intervention (1975)", treatment == 1) %>% summarize(m = median(income)) %>% pull(m)

ghist_upsidedown_pre <- ghist_upsidedown_pre + annotate("segment", x = m0, xend = m0, y = 0, yend = Inf, colour = '#56B4E9') +
  annotate("segment", x = m1, xend = m1, y = 0, yend = -Inf, colour = '#CC79A7')
ghist_upsidedown_pre <- ghist_upsidedown_pre + xlab('income (in US dollars)')
ghist_upsidedown_pre
ggsave(filename = "inst/causal/causal_inference_income_histogram_pre.pdf", ghist_upsidedown_pre, width = 8, height = 5)


ghist_upsidedown_post <- ggplot() + geom_histogram(data=df_histogram %>% filter(year == "post-intervention (1978)", treatment == 0), aes(x = income, y = ..density.., fill = treatment), alpha=0.4, bins = 50) +
  geom_histogram(data=df_histogram %>% filter(year == "post-intervention (1978)", treatment == 1), aes(x = income, y = -..density.., fill = treatment), alpha=0.4, bins = 50) +
  ylim(-2.5e-4, 2.5e-4) +
  scale_fill_manual(name = "", values = c( '0' = '#56B4E9',
                                '1' = '#CC79A7'), labels = c('Control', 'Treated'))
m0_post <- df_histogram %>% filter(year == "post-intervention (1978)", treatment == 0) %>% summarize(m = median(income)) %>% pull(m)
m1_post <- df_histogram %>% filter(year == "post-intervention (1978)", treatment == 1) %>% summarize(m = median(income)) %>% pull(m)
ghist_upsidedown_post <- ghist_upsidedown_post + annotate("segment", x = m0_post, xend = m0_post, y = 0, yend = Inf, colour = '#56B4E9') +
  annotate("segment", x = m1_post, xend = m1_post, y = 0, yend = -Inf, colour = '#CC79A7')
ghist_upsidedown_post <- ghist_upsidedown_post + xlab('income (in US dollars)')
ghist_upsidedown_post
ggsave(filename = "inst/causal/causal_inference_income_histogram_post.pdf", ghist_upsidedown_post, width = 8, height = 5)


load('inst/causal/causal_twostage_results.RData')

## Standardize the data
standardize <- function(x) {
  (x-mean(x))/sd(x)
}
lalonde$race <- factor(lalonde$race, levels = c('white', 'hispan', 'black'))
lalonde$age <- standardize(lalonde$age)
lalonde$educ <- standardize(lalonde$educ)
lalonde$re74 <- standardize(lalonde$re74)
lalonde$re75 <- standardize(lalonde$re75)
lalonde$age_sq <- standardize((lalonde$age)^2)

stage_1_model_mat <- model.matrix(treat ~ age + age_sq + educ + race + married + nodegree + re74 + re75,
                                  data = lalonde)


propscores <- stage_1_model_mat %*% theta1hat

## histogram of propscores for treated and control groups
df_propscores <- data.frame(propscore = as.vector(propscores),
                            treat = lalonde$treat)
g_propscorebytreat <- ggplot(df_propscores, aes(x = propscore, fill = factor(treat))) +
  geom_histogram(bins = 100, position = 'identity', alpha = 0.5) +
  labs(x = 'propensity score',
       y = 'count',
       fill = 'treatment') +
  scale_fill_manual(name = "", labels = c('Control', 'Treated'),
                    values = c('#56B4E9', '#CC79A7'))
g_propscorebytreat
ggsave(filename = "inst/causal/causal_inference_propscores.pdf", g_propscorebytreat, width = 8, height = 5)


## Causal inference results plot comparing the treatment effect under cut Bayesian and PB
# causal_df <- readRDS(df, file ='inst/causal/causal_twostage_results.rds')
# head(causal_df)
## rename methods
df_twostage$method <- dplyr::recode(df_twostage$method,
                                  'PBMI' = 'PBMI',
                                  'PB refresh' = 'PBMI-refresh',
                                  'cut Bayesian' = 'Cut')
#
my_palette <- c("black", rgb(1,0.7,0.7), rgb(0.8,0.5,0.5), rgb(0.7,0.7, 1), rgb(0.5, 0.5, 0.8), rgb(0, 0.3, 0.6), rgb(0, 0.6, 0.3))
#
#
# load('inst/causal/causal_twostep_laplace_results.RData')
# laplace_cut_samples <- fast_rmvnorm(1e5, laplace_cut_results$thetahat, laplace_cut_results$asympvar)
# laplace_cut_samples <- laplace_cut_samples[,11:16]
# causal_df_laplace <- data.frame(laplace_cut_samples, method = 'Cut-Laplace')
# colnames(causal_df_laplace)[1:6] <- colnames(causal_df)[1:6]
# causal_df <- rbind(causal_df, causal_df_laplace)

df_twostage %>% group_by(method) %>% summarise(mean_treat = mean(V1), sd_treat = sd(V1))

df_twostage$method %>% unique()

causal_plot <- ggplot(df_twostage, aes(x = V1, color = method, linetype = method)) + geom_density(size=1) +
  xlab('treatment') +
  scale_color_manual(name = "Distribution", values = c(my_palette[2],  my_palette[6], my_palette[7]))  +
  scale_linetype_manual(name = "Distribution", values = c(1, 4, 4)) +
  scale_x_continuous(breaks = seq(-.3, 0.9, 0.3))
causal_plot
#
# ggsave('inst/causal/causal_result.pdf', causal_plot, width = 8, height = 5)
ggsave('inst/causal/causal_treatcoef.pdf', causal_plot, width = 8, height = 5)


