# create figure 5: multi-panel plot of results from covariate models 
# requires running top ranked models to obtain results to plot

#### packages ####
library(tidyverse)
library(cowplot)
library(MCMCvis)
library(tidybayes)
library(bayesplot)
library(ggdist)
library(HDInterval)
library(bayestestR)
library(scales)
library(easystats)

#### load in model results - juvenile survival ####
# load in from where they are saved
load("age-model-v8.RData") # time dependent model
load("phijuv_apc1_ssst_phi3c.RData") # covariate model

# beta[2] - autumn storms
juv_beta2 <- mcmc.phijuv.apc1.ssst.phi3c$samples %>% 
  gather_draws(`beta[2]`) # spread_draws from tidybayes - get mcmc iterations into an object

# beta[3] - sst
juv_beta3 <- mcmc.phijuv.apc1.ssst.phi3c$samples %>% 
  gather_draws(`beta[3]`) # spread_draws from tidybayes - get mcmc iterations into an object

# create object with all mcmc iterations of autumn storms and ssst
juv_betas <- mcmc.phijuv.apc1.ssst.phi3c$samples %>%
  gather_draws(`beta[2]`, `beta[3]`)

juv_betas$.variable <- recode(juv_betas$.variable, `beta[2]` = "Autumn storms", `beta[3]` = "Summer SST ºC")

# create density plots showing posterior distribution and proportion of posterior 
# with the same sign as the mean

# autumn storms
juv_beta2$.variable <- recode(juv_beta2$.variable, `beta[2]` = "Autumn storms")

(dens_apc1 <- juv_beta2 %>%  
    ggplot(aes(x = .value, y = forcats::fct_rev(.variable), fill = after_stat(x < 0))) +
    stat_halfeye(.width = c(0.50, 0.95), point_interval = mean_hdi, alpha = 0.5, size = 0.5) +
    xlab("Estimate") + ylab("") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "gray50") +
    scale_y_discrete(expand = expansion(mult = c(0.4, 1.4))) +
    scale_fill_manual(values = c("gray80", "#4B0055"), guide = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6, colour = "black"),
          axis.text.y = element_text(size = 7, colour = "black", angle = 90, hjust = -0.15),
          axis.title = element_text(size = 8), 
          axis.ticks.y = element_blank(), 
          axis.ticks.x = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2)))

# summer SST
juv_beta3$.variable <- recode(juv_beta3$.variable, `beta[3]` = "Summer SST ºC")

(dens_ssst <- juv_beta3 %>%  
    ggplot(aes(x = .value, y = forcats::fct_rev(.variable), fill = after_stat(x < 0))) +
    stat_halfeye(.width = c(0.50, 0.95), point_interval = mean_hdi, alpha = 0.5, size = 0.5) +
    xlab("Estimate") + ylab("") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "gray50") +
    scale_y_discrete(expand = expansion(mult = c(0.4, 1.4))) +
    scale_fill_manual(values = c("gray80", "#4B0055"), guide = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6, colour = "black"),
          axis.text.y = element_text(size = 7, colour = "black", angle = 90, hjust = -0.15),
          axis.title = element_text(size = 8), 
          axis.ticks.y = element_blank(), 
          axis.ticks.x = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2)))

### plot logistic regression

# time-depdendent estimates of juvenile survival 
guille <- data.frame(MCMCsummary(mcmc.guille.agem8$samples, params = "phi.1"))
guille <- rename(guille, q2.5 = X2.5., q50 = X50., q97.5 = X97.5.)                    

guille$year <- seq(1995, 2020, 1)

cov <- read.csv("standardised_covariates_for_models_1995-2020.csv")
cov2 <- read.csv("covariates_for_models_1995-2020.csv")

# add covariates 
guille$apc1 <- cov$son_pc1   # autumn storms
guille$ssst <- cov$ssst      # standardised summer sst
guille$ssst_n <- cov2$ssst   # summer sst on normal scale

# plot autumn storms effect 
beta1 <- c(mcmc.phijuv.apc1.ssst.phi3c$samples$chain1[,'beta[1]'], # intercept, chain 1
           mcmc.phijuv.apc1.ssst.phi3c$samples$chain2[,'beta[1]'], 
           mcmc.phijuv.apc1.ssst.phi3c$samples$chain3[,'beta[1]']) # intercept, chain 3
beta2 <- c(mcmc.phijuv.apc1.ssst.phi3c$samples$chain1[,'beta[2]'], # slope, chain 1
           mcmc.phijuv.apc1.ssst.phi3c$samples$chain2[,'beta[2]'], 
           mcmc.phijuv.apc1.ssst.phi3c$samples$chain3[,'beta[2]']) # slope, chain 2
beta3 <- c(mcmc.phijuv.apc1.ssst.phi3c$samples$chain1[,'beta[3]'], # slope, chain 1
           mcmc.phijuv.apc1.ssst.phi3c$samples$chain2[,'beta[3]'], 
           mcmc.phijuv.apc1.ssst.phi3c$samples$chain3[,'beta[3]']) # slope, chain 2

# predict survival for each iteration
predicted_survival <- matrix(NA, 
                             nrow = length(beta1), 
                             ncol = length(guille$apc1))
for (i in 1:length(beta1)){
  for (j in 1:length(guille$apc1)){
    predicted_survival[i,j] <- plogis(beta1[i] + 
                                        beta2[i] * guille$apc1[j])
  }
}

# calculate posterior mean and CI
mean_survival <- apply(predicted_survival, 2, mean)
lci <- apply(predicted_survival, 2, quantile, prob = 2.5/100)
uci <- apply(predicted_survival, 2, quantile, prob = 97.5/100)
ord <- order(guille$apc1)
df <- data.frame(apc1 = guille$apc1[ord],
                 survival = mean_survival[ord],
                 lci = lci[ord],
                 uci = uci[ord])

# ggplot
(apc1_plot <- df %>%
    ggplot() + 
    aes(x = apc1, y = survival) + 
    geom_line(col = "#4B0055") + 
    geom_ribbon(aes(ymin = lci, ymax = uci), 
                fill = "#4B0055", 
                alpha = 0.5) + 
    geom_point(data = guille, aes(x = apc1, y = mean), col = "gray40", size = 0.5) +
    ylim(0,1) + 
    labs(x = "Autumn storms", y = "Juvenile survival") +
    theme_bw() + 
    annotate(geom = "text", x = -3, y = 0.05, label = "less\nstormy", size = 2) +
    annotate(geom = "text", x = 1.8, y = 0.05, label = "more\nstormy", size = 2) +
    theme(axis.text.x = element_text(size = 6, colour = "black"), 
          axis.text.y = element_text(size = 6, colour = "black", angle = 90), 
          axis.title = element_text(size = 8), 
          axis.ticks = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2))) 

### sst impact
# predict survival for each iteration
ssst_survival <- matrix(NA, 
                        nrow = length(beta1), 
                        ncol = length(guille$ssst))
for (i in 1:length(beta1)){
  for (j in 1:length(guille$ssst)){
    ssst_survival[i,j] <- plogis(beta1[i] + 
                                   beta3[i] * guille$ssst[j])
  }
}

# calculate posterior mean and CI
mean_ssst_survival <- apply(ssst_survival, 2, mean)
lci_s <- apply(ssst_survival, 2, quantile, prob = 2.5/100)
uci_s <- apply(ssst_survival, 2, quantile, prob = 97.5/100)
ord <- order(guille$ssst)
df2 <- data.frame(ssst = guille$ssst[ord],
                  survival = mean_ssst_survival[ord],
                  lci = lci_s[ord],
                  uci = uci_s[ord], 
                  ssst_n = guille$ssst_n[ord])

# plot
(ssst_plot <- df2 %>%
    ggplot() + 
    aes(x = ssst_n, y = survival) + 
    geom_line(col = "#4B0055") + 
    geom_ribbon(aes(ymin = lci, ymax = uci), 
                fill = "#4B0055", 
                alpha = 0.5) + 
    geom_point(data = guille, aes(x = ssst_n, y = mean), col = "gray40", size = 0.5) +
    ylim(0,1) + 
    labs(x = "Summer SST ºC", y = "Juvenile survival") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 6, colour = "black"), 
          axis.text.y = element_text(size = 6, colour = "black", angle = 90), 
          axis.title = element_text(size = 8), 
          axis.ticks = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2))) 

#### recruitment age 5 - lag winter sst ####
load("alpha5_lwsst.RData") # load in from wherever saved

# beta[2] - lag winter sst
psi5_beta2 <- mcmc.alpha5.lwsst$samples %>% 
  gather_draws(`beta[2]`) 

# time-dependent estimates of recruitment age5
alpha5 <- data.frame(MCMCsummary(mcmc.guille.agem8$samples, params = "alpha.5"))
alpha5 <- rename(alpha5, q2.5 = X2.5., q50 = X50., q97.5 = X97.5.)                    

alpha5$year <- seq(1995, 2020, 1)

# no juveniles in first four years of model
cov$lag_wsst[1] <- 0
cov$lag_wsst[2] <- 0
cov$lag_wsst[3] <- 0
cov$lag_wsst[4] <- 0

cov2$lag_wsst[1] <- 0
cov2$lag_wsst[2] <- 0
cov2$lag_wsst[3] <- 0
cov2$lag_wsst[4] <- 0

# add covariates to recruitment age 5 dataframe
alpha5$lag_wsst <- cov$lag_wsst
alpha5$lag_wsst_n <- cov2$lag_wsst

# plot winter storms effect 
beta1 <- c(mcmc.alpha5.lwsst$samples$chain1[,'beta[1]'], # intercept, chain 1
           mcmc.alpha5.lwsst$samples$chain2[,'beta[1]'], 
           mcmc.alpha5.lwsst$samples$chain3[,'beta[1]']) # intercept, chain 3
beta2 <- c(mcmc.alpha5.lwsst$samples$chain1[,'beta[2]'], # slope, chain 1
           mcmc.alpha5.lwsst$samples$chain2[,'beta[2]'], 
           mcmc.alpha5.lwsst$samples$chain3[,'beta[2]']) # slope, chain 2

# predict survival for each iteration
predicted_recruitment <- matrix(NA, 
                                nrow = length(beta1), 
                                ncol = length(alpha5$lag_wsst))
for (i in 1:length(beta1)){
  for (j in 1:length(alpha5$lag_wsst)){
    predicted_recruitment[i,j] <- plogis(beta1[i] + 
                                           beta2[i] * alpha5$lag_wsst[j])
  }
}

# calculate posterior mean and CI
mean_recruitment <- apply(predicted_recruitment, 2, mean)
lci <- apply(predicted_recruitment, 2, quantile, prob = 2.5/100)
uci <- apply(predicted_recruitment, 2, quantile, prob = 97.5/100)
ord <- order(alpha5$lag_wsst)
psi5_df <- data.frame(lag_wsst = alpha5$lag_wsst[ord],
                      recruitment = mean_recruitment[ord],
                      lci = lci[ord],
                      uci = uci[ord], 
                      lag_wsst_n = cov2$lag_wsst[ord])

psi5_df2 <- psi5_df %>% slice(1:11, 16:26)
alpha52 <- alpha5 %>% slice(5:26)

# plot logistic regression
(psi5_lagwsst <- psi5_df2 %>%
    ggplot() + 
    aes(x = lag_wsst_n, y = recruitment) + 
    geom_line(col = "#4B0055") + 
    geom_ribbon(aes(ymin = lci, ymax = uci), 
                fill = "#4B0055", 
                alpha = 0.5) + 
    geom_point(data = alpha52, aes(x = lag_wsst_n, y = mean), col = "gray40", size = 0.5) +
    ylim(0,0.5) + 
    labs(x = "Lag winter SST ºC", y = "Age 5 recruitment") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6, colour = "black"), 
          axis.text.y = element_text(size = 6, colour = "black", angle = 90), 
          axis.title = element_text(size = 8), 
          axis.ticks = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2))) 

# density plot of posterior
psi5_beta2$.variable <- recode(psi5_beta2$.variable, `beta[2]` = "Lag winter SST ºC")

(psi5_dens <- psi5_beta2 %>%  
    ggplot(aes(x = .value, y = forcats::fct_rev(.variable), fill = after_stat(x < 0))) +
    stat_halfeye(.width = c(0.50, 0.95), point_interval = "mean_hdi", alpha = 0.5, size = 0.5) +
    xlab("Estimate") + ylab("") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "gray50", linewidth = 0.5) +
    scale_fill_manual(values = c("gray80", "#4B0055"), guide = "none") +
    scale_y_discrete(expand = expansion(mult = c(0.4, 1.4))) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6, colour = "black"),
          axis.text.y = element_text(size = 7, colour = "black", angle = 90, hjust = -0.15),
          axis.title = element_text(size = 8), 
          axis.ticks.y = element_blank(), 
          axis.ticks.x = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2)))

#### recruitment age 6 - winter storms * winter sst ####
load("alpha6_wpc1_wsst_int.RData")

# beta 2 - storms
psi6_beta2 <- mcmc.alpha6.wpc1.wsst.int$samples %>%
  gather_draws(`beta[2]`)

psi6_beta2$.variable <- recode(psi6_beta2$.variable, `beta[2]` = "Storms")

# beta 3 and 4 - wsst and interaction
psi6_beta34 <- mcmc.alpha6.wpc1.wsst.int$samples %>%
  gather_draws(`beta[3]`, `beta[4]`)

psi6_beta34$.variable <- recode(psi6_beta34$.variable, `beta[3]` = "SST ºC", `beta[4]` = "Interaction")

# time-dependent estimates of recruitment age 6
alpha6 <- data.frame(MCMCsummary(mcmc.guille.agem8$samples, params = "alpha.6"))
alpha6 <- rename(alpha6, q2.5 = X2.5., q50 = X50., q97.5 = X97.5.)                    

alpha6$year <- seq(1995, 2020, 1)

# no 6 year old immatures in first 5 years
cov$wsst[1] <- 0
cov$wsst[2] <- 0
cov$wsst[3] <- 0
cov$wsst[4] <- 0
cov$wsst[5] <- 0

alpha6$wsst <- cov$wsst
alpha6$wsst_n <- cov2$wsst
alpha6$wpc1 <- cov$djf_pc1

# winter storms and wsst interaction plot
beta1 <- c(mcmc.alpha6.wpc1.wsst.int$samples$chain1[, 'beta[1]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain2[, 'beta[1]'],
           mcmc.alpha6.wpc1.wsst.int$samples$chain3[, 'beta[1]'])
beta2 <- c(mcmc.alpha6.wpc1.wsst.int$samples$chain1[, 'beta[2]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain2[, 'beta[2]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain3[, 'beta[2]'])
beta3 <- c(mcmc.alpha6.wpc1.wsst.int$samples$chain1[, 'beta[3]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain2[, 'beta[3]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain3[, 'beta[3]'])
beta4 <- c(mcmc.alpha6.wpc1.wsst.int$samples$chain1[, 'beta[4]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain2[, 'beta[4]'], 
           mcmc.alpha6.wpc1.wsst.int$samples$chain3[, 'beta[4]'])

# calculate values for plotting interaction
plus_sd <- mean(alpha6$wsst) + sd(alpha6$wsst) # wsst at 1 SD higher than mean
minus_sd <- mean(alpha6$wsst) - sd(alpha6$wsst) # wsst at 1 SD lower than mean

# define grid of values for winter storms 
predicted_survival <- matrix(NA, 
                             nrow = length(beta1), 
                             ncol = length(cov2$djf_pc1))

# predict mean recruitment for each MCMC iteration
for (i in 1:length(beta1)) {
  for (j in 1:length(cov2$djf_pc1)) {
    predicted_survival[i,j] <- plogis(beta1[i] + 
                                        beta2[i] * cov2$djf_pc1[j])
  }
}

# grid for recruitment +1 SD
predicted_plus <- matrix(NA, 
                         nrow = length(beta1), 
                         ncol = length(cov2$djf_pc1))

# predict +1 SD recruitment for each MCMC iteration
for (i in 1:length(beta1)) {
  for (j in 1:length(cov2$djf_pc1)) {
    predicted_plus[i,j] <- plogis(beta1[i] + 
                                    beta2[i] * cov2$djf_pc1[j] +
                                    beta3[i] * plus_sd +
                                    (beta4[i] * plus_sd * cov2$djf_pc1[j]))
  }
}

# grid for recruitment -1 SD
predicted_minus <- matrix(NA, 
                          nrow = length(beta1), 
                          ncol = length(cov2$djf_pc1))

# predict -1 SD recruitment for each MCMC iteration
for (i in 1:length(beta1)) {
  for (j in 1:length(cov2$djf_pc1)) {
    predicted_minus[i,j] <- plogis(beta1[i] + 
                                     beta2[i] * cov2$djf_pc1[j] +
                                     beta3[i] * minus_sd +
                                     (beta4[i] * minus_sd * cov2$djf_pc1[j]))
  }
}

# calculate posterior mean and credible interval
mean_plus <- apply(predicted_plus, 2, mean) # impact of winter storms on recruitment when wsst is 1 SD higher
mean_minus <- apply(predicted_minus, 2, mean) # impact of winter storms on recruitment when wsst is 1 SD lower
mean_survival <- apply(predicted_survival, 2, mean) # impact of winter storms on recruitment at mean wsst
lci <- apply(predicted_survival, 2, quantile, prob = 2.5/100) # mean lower CI
uci <- apply(predicted_survival, 2, quantile, prob = 97.5/100) # mean upper CI
plci <- apply(predicted_plus, 2, quantile, prob = 2.5/100) # plus 1 SD lower CI
puci <- apply(predicted_plus, 2, quantile, prob = 97.5/100) # plus 1 SD upper CI
mlci <- apply(predicted_minus, 2, quantile, prob = 2.5/100) # minus 1 SD lower CI
muci <- apply(predicted_minus, 2, quantile, prob = 97.5/100) # minus 1 SD Upper CI

ord <- order(cov2$djf_pc1)
dfi_p <- data.frame(c1 = cov2$djf_pc1[ord],                     # covariate - winter storms
                    survival = mean_survival[ord],               # mean recruitment
                    plus = mean_plus[ord],                       # recruitment +1 SD
                    minus = mean_minus[ord],                     # recruitment -1 SD
                    td = alpha6$mean[ord],                       # time-dependent recruitment estimates
                    lci = lci[ord],                              # mean lower CI
                    uci = uci[ord],                              # mean upper CI 
                    plci = plci[ord],                            # +1 SD lower CI
                    puci = puci[ord],                            # +1 SD upper CI
                    mlci = mlci[ord],                            # -1 SD lower CI
                    muci = muci[ord])                            # -1 SD upper CI

# rearrange data into long format to plot with legend
ymeans <- dfi_p %>% select(1:4)
ylower <- dfi_p %>% select(1, 6, 8, 10)
yupper <- dfi_p %>% select(1, 7 , 9, 11)

ymeanslong <- ymeans %>% pivot_longer(!c1, names_to = "group", values_to = "ys")
ymeanslong$group <- recode(ymeanslong$group, survival = "rmean")

ylowerlong <- ylower %>% pivot_longer(!c1, names_to = "group", values_to = "lci")
ylowerlong$group <- recode(ylowerlong$group, lci = "rmean", plci = "plus", mlci = "minus")

yupperlong <- yupper %>% pivot_longer(!c1, names_to = "group", values_to = "uci")
yupperlong$group <- recode(yupperlong$group, uci = "rmean", puci = "plus", muci = "minus")

y2 <- cbind(ymeanslong, ylowerlong, yupperlong)

y3 <- y2 %>% select(1,2,3,6,9)
colnames(y3)[1] <- "wpc1"

y3$group <- as.factor(y3$group)
levels(y3$group)

y3$group <- factor(y3$group, levels = c("rmean", "plus", "minus"))

alpha62 <- alpha6 %>% slice(6:26)
alpha62$group = "rmean"

# ggplot with legend
(a6i <- ggplot(data = y3, aes(x = wpc1, y = ys, fill = group, colour = group)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.4, colour = NA) + 
    geom_point(data = alpha62, aes(x = wpc1, y = mean), colour = "grey40", show.legend = F, size = 0.5) + 
    scale_fill_manual(values = c("#4B0055", "#009B95", "#FDE333"), labels = c("mean", "+1 SD", "-1 SD")) +
    scale_colour_manual(values = c("#4B0055", "#009B95", "#FDE333"),labels = c("mean", "+1 SD", "-1 SD")) +
    labs(x = "Winter storms", y = "Age 6 recruitment") +
    annotate(geom = "text", x = -3.7, y = 0.05, label = "less\nstormy", size = 2) +
    annotate(geom = "text", x = 1.9, y = 0.05, label = "more\nstormy", size = 2) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 6, colour = "black"), 
          axis.text.y = element_text(size = 6, colour = "black", angle = 90), 
          axis.title = element_text(size = 8), 
          axis.ticks = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2), 
          legend.title = element_blank(), 
          legend.position = "inside",
          legend.position.inside = c(0.5, 0.9), 
          legend.direction = "horizontal", 
          legend.text = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          legend.background = element_blank(), 
          legend.key.spacing = unit(0.2, "cm")))

# density plot
(psi6_dens2 <- ggplot() +
    stat_halfeye(data = psi6_beta2, aes(x = .value, y = .variable, fill = after_stat(x > 0)),
                 .width = c(0.50, 0.95), point_interval = "mean_hdi", alpha = 0.5, size = 0.5) +
    stat_halfeye(data = psi6_beta34, aes(x = .value, y = .variable, fill = after_stat(x < 0)),
                 .width = c(0.50, 0.95), point_interval = "mean_hdi", alpha = 0.5, size = 0.5) +xlab("Estimate") + ylab("") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "gray50", linewidth = 0.5) +
    scale_fill_manual(values = c("gray80", "#4B0055"), guide = "none") +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.6))) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6, colour = "black"),
          axis.text.y = element_text(size = 6, colour = "black", angle = 90, hjust = -0.15),
          axis.title = element_text(size = 8), 
          axis.ticks.y = element_blank(), 
          axis.ticks.x = element_line(size = 0.1), 
          panel.grid = element_line(size = 0.2)))

#### combine into one big plot ####
(fig <- plot_grid(apc1_plot, NA, ssst_plot, NA, psi5_lagwsst, NA, a6i, 
                  dens_apc1, NA, dens_ssst, NA, psi5_dens, NA, psi6_dens2, nrow = 2, ncol = 7,
                  rel_widths = c(1, -0.03, 1, -0.03, 1, -0.03, 1, 1, -0.01, 1, -0.01, 1, -0.01, 1), 
                  labels = c("A", NA, "B", NA, "C", NA, "D", "E", NA, "F", NA, "G", NA, "H"), 
                  label_size = 8, label_x = 0.02))


(fig_final <- ggdraw() +
    draw_plot(fig) +
    theme(panel.background = element_rect(fill = "white")) +
    panel_border(colour = "white"))

save_plot("model_results.pdf", fig_final, base_height = 4, base_width = 7, dpi = 600)
save_plot("model_results.png", fig_final, base_height = 4, base_width = 7, dpi = 600)
