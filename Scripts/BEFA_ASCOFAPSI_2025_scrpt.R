# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : BEFA_ASCOFAPSI_2025_scrpt.R                                ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 24-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Este script reproduce todos los análisis presentados durante el seminario
# "Introducción a la Psicometría Bayesiana: Aplicaciones al Análisis Factorial"

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

# Es necesario instalar Stan antes de utilizar blavaan!
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html

library(lavaan)
library(blavaan)
library(bayesplot)
library(semTools)
library(ggplot2)
library(cowplot)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Model specification and prior predictive check
# ─────────────────────────────────────────────────────────────────────────────

# Sintaxis del modelo nulo
null.model <- paste0("x", 1:9, "~~", "x", 1:9, collapse = "\n")

# Sintaxis del modelo: un factor común
model.2f <- ' visual_text =~ x1 + x2 + x3 + x4 + x5 + x6 
              speed =~ x7 + x8 + x9
'

# Sintaxis del modelo: dos factores relacionados
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

# Standardize observed variables
df <- HolzingerSwineford1939
df[,paste0("x", 1:9)] <- apply(df[,paste0("x", 1:9)], 2, scale)

# Prior predictive check
dpriors()

# Check prior distribution of lambda
curve(dnorm(x, mean = 0, sd = 10), from = -50, to = 50, 
      main = "Distribución a priori de Lambda",
      xlab = "Valores de Lambda", ylab = "Densidad de probabilidad", lwd = 2)

# Prior predictive check inside {blavaan}
# It will be warnings. Don't worry!
prior_predictive <- bcfa(model.2f, data = df, std.lv = TRUE, n.chains = 3,
                         meanstructure = TRUE, test = "none", sample = 20000, 
                         prisamp = TRUE, bcontrol = list(cores = 3)) 

# Save prior distribution
prior.draws <- blavInspect(prior_predictive, "mcmc")

# Prior predictive distribution: item intercepts
color_scheme_set("orange")
mcmc_hist(x = prior.draws, regex_pars = paste0("x", 1:6, "~1"), alpha = .5)

# Prior predictive distribution: factor loadings
color_scheme_set("blue")
mcmc_hist(x = prior.draws, regex_pars = "visual_text=~", alpha = .5)

# Prior predictive distribution: latent correlation matrix
color_scheme_set("teal")
mcmc_hist(x = prior.draws, pars = "visual_text~~speed", alpha = .5)

# Prior predictive distribution: residual variances
color_scheme_set("red")
mcmc_hist(x = prior.draws, pars = paste0("x", 1:6, "~~", "x", 1:6), alpha = .5)

# Change priors
new_priors <- dpriors(nu = "normal(0,1)", 
                      lambda = "normal(0,1)", 
                      rho = "beta(2,2)", 
                      psi = "gamma(1,1)")

# Sampling from the new prior
# It will be warnings. Don't worry!
new_prior_pred <- bcfa(model.2f, data = df, std.lv = TRUE, n.chains = 3,
                       meanstructure = TRUE, test = "none",  sample = 20000, 
                       dp = new_priors,
                       prisamp = TRUE, bcontrol = list(cores = 3)) 

# Save prior distribution
new.prior.draws <- blavInspect(new_prior_pred, "mcmc")

# Prior predictive distribution: item intercepts
color_scheme_set("orange")
mcmc_hist(x = new.prior.draws, regex_pars = paste0("x", 1:6, "~1"), alpha = .5)

# Prior predictive distribution: factor loadings
color_scheme_set("blue")
mcmc_hist(x = new.prior.draws, regex_pars = "visual_text=~", alpha = .5)

# Prior predictive distribution: latent correlation matrix
color_scheme_set("teal")
mcmc_hist(x = new.prior.draws, pars = "visual_text~~speed", alpha = .5)

# Prior predictive distribution: residual variances
color_scheme_set("red")
mcmc_hist(x = new.prior.draws, pars = paste0("x", 1:6, "~~", "x", 1:6), alpha = .5)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Bayesian model estimation
# ─────────────────────────────────────────────────────────────────────────────

# Estimate all the Bayesian models
blavaan.null.fit <- bcfa(null.model, data = df, burnin = 500, sample = 1000, 
                         n.chains = 3, meanstructure = TRUE,std.lv = TRUE, 
                         bcontrol = list(cores = 3))
blavaan.2f.fit <- bcfa(model.2f, data = df, burnin = 500, sample = 1000, 
                       n.chains = 3, meanstructure = TRUE,std.lv = TRUE, 
                       bcontrol = list(cores = 3))
blavaan.HS.fit <- bcfa(HS.model, data = df, burnin = 500, sample = 1000, 
                       n.chains = 3, meanstructure = TRUE,std.lv = TRUE, 
                       bcontrol = list(cores = 3))

# Fit also the frequentist CFA models with Maximum Likelihood
lavaan.null.fit <- cfa(null.model, data = df, std.lv = TRUE)
lavaan.2f.fit   <- cfa(model.2f,   data = df, std.lv = TRUE)
lavaan.HS.fit   <- cfa(HS.model,   data = df, std.lv = TRUE)

# Save blavaan models
# save(blavaan.null.fit, file = "blavaan objects/blavaan_null_fit.rdata")
# save(blavaan.2f.fit,   file = "blavaan objects/blavaan_2f_fit.rdata")
# save(blavaan.HS.fit,   file = "blavaan objects/blavaan_HS_fit.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Bayesian Convergence and Efficiency Assesment
# ─────────────────────────────────────────────────────────────────────────────

# Convergence: Potential Scale Reduction Factor (PSRF) must be below 1.05
which(blavInspect(blavaan.2f.fit, what = "rhat") > 1.05)
which(blavInspect(blavaan.HS.fit, what = "rhat") > 1.05)

# Efficiency: Effective Sample Sizes (at least 100 ESS per chain)
min(blavInspect(blavaan.2f.fit, "neff"))
min(blavInspect(blavaan.HS.fit, "neff"))

# Traceplots for parameters
posterior.2f <- blavInspect(blavaan.2f.fit, "mcmc")
posterior.HS <- blavInspect(blavaan.HS.fit, "mcmc")

color_scheme_set("mix-blue-red")
mcmc_trace(x = posterior.2f, regex_pars = "visual_text=~")
mcmc_trace(x = posterior.HS, regex_pars = c("visual=~", "textual=~"))

# If one chain goes wrong, you can highlight it
color_scheme_set("viridis")
mcmc_trace_highlight(x = posterior.HS, regex_pars = "visual~~textual", highlight = 3)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Bayesian Approximate Fit Measures
# ─────────────────────────────────────────────────────────────────────────────

# Bayesian fit indices: Two factors
(blav_fit_2f <- blavFitIndices(object = blavaan.2f.fit, baseline.model = blavaan.null.fit))

# Frequentist fit measures
fitmeasures(lavaan.2f.fit)[c("rmsea", "cfi", "tli")]

# Posterior distribution of bayesian fit indices
color_scheme_set("green")
mcmc_hist(data.frame(blav_fit_2f@indices), alpha = 0.5,
          pars = c("BRMSEA", "BCFI", "BTLI", "BGammaHat"))

# Bayesian fit indices: HS model
(blav_fit_HS <- blavFitIndices(object = blavaan.HS.fit, baseline.model = blavaan.null.fit))

# Frequentist fit measures
fitmeasures(lavaan.HS.fit)[c("rmsea", "cfi", "tli")]

# Posterior distribution of bayesian fit indices
color_scheme_set("brightblue")
mcmc_hist(data.frame(blav_fit_HS@indices), alpha = 0.5,
          pars = c("BRMSEA", "BCFI", "BTLI", "BGammaHat"))

# We can compute the posterior (empirical) cumulative density function of each statistic
tli_pdist <- ecdf(data.frame(blav_fit_HS@indices)[,"BTLI"])

# Probability that TLI value it's equal to or above 0.90
1 - tli_pdist(v = .90)

# We can compute the posterior (empirical) cumulative density function of each statistic
cfi_pdist <- ecdf(data.frame(blav_fit_HS@indices)[,"BCFI"])

# Probability that CFI value it's equal to or BELOW 0.90
cfi_pdist(v = .90)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Bayesian Model Comparison
# ─────────────────────────────────────────────────────────────────────────────

# Bayesian model comparison with several statistics
blav_com_2f_HS <- blavCompare(object1 = blavaan.2f.fit, 
                              object2 = blavaan.HS.fit)

# Watanabe-Akaike Information Criteria
blav_com_2f_HS$waic[[1]]    # Object 1 (i.e., blavaan.2f.fit) WAIC
blav_com_2f_HS$waic[[2]]    # Object 2 (i.e., blavaan.HS.fit) WAIC

# WAIC comparison: best model in the first row (here, model 2)
blav_com_2f_HS$diff_waic

# Leave-One-Out Cross-validation (better than WAIC)
blav_com_2f_HS$loo[[1]]    # Object 1 (i.e., blavaan.2f.fit) WAIC
blav_com_2f_HS$loo[[2]]    # Object 2 (i.e., blavaan.HS.fit) WAIC

# WAIC comparison: best model in the first row (here, model 2)
blav_com_2f_HS$diff_loo

# Log-bayes factor via Laplace-Metropolis Approximation
# Positive values favor object 1 (i.e., blavaan.2f.fit)
blav_com_2f_HS$bf # mll are marginal log-likelihoods

# Just exponentiate it see BF12
exp(blav_com_2f_HS$bf)[1] # Close-to zero.

# Now, BF21 is the inverse of BF12
1/exp(blav_com_2f_HS$bf)[1] 

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Plots and model summaries
# ─────────────────────────────────────────────────────────────────────────────

# Model 1 summary
summary(blavaan.2f.fit, standardized = TRUE, rsquare = TRUE)

# Model 2 summary
summary(blavaan.HS.fit, standardized = TRUE, rsquare = TRUE)

# Posterior distribution: histograms
color_scheme_set("orange")
mcmc_hist(x = posterior.HS, regex_pars = "~1", alpha = .5)

# Posterior distribution: density plots
color_scheme_set("teal")
mcmc_hist(x = posterior.HS, regex_pars = "textual~~|visual~~", alpha = .5)

# Posterior distribution: pairs plots
color_scheme_set("gray")
mcmc_pairs(posterior.HS, regex_pars = "textual=~", diag_fun = "hist", 
           off_diag_fun = "scatter", off_diag_args = list(alpha = 0.5))

# Posterior distribuition: uncertainty intervals
color_scheme_set("pink")
mcmc_intervals(posterior.HS, pars = paste0("x", 1:6, "~~", "x", 1:6))

# Posterior distribuition: uncertainty intervals with areas
color_scheme_set("blue")
mcmc_areas(posterior.HS, regex_pars = c("visual=~", "textual=~"),
           prob = 0.8, # 80% intervals
           prob_outer = 0.95, # 95%
           point_est = "median"
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Posterior Predictive Model Checks
# Simulate item responses 
# ─────────────────────────────────────────────────────────────────────────────

# Just simulate 100 datasets 
yrep <- sampleData(blavaan.HS.fit, nrep = 100, simplify = TRUE)

# Dataset with observed items
Y <- df[,7:15]

# Observed vs prior predicted distribution for one item
posterior_plots_list <- vector(mode = "list", length = ncol(Y))
for(i in 1:ncol(Y)){
  posterior_plots_list[[i]] <- ppc_dens_overlay(
    y = Y[,i], yrep = t(sapply(yrep, function(x) x[,i]))
  ) + 
    labs(title = paste("Item", i))
}

# All posterior plots
plot_grid(plotlist = posterior_plots_list)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Posterior Predictive Model Ckecks
# Comparing chi-square and SRMR
# ─────────────────────────────────────────────────────────────────────────────

# This will take a few minutes!!
posterior_checks_2f <- ppmc(object = blavaan.2f.fit, thin = 1, fit.measures = c("srmr","chisq"))
posterior_checks_HS <- ppmc(object = blavaan.HS.fit, thin = 1, fit.measures = c("srmr","chisq"))

# PPMC: Likelihood-Ratio test (i.e. chisq)
plot(posterior_checks_HS, element = "chisq")
plot(posterior_checks_2f, element = "chisq")
hist(posterior_checks_HS, element = "chisq")
hist(posterior_checks_2f, element = "chisq")

# PPMC: Standardized Root Mean Residuals
plot(posterior_checks_HS, element = "srmr")
plot(posterior_checks_2f, element = "srmr")
hist(posterior_checks_HS, element = "srmr")
hist(posterior_checks_2f, element = "srmr")

# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Posterior Predictive Model Checks
# Residual correlation matrix (i.e., observed - model-implied correlation)
# ─────────────────────────────────────────────────────────────────────────────

# Changing "discrepancy function": model-implied correlation-residual matrix
res_cors <- function(fit){ lavaan::lavResiduals(fit, zstat = FALSE, summary = FALSE)$cov } 
posterior_cor_residuals_2f <- ppmc(object = blavaan.2f.fit, thin = 1, discFUN = res_cors)
posterior_cor_residuals_HS <- ppmc(object = blavaan.HS.fit, thin = 1, discFUN = res_cors)

# Average residual correlation matrix
summary(posterior_cor_residuals_2f)$EAP
summary(posterior_cor_residuals_HS)$EAP

# How many residual correlations include zero inside their 95% CI?
lb_2f <- summary(posterior_cor_residuals_2f)$lower[lower.tri(summary(posterior_cor_residuals_2f)$lower)]
ub_2f <- summary(posterior_cor_residuals_2f)$upper[lower.tri(summary(posterior_cor_residuals_2f)$upper)]
sum(lb_2f <= 0 & 0 <= ub_2f) # Nº of correlations well-recovered
sum(lb_2f <= 0 & 0 <= ub_2f) / length(lb_2f) # proportion of correlations well-recovered

# How many residual correlations include zero inside their 95% CI?
lb_HS <- summary(posterior_cor_residuals_HS)$lower[lower.tri(summary(posterior_cor_residuals_HS)$lower)]
ub_HS <- summary(posterior_cor_residuals_HS)$upper[lower.tri(summary(posterior_cor_residuals_HS)$upper)]
sum(lb_HS <= 0 & 0 <= ub_HS) # Nº of correlations well-recovered
sum(lb_HS <= 0 & 0 <= ub_HS) / length(lb_HS) # proportion of correlations well-recovered

# Plot Posterior Predictive Checks for residual correlation:
plot(posterior_cor_residuals_2f, element = c("x1", "x2")) # Bad recovered
plot(posterior_cor_residuals_2f, element = c("x4", "x5")) # Well recovered

# PPMC: All correlations
unique_cors <- cbind(paste0("x", combn(1:9,2)[1,]), paste0("x", combn(1:9,2)[2,]))

# Two-factor model
par(mfrow=c(6,6))
for(i in 1:nrow(unique_cors)){
  plot(posterior_cor_residuals_2f, element = c(unique_cors[i,1], unique_cors[i,2]))
}
for(i in 1:nrow(unique_cors)){
  hist(posterior_cor_residuals_2f, element = c(unique_cors[i,1], unique_cors[i,2]))
}
par(mfrow=c(1,1))

# Three-factor model
par(mfrow=c(6,6))
for(i in 1:nrow(unique_cors)){
  plot(posterior_cor_residuals_HS, element = c(unique_cors[i,1], unique_cors[i,2]))
}
for(i in 1:nrow(unique_cors)){
  hist(posterior_cor_residuals_HS, element = c(unique_cors[i,1], unique_cors[i,2]))
}
par(mfrow=c(1,1))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 12: Posterior Predictive Model Checks
# Posterior distribution of reliability
# ─────────────────────────────────────────────────────────────────────────────

# Compote PPMC omega statistic
posterior_omega_2f <- ppmc(object = blavaan.2f.fit, thin = 1, discFUN = semTools::compRelSEM)
posterior_omega_HS <- ppmc(object = blavaan.HS.fit, thin = 1, discFUN = semTools::compRelSEM)

# Plot reliability PPMC: Two latent factors
par(mfrow=c(1,2))
plot(posterior_omega_2f, element = "visual_text")
plot(posterior_omega_2f, element = "speed")
par(mfrow=c(1,1))

# Plot reliability PPMC: Two latent factors
par(mfrow=c(1,2))
hist(posterior_omega_2f, element = "visual_text")
hist(posterior_omega_2f, element = "speed")
par(mfrow=c(1,1))

# Posterior reliability distribution: Two common factors
reliab.2f <- as.data.frame(do.call(rbind, posterior_omega_2f@obsDist$discFUN1))

# Posterior reliability distribution: Two common factors
par(mfrow=c(1,2))
hist(reliab.2f$visual_text, breaks = 17, 
     xlab = "Reliability: F1", ylab = "Density", 
     main = "Posterior distribution: reliability", 
     col = "peachpuff")
hist(reliab.2f$speed, breaks = 17, 
     xlab = "Reliability: F2", ylab = " ", 
     main = " ", 
     col = "lightblue1")
par(mfrow=c(1,1))

# Mean and standard deviation
colMeans(reliab.2f)
apply(reliab.2f, 2, sd)

# Posterior reliability distribution: Three common factors
reliab.HS <- do.call(rbind, posterior_omega_HS@obsDist$discFUN1)
par(mfrow=c(1,3))
hist(reliab.HS[,1], breaks = 30, 
     xlab = "Fiabilidad", ylab = "Densidad", 
     main = "Fiabilidad: visual", 
     col = "peachpuff")
hist(reliab.HS[,2], breaks = 30, 
     xlab = "Fiabilidad", ylab = "Densidad", 
     main = "Fiabilidad: textual", 
     col = "lightblue")
hist(reliab.HS[,3], breaks = 30, 
     xlab = "Fiabilidad", ylab = " ", 
     main = "Fiabilidad: speed", 
     col = "lightgreen")
par(mfrow=c(1,1))

# Mean and standard deviation
colMeans(reliab.HS)
apply(reliab.HS, 2, sd)

# ─────────────────────────────────────────────────────────────────────────────
