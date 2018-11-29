#'---
#' author: John Flournoy
#'---

#+echo = T, warning = F, error = F, message = F
library(simsem)
library(nlme)
library(lmerTest)
library(lme4)
library(semPlot)
library(dplyr)
library(tidyr)


#+echo = F, warning = F, error = F, message = F
var_at_wave <- lapply(c('x', 'y'), paste0, 1:12)

stabilityX <- paste(unlist(
  lapply(2:length(var_at_wave[[1]]), 
         function(i){
           paste0(var_at_wave[[1]][i], ' ~ .7*',var_at_wave[[1]][i-1])
         })),
  collapse = '\n')

stabilityY <- paste(unlist(
  lapply(2:length(var_at_wave[[1]]), 
         function(i){
           paste0(var_at_wave[[2]][i], ' ~ .7*', var_at_wave[[2]][i-1])
         })),
  collapse = '\n')


residuals <- paste(unlist(
  lapply(var_at_wave, 
         function(thing){
           paste0(thing, ' ~~ 1*', thing)
         })),
  collapse = '\n')

correlation <- paste(unlist(
  lapply(1:length(var_at_wave[[1]]), 
         function(i){
           paste0(var_at_wave[[1]][i], ' ~~ .3*', var_at_wave[[2]][i])
         })),
  collapse = '\n')

regressionsX <- paste(unlist(
  lapply(2:length(var_at_wave[[1]]), 
         function(i){
           paste0(var_at_wave[[1]][i], ' ~ .0*', var_at_wave[[2]][i-1])
         })),
  collapse = '\n')

regressionsY <- paste(unlist(
  lapply(2:length(var_at_wave[[1]]), 
         function(i){
           paste0(var_at_wave[[2]][i], ' ~ .0*', var_at_wave[[1]][i-1])
         })),
  collapse = '\n')


generatingModel <- paste(stabilityX,
                         stabilityY,
                          residuals, 
                          correlation,
                          regressionsX,
                          regressionsY, sep = '\n')

#'
#' We want to sample data from a null model where _x_ and _y_ are correlated over time
#' (due to some shared influences), but are _not_ causally related. This is the null model
#' we're implicitly testing against, so we need to know that our procedure's error rate performance
#' is as expected (that is, that if we set $\alpha=.05$ we're actually controlling or error rate at
#' 5%). I've created a model that has _x_ and _y_ observed at 12 occasions, with x and y having high stability
#' (.7) and being somewhat correlated (.3).  
#'
#+echo = T
someData <- simulateData(model=generatingModel, sample.nobs=250, empirical=T)
fit.DGM <- sem(generatingModel, someData, fixed.x=F)
semPaths(fit.DGM, what='est', rotation=2, exoCov=T, exoVar=T, layout = 'tree2', residuals = F)

someData_long <- someData %>% 
  mutate(id = 1:n()) %>%
  gather(key, value, -id) %>%
  extract(key, c('var', 'wave'), '(x|y)(\\d+)') %>%
  spread(var, value) %>%
  mutate(wave = as.numeric(wave)) %>%
  arrange(id, wave) %>%
  group_by(id) %>%
  mutate(x_lag = lag(x),
         y_lag = lag(y))

#'
#' Now we can see what happens if we predict y_{t} from x_{t-1} without accounting for previous
#' levels of y.
#' 

lag1.misspec <- lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
                 na.action = na.omit)
summary(lag1.misspec)

#' 
#' We see a very robust effect of lagged x on y, even though no causal effect exists, by construction.
#' This is due to the fact the x_lag is correlated with y_lag, and y_lag is correlated with y.
#' 

lag1.ylag <- lme(y ~ 1 + y_lag + x_lag, random = ~ 1 | id, data = someData_long,
                 na.action = na.omit)
summary(lag1.ylag)

#'
#' Once we control for lagged y, which is the true causal effect, we see that lagged x does not have an influence on y
#' 

lag1.ylag.xcor <- lme(y ~ 1 + y_lag + x_lag + x, random = ~ 1 | id, data = someData_long,
                 na.action = na.omit)
summary(lag1.ylag.xcor)

#'
#' Interestingly, when we add in contemporaneous values for x, we see the true 
#' correlation from the generating model, but we also see an induced negative
#' effect of x_lag on y. This is because we are not accounting for the fact that 
#' x and x_lag are correlated.
#' 
#' Can we use the residual autocorrelation in place of y_lag?
#'

lag1.yAR <-lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
               correlation = corAR1(form = ~ wave | id),
               na.action = na.omit)
summary(lag1.yAR) 

#'
#' Including the residual AR stucture helps, but we still see a positive
#' estimate for the effect of x_lag. Since we're working with data that reflects
#' exaclty the data generating model, this is cause for concern. How
#' much does this inflate our false positive rate if we sample
#' randomly from the data generating model?
#' 

generate_data <- function(n = 250, gen_mod, wcen = FALSE){
  someData <- simulateData(model=gen_mod, sample.nobs=n, empirical=F)
  
  someData_long <- someData %>% 
    mutate(id = 1:n()) %>%
    gather(key, value, -id) %>%
    extract(key, c('var', 'wave'), '(x|y)(\\d+)') %>%
    spread(var, value) %>%
    mutate(wave = as.numeric(wave)) %>%
    arrange(id, wave) %>%
    group_by(id) %>%
    mutate(x_lag = lag(x),
           y_lag = lag(y))
  
  if(wcen){
    someData_long$gmean_x <- mean(someData_long$x)
    someData_long$gmean_y <- mean(someData_long$y)
    someData_long <- someData_long %>%
      group_by(id) %>%
      mutate(id_mean_x = mean(x),
             wcen_x = x - id_mean_x,
             gcen_x = id_mean_x - gmean_x,
             wcen_x_lag = lag(wcen_x),
             id_mean_y = mean(y),
             wcen_y = y - id_mean_y,
             gcen_y = id_mean_y - gmean_y,
             wcen_y_lag = lag(wcen_y))
  }
  
  return(someData_long)
}

run_model <- function(n = 250, gen_mod, wcen = FALSE, use_lag_y = FALSE){
  someData_long <- generate_data(n = n, gen_mod = gen_mod, wcen = wcen)

  if(wcen){
    if(use_lag_y){
      lag1.yAR <-lme(y ~ 1 + wcen_y_lag + gcen_y + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
                     na.action = na.omit)  
    } else {
      lag1.yAR <-lme(y ~ 1 + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
                     correlation = corAR1(form = ~ wave | id),
                     na.action = na.omit)
    }
    
    asum <- summary(lag1.yAR)
    x_lag_p <- asum$tTable['wcen_x_lag', 'p-value']
    x_lag_bias <- asum$tTable['wcen_x_lag', 'Value']
  } else {
    if(use_lag_y){
      lag1.yAR <-lme(y ~ 1 + y_lag + x_lag, random = ~ 1 | id, data = someData_long,
                     na.action = na.omit)
    } else {
      lag1.yAR <-lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
                     correlation = corAR1(form = ~ wave | id),
                     na.action = na.omit)  
    }
    
    asum <- summary(lag1.yAR)
    x_lag_p <- asum$tTable['x_lag', 'p-value']
    x_lag_bias <- asum$tTable['x_lag', 'Value'] 
  }
  return(data.frame(n = n, x_lag_p = x_lag_p, x_lag_bias = x_lag_bias))
}

reps_per_n = 400
rerun = FALSE

if(rerun){
  library(parallel)
  system.time(no_wcen_ar_reps <- mclapply(seq(30, 250, 40),
                                       function(n) {
                                         somereps <- dplyr::bind_rows(
                                           replicate(reps_per_n,
                                                     run_model(n, generatingModel, wcen = F),
                                                     simplify = F))
                                         return(somereps)
                                       }, mc.cores = detectCores()))
  saveRDS(no_wcen_ar_reps, '~/code/lagged_mlm/no_wcen_ar.RDS')
  
  system.time(wcen_ar_reps <- mclapply(seq(30, 250, 40),
                                     function(n) {
                                       somereps <- dplyr::bind_rows(
                                         replicate(reps_per_n,
                                                   run_model(n, generatingModel, wcen = T),
                                                   simplify = F))
                                       return(somereps)
                                     }, mc.cores = detectCores()))
  saveRDS(wcen_ar_reps, '~/code/lagged_mlm/wcen_ar.RDS')
  
  system.time(no_wcen_ylag_reps <- mclapply(seq(30, 250, 40),
                                          function(n) {
                                            somereps <- dplyr::bind_rows(
                                              replicate(reps_per_n,
                                                        run_model(n, generatingModel, wcen = F, use_lag_y = T),
                                                        simplify = F))
                                            return(somereps)
                                          }, mc.cores = detectCores()))
  saveRDS(no_wcen_ylag_reps, '~/code/lagged_mlm/no_wcen_ylag_reps.RDS')
  
  system.time(wcen_ylag_reps <- mclapply(seq(30, 250, 40),
                                       function(n) {
                                         somereps <- dplyr::bind_rows(
                                           replicate(reps_per_n,
                                                     run_model(n, generatingModel, wcen = T, use_lag_y = T),
                                                     simplify = F))
                                         return(somereps)
                                       }, mc.cores = detectCores()))
  saveRDS(wcen_ylag_reps, '~/code/lagged_mlm/wcen_ylag_reps.RDS')
} else {
  no_wcen_ar_reps <- readRDS('~/code/lagged_mlm/no_wcen_ar.RDS')
  wcen_ar_reps <- readRDS('~/code/lagged_mlm/wcen_ar.RDS')
  no_wcen_ylag_reps <- readRDS('~/code/lagged_mlm/no_wcen_ylag_reps.RDS')
  wcen_ylag_reps <- readRDS('~/code/lagged_mlm/wcen_ylag_reps.RDS')
}

#'
#' # Using AR residual structure to account for stability
#'

library(ggplot2)
no_wcen_ar_reps_df <- dplyr::bind_rows(no_wcen_ar_reps) %>%
  group_by(n) %>%
  mutate(fp = as.numeric(x_lag_p < .05))

no_wcen_ar_reps_df_prop_sum <- no_wcen_ar_reps_df %>%
  group_by(n) %>%
  mutate(fp = mean(fp),
         fp.se = sqrt((fp * (1-fp)) / n()),
         fp.u = fp+1.96*fp.se,
         fp.l = fp-1.96*fp.se)

ggplot(no_wcen_ar_reps_df,
       aes(x = n, y = x_lag_bias)) + 
  geom_point(alpha = .4) +
  geom_smooth(method = 'gam') + 
  labs(x = 'Sample size', y = 'Estimate of effect of lagged x\n(should be 0, on average)') + 
  theme_minimal()

ggplot(no_wcen_ar_reps_df_prop_sum,
       aes(x = n, y = fp)) +
  geom_errorbar(aes(ymin = fp.l, ymax = fp.u), width = 0) +
  geom_point() + 
  geom_hline(yintercept = .05) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 3, bs = 'tp', fx = FALSE),
              se = FALSE) +
  coord_cartesian(ylim = c(0, .50)) +
  scale_x_continuous(breaks = seq(30, 250, 40)) + 
  labs(x = 'Sample size', y = 'Proportion of false positives\n(should be < .05 for all N)') + 
  theme_minimal()

#'
#' What we see here is a consistent bias in the estimates of the effect of lagged x on y. As our sample size
#' increases, we are better and better powered, which means we erroneously reject the null more often due
#' to this bias. In other words, if we are powered to detect an effect of lagged x and y, we are powered to
#' mistake this bias for a causal effect, even using the AR structure
#' to adjust for the stability of y.
#'

#'
#' # Using lagged DV to account for stability
#'


no_wcen_ylag_reps_df <- dplyr::bind_rows(no_wcen_ylag_reps) %>%
  group_by(n) %>%
  mutate(fp = as.numeric(x_lag_p < .05))

no_wcen_ylag_reps_df_prop_sum <- no_wcen_ylag_reps_df %>%
  group_by(n) %>%
  mutate(fp = mean(fp),
         fp.se = sqrt((fp * (1-fp)) / n()),
         fp.u = fp+1.96*fp.se,
         fp.l = fp-1.96*fp.se)

ggplot(no_wcen_ylag_reps_df,
       aes(x = n, y = x_lag_bias)) + 
  geom_point(alpha = .4) +
  geom_smooth(method = 'gam') + 
  labs(x = 'Sample size', y = 'Estimate of effect of lagged x\n(should be 0, on average)') + 
  theme_minimal()

ggplot(no_wcen_ylag_reps_df_prop_sum,
       aes(x = n, y = fp)) +
  geom_errorbar(aes(ymin = fp.l, ymax = fp.u), width = 0) +
  geom_point() + 
  geom_hline(yintercept = .05) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 3, bs = 'tp', fx = FALSE),
              se = FALSE) +
  coord_cartesian(ylim = c(0, .50)) +
  scale_x_continuous(breaks = seq(30, 250, 40)) + 
  labs(x = 'Sample size', y = 'Proportion of false positives\n(should be < .05 for all N)') + 
  theme_minimal()

#'
#' Here we see no bias in the estimate of the effect of lagged x on y, and correspondingly our error
#' rate is controlled appropriately.
#'

#'
#' # AR with within-person-centering
#'

wcen_ar_reps_df <- dplyr::bind_rows(wcen_ar_reps) %>%
  group_by(n) %>%
  mutate(fp = as.numeric(x_lag_p < .05))

wcen_ar_reps_df_prop_sum <- wcen_ar_reps_df %>%
  group_by(n) %>%
  mutate(fp = mean(fp),
         fp.se = sqrt((fp * (1-fp)) / n()),
         fp.u = fp+1.96*fp.se,
         fp.l = fp-1.96*fp.se)

ggplot(wcen_ar_reps_df,
       aes(x = n, y = x_lag_bias)) + 
  geom_point(alpha = .4) +
  geom_smooth(method = 'gam') + 
  labs(x = 'Sample size', y = 'Estimate of effect of lagged x\n(should be 0, on average)') + 
  theme_minimal()

ggplot(wcen_ar_reps_df_prop_sum,
       aes(x = n, y = fp)) +
  geom_errorbar(aes(ymin = fp.l, ymax = fp.u), width = 0) +
  geom_point() + 
  geom_hline(yintercept = .05) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 3, bs = 'tp', fx = FALSE),
              se = FALSE) +
  coord_cartesian(ylim = c(0, .50)) +
  scale_x_continuous(breaks = seq(30, 250, 40)) + 
  labs(x = 'Sample size', y = 'Proportion of false positives\n(should be < .05 for all N)') + 
  theme_minimal()

#'
#' Interstingly, when we within-person center the lagged x variable, we see a negative bias. This is also
#' picked up in the p-values as our sample size increases.
#'
#' # Lagged y with within-person-centering
#' 
#' Note that when we within-person cetner the IV, x, we also within-person center the lagged DV, y, and include
#' both the within- and between- person variables. If you don't do this, wow, does your error rate go up. 
#' The error rate still seems to be a little high, but since bias is 0, this may be due to a problem with the
#' standard errors, or the DF being used. Another possibility is that the random effects are mispecified. There is probably
#' an answer to this in the literature but I don't know what it is. I may see how this performs using
#' the Satterthwaite correction for degrees of freedom.
#'

wcen_ylag_reps_df <- dplyr::bind_rows(wcen_ylag_reps) %>%
  group_by(n) %>%
  mutate(fp = as.numeric(x_lag_p < .05))

wcen_ylag_reps_df_prop_sum <- wcen_ylag_reps_df %>%
  group_by(n) %>%
  mutate(fp = mean(fp),
         fp.se = sqrt((fp * (1-fp)) / n()),
         fp.u = fp+1.96*fp.se,
         fp.l = fp-1.96*fp.se)

ggplot(wcen_ylag_reps_df,
       aes(x = n, y = x_lag_bias)) + 
  geom_point(alpha = .4) +
  geom_smooth(method = 'gam') + 
  labs(x = 'Sample size', y = 'Estimate of effect of lagged x\n(should be 0, on average)') + 
  theme_minimal()

ggplot(wcen_ylag_reps_df_prop_sum,
       aes(x = n, y = fp)) +
  geom_errorbar(aes(ymin = fp.l, ymax = fp.u), width = 0) +
  geom_point() + 
  geom_hline(yintercept = .05) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 3, bs = 'tp', fx = FALSE),
              se = FALSE) +
  coord_cartesian(ylim = c(0, .50)) +
  scale_x_continuous(breaks = seq(30, 250, 40)) + 
  labs(x = 'Sample size', y = 'Proportion of false positives\n(should be < .05 for all N)') + 
  theme_minimal()
