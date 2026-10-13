#'---
#' author: John Flournoy
#' output: 
#'   html_document:
#'     toc: true
#'     toc_float: true
#'---

#' In 2018 I was helping a lab fit multilevel models to monthly data, asking whether one variable
#' this month predicts another next month. I suspected the usual way of setting that up finds effects
#' that aren't there, so I simulated data where x and y move together over time but neither causes the
#' other, and checked how often different models say otherwise. The short version: leave out last
#' month's y and you get a false effect nearly every time, and when people also differ in stable ways,
#' the usual fixes don't fully work either. The details, and what I'd do instead, are below.

#+echo = T, warning = F, error = F, message = F
library(simsem)
library(nlme)
library(lmerTest)
library(lme4)
library(semPlot)
library(dplyr)
library(tidyr)

reps_per_n = 1e4
rerun = FALSE
N_seq = c(seq(30, 150, 40), seq(300, 1050, 150))
cores = 10

# install.packages(c('simsem', 'nlme', 'lmerTest', 'lme4', 'semPlot', 'dplyr', 'tidyr'))

#+echo = F, warning = F, error = F, message = F


make_clpm_lavaan <- function(generate = T, generating_params = list(s = .7,cl = 0,r = 1,co = .5)){
  if(generate){
    b <- lapply(generating_params, paste0, '*')
  } else {
    b <- as.list(rep('', length(generating_params)))
    names(b) <- names(generating_params)
  }
  
  var_at_wave <- lapply(c('x', 'y'), paste0, 1:12)
  
  stabilityX <- paste(unlist(
    lapply(2:length(var_at_wave[[1]]), 
           function(i){
             paste0(var_at_wave[[1]][i], ' ~ ', b$s, var_at_wave[[1]][i-1])
           })),
    collapse = '\n')
  
  stabilityY <- paste(unlist(
    lapply(2:length(var_at_wave[[1]]), 
           function(i){
             paste0(var_at_wave[[2]][i], ' ~ ', b$s, var_at_wave[[2]][i-1])
           })),
    collapse = '\n')
  
  
  residuals <- paste(unlist(
    lapply(var_at_wave, 
           function(thing){
             paste0(thing, ' ~~ ', b$r, thing)
           })),
    collapse = '\n')
  
  correlation <- paste(unlist(
    lapply(1:length(var_at_wave[[1]]), 
           function(i){
             paste0(var_at_wave[[1]][i], ' ~~ ', b$co, var_at_wave[[2]][i])
           })),
    collapse = '\n')
  
  regressionsX <- paste(unlist(
    lapply(2:length(var_at_wave[[1]]), 
           function(i){
             paste0(var_at_wave[[1]][i], ' ~ ', b$cl, var_at_wave[[2]][i-1])
           })),
    collapse = '\n')
  
  regressionsY <- paste(unlist(
    lapply(2:length(var_at_wave[[1]]), 
           function(i){
             paste0(var_at_wave[[2]][i], ' ~ ', b$cl, var_at_wave[[1]][i-1])
           })),
    collapse = '\n')
  
  
  generatingModel <- paste(stabilityX,
                           stabilityY,
                           residuals, 
                           correlation,
                           regressionsX,
                           regressionsY, sep = '\n')
  return(generatingModel)
}

#generate data based on RI-CLPM
make_riclpm_lavaan <- function(generate = T, generating_params = list(s = .7,cl = 0,r = 1,v = 1,cv = 0,co = .5)){
  if(generate){
    b <- lapply(generating_params, paste0, '*')
  } else {
    b <- as.list(rep('', length(generating_params)))
    names(b) <- names(generating_params)
  }
  
  latent_means <- paste(paste0(c('kappa =~ ', 'omega =~ '),
                               lapply(lapply(c('x', 'y'), paste0, 1:12), 
                                      function(vars){
                                        paste(paste0('1*', vars), collapse = ' + ')
                                      })), 
                        collapse = '\n')
  
  intercepts <- paste(unlist(lapply(list(c('x','mu'), c('y', 'pi')),
                                    function(pair){
                                      lapply(1:12, function(wave){
                                        paste0(pair[1], wave, ' ~ ', pair[2], wave, '*1')
                                      })
                                    })), collapse = '\n')
  
  latent_covar <- paste0('
kappa ~~ ',b$v,'kappa #variance
omega ~~ ',b$v,'omega #variance
kappa ~~ ',b$cv,'omega #covariance')
  
  latent_resids <- paste(unlist(lapply(list(c('p','x'), c('q', 'y')),
                                       function(pair){
                                         lapply(1:12, function(wave){
                                           paste0(pair[1], wave, ' =~ ', '1*', pair[2], wave)
                                         })
                                       })), collapse = '\n')
  
  regressions <- paste(
    unlist(lapply(list(c('p','q'), c('q', 'p')),
                  function(pair){
                    lapply(12:2, function(wave){
                      paste0(pair[1], wave, ' ~ ', b$s, pair[1], wave-1, ' + ', b$cl, pair[2], wave-1)
                    })
                  })), collapse = '\n')
  
  first_wave_varcovar <- paste0('
p1 ~~ p1
q1 ~~ q1
p1 ~~ ',b$co,'q1')
  
  residvar <- paste(unlist(lapply(c('p','q'),
                                  function(avar){
                                    lapply(2:12, function(wave){
                                      paste0(avar, wave, ' ~~ ', b$r, avar, wave)
                                    })
                                  })), collapse = '\n')
  
  contemporaneous <- paste(lapply(2:12, 
                                  function(wave){
                                    paste0('p', wave, ' ~~ ', b$co, 'q', wave)
                                  }), 
                           collapse = '\n')
  
  
  generatingModel <- paste(latent_means,
                           intercepts,
                           latent_covar,
                           latent_resids,
                           regressions,
                           first_wave_varcovar,
                           residvar,
                           contemporaneous, sep = '\n')
  return(generatingModel)
}

generatingModel <- make_clpm_lavaan()

#'
#' We want to sample data from a null model where _x_ and _y_ are correlated over time
#' (due to some shared influences), but are _not_ causally related. This is the null model
#' we're implicitly testing against, so we need to know that our procedure's error rate performance
#' is as expected (that is, that if we set $\alpha=.05$ we're actually controlling or error rate at
#' 5%). I've created a model that has _x_ and _y_ observed at 12 occasions, with x and y having high stability
#' (.7) and being somewhat correlated (.5).  
#' 
#'
#+echo = T
someData <- simulateData(model=generatingModel, sample.nobs=250, empirical=T,
                         auto.fix.first = FALSE, auto.var = FALSE, auto.fix.single = FALSE,
                         auto.cov.lv.x = FALSE, auto.cov.y = FALSE)
fit.DGM <- lavaan(make_clpm_lavaan(generate = F), someData, fixed.x=F)
semPaths(fit.DGM, what='est', rotation=2, exoCov=T, exoVar=T, layout = 'tree2',
         residuals = F, structural = F, sizeLat = 3, sizeMan = 3, sizeMan2 = 2, 
         intercepts = F, layoutSplit = F, bifactor = c('omega', 'kappa', paste0('q', 1:12)))
# summary(fit.DGM, standardize=T)

someData_long <- someData %>% 
  mutate(id = 1:n()) %>%
  gather(key, value, -id) %>%
  extract(key, c('var', 'wave'), '(x|y)(\\d+)') %>%
  spread(var, value) %>%
  mutate(wave = as.numeric(wave)) %>%
  arrange(id, wave)

someData_long$gmean_x <- mean(someData_long$x)
someData_long$gmean_y <- mean(someData_long$y)

someData_long <- someData_long %>%
  group_by(id) %>%
  mutate(x_lag = lag(x),
         y_lag = lag(y),
         id_mean_x = mean(x),
         wcen_x = x - id_mean_x,
         gcen_x = id_mean_x - gmean_x,
         wcen_x_lag = lag(wcen_x),
         id_mean_y = mean(y),
         wcen_y = y - id_mean_y,
         gcen_y = id_mean_y - gmean_y,
         wcen_y_lag = lag(wcen_y))

#'
#' Now we can see what happens if we predict $y_{t}$ from $x_{t-1}$ without accounting for previous
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
#' Can we use the residual autocorrelation structure in place of y_lag?
#'

lag1.yAR <-lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
               correlation = corAR1(form = ~ wave | id),
               na.action = na.omit)
summary(lag1.yAR) 

#'
#' Including the residual AR structure helps, but we still see a positive
#' estimate for the effect of x_lag. Since we're working with data that reflects
#' exactly the data generating model, we should get exactly the same numbers that
#' we put in -- so this positive estimate  is cause for concern. How
#' much does this inflate our false positive rate if we sample
#' randomly from the data generating model (as we do when we run a study)?
#' 

#'
#' If we're just estimating contemporaneous effects without allowing for possible lagged causes
#' what do we get?
#'

cont1.y <-lme(y ~ 1 + x, random = ~ 1 | id, data = someData_long,
              na.action = na.omit)
summary(cont1.y) 

cont1.yAR <-lme(y ~ 1 + x, random = ~ 1 | id, data = someData_long,
               correlation = corAR1(form = ~ wave | id),
               na.action = na.omit)
summary(cont1.yAR) 

#'
#' This looks right, 'cause this is very close to the generating model with no possible confounds.
#'

#' 
#' # Numerous simulations

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

run_model <- function(n = 250, gen_mod, wcen = FALSE, use_lag_y = FALSE, use_lme4 = FALSE, wcen_re = FALSE){
  someData_long <- generate_data(n = n, gen_mod = gen_mod, wcen = wcen)

  if(wcen){
    if(use_lag_y){
      if(use_lme4){
        if(wcen_re){
          lag1.yAR <-lmer(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x + (1 + wcen_x_lag || id), data = someData_long,
                          REML = TRUE,
                          na.action = na.omit)
        } else {
          lag1.yAR <-lmer(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x + (1 | id), data = someData_long,
                          REML = TRUE,
                          na.action = na.omit)
        }
      } else {
        lag1.yAR <-lme(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
                       na.action = na.omit)
      }
    } else {
      lag1.yAR <-lme(y ~ 1 + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
                     correlation = corAR1(form = ~ wave | id),
                     na.action = na.omit)
    }
    if(use_lme4){
      asum <- coef(summary(lag1.yAR, ddf = 'Satterthwaite'))
      x_lag_p <- asum['wcen_x_lag', 'Pr(>|t|)']
      x_lag_bias <- asum['wcen_x_lag', 'Estimate']
    } else {
      asum <- summary(lag1.yAR)
      x_lag_p <- asum$tTable['wcen_x_lag', 'p-value']
      x_lag_bias <- asum$tTable['wcen_x_lag', 'Value'] 
    }
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

if(rerun){
  library(parallel)
  system.time(no_wcen_ar_reps <- mclapply(N_seq,
                                       function(n) {
                                         somereps <- dplyr::bind_rows(
                                           replicate(reps_per_n,
                                                     run_model(n, generatingModel, wcen = F),
                                                     simplify = F))
                                         return(somereps)
                                       }, mc.cores = cores))
  saveRDS(no_wcen_ar_reps, 'no_wcen_ar.RDS')
  
  system.time(wcen_ar_reps <- mclapply(N_seq,
                                     function(n) {
                                       somereps <- dplyr::bind_rows(
                                         replicate(reps_per_n,
                                                   run_model(n, generatingModel, wcen = T),
                                                   simplify = F))
                                       return(somereps)
                                     }, mc.cores = cores))
  saveRDS(wcen_ar_reps, 'wcen_ar.RDS')
  
  system.time(no_wcen_ylag_reps <- mclapply(N_seq,
                                          function(n) {
                                            somereps <- dplyr::bind_rows(
                                              replicate(reps_per_n,
                                                        run_model(n, generatingModel, wcen = F, use_lag_y = T),
                                                        simplify = F))
                                            return(somereps)
                                          }, mc.cores = cores))
  saveRDS(no_wcen_ylag_reps, 'no_wcen_ylag_reps.RDS')
  
  system.time(wcen_ylag_reps <- mclapply(N_seq,
                                       function(n) {
                                         somereps <- dplyr::bind_rows(
                                           replicate(reps_per_n,
                                                     run_model(n, generatingModel, wcen = T, use_lag_y = T),
                                                     simplify = F))
                                         return(somereps)
                                       }, mc.cores = cores))
  saveRDS(wcen_ylag_reps, 'wcen_ylag_reps.RDS')
  
  system.time(wcen_ylag_l4_reps <- mclapply(N_seq,
                                            function(n) {
                                              somereps <- dplyr::bind_rows(
                                                replicate(reps_per_n,
                                                          run_model(n, generatingModel, wcen = T, use_lag_y = T,
                                                                    use_lme4 = T, wcen_re = F),
                                                          simplify = F))
                                              return(somereps)
                                            }, mc.cores = cores))
  saveRDS(wcen_ylag_l4_reps, 'wcen_ylag_l4_reps.RDS')
  
  system.time(wcen_ylag_l4_wcenre_reps <- mclapply(N_seq,
                                                   function(n) {
                                                     somereps <- dplyr::bind_rows(
                                                       replicate(reps_per_n,
                                                                 run_model(n, generatingModel, wcen = T, use_lag_y = T,
                                                                           use_lme4 = T, wcen_re = T),
                                                                 simplify = F))
                                                     return(somereps)
                                                   }, mc.cores = cores))
  saveRDS(wcen_ylag_l4_wcenre_reps, 'wcen_ylag_l4_wcenre_reps.RDS')
} else {
  no_wcen_ar_reps <- readRDS('no_wcen_ar.RDS')
  wcen_ar_reps <- readRDS('wcen_ar.RDS')
  no_wcen_ylag_reps <- readRDS('no_wcen_ylag_reps.RDS')
  wcen_ylag_reps <- readRDS('wcen_ylag_reps.RDS')
  wcen_ylag_l4_reps <- readRDS('wcen_ylag_l4_reps.RDS')
  wcen_ylag_l4_wcenre_reps <- readRDS('wcen_ylag_l4_wcenre_reps.RDS')
}

#'
#' # Using AR residual structure to account for stability
#'
#' ```
#' lag1.yAR <-lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
#'                correlation = corAR1(form = ~ wave | id),
#'                na.action = na.omit)
#' ```
#'

#+echo = F
plot_results <- function(somerez, max_fp = .20, plot_points = FALSE, which = c('both', 'bias', 'error')){
  library(ggplot2)
  somerez <- lapply(somerez, function(arez){
    if(inherits(arez, what = 'try-error')){
      return(NULL)
    } else {
      return(arez)
    }
  })
  somerez_df <- dplyr::bind_rows(somerez) %>%
    group_by(n) %>%
    mutate(fp = as.numeric(x_lag_p < .05))
  
  somerez_df_prop_sum <- somerez_df %>%
    group_by(n) %>%
    summarize(fp = mean(fp),
           fp.se = sqrt((fp * (1-fp)) / n()),
           fp.u = fp+1.96*fp.se,
           fp.l = fp-1.96*fp.se)
  
  if(which[[1]] %in% c('both', 'bias')){
    biasplot <- ggplot(somerez_df,
                       aes(x = n, y = x_lag_bias)) + 
      geom_violin(aes(group = n)) + 
      geom_hline(yintercept = 0)
    if(plot_points){
      biasplot <- biasplot +
        geom_point(alpha = .2, position = position_jitter(width = 2))  
    } else {
      biasplot <- biasplot +
        geom_boxplot(alpha = 0, size = .5, aes(group = factor(n)), width = 25)
    }
    biasplot <- biasplot + 
      geom_smooth(method = 'gam', formula = y ~ s(x, k = 3, bs = 'tp', fx = FALSE), size = .75) + 
      labs(x = 'Sample size', y = 'Estimate of effect of lagged x\n(should be 0, on average)') + 
      theme_minimal()
    
    print(biasplot)
  }
  if(which[[1]] %in% c('both', 'error')){
    print(ggplot(somerez_df_prop_sum,
                 aes(x = n, y = fp)) +
            geom_errorbar(aes(ymin = fp.l, ymax = fp.u), width = 0) +
            geom_point() + 
            geom_hline(yintercept = .05) + 
            geom_smooth(method = 'gam', formula = y ~ s(x, k = 3, bs = 'tp', fx = FALSE),
                        se = FALSE) +
            coord_cartesian(ylim = c(0, max_fp)) +
            scale_x_continuous(breaks = N_seq) + 
            labs(x = 'Sample size', y = 'Proportion of false positives\n(should be < .05 for all N)') + 
            theme_minimal())
  }
}

#+echo = T
plot_results(no_wcen_ar_reps, max_fp = 1)

#'
#' What we see here is a consistent bias in the estimates of the effect of lagged x on y. As our sample size
#' increases, we are better and better powered, which means we erroneously reject the null more often due
#' to this bias. In other words, if we are powered to detect an effect of lagged x on y, we are powered to
#' mistake this bias for a causal effect, even using the AR structure
#' to adjust for the stability of y.
#'

#'
#' # Using lagged DV to account for stability
#'
#' ```
#' lag1.yAR <-lme(y ~ 1 + y_lag + x_lag, random = ~ 1 | id, data = someData_long,
#'                na.action = na.omit)
#' ```
#' 

plot_results(no_wcen_ylag_reps)

#'
#' Here we see no bias in the estimate of the effect of lagged x on y, and correspondingly our error
#' rate is controlled appropriately.
#'
#'
#'
#' # AR with within-person-centering
#' 
#' We might want to create separate IVs for the within-person versus the between-person effect. In this
#' example, I've created a variable that is just each person's mean across all waves, centered at the 
#' group mean (gcen) and a within-person variable that is each person's score at a wave minus their mean.
#' 
#' This is in line with the advice in:
#' 
#' >Wang, L. (Peggy), & Maxwell, S. E. (2015). On disaggregating between-person and within-person effects 
#' with longitudinal data using multilevel models. Psychological Methods, 20(1), 63–83. 
#' https://doi.org/10.1037/met0000030
#'
#' ```
#' lag1.yAR <-lme(y ~ 1 + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
#'                correlation = corAR1(form = ~ wave | id),
#'                na.action = na.omit)
#' ```
#'

plot_results(wcen_ar_reps)

#'
#' Interestingly, when we within-person center the lagged x variable, we see a negative bias. This is also
#' picked up in the p-values as our sample size increases.
#'
#' # Lagged y with within-person-centering
#' 
#' Note that when we within-person center the IV, x, we also within-person center the lagged DV, y, and include that (the 
#' between-person component of y is redundant with the random intercept). If you don't do this, wow, the false positive 
#' error rate goes way up. I'm not quite sure why this is (simulations not shown).
#' 
#' ```
#' lag1.yAR <-lme(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
#'                na.action = na.omit) 
#' ```
#'

plot_results(wcen_ylag_reps)

#'
#' The error rate in this case is still higher than the nominal rate, but since bias is 0, this may be due to a problem with the
#' standard errors, or the DF being used. Perhaps another possibility is that the random effects are misspecified. I'm not
#' sure why this is, but we can test whether using the Satterthwaite correction, or including a random effect of the within-person
#' predictor helps.
#'
#'
#' # Lagged y with within-person-centering, Satterthwaite correction
#'
#' Perhaps using the Satterthwaite correction will appropriately adjust the _p_-values.
#' 
#' ```
#' lag1.yAR <-lmer(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x + (1 | id), data = someData_long,
#'                REML = TRUE,
#'                na.action = na.omit)
#' summary(lag1.yAR, ddf = 'Satterthwaite')
#' ```
#' 

plot_results(wcen_ylag_l4_reps)

#'
#' Bummer. Still a higher than expected error rate.
#'
#' # Lagged y with within-person-centering, Satterthwaite correction + WCEN RE
#'
#' The last thing I'll try is to add the within-person centered lagged predictor as random effect.
#'
#' ```
#' lag1.yAR <-lmer(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x + (1 + wcen_x_lag || id), data = someData_long,
#'                REML = TRUE,
#'                na.action = na.omit)
#' summary(lag1.yAR, ddf = 'Satterthwaite')
#' ```
#' 

plot_results(wcen_ylag_l4_wcenre_reps)

#'
#' That looks like it helped a little bit, but it also seems to induce a small negative bias.
#' 
#' # What's going on?
#' 
#' The naive model (y on lagged x, with a random intercept and nothing else) finds an effect that
#' isn't there, nearly every time: the average estimate is about .28 on the CLPM data and .25 to .26 when
#' people also have stable traits, with p < .05 in at least 97% of replications at every N. It's an
#' ordinary omitted-variable problem, where the omitted variable is the outcome one wave earlier:
#' $y_{t-1}$ causes $y_{t}$ and is correlated with $x_{t-1}$, so a regression of $y_{t}$ on
#' $x_{t-1}$ alone reads that path as an effect. Stable traits add a second path, through
#' between-person differences in people's typical levels of x and y.
#' 
#' On the 2018 CLPM data, adding lagged y fixed it (false-positive rates of .049 to .055 where the
#' runs finished, average estimates within .001 of zero). AR(1) residuals left a constant bias of
#' about .03, so the rate climbed with N (.115 at N = 30, .721 at N = 750) as the standard error
#' shrank around an off-center estimate. Person-mean centering plus lagged y was unbiased but ran a
#' little hot (.074 to .083). Note well that lagged y is the true model for CLPM data, so that was
#' its best case.
#' 

#+echo = F, results = 'asis'
rep_sets <- list(
  'y ~ x_lag, AR(1)' = no_wcen_ar_reps,
  'y ~ y_lag + x_lag' = no_wcen_ylag_reps,
  'centered x_lag, AR(1)' = wcen_ar_reps,
  'centered y_lag + x_lag (nlme)' = wcen_ylag_reps,
  'centered y_lag + x_lag (lme4, Satterthwaite)' = wcen_ylag_l4_reps,
  'centered y_lag + x_lag (lme4, Satterthwaite, random slope)' = wcen_ylag_l4_wcenre_reps)

rep_summary <- dplyr::bind_rows(lapply(names(rep_sets), function(label){
  somerez <- rep_sets[[label]]
  if(length(somerez) != length(N_seq)){
    stop('Expected one result per element of N_seq for "', label, '"; got ', length(somerez))
  }
  failed <- sapply(somerez, inherits, what = 'try-error')
  dplyr::bind_rows(
    dplyr::bind_rows(somerez[!failed]) %>%
      group_by(n) %>%
      summarize(fp = mean(x_lag_p < .05), est = mean(x_lag_bias)) %>%
      mutate(model = label, failed = FALSE),
    data.frame(n = N_seq[failed], fp = rep(NA_real_, sum(failed)), est = rep(NA_real_, sum(failed)), model = rep(label, sum(failed)),
               failed = rep(TRUE, sum(failed)), stringsAsFactors = FALSE))
}))

wide_table <- function(value, digits){
  rep_summary %>%
    mutate(cell = ifelse(failed, 'failed', formatC(.data[[value]], format = 'f', digits = digits)),
           model = factor(model, levels = names(rep_sets))) %>%
    select(n, model, cell) %>%
    spread(model, cell)
}

print(knitr::kable(wide_table('fp', 3),
                   caption = 'CLPM data (2018 runs, 10,000 replications per cell): proportion of replications with p < .05 for lagged x (nominal rate .05).'))
print(knitr::kable(wide_table('est', 4),
                   caption = 'CLPM data (2018 runs): average estimate of the effect of lagged x (true value 0).'))

failed_cells <- rep_summary %>% filter(failed)
if(nrow(failed_cells) > 0){
  cat('\n\nCells where the 2018 simulation run failed (no results saved):\n\n')
  cat(paste0('- ', failed_cells$model, ', N = ', failed_cells$n, collapse = '\n'), '\n')
}

#'
#' ## Data with stable traits
#' 
#' To see what a stable trait does, I generated data from the RI-CLPM, with trait variance .5 or 1
#' and everything else as before, and also fit the RI-CLPM itself (with riclpmr).^[As written,
#' `make_riclpm_lavaan()` leaves the observed residual variances free, so `simulateData()` adds
#' measurement error the model never mentions; `sim_riclpm.R` fixes them to 0, so the data are
#' exactly an RI-CLPM.] With 1,000 replications per cell, the Monte Carlo SE of a rate near .05 is
#' about .007.
#' 
#' Lagged y without centering is now biased negative (about −.03 at trait variance .5, −.04 at 1),
#' and by N = 750 it rejects in nearly every replication, which is probably the well-known problem of a
#' lagged outcome that is correlated with the random intercept. Centered x with AR(1) is biased by
#' about −.02, with rates mostly rising with N, to about .55.^[The N = 300 cells are worse (.54 and .64). A small
#' side check (60 replications at each of three N) suggests `lmeControl(opt = "optim")` stopping
#' short of the maximum there; the default optimizer gave about −.023. I can't rule out milder
#' versions of this at other N.] Centered x with centered lagged y stays unbiased and a little hot
#' (.063 to .097 without a random slope; .053 to .085 with one, plus a small negative bias). The
#' RI-CLPM is unbiased (within .002 of zero) and close to nominal (.036 to .074, the only two cells
#' more than about two Monte Carlo SEs from .05), though at N = 30 its fits raised warnings in 19%
#' (trait variance .5) and 7% (trait variance 1) of replications (table below).
#' 
#' One surprise: AR(1) without centering does better with a trait than without one. Its average
#' estimate is about .008 at trait variance .5 (the rate creeps from .071 to .177) and within .001
#' of zero at 1 (rates of .058 to .073). I don't know why. My guess is that its positive bias on CLPM
#' data is being offset by a negative one from the trait, like the one that hits lagged y, and that
#' at trait variance 1 they roughly cancel. I wouldn't count on it.
#' 

#+echo = F, results = 'asis'
riclpm_summary_file <- 'riclpm_sims_r1000/riclpm_sims_summary.RDS'
if(!file.exists(riclpm_summary_file)){
  stop('Expected the RI-CLPM simulation summary at ', riclpm_summary_file, '; run sim_riclpm.R with 1000 reps to create it.')
}
riclpm_sum <- as.data.frame(readRDS(riclpm_summary_file))

spec_labels <- c(
  naive = 'y ~ x_lag (naive)',
  no_wcen_ar = 'y ~ x_lag, AR(1)',
  no_wcen_ylag = 'y ~ y_lag + x_lag',
  wcen_ar = 'centered x_lag, AR(1)',
  wcen_ylag = 'centered y_lag + x_lag (nlme)',
  wcen_ylag_l4 = 'centered y_lag + x_lag (lme4, Satterthwaite)',
  wcen_ylag_l4_wcenre = 'centered y_lag + x_lag (lme4, Satterthwaite, random slope)',
  riclpm = 'RI-CLPM')

expected_cells <- rbind(
  expand.grid(trait_var = 0, spec = 'naive', n = N_seq, stringsAsFactors = FALSE),
  expand.grid(trait_var = c(.5, 1), spec = names(spec_labels), n = N_seq, stringsAsFactors = FALSE))
missing_cells <- dplyr::anti_join(expected_cells, riclpm_sum, by = c('trait_var', 'spec', 'n'))
if(nrow(missing_cells) > 0 || nrow(riclpm_sum) != nrow(expected_cells)){
  stop('RI-CLPM simulation summary has ', nrow(riclpm_sum), ' cells; expected ', nrow(expected_cells),
       '. Missing: ', paste(apply(missing_cells, 1, paste, collapse = '/'), collapse = ', '))
}
if(any(riclpm_sum$reps != 1000)){
  stop('Every RI-CLPM simulation cell should have 1000 replications; these do not: ',
       paste(with(riclpm_sum[riclpm_sum$reps != 1000, ], paste(trait_var, spec, n, sep = '/')), collapse = ', '))
}
if(any(riclpm_sum$n_valid + riclpm_sum$failed != riclpm_sum$reps)){
  stop('Some RI-CLPM simulation cells have replications that are neither valid nor recorded as failed.')
}

table_n <- c(30, 150, 300, 1050)
riclpm_wide <- function(value, digits){
  riclpm_sum %>%
    filter(trait_var > 0, n %in% table_n) %>%
    mutate(cell = formatC(.data[[value]], format = 'f', digits = digits),
           col = factor(paste0('trait ', trait_var, ', N = ', n),
                        levels = paste0('trait ', rep(c(.5, 1), each = length(table_n)), ', N = ', table_n)),
           model = factor(spec_labels[spec], levels = spec_labels)) %>%
    select(model, col, cell) %>%
    spread(col, cell) %>%
    arrange(model)
}

print(knitr::kable(riclpm_wide('fp_rate', 3),
                   caption = paste0('RI-CLPM data (1,000 replications per cell): proportion of replications with p < .05 for lagged x (nominal rate .05). ',
                                    'The centered AR(1) cells at N = 300 are probably affected by the optimizer problem described in the footnote.')))
print(knitr::kable(riclpm_wide('mean_est', 3),
                   caption = 'RI-CLPM data: average estimate of the effect of lagged x (true value 0).'))

warn_cells <- riclpm_sum %>%
  filter(reps_with_warnings > 0) %>%
  mutate(model = factor(spec_labels[spec], levels = spec_labels)) %>%
  arrange(model, trait_var, n)
warn_table <- warn_cells %>%
  group_by(model, trait_var) %>%
  summarize(total = sum(reps_with_warnings),
            `by N` = paste0(n, ': ', reps_with_warnings, collapse = '; '),
            .groups = 'drop')
print(knitr::kable(warn_table %>% rename(`trait variance` = trait_var),
                   caption = 'Replications (out of 1,000 per cell) that raised at least one warning; their estimates are kept in the rates above.'))

failed_new <- riclpm_sum %>% filter(failed > 0)
if(nrow(failed_new) > 0){
  cat('\n\nReplications that failed to fit (left out of the rates above):\n\n')
  cat(paste0('- ', spec_labels[failed_new$spec], ', trait variance ', failed_new$trait_var,
             ', N = ', failed_new$n, ': ', failed_new$failed, collapse = '\n'), '\n\n')
}

#+echo = F, fig.width = 9, fig.height = 5.5
library(ggplot2)
clpm_2018_fp <- rep_summary %>%
  filter(!failed) %>%
  mutate(spec = names(spec_labels)[match(model, spec_labels)], trait_var = 0) %>%
  select(trait_var, spec, n, fp_rate = fp)
if(any(is.na(clpm_2018_fp$spec))) stop('A 2018 model label does not match any label in spec_labels.')
fp_plot_df <- dplyr::bind_rows(clpm_2018_fp, riclpm_sum %>% select(trait_var, spec, n, fp_rate)) %>%
  mutate(model = factor(spec_labels[spec], levels = spec_labels),
         trait = factor(trait_var, levels = c(0, .5, 1), labels = c('0 (CLPM)', '.5', '1')))
print(ggplot(fp_plot_df, aes(x = n, y = fp_rate, color = trait, shape = trait)) +
        geom_hline(yintercept = .05, color = 'gray40') +
        geom_line(linewidth = .6) +
        geom_point(size = 2) +
        scale_y_log10(breaks = c(.025, .05, .1, .25, .5, 1), labels = c('.025', '.05', '.1', '.25', '.5', '1')) +
        scale_color_manual(values = c('#8a8a8a', '#2a78d6', '#eb6834')) +
        facet_wrap(~ model, ncol = 4, labeller = label_wrap_gen(width = 28)) +
        labs(x = 'Sample size', y = 'Proportion of false positives (log scale)',
             color = 'Trait variance', shape = 'Trait variance') +
        theme_minimal() +
        theme(legend.position = 'bottom'))

#'
#' So what should you do? Separate within-person from between-person variation at the latent level,
#' with the RI-CLPM or with latent person-mean centering (as in DSEM), though only the RI-CLPM was
#' simulated here. Observed person-mean centering plus centered lagged y comes close (unbiased, with
#' rates of roughly .05 to .10), which may be fine as long as you know it runs a little hot.
#' 
#' There's one generating model and one set of parameters here (12 waves, stability .7, a
#' within-wave correlation of .5, equal and uncorrelated traits), so I can't say how any of this
#' scales. Only the null was simulated, so nothing here speaks to power. And the new runs have
#' 1,000 replications per cell, not 10,000.
#' 
