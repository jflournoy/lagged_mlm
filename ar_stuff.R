#'---
#' author: John Flournoy
#' output: 
#'   html_document:
#'     toc: true
#'     toc_float: true
#'---

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
#' Including the residual AR stucture helps, but we still see a positive
#' estimate for the effect of x_lag. Since we're working with data that reflects
#' exaclty the data generating model, we should get exactly the same numbers that
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
                       na.action = na.omit,
                       control=lmeControl(opt = "optim")) 
      }
    } else {
      lag1.yAR <-lme(y ~ 1 + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
                     correlation = corAR1(form = ~ wave | id),
                     na.action = na.omit,
                     control=lmeControl(opt = "optim"))
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
      geom_smooth(method = 'gam', size = .75) + 
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
#' Interstingly, when we within-person center the lagged x variable, we see a negative bias. This is also
#' picked up in the p-values as our sample size increases.
#'
#' # Lagged y with within-person-centering
#' 
#' Note that when we within-person center the IV, x, we also within-person center the lagged DV, y, and include that (the 
#' between-person component of y is redundant with the random intercept. If you don't do this, wow, the false positive 
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
#' standard errors, or the DF being used. Perhaps another possibility is that the random effects are mispecified. I'm not
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
#' 
