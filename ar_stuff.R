#'---
#' title: "Do we need to control for the outcome's t-1 score?"
#' author: John Flournoy
#' output: 
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: hide
#'---

#' In 2018 I started a postdoc and was collaborating on decisions about how to model monthly data over the course of a year. The goal was to see if some variable 
#' measured the *previous* month ($x_{t-1}$) predicts the outcome *this* month ($y_{t}$). Having worked on several longitudinal datasets in graduate school, I had learned that it is crucial to 
#' control for the *previous* month's level of the outcome variable, or else risk spurious effects driven by non-causal correlations between the constructs being measured. I wanted to see if the way collaborators were proposing to run these analyses was really problematic, and I wanted to explore any other ways to correct for these kinds of confounds. To do so, I simulated data where x and y move together over time but neither causes the
#' other, and checked how often we conclude erroneously that $x_{t-1} \rightarrow y_{t}$. The short version: leave out last
#' month's y and you get a high false positive rate, and when people also differ in stable ways on these constructs,
#' some fixes that work in some situations don't fully work.

#+echo = T, warning = F, error = F, message = F
library(simsem)
library(nlme)
library(lmerTest)
library(lme4)
library(semPlot)
library(dplyr)
library(tidyr)
library(ggplot2)

N_seq = c(seq(30, 150, 40), seq(300, 1050, 150))

# Simulation results come from sim_riclpm.R (see "How the simulations work").
# LAGGED_RESULTS_DIR and LAGGED_EXPECTED_REPS exist only so a pilot run can be
# rendered without editing this file; the published page uses the defaults.
results_dir <- Sys.getenv("LAGGED_RESULTS_DIR", "riclpm_sims_r10000")
expected_reps <- as.integer(Sys.getenv("LAGGED_EXPECTED_REPS", "10000"))

# install.packages(c('simsem', 'nlme', 'lmerTest', 'lme4', 'semPlot', 'dplyr', 'tidyr', 'ggplot2'))

#+echo = T
if(is.na(expected_reps) || expected_reps < 1){
  stop('LAGGED_EXPECTED_REPS must be a positive integer; got "', Sys.getenv("LAGGED_EXPECTED_REPS"), '".')
}
if(!dir.exists(results_dir)){
  stop('Simulation results directory "', results_dir, '" does not exist. Run sim_riclpm.R ',
       '(Rscript sim_riclpm.R ', expected_reps, ' <cores>) or set LAGGED_RESULTS_DIR.')
}
summary_file <- file.path(results_dir, 'riclpm_sims_summary.RDS')
if(!file.exists(summary_file)){
  stop('Expected the simulation summary at ', summary_file, '; sim_riclpm.R writes it when every chunk is done.')
}
sim_sum <- as.data.frame(readRDS(summary_file))

needed_cols <- c('trait_var', 'spec', 'n', 'reps', 'failed', 'n_valid', 'reps_with_warnings',
                 'fp_rate', 'mc_se', 'mean_est', 'sd_est')
if(!all(needed_cols %in% names(sim_sum))){
  stop('Simulation summary is missing columns: ', paste(setdiff(needed_cols, names(sim_sum)), collapse = ', '))
}

spec_labels <- c(
  naive = 'y ~ x_lag (naive)',
  no_wcen_ar = 'y ~ x_lag, AR(1)',
  no_wcen_ylag = 'y ~ y_lag + x_lag',
  wcen_ar = 'centered x_lag, AR(1)',
  wcen_ylag = 'centered y_lag + x_lag (nlme)',
  wcen_ylag_l4 = 'centered y_lag + x_lag (lme4, Satterthwaite)',
  wcen_ylag_l4_wcenre = 'centered y_lag + x_lag (lme4, Satterthwaite, random slope)',
  riclpm = 'RI-CLPM')

# Every spec on all three generating models, except the RI-CLPM fit, which is
# not run on CLPM data (its trait variance would sit at the boundary of zero).
expected_cells <- rbind(
  expand.grid(trait_var = c(0, .5, 1), spec = setdiff(names(spec_labels), 'riclpm'), n = N_seq,
              stringsAsFactors = FALSE),
  expand.grid(trait_var = c(.5, 1), spec = 'riclpm', n = N_seq, stringsAsFactors = FALSE))
if(nrow(expected_cells) != 230) stop('Internal error: the expected design should have 230 cells, not ', nrow(expected_cells))

cell_key <- function(d) paste(d$trait_var, d$spec, d$n, sep = '/')
if(anyDuplicated(cell_key(sim_sum))){
  stop('Simulation summary has duplicated cells: ', paste(unique(cell_key(sim_sum)[duplicated(cell_key(sim_sum))]), collapse = ', '))
}
missing_cells <- setdiff(cell_key(expected_cells), cell_key(sim_sum))
extra_cells <- setdiff(cell_key(sim_sum), cell_key(expected_cells))
if(length(missing_cells) > 0 || length(extra_cells) > 0 || nrow(sim_sum) != nrow(expected_cells)){
  stop('Simulation summary has ', nrow(sim_sum), ' cells; expected ', nrow(expected_cells), '. Missing: ',
       paste(missing_cells, collapse = ', '), '. Unexpected: ', paste(extra_cells, collapse = ', '))
}
if(!identical(sort(unique(sim_sum$n)), sort(N_seq))){
  stop('Simulation summary sample sizes do not match N_seq.')
}
if(any(sim_sum$reps != expected_reps)){
  stop('Every simulation cell should have ', expected_reps, ' replications; these do not: ',
       paste(cell_key(sim_sum[sim_sum$reps != expected_reps, ]), collapse = ', '))
}
if(any(sim_sum$n_valid + sim_sum$failed != sim_sum$reps)){
  stop('Some simulation cells have replications that are neither valid nor recorded as failed: ',
       paste(cell_key(sim_sum[sim_sum$n_valid + sim_sum$failed != sim_sum$reps, ]), collapse = ', '))
}
if(anyNA(sim_sum[needed_cols])){
  stop('Simulation summary has missing values in: ', paste(needed_cols[sapply(sim_sum[needed_cols], anyNA)], collapse = ', '))
}
if(any(sim_sum$fp_rate < 0 | sim_sum$fp_rate > 1) || any(sim_sum$reps_with_warnings < 0 | sim_sum$reps_with_warnings > sim_sum$reps)){
  stop('Simulation summary has rates outside [0, 1] or warning counts outside [0, reps].')
}
if(any(sim_sum$n_valid < 2)){
  stop('Some simulation cells have fewer than 2 valid replications: ',
       paste(cell_key(sim_sum[sim_sum$n_valid < 2, ]), collapse = ', '))
}
trait_levels <- c('0 (CLPM)', '.5', '1')
trait_label <- function(tv){
  out <- c('0' = '0 (CLPM)', '0.5' = '.5', '1' = '1')[as.character(tv)]
  if(any(is.na(out))) stop('Unexpected trait variance: ', paste(tv[is.na(out)], collapse = ', '))
  unname(out)
}
trait_colors <- c('0 (CLPM)' = '#7a7a7a', '.5' = '#2a78d6', '1' = '#d4572a')
trait_shapes <- c('0 (CLPM)' = 16, '.5' = 17, '1' = 15)
fp_floor <- .01

# Leading-zero-free numbers, as in ".05".
f3 <- function(x, digits = 3) sub('^(-?)0\\.', '\\1.', formatC(round(x, digits) + 0, format = 'f', digits = digits))

sim_theme <- function(){
  theme_minimal(base_size = 12) +
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
n_scale <- function() scale_x_sqrt(breaks = N_seq, name = 'Sample size (N, square-root scale)')
trait_scales <- function(){
  list(scale_color_manual(values = trait_colors, limits = trait_levels, drop = FALSE, name = 'Trait variance'),
       scale_shape_manual(values = trait_shapes, limits = trait_levels, drop = FALSE, name = 'Trait variance'))
}

# One plot pair per model specification, from the summary: the false-positive
# rate on a log scale (reference line at .05, bars from a Wilson interval with
# z = 2, about +/-2 Monte Carlo SE), and
# the average estimate (reference line at 0, bars at +/-2 SE of the mean).
plot_spec <- function(spec_name){
  d <- sim_sum[sim_sum$spec == spec_name, ]
  if(nrow(d) == 0) stop('No simulation results for spec "', spec_name, '".')
  d$trait <- factor(trait_label(d$trait_var), levels = trait_levels)
  d$fp_plot <- pmax(d$fp_rate, fp_floor)
  # Wilson interval with z = 2, which (unlike rate +/- 2 SE) doesn't collapse
  # to a point at rates of 0 or 1.
  z <- 2
  centre <- (d$fp_rate + z^2 / (2 * d$n_valid)) / (1 + z^2 / d$n_valid)
  half <- z * sqrt(d$fp_rate * (1 - d$fp_rate) / d$n_valid + z^2 / (4 * d$n_valid^2)) / (1 + z^2 / d$n_valid)
  d$fp_lo <- pmax(centre - half, fp_floor)
  d$fp_hi <- pmax(pmin(centre + half, 1), fp_floor)
  d$est_se <- d$sd_est / sqrt(d$n_valid)
  dodge <- position_dodge(width = .8)
  n_floored <- sum(d$fp_rate < fp_floor)
  floor_note <- if(n_floored > 0) sprintf('%d rate(s) below %s drawn at %s', n_floored, f3(fp_floor, 2), f3(fp_floor, 2)) else NULL

  fp_plot <- ggplot(d, aes(x = n, y = fp_plot, color = trait, shape = trait)) +
    geom_hline(yintercept = .05, color = 'gray30', linetype = 'dashed') +
    geom_linerange(aes(ymin = fp_lo, ymax = fp_hi), position = dodge) +
    geom_line(position = dodge, linewidth = .5) +
    geom_point(position = dodge, size = 2, show.legend = TRUE) +
    scale_y_log10(limits = c(fp_floor, 1), breaks = c(.01, .025, .05, .1, .25, .5, 1),
                  labels = c('.01', '.025', '.05', '.1', '.25', '.5', '1')) +
    n_scale() + trait_scales() +
    labs(y = 'Proportion with p < .05 (log scale)', title = 'False-positive rate', caption = floor_note) +
    sim_theme()

  est_plot <- ggplot(d, aes(x = n, y = mean_est, color = trait, shape = trait)) +
    geom_hline(yintercept = 0, color = 'gray30', linetype = 'dashed') +
    geom_linerange(aes(ymin = mean_est - 2 * est_se, ymax = mean_est + 2 * est_se), position = dodge) +
    geom_line(position = dodge, linewidth = .5) +
    geom_point(position = dodge, size = 2, show.legend = TRUE) +
    n_scale() + trait_scales() +
    labs(y = 'Average estimate of lagged x (true value 0)', title = 'Average estimate') +
    sim_theme()

  print(fp_plot)
  print(est_plot)
}

# Lookups for numbers quoted in the text; each stops if a cell is missing.
cells <- function(spec, tv, n = N_seq){
  d <- sim_sum[sim_sum$spec %in% spec & sim_sum$trait_var %in% tv & sim_sum$n %in% n, ]
  if(nrow(d) != length(spec) * length(tv) * length(n)){
    stop('Expected ', length(spec) * length(tv) * length(n), ' cells for ', paste(spec, collapse = '+'),
         ' at trait variance ', paste(tv, collapse = '/'), '; found ', nrow(d))
  }
  d
}
avg_est <- function(spec, tv) mean(cells(spec, tv)$mean_est)
rng <- function(x, digits = 3) paste0(f3(min(x), digits), ' to ', f3(max(x), digits))
max_abs <- function(x, digits = 3) f3(max(abs(x)), digits)
# Cells more than two Monte Carlo SEs from .05, using the SE a rate of .05
# would have (the SE at the observed rate is 0 when the rate is 0 or 1).
null_se <- function(d) sqrt(.05 * .95 / d$n_valid)
n_hot <- function(spec, tv){ d <- cells(spec, tv); sum(d$fp_rate - .05 > 2 * null_se(d)) }
n_off <- function(spec, tv){ d <- cells(spec, tv); sum(abs(d$fp_rate - .05) > 2 * null_se(d)) }
n_cells <- function(spec, tv) nrow(cells(spec, tv))
max_n <- max(N_seq)


#+echo = T
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
#' The CLPM has no stable differences between people: everyone's x and y bounce around the same
#' average. To see what a stable trait does, I also generate data from the random-intercept CLPM
#' (RI-CLPM), which gives each person their own typical level of x and of y (trait variance .5 or
#' 1, uncorrelated), with the same within-person process as above (stability .7, within-wave
#' correlation .5, no cross-lagged effects). As written, `make_riclpm_lavaan()` leaves the observed
#' residual variances free, so `simulateData()` adds measurement error the model never mentions;
#' `sim_riclpm.R` fixes them to 0, so the data are exactly an RI-CLPM.
#'
#' # One simulated dataset
#'
#' First, a single dataset of 250 people drawn from the CLPM (with `empirical = T`, so its covariances
#' match the model exactly).
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
#' # Simulations
#' 
#' ## How the simulations work
#' 
#' Every simulation has a null cross-lag: over 12 waves, x and y are each stable (.7) and correlated
#' within each wave (.5), but neither affects the other. For each generating model (the CLPM, or
#' the RI-CLPM with trait variance .5 or 1) and each sample size (N = `r paste(N_seq, collapse = ', ')`),
#' `sim_riclpm.R` draws `r format(expected_reps, big.mark = ',')` datasets, fits each model below,
#' and records the estimate and _p_-value for lagged x. A fit that errors or doesn't converge counts
#' as a failure and is left out of that cell (the table at the end shows how many); fits that only
#' raise warnings are kept. With `r format(expected_reps, big.mark = ',')` replications, the Monte
#' Carlo SE of a rate near .05 is about `r f3(sqrt(.05 * .95 / expected_reps))`; the error bars in
#' the plots span about ±2 SE. The 2018 runs used `lmeControl(opt = "optim")`, which stopped short of the
#' maximum around N = 300 and made a spurious spike there for centered x with AR(1); these runs use
#' nlme's default optimizer, which fixes that.
#' 

#+echo = T, results = 'asis'
if(results_dir != 'riclpm_sims_r10000'){
  cat('\n\n**Note: this page was rendered from `', results_dir, '`, not the published results in `riclpm_sims_r10000`.**\n\n', sep = '')
}

#' The data generation and model fitting code is below. `sim_riclpm.R` reuses these definitions
#' (and the generators above), and adds the naive model and the RI-CLPM fit. Running it
#' (`Rscript sim_riclpm.R 10000 <cores>`) writes the summary this page reads,
#' `riclpm_sims_summary.RDS`, into the folder named in `results_dir` at the top.
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

#'
#' ## The naive model
#'
#' ```
#' lag1.naive <- lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
#'                   na.action = na.omit)
#' ```
#'
#' This is the first model in the demonstration above, with nothing to account for the stability
#' of y.
#'

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('naive')

#'
#' ## Using AR residual structure to account for stability
#'
#' ```
#' lag1.yAR <-lme(y ~ 1 + x_lag, random = ~ 1 | id, data = someData_long,
#'                correlation = corAR1(form = ~ wave | id),
#'                na.action = na.omit)
#' ```
#'

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('no_wcen_ar')

#'
#' What we see here on the CLPM data (gray) is a consistent bias in the estimates of the effect of lagged x on y. As our sample size
#' increases, we are better and better powered, which means we erroneously reject the null more often due
#' to this bias. In other words, if we are powered to detect an effect of lagged x on y, we are powered to
#' mistake this bias for a causal effect, even using the AR structure
#' to adjust for the stability of y. (The trait data behave differently; see "What's going on?" below.)
#'

#'
#' ## Using lagged DV to account for stability
#'
#' ```
#' lag1.yAR <-lme(y ~ 1 + y_lag + x_lag, random = ~ 1 | id, data = someData_long,
#'                na.action = na.omit)
#' ```
#' 

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('no_wcen_ylag')

#'
#' On the CLPM data, we see no bias in the estimate of the effect of lagged x on y, and correspondingly our error
#' rate is controlled appropriately. With stable traits, though, the estimate is biased negative, and the
#' error rate climbs with N.
#'
#'
#'
#' ## AR with within-person-centering
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

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('wcen_ar')

#'
#' Interestingly, when we within-person center the lagged x variable, we see a negative bias. This is also
#' picked up in the p-values as our sample size increases.
#'
#' ## Lagged y with within-person-centering
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

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('wcen_ylag')

#'
#' The error rate in this case is still higher than the nominal rate, but since bias is 0, this may be due to a problem with the
#' standard errors, or the DF being used. Perhaps another possibility is that the random effects are misspecified. I'm not
#' sure why this is, but we can test whether using the Satterthwaite correction, or including a random effect of the within-person
#' predictor helps.
#'
#'
#' ## Lagged y with within-person-centering, Satterthwaite correction
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

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('wcen_ylag_l4')

#'
#' Bummer. Still a higher than expected error rate.
#'
#' ## Lagged y with within-person-centering, Satterthwaite correction + WCEN RE
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

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('wcen_ylag_l4_wcenre')

#'
#' That looks like it helped a little bit, but it also seems to induce a small negative bias.
#' 
#' ## The RI-CLPM
#' 
#' ```
#' riclpm_fit_text <- riclpmr::riclpm_text(list(x = paste0("x", 1:12), y = paste0("y", 1:12)),
#'                                         constrain_over_waves = TRUE)
#' fit <- riclpmr::lavriclpm(riclpmModel = riclpm_fit_text, data = someData)
#' ```
#' 
#' Last, the RI-CLPM itself, fit to the wide data with riclpmr, with each path constrained equal
#' over waves. The estimate is the path from x at _t_ − 1 to y at _t_. It's the true model for the
#' trait data, so it's the benchmark here. It isn't fit to the CLPM data, where its trait variances
#' would sit at their boundary of zero.
#' 

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('riclpm')

#' 
#' ## The results in two tables
#' 

#+echo = T, results = 'asis'
table_n <- c(30, 300, max_n)
tab_cells <- sim_sum[sim_sum$n %in% table_n, ]
tab_cells$cell_fp <- f3(tab_cells$fp_rate)
tab_cells$cell_est <- f3(tab_cells$mean_est)
rate_wide <- tab_cells %>%
  mutate(col = paste0('rate, N = ', n)) %>%
  select(spec, trait_var, col, cell_fp) %>%
  spread(col, cell_fp)
est_wide <- tab_cells %>%
  mutate(col = paste0('estimate, N = ', n)) %>%
  select(spec, trait_var, col, cell_est) %>%
  spread(col, cell_est)
results_table <- dplyr::inner_join(rate_wide, est_wide, by = c('spec', 'trait_var'))
if(nrow(results_table) != nrow(expected_cells) / length(N_seq)) stop('Results table lost or gained rows in the join.')
results_table <- results_table %>%
  mutate(model = factor(spec_labels[spec], levels = spec_labels),
         `trait variance` = factor(trait_label(trait_var), levels = trait_levels)) %>%
  arrange(model, `trait variance`) %>%
  select(model, `trait variance`, paste0('rate, N = ', table_n), paste0('estimate, N = ', table_n))
print(knitr::kable(results_table, align = c('l', 'l', rep('r', 2 * length(table_n))),
                   caption = paste0('Proportion of replications with p < .05 for lagged x (nominal .05), and the average estimate ',
                                    '(true value 0), at three sample sizes. ', format(expected_reps, big.mark = ','),
                                    ' replications per cell.')))

problem_table <- sim_sum %>%
  group_by(spec, trait_var) %>%
  summarize(failed_total = sum(failed),
            failed_max = max(failed),
            failed_max_n = paste(n[failed == max(failed)], collapse = ', '),
            warn_total = sum(reps_with_warnings),
            warn_max = max(reps_with_warnings),
            warn_max_n = paste(n[reps_with_warnings == max(reps_with_warnings)], collapse = ', '),
            .groups = 'drop') %>%
  filter(failed_total > 0 | warn_total > 0) %>%
  mutate(model = factor(spec_labels[spec], levels = spec_labels),
         `trait variance` = factor(trait_label(trait_var), levels = trait_levels),
         `failed fits` = failed_total,
         `most failures in one cell (N)` = ifelse(failed_total > 0, paste0(failed_max, ' (', failed_max_n, ')'), ''),
         `fits with warnings` = warn_total,
         `most warnings in one cell (N)` = ifelse(warn_total > 0, paste0(warn_max, ' (', warn_max_n, ')'), '')) %>%
  arrange(model, `trait variance`) %>%
  select(model, `trait variance`, `failed fits`, `most failures in one cell (N)`, `fits with warnings`, `most warnings in one cell (N)`)
if(nrow(problem_table) == 0){
  cat('\n\nNo fit failed or raised a warning in any cell.\n\n')
} else {
  print(knitr::kable(problem_table, align = c('l', 'l', 'r', 'r', 'r', 'r'),
                     caption = paste0('Fits that failed (errors or non-convergence; left out of that cell) and fits that raised ',
                                      'at least one warning (kept), summed over the ', length(N_seq), ' sample sizes, out of ',
                                      format(expected_reps * length(N_seq), big.mark = ','),
                                      ' per row. Rows with neither are left out.')))
}

#+echo = T
# The text below states some results in words. Each such statement is checked
# here, so if new results contradict it the page stops instead of publishing
# prose that no longer matches the numbers.
check_claim <- function(ok, claim){
  if(!isTRUE(ok)) stop('The text claims ', claim, ', but these results do not show it. Revise the sentence that says so (in "What\'s going on?" or that model\'s section).')
}
pct_floor <- function(x) paste0(floor(100 * x), '%')
naive <- cells('naive', c(0, .5, 1))
ylag0 <- cells('no_wcen_ylag', 0)
ar0 <- cells('no_wcen_ar', 0)
wcen_ar_tr <- cells('wcen_ar', c(.5, 1))
wcen_ylag_all <- cells('wcen_ylag', c(0, .5, 1))
ric <- cells('riclpm', c(.5, 1))
ric_large <- ric[ric$n >= 150, ]
ric_30 <- ric[ric$n == 30, ]
rate_at <- function(spec, tv, n) cells(spec, tv, n)$fp_rate

check_claim(min(naive$fp_rate) > .9, 'the naive model rejects in nearly every replication')
check_claim(max(abs(ylag0$mean_est)) < .01, 'lagged y is unbiased on the CLPM data')
check_claim(avg_est('no_wcen_ylag', .5) < 0 && avg_est('no_wcen_ylag', 1) < 0, 'lagged y is biased negative with traits')
check_claim(all(rate_at('no_wcen_ylag', c(.5, 1), max_n) > rate_at('no_wcen_ylag', c(.5, 1), 30)),
            'the lagged y false-positive rate climbs with N when there are traits')
check_claim(avg_est('no_wcen_ar', 0) > 0, 'AR(1) without centering has a positive bias on the CLPM data')
check_claim(rate_at('no_wcen_ar', 0, max_n) > rate_at('no_wcen_ar', 0, 30), 'the AR(1) rate on CLPM data climbs with N')
check_claim(avg_est('wcen_ar', .5) < 0 && avg_est('wcen_ar', 1) < 0, 'centered x with AR(1) is biased negative with traits')
check_claim(max(abs(wcen_ylag_all$mean_est)) < .02, 'centered x with centered lagged y is close to unbiased')
check_claim(mean(wcen_ylag_all$fp_rate) > .05 && mean(wcen_ylag_all$fp_rate) < .10 && n_hot('wcen_ylag', c(0, .5, 1)) > 0,
            'centered x with centered lagged y runs a little hot (pooled rate between .05 and .10, some cells clearly above .05)')
# "Calibrated" here means every rate within .01 of .05, beyond Monte Carlo noise.
calibrated <- function(d) all(abs(d$fp_rate - .05) - 2 * null_se(d) < .01)
check_claim(calibrated(ylag0), 'lagged y controls the error rate on the CLPM data')
check_claim(calibrated(ric_large), 'the RI-CLPM error rate is close to nominal from N = 150 up')
check_claim(all(ric_30$fp_rate - .05 > 2 * null_se(ric_30)), 'the RI-CLPM runs hot at N = 30')
check_claim(sum(ric$reps_with_warnings[ric$n < 150]) >= .9 * sum(ric$reps_with_warnings),
            'nearly all RI-CLPM fit warnings are at N < 150')
check_claim(avg_est('wcen_ar', 0) < 0 && rate_at('wcen_ar', 0, max_n) > rate_at('wcen_ar', 0, 30),
            'centered x with AR(1) is biased negative on the CLPM data, with the rate rising with N (2018 remark)')
check_claim(mean(cells('wcen_ylag_l4', c(0, .5, 1))$fp_rate) > .05, 'the Satterthwaite version still runs above .05 (2018 remark: "Bummer")')
check_claim(n_hot('wcen_ylag_l4_wcenre', c(0, .5, 1)) <= n_hot('wcen_ylag_l4', c(0, .5, 1)) &&
              mean(cells('wcen_ylag_l4_wcenre', c(0, .5, 1))$mean_est) < 0,
            'the random slope helps a little but adds a small negative bias (2018 remark)')
check_claim(max(abs(ric$mean_est)) < .02, 'the RI-CLPM is unbiased')
check_claim(abs(avg_est('no_wcen_ar', .5)) < abs(avg_est('no_wcen_ar', 0)) && abs(avg_est('no_wcen_ar', 1)) < abs(avg_est('no_wcen_ar', 0)),
            'the AR(1) bias is smaller with a trait than on CLPM data')

#' 
#' # What's going on?
#' 
#' The naive model (y on lagged x, with a random intercept and nothing else) finds an effect that
#' isn't there nearly every time: _p_ < .05 in at least `r pct_floor(min(naive$fp_rate))` of
#' replications in every cell, with average estimates from `r rng(naive$mean_est)`. It's an
#' ordinary omitted-variable problem, where the omitted variable is the outcome one wave earlier:
#' $y_{t-1}$ causes $y_{t}$ and is correlated with $x_{t-1}$, so a regression of $y_{t}$ on
#' $x_{t-1}$ alone reads that path as an effect. Stable traits add a second path, through
#' between-person differences in people's typical levels of x and y.
#' 
#' Adding lagged y fixes it on the CLPM data (rates from `r rng(ylag0$fp_rate)`, average estimates
#' within `r max_abs(ylag0$mean_est)` of zero). But lagged y is the true model for CLPM data, so
#' that's its best case. With stable traits it's biased negative (averaged over N, 
#' `r f3(avg_est('no_wcen_ylag', .5))` at trait variance .5 and `r f3(avg_est('no_wcen_ylag', 1))` at 1).
#' The bias stays put while the standard error shrinks, so the rate climbs with N, to 
#' `r f3(rate_at('no_wcen_ylag', .5, max_n))` and `r f3(rate_at('no_wcen_ylag', 1, max_n))` at
#' N = `r max_n`. This is probably the well-known problem of a lagged outcome that is correlated with
#' the random intercept.
#' 
#' AR(1) residuals have the same kind of trouble. On the CLPM data the average estimate is about
#' `r f3(avg_est('no_wcen_ar', 0))`, so the rate goes from `r f3(rate_at('no_wcen_ar', 0, 30))` at
#' N = 30 to `r f3(rate_at('no_wcen_ar', 0, max_n))` at N = `r max_n`. Centering x and keeping the
#' AR(1) gives a negative bias when there's a trait (about `r f3(avg_est('wcen_ar', .5))` at trait
#' variance .5 and `r f3(avg_est('wcen_ar', 1))` at 1), with rates up to `r f3(max(wcen_ar_tr$fp_rate))`.
#' 
#' Centered x with centered lagged y stays close to unbiased everywhere (every average estimate
#' within `r max_abs(wcen_ylag_all$mean_est)` of zero) but runs a little hot: rates from
#' `r rng(wcen_ylag_all$fp_rate)`, with `r n_hot('wcen_ylag', c(0, .5, 1))` of its `r n_cells('wcen_ylag', c(0, .5, 1))` cells more
#' than two Monte Carlo SEs above .05. With Satterthwaite degrees of freedom that count is
#' `r n_hot('wcen_ylag_l4', c(0, .5, 1))`, and with a random slope as well it's
#' `r n_hot('wcen_ylag_l4_wcenre', c(0, .5, 1))`. The RI-CLPM is unbiased at every N (within
#' `r max_abs(ric$mean_est)` of zero), and from N = 150 up its rates are close to .05
#' (`r rng(ric_large$fp_rate)`). Below that it runs hot (`r rng(ric_30$fp_rate)` at N = 30), which is
#' about what you'd expect from a large-sample z-test with 24 observed variables and only a few dozen
#' people, and it's also where nearly all of its fit warnings are.
#' 
#' One thing I didn't expect: AR(1) without centering is less biased with a trait than without one.
#' Averaged over N, its estimate is `r f3(avg_est('no_wcen_ar', 0))` on the CLPM data,
#' `r f3(avg_est('no_wcen_ar', .5))` at trait variance .5, and `r f3(avg_est('no_wcen_ar', 1))` at 1.
#' I don't know why. My guess is that its positive bias on CLPM data is being offset by a negative
#' one from the trait, like the one that hits lagged y. I wouldn't count on it.
#' 
#' So what should you do? Separate within-person from between-person variation at the latent level,
#' with the RI-CLPM or with latent person-mean centering (as in DSEM), though only the RI-CLPM was
#' simulated here. Observed person-mean centering plus centered lagged y comes close, which may be
#' fine as long as you know it runs a little hot.
#' 
#' # Limits
#' 
#' There are two generating models and one set of parameters here (12 waves, stability .7, a
#' within-wave correlation of .5, equal and uncorrelated traits), so I can't say how any of this
#' scales. Only the null was simulated, so nothing here speaks to power.
#' 
