#'---
#' title: "False positives in multilevel models for longitudinal data"
#' subtitle: "Do we need to control for the outcome's t-1 score? Lagged and concurrent effects, prior-outcome control, person-mean centering, AR(1) residuals, and the RI-CLPM, by simulation"
#' author: John Flournoy
#' output: 
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: hide
#'     includes:
#'       in_header: [false_positives_mlm_head.html, lightbox.html]
#' bibliography: refs.bib
#' link-citations: true
#'---

#' In 2018 I started a postdoc and was collaborating on decisions about how to model monthly data over the course of a year. The goal was to see if some variable 
#' measured the *previous* month ($x_{t-1}$) predicts the outcome *this* month ($y_{t}$). Having worked on several longitudinal datasets in graduate school, I had learned that it is crucial to 
#' control for the *previous* month's level of the outcome variable, or else risk spurious effects driven by non-causal correlations between the constructs being measured. I wanted to see if the way collaborators were proposing to run these analyses was really problematic, and I wanted to explore any other ways to correct for these kinds of confounds. To do so, I simulated data where x and y move together over time but neither causes the
#' other, and checked how often we conclude erroneously that $x_{t-1} \rightarrow y_{t}$. The short version: leave out last
#' month's y and you get a high false positive rate, and when people also differ in stable ways on these constructs,
#' some fixes that work in some situations don't fully work. None of these problems is new (see the
#' references); this page is a tutorial-style walk-through that puts them side by side for monthly data, with the code to reproduce each step. Varying the design (near the end) adds that
#' centering on observed person means brings a bias of its own when y is less stable than x and the two
#' move together within a wave.
#' 
#' _Keywords: `r paste(jsonlite::fromJSON('.zenodo.json')$keywords, collapse = ', ')`._
#' 
#' [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23383372.svg){data-external="1"}](https://doi.org/10.5281/zenodo.23383372) Cite as: Flournoy, J. C. (2026). _False positives in
#' multilevel models for longitudinal data_ (Version 1.0.0) [Computer software]. Zenodo.
#' <https://doi.org/10.5281/zenodo.23383372>
#' 
#' Related: [my RI-CLPM tutorial](https://johnflournoy.science/riclpm-lavaan-demo.html), which walks through
#' the model this page uses as its benchmark, and [riclpmr](https://johnflournoy.science/riclpmr/)
#' ([GitHub](https://github.com/jflournoy/riclpmr)), the R package that writes and fits the RI-CLPMs here.

#+echo = T
# The keywords live in .zenodo.json; the page's meta tags (false_positives_mlm_head.html)
# must list the same ones, so the two can't drift apart.
page_keywords <- jsonlite::fromJSON('.zenodo.json')$keywords
head_keywords <- sub('.*<meta name="keywords" content="([^"]*)".*', '\\1',
                     paste(readLines('false_positives_mlm_head.html'), collapse = ' '))
if (!identical(strsplit(head_keywords, ', ')[[1]], page_keywords)) {
  stop('false_positives_mlm_head.html keywords differ from .zenodo.json; update the meta tag to match.')
}

#+echo = T, warning = F, error = F, message = F
library(simsem)
library(nlme)
library(lmerTest)
library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)

N_seq = c(seq(30, 150, 40), seq(300, 1050, 150))

# Simulation results come from sim_riclpm.R (see "How the simulations work").
# LAGGED_RESULTS_DIR and LAGGED_EXPECTED_REPS exist only so a pilot run can be
# rendered without editing this file; the published page uses the defaults.
results_dir <- Sys.getenv("LAGGED_RESULTS_DIR", "riclpm_sims_r10000")
expected_reps <- as.integer(Sys.getenv("LAGGED_EXPECTED_REPS", "10000"))

# install.packages(c('simsem', 'nlme', 'lmerTest', 'lme4', 'dplyr', 'tidyr', 'ggplot2'))

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
  sea = 'centered x_lag + person mean + time (within/between)',
  no_wcen_ar = 'y ~ x_lag, AR(1)',
  no_wcen_ylag = 'y ~ y_lag + x_lag',
  wcen_ar = 'centered x_lag, AR(1)',
  wcen_ylag = 'centered y_lag + x_lag (nlme)',
  wcen_x_raw_ylag = 'raw y_lag + centered x_lag',
  wcen_ylag_l4 = 'centered y_lag + x_lag (lme4, Satterthwaite)',
  wcen_ylag_l4_wcenre = 'centered y_lag + x_lag (lme4, Satterthwaite, random slope)',
  riclpm = 'RI-CLPM')

# Every spec on all three generating models, except the RI-CLPM fit, which is
# not run on CLPM data (its trait variance would sit at the boundary of zero).
expected_cells <- rbind(
  expand.grid(trait_var = c(0, .5, 1), spec = setdiff(names(spec_labels), 'riclpm'), n = N_seq,
              stringsAsFactors = FALSE),
  expand.grid(trait_var = c(.5, 1), spec = 'riclpm', n = N_seq, stringsAsFactors = FALSE))
if(nrow(expected_cells) != 290) stop('Internal error: the expected design should have 290 cells, not ', nrow(expected_cells))

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
  list(scale_color_manual(values = trait_colors, limits = trait_levels, drop = FALSE, name = 'Data-generating model\ntrait variance'),
       scale_shape_manual(values = trait_shapes, limits = trait_levels, drop = FALSE, name = 'Data-generating model\ntrait variance'))
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

#' # The data-generating model
#' 
#' All the data here come from one model, drawn below from the same lavaan syntax the simulations
#' use. Each person's x and y at each wave is their usual level (a trait) plus a within-person part
#' (the starred circles); the within-person parts carry over from wave to wave and are correlated
#' within a wave. That's the RI-CLPM, and with no trait variance the traits drop out and it's the
#' CLPM. I think of it as the realistic case for monthly data: people differ in their usual levels,
#' and both x and y carry over from month to month and move together within a month. The first set of
#' simulations holds it at one setting (both stabilities .7, within-wave correlation .5, no cross-lag,
#' trait variance 0, .5 or 1 with uncorrelated traits); the factorial near the end varies the orange paths.
#' 

#+echo = T, fig.width = 21, fig.height = 7, out.width = '100%', fig.cap = "The data-generating model over 12 waves. Squares are observed scores, starred circles their within-person parts, and the ovals each person's traits. Orange paths are the ones the factorial varies; the key gives each path's value in the main simulations and its values in the factorial. With trait variance 0 the model is the CLPM."
source('sim_factorial_lib.R')
source('dgp_diagram.R')
dgp_diagram(12)

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
#' within each wave (.5), but neither affects the other. The data come from three data-generating
#' models, which are the three colors in every plot:
#'
#' - the CLPM (gray), with no stable differences between people (trait variance 0; the same
#'   generator as the 2018 runs);
#' - the RI-CLPM with trait variance .5 (blue);
#' - the RI-CLPM with trait variance 1 (orange).
#'
#' Every model below is fit to data from all three, except the RI-CLPM itself, which is only fit
#' to the RI-CLPM data (on CLPM data its trait variances would sit at zero). For each
#' data-generating model and each sample size (N = `r paste(N_seq, collapse = ', ')`),
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
#' (and the generators above), and adds the naive model, the within/between model, and the RI-CLPM fit. Running it
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
#' ## The within/between model
#'
#' ```
#' lag1.sea <- lme(y ~ 1 + wcen_x_lag + gcen_x + time_c, random = ~ 1 | id, data = someData_long,
#'                 na.action = na.omit)
#' ```
#'
#' This is a common way to set up month-level lagged analyses: x is split into each person's
#' deviation from their own mean (lagged one month) and their mean itself (grand-mean centered), with
#' a time term and a random intercept [@curran2011disaggregation; @wang2015disaggregating]. It's often fit in brms with weak priors; here it's fit with lme,
#' which should give essentially the same estimates. Like the naive model, it has nothing for the prior
#' month's y.
#'

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('sea')

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
#' error rate goes way up. I'm not quite sure why this is (see next sim).
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
#' ## Centered x with raw lagged y
#'
#' ```
#' lag1.rawy <- lme(y ~ 1 + y_lag + wcen_x_lag + gcen_x, random = ~ 1 | id, data = someData_long,
#'                  na.action = na.omit)
#' ```
#'
#' This is what happens if you center x but leave the lagged outcome raw.
#'

#+echo = T, fig.show = 'hold', out.width = '50%', fig.width = 5, fig.height = 4.5
plot_spec('wcen_x_raw_ylag')

#'
#' My best guess at why: each person's mean of x is taken over all 12 waves, including this month,
#' so centered lagged x carries a little of this month's shock to x (with a minus sign), and this
#' month's shock to x is correlated with this month's shock to y. Centered lagged y carries the same
#' kind of piece of this month's y, which seems to soak most of that up; raw lagged y doesn't, so
#' lagged x picks it up as a negative effect. I haven't checked this beyond these simulations.
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
#' Last, the RI-CLPM itself, fit to the wide data with [riclpmr](https://github.com/jflournoy/riclpmr),
#' with each path constrained equal over waves (for a walk-through of the model, see
#' [my RI-CLPM tutorial](https://johnflournoy.science/riclpm-lavaan-demo.html)). The estimate is the path from x at _t_ − 1 to y at _t_. It's the true model for the
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
         `data-generating trait variance` = factor(trait_label(trait_var), levels = trait_levels)) %>%
  arrange(model, `data-generating trait variance`) %>%
  select(model, `data-generating trait variance`, paste0('rate, N = ', table_n), paste0('estimate, N = ', table_n))
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
         `data-generating trait variance` = factor(trait_label(trait_var), levels = trait_levels),
         `failed fits` = failed_total,
         `most failures in one cell (N)` = ifelse(failed_total > 0, paste0(failed_max, ' (', failed_max_n, ')'), ''),
         `fits with warnings` = warn_total,
         `most warnings in one cell (N)` = ifelse(warn_total > 0, paste0(warn_max, ' (', warn_max_n, ')'), '')) %>%
  arrange(model, `data-generating trait variance`) %>%
  select(model, `data-generating trait variance`, `failed fits`, `most failures in one cell (N)`, `fits with warnings`, `most warnings in one cell (N)`)
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
rawy <- cells('wcen_x_raw_ylag', c(0, .5, 1))
check_claim(mean(rawy$fp_rate) > mean(cells('wcen_ylag', c(0, .5, 1))$fp_rate) + .2,
            'leaving lagged y raw when x is centered makes the error rate go way up (2018 remark)')
check_claim(all(rawy$mean_est < 0), 'with raw lagged y, centered lagged x picks up a negative effect')
check_claim(max(abs(ylag0$mean_est)) < .01, 'lagged y is unbiased on the CLPM data')
check_claim(avg_est('no_wcen_ylag', .5) < 0 && avg_est('no_wcen_ylag', 1) < 0, 'lagged y is biased negative with traits')
check_claim(all(rate_at('no_wcen_ylag', c(.5, 1), max_n) > rate_at('no_wcen_ylag', c(.5, 1), 30)),
            'the lagged y false-positive rate climbs with N when there are traits')
check_claim(avg_est('no_wcen_ar', 0) > 0, 'AR(1) without centering has a positive bias on the CLPM data')
check_claim(rate_at('no_wcen_ar', 0, max_n) > rate_at('no_wcen_ar', 0, 30), 'the AR(1) rate on CLPM data climbs with N')
check_claim(avg_est('wcen_ar', .5) < 0 && avg_est('wcen_ar', 1) < 0, 'centered x with AR(1) is biased negative with traits')
check_claim(max(abs(wcen_ylag_all$mean_est)) < .02, 'centered x with centered lagged y is close to unbiased')
check_claim(mean(wcen_ylag_all$fp_rate) > .05 && mean(wcen_ylag_all$fp_rate) < .10 && n_hot('wcen_ylag', c(0, .5, 1)) > 0,
            'centered x with centered lagged y has a somewhat inflated error rate (pooled rate between .05 and .10, some cells clearly above .05)')
# "Calibrated" here means every rate within .01 of .05, beyond Monte Carlo noise.
calibrated <- function(d) all(abs(d$fp_rate - .05) - 2 * null_se(d) < .01)
check_claim(calibrated(ylag0), 'lagged y controls the error rate on the CLPM data')
check_claim(calibrated(ric_large), 'the RI-CLPM error rate is close to nominal from N = 150 up')
check_claim(all(ric_30$fp_rate - .05 > 2 * null_se(ric_30)), 'the RI-CLPM error rate is inflated at N = 30')
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
#' between-person differences in people's typical levels of x and y [@hamaker2015critique].
#' 
#' When is leaving out $y_{t-1}$ not a problem? Write a person's departure from their own usual level
#' of y as $y_{t} = \phi y_{t-1} + e_{t}$. A model without $y_{t-1}$ then picks up $\phi$ times the
#' within-person covariance of $x_{t-1}$ and $y_{t-1}$. So it's only safe in two cases: y has no
#' within-person stability ($\phi = 0$, so month-to-month departures from a person's usual level are
#' just noise), or x and y aren't correlated within person at the same time point (a month that's high
#' on x for someone isn't also high on y). Neither seems likely for most things we'd measure monthly,
#' and if you want to lean on one, you'd need to show it rather than assume it. Without $y_{t-1}$, the
#' model still answers a real question (does last month's x go along with this month's y?), just not
#' whether x adds anything beyond y's own momentum. (Varying the design, below, adds caveats: that
#' reasoning is for a plain regression of within-person departures, and the random-intercept models
#' people actually fit can be biased even when y has no stability; and when x and y are both stable but
#' otherwise unconnected, the estimate is unbiased but the error rate is still inflated.)
#' 
#' Adding lagged y fixes it on the CLPM data (rates from `r rng(ylag0$fp_rate)`, average estimates
#' within `r max_abs(ylag0$mean_est)` of zero). But lagged y is the true model for CLPM data, so
#' that's its best case. With stable traits it's biased negative (averaged over N, 
#' `r f3(avg_est('no_wcen_ylag', .5))` at trait variance .5 and `r f3(avg_est('no_wcen_ylag', 1))` at 1).
#' The bias stays put while the standard error shrinks, so the rate climbs with N, to 
#' `r f3(rate_at('no_wcen_ylag', .5, max_n))` and `r f3(rate_at('no_wcen_ylag', 1, max_n))` at
#' N = `r max_n`. This is probably the well-known problem of a lagged outcome that is correlated with
#' the random intercept [@falkenstrom2017dynamic; @allison2017maximum; for the fixed-effects version,
#' @nickell1981biases].
#' 
#' AR(1) residuals have the same kind of trouble. On the CLPM data the average estimate is about
#' `r f3(avg_est('no_wcen_ar', 0))`, so the rate goes from `r f3(rate_at('no_wcen_ar', 0, 30))` at
#' N = 30 to `r f3(rate_at('no_wcen_ar', 0, max_n))` at N = `r max_n`. Centering x and keeping the
#' AR(1) gives a negative bias when there's a trait (about `r f3(avg_est('wcen_ar', .5))` at trait
#' variance .5 and `r f3(avg_est('wcen_ar', 1))` at 1), with rates up to `r f3(max(wcen_ar_tr$fp_rate))`.
#' 
#' Centered x with centered lagged y stays close to unbiased in these simulations (every average estimate
#' within `r max_abs(wcen_ylag_all$mean_est)` of zero) but has a somewhat inflated error rate: rates from
#' `r rng(wcen_ylag_all$fp_rate)`, with `r n_hot('wcen_ylag', c(0, .5, 1))` of its `r n_cells('wcen_ylag', c(0, .5, 1))` cells more
#' than two Monte Carlo SEs above .05. With Satterthwaite degrees of freedom that count is
#' `r n_hot('wcen_ylag_l4', c(0, .5, 1))`, and with a random slope as well it's
#' `r n_hot('wcen_ylag_l4_wcenre', c(0, .5, 1))`. The RI-CLPM is unbiased at every N (within
#' `r max_abs(ric$mean_est)` of zero), and from N = 150 up its rates are close to .05
#' (`r rng(ric_large$fp_rate)`). Below that its error rate is inflated (`r rng(ric_30$fp_rate)` at N = 30), which is
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
#' with the RI-CLPM [@hamaker2015critique] or with latent person-mean centering [as in DSEM,
#' @asparouhov2018dynamic; see also @ludtke2008multilevel], though only the RI-CLPM was
#' simulated here. Observed person-mean centering plus centered lagged y comes close when y is about as
#' stable as x, which may be fine as long as you know its error rate is somewhat inflated; when y is less stable than x
#' and the two are correlated within a wave, it's biased (see Varying the design).
#' 
#' # Varying the design
#' 
#+echo = T
fac_dir <- 'factorial_sims_r2500'
fac_reps <- 2500L
fac_file <- file.path(fac_dir, 'summary_cells.RDS')
if(!file.exists(fac_file)) stop('Expected the factorial summary at ', fac_file, '; run summarize_factorial.R on ', fac_dir, '.')
fac <- readRDS(fac_file)
fac_manifest <- jsonlite::fromJSON(file.path(fac_dir, 'manifest.json'))
if(fac_manifest$reps != fac_reps) stop('The factorial manifest says ', fac_manifest$reps, ' reps per cell, not ', fac_reps, '.')
if(!identical(unique(fac$partial), FALSE)) stop('The factorial summary is partial; rerun summarize_factorial.R without --partial.')
if(nrow(fac) != 3552 || any(fac$n_reps != fac_reps)) stop('The factorial summary should have 3552 rows of ', fac_reps, ' reps each.')
if(anyNA(fac[c('reject_rate', 'mean_est', 'sd_over_se')])) stop('The factorial summary has missing statistics.')
fac_conds <- unique(fac[c('n', 'waves', 'sy', 'co', 'clxy', 'v', 'tcor')])
fac$trait <- ifelse(fac$v == 0, 'none', ifelse(fac$tcor == 0, 'uncorrelated', 'correlated (.3)'))

fac_labels <- c(naive = 'y ~ x_lag (naive)', sea = 'centered x_lag + person mean + time (within/between)',
                no_wcen_ylag = 'y ~ y_lag + x_lag', wcen_x_raw_ylag = 'raw y_lag + centered x_lag',
                wcen_ylag = 'centered y_lag + x_lag', wcen_ylag_wcenre = 'centered y_lag + x_lag, random slope',
                no_wcen_ar = 'y ~ x_lag, AR(1)', wcen_ar = 'centered x_lag, AR(1)', riclpm = 'RI-CLPM',
                sea_concurrent = 'concurrent within/between model', sea_concurrent_ylag = '+ raw y_lag',
                sea_concurrent_wcen_ylag = '+ centered y_lag')

# Mean over the matching cells of a statistic for one spec and term role;
# stops if no cell matches, so a typo can't turn into an empty average.
fac_cells <- function(spec, role, ...){
  keep <- fac$spec %in% spec & fac$role == role
  for(cond in list(...)) keep <- keep & cond(fac)
  d <- fac[keep, ]
  if(nrow(d) == 0) stop('No factorial cells for ', paste(spec, collapse = '+'), ' / ', role)
  d
}
fac_mean <- function(stat, spec, role, ...) mean(fac_cells(spec, role, ...)[[stat]])
is_ <- function(col, val) function(d) d[[col]] %in% val

#' Everything above uses one set of parameters. To see how far it generalizes, `sim_factorial.R`
#' crosses y's stability (0, .3, .7; x's stays at .7), the within-wave correlation (0, .5), the
#' number of waves (6, 12), a true effect of lagged x on y (0, .1), and the trait variance (0, .5,
#' 1) with, when there are traits, a trait correlation (0, .3), at N = 50 and 150. That's
#' `r nrow(fac_conds)` cells with `r format(fac_reps, big.mark = ',')` replications each, and every
#' model is fit to the same simulated datasets within a cell. At the original settings and N = 150 (the
#' one sample size both runs have), it reproduces the 10,000-replication results above within Monte
#' Carlo error, model by model. The models are the
#' ones above, except that those fit with nlme without an AR structure are now fit with lme4
#' (Satterthwaite degrees of freedom), and there are three new ones for the concurrent question below.
#' 


#+echo = T

# Rates of null cells for the lagged models, one panel per model and one line
# per N x waves x trait condition (trait variance .5 and 1 averaged).
lag_specs <- c('naive', 'sea', 'no_wcen_ylag', 'wcen_x_raw_ylag', 'wcen_ylag', 'wcen_ylag_wcenre',
               'no_wcen_ar', 'wcen_ar', 'riclpm')
d_lag <- fac[fac$role == 'lag' & fac$null & fac$spec %in% lag_specs, ]
d_lag <- aggregate(reject_rate ~ spec + sy + co + trait + n + waves, data = d_lag, FUN = mean)
d_lag$N <- factor(d_lag$n, levels = c(50, 150))
d_lag$Waves <- factor(d_lag$waves, levels = c(6, 12))
d_lag$model <- factor(fac_labels[d_lag$spec], levels = fac_labels[lag_specs])
d_lag$co <- factor(paste('within-wave r =', d_lag$co))
d_lag$trait <- factor(d_lag$trait, levels = c('none', 'uncorrelated', 'correlated (.3)'))

#+echo = T, fig.width = 9, fig.height = 22, out.width = '85%'
ggplot(d_lag, aes(x = sy, y = pmax(reject_rate, fp_floor), color = trait, shape = trait, alpha = N, linetype = Waves,
                  group = interaction(trait, N, Waves))) +
  geom_hline(yintercept = .05, color = 'gray30', linetype = 'dotted') +
  geom_line(linewidth = .6) + geom_point(size = 1.8) +
  scale_alpha_manual(values = c('50' = .5, '150' = 1), name = 'N') +
  scale_linetype_manual(values = c('6' = 'dashed', '12' = 'solid'), name = 'Waves') +
  facet_grid(model ~ co, labeller = label_wrap_gen(22)) +
  scale_y_log10(limits = c(fp_floor, 1), breaks = c(.01, .05, .1, .25, .5, 1), labels = c('.01', '.05', '.1', '.25', '.5', '1')) +
  scale_x_continuous(breaks = c(0, .3, .7), labels = c('0', '.3', '.7'), name = "y's stability (x's is .7)") +
  scale_color_manual(values = c(none = '#7a7a7a', uncorrelated = '#2a78d6', 'correlated (.3)' = '#d4572a'), name = 'Traits') +
  scale_shape_manual(values = c(none = 16, uncorrelated = 17, 'correlated (.3)' = 15), name = 'Traits') +
  labs(y = 'Proportion with p < .05 (log scale)',
       caption = 'No true effect. With traits, each point averages trait variance .5 and 1.') +
  sim_theme() + theme(axis.text.x = element_text(angle = 0, hjust = .5),
                      strip.text.y = element_text(angle = 0, hjust = 0),
                      panel.spacing.y = unit(1.1, 'lines'), panel.spacing.x = unit(1.5, 'lines'),
                      legend.box = 'vertical')

#+echo = T
# Claims made in the text below.
both_stable <- function(d) d$sy == .7 & d$co == 0 & d$v == 0
check_claim(fac_mean('mean_est', c('naive', 'sea'), 'lag', is_('clxy', 0), both_stable) < .005 &&
              fac_mean('reject_rate', c('naive', 'sea'), 'lag', is_('clxy', 0), both_stable) > .1 &&
              fac_mean('reject_rate', c('no_wcen_ar', 'no_wcen_ylag'), 'lag', is_('clxy', 0), both_stable) < .06,
            'with both series stable and nothing else connecting them, naive and within/between are unbiased but reject > .1, and AR(1) and lagged y are near .05')
check_claim(fac_mean('reject_rate', 'naive', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0)) < .08 &&
              fac_mean('reject_rate', 'sea', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0)) > .5 &&
              fac_mean('mean_est', 'sea', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0)) < -.1,
            'with no y stability, the naive model is close to .05 but the within/between model is biased negative and rejects most of the time')
pair <- function(spec, tc) fac_mean('reject_rate', spec, 'lag', is_('clxy', 0), is_('co', 0), function(d) d$v > 0, is_('tcor', tc))
check_claim(all(sapply(c('naive', 'no_wcen_ylag', 'no_wcen_ar'), function(s) pair(s, .3) > pair(s, 0) + .04)) &&
              all(sapply(c('sea', 'wcen_x_raw_ylag', 'wcen_ylag', 'wcen_ylag_wcenre', 'wcen_ar', 'riclpm'), function(s) abs(pair(s, .3) - pair(s, 0)) < .01)),
            'a trait correlation raises the rate only for the models that do not center x')
wcy <- function(stat, s, w) fac_mean(stat, 'wcen_ylag', 'lag', is_('clxy', 0), is_('co', .5), is_('sy', s), is_('waves', w))
check_claim(wcy('mean_est', 0, 6) < -.1 && wcy('mean_est', 0, 12) < -.05 &&
              abs(wcy('mean_est', 0, 12)) < .6 * abs(wcy('mean_est', 0, 6)) &&
              abs(wcy('mean_est', .7, 6)) < .005 && abs(wcy('mean_est', .7, 12)) < .005,
            'centered lagged y is biased negative when y is less stable than x, the bias about halving from 6 to 12 waves, and unbiased at equal stabilities')
bs <- function(stat, w) fac_mean(stat, c('naive', 'sea'), 'lag', is_('clxy', 0), both_stable, is_('waves', w))
check_claim(bs('reject_rate', 12) > bs('reject_rate', 6) && bs('sd_over_se', 6) > 1.05 &&
              max(fac_cells(c('no_wcen_ar', 'no_wcen_ylag'), 'lag', is_('clxy', 0), both_stable)$reject_rate) < .06,
            'both stable: naive and within/between error rates are more inflated with 12 waves than 6, SEs too small; uncentered lagged y and AR(1) every cell below .06')
check_claim(fac_mean('mean_est', 'naive', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), function(d) d$v > 0) < -.05 &&
              fac_mean('reject_rate', 'naive', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), function(d) d$v > 0) > .3 &&
              all(fac_cells('sea', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0))$mean_est < -.05),
            'with no y stability, the naive model is biased once there are traits, and the within/between model is biased without them')
check_claim(max(abs(fac_cells('wcen_ylag', 'lag', is_('clxy', 0), is_('co', 0))$mean_est)) < .01,
            'centered x with centered lagged y is unbiased whenever there is no within-wave correlation')
tc_diff <- function(spec, co) {
  a <- fac_cells(spec, 'lag', is_('clxy', 0), is_('co', co), function(d) d$v > 0)
  key <- paste(a$n, a$waves, a$sy, a$v)
  r3 <- a$reject_rate[a$tcor == .3][order(key[a$tcor == .3])]
  r0 <- a$reject_rate[a$tcor == 0][order(key[a$tcor == 0])]
  r3 - r0
}
check_claim(all(sapply(c('sea', 'wcen_x_raw_ylag', 'wcen_ylag', 'wcen_ylag_wcenre', 'wcen_ar', 'riclpm'), function(s) max(abs(tc_diff(s, .5))) < .04)) &&
              any(tc_diff('no_wcen_ylag', .5) < -.1) && any(tc_diff('no_wcen_ar', .5) > .1) && any(tc_diff('no_wcen_ar', .5) < -.1),
            'with a within-wave correlation, the trait correlation leaves centered models and the RI-CLPM alone (within .04) but moves uncentered models both ways')
fac_old <- fac[fac$sy == .7 & fac$co == .5 & fac$waves == 12 & fac$clxy == 0 & fac$tcor == 0 & fac$n == 150 & fac$role == 'lag', ]
old_name <- c(naive = 'naive', sea = 'sea', no_wcen_ylag = 'no_wcen_ylag', wcen_x_raw_ylag = 'wcen_x_raw_ylag',
              wcen_ylag = 'wcen_ylag_l4', wcen_ylag_wcenre = 'wcen_ylag_l4_wcenre', no_wcen_ar = 'no_wcen_ar',
              wcen_ar = 'wcen_ar', riclpm = 'riclpm')
# Rates of exactly 1 in both runs have an SE of 0 (z is NaN); those are
# compared directly. Estimates are compared too.
repro <- do.call(rbind, mapply(function(s, v, r, se, est, sd, nv){
  o <- cells(old_name[[s]], v, 150)
  data.frame(z_rate = (r - o$fp_rate) / sqrt(se^2 + o$mc_se^2), both_one = r == 1 && o$fp_rate == 1,
             z_est = (est - o$mean_est) / sqrt(sd^2 / nv + o$sd_est^2 / o$n_valid))
}, fac_old$spec, fac_old$v, fac_old$reject_rate, fac_old$mc_se, fac_old$mean_est, fac_old$sd_est, fac_old$n_valid, SIMPLIFY = FALSE))
check_claim(nrow(fac_old) == 26 && all(ifelse(repro$both_one, TRUE, abs(repro$z_rate) < 3)) && all(abs(repro$z_est) < 3.5),
            'at the original settings and N = 150, the factorial reproduces the 10,000-replication rates and estimates within Monte Carlo error')

#' The first thing the factorial shows is that the "safe cases" above are only half safe. With y as
#' stable as x and nothing else connecting them (no within-wave correlation, no traits), the naive and
#' within/between models aren't biased (average estimates within
#' `r max_abs(fac_cells(c('naive', 'sea'), 'lag', is_('clxy', 0), both_stable)$mean_est)` of zero),
#' but they still reject `r f3(bs('reject_rate', 6))` of the time with 6 waves and
#' `r f3(bs('reject_rate', 12))` with 12. The estimates vary `r f3(bs('sd_over_se', 6), 2)` and
#' `r f3(bs('sd_over_se', 12), 2)` times as much as the standard errors say, probably because the
#' residuals of a model without $y_{t-1}$ are autocorrelated and the standard errors assume they aren't
#' [the old problem of relating two autocorrelated series, @bartlett1935some; @granger1974spurious].
#' Uncentered lagged y and AR(1) residuals are fine there. The other safe case, a y with no stability,
#' holds for the naive model only when people don't differ in their usual levels
#' (`r f3(fac_mean('reject_rate', 'naive', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0)))`
#' with no traits); with traits its average estimate is
#' `r f3(fac_mean('mean_est', 'naive', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), function(d) d$v > 0))`
#' and it rejects `r f3(fac_mean('reject_rate', 'naive', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), function(d) d$v > 0))`
#' of the time. The within/between model is biased even without traits: average estimates of
#' `r f3(fac_mean('mean_est', 'sea', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0), is_('waves', 6)))`
#' with 6 waves and
#' `r f3(fac_mean('mean_est', 'sea', 'lag', is_('clxy', 0), is_('sy', 0), is_('co', .5), is_('v', 0), is_('waves', 12)))`
#' with 12.
#' 
#' That negative bias seems to come from centering on the observed person mean (a random-intercept model
#' partly centers x too, which may be why the naive model picks it up once there are traits), and it also
#' shows up in centered x with centered lagged y once y is less stable than x and the two are correlated
#' within a wave. With no y stability, its average estimate is `r f3(wcy('mean_est', 0, 6))` with 6 waves and
#' `r f3(wcy('mean_est', 0, 12))` with 12, and the rate is
#' `r f3(wcy('reject_rate', 0, 6))` and `r f3(wcy('reject_rate', 0, 12))`. It's unbiased when y is as
#' stable as x (the only case the earlier simulations looked at) and whenever there's no within-wave
#' correlation. It looks like a cousin of the bias @nickell1981biases described for a lagged outcome with
#' fixed effects, which is negative, shrinks roughly as one over the number of waves, and also shows up in
#' person-mean-centered multilevel models [@hamaker2015center; @falkenstrom2017dynamic]. Fixed-effects
#' estimates are known to be biased whenever a predictor isn't strictly exogenous [@allison2017maximum],
#' and here lagged x isn't, because x moves with y within a wave. My guess at the mechanism: the person
#' mean includes this month's x, which moves with this month's y, so subtracting it pulls the lagged
#' coefficient down by something like one over the number of waves. That fits the bias about halving
#' from 6 to 12 waves (see the figure below). I don't know why it cancels when the stabilities are
#' equal, and y was never more stable than x here.
#' 
#' Without a within-wave correlation, a trait correlation hurts only the models that don't center x:
#' going from uncorrelated to correlated traits takes the rate for lagged y from `r f3(pair('no_wcen_ylag', 0))`
#' to `r f3(pair('no_wcen_ylag', .3))`, AR(1) from `r f3(pair('no_wcen_ar', 0))` to
#' `r f3(pair('no_wcen_ar', .3))`, and the naive model from `r f3(pair('naive', 0))` to
#' `r f3(pair('naive', .3))`, and the centered models and the RI-CLPM stay where they were. With a
#' within-wave correlation, it still leaves the centered models and the RI-CLPM alone, but it pushes the
#' uncentered models' rates around in both directions, sometimes down, by offsetting their other biases.
#' 

#+echo = T, fig.width = 8, fig.height = 4
d_bias <- fac[fac$role == 'lag' & fac$clxy == 0 & fac$co == .5 & fac$spec %in% c('sea', 'wcen_ylag', 'wcen_x_raw_ylag', 'wcen_ar'), ]
d_bias <- aggregate(mean_est ~ spec + sy + waves, data = d_bias, FUN = mean)
d_bias$model <- factor(fac_labels[d_bias$spec], levels = fac_labels[c('sea', 'wcen_ar', 'wcen_x_raw_ylag', 'wcen_ylag')])
d_bias$waves <- factor(d_bias$waves)
ggplot(d_bias, aes(x = sy, y = mean_est, color = waves, shape = waves)) +
  geom_hline(yintercept = 0, color = 'gray30', linetype = 'dashed') +
  geom_line(linewidth = .5) + geom_point(size = 2) +
  facet_wrap(~ model, nrow = 1, labeller = label_wrap_gen(18)) +
  scale_x_continuous(breaks = c(0, .3, .7), labels = c('0', '.3', '.7'), name = "y's stability (x's is .7)") +
  scale_color_manual(values = c('6' = '#2a78d6', '12' = '#d4572a'), name = 'Waves') +
  scale_shape_manual(values = c('6' = 16, '12' = 17), name = 'Waves') +
  scale_y_continuous(labels = function(x) f3(x, 1)) +
  labs(y = 'Average estimate (true value 0)',
       caption = 'Within-wave correlation .5, no true effect; averaged over N and the trait conditions.') +
  sim_theme() + theme(axis.text.x = element_text(angle = 0, hjust = .5))

#' ## The concurrent model
#' 
#+echo = T
conc_specs <- c('sea_concurrent', 'sea_concurrent_ylag', 'sea_concurrent_wcen_ylag')
cw <- function(stat, spec, ...) fac_mean(stat, spec, 'within', function(d) d$null, ...)
check_claim(max(abs(fac_cells(conc_specs, 'within', function(d) d$null)$mean_est)) < .01 &&
              cw('reject_rate', 'sea_concurrent', is_('sy', .7)) > .1 && cw('reject_rate', 'sea_concurrent', is_('sy', 0)) < .06 &&
              cw('reject_rate', 'sea_concurrent', is_('sy', .7), is_('waves', 12)) > cw('reject_rate', 'sea_concurrent', is_('sy', .7), is_('waves', 6)),
            'the concurrent within-person term is unbiased but its error rate is inflated as y stability rises, more so with more waves')
check_claim(abs(cw('reject_rate', 'sea_concurrent_ylag') - .05) < .01 &&
              cw('reject_rate', 'sea_concurrent_wcen_ylag', is_('sy', .7)) < cw('reject_rate', 'sea_concurrent', is_('sy', .7)) &&
              cw('reject_rate', 'sea_concurrent_wcen_ylag', is_('sy', 0)) > cw('reject_rate', 'sea_concurrent', is_('sy', 0)) &&
              min(fac_cells('sea_concurrent_wcen_ylag', 'within', function(d) d$null)$reject_rate) > .05,
            'raw lagged y fixes the concurrent term; centered lagged y helps at high y stability but leaves the error rate inflated everywhere, worse than the bare model at sy = 0')
check_claim(abs(cw('reject_rate', 'sea_concurrent', is_('tcor', .3), function(d) d$v > 0) -
                  cw('reject_rate', 'sea_concurrent', is_('tcor', 0), function(d) d$v > 0)) < .01 &&
              fac_mean('reject_rate', 'sea_concurrent', 'between', is_('co', 0), is_('clxy', 0), is_('tcor', .3), function(d) d$v > 0) > .2,
            'the trait correlation shows up in the between-person term, not the within-person one')
# RI-CLPM rejection rate under the .1 effect, for the next section.
rp <- function(n, w) fac_mean('reject_rate', 'riclpm', 'lag', is_('clxy', .1), is_('n', n), is_('waves', w))
# Person-mean (between) term with no within-wave correlation and no cross-lag,
# by trait correlation; its true value is 0 only when tcor = 0.
cb <- function(spec, tc) fac_mean('reject_rate', spec, 'between', is_('co', 0), is_('clxy', 0), function(d) d$v > 0, is_('tcor', tc))
check_claim(all(sapply(conc_specs, function(s) cb(s, .3) > .2)) && cb('sea_concurrent_ylag', .3) < cb('sea_concurrent', .3) - .03 &&
              max(abs(fac_cells(conc_specs, 'between', function(d) d$null)$reject_rate - .05)) < .025 &&
              all(abs(sapply(c(0, .3, .7), function(y) cw('reject_rate', 'sea_concurrent', is_('sy', y), is_('tcor', .3), function(d) d$v > 0) -
                                                     cw('reject_rate', 'sea_concurrent', is_('sy', y), is_('tcor', 0), function(d) d$v > 0))) < .02),
            'the person-mean term picks up the trait correlation (less so with raw lagged y), holds its error rate otherwise, and the within-person term is unaffected by it at every y stability')
# Both terms, no within-wave correlation and no cross-lag, every trait condition.
d_conc <- fac[fac$spec %in% conc_specs & fac$co == 0 & fac$clxy == 0, ]
d_conc <- aggregate(reject_rate ~ spec + role + sy + trait + n + waves, data = d_conc, FUN = mean)
d_conc$model <- factor(fac_labels[d_conc$spec], levels = fac_labels[conc_specs])
d_conc$term <- factor(ifelse(d_conc$role == 'within', 'within-person x\n(true value 0)', 'person mean of x\n(not 0 when traits correlate)'),
                      levels = c('within-person x\n(true value 0)', 'person mean of x\n(not 0 when traits correlate)'))
d_conc$trait <- factor(d_conc$trait, levels = c('none', 'uncorrelated', 'correlated (.3)'))
d_conc$N <- factor(d_conc$n, levels = c(50, 150))
d_conc$Waves <- factor(d_conc$waves, levels = c(6, 12))

#' The within/between model is often used concurrently, too: this month's y on this month's centered x,
#' the person mean of x, and time. In Sanjay Srivastava's lab we'd have added last month's y, to ask whether x goes
#' with y over and above where y already was. Under the null for the concurrent within-person term (no
#' within-wave correlation and no cross-lag), the estimate is unbiased, but the error rate is inflated as
#' y gets more stable, probably for the same reason as above: x and y are each autocorrelated, so the model
#' sees more independent information than there is. Adding raw lagged y
#' fixes it. Centered lagged y helps when y is quite stable but leaves the error rate somewhat inflated everywhere, even where
#' the bare model was fine (see the figure). The trait correlation shows up where it should, in the
#' person-mean term (which finds it `r f3(cb('sea_concurrent', .3))` of the time, averaged over N and waves),
#' and doesn't leak into the within-person one at any level of y's stability, including none
#' (`r f3(cw('reject_rate', 'sea_concurrent', is_('sy', 0), is_('tcor', .3), function(d) d$v > 0))` with
#' correlated traits and no y stability). Adding raw lagged y takes some of that between-person association
#' with it (`r f3(cb('sea_concurrent_ylag', .3))`), probably because last month's raw y carries each person's
#' usual level; the person-mean terms hold their error rate when the traits are uncorrelated.
#' 



#+echo = T, fig.width = 9, fig.height = 7
ggplot(d_conc, aes(x = sy, y = reject_rate, color = trait, shape = trait, alpha = N, linetype = Waves,
                   group = interaction(trait, N, Waves))) +
  geom_hline(yintercept = .05, color = 'gray30', linetype = 'dotted') +
  geom_line(linewidth = .6) + geom_point(size = 1.8) +
  facet_grid(term ~ model) +
  scale_y_log10(breaks = c(.025, .05, .1, .25, .5, 1), labels = c('.025', '.05', '.1', '.25', '.5', '1'),
                name = 'Proportion with p < .05 (log scale)') +
  scale_x_continuous(breaks = c(0, .3, .7), labels = c('0', '.3', '.7'), name = "y's stability (x's is .7)") +
  scale_color_manual(values = c(none = '#7a7a7a', uncorrelated = '#2a78d6', 'correlated (.3)' = '#d4572a'), name = 'Traits') +
  scale_shape_manual(values = c(none = 16, uncorrelated = 17, 'correlated (.3)' = 15), name = 'Traits') +
  scale_alpha_manual(values = c('50' = .5, '150' = 1), name = 'N') +
  scale_linetype_manual(values = c('6' = 'dashed', '12' = 'solid'), name = 'Waves') +
  labs(caption = 'No within-wave correlation and no cross-lag. With traits, each point averages trait variance .5 and 1.') +
  sim_theme() + theme(axis.text.x = element_text(angle = 0, hjust = .5))

#' ## With a true effect
#' 
#' With a lagged effect of .1, power is only meaningful for models that hold their error rate. The
#' RI-CLPM comes closest: its null rates are `r rng(fac_cells('riclpm', 'lag', is_('clxy', 0))$reject_rate)`,
#' slightly inflated at N = 50 with 6 waves, which is also where most of its fit warnings are (and it's the
#' generating model here, so this is its best case). It recovers the effect (average estimate
#' `r f3(fac_mean('mean_est', 'riclpm', 'lag', is_('clxy', .1)))`), and the number of waves matters more
#' for finding it than N does: `r f3(rp(50, 6))` at N = 50 with 6 waves, `r f3(rp(50, 12))` with 12,
#' and `r f3(rp(150, 6))` and `r f3(rp(150, 12))` at N = 150. The other models' rates under the alternative mix the
#' effect with their biases, so I haven't summarized them here; they're in `summary_cells.csv`.
#' 

#+echo = T
ric_null <- fac_cells('riclpm', 'lag', is_('clxy', 0))
hot_cell <- function(d) d$n == 50 & d$waves == 6
check_claim(abs(fac_mean('mean_est', 'riclpm', 'lag', is_('clxy', .1)) - .1) < .005 &&
              max(ric_null$reject_rate) < .08 &&
              mean(ric_null$reject_rate[hot_cell(ric_null)]) == max(tapply(ric_null$reject_rate, paste(ric_null$n, ric_null$waves), mean)) &&
              mean(ric_null$warn_rate[hot_cell(ric_null)]) == max(tapply(ric_null$warn_rate, paste(ric_null$n, ric_null$waves), mean)) &&
              rp(50, 12) - rp(50, 6) > rp(150, 6) - rp(50, 6),
            'the RI-CLPM recovers .1; its null rates are below .08, most inflated at N = 50 / 6 waves where warnings are most common; waves matter more than N for power')

#' # Limits
#' 
#' The factorial only runs N = 50 and 150, with x's stability fixed at .7 and y's never above it, one
#' nonzero level each of the within-wave correlation, the trait correlation and the true effect, no
#' effect of y on x, equal trait variances for x and y, and complete, balanced, normal data. The RI-CLPM
#' is also the model that generates the data with traits, and it isn't fit to the data without them. N = 600 and 24 waves (most of the computing
#' cost) are the obvious next additions; they'd say whether the centering bias keeps
#' shrinking with more waves and whether the biases that grow with N show up the same way here.
#' 
#' 
#' # References
#' 
