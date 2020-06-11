## Model_runs
# Follows scripts "data_preparation.R" and "stan_models.R"
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

# Auxiliary function to extract parameters from model fit
ext_pars <- function(mod_fit, startval = FALSE){
  pars <- paste0(c("alphaB", "betaB"), startval + 1)
  sel_pars <- extract(mod_fit, pars = pars)
  calc_p <- function(vec) 
    round(c(Mean = mean(vec), quantile(vec, c(0.025, 0.975))), 4)
  rbind(
    Diff = calc_p(sel_pars[[pars[2]]]),
    Herbs = calc_p(sel_pars[[pars[1]]]),
    Woody = calc_p(sel_pars[[pars[1]]] + sel_pars[[pars[2]]])
  )
}

## Winter freezing
# Mortality
dat <- prep_dat(dat_wf, "alive", tree, spec_codes, startval = FALSE, 
  t_log = FALSE)
mod_fit <- stan(model_code = mod_code_mort, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/wf_alive")

# Plain view
dat <- prep_dat(dat_wf, "plainw", tree, spec_codes, startval = FALSE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/wf_plainw")

# Height
dat <- prep_dat(dat_wf, "length", tree, spec_codes, startval = FALSE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/wf_length")

# Number of leaves
dat <- prep_dat(dat_wf, "nleaves", tree, spec_codes, startval = FALSE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/wf_nleaves")

# Total biomass
dat <- prep_dat(dat_wf, "bio", tree, spec_codes, startval = FALSE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/wf_bio")

# New biomass
dat <- prep_dat(dat_wf, "bion", tree, spec_codes, startval = FALSE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/wf_bion")

## Spring freezing
# Mortality
dat <- prep_dat(dat_sf, "alive", tree, spec_codes, startval = FALSE, 
  t_log = FALSE)
mod_fit <- stan(model_code = mod_code_mort, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/sf_alive")

#Plain view
dat <- prep_dat(dat_sf, "plainw", tree, spec_codes, startval = TRUE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code_startval, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit, startval = TRUE)
# save(mod_fit, file = "models/sf_plainw")

# Height
dat <- prep_dat(dat_sf, "length", tree, spec_codes, startval = TRUE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code_startval, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit, startval = TRUE)
# save(mod_fit, file = "models/sf_length")

# Number of leaves
dat <- prep_dat(dat_sf, "nleaves", tree, spec_codes, startval = TRUE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code_startval, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit, startval = TRUE)
# save(mod_fit, file = "models/sf_nleaves")

## Herbivory
# Mortality
dat <- prep_dat(dat_he, "alive", tree, spec_codes, startval = FALSE, 
  t_log = FALSE)
mod_fit <- stan(model_code = mod_code_mort, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit)
# save(mod_fit, file = "models/he_alive")

#Plain view
dat <- prep_dat(dat_he, "plainw", tree, spec_codes, startval = TRUE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code_startval, data = dat, iter = 4000, 
  control=list(adapt_delta=0.999, max_treedepth=12), seed = 105)
ext_pars(mod_fit, startval = TRUE)
# save(mod_fit, file = "models/he_plainw")

#Height
dat <- prep_dat(dat_he, "length", tree, spec_codes, startval = TRUE, 
  t_log = TRUE)
mod_fit <- stan(model_code = mod_code_startval, data = dat, iter = 4000, 
  control = list(adapt_delta = 0.98), seed = 105)
ext_pars(mod_fit, startval = TRUE)
# save(mod_fit, file = "models/he_length")
