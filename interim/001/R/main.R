library(rstan)
library(data.table)
library(OrvacRCT)
library(brms)

# documentation for stan functions via lookup
rstan::lookup("inv_logit")

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())


"
# Notes


## Reminders

renaming columns in data.table :
  rename serot3 setnames(d, \"serot3\", \"y\")

brms style priors:
stan_dat <- brms::make_standata(formula = as.formula(y~trt), d)
  stan_priors <- c(
  brms::prior(student_t(3, 0, 2.5), class = \"Intercept\"),
  brms::prior(student_t(3, 0, 2.5), class = \"b\")
)

brms make stancode:
mod_defn <- brms::make_stancode(formula = as.formula(y~trt), d,
                                  family = \"bernoulli\",
                                  prior = stan_priors,
                                  save_model = myf)

rstan sampling reminder:
fit <- rstan::sampling(object  = mod,
                         data    = ld,
                         #pars    = pars,
                         #init    = gen_init,
                         chains  = 1,
                         iter    = 1000 #,
                         #warmup  = nwarmup,
                         #thin    = nthin,
                         # control = list(adapt_delta   = adapt_delta,
                         #                stepsize      = stepsize,
                         #                max_treedepth = max_treedepth)
                         )
"
main <- function(){
  
  library(rstan)
  library(data.table)
  library(OrvacRCT)
  library(brms)
  
  # documentation for stan functions via lookup
  rstan::lookup("inv_logit")
  
  # ensure that the compiled version of the stan model
  # is written to the hard disk in the same dir
  # as the .stan file.
  rstan_options(auto_write = TRUE)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
  options(mc.cores = parallel::detectCores())
  
  # temp simulated data
  mydir <- file.path(getwd(), "interim", "001", "etc")
  cfg <- OrvacRCT::get_cfg(cfgdir = mydir, print = F)
  d <- data.table(OrvacRCT::get_data(cfg))
  names(d) <- cfg$dnames
  
  # only have 70 observations
  N <- cfg$n_start
  d <- d[1:N,]

  ld <- list(
    N = cfg$n_start,
    N_pred = cfg$n_max_sero - cfg$n_start,
    K = 1,
    X = array(d$trt, dim = c(cfg$n_start, 1)),
    # pass X_new in but only used if predict is T
    X_new = array(rep(0:1, len = cfg$n_max_sero - cfg$n_start), 
                  dim = c(cfg$n_max_sero-nrow(d), 1)),
    y = d$serot3
  )

  # compile code
  myf <- file.path(getwd(), "interim", "001", "stan", "logistic001.stan")
  mod <- rstan::stan_model(myf ,verbose = F)
  
  fit <- rstan::sampling(object  = mod,
                         data    = ld,
                         chains  = 1,
                         iter    = 10000 
                         )
  
  
  # report odds ratio and probabilities - 
  # compute prob and diff in generated data block
  
  
  # assess futility
  res <- rstan::extract(fit)
  y_rep <- res$y_rep
  trtnew <- c(d$trt, rep(0:1, len = cfg$n_max_sero - cfg$n_start))
  table(trtnew)
  dim(y_rep)
  
  # set up new data
  ld <- list(
    N = cfg$n_max_sero,
    N_pred = 1,
    K = 1,
    X = array(trtnew, dim = c(cfg$n_max_sero, 1)),
    # pass X_new in but only used if predict is T
    X_new = array(0, dim = c(1, 1)),
    y = NA
  )
  
  nsim <- 1000
  win <- numeric(nsim)
  for(i in 1:nsim){
    
    ld$y <- c(d$serot3, y_rep[i, ])
    fit_max <- rstan::sampling(object  = mod,
                           data    = ld,
                           pars    = c("alpha", "b"),
                           chains  = 1,
                           iter    = 1000 
    )
    
    m <- as.matrix(fit_max)
    # hist(m[, "b[1]"])
    win[i] <- mean(m[, "b[1]"]>0) > cfg$thresh_p_sup
  }
  
  futile <- mean(win) < cfg$thresh_p_fut
  
  message(paste0("Of the ", nsim, " simulations, there were ", sum(win),
                 " trials where the probability of a",
                 " treatment effect was > ", cfg$thresh_p_sup, "."))
  
  msg <- ifelse(futile, "futile and should be abandoned.", 
                " not futile and enrollments should continue.")
  
  message(paste0("Using a futility threshold of  ", cfg$thresh_p_fut,
                 " this implies that the trial is ",
                 msg))
  
  if(futile) stop("Trial abandoned on basis of futility. Stop.")
  
  logoddtrt <- as.numeric(res$b)
  ptrt <- mean(logoddtrt>0)
  stopvsamp <- ptrt > cfg$thresh_pp_es
  
  message(paste0("The posterior probability that the log-odds of",
                 " the treatment effect was greater than zero is ",
                 round(ptrt, 2), "."))
  
  msg <- ifelse(stopvsamp, "cease.", "continue.")
  
  message(paste0("Using a decision threshold of  ", cfg$thresh_pp_es,
                 " this implies that venous sampling should ",
                 msg))
  
}


main()
