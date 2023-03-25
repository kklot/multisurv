source("~/Documents/libs.r")
# devtools::install('~/git/adcomp/TMB')

MAX_AGE = 50
N_PAR = 8

library(TMB)
TMB::compile("code/model.cpp")
base::dyn.load(TMB::dynlib("code/model"))
invisible(TMB::config(tape.parallel = 0, DLL = "model"))
TMB::openmp(1)

data <- readRDS(here("data/prep.rds")) %>%
    mutate(across(everything(), haven::zap_label)) %>%
    mutate(across(everything(), as.integer)) %>%
    slice_sample(n = 100) %>%
    as.list()

expect_true(max(data$age) <= MAX_AGE)
expect_true(all((data$afs + data$n_N) <= data$age))
expect_true(all((data$aam + data$n_N) <= data$age))

data$sd_b <- c(0, 0.1)
data$sd_q <- c(log(0.1), 1) # log normal mean and sd
str(data)

init <- list(
    betav = rnorm((MAX_AGE + 1) * N_PAR, 0, 0.1),
    lqv = rnorm(N_PAR, log(0.1), 1),
    pacf = c(0, 1)
)
str(init)

obj <- TMB::MakeADFun(data, init, DLL = "model", silent = FALSE)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

