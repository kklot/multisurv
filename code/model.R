source("~/Documents/libs.r")
# devtools::install('~/git/adcomp/TMB')

library(TMB)
TMB::compile("code/model.cpp")
base::dyn.load(TMB::dynlib("code/model"))
invisible(TMB::config(tape.parallel = 0, DLL = "model"))

init <- list(
    betav = rnorm(2, 0, 0.1),
    lqv = rnorm(8, log(0.1), 1)
)

data <- readRDS(here("data/prep.rds")) %>%
    mutate(across(everything(), haven::zap_label)) %>%
    mutate(across(everything(), as.integer)) %>%
    slice_sample(n = 100) %>%
    as.list()

data$sd_b <- c(0, 0.1)
data$sd_q <- c(log(0.1), 1) # log normal mean and sd
str(data)

TMB::openmp(1)
obj <- TMB::MakeADFun(data, init, DLL = "model", silent = FALSE)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

