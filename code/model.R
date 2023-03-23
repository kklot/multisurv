library(TMB)
TMB::compile("code/model.cpp")
base::dyn.load(TMB::dynlib("code/model"))
TMB::config(tape.parallel = 0, DLL = "model")

init <- list(
    betav = rnorm(2, 0, 0.1),
    lqv = rnorm(11, 0, 0.1)
)

data <- readRDS(here("data/prep.rds")) %>%
    mutate(across(everything(), as.integer)) %>%
    as.list()

data$sd_b <- data$sd_q <- c(0, .3)
str(data)

TMB::openmp(2)
obj <- TMB::MakeADFun(data, init, DLL = "model")
fit <- nlminb(obj$par, obj$fn, obj$gr)
