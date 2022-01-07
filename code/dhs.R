#' Pull with `rdhs`
library(rdhs)
library(tidyverse)
library(ktools)
#' get one sex one country first
datasets <- dhs_datasets("MW", fileFormat = "flat", fileType = c("IR")) %>%
    filter(SurveyYear == 2015)  %>%
    mutate(sex = if_else(FileType == "Individual Recode", "female", "male"))
datasets
downloads <- get_datasets(datasets$FileName)
names(downloads)

#' Extract 
x <- 1 # adapt from multiple datasets
o <- readRDS(downloads[[x]]) %>% as_tibble %>% mutate(lab = datasets$SurveyId[x])

dta <- o %>% 
select(
    psu = v021, strata = v023, weights = v005, dob = v011, doi = v008, afs = v531, lab, 
    marriage_age = v511, 
    marital_status = v501
) %>%
mutate(
    iso = substr(lab, 1, 2),
    yob = cmc_to_year(dob), svy = cmc_to_year(doi), 
    age = svy - yob, sex = datasets$sex[x]
)

#' cleaning
dta %<>% filter(afs <= age, marriage_age <= age)

#' states indicators
dta %>% pull(marital_status) %>% table(useNA='a')

dta %<>%
mutate(
    eversexi = afs != 0, 
    marriagei = marital_status %in% 1, 
    # partneri = marital_status %in% 1:2, 
    widowedi = marital_status %in% 3, 
    divorcei = marital_status %in% 4 # how to treat separated (2)
) %>%
select(yob, age, afs, marriage_age, eversexi, marriagei, widowedi, divorcei) %>%
mutate(across(where(is.logical), as.integer))

#' multi stage format
dta %<>%
mutate(
    age = as.integer(age),
    eversex = if_else(eversexi==1, afs, age),
    marriage = if_else(marriagei==1, marriage_age, age),
    widowed = if_else(widowedi==1, age, age), 
    divorce = if_else(widowedi==1, age, age)
)
library("mstate")
tmat <- mstate::transMat(
    list(c(2), c(3), c(4, 5), c(), c()), 
    c("born", "eversex", "marriage", "widowed", "divorce"))
tmat

#' number of transition
dta$oid <- 1:nrow(dta)
twotrans_id <- dta %>% filter(eversexi==1, marriagei==1) %>% select(oid)

dta %>% filter(eversexi==1, marriagei==1)
dta %>% filter(eversexi==1, divorcei==1)
    dta %>% filter(eversexi==1, marriagei==1, widowedi==1)

msebmt <- msprep(
    data = dta, trans = tmat, 
    time = c(NA, "eversex", "marriage","widowed", "divorce"), 
    status = c(NA, "eversexi", "marriagei", "widowedi", "divorcei")
    )

library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE)
fits_wei <- vector(mode = "list", length = n_trans)


for (i in 1:n_trans){
  fits_wei[[i]] <- flexsurvreg(
      Surv(time, status) ~ 1,
      data = subset(msebmt, trans == i), 
      dist = "lnorm")
}

pat_2 <- data.frame(msebmt[msebmt$id == 1, ][1, ])
yr_grid <- seq(0, 30, .1)
cumhaz_grid <- seq(0, max(msebmt$time), .01)
cumhaz_pat_2 <- msfit.flexsurvreg(
    fits_wei, trans = tmat, 
    t = cumhaz_grid,
    newdata = pat_2,
    variance = FALSE)

cumhaz_pat_2$Haz %>%
ggplot() +
geom_line(aes(time, Haz, color = factor(trans))) +
labs()
