library(rdhs)
library(survival)
library(dplyr)
library(haven)
library(tidyr)

#' Malawi 2015 DHS

ir <- readRDS(get_datasets("MWIR7AFL.ZIP")[[1]])

df <- ir %>%
  transmute(
    caseid, 
    psu = v021,
    strata = v023,
    weights = v005,
    dob_cmc = v011,
    doi_cmc = v008,
    afs = v525,
    afs_recode = v531,
    afs_flag = as_factor(v532),
    first_birth_cmc = v211,
    marital_status = as_factor(v501),
    marriage_age = v511,
    first_union_cmc = v509,
    n_union = v503
  )

#' Minimum and maximum first sex month based on AFS
#'
#' Minumum defined by:
#' 
#' * start of AFS year (dob_cmc + afs * 12)
#'
#' Maximum defined by minimum of:
#'
#' * end of AFS year (dob_cmc + afs*12 + 12)  [Need to think about if this should be + 12 or +11]
#' * month of interview
#' * month of first union
#' * month of first birth minus 9
#'
#' There are going to be a LOT of inconsistent cases (max < min).  For now just take
#' these to be the maximum; but ponder more...


df <- df %>%
  mutate(
    eversex = afs > 0,
    afs = if_else(afs == 0, NA_real_, as.double(afs)),
    first_sex_cmc_min = if_else(eversex, dob_cmc + afs * 12, NA_real_),
    first_sex_cmc_max = if_else(eversex,
                               pmin(dob_cmc + afs * 12 + 12,
                                    first_union_cmc,
                                    first_birth_cmc - 9,
                                    doi_cmc,
                                    na.rm = TRUE),
                               NA_real_),
    first_sex_cmc_min_adj = pmin(first_sex_cmc_min, first_sex_cmc_max)
  )

#' Inconsistent for about 23% of cases! (mostly based on union before first sex)
df %>%
  count(first_sex_cmc_min > first_sex_cmc_max)

#' Impute month of first sex

df <- df %>%
  mutate(
    first_sex_cmc = round(runif(n(), first_sex_cmc_min_adj, first_sex_cmc_max))
  )


#' Minimal dataset; no faffing with weights, etc.
marriage_epis <- df %>%
  select(caseid,
         dob_cmc,
         doi_cmc,
         eversex,
         first_sex_cmc,
         first_union_cmc)

#' Define start and end of calendar periods, duration periods and age periods
marriage_epis <- marriage_epis %>%
  filter(eversex) %>%
  mutate(
    time_start = first_sex_cmc,
    time_end = pmin(first_union_cmc, doi_cmc, na.rm = TRUE),
    event = as.integer(!is.na(first_union_cmc) & first_union_cmc == time_end),
    duration_start = 0,
    duration_end = time_end - time_start,
    age_start = time_start - dob_cmc,
    age_end = time_end - dob_cmc
  )

#' survSplit() can't handle cases where tstart == tstop (first sex coincides with first union)
#' Drop these.

marriage_epis <- filter(marriage_epis, time_end > time_start)

#' Split episodes by year up to 5+ years

marriage_epis_durspl <- survSplit(Surv(duration_start, duration_end, event) ~ .,
                                  data = marriage_epis, cut = c(12, 24, 36, 48, 60),
                                  start = "duration_start", end = "duration_end",
                                  episode = "duration_cat", event = "event")

marriage_epis_durspl$duration_cat <- factor(marriage_epis_durspl$duration_cat, 1:6,
                                            c("<1 year", "1 year", "2 years", "3 years",
                                              "4 years", "5+ years"))

#' Need to update start/stop of time and age episodes

marriage_epis_durspl <- marriage_epis_durspl %>%
  mutate(
    time_start = first_sex_cmc + duration_start,
    time_end = first_sex_cmc + duration_end,
    age_start = time_start - dob_cmc,
    age_end = time_end - dob_cmc
  )

#' Verify that same person-months for duration, time, and age episodes

marriage_epis_durspl %>%
  summarise(
    time_pms = sum(time_end - time_start),
    duration_pms = sum(duration_end - duration_start),
    age_pms = sum(age_end - age_start),
    events = sum(event)
  )

#' Split episodes by five-year age groups from <10 to >35


marriage_epis_agespl <- survSplit(Surv(age_start, age_end, event) ~ .,
                                  data = marriage_epis_durspl, cut = 12 * seq(10, 35, 5),
                                  start = "age_start", end = "age_end",
                                  episode = "age_cat", event = "event")

marriage_epis_agespl$age_cat <- factor(marriage_epis_agespl$age_cat, 1:7,
                                       c("<10 years", "10-14", "15-19", "20-24", "25-29",
                                         "30-34", ">=35 years"))

#' Update start/stop of duration and time episodes

marriage_epis_agespl <- marriage_epis_agespl %>%
  mutate(
    time_start = dob_cmc + age_start,
    time_end = dob_cmc + age_end,
    duration_start = time_start - first_sex_cmc,
    duration_end = time_end - first_sex_cmc
  )

#' Verify that same person-months for duration, time, and age episodes

marriage_epis_agespl %>%
  summarise(
    time_pms = sum(time_end - time_start),
    duration_pms = sum(duration_end - duration_start),
    age_pms = sum(age_end - age_start),
    events = sum(event)
  )


marriage_epis_agespl$age_cat <- relevel(marriage_epis_agespl$age_cat, "20-24")

fitA <- glm(event ~ duration_cat + age_cat, family = "poisson",
            data = marriage_epis_agespl, offset = log((time_end - time_start) / 12))

summary(fitA)



#' ## Use pyears() function
#'
#' For Poisson regression, we only need to aggregate the counts of person-time
#' and events. We can use the pyears() function with tcut() as a shortcut to
#' do cuts on all three dimensions at once.
#'
#' For some reason, tcut() doesn't like open-ended intervals (ending in Inf)...

marriage_epis <- marriage_epis %>%
  mutate(
    duration_cat = tcut(time_start - first_sex_cmc, 12 * c(0:5, 100),
                        c("<1 year", "1 year", "2 years", "3 years", "4 years", "5+ years")),
    age_cat = tcut(time_start - dob_cmc, 12 * c(0, seq(10, 35, 5), 100),
                   c("<10 years", "10-14", "15-19", "20-24", "25-29", "30-34", ">=35 years"))
  )
  

#' Argument `scale=` converts to years. To keep person-time in months, specify `scale = 1`.
#' To convert to years, specify `scale = 12`.

marriage_pyears <- pyears(Surv(time_end - time_start, event) ~ duration_cat + age_cat,
                          data = marriage_epis, scale = 12, data.frame = TRUE)
                          

marriage_pyears$data

marriage_pyears$data$age_cat <- relevel(marriage_pyears$data$age_cat, "20-24")

fitB <- glm(event ~ duration_cat + age_cat, family = "poisson",
            data = marriage_pyears$data, offset = log(pyears))

summary(fitB)

summary(fitA)