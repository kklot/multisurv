#' ---
#' title: "Married, divorce, widowed in women - Malawi 2015"
#' bibliography: "../zotero_biblatex.bib"
#' output:
#'   pdf_document:
#'     toc: true
#'     number_sections: true
#'     keep_tex: false
#'     include:
#'       in_header: "~/templates/floatHforRmarkdown.tex"
#' documentclass: article
#' linkcolor: blue
#' urlcolor: blue
#' citecolor: red
#' ---
#'

#+ packages and config, include=FALSE
library(rdhs)
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(janitor)
tabyl <- function(dat, ...) janitor::tabyl(dat, ..., show_missing_level = FALSE)
devtools::load_all("~/Code/R/ktools/")
library(knitr)
opts_chunk$set(echo = FALSE, cache = FALSE, out.extra = "")

# Pull with `rdhs`
# set_rdhs_config(email = "ath19@ic.ac.uk", project = "Statistics and Machine Learning for HIV")
# datasets <- dhs_datasets("MW", fileFormat = "flat", fileType = c("IR")) %>%
#     filter(SurveyYear == 2015)  %>%
#     mutate(sex = if_else(FileType == "Individual Recode", "female", "male"))
# datasets
# downloads <- get_datasets(datasets$FileName)
# names(downloads)

x <- 1 # adapt from multiple datasets
o <- readRDS(downloads[[x]]) %>%
    as_tibble() %>%
    mutate(lab = datasets$SurveyId[x])

# Extract and converting date time 
dta <- o %>%
    select(
        psu = v021, strata = v023, weights = v005, dob = v011, doi = v008, afs = v531, lab,
        marriage_age = v511,
        n_union = v503,
        month_1st_union = v507,
        year_1st_union = v508,
        cmc_uninon = v509,
        marital_status = v501
    ) %>%
    mutate(
        marriage_age =  if_else(is.na(marriage_age), 0L, marriage_age),
        marriage_age2 = (cmc_uninon - dob) / 12,
        iso = substr(lab, 1, 2),
        yob = cmc_to_year(dob), svy = cmc_to_year(doi),
        age = svy - yob, sex = datasets$sex[x]
    )

#' cleaning
#' # Data cleaning
#'
#' Remove those
#'
#' - afs or marriage age greater than age
#' - married but no afs

dta %<>% filter(!(afs > age | marriage_age > age))
dta %<>% filter(!(marital_status != 0 & afs == 0))
# dta %<>% filter(afs <= marriage_age)
dta <- dta |>
    mutate(across(c(marital_status, n_union), as_factor)) %>%
    mutate(marital_status = fct_recode(
        marital_status,
        separated = "no longer living together/separated",
        union = "living with partner"
    ))

#' # Descriptive statistics
#'
#' ## Key variables
#'
#' All married individuals responded to number of union question. Among those
#' with more than once unions, more than 19% had more than one union in any of
#' the marital statuses. The following differences are needed to be considered
#' when modelling:
#'
#' - The number of currently married with more than one unions reflects
#'   *partially* the number of remarried. Because remarried can also be
#'   currently in other states, including living with partner, widowed,
#'   divorced, and separated. Using this as an estimate of remarried rate would
#'   underestimate it.
#' - Among more than once union respondents, how was the 1st union ended is
#'   unknown in any of the states.

#+ echo=F, results="asis"
dta %>% tabyl(marital_status) |> adorn_pct_formatting() |> kable(caption = "Marital distribution")

dta %>% tabyl(n_union) |> adorn_pct_formatting() |> kable(caption = "Number of union")

dta %>%
    tabyl(marital_status, n_union) %>%
    adorn_percentages() %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    kable(caption = "Percent rowwise")

dta %>%
    tabyl(marital_status, n_union) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    kable(caption = "Percent colwise")

#' ## Union - time since sexual debut to first union
#'
#' ### Delay sexual intercourse in child marriage
#' 
dta %<>%
    mutate(
        time_since_debut_c = marriage_age2 - afs,
        marriage_age_d = findInterval2(marriage_age2, seq(15, 30, 5)),
        afs_d = findInterval2(afs, seq(15, 30, 5)),
        time_since_debut_d = findInterval2(time_since_debut_c, c(-2, -1, 0, 1, 2, 5, 10)), 
        time_since_debut_d = fct_recode(time_since_debut_d, '<-3'='-12--3')
    )
#' Very young age of marriage has larger proportion of having *sexual debut
#' after marriage* (\@ref(fig:time_since_db_fig)), up to 50% in those married
#' under 14 (\@ref(tab:time_since_db_tb)). This period is mostly in a year or
#' two; three or more years delay of sexual intercourse after marriage is
#' visible only in those married at age under 14. We might assume the delay is
#' correctly reported.
#'
#' This implies risk of STDs transmission in very young age should not be
#' inferred based on married status but AFS. Also, this group is already in a
#' stable relationship. These suggest excluding this group from the analyses of
#' union formation is reasonable and standard survival model can be used.
#'
#+ time_since_db_fig, fig.cap="Time since debut to married by age of marriage"
dta %>%
    group_by(time_since_debut_d, marriage_age_d) %>%
    count() %>%
    filter(marriage_age_d != "NA-NA") %>%
    ggplot() +
    geom_col(aes(marriage_age_d, n, fill = time_since_debut_d), position = position_fill()) +
    theme(axis.text.x = element_text(angle = 30)) +
    scale_fill_manual(values = ktools::gen_colors(okabe, 8)) +
    labs(x = "Age at marriage", title = "Time since AFS to first union")

#+ fig.cap="Time since debut to married by age of marriage - filtered"
dta %>%
    filter(time_since_debut_c >= 0) %>%
    group_by(time_since_debut_d, marriage_age_d) %>%
    count() %>%
    filter(marriage_age_d != "NA-NA") %>%
    ggplot() +
    geom_col(aes(marriage_age_d, n, fill = time_since_debut_d), position = position_fill()) +
    theme(axis.text.x = element_text(angle = 30)) +
    scale_fill_manual(values = ktools::gen_colors(okabe, 8)) +
    labs(x = "Age at marriage", title="Time since AFS to first union - exclude delay sex group")

#' ### Late sexual debut - short time to marriage
#' 
#' Plot time since debut by AFS shows late sexual debut shorten the time to
#' marriage. Might need to include this in the model as covariate.
#'
#+ fig.cap="Time since debut to married by age at first sex"
dta %>%
    filter(time_since_debut_c >= 0) %>%
    filter(afs != 0) %>%
    group_by(time_since_debut_d, afs_d) %>%
    count() %>%
    filter(time_since_debut_d != "NA-NA") %>%
    ggplot() +
    geom_col(aes(afs_d, n, fill = time_since_debut_d), position = position_fill()) +
    theme(axis.text.x = element_text(angle = 30)) +
    scale_fill_manual(values = ktools::gen_colors(okabe, 9)) +
    labs(x = "Age at first sex", title = "Time since AFS to first union")

dta %<>%
    rownames_to_column("id") %>%
    mutate(afs = as.numeric(haven::zap_labels(afs)), marriage_age = as.numeric(marriage_age)) %>%
    mutate(
        t_virgin = if_else(afs == 0, age, afs),
        t_debut = if_else(afs != 0, afs, NA_real_),
        t_union = if_else(marriage_age != 0, marriage_age, NA_real_),
        t_union = if_else(marital_status == 2, age, t_union), 
        t_widowed = if_else(marital_status == 3, age, NA_real_), 
        t_separate = if_else(marital_status %in% 4:5, age, NA_real_),
        censored = marital_status %in% 2:5
    )
dta

library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE)
#' plot raw data
library(ggpubr, help, pos = 2, lib.loc = NULL)

four_state_raw <- ggarrange(
    dta %>% ggplot() +
        geom_density(aes(t_debut)) +
        labs(title = "t_debut", x = "Time"),
    dta %>%
        mutate(diff = if_else(censored, t_union, t_union - t_debut)) %>%
        filter(diff >= 0) %>% 
        ggplot() +
        geom_density(aes(diff, fill = factor(censored)), alpha=.7) +
        labs(title = "Debut before Married", fill = "Censored", x = "Time"),
    dta %>%
        mutate(diff = if_else(censored, t_union, t_union - t_debut)) %>%
        filter(diff < 0) %>% 
        ggplot() +
        geom_density(aes(diff, fill = factor(censored)), alpha=.7) +
        labs(title = "Debut after Married", fill = "Censored", x = "Time"),
    
    dta %>%
        mutate(
            diff = if_else(marriage_age != 0, t_separate - t_union, t_separate),
            know_aam = marriage_age != 0
        ) %>%
        ggplot() +
            geom_density(aes(diff, fill = factor(know_aam)), alpha = .7) +
            labs(title = "Separated", fill = "Known age at married?", x = "Time"),
    
    dta %>%
        mutate(
            diff = if_else(marriage_age != 0, t_widowed - t_union, t_widowed),
            know_aam = marriage_age != 0
        ) %>%
        ggplot() +
            geom_density(aes(diff, fill = factor(know_aam)), alpha = .7) +
            labs(title = "Widowed", fill = "Known age at married?", x = "Time"),
   
    ncol = 3, nrow = 2
)
ggsave("img/four_state_raw.png")

#' Prepare msm format
msdta <-
dta %>% 
    select(age, marital_status, marriage_age, starts_with("t_")) %>%
    filter(t_union >= t_debut) %>%
    mutate(
        pid = 1:n(),
        t_debut = if_else(t_debut == t_virgin, t_virgin + 1 / 12, t_debut),
        t_union = if_else(t_union <= t_debut, t_debut + 1 / 12, t_union),
        t_separate = if_else(t_separate <= t_union, t_union + 1 / 12, t_separate),
        t_widowed = if_else(t_widowed <= t_union, t_union + 1 / 12, t_widowed),
    ) %>%
    pivot_longer(starts_with("t_"), names_to = "state", values_to = "time", names_prefix = "t_") %>%
    mutate(
        state = factor(state, levels = char(virgin, debut, union, separate, widowed)),
        state_numeric = as.numeric(state)
    ) %>% drop_na()
    
msdta

saveRDS(msdta, "data/mw2015.rds")

#' 
#' ## Dissolution - time since married
#' 
dta %<>%
    mutate(
        time_since_married_c = (doi - cmc_uninon) / 12,
        time_since_married_d = findInterval2(time_since_married_c, seq(0, 20, 5))
    )
#' ### Once union - competing risk model with right censored data
#'
#' **Among those with only once union**, a majority of the population have
#' remained in the same marital status for more than 5 years. ~~Nearly 95% of
#' the widowed respondents has not remarried for 5 to 20+ year after the first
#' union~~. But since we do not know the time of the event, e.g., time of
#' partner's death, the above interpretation is only correct for the married and
#' partnered group. The others groups, widowed, divorced, and separated, the
#' numbers would more likely reflect the varying in the time of the events.
#'
#' Limited to this set of data, a survival model using the time since married
#' would treat separated, divorce, and widowed as competing events while those
#' in married or partnered states will be right-censored.

#+ fig.cap="Time since married"
dta %>%
    group_by(marital_status, time_since_married_d, n_union) %>%
    count() %>%
    drop_na() %>%
    ggplot() +
    geom_col(aes(marital_status, n, fill = time_since_married_d), position = position_fill()) +
    facet_wrap(~n_union) +
    theme(axis.text.x = element_text(angle = 30))

#' All those in union have records of time of union
#+ union_by_marriage_age, results = "asis"
dta %>%
    filter(n_union == "once") %>%
    tabyl(marriage_age_d, marital_status) %>%
    kable(caption = "Time since union by marital statue")




#' 
#' ## More than once union - left-censoring of both event and time
#' 
#' **Among those with more than once unions**, it is unknown what states a
#' respondent has passed through before reaching the current state. Inexact
#' number of unions also does not allow to narrow down the possible pathways. If
#' we are interested in hazard of union dissolution from the first union, where
#' union dissolution is defined as either widowed, divorce, or separate, we can
#' treat all the states as left-censored at the current time since married.
#'
#' However, if we are interested in in specific events of widowed or divorce,
#' how to specify which state a respondent is censored? Could a probability of
#' the first union ended in widowed or divorce be used as a weight in specifying
#' one of the censored events (likelihood contribution). Event if this works,
#' there are still cases where both states widowed and divorce has been passed
#' through. So we could only limit the analyses at best to hazard of an event
#' *after the first union*.
#' 

saveRDS(msdta, "data/mw2015.rds")
