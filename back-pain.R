

# ---- libs ----

suppressPackageStartupMessages(suppressWarnings({
  
  library("dplyr") # tidyverse
  library("tidyr")
  library("readr")
  library("forcats")
  library("ggplot2")
  
  library("GGally") # additional ggplot-type plotting
  
  library("compositions")
  library("zCompositions") # this one for lr_EM
  
  library("performance") # model checking
  library("mice")        # missing data functions
  library("car")         # Anova() for comparing models
  library("knitr")       # kable() for pretty printing
  library("foreach")     # powerful looping
  
  library("boot")        # bootstrap confidence intervals
  library("tictoc")      # check time between tic() and toc()
  
}))
  
# ---- notes ----

# Introduction

# I want to perform compositional isotemporal substitution based on
# binary regression. I would like to ask for a review of this R script.
# I want to observe results expressed in odds ratios (OR).

# Data

# Exposture:
#   - Time_Sleep: in minutes per day
# - Time_Sedentary: in minutes per day
# - Time_LPA: in minutes per day
# - Time_MVPA: in minutes per day

# There were 4 low back pain questions:
#   
#   Question 1: On how many days have you experienced low back pain in the last 12 months?
#   Response options: 0 days / 1-7 days / 8-30 days / 31-90 days / more than 90 days, but not every day / every day
# 
# Question 2: How would you rate the average intensity of your low back pain during the last 12 months (average pain intensity on days when you experienced pain)?
#   Response options: 0 - 100 (ranging from no pain to worst pain imageable)
# 
# Question 3: How intense was the worst low back pain that you experienced in the last month?
#   Response options: 0 - 100 (ranging from no pain to worst pain imageable)
# 
# Question 4: How intense was the worst low back pain that you experienced in the last week?
#   Response options: 0 - 100 (ranging from no pain to worst pain imageable)
# 
# Note: Only participants who reported to experience pain in the past 12 months (question 1) were asked to complete the questions 2, 3, and 4.
# 
# 
# We want to explore the effect of reallocating behaviours on low back pain:
#   - occurance (by including all participants into the analysis)
# - frequency (by including only low back pain sufferers into the analysis (ie, those who reported to experience pain in the past 12 months))
# - intensity (by including only low back pain sufferers into the analysis)
# 
# 
# In case models will not fit, we can compute the outcomes to become a binary outcomes:
#   
#   Pain intensity could be categorised as (Boonstra et al., 2014):
#   - no pain (0)
# - mild pain (1-38)
# - moderate pain (39-57)
# - severe pain (58-100)
# 
# Those cut-offs could be used to compute:
#   - binary pain occurance outcome (no pain / mild+moderate+severe pain; no pain+mild pain / moderate+severe pain; no pain+mild+moderate pain / severe pain) by including all participants into the analysis
# - binary pain intensity outcome (no pain / mild+moderate+severe pain; no pain+mild pain / moderate+severe pain; no pain+mild+moderate pain / severe pain) by including only low back pain sufferers into the analysis
# 
# Binary frequency outcome could be: 1-90 days / more than 90 days (this approximately correspond to: non-chronic pain / chronic pain) by including only low back pain sufferers into the analysis.

# Covariates:
#   - sex: sex of the participant. Categorical variable taking one of two categories: female, male.
# - age: age of the participant. Categorical variable taking one of three categories: younger, middle, older.
# - bmi: bmi of the participant. Categorical variable taking one of two categories: normal, overweight.
# - smoking: smoking status of the participant. Categorical variable taking one of two categories: nonsmoker, smoker.
# - stress: stress status of the participant. Categorical variable taking one of two categories: normal, stressed.
# - education: education level of the participant. Categorical variable taking one of two categories: lower, higher.
# - ses: socioeconomic status of the participant. Categorical variable taking one of three categories: lower, middle, higher.


# ---- consts ----

pred_comps <- c("Time_Sleep", "Time_Sedentary", "Time_LPA", "Time_MVPA")
(D <- length(pred_comps))
pred_covs <- c("age", "sex", "bmi", "stress", "smoking", "education", "ses")
outcs <- 
  c(
    "LBP_frequency_year", "LBP_intensity_year", 
    "LBP_intensity_month", "LBP_intensity_week"
  )

# default RHS of model formulas
# (rhs_formula <- paste(c(paste(pred_covs, collapse = " + "), "ilr"), collapse = " + "))
(rhs_formula <- paste(pred_covs, collapse = " + "))


# this is the sequential binary partition matrix to be used for ilr creation
sbp1 <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)

# a way of creating ilr names automatically from SBP matrix
create_ilr_names <- function(sbp_matrix) {
  ilr_sbp_nms <- apply(sbp_matrix, 1, paste, collapse = "")
  ilr_sbp_nms <- gsub("-1", "-", ilr_sbp_nms)
  ilr_sbp_nms <- gsub("1", "+", ilr_sbp_nms)
  ilr_sbp_nms <- gsub("0", ".", ilr_sbp_nms)
  return(paste0("ilr(", ilr_sbp_nms, ")"))
}
create_ilr_names(sbp1)

do_closure <- function(x, clo_val = 1) {
  return(clo_val * x / sum(x))
}
calc_comp_mean <- function(x, clo_val = 1) {
  unclose_mean <- NULL
  if (is.null(dim(x))) {
    return(x)
  } else if (ncol(x) == 1) { # column matrix
    return(as.numeric(x))
  } else {
    unclose_mean <- apply(x, 2, function(x) exp(mean(log(x))))
  }
  
  return(do_closure(unclose_mean, clo_val = clo_val))
  
}


# ---- read ----


bpd_col_spec <- 
  cols(
    Time_Sleep = col_double(),
    Time_Sedentary = col_double(),
    Time_LPA = col_double(),
    Time_MVPA = col_double(),
    age = col_double(),
    sex = col_character(),
    bmi = col_double(),
    stress = col_character(),
    smoking = col_character(),
    education = col_character(),
    ses = col_character(),
    LBP_sufferer = col_character(),
    LBP_frequency_year = col_character(),
    LBP_intensity_year = col_double(),
    LBP_intensity_month = col_double(),
    LBP_intensity_week = col_double()
  )


bpd <- read_csv("dat/bpd.csv", col_types = bpd_col_spec)
# head(bpd)

summary(bpd)


# ---- tidy ----

# relevel categories in LBP_frequency_year
y_lab <- "LBP_frequency_year"
sort(unique(bpd[[y_lab]]))

bpd[[y_lab]] <-
  if_else(
    bpd[[y_lab]] == "more_than90days_but_not_everyday",
    "91+_not_evday",
    bpd[[y_lab]]
  )

# check
table(bpd[[y_lab]], useNA = "ifany")



bpd[[y_lab]] <- factor(bpd[[y_lab]])
levels(bpd[[y_lab]])
lvls_ord <- c(1, 2, 4, 3, 5, 6)
# right order?
# levels(bpd[[y_lab]])[lvls_ord]
bpd[[y_lab]] <- lvls_reorder(bpd[[y_lab]], lvls_ord)
### right order?
levels(bpd[[y_lab]])


with(bpd, table(LBP_frequency_year, LBP_sufferer, useNA = "ifany"))




### comment these lines for sensitivity analysis
bpd$age <- 
  cut(
    bpd$age, 
    breaks = c(17, 44, 64, 100), 
    labels = c("1_younger", "2_middle", "3_older")
  )
table(bpd$age)

bpd$bmi <- 
  cut(
    bpd$bmi, 
    breaks = c(15, 18.5, 25, 70), 
    right = FALSE,
    labels = c("1_underweight", "2_normal", "3_overweight")
  )
table(bpd$bmi)




# ---- impute ----



# Do I have zero values in my composition? (yes in MVPA)
### See: summary(bpd)

# We need to make compositions before we do the lrEM method. The most straightforward way is to create separate datasets and adjust them accordingly.
comp1 <- bpd[, pred_comps]

# How much participants have zero MVPA? 159 participants (6.8% of the sample)
missingSummary(comp1)
sum(rowSums(is.na(comp1) | (comp1 < 0.1)), na.rm = FALSE)
sum(which_0 <- as.logical(rowSums(is.na(comp1) | (comp1 < 0.1))))

# these are 0 vals anywhere in composition (or NA)
bpd[which_0, pred_comps]

# I have zeroes in MVPA - lrEM function will be applied

# ?lrEM
# what is the smallest time-use value above 0? [in minutes]
min(comp1[comp1 > 0])

thresh_detect <- 10 / 1440 
# thresh_detect <- 0.01

comp1.a <- comp1 / 1440 # Create % based composition
dl <- c(rep(thresh_detect, times = D)) # threshold limit for the replacement            

comp1.zr <- lrEM(comp1.a, label = 0, dl = dl) # conduct the lrEM Zero Replacement
comp1.zr <- as_tibble(comp1.zr * 1440)
# composition is larger than 1440 for those who have imputated MVPA 
# (all behaviours will be proportionally downscaled to fit 1440 min when constructing the composition)

# look at imputed values
comp1.zr[which_0, ]

# build dataset that contain imputed values for our 24-h composition

# add new compositions to other noncompositional data
# remove 24-h data from the datase
# bpd <- subset(bpd, select = -c(id, Time_Sleep, Time_Sedentary, Time_LPA, Time_MVPA)) 
head(bpd[which_0, pred_comps])
bpd <- bpd[, !(colnames(bpd) %in% pred_comps)] # remove ori time-use cols
bpd <- bind_cols(comp1.zr, bpd) # add imputed 24-h data
head(bpd[which_0, pred_comps])




# ---- ilr_create ----


add_ilrs_to_data <- function(dataset, comp_vars = pred_comps, sbp_matrix = sbp1) {

  # the time-use composition
  comp <- dataset[, comp_vars]
  comp <- acomp(comp) # designate it as a compositional variable
  
  # define sequential binary partition (SBP)
  psi1 <- gsi.buildilrBase(t(sbp_matrix)) #  The orthonormal matrix

  # find the mean composition
  (m <- mean(comp)) # comp has been designated as acomp, therefore R knows it’s a composition and returns the compositional mean.
  # cat(
  #   "\nThis is the compositional mean [in mins] of the columns (", 
  #   paste(comp_vars, collapse = ", "), 
  #   ")\n\n", 
  #   sep = ""
  # )
  # print(clo(m, total = 1440)) # to look at the mean in minutes/day.
  # cat("\n\n")
  
  # create isometric log ratios (ilr.1) using the above SBP and orthonormal b asis V=psi1. Put ilrs into a data.frame, (X).
  ilrs_from_comp <- ilr(comp, V = psi1)
  colnames(ilrs_from_comp) <- create_ilr_names(sbp_matrix)
  # colnames(ilrs_from_comp) <- paste0("coord", 1:(length(comp_vars) - 1))
  dataset$ilr <- ilrs_from_comp

  return(dataset)
  
}

# use function: creates the ilr columns nested in the single column "ilr"
bpd <- add_ilrs_to_data(bpd)
# check
bpd[, c("ilr", pred_comps)]

# also create version of data without the nested ilrs
bpd_clean <- as.data.frame(bpd)
bpd_clean$ilr <- NULL # remove nested cols
bpd_clean <- cbind(bpd_clean, as.data.frame(bpd$ilr))


# ---- explore1 ----

### Missing data summary for LBP suffers
bpd_clean %>%
  dplyr::filter(LBP_sufferer == "yes") %>%
  md.pattern(., rotate.names = TRUE)


### Missing data summary for _non_ LBP suffers
bpd_clean %>%
  dplyr::filter(LBP_sufferer == "no") %>%
  md.pattern(., rotate.names = TRUE)

### ===> data doesn't have mistiness for analysis


# ---- explore2 ----


### plot pairwise comparisons of time-use and outcomes
if (FALSE) { # takes 30 sec
  suppressWarnings({
    bpd_clean %>%
      dplyr::select(all_of(pred_comps), all_of(outcs)) %>%
      ggpairs(
        ., 
        progress = FALSE, 
        ggplot2::aes(
          colour = LBP_frequency_year, 
          fill = LBP_frequency_year, 
          alpha = 0.25
        )
      ) + 
      theme_bw() +
      scale_colour_viridis_d()+
      scale_fill_viridis_d()
  })
}

### plot pairwise comparisons of _ilrs_ and outcomes
if (TRUE) { # takes 30 sec
  suppressWarnings({
    bpd_clean %>%
      dplyr::select(starts_with("ilr"), all_of(outcs)) %>%
      ggpairs(
        ., 
        progress = FALSE, 
        ggplot2::aes(
          colour = LBP_frequency_year, 
          fill = LBP_frequency_year, 
          alpha = 0.25
        )
      ) + 
      theme_bw() +
      scale_colour_viridis_d()+
      scale_fill_viridis_d()
  })
}



# ---- outcome0 ----

bpd <-
  bpd %>% 
  mutate(lbp_occurr = as.integer(LBP_sufferer == "yes"))

(this_outcome <- "lbp_occurr")

# (mod_form_null <-as.formula(paste0(this_outcome, " ~ ", rhs_formula)))
(mod_form_ilrs <-as.formula(paste0(this_outcome, " ~ ", rhs_formula, " + ilr")))

table(bpd[, this_outcome], useNA = "ifany")





# logistic regression model __with__ ilrs
bpd_occurr_ilrs <- glm(mod_form_ilrs, data = bpd, family = binomial())
summary(bpd_occurr_ilrs)
Anova(bpd_occurr_ilrs)





# ---- outcome0_diag ----

### check binned residuals are acceptable
# From the help file:
# Binned residual plots are achieved by “dividing the data into categories 
# (bins) based on their fitted values, and then plotting the average residual 
# versus the average fitted value for each bin.” (Gelman, Hill 2007: 97). 
# If the model were true, one would expect about 95% of the residuals to 
# fall inside the error bounds. 
bin_res_overall <- binned_residuals(bpd_occurr_ilrs)
bin_res_overall
plot(bin_res_overall)






# ---- outcome0_pred ----


# create dataset for predictions
newdata <- 
  bpd %>%
  dplyr::select(all_of(pred_covs), ilr) %>%
  distinct(pick(all_of(pred_covs)), .keep_all = TRUE) %>%
  arrange(pick(all_of(pred_covs)))

mean_ilr <- mean(bpd$ilr)
dev_null <- foreach(i = 1:nrow(newdata)) %do% {
  newdata$ilr[i, ] <- mean_ilr
}

# make preds and then put in long format for ggplot
predictions_probs <-
  cbind(
    `P(LBP)` = predict(bpd_occurr_ilrs, newdata, type = "response"),
    newdata
  ) %>%
  dplyr::select(-ilr) 

# predictions_probs


## model predictions for specific values
predictions_probs %>%
  dplyr::filter(
    # sex == "1_female", 
    stress == "1_normal", 
    smoking == "2_nonsmoker", 
    education == "2_higher", 
    # ses == "2_middle"
  ) %>%
  ggplot(., aes(age, `P(LBP)`, group = bmi)) +
  geom_line(aes(colour = bmi), linewidth = 1) +
  geom_point(aes(colour = bmi), size = 2) +
  facet_grid(sex~ ses, labeller = label_both) +
  labs(x = "age", y = "estimated P(lower back pain sufferer)") +
  theme_bw() +
  scale_color_manual(values = c("orange2", "turquoise", "purple")) +
  theme(axis.text.x = element_text(angle = 45))




# create a RHS of regression equation dataset for time-reallocation
predict_basis <-
  bpd %>%
  dplyr::select(all_of(pred_covs), all_of(pred_comps)) %>%
  dplyr::filter(
    age == "2_middle",
    sex == "1_female", 
    stress == "1_normal", 
    smoking == "2_nonsmoker", 
    education == "2_higher", 
    ses == "2_middle",
    bmi == "2_normal"
  ) 
  
### continuous scenario  
# predict_basis$age <- mean(predict_basis$age)
  
(predict_basis <-
  predict_basis %>%
  distinct(across(all_of(pred_covs)), .keep_all = TRUE) %>%
  as.data.frame())




# compositional mean: geometric mean to closure
# (comp_mean <- mean(acomp(bpd[, pred_comps])))
(comp_mean <- calc_comp_mean(bpd[, pred_comps], clo_val = 1440))
predict_basis0 <- predict_basis
predict_basis0[, pred_comps] <- comp_mean

predict_basis0

# +15 minutes to Time_MVPA and -15 minutes from Time_Sedentary 
comp_mean_changed <- comp_mean
comp_mean_changed["Time_MVPA"] <- comp_mean_changed["Time_MVPA"] + 15
comp_mean_changed["Time_Sedentary"] <- comp_mean_changed["Time_Sedentary"] - 15
# check
comp_mean_changed - comp_mean

predict_basis1 <- predict_basis
predict_basis1[, pred_comps] <- comp_mean_changed

pred_df <- rbind(predict_basis0, predict_basis1)
pred_df <- add_ilrs_to_data(pred_df, comp_vars = pred_comps, sbp_matrix = sbp1)
pred_df
predict(bpd_occurr_ilrs, pred_df, type = "link")
# ratio of odds ratios
exp(diff(predict(bpd_occurr_ilrs, pred_df, type = "link")))


get_pred_diff <- function(mod, new_dat) {
  log_odds_pred <- predict(mod, new_dat, type = "link")
  odds_ratio_ratio <- exp(log_odds_pred[2] - log_odds_pred[1])
  return(odds_ratio_ratio)
}
(est_v1 <- get_pred_diff(bpd_occurr_ilrs, pred_df))

fit_mod_boot <- function(data, i, pred_dat) {
  
  this_dat <- data[i, ]
  this_logis <- glm(mod_form_ilrs, data = this_dat, family = binomial())
  est <- get_pred_diff(this_logis, new_dat = pred_dat)
  return(est)
  
}

### CI method #1 (bootstrapping):
alpha <- 0.05
(ci_v1 <- 
  c(
    est = est_v1,
    quantile(
      boot(bpd, fit_mod_boot, R = 100, pred_dat = pred_df)$t,
      c(alpha / 2, 1 - alpha / 2)
    )))


### alternative CI method #2 (Wald approximation - re-transformed):
pred_df[, "ilr"]
diff(pred_df[, "ilr"])
x_0_red <- matrix(as.numeric(diff(pred_df[, "ilr"])), nrow = 1)
x_0_red
betas <- coef(bpd_occurr_ilrs)
nms_kp <- grepl("^ilr", names(betas))
betas_red <- as.matrix(betas[nms_kp])
Sigma <- stats::vcov(bpd_occurr_ilrs)
nms_kp <- grepl("^ilr", colnames(Sigma))
sigma_red <- Sigma[nms_kp, nms_kp]
sigma_red
est_red <- x_0_red %*% betas_red
se_red <- sqrt(x_0_red %*% sigma_red %*% t(x_0_red))
z_star <- qnorm(0.975)
(ci_v2 <- 
  exp(c(
    est = est_red, 
    lo = est_red - z_star * se_red, 
    hi = est_red + z_star * se_red
  )))

### alternative CI method #3 (delta method)
# (first order approximation, although still linear combin of param ests):
approx_ci <-
  deltaMethod(
    bpd_occurr_ilrs, 
    "-0.2418864 * `ilrilr(++--)` + 0.02454151 * `ilrilr(+-..)` + -0.3175375 * `ilrilr(..+-)`"
  )
(ci_v3 <- 
  exp(c(
    est = approx_ci[["Estimate"]], 
    lo = approx_ci[["2.5 %"]], 
    hi = approx_ci[["97.5 %"]]
  )))


### compare CIs
kable(rbind(ci_v1, ci_v2, ci_v3))


do_multi_realloc <- function(mod, basis_data, timeusenames, time_changes, sbp_matrix = sbp1) {
  
  
  x0 <- basis_data
  
  plot_dat <- 
    foreach(i = 1:length(timeusenames), .combine = bind_rows) %do% {
      print(paste("i: ", i))
      foreach(j = 1:length(timeusenames), .combine = bind_rows) %do% {
        print(paste("  j: ", j))
        foreach(d = 1:length(time_changes), .combine = bind_rows) %do% {
          print(paste("    d: ", d))
          
          timeuse_to <- timeusenames[i]
          timeuse_from <- timeusenames[j]
          change_time <- time_changes[d]
          
          proposed_change_1 <- x0[timeuse_to] + change_time
          proposed_change_2 <- x0[timeuse_from] - change_time
          
          if (timeuse_to == timeuse_from) {
            NULL # reallocation exceeds 0 or max time
          } else if ((proposed_change_1 < 0) | (proposed_change_1 > 1440)) {
            NULL # reallocation exceeds 0 or max time
          } else if ((proposed_change_2 < 0) | (proposed_change_2 > 1440)) {
            NULL # reallocation exceeds 0 or max time
          } else {
              
            x1 <- x0
            x1[timeuse_to] <- x1[timeuse_to] + change_time
            x1[timeuse_from] <- x1[timeuse_from] - change_time
            
            pred_df <- rbind(x0, x1)
            pred_df <- add_ilrs_to_data(pred_df, comp_vars = timeusenames, sbp_matrix = sbp_matrix)
            
            ratio_of_odds_ratios <- get_pred_diff(mod, pred_df)
            
            bootstrapped_ests <- boot(bpd, fit_mod_boot, R = 1000, pred_dat = pred_df)$t
            ci_est <- quantile(as.numeric(bootstrapped_ests), c(alpha / 2, 1 - alpha / 2))
              
            tibble(
              to = timeuse_to, 
              from = timeuse_from, 
              change_time = change_time,
              ratio_of_odds_ratios = ratio_of_odds_ratios,
              ci_lo = ci_est[1],
              ci_hi = ci_est[2]
            )
          }
        }
      }
    }
  
  plot_dat$to <- factor(plot_dat$to, levels = timeusenames)
  plot_dat$from <- factor(plot_dat$from, levels = timeusenames)
  
  return(plot_dat)
  
}




# takes ~25 min (single core)

### Uncomment to generate bootstrapping

# tic()
# set.seed(1234)
# realloc_plot_data <-
#   do_multi_realloc(
#     bpd_occurr_ilrs,
#     predict_basis0,
#     pred_comps,
#     seq(-30, 30, by = 10)
#   )
# saveRDS(realloc_plot_data, file = "res/logistic_realloc_boot_res.rda")
# toc()


realloc_plot_data <- readRDS(file = "res/logistic_realloc_boot_res.rda")


levels(realloc_plot_data$to) <- paste0(levels(realloc_plot_data$to), "+Delta")
levels(realloc_plot_data$from) <- paste0(levels(realloc_plot_data$from), "-Delta")

ggplot(realloc_plot_data) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 1, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3) +
  geom_line(aes(x = change_time , y = ratio_of_odds_ratios, col = to)) +
  geom_point(aes(x = change_time , y = ratio_of_odds_ratios, col = to), size = 1) +
  facet_grid(from ~ to, labeller = label_parsed) +
  theme_bw() +
  scale_colour_manual(values = c("darkorange", "purple", "cyan4", "dodgerblue")) +
  scale_fill_manual(values = c("darkorange", "purple", "cyan4", "dodgerblue")) +
  labs(
    x = paste0("Change/delta in composition (mins)"),
    y = paste0("Ratio of odds-ratios (after reallocation:before reallocation)")
  ) +
  theme(legend.position = "none")


ggsave(
  filename = "fig/lbp_occur_logistic_odds.png",
  dpi = 600, # print quality
  width = 10,
  height = 10
)


# ---- update_data ----

bpd_yes <- bpd %>% dplyr::filter(LBP_sufferer == "yes")
nrow(bpd)
nrow(bpd_yes)


bpd_clean_yes <- bpd_clean %>% dplyr::filter(LBP_sufferer == "yes")




# ---- outcome1 ----

(this_outcome <- outcs[1])

# (mod_form_null <-as.formula(paste0(this_outcome, " ~ ", rhs_formula)))
(mod_form_ilrs <-as.formula(paste0(this_outcome, " ~ ", rhs_formula, " + ilr")))

table(bpd_yes[, this_outcome], useNA = "ifany")

bpd_yes[[this_outcome]] <- fct_drop(bpd_yes[[this_outcome]])
table(bpd_yes[, this_outcome], useNA = "ifany")


## model without ilrs
# bpd_ordinal_null <- polr(mod_form_null, data = bpd, Hess = TRUE, method = "logistic")
# summary(bpd_ordinal_null)


## model __with__ ilrs
bpd_ordinal_ilrs <- polr(mod_form_ilrs, data = bpd_yes, Hess = TRUE, method = "logistic")
summary(bpd_ordinal_ilrs)
Anova(bpd_ordinal_ilrs)

# pr <- profile(bpd_ordinal_ilrs)
# confint(pr)
# plot(pr)
# pairs(pr)


est_ci_df <- cbind(est = coef(bpd_ordinal_ilrs), confint(bpd_ordinal_ilrs)) # profiled CIs
kable(est_ci_df, digits = 3)      # these are the log-odds scale estimates (and CI)
kable(exp(est_ci_df), digits = 3) # these are the odds ratios (and approx CIs)


# ---- outcome1_diag ----


# deviance test 
g2 <- deviance(bpd_ordinal_ilrs)
df <- df.residual(bpd_ordinal_ilrs)
1 - pchisq(g2, df)


with(bpd_yes, 
  table(
    LBP_frequency_year, 
    as.numeric(LBP_frequency_year), 
    useNA = "ifany"
  )
)


## checking parallel slopes assumptions can be done by fitting successive logistic regressions
## while creating a binary outcome using different thresholds of the ordinal outcome
### (note the rhs/linear predictor is negative so coefs should be approx same
### as main model except negative)

# e.g. this is a single logistic regression
coef(glm(
  I(as.numeric(LBP_frequency_year) <= 1) ~ 
    age + sex + bmi + stress + smoking + education + ses + ilr, 
  family = "binomial", 
  data = bpd_yes
))

# this is running multiple logistic regressions
## we want to see the coefficients to be roughly the same (intercepts and negative coefs - see above)

### note that the below shows there may be reason to include an age variable that has
### non-constant coefficient for each level of the outcome (or subgroup analyses for each age cohort)
### we can see this because the age coefs increase/decrease monotonically

foreach(i = 1:(length(levels(bpd_yes$LBP_frequency_year)) - 1), .combine = cbind) %do% {
  log_coefs <-
    coef(glm(
      I(as.numeric(LBP_frequency_year) <= i) ~ 
         age + sex + bmi + stress + smoking + education + ses + ilr, 
      family = "binomial", 
      data = bpd_yes
    ))
  log_coefs <- as.data.frame(log_coefs)
  colnames(log_coefs) <- paste0("logit(P(Y<=", i, "))")
  log_coefs
} %>%
  kable(., digits = 2)





# ---- outcome1_pred_a ----


# create dataset for predictions
newdata <- 
  bpd_yes %>%
  dplyr::select(all_of(pred_covs), ilr) %>%
  distinct(pick(all_of(pred_covs)), .keep_all = TRUE) %>%
  arrange(pick(all_of(pred_covs)))

mean_ilr <- mean(bpd_yes$ilr)
dev_null <- foreach(i = 1:nrow(newdata)) %do% {
  newdata$ilr[i, ] <- mean_ilr
}

# make preds and then put in long format for ggplot
predictions_probs <-
  cbind(
    predict(bpd_ordinal_ilrs, newdata, type = "probs"),
    newdata
  ) %>%
  dplyr::select(-ilr) %>%
  pivot_longer(
    cols = -all_of(pred_covs),
    names_to = "outcome",
    values_to = "P(outc)"
  )

predictions_probs


## model predictions for specific values
predictions_probs %>%
  dplyr::filter(
    # sex == "1_female", 
    stress == "1_normal", 
    smoking == "2_nonsmoker", 
    education == "2_higher", 
    # ses == "2_middle"
  ) %>%
  ggplot(., aes(outcome, `P(outc)`)) +
  geom_line(aes(colour =  age, group =  age), linewidth = 1) +
  geom_point(aes(colour =  age), size = 2) +
  facet_grid(sex * bmi ~ ses, labeller = label_both) +
  labs(x = "Outcome: back pain freq", y = "estimated P(outcome)") +
  theme_bw() +
  scale_colour_viridis_d()
  
  




# ---- outcome1_pred_b ----


# create a RHS of regression equation dataset for time-reallocation
predict_basis <-
   bpd_yes %>%
   dplyr::select(all_of(pred_covs), all_of(pred_comps)) %>%
   dplyr::filter(
     age == "2_middle",
     sex == "1_female", 
     stress == "1_normal", 
     smoking == "2_nonsmoker", 
     education == "2_higher", 
     ses == "2_middle",
     bmi == "2_normal"
   ) 
 
### continuous situation
# predict_basis$age <- mean(predict_basis$age)
 
(predict_basis <-
  predict_basis %>%
   distinct(across(all_of(pred_covs)), .keep_all = TRUE) %>%
   as.data.frame())



# compositional mean: geometric mean to closure
# (comp_mean <- mean(acomp(bpd_yes[, pred_comps])))
(comp_mean <- calc_comp_mean(bpd_yes[, pred_comps], clo_val = 1440))
predict_basis0 <- predict_basis
predict_basis0[, pred_comps] <- comp_mean

predict_basis0

# +15 minutes to Time_MVPA and -15 minutes from Time_Sedentary 
comp_mean_changed <- comp_mean
comp_mean_changed["Time_MVPA"] <- comp_mean_changed["Time_MVPA"] + 15
comp_mean_changed["Time_Sedentary"] <- comp_mean_changed["Time_Sedentary"] - 15
# check
comp_mean_changed - comp_mean

predict_basis1 <- predict_basis
predict_basis1[, pred_comps] <- comp_mean_changed

pred_df <- rbind(predict_basis0, predict_basis1)
pred_df <- add_ilrs_to_data(pred_df, comp_vars = pred_comps, sbp_matrix = sbp1)
pred_df <- pred_df[, !(colnames(pred_df) %in% pred_comps)] # get rid of compositions
pred_df

# model.matrix(formula(bpd_ordinal_ilrs), data = cbind(LBP_frequency_year = 0, pred_df))


df <- bpd_yes[, attr(formula(bpd_ordinal_ilrs), "term.labels")]
# this is a list of levels for each factor in the original df (after applying factor funciton)
xlevs <- lapply(df[,sapply(df, is.character), drop = F], function(j) {
  levels(factor(j))
})

# calling "xlev = " builds out a model.matrix with identical levels as the original df
mm_new <- model.matrix( ~ ., data = pred_df, xlev = xlevs)
colnames(mm_new)
mm_new <- mm_new[, -1] # remove intercept


colnames(mm_new)[grepl("^ilr", colnames(mm_new))] <- paste0("ilr", create_ilr_names(sbp1))
# colnames(mm_new)
# don't need intercept  # c("(Intercept)" = 1, coef(bpd_ordinal_ilrs))
betas <- as.matrix(coef(bpd_ordinal_ilrs)) # should be col matrix
# rownames(betas)

if (!all(colnames(mm_new) == rownames(betas))) {
  stop("design and parameter est matrices non-conform")
}

# note as linear predictor is taken from the K intercepts the ratio of odds ratios is flipped
# i.e. after:before of odds is calculated as exp(before_log_odds / after_log_odds)
preds <- mm_new %*% betas
exp(preds[1] - preds[2])

# check manual calcs agree with model
mm_old <- model.matrix( ~ ., data = df, xlev = xlevs)
mm_old <- mm_old[, -1] # remove intercept
# colnames(mm_old)

# model and manual calcs agree?
# note that bpd_ordinal_ilrs$lp are the eta/linear predictor that is taken 
# away from the xi_k intercept
all(abs(as.numeric(mm_old %*% betas) - bpd_ordinal_ilrs$lp) < 1e-9)

# bpd_ordinal_ilrs$lp # linear predictor

get_pred_diff <- function(mod, new_dat) {
  betas_ <- as.matrix(coef(mod))
  if (!all(colnames(new_dat) == rownames(betas_))) {
    print(paste(paste(colnames(new_dat), collapse = "|"), "vs", paste(rownames(betas_), collapse = "|")))
    stop("design and parameter est matrices non-conform")
  }
  log_odds_pred <- as.numeric(new_dat %*% betas_)
  # note reversal of order (see above)
  odds_ratio_ratio <- exp(log_odds_pred[1] - log_odds_pred[2])
  return(odds_ratio_ratio)
}
(est_v1 <- get_pred_diff(bpd_ordinal_ilrs, mm_new))

fit_mod_boot <- function(data, i, pred_dat) {
  
  this_dat <- data[i, ]
  this_ordinal <- polr(mod_form_ilrs, data = this_dat, Hess = TRUE, method = "logistic")

  df <- this_dat[, attr(formula(this_ordinal), "term.labels")]
  # this is a list of levels for each factor in the original df (after applying factor funciton)
  xlevs <- lapply(df[,sapply(df, is.character), drop = F], function(j) {
    levels(factor(j))
  })
  
  
  # calling "xlev = " builds out a model.matrix with identical levels as the original df
  mm_new <- model.matrix( ~ ., data = pred_dat, xlev = xlevs)
  mm_new <- mm_new[, -1] # remove intercept
  # make sure ilr colnames are legit/match coeffs
  colnames(mm_new)[grepl("^ilr", colnames(mm_new))] <- paste0("ilr", create_ilr_names(sbp1))
  colnames(mm_new)
  
  est <- get_pred_diff(this_ordinal, new_dat = mm_new)
  return(est)
  
}

### CI method #1 (bootstrapping):
alpha <- 0.05
(ci_v1 <- 
    c(
      est = est_v1,
      quantile(
        boot(bpd_yes, fit_mod_boot, R = 100, pred_dat = pred_df)$t,
        c(alpha / 2, 1 - alpha / 2)
      )))


### alternative CI method #2 (Wald approximation - re-transformed):
pred_df[, "ilr"]
diff(pred_df[, "ilr"])
# x_0_red <- matrix(- as.numeric(diff(pred_df[, "ilr"])), nrow = 1)
x_0_red <- matrix(as.numeric(pred_df[1, "ilr"] - pred_df[2, "ilr"]), nrow = 1)
x_0_red
betas <- coef(bpd_ordinal_ilrs)
nms_kp <- grepl("^ilr", names(betas))
betas_red <- as.matrix(betas[nms_kp])
Sigma <- stats::vcov(bpd_ordinal_ilrs)
nms_kp <- grepl("^ilr", colnames(Sigma))
sigma_red <- Sigma[nms_kp, nms_kp]
sigma_red
est_red <- x_0_red %*% betas_red
se_red <- sqrt(x_0_red %*% sigma_red %*% t(x_0_red))
z_star <- qnorm(0.975)
(ci_v2 <- 
    exp(c(
      est = est_red, 
      lo = est_red - z_star * se_red, 
      hi = est_red + z_star * se_red
    )))

### alternative CI method #3 (delta method)
# (first order approximation, although still linear combin of param ests):
as.numeric(x_0_red)
(g_form <- paste(
  paste(
    as.numeric(x_0_red), 
    "*", 
    c("`ilrilr(++--)`", "`ilrilr(+-..)`", "`ilrilr(..+-)`")
  ), 
  collapse = " + "
))
approx_ci <-deltaMethod(bpd_ordinal_ilrs, g_form)
(ci_v3 <- 
    exp(c(
      est = approx_ci[["Estimate"]], 
      lo = approx_ci[["2.5 %"]], 
      hi = approx_ci[["97.5 %"]]
    )))


### compare CIs
kable(rbind(ci_v1, ci_v2, ci_v3))


do_multi_realloc <- function(mod, basis_data, timeusenames, time_changes, sbp_matrix = sbp1) {
  
  
  x0 <- basis_data
  
  plot_dat <- 
    foreach(i = 1:length(timeusenames), .combine = bind_rows) %do% {
      print(paste("i: ", i))
      foreach(j = 1:length(timeusenames), .combine = bind_rows) %do% {
        print(paste("  j: ", j))
        foreach(d = 1:length(time_changes), .combine = bind_rows) %do% { # %dopar%
          print(paste("    d: ", d))
          
          timeuse_to <- timeusenames[i]
          timeuse_from <- timeusenames[j]
          change_time <- time_changes[d]
          
          proposed_change_1 <- x0[timeuse_to] + change_time
          proposed_change_2 <- x0[timeuse_from] - change_time
          
          if (timeuse_to == timeuse_from) {
            NULL # reallocation exceeds 0 or max time
          } else if ((proposed_change_1 < 0) | (proposed_change_1 > 1440)) {
            NULL # reallocation exceeds 0 or max time
          } else if ((proposed_change_2 < 0) | (proposed_change_2 > 1440)) {
            NULL # reallocation exceeds 0 or max time
          } else {
            
            x1 <- x0
            x1[timeuse_to] <- x1[timeuse_to] + change_time
            x1[timeuse_from] <- x1[timeuse_from] - change_time
            
            pred_df <- rbind(x0, x1)
            pred_df <- add_ilrs_to_data(pred_df, comp_vars = timeusenames, sbp_matrix = sbp_matrix)
            
            ### alternative CI method #3 (delta method)
            # x_0_red <- -as.numeric(diff(pred_df[, "ilr"]))
            x_0_red <- as.numeric(pred_df[1, "ilr"] - pred_df[2, "ilr"])
            # (first order approximation, although still linear combin of param ests):
            (g_form <- paste(
              paste(
                x_0_red, 
                "*", 
                c("`ilrilr(++--)`", "`ilrilr(+-..)`", "`ilrilr(..+-)`")
              ), 
              collapse = " + "
            ))
            approx_ci <-deltaMethod(bpd_ordinal_ilrs, g_form)
            this_ci <- 
                exp(c(
                  est = approx_ci[["Estimate"]], 
                  lo = approx_ci[["2.5 %"]], 
                  hi = approx_ci[["97.5 %"]]
                ))
            
            
            ### bootstrapping takes too long
            # pred_df <- pred_df[, !(colnames(pred_df) %in% timeusenames)] # get rid of compositions
            # 
            # df <- bpd_yes[, attr(formula(bpd_ordinal_ilrs), "term.labels")]
            # # this is a list of levels for each factor in the original df (after applying factor funciton)
            # xlevs <- lapply(df[,sapply(df, is.character), drop = F], function(j) {
            #   levels(factor(j))
            # })
            # 
            # # calling "xlev = " builds out a model.matrix with identical levels as the original df
            # mm_new <- model.matrix( ~ ., data = pred_df, xlev = xlevs)
            # mm_new <- mm_new[, -1] # remove intercept
            # # make sure ilr colnames are legit/match coeffs
            # colnames(mm_new)[grepl("^ilr", colnames(mm_new))] <- paste0("ilr", create_ilr_names(sbp1))
            # 
            # ratio_of_odds_ratios <- get_pred_diff(mod, new_dat = mm_new)
            # bootstrapped_ests <- boot(bpd_yes, fit_mod_boot, R = 10, pred_dat = pred_df)$t
            # ci_est <- quantile(as.numeric(bootstrapped_ests), c(alpha / 2, 1 - alpha / 2))
            
            tibble(
              to = timeuse_to, 
              from = timeuse_from, 
              change_time = change_time,
              ratio_of_odds_ratios = this_ci["est"],
              ci_lo = this_ci["lo"],
              ci_hi = this_ci["hi"]
            )
          }
        }
      }
    }
  
  plot_dat$to <- factor(plot_dat$to, levels = timeusenames)
  plot_dat$from <- factor(plot_dat$from, levels = timeusenames)
  
  return(plot_dat)
  
}


# takes ~ 3h (single core) for bootstrapping
# takes ~ 4sec (single core) for delta/wald method



### Uncomment to generate bootstrapping

# tic()
# set.seed(1234)

  # # library("doParallel")
  # # no_cores <- detectCores() - 1 # Calculate the number of cores (leave one free)
  # # cl <- makeCluster(no_cores) # Create clusters
  # # registerDoParallel(cl) # and register

# realloc_plot_data <-
#   do_multi_realloc(
#     bpd_ordinal_ilrs,
#     predict_basis0,
#     pred_comps,
#     seq(-30, 30, by = 10)
#   )

  # # # close para comp
  # # stopCluster(cl) 

# saveRDS(realloc_plot_data, file = "res/ordinal_realloc_wald_res.rda")
# toc()

# saveRDS(realloc_plot_data, file = "res/ordinal_realloc_boot_res.rda")
# realloc_plot_data <- readRDS(file = "res/ordinal_realloc_boot_res.rda")


realloc_plot_data <- readRDS(file = "res/ordinal_realloc_wald_res.rda")


levels(realloc_plot_data$to) <- paste0(levels(realloc_plot_data$to), "+Delta")
levels(realloc_plot_data$from) <- paste0(levels(realloc_plot_data$from), "-Delta")

ggplot(realloc_plot_data) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 1, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3) +
  geom_line(aes(x = change_time , y = ratio_of_odds_ratios, col = to)) +
  geom_point(aes(x = change_time , y = ratio_of_odds_ratios, col = to), size = 1) +
  facet_grid(from ~ to, labeller = label_parsed) +
  theme_bw() +
  scale_colour_manual(values = c("darkorange", "purple", "cyan4", "dodgerblue")) +
  scale_fill_manual(values = c("darkorange", "purple", "cyan4", "dodgerblue")) +
  labs(
    x = paste0("Change/delta in composition (mins)"),
    y = paste0("Ratio of odds-ratios (after reallocation:before reallocation)"),
    subtitle = "Note that odds ratios relate to the probability of having _decreased_ frequency (per year) of pain"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = "fig/lbp_freq_ordinal_odds_v1.png",
  dpi = 600, # print quality
  width = 10,
  height = 10
)





time_lvls <- gsub("Time_", "", pred_comps)


rep_char <- function(n, char = " ") paste(rep(char, n), collapse = "")
rep_char(3)
rep_char(0)

rep_char <- Vectorize(rep_char, vectorize.args = "n")
rep_char(0:7)


pd2 <-
  realloc_plot_data %>%
  mutate(
    to = gsub("Time_", "", to),
    from = gsub("Time_", "", from),
    to = gsub("+Delta", "", to, fixed = TRUE),
    from = gsub("-Delta", "", from, fixed = TRUE),
    to_len = nchar(to),
    to_max = max(to_len),
    from_len = nchar(from),
    from_max = max(from_len),
    to_pad =  rep_char(pmax(0, from_max - to_len)),
    from_pad = rep_char(pmax(0, to_max - from_len)),
    to = factor(to, levels = time_lvls),
    from = factor(from, levels = time_lvls),
    to_num = as.numeric(to),
    from_num = as.numeric(from)
  ) %>%
  dplyr::filter(to_num > from_num) %>%
  mutate(
    ratio_of_odds_ratios = 1 / ratio_of_odds_ratios,
    tmp = 1 / ci_lo,
    ci_lo = 1 / ci_hi,
    ci_hi = tmp,
    # from_to = paste0("     ", "+", from, rep_char(10), from_pad, "\u2194", to_pad, rep_char(10), "+", to, "     ")
    from_to = paste0("+", from, rep_char(13), from_pad, "", to_pad, rep_char(13), "+", to)
  ) %>%
  arrange(from, to)

unique(pd2$from_to)
pd2$from_to <- factor(pd2$from_to, levels = unique(pd2$from_to))

this_breaks <- seq(-30, 30, 10)
this_labs <- sprintf("+%2.0f", abs(seq(-30, 30, 10)))
this_labs[this_labs == "+ 0"] <- ""
this_labs

ggplot(pd2) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 1, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3, col = NA, fill = "cyan4") +
  geom_line(aes(x = change_time , y = ratio_of_odds_ratios, col = to), col = "cyan4") +
  geom_point(aes(x = change_time , y = ratio_of_odds_ratios, col = to), size = 1, col = "cyan4") +
  facet_wrap(~ from_to, labeller = label_bquote(.(from_to))) +
  theme_bw() +
  scale_x_continuous(breaks = this_breaks, labels = this_labs) +
  labs(
    x = paste0("Reallocation between pair of compositional parts (minutes)"),
    y = paste0("Odds-ratio of increased LBP frequency (after reallocation:before reallocation)")
    # subtitle = "Note that odds ratios relate to the probability of having _increased_ frequency (per year) of pain"
  ) +
  theme(
    legend.position = "none",
    text = element_text(family = "serif"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text =  element_text(size = 10),
    axis.title =  element_text(size = 12)
  )


ggsave(filename = "fig/lbp_freq_ordinal_odds_v2.png", width = 14, height = 9, dpi = 600)
# ggsave(filename = "fig/lbp_freq_ordinal_odds.pdf", width = 10, height = 8)



 # ---- outcome1_pred_not_use ----



# logitP(Y≤k∣x)=ζ_k−η
# zeta_{1-7days|8-30days} = -0.2910
# eta = 0.3184 + -0.1786 + -0.1110 + -0.3463 + -0.3017
coef(bpd_ordinal_ilrs)
# summary(bpd_ordinal_ilrs)
bpd_ordinal_ilrs$zeta


# bpd_ordinal_ilrs$lp
# p_0 <- predict(bpd_ordinal_ilrs, pred_df, type = "prob")
# (lodr <- log(p_0 / (1-p_0)))
# # ratio of odds ratios
# exp(apply(lodr, 2, diff))
# # predicted class  argmin_k{abs(zeta_k - eta)}?
# predict(bpd_ordinal_ilrs, type = "class")[1:3]
# p_m <- matrix(rep(bpd_ordinal_ilrs$lp, 5), ncol = 5)
# co_m <- matrix(rep(c(bpd_ordinal_ilrs$zeta, 0), nrow(p_m)), ncol = 5, byrow = TRUE)
# apply(abs(p_m - co_m), 1, which.min)[1:3]
# 
# table(
#   predict(bpd_ordinal_ilrs, type = "class"), 
#   apply(abs(p_m - co_m), 1, which.min) # + 1) %% 5
# )


# ---- outcome2 ----

(this_outcome <- outcs[2])

# (mod_form_null <-as.formula(paste0(this_outcome, " ~ ", rhs_formula)))
(mod_form_ilrs <- as.formula(paste0(this_outcome, " ~ ", rhs_formula, " + ilr")))

hist(bpd_yes[[this_outcome]])





lbp_intensity_lm <- lm(mod_form_ilrs, data = bpd_yes)
summary(lbp_intensity_lm)
car::Anova(lbp_intensity_lm)

### THis is logodds transform of outcome, not a good fit
# # move extreme values off boundary
# bpd_yes$intensity <- bpd_yes$LBP_intensity_year
# bpd_yes$intensity[bpd_yes$LBP_intensity_year < 0.5] <- 0.5
# bpd_yes$intensity[bpd_yes$LBP_intensity_year > (100 - 0.5)] <- 100 - 0.5
# bpd_yes$logodds_intensity <- with(bpd_yes, log((intensity / 100) / (1 - intensity / 100)))
# bpd_yes$intensity <- NULL
# (mod_form_logodds_ilrs <- as.formula(paste0("logodds_intensity ~ ", rhs_formula, " + ilr")))
# lbp_intensity_logodds_lm <- lm(mod_form_logodds_ilrs, data = bpd_yes)
# summary(lbp_intensity_logodds_lm)
# car::Anova(lbp_intensity_logodds_lm)
# check_model(lbp_intensity_logodds_lm)

lbp_intensity_pois <- glm(mod_form_ilrs, family = "poisson", data = bpd_yes)
summary(lbp_intensity_pois)

# check the goodness of fit test not significant, 
# p > 0.05: indicates model fit the data
# p < 0.05: indicates model DOES NOIT fit the data
with(
  lbp_intensity_pois, 
  cbind(
    res.deviance = deviance, 
    df = df.residual,
    p = pchisq(deviance, df.residual, lower.tail = FALSE)
  )
)



lbp_intensity_nb <- glm.nb(mod_form_ilrs, data = bpd_yes)
summary(lbp_intensity_nb)
car::Anova(lbp_intensity_nb)

(est <- cbind(Estimate = coef(lbp_intensity_nb), confint(lbp_intensity_nb)))
exp(est) 

# likelihood ratio test 
# dispersion parameter check is it equal to zero (if not, NB mod preferred)
pchisq(
  2 * (logLik(lbp_intensity_nb) - logLik(lbp_intensity_pois)), 
  df = 1, 
  lower.tail = FALSE
)

par(mfrow = c(2, 1))
hist(bpd_yes$LBP_intensity_year, xlim = c(0, 160), breaks = 10, 
     main = "observed lower back pain intensity (0-100)")
hist(simulate(lbp_intensity_nb)$sim_1, xlim = c(0, 160), breaks = 20, 
     main = "neg binomial predicted values (0-inf)")
par(mfrow = c(1, 1))


# Pain intensity could be categorised as (Boonstra et al., 2014):
# - no pain (0)
# - mild pain (1-38)
# - moderate pain (39-57)
# - severe pain (58-100)

bpd_yes$intens_ord <- 
  cut(
    bpd_yes[[this_outcome]], 
    breaks = c(-1, 0, 38, 57, 101),
    labels = c(
      "no pain (0)", "mild pain (1-38)", 
      "moderate pain (39-57)", "severe pain (58-100)"
    )
  )
class(bpd_yes$intens_ord)
table(
  bpd_yes$intens_ord, 
  cut(bpd_yes[[this_outcome]], breaks = c(-1, 0, 38, 57, 101)), 
  useNA = "ifany"
)

(mod_form_ord_ilrs <-as.formula(paste0("intens_ord ~ ", rhs_formula, " + ilr")))


## model __with__ ilrs
bpd_intens_ord_ilrs <- polr(mod_form_ord_ilrs, data = bpd_yes, Hess = TRUE, method = "logistic")
summary(bpd_intens_ord_ilrs)
Anova(bpd_intens_ord_ilrs)

# profiled CIs
est_ci_df <- cbind(est = coef(bpd_intens_ord_ilrs), confint(bpd_intens_ord_ilrs)) 
kable(est_ci_df, digits = 3)      # these are the log-odds scale estimates (and CI)
kable(exp(est_ci_df), digits = 3) # these are the odds ratios (and approx CIs)


# ---- outcome2_diag ----

## plain linear model
check_model(lbp_intensity_lm) # acceptable?

## Poisson regression (bad)
check_model(lbp_intensity_pois) # horrible

## Negative Binomial regression
check_model(
  lbp_intensity_nb,
  check = c("pp_check", "homogeneity", "outliers", "vif")
)



plot(
  predict(lbp_intensity_nb, type = "link"),
  residuals(lbp_intensity_nb, type = "deviance")
)


deviance(lbp_intensity_nb)/lbp_intensity_nb$df.residual





## Ordinal logistic regression (looks ok)
# this is running multiple logistic regressions
## we want to see the coefficients to be roughly the same EXCEPT for the 
## (intercept) values
foreach(i = 2:length(levels(bpd_yes$intens_ord)), .combine = cbind) %do% {
  log_coefs <-
    coef(glm(
      I(as.numeric(intens_ord) >= i) ~ 
        age + sex + bmi + stress + smoking + education + ses + ilr, 
      family = "binomial", 
      data = bpd_yes
    ))
  log_coefs <- as.data.frame(log_coefs)
  colnames(log_coefs) <- paste0("logit(P(Y>=", i, "))")
  log_coefs
} %>%
  kable(., digits = 2)


# ---- outcome2_pred_a ----

# create dataset for predictions
newdata <- 
  bpd_yes %>%
  dplyr::select(all_of(pred_covs), ilr) %>%
  distinct(pick(all_of(pred_covs)), .keep_all = TRUE) %>%
  arrange(pick(all_of(pred_covs)))

(mean_ilr <- mean(bpd_yes$ilr))
dev_null <- foreach(i = 1:nrow(newdata)) %do% {
  newdata$ilr[i, ] <- mean_ilr
}

# make preds and then put in long format for ggplot
predictions_intens <-
  cbind(
    pain_intens = predict(lbp_intensity_nb, newdata, type = "response"),
    newdata
  ) %>%
  dplyr::select(-ilr) 

head(predictions_intens)



# newdata2 <- cbind(newdata2, predict(lbp_nb, newdata2, type = "link", se.fit=TRUE))
# newdata2 <- within(newdata2, {
#   lbp_pred <- exp(fit)
#   LL <- exp(fit - 1.96 * se.fit)
#   UL <- exp(fit + 1.96 * se.fit)
# })



predictions_intens %>%
  dplyr::filter(
    # sex == "1_female", 
    stress == "1_normal", 
    smoking == "2_nonsmoker", 
    education == "2_higher", 
    # ses == "2_middle"
  ) %>%
  ggplot(., aes(age, pain_intens, group = bmi)) +
  geom_line(aes(colour = bmi), linewidth = 1) +
  geom_point(aes(colour = bmi), size = 2) +
  facet_grid(sex~ ses, labeller = label_both) +
  labs(x = "age", y = "estimated back pain intesity") +
  theme_bw() +
  scale_color_manual(values = c("darkorange", "turquoise", "purple")) +
  theme(axis.text.x = element_text(angle = 45))





# ---- outcome2_pred_b_abs ----



# create a RHS of regression equation dataset for time-reallocation
(predict_basis <-
    bpd_yes %>%
    dplyr::select(all_of(pred_covs), all_of(pred_comps)) %>%
    dplyr::filter(
      age == "2_middle",
      sex == "1_female", 
      stress == "1_normal", 
      smoking == "2_nonsmoker", 
      education == "2_higher", 
      ses == "2_middle",
      bmi == "2_normal"
    ) %>%
    distinct(across(all_of(pred_covs)), .keep_all = TRUE) %>%
    as.data.frame())





# compositional mean: geometric mean to closure
# (comp_mean <- mean(acomp(bpd_yes[, pred_comps])))
(comp_mean <- calc_comp_mean(bpd_yes[, pred_comps], clo_val = 1440))
predict_basis0 <- predict_basis
predict_basis0[, pred_comps] <- comp_mean

predict_basis0

# +15 minutes to Time_MVPA and -15 minutes from Time_Sedentary 
comp_mean_changed <- comp_mean
comp_mean_changed["Time_MVPA"] <- comp_mean_changed["Time_MVPA"] + 15
comp_mean_changed["Time_Sedentary"] <- comp_mean_changed["Time_Sedentary"] - 15
# check
comp_mean_changed - comp_mean

predict_basis1 <- predict_basis
predict_basis1[, pred_comps] <- comp_mean_changed

pred_df <- rbind(predict_basis0, predict_basis1)
pred_df <- add_ilrs_to_data(pred_df, comp_vars = pred_comps, sbp_matrix = sbp1)
pred_df
predict(lbp_intensity_nb, pred_df, type = "link")
# exponentiate difference in the log back pain intensity (ratio of back pain preds)
exp(diff(predict(lbp_intensity_nb, pred_df, type = "link")))
# abs difference in the mean back pain intensity
diff(predict(lbp_intensity_nb, pred_df, type = "response"))
(p_0 <- predict(lbp_intensity_nb, pred_df, type = "response"))
# % increase in pain intensity
(p_0[2]  - p_0[1]) / p_0[1]

# ratio version
get_pred_diff_rat <- function(mod, new_dat) {
  log_ratio_pred <- predict(mod, new_dat, type = "link")
  ratio_outc <- exp(log_ratio_pred[2] - log_ratio_pred[1])
  return(ratio_outc)
}
get_pred_diff_rat(lbp_intensity_nb, pred_df)

# absolute difference version
get_pred_diff_abs <- function(mod, new_dat) {
  log_ratio_pred <- predict(mod, new_dat, type = "response")
  ratio_outc <- log_ratio_pred[2] - log_ratio_pred[1]
  return(ratio_outc)
}
get_pred_diff_abs(lbp_intensity_nb, pred_df)

# wrapper:
get_pred_diff <- function(mod, new_dat, type = "abs") {
  if (type == "abs") {
    return(get_pred_diff_abs(mod = mod, new_dat = new_dat))
  } else if (type == "rat") {
    return(get_pred_diff_rat(mod = mod, new_dat = new_dat))
  } else {
    stop("'type' must be 'abs' (absolute differnce) or 'rat' (ratio)")
  }
}
get_pred_diff(lbp_intensity_nb, pred_df, type = "abs")
get_pred_diff(lbp_intensity_nb, pred_df, type = "rat")



fit_mod_boot <- function(data, i, pred_dat, type = "abs") {
  
  this_dat <- data[i, ]
  this_nbr <- glm.nb(mod_form_ilrs, data = this_dat)
  est <- get_pred_diff(this_nbr, new_dat = pred_dat, type = type)
  return(est)
  
}
alpha <- 0.05
quantile(boot(bpd_yes, fit_mod_boot, R = 10, pred_dat = pred_df)$t, c(alpha / 2, 1 - alpha / 2))


do_multi_realloc <- function(mod, basis_data, timeusenames, time_changes, sbp_matrix = sbp1) {
  
  
  x0 <- basis_data
  
  plot_dat <- 
    foreach(i = 1:length(timeusenames), .combine = bind_rows) %do% {
      print(paste("i: ", i))
      foreach(j = 1:length(timeusenames), .combine = bind_rows) %do% {
        print(paste("  j: ", j))
        foreach(d = 1:length(time_changes), .combine = bind_rows) %do% {
          print(paste("    d: ", d))
          
          timeuse_to <- timeusenames[i]
          timeuse_from <- timeusenames[j]
          change_time <- time_changes[d]
          
          proposed_change_1 <- x0[timeuse_to] + change_time
          proposed_change_2 <- x0[timeuse_from] - change_time
          
          if (timeuse_to == timeuse_from) {
            NULL # reallocation exceeds 0 or max time
          } else if ((proposed_change_1 < 0) | (proposed_change_1 > 1440)) {
            NULL # reallocation exceeds 0 or max time
          } else if ((proposed_change_2 < 0) | (proposed_change_2 > 1440)) {
            NULL # reallocation exceeds 0 or max time
          } else {
            
            x1 <- x0
            x1[timeuse_to] <- x1[timeuse_to] + change_time
            x1[timeuse_from] <- x1[timeuse_from] - change_time
            
            pred_df <- rbind(x0, x1)
            pred_df <- add_ilrs_to_data(pred_df, comp_vars = timeusenames, sbp_matrix = sbp_matrix)
            
            outc_ratio <- get_pred_diff(mod, pred_df)
            
            bootstrapped_ests <- boot(bpd_yes, fit_mod_boot, R = 1000, pred_dat = pred_df)$t
            ci_est <- quantile(as.numeric(bootstrapped_ests), c(alpha / 2, 1 - alpha / 2))
            
            tibble(
              to = timeuse_to, 
              from = timeuse_from, 
              change_time = change_time,
              outc_ratio = outc_ratio,
              ci_lo = ci_est[1],
              ci_hi = ci_est[2]
            )
          }
        }
      }
    }
  
  plot_dat$to <- factor(plot_dat$to, levels = timeusenames)
  plot_dat$from <- factor(plot_dat$from, levels = timeusenames)
  
  return(plot_dat)
  
}

set.seed(1234)


# takes ~60 min (single core) for bootstrapped CIs (R = 1000)
# takes ~ 6 min (single core) for bootstrapped CIs (R = 100)

### Uncomment to generate bootstrapping
# tic()
# realloc_plot_data <-
#   do_multi_realloc(
#     lbp_intensity_nb,
#     predict_basis0,
#     pred_comps,
#     seq(-30, 30, by = 10)
#   )
# saveRDS(realloc_plot_data, file = "res/negbin_realloc_boot_res(abs).rda")
# toc()



realloc_plot_data <- readRDS(file = "res/negbin_realloc_boot_res(abs).rda")


levels(realloc_plot_data$to) <- paste0(levels(realloc_plot_data$to), "+Delta")
levels(realloc_plot_data$from) <- paste0(levels(realloc_plot_data$from), "-Delta")

ggplot(realloc_plot_data) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 0, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3) +
  geom_line(aes(x = change_time , y = outc_ratio, col = to)) +
  geom_point(aes(x = change_time , y = outc_ratio, col = to), size = 1) +
  facet_grid(from ~ to, labeller = label_parsed) +
  theme_bw() +
  scale_colour_manual(values = c("darkorange","purple","cyan4", "dodgerblue")) +
  scale_fill_manual(values = c("darkorange","purple","cyan4", "dodgerblue")) +
  labs(
    x = paste0("Change/delta in composition (mins)"),
    y = paste0(
      "Absolute difference in mean LBP intensity (on % outcome scale, ",
      "after reallocation - before reallocation)")
  ) +
  theme(legend.position = "none")

ggsave(
  filename = "fig/lbp_intens_negbin_abs_v1.png",
  dpi = 600, # print quality
  width = 10,
  height = 10
)



pd2 <-
  realloc_plot_data %>%
  mutate(
    to = gsub("Time_", "", to),
    from = gsub("Time_", "", from),
    to = gsub("+Delta", "", to, fixed = TRUE),
    from = gsub("-Delta", "", from, fixed = TRUE),
    to_len = nchar(to),
    to_max = max(to_len),
    from_len = nchar(from),
    from_max = max(from_len),
    to_pad =  rep_char(pmax(0, from_max - to_len)),
    from_pad = rep_char(pmax(0, to_max - from_len)),
    to = factor(to, levels = time_lvls),
    from = factor(from, levels = time_lvls),
    to_num = as.numeric(to),
    from_num = as.numeric(from)
  ) %>%
  dplyr::filter(to_num > from_num) %>%
  mutate(
    # from_to = paste0("     ", "+", from, rep_char(10), from_pad, "\u2194", to_pad, rep_char(10), "+", to, "     ")
    from_to = paste0("+", from, rep_char(13), from_pad, "", to_pad, rep_char(13), "+", to)
  ) %>%
  arrange(from, to)

unique(pd2$from_to)
pd2$from_to <- factor(pd2$from_to, levels = unique(pd2$from_to))

this_breaks <- seq(-30, 30, 10)
this_labs <- sprintf("+%2.0f", abs(seq(-30, 30, 10)))
this_labs[this_labs == "+ 0"] <- ""
this_labs

ggplot(pd2) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 0, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3, col = NA, fill = "cyan4") +
  geom_line(aes(x = change_time , y = outc_ratio, col = to), col = "cyan4") +
  geom_point(aes(x = change_time , y = outc_ratio, col = to), size = 1, col = "cyan4") +
  facet_wrap(~ from_to, labeller = label_bquote(.(from_to))) +
  theme_bw() +
  scale_x_continuous(breaks = this_breaks, labels = this_labs) +
  labs(
    x = paste0("Reallocation between pair of compositional parts (minutes)"),
    y = paste0(
      "Absolute difference in mean LBP intensity (on % outcome scale, ",
      "after reallocation - before reallocation)"
    )
    # subtitle = "Note that odds ratios relate to the probability of having _increased_ frequency (per year) of pain"
  ) +
  theme(
    legend.position = "none",
    text = element_text(family = "serif"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text =  element_text(size = 10),
    axis.title =  element_text(size = 12)
  )


ggsave(filename = "fig/lbp_intens_negbin_abs_v2.png", width = 14, height = 9, dpi = 600)





# ---- outcome2_pred_b_rat ----



# wrapper:
get_pred_diff <- function(mod, new_dat, type = "rat") {
  if (type == "abs") {
    return(get_pred_diff_abs(mod = mod, new_dat = new_dat))
  } else if (type == "rat") {
    return(get_pred_diff_rat(mod = mod, new_dat = new_dat))
  } else {
    stop("'type' must be 'abs' (absolute differnce) or 'rat' (ratio)")
  }
}
get_pred_diff(lbp_intensity_nb, pred_df, type = "abs")
get_pred_diff(lbp_intensity_nb, pred_df, type = "rat")
get_pred_diff(lbp_intensity_nb, pred_df)



fit_mod_boot <- function(data, i, pred_dat, type = "rat") {
  
  this_dat <- data[i, ]
  this_nbr <- glm.nb(mod_form_ilrs, data = this_dat)
  est <- get_pred_diff(this_nbr, new_dat = pred_dat, type = type)
  return(est)
  
}
alpha <- 0.05
quantile(boot(bpd_yes, fit_mod_boot, R = 10, pred_dat = pred_df)$t, c(alpha / 2, 1 - alpha / 2))





# takes ~60 min (single core) for bootstrapped CIs (R = 1000)
# takes ~ 6 min (single core) for bootstrapped CIs (R = 100)


### Uncomment to generate bootstrapping
# set.seed(1234)
# tic()
# realloc_plot_data <-
#   do_multi_realloc(
#     lbp_intensity_nb,
#     predict_basis0,
#     pred_comps,
#     seq(-30, 30, by = 10)
#   )
# saveRDS(realloc_plot_data, file = "res/negbin_realloc_boot_res(rat).rda")
# toc()




realloc_plot_data <- readRDS(file = "res/negbin_realloc_boot_res(rat).rda")


levels(realloc_plot_data$to) <- paste0(levels(realloc_plot_data$to), "+Delta")
levels(realloc_plot_data$from) <- paste0(levels(realloc_plot_data$from), "-Delta")

ggplot(realloc_plot_data) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 1, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3) +
  geom_line(aes(x = change_time , y = outc_ratio, col = to)) +
  geom_point(aes(x = change_time , y = outc_ratio, col = to), size = 1) +
  facet_grid(from ~ to, labeller = label_parsed) +
  theme_bw() +
  scale_colour_manual(values = c("darkorange","purple","cyan4", "dodgerblue")) +
  scale_fill_manual(values = c("darkorange","purple","cyan4", "dodgerblue")) +
  labs(
    x = paste0("Change/delta in composition (mins)"),
    y = paste0("Ratio of back pain intensity (after reallocation:before reallocation)")
  ) +
  theme(legend.position = "none")

ggsave(
  filename = "fig/lbp_intens_negbin_rat_v1.png",
  dpi = 600, # print quality
  width = 10,
  height = 10
)



pd2 <-
  realloc_plot_data %>%
  mutate(
    to = gsub("Time_", "", to),
    from = gsub("Time_", "", from),
    to = gsub("+Delta", "", to, fixed = TRUE),
    from = gsub("-Delta", "", from, fixed = TRUE),
    to_len = nchar(to),
    to_max = max(to_len),
    from_len = nchar(from),
    from_max = max(from_len),
    to_pad =  rep_char(pmax(0, from_max - to_len)),
    from_pad = rep_char(pmax(0, to_max - from_len)),
    to = factor(to, levels = time_lvls),
    from = factor(from, levels = time_lvls),
    to_num = as.numeric(to),
    from_num = as.numeric(from)
  ) %>%
  dplyr::filter(to_num > from_num) %>%
  mutate(
    # from_to = paste0("     ", "+", from, rep_char(10), from_pad, "\u2194", to_pad, rep_char(10), "+", to, "     ")
    from_to = paste0("+", from, rep_char(13), from_pad, "", to_pad, rep_char(13), "+", to)
  ) %>%
  arrange(from, to)

unique(pd2$from_to)
pd2$from_to <- factor(pd2$from_to, levels = unique(pd2$from_to))

this_breaks <- seq(-30, 30, 10)
this_labs <- sprintf("+%2.0f", abs(seq(-30, 30, 10)))
this_labs[this_labs == "+ 0"] <- ""
this_labs

ggplot(pd2) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_hline(yintercept = 1, col = "grey60") +
  geom_ribbon(aes(x = change_time, ymin = ci_lo, ymax = ci_hi, fill = to), alpha = 0.3, col = NA, fill = "cyan4") +
  geom_line(aes(x = change_time , y = outc_ratio, col = to), col = "cyan4") +
  geom_point(aes(x = change_time , y = outc_ratio, col = to), size = 1, col = "cyan4") +
  facet_wrap(~ from_to, labeller = label_bquote(.(from_to))) +
  theme_bw() +
  scale_x_continuous(breaks = this_breaks, labels = this_labs) +
  labs(
    x = paste0("Reallocation between pair of compositional parts (minutes)"),
    y = paste0("Ratio of mean LBP intensity (after reallocation:before reallocation)")
    # subtitle = "Note that odds ratios relate to the probability of having _increased_ frequency (per year) of pain"
  ) +
  theme(
    legend.position = "none",
    text = element_text(family = "serif"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text =  element_text(size = 10),
    axis.title =  element_text(size = 12)
  )


ggsave(filename = "fig/lbp_intens_negbin_rat_v2.png", width = 14, height = 9, dpi = 600)










