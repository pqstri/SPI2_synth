# Import required libreries ----------------------------------------------------
require(tidyverse)
require(synthpop)
require(sdcMicro)
require(ggpubr)
require(broom.mixed)

# Load original data -----------------------------------------------------------
original.data <- "load original SPI2 data here"

# Synthetic data generation ----------------------------------------------------
## 1. Noise Addition  ####
# set a seed to make the synthesis reproducible
set.seed(1)

# set the desired amount of noise (in percentages)
noise.level <- 150 

# Add 10 times noise to bootstraped samples from the original sample
nadd.set <- lapply(1:10, function(i) {
  
  db_sampled <- db[sample(1:nrow(db), size = nrow(db), replace = TRUE),]
  
  fdb <- sdcMicro::addNoise(
    obj = as.data.frame(data.matrix(db_sampled)), 
    noise = noise.level, 
    method = "correlated")$xm
  
  # convert factor variable with noise back to the original scale
  handle_factors <- function(original, synth) {
    ref <- levels(original)[2]
    thrq <- (1-mean(original == ref))*100
    thr <- min(synth[ntile(synth, 100) > thrq])
    synth <- synth > thr
    synth <- factor(synth, labels = levels(original))
    return(synth)
  }
  
  # apply the conversion to the variables TREATMENT, MSFRM, SEX, and REGION
  binary.factor.vars <- c("TREATMENT", "MSFRM", "SEX", "REGION")
  for (v in binary.factor.vars) {
    fdb[,'v'] <- handle_factors(original.data[,'v'], fdb[,'v'])
  }
  
  return(fdb)
})

## 2. Chains of Conditional Distributions  #####
# set a seed to make the synthesis reproducible
set.seed(1)

# train the synthesizer and sample 10 times the original sample size
chcd.set <- synthpop::syn(db, m = 10)$syn

## 3. Multivariate modeling via Gaussian Copulas  #####
# Synthetic data are produced by GaussianCopulaSynthesizer of Python sdv library
# and exported from python as 10 csv files. The code to replicate this process 
# is provided in the gcs.py file, starting from the "original.data".
mvtn.set <- lapply(0:9, \(f) read.csv(sprintf("~/Downloads/temp%.csv", f)))


## 4. Generative Adversarial Networks ####
# Synthetic data are produced by CTGANSynthesizer of Python sdv library
# and exported from python as 10 csv files. The code to replicate this process 
# is provided in the ctgan.py file, starting from the "original.data".
gans.set <- lapply(0:9, \(f) read.csv(sprintf("~/Downloads/temp_ctgan%s.csv", f)))

# Structure data sets ----------------------------------------------------------
d <- list(
  "Original sample"                     = list(original.data),
  "Noise Addition"                      = nadd.set,
  "Chains of Conditional Distributions" = chcd.set, 
  "Multivariate modeling"               = mvtn.set, 
  "Generative Adversarial Networks"     = gans.set)

# Evaluation metrics -----------------------------------------------------------

## Proportion of Confidence Interval Overlap #####
# Function to Calculate Proportion of Confidence Interval Overlap
# Arguments:
# - comparison: A data frame or list containing the confidence interval 
# (conf.low and conf.high) of the comparison group.
# - reference: A data frame or list containing the confidence interval 
# (conf.low and conf.high) of the reference group.
# Returns:
# - The proportion of overlap between the two confidence intervals, ranging 
# from 0 to 1.
prop.ci.overlap <- function(comparison, reference) {
  
  # Calculate the absolute difference between lower bounds of confidence intervals
  exciding.side1 <- abs(reference$conf.low - comparison$conf.low)
  
  # Calculate the absolute difference between upper bounds of confidence intervals
  exciding.side2 <- abs(reference$conf.high - comparison$conf.high)
  
  # Calculate the total mismatch between confidence intervals
  mismatch <- exciding.side1 + exciding.side2
  
  # Calculate the range of the comparison group's confidence interval
  comparison.range <- abs(comparison$conf.low - comparison$conf.high)
  
  # Calculate the range of the reference group's confidence interval
  reference.range <- abs(reference$conf.low - reference$conf.high)
  
  # Calculate the total span of both confidence intervals
  total.span <- comparison.range + reference.range
  
  # Calculate the proportion of overlap or distance
  overlap_or_distance <- 1 - mismatch / total.span
  
  # Ensure overlap is non-negative
  overlap <- ifelse(overlap_or_distance < 0, 0, overlap_or_distance)
  
  # Return the calculated overlap proportion
  return(overlap)
}

## Standardized Mean Difference #####
# Function to Calculate Standardized Mean Difference (SMD)
# Arguments:
# - comparison: A data frame or list containing the estimate and confidence 
# interval (conf.low and conf.high) of the comparison group.
# - reference: A data frame or list containing the estimate and confidence
# interval (conf.low and conf.high) of the reference group.
# Returns:
# - The calculated Standardized Mean Difference.
smd <- function(comparison, reference) {
  
  # Calculate the standard error of the reference group
  reference.se <- (reference$conf.high - reference$conf.low) / (2 * 1.96)
  
  # Calculate the absolute distance between estimates
  absolute.distance <- abs(comparison$estimate - reference$estimate)
  
  # Calculate and return the Standardized Mean Difference
  return(absolute.distance / mean(reference.se))
}

# Summarizing tool -------------------------------------------------------------
summarise.q <- function(dat, ...) {
  dat %>% 
    group_by(model, term, ...) %>% 
    summarise_all(mean)
}

# Research scenarios -----------------------------------------------------------
## Q1 ####
q1.formula <- as.formula("edss_improved ~ TREATMENT + MSFRM + REGION")

q1.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    model <- glm(q1.formula,data = sample, family = "binomial")
    model.res <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
    model.res[-1, c("term", "estimate", "conf.low", "conf.high")]
    })
  })

q1.results.main <- q1.results[q1.results$term == "TREATMENTMD1003",]

## Q2 ####
q2.stat <- function(dat) {
  tst <- cor.test(dat$TW25BL, dat$EDSSBL)
  data.frame(
    term = "",
    conf.low = mean(tst$conf.int[1]),
    estimate = mean(tst$estimate),
    conf.high = mean(tst$conf.int[2])
  )
}

q2.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q2.stat(sample)
  })
})

## Q3 ----
q3.stat <- function(dat) {
  tst <- cor.test(dat$TW25M15 - dat$TW25BL, 
                  dat$EDSSM15 - dat$EDSSBL)
  data.frame(
    term = "",
    conf.low = mean(tst$conf.int[1]),
    estimate = mean(tst$estimate),
    conf.high = mean(tst$conf.int[2])
  )
}

q3.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q3.stat(sample)
  })
})

## Q4 ----
q4.stat <- function(dat) {
  m <- t.test((dat$TW25M15 - dat$TW25BL) ~ dat$MSFRM)
  data.frame(
    term = "",
    conf.low = m$conf.int[1],
    estimate = unname(-diff(m$estimate)),
    conf.high = m$conf.int[2]
  )
}

q4.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q4.stat(sample)
  })
})

## Q5 ----
lonitudify_tw25 <- function(dat) {
  dat %>% 
    select(contains("tw25"), -contains("imp")) %>% 
    rownames_to_column("id") %>% 
    gather("time", "TW25", -id) %>% 
    mutate(time = case_when(
      time == "TW25MSCR" ~ -2,
      time == "TW25BL"   ~ -1,
      time == "TW25M0"   ~ 0,
      time == "TW25M3"   ~ 3,
      time == "TW25M6"   ~ 6,
      time == "TW25M9"   ~ 9,
      time == "TW25M12"  ~ 12,
      time == "TW25M15"  ~ 15))
}

q5.stat <- function(dat) {
  m <- lme4::lmer(TW25 ~ time + (1|id), data = lonitudify_tw25(dat))
  b <- broom.mixed::tidy(m, conf.int=T)
  b <- b[b$term == "time",]
  data.frame(
    term = "",
    conf.low = b$conf.low,
    estimate = b$estimate,
    conf.high = b$conf.high
  )
}

q5.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q5.stat(sample)
  })
})

## Q6 ####
q6.stat.generator <- function(formula, extractor) {
  return(function(dat) {
    m <- lm(formula, data = dat)
    b <- broom::tidy(m, conf.int = TRUE)
    b <- b[b$term == extractor,]
    data.frame(
      term = "",
      conf.low = b$conf.low,
      estimate = b$estimate,
      conf.high = b$conf.high
    )
  })
}

### Step 1 ----
q6_1.stat <- q6.stat.generator(
  formula = as.formula("I(TW25M0 - TW25M6) ~ BMI"),
  extractor = "BMI"
)

q6_1.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q6_1.stat(sample)
  })
})

### Step 2 ----
q6_2.stat <- q6.stat.generator(
  formula = as.formula("I(TW25M0 - TW25M6) ~ SEX + AGE + BMI"),
  extractor = "BMI"
)

q6_2.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q6_2.stat(sample)
  })
})

### Step 3 ----
q6_3.stat <- q6.stat.generator(
  formula = as.formula("I(TW25M0 - TW25M6) ~ SEX + AGE + BMI * ntile(EDSSM0, 2)"),
  extractor = "BMI:ntile(EDSSM0, 2)"
)

q6_3.results <- map_df(d, .id = "model", function(set) {
  set.results <- map_df(set, function(sample) {
    q6_3.stat(sample)
  })
})

# Figure 1 ---------------------------------------------------------------------
Qs <- list(
  "Q1" = q1.results.main, 
  "Q2" = q2.results, 
  "Q3" = q3.results, 
  "Q4" = q4.results, 
  "Q5" = q5.results, 
  "Q6" = q6_3.results) %>% 
  bind_rows(.id = "q") %>% 
  group_by(q) %>% 
  summarise.q(q) %>% 
  mutate(ref = ifelse(model == "Original sample", estimate, NA)) %>% 
  mutate(refl = ifelse(model == "Original sample", conf.low, NA)) %>% 
  mutate(refh = ifelse(model == "Original sample", conf.high, NA)) %>% 
  ungroup()

Qs <- Qs %>% 
  mutate(model = factor(case_when(
    model == "Original sample"                     ~ "Original sample",
    model == "Noise Addition"                      ~ "Noise addition",
    model == "Chains of Conditional Distributions" ~ "Conditional distributions",
    model == "Multivariate modeling"               ~ "Multivariate modeling",
    model == "Generative Adversarial Networks"     ~ "Generative network"), 
    ordered = TRUE, 
    levels = c("Original sample", 
               "Noise addition", 
               "Conditional distributions", 
               "Multivariate modeling", 
               "Generative network")))

ggplot(Qs, aes(estimate, model, xmin = conf.low, xmax = conf.high)) + 
  facet_wrap(~q, scales = "free_x") +
  geom_rect(aes(fill = model), ymin = -1, ymax = 6, alpha = 0.1) +
  geom_vline(aes(xintercept = ifelse(q == "Q1", 1, 0)), lty = 1) +
  geom_vline(aes(xintercept = ref), lty = 2, col = 2) +
  geom_pointrange(aes(shape = model, col = model), size = 0.8) +
  ggpubr::theme_pubr(base_size = 8, border = TRUE) +
  theme(plot.margin = margin(t = 0.5, r = 1.5, b = 0.5, l = 0.5, unit = "cm"),
        strip.text = element_text(face = "bold")) +
  scale_shape_manual(guide = "none", values = c(18, rep(20, 4))) +
  scale_color_manual(guide = "none", values = c(2, rep(1, 4))) +
  scale_fill_manual(guide = "none", values = c(2, rep("transparent", 4))) +
  labs(x = "\nEstimate", y = "")

# ggsave("~/fig1.jpg", dpi = 400, width = 10, height = 6, units = "cm", scale = 2)